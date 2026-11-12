"""Construct the trial-scoped cross-segment identity map from idmatcher.ai output.

See plan section 10 for the full integration rules.

The canonical input is ``datasets/idmatcher_annotated.csv`` produced by
``01_5_annotate_idmatcher_GOT.R``. Each row links a Segment-2 anchor fish to a
candidate Segment-3 fish with a ``match_quality`` label and a boolean
``match_accepted`` flag. From accepted rows we build:

* a long table keyed on ``(trial, segment, fish_id_segment)`` carrying
  ``fish_id_matched_trial`` of the form ``T<trial>_MF<n>`` (matching the R
  helper's convention), plus ``match_quality`` and ``match_accepted``.

Segment 1 is never linked to a matched identity (idmatcher only matches
seg2 <-> seg3). Trials absent from the idmatcher CSV propagate as NA matched
IDs throughout the pipeline -- never silently filled.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


IDMATCHER_REQUIRED_COLUMNS = {
    "trial", "fish_id_master", "master_segment", "segment",
    "fish_id", "direct_score", "indirect_score",
    "match_quality", "match_accepted",
}


def load_idmatcher_annotated(path: Path) -> pd.DataFrame:
    """Read ``idmatcher_annotated.csv`` with strict column checking."""
    if not path.exists():
        raise FileNotFoundError(
            f"IDMATCHER_ANNOTATED_CSV not found at {path}. "
            "Run 01_5_annotate_idmatcher_GOT.R, or set "
            "config.IDMATCHER_ANNOTATED_CSV = None to disable the cross-segment join."
        )
    df = pd.read_csv(path)
    missing = IDMATCHER_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {sorted(missing)}"
        )
    # Normalise types
    df["trial"] = df["trial"].astype(int)
    df["fish_id_master"] = df["fish_id_master"].astype(int)
    df["master_segment"] = df["master_segment"].astype(int)
    df["segment"] = df["segment"].astype(int)
    df["fish_id"] = df["fish_id"].astype(int)
    df["match_accepted"] = df["match_accepted"].astype(bool)
    df["match_quality"] = df["match_quality"].astype(str)
    return df


def build_idmatcher_map(annotated: pd.DataFrame) -> pd.DataFrame:
    """Produce the ``(trial, segment, fish_id_segment) -> matched_trial`` map.

    Returns
    -------
    pd.DataFrame
        Columns: ``trial`` (int), ``segment`` (int),
        ``fish_id_segment`` (str), ``fish_id_matched_trial`` (str),
        ``match_quality`` (str), ``match_accepted`` (bool).

    Notes
    -----
    * Seg2 anchors get ``fish_id_matched_trial`` regardless of acceptance --
      the anchor identity always exists.
    * Seg3 candidates get ``fish_id_matched_trial`` only when the row was
      accepted.
    * Rejected seg3 rows are still emitted with ``fish_id_matched_trial = NA``
      and the original ``match_quality`` (e.g. ``rejected_conflict``).
    """
    rows: list[dict[str, Any]] = []
    seg2_seen: set[tuple[int, int]] = set()

    for _, r in annotated.iterrows():
        trial = int(r["trial"])
        master_seg = int(r["master_segment"])
        seg = int(r["segment"])
        fid_master = int(r["fish_id_master"])
        fid_seg = int(r["fish_id"])
        accepted = bool(r["match_accepted"])
        quality = str(r["match_quality"])
        matched_id = f"T{trial}_MF{fid_master}"

        # Seg2 anchor: always known
        key = (trial, master_seg, fid_master)
        if key not in seg2_seen:
            rows.append({
                "trial": trial,
                "segment": master_seg,
                "fish_id_segment": str(fid_master),
                "fish_id_matched_trial": matched_id,
                "match_quality": "anchor",
                "match_accepted": True,
            })
            seg2_seen.add(key)

        # Seg3 candidate: matched only if accepted
        rows.append({
            "trial": trial,
            "segment": seg,
            "fish_id_segment": str(fid_seg),
            "fish_id_matched_trial": matched_id if accepted else pd.NA,
            "match_quality": quality,
            "match_accepted": accepted,
        })

    df = pd.DataFrame(rows)

    # Detect ambiguous many-to-one: same (trial, segment, fish_id_segment)
    # carrying two different non-NA matched IDs.
    dup = (
        df.dropna(subset=["fish_id_matched_trial"])
          .groupby(["trial", "segment", "fish_id_segment"])
          ["fish_id_matched_trial"].nunique()
    )
    bad = dup[dup > 1]
    if not bad.empty:
        raise ValueError(
            "idmatcher map is ambiguous (same segment-local fish mapped "
            f"to multiple matched IDs):\n{bad}"
        )

    # If both an accepted and a rejected row exist for the same key, keep the
    # accepted one (matched_trial wins over NA).
    df = (
        df.sort_values(
            by=["trial", "segment", "fish_id_segment", "match_accepted"],
            ascending=[True, True, True, False],
        )
        .drop_duplicates(
            subset=["trial", "segment", "fish_id_segment"], keep="first"
        )
        .reset_index(drop=True)
    )
    return df


def covered_trials(annotated: pd.DataFrame) -> list[int]:
    """Return sorted list of trial numbers present in the idmatcher CSV."""
    return sorted(int(t) for t in annotated["trial"].unique())


# ---------------------------------------------------------------------------
# Body-length diagnostics (optional)
# ---------------------------------------------------------------------------

BL_DIAGNOSTICS_OPTIONAL_COLUMNS = (
    "trial", "fish_id_master",
    "body_length_seg2_cm", "body_length_seg3_cm",
    "body_length_mean_cm", "body_length_source",
)


def load_body_length_diagnostics(path: Path | None) -> pd.DataFrame | None:
    """Load optional ``body_length_diagnostics_*.csv`` from the R STEP 1.5.

    Accepts either of two schemas:

    * Legacy: columns ``trial``, ``fish_id_master``, ``body_length_mean_cm``.
    * Current (May 2026 helper): columns ``trial_id``, ``fish_id_matched_trial``
      (e.g. ``T5_MF1``), ``body_length_seg2_cm``, ``body_length_seg3_cm``,
      ``body_length_mean_cm`` (NaN where stability checks failed).

    For the current schema this function (a) renames ``trial_id`` -> ``trial``,
    (b) parses ``fish_id_matched_trial`` -> ``fish_id_master`` (the integer
    after ``_MF``), and (c) when ``body_length_mean_cm`` is NaN, fills it with
    the mean of the available per-segment values so downstream events get a
    body length even for fish flagged as unstable.
    """
    if path is None:
        return None
    if not path.exists():
        logger.info(
            "BODY_LENGTH_DIAGNOSTICS_CSV not found at %s; "
            "body-length-mean-cm columns will be NA", path,
        )
        return None
    df = pd.read_csv(path)

    # Adapt the current-schema column names to what the rest of the pipeline
    # expects. The legacy schema is left untouched.
    if "trial" not in df.columns and "trial_id" in df.columns:
        df = df.rename(columns={"trial_id": "trial"})
    if "fish_id_master" not in df.columns and "fish_id_matched_trial" in df.columns:
        df["fish_id_master"] = (
            df["fish_id_matched_trial"]
            .astype(str)
            .str.extract(r"_MF(\d+)$")[0]
        )

    # Fallback: if body_length_mean_cm is NaN but per-segment values exist,
    # use the mean of whatever is available. The diagnostic was leaving the
    # mean blank for fish flagged as unstable, but the per-segment values
    # are still informative for body-length normalisation.
    if "body_length_mean_cm" in df.columns:
        seg_cols = [c for c in ("body_length_seg2_cm", "body_length_seg3_cm")
                    if c in df.columns]
        if seg_cols:
            fallback = df[seg_cols].mean(axis=1, skipna=True)
            mask = df["body_length_mean_cm"].isna() & fallback.notna()
            n_filled = int(mask.sum())
            if n_filled:
                df.loc[mask, "body_length_mean_cm"] = fallback.loc[mask]
                logger.info(
                    "body_length_mean_cm: filled %d NaN rows via mean of "
                    "available per-segment values.", n_filled,
                )

    have = set(df.columns)
    missing = [c for c in BL_DIAGNOSTICS_OPTIONAL_COLUMNS if c not in have]
    if missing:
        logger.warning(
            "%s is missing optional columns %r; "
            "body-length enrichment will be partial", path, missing,
        )
    df["trial"] = df["trial"].astype(int)
    if "fish_id_master" in df.columns:
        df["fish_id_master"] = df["fish_id_master"].astype(int)
    return df


def attach_body_length_cm(
    events: pd.DataFrame,
    bl_diag: pd.DataFrame | None,
) -> pd.DataFrame:
    """Add attacker_/attacked_body_length_mean_cm columns to ``events``.

    Primary join: exact match on (trial, fish_number) extracted from the
    matched-trial ID string (e.g. 'T15_MF1' → fish_number '1').

    Fallback: fish that have a matched trial ID but are absent from the body-
    length CSV (e.g. too unstable for measurement in any segment) receive the
    mean body length of all measured fish in the same trial.  This imputed
    value is used only when no individual estimate is available.
    """
    out = events.copy()
    if bl_diag is None or "body_length_mean_cm" not in bl_diag.columns:
        out["attacker_body_length_mean_cm"] = pd.NA
        out["attacked_body_length_mean_cm"] = pd.NA
        return out

    lookup = bl_diag[["trial", "fish_id_master", "body_length_mean_cm"]].copy()
    lookup["fish_id_master"] = lookup["fish_id_master"].astype(str)
    lookup = lookup.rename(columns={
        "fish_id_master": "matched_fish_no",
        "body_length_mean_cm": "_bl_mean_cm",
    })

    # Trial-mean fallback for fish absent from the body-length CSV.
    trial_mean_bl = (
        lookup.groupby("trial")["_bl_mean_cm"]
        .mean()
        .rename("_trial_mean_bl")
        .reset_index()
    )

    # The matched ID format is T<trial>_MF<n>; extract <n>
    def _fish_no(s: Any) -> Any:
        if s is None or s is pd.NA or (isinstance(s, float) and pd.isna(s)):
            return pd.NA
        s = str(s)
        if "_MF" not in s:
            return pd.NA
        return s.split("_MF", 1)[1]

    out["_atk_n"] = out["attacker_fish_id_matched_trial"].map(_fish_no)
    out["_vct_n"] = out["attacked_fish_id_matched_trial"].map(_fish_no)

    out = out.merge(
        lookup.rename(columns={
            "matched_fish_no": "_atk_n",
            "_bl_mean_cm": "attacker_body_length_mean_cm",
        }),
        on=["trial", "_atk_n"], how="left",
    )
    out = out.merge(
        lookup.rename(columns={
            "matched_fish_no": "_vct_n",
            "_bl_mean_cm": "attacked_body_length_mean_cm",
        }),
        on=["trial", "_vct_n"], how="left",
    )

    # Fill remaining NAs (fish with matched ID but absent from CSV) with
    # the trial mean body length.
    out = out.merge(trial_mean_bl, on="trial", how="left")
    for role, id_col, bl_col in (
        ("attacker", "attacker_fish_id_matched_trial", "attacker_body_length_mean_cm"),
        ("attacked", "attacked_fish_id_matched_trial", "attacked_body_length_mean_cm"),
    ):
        fill_mask = out[id_col].notna() & out[bl_col].isna()
        n_filled = int(fill_mask.sum())
        if n_filled:
            out.loc[fill_mask, bl_col] = out.loc[fill_mask, "_trial_mean_bl"]
            logger.info(
                "attach_body_length_cm: %d %s rows absent from body-length CSV; "
                "filled with trial-mean body length.",
                n_filled, role,
            )

    return out.drop(columns=["_atk_n", "_vct_n", "_trial_mean_bl"])
