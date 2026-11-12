"""Master per-event merge + descriptive summaries (no inferential statistics).

Combines every ``trajectories_trial{N}_Segment{M}_events.csv`` produced by
the ``assign`` step into a single tidy table and emits two flavours of
per-fish summary:

* ``Attacker_Summary_segment`` / ``Attacked_Summary_segment`` -- keyed on
  ``(trial, segment, fish_id_segment)``. Always valid.
* ``Attacker_Summary_matched`` / ``Attacked_Summary_matched`` -- keyed on
  ``fish_id_matched_trial``. Rows where the matched ID is NA are excluded
  from the matched-key summaries; the excluded count is logged.

Statistical tests are deliberately omitted (see plan section 6).
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def merge_event_csvs(results_dir: Path) -> pd.DataFrame:
    """Concatenate every ``*_events.csv`` under ``results_dir``."""
    if not results_dir.exists():
        raise FileNotFoundError(f"RESULTS_DIR does not exist: {results_dir}")
    csvs = sorted(results_dir.glob("*_events.csv"))
    if not csvs:
        raise FileNotFoundError(f"No *_events.csv files found in {results_dir}")
    parts = []
    for p in csvs:
        try:
            chunk = pd.read_csv(p)
            if not chunk.empty:
                parts.append(chunk)
        except Exception:
            # Skip truly-empty or malformed files (e.g. written by the assign
            # step when a trial has no events after scope filtering).
            pass
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    logger.info("Merged %d files, %d events total", len(csvs), len(df))
    return df


def make_segment_keyed_summary(
    df: pd.DataFrame, role: str
) -> pd.DataFrame:
    """Aggregate events by ``(trial, segment, <role>_id)``."""
    role_id = f"{role}_id"
    if role_id not in df.columns or df.empty:
        return pd.DataFrame()
    grouped = (
        df.groupby(["trial", "segment", role_id], dropna=False)
          .agg(
              n_events=(role_id, "size"),
              mean_attacker_accel_px_per_frame2=(
                  "attacker_max_acceleration_px_per_frame2", "mean"
              ),
              sd_attacker_accel_px_per_frame2=(
                  "attacker_max_acceleration_px_per_frame2", "std"
              ),
              mean_attacker_accel_BL_per_s2=(
                  "attacker_max_acceleration_BL_per_s2", "mean"
              ),
              sd_attacker_accel_BL_per_s2=(
                  "attacker_max_acceleration_BL_per_s2", "std"
              ),
              mean_distance_px=("distance_at_endpoint_px", "mean"),
              mean_distance_BL=("distance_at_endpoint_BL", "mean"),
          )
          .reset_index()
          .rename(columns={role_id: "fish_id_segment"})
    )
    return grouped


def make_matched_keyed_summary(
    df: pd.DataFrame, role: str
) -> tuple[pd.DataFrame, int]:
    """Aggregate events by ``<role>_fish_id_matched_trial``.

    Returns
    -------
    tuple
        ``(summary_df, n_excluded)`` where ``n_excluded`` is the number of
        events whose matched ID is NA (and therefore not in the summary).
    """
    key = f"{role}_fish_id_matched_trial"
    if key not in df.columns or df.empty:
        return pd.DataFrame(), 0
    mask = df[key].notna()
    n_excluded = int((~mask).sum())
    sub = df.loc[mask].copy()
    if sub.empty:
        return pd.DataFrame(), n_excluded
    grouped = (
        sub.groupby([key], dropna=False)
           .agg(
               trial=("trial", "first"),
               n_events=(key, "size"),
               mean_attacker_accel_BL_per_s2=(
                   "attacker_max_acceleration_BL_per_s2", "mean"
               ),
               sd_attacker_accel_BL_per_s2=(
                   "attacker_max_acceleration_BL_per_s2", "std"
               ),
               mean_distance_BL=("distance_at_endpoint_BL", "mean"),
           )
           .reset_index()
           .rename(columns={key: "fish_id_matched_trial"})
    )
    return grouped, n_excluded


def _build_side_coverage(
    df: pd.DataFrame,
    side_diag: pd.DataFrame | None,
    skipped_by_segment: pd.DataFrame | None,
) -> pd.DataFrame:
    """Per-trial table summarising side-annotation coverage and skip reasons."""
    rows: list[dict] = []

    # Per-trial counts of ok / missing / ambiguous from side_diag.
    diag_by_trial: dict[int, dict[str, int]] = {}
    if side_diag is not None and not side_diag.empty:
        # Derive trial from the observation id (trialN_...).
        def _trial_from_obs(obs: str) -> int | None:
            try:
                if isinstance(obs, str) and obs.lower().startswith("trial"):
                    return int(obs.split("_", 1)[0][5:])
            except Exception:
                return None
            return None

        side_diag = side_diag.copy()
        side_diag["_trial"] = side_diag["Observation id"].map(_trial_from_obs)
        for trial, sub in side_diag.groupby("_trial"):
            if pd.isna(trial):
                continue
            diag_by_trial[int(trial)] = {
                "n_side_ok": int((sub["side_status"] == "ok").sum()),
                "n_side_missing": int((sub["side_status"] == "missing").sum()),
                "n_side_ambiguous": int((sub["side_status"] == "ambiguous").sum()),
            }

    # Per-trial counts of side-related skip reasons from skipped_by_segment.
    skip_by_trial: dict[int, dict[str, int]] = {}
    if skipped_by_segment is not None and not skipped_by_segment.empty:
        side_cols = [
            c for c in (
                "side_missing", "side_ambiguous", "side_no_candidates",
                "side_one_candidate", "attacker_position_nan_at_event",
            ) if c in skipped_by_segment.columns
        ]
        if side_cols:
            grouped = (
                skipped_by_segment.groupby("trial")[side_cols].sum().astype(int)
            )
            for trial, row in grouped.iterrows():
                skip_by_trial[int(trial)] = {c: int(row[c]) for c in side_cols}

    # Per-trial counts of processed events grouped by attacker side.
    proc_by_trial: dict[int, dict[str, int]] = {}
    if not df.empty and {"trial", "annotated_side"}.issubset(df.columns):
        for trial, sub in df.groupby("trial"):
            n_legacy = (
                int((~sub["side_filter_applied"].fillna(True)).sum())
                if "side_filter_applied" in sub.columns else 0
            )
            proc_by_trial[int(trial)] = {
                "n_processed_total": len(sub),
                "n_processed_flow": int((sub["annotated_side"] == "flow").sum()),
                "n_processed_calm": int((sub["annotated_side"] == "calm").sum()),
                "n_processed_legacy_fallback": n_legacy,
            }

    all_trials = sorted(
        set(diag_by_trial) | set(skip_by_trial) | set(proc_by_trial)
    )
    for trial in all_trials:
        row = {"trial": trial}
        row.update(diag_by_trial.get(trial, {}))
        row.update(proc_by_trial.get(trial, {}))
        row.update(skip_by_trial.get(trial, {}))
        rows.append(row)

    return pd.DataFrame(rows)


def _build_annotated_vs_trajectory_side(df: pd.DataFrame) -> pd.DataFrame:
    """Crosstab annotated_side vs attacker_side_at_event (should be diagonal).

    Restricted to events where the side filter was actually applied
    (``side_filter_applied == True``).  Legacy-fallback events have
    ``annotated_side`` = NA and are excluded so the 100 % diagonal invariant
    remains meaningful.
    """
    if df.empty or not {"annotated_side", "attacker_side_at_event"}.issubset(df.columns):
        return pd.DataFrame()
    # Exclude legacy-fallback events (annotated_side is NA for those)
    sub = df[df["annotated_side"].notna()].copy()
    if sub.empty:
        return pd.DataFrame()
    return (
        pd.crosstab(
            sub["annotated_side"], sub["attacker_side_at_event"], dropna=False,
        )
        .reset_index()
    )


def write_master_xlsx(
    df: pd.DataFrame,
    out_path: Path,
    scope_banner: str,
    side_diag: pd.DataFrame | None = None,
    skipped_by_segment: pd.DataFrame | None = None,
) -> dict[str, int]:
    """Write the master Excel with descriptive summaries only.

    The side-aware variant adds two extra sheets:
      * ``Side_Coverage`` -- per-trial side-annotation OK / missing / ambiguous
        counts plus side-related skip counts.
      * ``Annotated_vs_Trajectory_Side`` -- crosstab; expected 100 % diagonal.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    info: dict[str, int] = {}

    with pd.ExcelWriter(out_path, engine="openpyxl") as xl:
        df.to_excel(xl, sheet_name="Aggression_Events", index=False)

        meta = pd.DataFrame({
            "key": [
                "generated_at", "scope_banner",
                "n_events", "n_trials", "n_segments",
            ],
            "value": [
                datetime.now().isoformat(timespec="seconds"),
                scope_banner,
                len(df),
                int(df["trial"].nunique()) if "trial" in df.columns else 0,
                int(df["segment"].nunique()) if "segment" in df.columns else 0,
            ],
        })
        meta.to_excel(xl, sheet_name="Run_Info", index=False)

        # Counts
        if "interaction_type" in df.columns:
            (
                df["interaction_type"]
                .value_counts(dropna=False)
                .rename_axis("interaction_type")
                .reset_index(name="count")
                .to_excel(xl, sheet_name="Interaction_Types", index=False)
            )

        for role in ("attacker", "attacked"):
            seg_summary = make_segment_keyed_summary(df, role)
            if not seg_summary.empty:
                seg_summary.to_excel(
                    xl,
                    sheet_name=f"{role.capitalize()}_Summary_segment",
                    index=False,
                )

            matched_summary, n_excluded = make_matched_keyed_summary(df, role)
            info[f"{role}_n_excluded_from_matched"] = n_excluded
            if not matched_summary.empty:
                matched_summary.to_excel(
                    xl,
                    sheet_name=f"{role.capitalize()}_Summary_matched",
                    index=False,
                )

        # Crosstabs
        if {"treatment", "condition"}.issubset(df.columns):
            crosstab = (
                pd.crosstab(df["treatment"], df["condition"], dropna=False)
                .reset_index()
            )
            crosstab.to_excel(xl, sheet_name="Treatment_x_Condition", index=False)
        if {"treatment", "interaction_type"}.issubset(df.columns):
            crosstab = (
                pd.crosstab(df["treatment"], df["interaction_type"], dropna=False)
                .reset_index()
            )
            crosstab.to_excel(xl, sheet_name="Type_by_Treatment", index=False)

        # Side-aware diagnostic sheets
        side_cov = _build_side_coverage(df, side_diag, skipped_by_segment)
        if not side_cov.empty:
            side_cov.to_excel(xl, sheet_name="Side_Coverage", index=False)
        cross_side = _build_annotated_vs_trajectory_side(df)
        if not cross_side.empty:
            cross_side.to_excel(
                xl, sheet_name="Annotated_vs_Trajectory_Side", index=False,
            )

    info["n_events"] = len(df)
    return info


def derive_master_filename(
    base_name: str,
    tp1_switch: bool,
    drop_pharma: bool,
    side_aware: bool = False,
    hybrid: bool = False,
) -> str:
    suffix = ""
    if tp1_switch:
        suffix += "_noTP1"
    if drop_pharma:
        suffix += "_noPharma"
    if side_aware:
        suffix += "_side"
    if hybrid:
        suffix += "_hybrid"
    return f"{base_name}{suffix}.xlsx"
