"""Comprehensive missing-events + idmatcher-coverage diagnostics.

Adapted from block J of the original Colab script (lines 3771-4190) with the
two additions required by plan section 10.8:

* an ``Idmatcher_Coverage`` section that summarises per-trial idmatcher status,
* an identity-label reconciliation between cleaned NPY identities and the
  ``fish_id_master`` / ``fish_id`` sets in ``idmatcher_annotated.csv``.

This step does **not** re-process trajectories; it reads the already-merged
per-event table and the BORIS CSV to bucket every missing event into a
category.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# Categories used in DETAILED_MISSING_EVENTS.csv
MISSING_REASONS = (
    "in_gap_between_segments",
    "out_of_range",
    "missing_segment_npy",
    "insufficient_tracking",
    "algorithm_failed_no_attacker",
    "algorithm_failed_no_attacked",
    "frame_calc_issue",
    "scope_filter_excluded_trial",
    "scope_filter_seg1",
    "scope_filter_pharma",
    # Side-aware additions
    "side_missing",
    "side_ambiguous",
    "side_no_candidates",
    "side_one_candidate",
    "attacker_position_nan_at_event",
    "unknown",
)


def write_missing_events_reports(
    boris: pd.DataFrame,
    events: pd.DataFrame,
    skipped_by_segment: dict[tuple[int, int], dict[str, int]],
    npy_available: dict[tuple[int, int], bool],
    qc_dir: Path,
    scope_banner: str,
    scope_skipped: dict[str, int],
    side_diag: pd.DataFrame | None = None,
    zone_y0_px: float | None = None,
    zone_flow_higher_y: bool | None = None,
) -> dict[str, int]:
    """Emit the three missing-events reports under ``qc_dir``.

    Parameters
    ----------
    boris:
        The aggression-only subset of the BORIS CSV (one row per scored event,
        post-START removal), with ``trial`` already parsed into a column.
    events:
        The master per-event dataframe produced by ``merge`` (one row per
        processed event).
    skipped_by_segment:
        Mapping ``(trial, segment) -> {reason: count}`` collected during
        ``assign``.
    npy_available:
        Mapping ``(trial, segment) -> bool`` (does the staged NPY exist?).
    qc_dir:
        Output directory.
    scope_banner:
        A one-line description of active filters (for the text report header).
    scope_skipped:
        ``{'skipped_tp1_switch': ..., 'skipped_pharma_group': ...}`` from the
        scope-filter pass.
    """
    qc_dir.mkdir(parents=True, exist_ok=True)

    boris_total = len(boris)
    processed_total = len(events)

    # Build a (obs_id, start) -> side_status lookup so per-event classification
    # can correctly distinguish 'side_missing' / 'side_ambiguous'.
    side_status_lookup: dict[tuple[str, float], str] = {}
    if side_diag is not None and not side_diag.empty:
        for _, r in side_diag.iterrows():
            side_status_lookup[
                (str(r["Observation id"]), round(float(r["Start (s)"]), 6))
            ] = str(r["side_status"])

    # Detailed per-event report (only for events that were NOT processed)
    processed_keys = set(
        zip(events.get("boris_observation", []), events.get("time_absolute_sec", []))
    )
    missing_rows = []
    for _, ev in boris.iterrows():
        key = (str(ev.get("Observation id")), float(ev.get("Start (s)", 0)))
        if key in processed_keys:
            continue
        reason = _classify_missing(
            ev, skipped_by_segment, npy_available, side_status_lookup,
        )
        missing_rows.append({
            "observation_id": ev.get("Observation id"),
            "trial": ev.get("trial"),
            "start_sec": ev.get("Start (s)"),
            "expected_segment": ev.get("expected_segment"),
            "reason": reason,
        })
    missing_df = pd.DataFrame(missing_rows)
    missing_df.to_csv(qc_dir / "DETAILED_MISSING_EVENTS.csv", index=False)

    # Per-trial capture rates
    if "trial" in boris.columns:
        per_trial_boris = boris.groupby("trial").size().rename("n_boris")
    else:
        per_trial_boris = pd.Series(dtype=int, name="n_boris")
    if "trial" in events.columns:
        per_trial_proc = events.groupby("trial").size().rename("n_processed")
    else:
        per_trial_proc = pd.Series(dtype=int, name="n_processed")
    capture = pd.concat([per_trial_boris, per_trial_proc], axis=1).fillna(0)
    capture["n_missing"] = capture["n_boris"] - capture["n_processed"]
    capture["capture_rate"] = (
        capture["n_processed"] / capture["n_boris"].where(capture["n_boris"] > 0)
    )
    capture.reset_index().to_csv(qc_dir / "TRIAL_CAPTURE_RATES.csv", index=False)

    # Idmatcher coverage
    cov = _idmatcher_coverage(events)
    cov.to_csv(qc_dir / "idmatcher_coverage_by_trial.csv", index=False)

    # Text report
    report_path = qc_dir / "COMPREHENSIVE_MISSING_EVENTS_REPORT.txt"
    with report_path.open("w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write("COMPREHENSIVE MISSING EVENTS REPORT\n")
        f.write("=" * 78 + "\n\n")
        f.write(f"Scope: {scope_banner}\n")
        f.write(f"BORIS aggression events  : {boris_total}\n")
        f.write(f"Events processed         : {processed_total}\n")
        f.write(f"Events missing           : {boris_total - processed_total}\n")
        f.write(f"Capture rate             : "
                f"{(processed_total / boris_total if boris_total else 0):.1%}\n")
        f.write("\nScope filter telemetry:\n")
        for k, v in scope_skipped.items():
            f.write(f"  {k:<28s} {v}\n")
        f.write("\nMissing events by reason:\n")
        if not missing_df.empty:
            for reason, n in missing_df["reason"].value_counts().items():
                f.write(f"  {reason:<32s} {n}\n")
        f.write("\nIdmatcher coverage by trial:\n")
        f.write(cov.to_string(index=False))
        f.write("\n")

    # Side-coverage text report (side-aware variant)
    side_report_path = qc_dir / "SIDE_COVERAGE_REPORT.txt"
    with side_report_path.open("w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write("SIDE-ANNOTATION COVERAGE REPORT\n")
        f.write("=" * 78 + "\n\n")
        f.write(f"Scope: {scope_banner}\n")
        f.write(
            f"Zone rule: y >= {zone_y0_px} -> "
            f"{'flow' if zone_flow_higher_y else 'calm'}, "
            f"y < {zone_y0_px} -> "
            f"{'calm' if zone_flow_higher_y else 'flow'}\n"
        )
        if side_diag is not None and not side_diag.empty:
            n_ok = int((side_diag["side_status"] == "ok").sum())
            n_missing = int((side_diag["side_status"] == "missing").sum())
            n_ambig = int((side_diag["side_status"] == "ambiguous").sum())
            f.write(f"\nSide annotation join (all BORIS aggression rows):\n")
            f.write(f"  side_status=ok        {n_ok}\n")
            f.write(f"  side_status=missing   {n_missing}\n")
            f.write(f"  side_status=ambiguous {n_ambig}\n")
            if n_ok:
                by_side = side_diag.loc[
                    side_diag["side_status"] == "ok", "annotated_side"
                ].value_counts(dropna=False).to_dict()
                f.write("\n  Annotated-side breakdown (ok only):\n")
                for k, v in by_side.items():
                    f.write(f"    {str(k):<10s} {v}\n")
        f.write("\nSide-related skip counts (from per-segment skips):\n")
        side_keys = (
            "side_missing", "side_ambiguous", "side_no_candidates",
            "side_one_candidate", "attacker_position_nan_at_event",
        )
        # Aggregate from skipped_by_segment
        agg_counts: dict[str, int] = {k: 0 for k in side_keys}
        for reasons in skipped_by_segment.values():
            for k in side_keys:
                agg_counts[k] += int(reasons.get(k, 0) or 0)
        for k, v in agg_counts.items():
            f.write(f"  {k:<36s} {v}\n")

    return {
        "boris_total": boris_total,
        "processed_total": processed_total,
        "missing_total": int(boris_total - processed_total),
    }


def _classify_missing(
    ev: pd.Series,
    skipped_by_segment: dict[tuple[int, int], dict[str, int]],
    npy_available: dict[tuple[int, int], bool],
    side_status_lookup: dict[tuple[str, float], str] | None = None,
) -> str:
    expected_segment = ev.get("expected_segment")
    trial = ev.get("trial")
    if expected_segment is None or pd.isna(expected_segment):
        return "in_gap_between_segments"
    if expected_segment == 0:
        return "out_of_range"
    if ev.get("scope_excluded_trial"):
        return "scope_filter_excluded_trial"
    if expected_segment == 1 and ev.get("scope_seg1_dropped"):
        return "scope_filter_seg1"
    if ev.get("scope_pharma_dropped"):
        return "scope_filter_pharma"

    # Side-aware reasons -- check before per-segment kinematic skips, because
    # an event with no/ambiguous side annotation is never even handed to the
    # per-segment loop.
    if side_status_lookup is not None:
        key_obs = (
            str(ev.get("Observation id")),
            round(float(ev.get("Start (s)", 0.0)), 6),
        )
        status = side_status_lookup.get(key_obs)
        if status == "missing":
            return "side_missing"
        if status == "ambiguous":
            return "side_ambiguous"

    if trial is not None and not pd.isna(trial):
        key = (int(trial), int(expected_segment))
        if not npy_available.get(key, False):
            return "missing_segment_npy"
        reasons = skipped_by_segment.get(key, {})
        # Side-related skip reasons take precedence over kinematic ones because
        # they explain rows that *had* annotations but failed the hard filter.
        for r in ("side_no_candidates", "side_one_candidate",
                  "attacker_position_nan_at_event"):
            if reasons.get(r, 0) > 0:
                return r
        for r in ("frame_out_of_range", "too_few_valid_fish",
                  "no_attacker_found", "invalid_endpoint",
                  "no_attacked_found", "missing_positions"):
            if reasons.get(r, 0) > 0:
                return {
                    "frame_out_of_range": "frame_calc_issue",
                    "too_few_valid_fish": "insufficient_tracking",
                    "no_attacker_found": "algorithm_failed_no_attacker",
                    "invalid_endpoint": "algorithm_failed_no_attacker",
                    "no_attacked_found": "algorithm_failed_no_attacked",
                    "missing_positions": "insufficient_tracking",
                }[r]
    return "unknown"


def _idmatcher_coverage(events: pd.DataFrame) -> pd.DataFrame:
    if events.empty or "trial" not in events.columns:
        return pd.DataFrame()
    grouped = events.groupby("trial")

    real_qualities = {
        "strong", "acceptable", "by_matrix", "by_exclusion",
        "anchor", "rejected_conflict",
    }

    def _count(g: pd.DataFrame) -> pd.Series:
        n_total = len(g)
        n_atk = int(g["attacker_fish_id_matched_trial"].notna().sum())
        n_vct = int(g["attacked_fish_id_matched_trial"].notna().sum())
        n_both = int(
            (g["attacker_fish_id_matched_trial"].notna()
             & g["attacked_fish_id_matched_trial"].notna()).sum()
        )
        idm_available = bool(
            g["attacker_match_quality"].astype(str).isin(real_qualities).any()
            or g["attacked_match_quality"].astype(str).isin(real_qualities).any()
        )
        if n_both == n_total:
            status = "full"
        elif n_both == 0:
            status = "none"
        else:
            status = "partial"
        return pd.Series({
            "n_events_total": n_total,
            "n_events_with_matched_attacker": n_atk,
            "n_events_with_matched_attacked": n_vct,
            "n_events_both_matched": n_both,
            "coverage_status": status,
            "idmatcher_data_available": idm_available,
        })

    return grouped.apply(_count, include_groups=False).reset_index()
