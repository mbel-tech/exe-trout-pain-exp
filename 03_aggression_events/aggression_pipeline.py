"""CLI entry point for the cleaned BORIS-to-tracking aggression pipeline.

Subcommands:

    extract       Walk SESSIONS_DIR -> stage canonical NPYs in NPY_RAW_DIR.
                  Applies identity-label cleanup; honours TP1_SWITCH.

    interpolate   Apply jump-protection + PCHIP gap fill -> NPY_INTERP_DIR.

    assign        Run the v2 aggressor detector on each segment, attach the
                  idmatcher cross-segment ID, write one CSV per session to
                  RESULTS_DIR.

    merge         Concatenate all per-segment CSVs into a master XLSX with
                  descriptive summaries (no inferential statistics).

    qc            Emit the comprehensive missing-events + idmatcher-coverage
                  reports to QC_DIR.

    all           Run extract -> interpolate -> assign -> merge -> qc in order.

Each subcommand reads ``config.py`` for paths and thresholds; override at the
command line with ``--config-key VALUE`` (limited; see ``--help``).
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import config
from trajectory_io import (
    discover_session_files,
    stage_session,
    load_npy_payload,
    save_npy_payload,
    canonical_npy_name,
)
from interpolation import interpolate_trajectory, summarize_interpolation
from idmatcher_join import (
    load_idmatcher_annotated,
    build_idmatcher_map,
    covered_trials as idmatcher_covered_trials,
    load_body_length_diagnostics,
)
from boris_matching import (
    load_boris,
    get_observation_start,
    extract_trial_info,
    assign_events_for_segment,
    attach_idmatcher_columns,
    is_pharma,
    in_segment1,
    build_side_annotations,
)
from merge_results import (
    merge_event_csvs,
    write_master_xlsx,
    derive_master_filename,
)
from qc_missing_events import write_missing_events_reports


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

def setup_logging(log_dir: Path, name: str) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = log_dir / f"{name}_{ts}.log"
    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    logging.basicConfig(
        level=getattr(logging, config.LOG_LEVEL.upper(), logging.INFO),
        format=fmt,
        handlers=[
            logging.FileHandler(path, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )
    logging.getLogger(__name__).info("Logging to %s", path)
    return path


# ---------------------------------------------------------------------------
# Scope banner
# ---------------------------------------------------------------------------

def scope_banner() -> str:
    bits = [
        f"side_aware=True+hybrid_fallback (y0={config.ZONE_Y0_PX}, "
        f"flow_higher_y={config.ZONE_FLOW_IS_HIGHER_Y})"
    ]
    if config.EXCLUDED_TRIALS:
        bits.append(f"EXCLUDED_TRIALS={config.EXCLUDED_TRIALS}")
    if config.TP1_SWITCH:
        bits.append("TP1_SWITCH=True (seg1 dropped)")
    if config.DROP_PHARMA_GROUPS:
        bits.append("DROP_PHARMA_GROUPS=True (naloxone/morphine dropped)")
    return "; ".join(bits)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_extract(args: argparse.Namespace) -> int:
    log = logging.getLogger("extract")
    config.NPY_RAW_DIR.mkdir(parents=True, exist_ok=True)
    entries = discover_session_files(config.SESSIONS_DIR, excluded_trials=config.EXCLUDED_TRIALS)
    log.info("Discovered %d sessions under %s", len(entries), config.SESSIONS_DIR)
    if config.EXCLUDED_TRIALS:
        log.info("EXCLUDED_TRIALS=%s (applied at discovery)", config.EXCLUDED_TRIALS)

    if config.TP1_SWITCH:
        before = len(entries)
        entries = [e for e in entries if e["segment"] != 1]
        log.info("TP1_SWITCH=True -> dropping %d Segment1 sessions", before - len(entries))

    cleanup_log: list[dict] = []
    for entry in entries:
        try:
            out, info = stage_session(
                entry,
                out_dir=config.NPY_RAW_DIR,
                fps_default=config.FPS_DEFAULT,
                body_len_default_px=config.BODY_LEN_DEFAULT_PX,
                overwrite=bool(args.overwrite),
            )
            cleanup_log.append({
                "trial": entry["trial"],
                "segment": entry["segment"],
                "out_path": str(out),
                "original_labels": info["original_labels"],
                "dropped": info["dropped"],
                "reason": info["reason"],
            })
        except Exception as exc:
            log.error("Failed staging %s: %s", entry["session_dir"], exc)
            cleanup_log.append({
                "trial": entry["trial"],
                "segment": entry["segment"],
                "out_path": "",
                "original_labels": "",
                "dropped": "",
                "reason": f"ERROR: {exc}",
            })

    pd.DataFrame(cleanup_log).to_csv(
        config.NPY_RAW_DIR / "extraction_log.csv", index=False
    )
    log.info("Wrote %d entries to extraction_log.csv", len(cleanup_log))
    return 0


def cmd_interpolate(args: argparse.Namespace) -> int:
    log = logging.getLogger("interpolate")
    config.NPY_INTERP_DIR.mkdir(parents=True, exist_ok=True)
    npys = sorted(config.NPY_RAW_DIR.glob("trajectories_trial*_Segment*.npy"))
    if not npys:
        log.error("No staged NPYs found under %s. Run extract first.",
                  config.NPY_RAW_DIR)
        return 1

    summary_log: list[dict] = []
    for src in npys:
        out_path = config.NPY_INTERP_DIR / src.name
        if out_path.exists() and not args.overwrite:
            log.info("[%s] already interpolated; skipping", src.name)
            continue
        payload = load_npy_payload(src)
        before = payload["trajectories"]
        after = interpolate_trajectory(
            before,
            max_gap=config.INTERP_MAX_GAP_FRAMES,
            remove_jumps=True,
            jump_threshold=config.JUMP_THRESHOLD_PX,
            jump_window=config.JUMP_WINDOW,
        )
        payload["trajectories"] = after
        payload["interpolation"] = summarize_interpolation(before, after)
        save_npy_payload(out_path, payload)
        log.info("[%s] %s", src.name, payload["interpolation"])
        summary_log.append({
            "file": src.name,
            **payload["interpolation"],
        })

    if summary_log:
        pd.DataFrame(summary_log).to_csv(
            config.NPY_INTERP_DIR / "interpolation_log.csv", index=False
        )
    return 0


def cmd_assign(args: argparse.Namespace) -> int:
    log = logging.getLogger("assign")
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    boris_full = load_boris(
        config.BORIS_FILE,
        sheet_name=getattr(config, "BORIS_SHEET_NAME", None),
        start_column=getattr(config, "BORIS_START_COLUMN", "Start (s)"),
    )
    log.info("Loaded BORIS file %s: %d rows", config.BORIS_FILE, len(boris_full))

    # Build side-annotation lookup (side-aware variant)
    side_lookup, side_diag = build_side_annotations(boris_full)
    n_ok = int((side_diag["side_status"] == "ok").sum())
    n_missing = int((side_diag["side_status"] == "missing").sum())
    n_ambig = int((side_diag["side_status"] == "ambiguous").sum())
    log.info(
        "Side annotations: %d ok, %d missing, %d ambiguous "
        "(zone y0=%s px, flow_higher_y=%s)",
        n_ok, n_missing, n_ambig,
        config.ZONE_Y0_PX, config.ZONE_FLOW_IS_HIGHER_Y,
    )

    # Build idmatcher map
    if config.IDMATCHER_ANNOTATED_CSV is None:
        log.warning("IDMATCHER_ANNOTATED_CSV is None; "
                    "fish_id_matched_trial will be NA for all events.")
        idmatcher_map = None
        covered_set: set[int] = set()
    else:
        annotated = load_idmatcher_annotated(config.IDMATCHER_ANNOTATED_CSV)
        idmatcher_map = build_idmatcher_map(annotated)
        covered_set = set(idmatcher_covered_trials(annotated))
        log.info("Built idmatcher map: %d rows; %d covered trials",
                 len(idmatcher_map), len(covered_set))

    bl_diag = load_body_length_diagnostics(config.BODY_LENGTH_DIAGNOSTICS_CSV)

    # Aggression-only subset (drop START rows from event-loop scope)
    agg = boris_full[boris_full["Behavior"].str.lower() == "aggression"].copy()
    log.info("Aggression events: %d (out of %d BORIS rows)", len(agg), len(boris_full))

    # Parse observation IDs once
    parsed = agg["Observation id"].map(extract_trial_info)
    agg["trial"] = parsed.map(lambda d: d.get("trial"))
    agg["treatment"] = parsed.map(lambda d: d.get("treatment"))
    agg["condition"] = parsed.map(lambda d: d.get("condition"))

    # Pre-fetch observation start times
    obs_start_lookup: dict[str, float] = {}
    for obs_id in agg["Observation id"].unique():
        st = get_observation_start(boris_full, obs_id)
        if st is None:
            log.warning("[%s] no START marker found; using 0.0", obs_id)
            st = 0.0
        obs_start_lookup[obs_id] = st

    # Determine expected segment per event (1, 2, 3, 0=out-of-range, NaN=in gap)
    def _expected_segment(row: pd.Series) -> int | float:
        obs_start = obs_start_lookup.get(row["Observation id"], 0.0)
        rel = float(row["Start (s)"]) - obs_start
        if rel < 0:
            return 0
        for i, (start, end) in enumerate(config.SEGMENT_SCHEDULE, start=1):
            if start <= rel < end:
                return i
        # past last segment
        if rel >= config.SEGMENT_SCHEDULE[-1][1]:
            return 0
        return float("nan")
    agg["expected_segment"] = agg.apply(_expected_segment, axis=1)

    # Scope filters
    agg["scope_seg1_dropped"] = False
    agg["scope_pharma_dropped"] = False
    agg["scope_excluded_trial"] = False

    if config.EXCLUDED_TRIALS:
        mask = agg["trial"].isin(config.EXCLUDED_TRIALS)
        agg.loc[mask, "scope_excluded_trial"] = True
        log.info("EXCLUDED_TRIALS=%s dropped %d events", config.EXCLUDED_TRIALS, int(mask.sum()))

    if config.DROP_PHARMA_GROUPS:
        mask = agg.apply(
            lambda r: is_pharma(r.get("treatment")) or is_pharma(r.get("condition")),
            axis=1,
        )
        agg.loc[mask, "scope_pharma_dropped"] = True
        log.info("DROP_PHARMA_GROUPS dropped %d events", int(mask.sum()))
    if config.TP1_SWITCH:
        mask = agg["expected_segment"] == 1
        agg.loc[mask, "scope_seg1_dropped"] = True
        log.info("TP1_SWITCH dropped %d events", int(mask.sum()))

    agg_scope = agg[~(agg["scope_seg1_dropped"] | agg["scope_pharma_dropped"] | agg["scope_excluded_trial"])].copy()

    # Iterate over staged segment files
    npys = sorted(config.NPY_INTERP_DIR.glob("trajectories_trial*_Segment*.npy"))
    if not npys:
        log.error("No interpolated NPYs found under %s. Run interpolate first.",
                  config.NPY_INTERP_DIR)
        return 1

    skipped_by_segment: dict[tuple[int, int], dict[str, int]] = {}
    npy_available: dict[tuple[int, int], bool] = {}

    # Map (trial, segment) -> npy
    seg_to_npy: dict[tuple[int, int], Path] = {}
    for p in npys:
        from trajectory_io import parse_session_name
        parsed_name = parse_session_name(p.stem)
        if parsed_name is None:
            log.warning("Cannot parse trial/segment from %s", p.name)
            continue
        seg_to_npy[parsed_name] = p
        npy_available[parsed_name] = True

    log.info("Processing %d (trial, segment) pairs", len(seg_to_npy))

    n_events_total = 0
    for (trial, segment), npy_path in sorted(seg_to_npy.items()):
        out_csv = config.RESULTS_DIR / f"trajectories_trial{trial}_Segment{segment}_events.csv"
        if config.TP1_SWITCH and segment == 1:
            # Overwrite any stale Segment-1 CSV so the merge step sees an
            # empty file rather than old pre-TP1_SWITCH results.
            pd.DataFrame().to_csv(out_csv, index=False)
            continue
        # Aggression events from this trial
        sub = agg_scope[agg_scope["trial"] == trial]
        if sub.empty:
            # Write/overwrite with an empty CSV so stale files from earlier
            # runs (without pharma filter, or without the `side` column) don't
            # pollute the merge step.
            pd.DataFrame().to_csv(out_csv, index=False)
            continue
        # Only feed events whose absolute time falls in this segment's window
        obs_first = sub.iloc[0]["Observation id"]
        obs_start = obs_start_lookup.get(obs_first, 0.0)
        treatment = sub.iloc[0].get("treatment")
        condition = sub.iloc[0].get("condition")

        events_df, skipped = assign_events_for_segment(
            npy_path=npy_path,
            boris_events=sub,
            obs_start_abs=obs_start,
            trial=trial,
            segment_idx=segment,
            segment_schedule=config.SEGMENT_SCHEDULE,
            fps_default=config.FPS_DEFAULT,
            body_len_default_px=config.BODY_LEN_DEFAULT_PX,
            window_s=config.WINDOW_S,
            accel_frames=config.ACCEL_FRAMES,
            endpoint_s=config.ENDPOINT_S,
            min_valid_fish=config.MIN_VALID_FISH,
            min_valid_frames_in_window=config.MIN_VALID_FRAMES_IN_WINDOW,
            interaction_thresholds=config.INTERACTION_THRESHOLDS,
            treatment=treatment,
            condition=condition,
            side_lookup=side_lookup,
            zone_y0_px=config.ZONE_Y0_PX,
            zone_flow_higher_y=config.ZONE_FLOW_IS_HIGHER_Y,
        )
        skipped_by_segment[(trial, segment)] = skipped
        if events_df.empty:
            log.info("[T%02d Seg%d] 0 events (skipped %s)", trial, segment, skipped)
            events_df.to_csv(out_csv, index=False)  # overwrite any stale CSV
            continue

        # Attach idmatcher columns
        events_df = attach_idmatcher_columns(
            events_df,
            idmatcher_map=idmatcher_map,
            covered_trial_set=covered_set,
            bl_diag=bl_diag,
        )

        events_df.to_csv(out_csv, index=False)
        n_events_total += len(events_df)
        log.info("[T%02d Seg%d] %d events written (skipped %s)",
                 trial, segment, len(events_df), skipped)

    # Persist auxiliary tables for the QC step
    (config.RESULTS_DIR / "_aux").mkdir(exist_ok=True)
    pd.DataFrame([
        {"trial": t, "segment": s, **reasons}
        for (t, s), reasons in skipped_by_segment.items()
    ]).to_csv(config.RESULTS_DIR / "_aux" / "skipped_by_segment.csv", index=False)
    pd.DataFrame([
        {"trial": t, "segment": s, "npy_available": a}
        for (t, s), a in npy_available.items()
    ]).to_csv(config.RESULTS_DIR / "_aux" / "npy_available.csv", index=False)
    agg.to_csv(config.RESULTS_DIR / "_aux" / "boris_aggression_parsed.csv", index=False)
    # Persist the side-annotation diagnostic so QC can report missing/ambiguous
    side_diag.to_csv(
        config.RESULTS_DIR / "_aux" / "side_annotations.csv", index=False,
    )

    log.info("Assign complete: %d total events across %d segment files",
             n_events_total, len(seg_to_npy))
    return 0


def _attach_trials_metadata(df: pd.DataFrame, log: logging.Logger) -> pd.DataFrame:
    """Left-join the trials-summary metadata onto the events dataframe.

    Adds ``group``, ``treatment_clean``, and ``date`` columns (from
    ``config.TRIALS_SUMMARY_XLSX``) keyed on ``trial`` × ``segment``
    (= Trial × timepoint in the summary file).  Columns are inserted
    immediately after ``segment`` for readability.  Skipped when
    ``config.TRIALS_SUMMARY_XLSX`` is None or the file doesn't exist.
    """
    src = getattr(config, "TRIALS_SUMMARY_XLSX", None)
    if src is None or not Path(src).exists():
        log.warning("TRIALS_SUMMARY_XLSX not set or missing; skipping metadata join.")
        return df

    try:
        meta = pd.read_excel(src)
    except Exception as exc:
        log.warning("Could not read TRIALS_SUMMARY_XLSX (%s): %s; skipping join.", src, exc)
        return df

    meta = meta.rename(columns={
        "Trial": "trial",
        "timepoint": "segment",
        "Group": "group",
        "Treatment": "treatment_clean",
        "Date": "date",
    })[["trial", "segment", "group", "treatment_clean", "date"]]

    before = len(df)
    merged = df.merge(meta, on=["trial", "segment"], how="left")
    n_unmatched = merged["group"].isna().sum()
    if n_unmatched:
        log.warning(
            "Trials-metadata join: %d / %d events have no match in summary "
            "(trial/segment pairs absent from TRIALS_SUMMARY_XLSX).",
            n_unmatched, before,
        )
    else:
        log.info("Trials-metadata join: all %d events matched (group/treatment_clean/date added).", before)

    # Reorder: insert group, treatment_clean, date right after 'segment'
    new_cols = ["group", "treatment_clean", "date"]
    other_cols = [c for c in merged.columns if c not in new_cols]
    insert_pos = other_cols.index("segment") + 1
    ordered = other_cols[:insert_pos] + new_cols + other_cols[insert_pos:]
    return merged[ordered]


def _write_skinny_table(df: pd.DataFrame, out_dir: Path, log: logging.Logger) -> None:
    """Emit a slim modelling-ready table.

    Keeps only columns relevant for interpretation and statistical modelling.
    Renames a handful of columns for brevity. Writes both CSV and XLSX next
    to the master output so analysts can grab whichever they prefer.
    """
    skinny_cols = [
        ("trial",                              "trial"),
        ("segment",                            "segment"),
        ("group",                              "group"),
        ("treatment_clean",                    "treatment"),
        ("date",                               "date"),
        ("time_within_segment_sec",            "time_sec"),
        ("attacker_fish_id_matched_trial",     "attacker_id"),
        ("attacked_fish_id_matched_trial",     "attacked_id"),
        ("side",                               "side"),
        ("individuals_in_calm",                "n_calm"),
        ("interaction_type",                   "intensity_class"),
        ("attacker_max_acceleration_BL_per_s2", "acceleration_BL_s2"),
        ("distance_at_endpoint_BL",            "distance_BL"),
        ("attacker_body_length_mean_cm",       "attacker_body_length_cm"),
        ("attacked_body_length_mean_cm",       "attacked_body_length_cm"),
    ]
    missing = [src for src, _ in skinny_cols if src not in df.columns]
    if missing:
        log.warning("Skinny table: missing source columns %r; skipping.", missing)
        return
    skinny = df[[src for src, _ in skinny_cols]].copy()
    skinny.columns = [dst for _, dst in skinny_cols]

    csv_path = out_dir / "master_events_skinny.csv"
    xlsx_path = out_dir / "master_events_skinny.xlsx"
    skinny.to_csv(csv_path, index=False)
    skinny.to_excel(xlsx_path, index=False, sheet_name="Events")
    log.info("Skinny table written: %s (rows=%d, cols=%d)",
             csv_path, len(skinny), len(skinny.columns))


def _coerce_bool(series: "pd.Series") -> "pd.Series":
    """Robustly coerce a possibly-string/object column to a plain bool Series.

    CSV round-trips render Python booleans as ``"True"``/``"False"`` strings and
    unmatched rows as NaN; treat anything that is not an explicit truthy token as
    False.
    """
    if series.dtype == bool:
        return series.fillna(False)
    return (
        series.astype(str).str.strip().str.lower()
        .isin(["true", "1", "t", "yes"])
    )


def _apply_side_override(df: "pd.DataFrame", log: logging.Logger) -> "pd.DataFrame":
    """Trajectory-primary / annotation-fallback resolution of the event ``side``.

    ``side`` becomes the trajectory-derived attacker side
    (``attacker_side_at_event``) whenever the attacker was clearly matched
    (``attacker_match_accepted`` is True). When the attacker was NOT clearly
    matched, fall back to the human ``annotated_side``; if that is also missing,
    keep the trajectory side so the event never loses its side label.
    """
    required = {
        "attacker_side_at_event", "annotated_side",
        "attacker_match_accepted", "side",
    }
    missing = required - set(df.columns)
    if missing or df.empty:
        log.warning("Side override skipped (missing columns: %r).", sorted(missing))
        return df

    accepted = _coerce_bool(df["attacker_match_accepted"])
    traj = df["attacker_side_at_event"]
    annot = df["annotated_side"]
    annot_present = annot.notna() & (annot.astype(str).str.strip() != "")

    # not-accepted -> annotation if present, else trajectory; accepted -> trajectory
    fallback = annot.where(annot_present, traj)
    new_side = traj.where(accepted, fallback)

    sentinel = "∅"  # ∅, so NaN==NaN compares equal for the change count
    n_changed = int((new_side.fillna(sentinel) != df["side"].fillna(sentinel)).sum())
    df = df.copy()
    df["side"] = new_side
    log.info(
        "Side override applied (trajectory primary, annotation fallback on "
        "unclear attacker): %d/%d event sides changed; %d rows fell back to "
        "annotation.",
        n_changed, len(df), int((~accepted).sum()),
    )
    return df


def cmd_merge(args: argparse.Namespace) -> int:
    log = logging.getLogger("merge")
    df = merge_event_csvs(config.RESULTS_DIR)
    df = _attach_trials_metadata(df, log)

    # ── Side override (June-2026) ──────────────────────────────────────────────
    # Make the trajectory-derived ATTACKER side the authoritative event `side`,
    # and only fall back to the human annotation when the attacker could NOT be
    # clearly matched from the trajectory (attacker_match_accepted == False). If
    # the annotation is also absent for such a row, keep the trajectory side so no
    # event loses its side. On the current data this reproduces the existing
    # `side` exactly (the trajectory side already always wins and never disagrees
    # with the annotation where present); the rule is encoded explicitly so it
    # remains correct for future data where the two could diverge.
    df = _apply_side_override(df, log)

    base = "master_events_v2"
    name = derive_master_filename(
        base,
        tp1_switch=config.TP1_SWITCH,
        drop_pharma=config.DROP_PHARMA_GROUPS,
        side_aware=True,
        hybrid=True,
    )
    out = config.RESULTS_DIR / name
    # Hand the side-annotation diagnostic + skip counts to merge so it can
    # write the Side_Coverage / Annotated_vs_Trajectory_Side sheets.
    side_diag_path = config.RESULTS_DIR / "_aux" / "side_annotations.csv"
    side_diag = pd.read_csv(side_diag_path) if side_diag_path.exists() else None
    skipped_path = config.RESULTS_DIR / "_aux" / "skipped_by_segment.csv"
    skipped_df = pd.read_csv(skipped_path) if skipped_path.exists() else None
    info = write_master_xlsx(
        df, out,
        scope_banner=scope_banner(),
        side_diag=side_diag,
        skipped_by_segment=skipped_df,
    )
    df.to_csv(out.with_suffix(".csv"), index=False)
    log.info("Master XLSX written: %s (rows=%d)", out, len(df))
    log.info("Summary info: %s", info)

    # R8: also write a canonical clean CSV at a fixed location (locked schema
    # for the future R consumer). Skip if CLEAN_CSV_PATH is None.
    clean_csv_path = getattr(config, "CLEAN_CSV_PATH", None)
    if clean_csv_path is not None:
        try:
            clean_csv_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(clean_csv_path, index=False)
            log.info(
                "R8: clean CSV written for R consumer: %s (rows=%d, cols=%d)",
                clean_csv_path, len(df), len(df.columns),
            )
        except Exception as exc:
            log.warning("R8: failed to write clean CSV at %s: %s",
                        clean_csv_path, exc)

    # Emit the slim modelling-ready table alongside the master.
    _write_skinny_table(df, config.RESULTS_DIR, log)
    return 0


def cmd_qc(args: argparse.Namespace) -> int:
    log = logging.getLogger("qc")
    events_csv = config.RESULTS_DIR / derive_master_filename(
        "master_events_v2",
        config.TP1_SWITCH,
        config.DROP_PHARMA_GROUPS,
        side_aware=True,
        hybrid=True,
    ).replace(".xlsx", ".csv")
    if not events_csv.exists():
        log.error("Master CSV not found at %s. Run `merge` first.", events_csv)
        return 1
    events = pd.read_csv(events_csv)

    aux = config.RESULTS_DIR / "_aux"
    skipped = {}
    if (aux / "skipped_by_segment.csv").exists():
        sk = pd.read_csv(aux / "skipped_by_segment.csv")
        for _, r in sk.iterrows():
            key = (int(r["trial"]), int(r["segment"]))
            skipped[key] = {
                c: int(r[c]) for c in r.index
                if c not in ("trial", "segment") and pd.notna(r[c])
            }
    npy_available: dict[tuple[int, int], bool] = {}
    if (aux / "npy_available.csv").exists():
        nv = pd.read_csv(aux / "npy_available.csv")
        for _, r in nv.iterrows():
            npy_available[(int(r["trial"]), int(r["segment"]))] = bool(r["npy_available"])

    boris_parsed = pd.read_csv(aux / "boris_aggression_parsed.csv")

    side_diag_path = aux / "side_annotations.csv"
    side_diag = pd.read_csv(side_diag_path) if side_diag_path.exists() else None

    scope_skipped = {
        "skipped_excluded_trial": int(boris_parsed["scope_excluded_trial"].sum()) if "scope_excluded_trial" in boris_parsed.columns else 0,
        "skipped_tp1_switch": int(boris_parsed["scope_seg1_dropped"].sum()) if "scope_seg1_dropped" in boris_parsed.columns else 0,
        "skipped_pharma_group": int(boris_parsed["scope_pharma_dropped"].sum()) if "scope_pharma_dropped" in boris_parsed.columns else 0,
    }
    if side_diag is not None:
        scope_skipped["side_missing"] = int((side_diag["side_status"] == "missing").sum())
        scope_skipped["side_ambiguous"] = int((side_diag["side_status"] == "ambiguous").sum())

    write_missing_events_reports(
        boris=boris_parsed,
        events=events,
        skipped_by_segment=skipped,
        npy_available=npy_available,
        qc_dir=config.QC_DIR,
        scope_banner=scope_banner(),
        scope_skipped=scope_skipped,
        side_diag=side_diag,
        zone_y0_px=config.ZONE_Y0_PX,
        zone_flow_higher_y=config.ZONE_FLOW_IS_HIGHER_Y,
    )
    log.info("QC reports written to %s", config.QC_DIR)
    return 0


def cmd_all(args: argparse.Namespace) -> int:
    for fn in (cmd_extract, cmd_interpolate, cmd_assign, cmd_merge, cmd_qc):
        code = fn(args)
        if code:
            return code
    return 0


# ---------------------------------------------------------------------------
# argparse plumbing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="BORIS-to-tracking aggression pipeline (v2)."
    )
    sub = p.add_subparsers(dest="command", required=True)

    for name, fn in [
        ("extract", cmd_extract),
        ("interpolate", cmd_interpolate),
        ("assign", cmd_assign),
        ("merge", cmd_merge),
        ("qc", cmd_qc),
        ("all", cmd_all),
    ]:
        sp = sub.add_parser(name, help=fn.__doc__)
        sp.add_argument(
            "--overwrite", action="store_true",
            help="Overwrite existing staged files / outputs",
        )
        sp.set_defaults(func=fn)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(config.LOG_DIR, name=args.command)
    log = logging.getLogger("main")
    log.info("Scope: %s", scope_banner())
    log.info("WORK_DIR=%s", config.WORK_DIR)
    try:
        return args.func(args)
    except Exception:
        log.error("Pipeline failed:\n%s", traceback.format_exc())
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
