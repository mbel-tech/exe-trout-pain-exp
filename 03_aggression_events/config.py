"""Paths and thresholds for the BORIS-to-trajectory aggression workflow
(Supplementary Material, "Quantification of aggression events from BORIS
annotations and idtracker.ai trajectories").

This is the configuration of the run used in the article (June 2026). All input
paths are relative to the repository root and point into ``data/``; outputs are
written to ``03_aggression_events/results/``. Run from this folder with
``python aggression_pipeline.py all``.
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"

# ---------------------------------------------------------------------------
# Input roots
# ---------------------------------------------------------------------------
# Root folder containing one subfolder per session:
#   session_trial{N}_Segment{M}/trajectories/trajectories.{h5,npy}
# (70 session folders, each carrying a trajectories.npy + trajectories.h5.)
SESSIONS_DIR = DATA_DIR / "tracking_sessions"

# BORIS per-event export (Aggregated events workbook). This is a full BORIS export (one row
# per scored event) carrying paired CALM SIDE / FLOW SIDE rows that share the
# same `Start (s)` value as each `aggression` row.
BORIS_FILE = DATA_DIR / "aggression" / "new_aggression_15.6.2026.xlsx"
BORIS_SHEET_NAME = "Aggregated events"
# BORIS's native column name for the event time.
BORIS_START_COLUMN = "Start (s)"

# idmatcher.ai annotated table written by 02_trajectory_pipeline/01_5_annotate_identity_matches.R
# (12 trials x 4 fish = 48 rows). Set to None to disable the cross-segment
# identity join entirely; the pipeline will then only produce segment-local IDs
# and flag this loudly in the log.
# Patched version: trial 8 Seg 3 rows corrected to matrix-optimal
# assignment (Seg2_1->Seg3_1, Seg2_2->Seg3_4, Seg2_3->Seg3_2, Seg2_4->Seg3_3),
# quality=by_matrix, match_accepted=True. Base: STEP1_5_20260524_160045.
IDMATCHER_ANNOTATED_CSV = DATA_DIR / "aggression" / "idmatcher_annotated_june2026_patched.csv"

# Optional body-length diagnostics CSV produced by the R helper. When provided,
# per-event rows gain attacker/attacked_body_length_mean_cm columns. When None,
# those columns are emitted as NA.
# Patched version: trial 8 rows added from STEP2 bodylength_by_fish.csv
# using the matrix-corrected identity mapping. Base: STEP1_5_20260524_160045.
BODY_LENGTH_DIAGNOSTICS_CSV: Path | None = DATA_DIR / "aggression" / "body_length_diagnostics_june2026_patched.csv"

# Trials metadata table (one row per trial × timepoint).
# Columns: Trial, Group, Treatment, Date, timepoint.
# Joined onto the master events table at merge time (Trial=trial,
# timepoint=segment) to add clean `group`, `treatment_clean`, `date` columns.
# Set to None to skip the join.
TRIALS_SUMMARY_XLSX = DATA_DIR / "trials summary GOT.xlsx"


# ---------------------------------------------------------------------------
# Output roots
# ---------------------------------------------------------------------------
# All outputs are written under 03_aggression_events/results/.
WORK_DIR = Path(__file__).parent / "results"
# Staged + interpolated NPYs, staged from the tracking sessions in SESSIONS_DIR.
NPY_RAW_DIR = WORK_DIR / "npy_raw"
NPY_INTERP_DIR = WORK_DIR / "npy_interpolated"
# Per-run output subfolders.
RESULTS_DIR = WORK_DIR / "boris_event_results_v2_side"
QC_DIR = WORK_DIR / "qc_missing_events_side"
LOG_DIR = WORK_DIR / "logs_side"

# ---------------------------------------------------------------------------
# Clean CSV export for the R statistics
# ---------------------------------------------------------------------------
# Single flat CSV, one row per scored aggression event. Set to None to skip.
CLEAN_CSV_PATH = WORK_DIR / "aggression_events_clean.csv"


# ---------------------------------------------------------------------------
# Kinematic defaults & detection thresholds (block H of the original script)
# ---------------------------------------------------------------------------
FPS_DEFAULT = 25.0                # used only when NPY metadata lacks fps
BODY_LEN_DEFAULT_PX = 100.0       # used only when NPY metadata lacks body_length

# Segment timing (seconds, relative to observation start = time of the 'START'
# BORIS event for each trial). Three 10-min windows separated by 10-min gaps.
SEGMENT_SCHEDULE = [(0.0, 600.0), (1200.0, 1800.0), (2400.0, 3000.0)]
SEGMENT_NAMES = ["Segment1", "Segment2", "Segment3"]

# Aggressor detection (v2)
WINDOW_S = 1.0                    # +/- window around BORIS event time
ACCEL_FRAMES = 5                  # n-frame stride for finite-difference accel
ENDPOINT_S = 0.5                  # forward projection horizon for attacker
MIN_VALID_FISH = 2                # need at least this many fish with positions
MIN_VALID_FRAMES_IN_WINDOW = 3    # per fish, frames with non-NaN position

# Interpolation (block C STEP 3)
INTERP_MAX_GAP_FRAMES = 24        # ~1 s @ 25 fps
JUMP_THRESHOLD_PX = 80.0          # rolling-median outlier threshold
JUMP_WINDOW = 11                  # rolling-median window length

# Interaction-type thresholds in px/frame^2 (kept identical to v2 for traceability)
INTERACTION_THRESHOLDS = {
    "high_intensity": 10.0,
    "medium_intensity": 5.0,
    "low_intensity": 0.0,
}


# ---------------------------------------------------------------------------
# Zone-assignment (calm / flow) -- ported from R main pipeline
# `01_idtracker_pipeline_GOT.R::zone_occupancy_y_boundary` (lines 1408-1443).
# Boundary is a single horizontal pixel-y threshold; trial- and segment-constant.
# ---------------------------------------------------------------------------
ZONE_Y0_PX = 480.0                # y-coordinate split between calm/flow zones
ZONE_FLOW_IS_HIGHER_Y = True      # True -> y >= y0 is "flow"; False -> reversed


# ---------------------------------------------------------------------------
# Scope-filter switches (plan section 11) -- unchanged from config.py
# ---------------------------------------------------------------------------
# Trial numbers to exclude entirely from the pipeline. Excludes any NPY files
# and any BORIS events from these trial IDs.
EXCLUDED_TRIALS = [1, 2, 3, 4]

# Drop Segment 1 from extraction and from BORIS event matching.
TP1_SWITCH = True

# Drop trials whose treatment/condition contains "naloxone" (or the typo
# "nloxone") or "morphine" (case-insensitive).
DROP_PHARMA_GROUPS = True

# Read BORIS at extract time and skip H5/NPY staging for pharma trials.
# Off by default to keep STEP 1 independent of STEP 4 inputs.
DROP_PHARMA_AT_EXTRACT = False


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
LOG_LEVEL = "INFO"
