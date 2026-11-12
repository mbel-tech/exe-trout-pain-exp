"""BORIS-event-to-trajectory matching with v2 aggressor assignment.

Implements block H of the original Colab script (the canonical aggressor
assignment) plus the idmatcher attachment from plan section 10 and the
scope-filter switches from plan section 11.

Per-event flow
--------------
For each BORIS aggression event:

1. Locate the segment of the observation start (the ``START`` BORIS row
   defines absolute ``t=0`` for that trial).
2. Map absolute event time to a segment index via ``SEGMENT_SCHEDULE``.
3. Convert to a frame index inside that segment's NPY at the segment's fps.
4. In a +/- ``WINDOW_S`` window, compute per-fish velocity and ``ACCEL_FRAMES``-stride
   finite-difference acceleration. Skip any fish with <
   ``MIN_VALID_FRAMES_IN_WINDOW`` valid (non-NaN) positions.
5. Select the fish with the highest acceleration in the window (regardless of
   sign) as the attacker; project it ``ENDPOINT_S`` forward; find the closest
   other valid fish to that endpoint as the attacked fish.
6. Classify intensity by raw acceleration.
7. Emit one row per event with attacker/attacked segment-local IDs and
   kinematics. Both raw px/frame^2 and BL-normalised BL/s^2 are emitted.
8. Attach idmatcher columns (attacker_/attacked_fish_id_matched_trial,
   match_quality, match_accepted) via left-merge.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean

from trajectory_io import load_npy_payload
from idmatcher_join import attach_body_length_cm
from zone_assignment import assign_zone, FLOW_LABEL, CALM_LABEL

logger = logging.getLogger(__name__)


BORIS_REQUIRED_COLUMNS = {"Observation id", "Behavior", "Start (s)"}

# Side-annotation behaviours emitted by the new BORIS export (May 2026 file).
SIDE_BEHAVIOR_TO_LABEL = {
    "CALM SIDE": CALM_LABEL,
    "FLOW SIDE": FLOW_LABEL,
}
# Rounding precision for matching aggression rows to side-annotation rows by
# their `Start` value. 1e-6 s = 1 microsecond, well below any frame interval.
SIDE_START_ROUND_DECIMALS = 6


# ---------------------------------------------------------------------------
# Observation-id parsing
# ---------------------------------------------------------------------------

# Format: trial{N}_{treatment}_{condition}_{DD.MM.YY}
_OBS_RE = re.compile(
    r"^trial(?P<trial>\d+)_"
    r"(?P<treatment>[^_]+)_"
    r"(?P<condition>[^_]+)_"
    r"(?P<date>\d{1,2}\.\d{1,2}\.\d{2,4})$",
    re.IGNORECASE,
)


def extract_trial_info(observation_id: str) -> dict[str, Any]:
    """Parse a BORIS observation id into (trial, treatment, condition, date).

    Falls back to a best-effort split on underscores when the canonical
    pattern doesn't match; unparseable IDs return ``None`` for the offending
    fields.
    """
    m = _OBS_RE.match(observation_id.strip())
    if m is not None:
        return {
            "trial": int(m.group("trial")),
            "treatment": m.group("treatment"),
            "condition": m.group("condition"),
            "date": m.group("date"),
            "raw": observation_id,
        }
    # Fallback: try to lift at least the trial number
    parts = observation_id.split("_")
    trial = None
    if parts and parts[0].lower().startswith("trial"):
        try:
            trial = int(parts[0][5:])
        except ValueError:
            trial = None
    return {
        "trial": trial,
        "treatment": parts[1] if len(parts) > 1 else None,
        "condition": parts[2] if len(parts) > 2 else None,
        "date": parts[3] if len(parts) > 3 else None,
        "raw": observation_id,
    }


# ---------------------------------------------------------------------------
# Scope-filter helpers (plan section 11)
# ---------------------------------------------------------------------------

_PHARMA_PATTERNS = ("naloxone", "nloxone", "morphine")


def is_pharma(s: Any) -> bool:
    """Case-insensitive substring match for naloxone / nloxone / morphine."""
    if s is None:
        return False
    try:
        sl = str(s).lower()
    except Exception:
        return False
    return any(p in sl for p in _PHARMA_PATTERNS)


def in_segment1(event_time_abs: float, obs_start_abs: float,
                segment_schedule: list[tuple[float, float]]) -> bool:
    """True if the absolute event time falls in segment-1's window."""
    if not segment_schedule:
        return False
    seg1_start, seg1_end = segment_schedule[0]
    rel = event_time_abs - obs_start_abs
    return seg1_start <= rel < seg1_end


# ---------------------------------------------------------------------------
# BORIS CSV loader
# ---------------------------------------------------------------------------

def load_boris(
    path: Path,
    sheet_name: str | None = None,
    start_column: str = "Start (s)",
) -> pd.DataFrame:
    """Read the BORIS per-event file (CSV or XLSX) with column validation.

    The May-2026 XLSX export uses ``Start`` (not ``Start (s)``) and a sheet
    named ``"Aggregated events"``. If ``start_column != "Start (s)"`` the
    column is renamed in place so the rest of the codebase sees the canonical
    name.

    Parameters
    ----------
    path
        Path to the BORIS export (``.csv`` or ``.xlsx``).
    sheet_name
        Required for ``.xlsx`` files. Ignored for CSV.
    start_column
        Name of the start-time column in the source file. Renamed to
        ``"Start (s)"`` after loading.
    """
    if not path.exists():
        raise FileNotFoundError(f"BORIS_FILE not found: {path}")
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        if sheet_name is None:
            raise ValueError(
                f"sheet_name is required when loading an Excel file ({path})"
            )
        df = pd.read_excel(path, sheet_name=sheet_name)
    else:
        df = pd.read_csv(path)

    # Normalise start-column name
    if start_column != "Start (s)":
        if start_column not in df.columns:
            raise ValueError(
                f"{path} is missing the configured start column "
                f"{start_column!r}. Columns present: {list(df.columns)}"
            )
        df = df.rename(columns={start_column: "Start (s)"})

    missing = BORIS_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required BORIS columns: {sorted(missing)}"
        )
    df["Behavior"] = df["Behavior"].astype(str)
    df["Observation id"] = df["Observation id"].astype(str)
    df["Start (s)"] = df["Start (s)"].astype(float)
    return df


def build_side_annotations(
    boris: pd.DataFrame,
) -> tuple[dict[tuple[str, float], str], pd.DataFrame]:
    """Build the (Observation id, Start) -> side lookup from a BORIS dataframe.

    Each aggression row in ``boris`` is paired with exactly one side-annotation
    row sharing the same ``(Observation id, Start (s))``. Returns:

    * a ``dict`` keyed by ``(observation_id, round(start, 6))`` mapping to
      ``"calm"`` / ``"flow"``;
    * a diagnostic DataFrame with one row per aggression event, columns
      ``Observation id``, ``Start (s)``, ``annotated_side``, ``side_status``
      where ``side_status`` is ``"ok"`` / ``"missing"`` / ``"ambiguous"``.

    The diagnostic DataFrame is what downstream code uses for QC counts; the
    dict is what the per-event loop consults.
    """
    behav = boris["Behavior"].astype(str)
    is_agg = behav.str.lower() == "aggression"
    is_side = behav.isin(SIDE_BEHAVIOR_TO_LABEL.keys())

    side_rows = boris.loc[is_side, ["Observation id", "Start (s)", "Behavior"]].copy()
    side_rows["_key_start"] = side_rows["Start (s)"].round(SIDE_START_ROUND_DECIMALS)
    side_rows["_side"] = side_rows["Behavior"].map(SIDE_BEHAVIOR_TO_LABEL)

    # Detect duplicates / multi-label rows at the same key
    grouped = side_rows.groupby(["Observation id", "_key_start"])["_side"].agg(
        lambda s: list(s.unique())
    )

    lookup: dict[tuple[str, float], str] = {}
    ambiguous_keys: set[tuple[str, float]] = set()
    for key, sides in grouped.items():
        if len(sides) == 1:
            lookup[key] = sides[0]
        else:
            ambiguous_keys.add(key)

    diag_rows: list[dict[str, Any]] = []
    agg_rows = boris.loc[is_agg, ["Observation id", "Start (s)"]].copy()
    for _, ev in agg_rows.iterrows():
        obs = str(ev["Observation id"])
        start_raw = float(ev["Start (s)"])
        key = (obs, round(start_raw, SIDE_START_ROUND_DECIMALS))
        if key in ambiguous_keys:
            status = "ambiguous"
            side = None
        elif key in lookup:
            status = "ok"
            side = lookup[key]
        else:
            status = "missing"
            side = None
        diag_rows.append({
            "Observation id": obs,
            "Start (s)": start_raw,
            "annotated_side": side,
            "side_status": status,
        })

    return lookup, pd.DataFrame(diag_rows)


def get_observation_start(boris: pd.DataFrame, observation_id: str) -> float | None:
    """Return the absolute time of the START marker for one observation."""
    sub = boris[
        (boris["Observation id"] == observation_id)
        & (boris["Behavior"].str.upper() == "START")
    ]
    if sub.empty:
        return None
    return float(sub.iloc[0]["Start (s)"])


# ---------------------------------------------------------------------------
# Kinematics
# ---------------------------------------------------------------------------

def calculate_velocity(positions: np.ndarray, fps: float) -> np.ndarray:
    """Per-frame speed magnitude (px/s). ``positions`` shape: (frames, 2)."""
    diffs = np.diff(positions, axis=0)
    speed = np.linalg.norm(diffs, axis=-1) * fps
    return np.concatenate([[np.nan], speed])  # align to frame index


def calculate_acceleration(velocity: np.ndarray, n_frames: int) -> np.ndarray:
    """Finite-difference acceleration with stride ``n_frames`` (px/frame^2-ish).

    Mirrors block H's helper: ``a[i] = v[i] - v[i - n_frames]``. Output uses
    NaN for indices < n_frames.
    """
    n = velocity.shape[0]
    out = np.full(n, np.nan, dtype=np.float64)
    out[n_frames:] = velocity[n_frames:] - velocity[:n - n_frames]
    return out


def get_trajectory_endpoint(
    positions: np.ndarray, current_frame: int, n_frames_forward: int
) -> np.ndarray | None:
    """Return the position ``n_frames_forward`` frames ahead, if available."""
    n = positions.shape[0]
    end = min(n - 1, current_frame + n_frames_forward)
    if end <= current_frame:
        return None
    candidate = positions[end]
    if np.isnan(candidate).any():
        # walk back until a valid position is found
        for f in range(end - 1, current_frame, -1):
            if not np.isnan(positions[f]).any():
                return positions[f]
        return None
    return candidate


def classify_intensity(max_accel: float, thresholds: dict[str, float]) -> str:
    if max_accel > thresholds["high_intensity"]:
        return "high_intensity"
    if max_accel > thresholds["medium_intensity"]:
        return "medium_intensity"
    if max_accel > thresholds["low_intensity"]:
        return "low_intensity"
    return "mutual_or_ambiguous"


# ---------------------------------------------------------------------------
# Per-segment processing
# ---------------------------------------------------------------------------

def assign_events_for_segment(
    npy_path: Path,
    boris_events: pd.DataFrame,
    obs_start_abs: float,
    trial: int,
    segment_idx: int,                  # 1-based
    segment_schedule: list[tuple[float, float]],
    fps_default: float,
    body_len_default_px: float,
    window_s: float,
    accel_frames: int,
    endpoint_s: float,
    min_valid_fish: int,
    min_valid_frames_in_window: int,
    interaction_thresholds: dict[str, float],
    treatment: str | None,
    condition: str | None,
    side_lookup: dict[tuple[str, float], str] | None = None,
    zone_y0_px: float = 480.0,
    zone_flow_higher_y: bool = True,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Process all BORIS aggression events that fall in this segment's window.

    Hybrid side-aware variant: fish whose position at the event frame matches
    the human-annotated side are preferred (side filter applied when ≥
    ``min_valid_fish`` same-side candidates exist).  When 0 or 1 candidates lie
    on the annotated side the algorithm falls back to the legacy kinematic rule
    applied to *all* valid fish — ``annotated_side`` is set to NA for those
    events and ``side_filter_applied`` is False.  Events with *missing* or
    *ambiguous* side annotations are still hard-skipped (no legacy fallback).
    Pass ``side_lookup=None`` to disable side filtering entirely; the function
    then behaves identically to the non-side pipeline.
    """
    payload = load_npy_payload(npy_path)
    trajectories: np.ndarray = payload["trajectories"]
    fps = float(payload.get("frames_per_second") or fps_default)
    body_length_px = float(payload.get("body_length") or body_len_default_px)
    identities = [str(x) for x in payload.get("identities_labels", [])]
    if not identities:
        identities = [str(i + 1) for i in range(trajectories.shape[1])]

    n_frames = trajectories.shape[0]
    window_frames = int(round(window_s * fps))
    endpoint_frames = int(round(endpoint_s * fps))

    seg_start_rel, seg_end_rel = segment_schedule[segment_idx - 1]
    seg_start_abs = obs_start_abs + seg_start_rel
    seg_end_abs = obs_start_abs + seg_end_rel

    # Pre-compute per-fish velocity + acceleration over the whole segment
    fish_velocity: dict[str, np.ndarray] = {}
    fish_acceleration: dict[str, np.ndarray] = {}
    for fish_idx, fish_id in enumerate(identities):
        positions = trajectories[:, fish_idx, :]
        velocity = calculate_velocity(positions, fps)
        fish_velocity[fish_id] = velocity
        fish_acceleration[fish_id] = calculate_acceleration(velocity, accel_frames)

    results: list[dict[str, Any]] = []
    skipped: dict[str, int] = {
        "out_of_segment": 0,
        "frame_out_of_range": 0,
        "too_few_valid_fish": 0,
        "no_attacker_found": 0,
        "invalid_endpoint": 0,
        "no_attacked_found": 0,
        "missing_positions": 0,
        # side-aware additions
        "side_missing": 0,
        "side_ambiguous": 0,
        "side_no_candidates": 0,
        "side_one_candidate": 0,
        "attacker_position_nan_at_event": 0,
    }

    bl_per_s2 = fps * fps / body_length_px
    side_aware = side_lookup is not None

    for _, ev in boris_events.iterrows():
        t_abs = float(ev["Start (s)"])
        if not (seg_start_abs <= t_abs < seg_end_abs):
            skipped["out_of_segment"] += 1
            continue
        t_within = t_abs - seg_start_abs
        event_frame = int(t_within * fps)
        if event_frame < 0 or event_frame >= n_frames:
            skipped["frame_out_of_range"] += 1
            continue

        # Side-annotation lookup -------------------------------------------
        obs_id = str(ev["Observation id"])
        annotated_side: str | None = None
        side_status = "ok"
        side_annotation_missing = False  # track if annotation was absent
        side_annotation_ambiguous = False  # track if annotation was ambiguous
        if side_aware:
            key = (obs_id, round(t_abs, SIDE_START_ROUND_DECIMALS))
            sentinel = side_lookup.get(key, "__missing__")
            if sentinel == "__missing__":
                # Annotation not found in lookup. Treat as fallback candidate
                # (not hard-skip). We'll process via legacy rule.
                side_annotation_missing = True
                annotated_side = None
                # Continue to candidate-building; will fall back to all valid fish
            else:
                annotated_side = sentinel

        # Per-fish side classification at event frame ----------------------
        fish_sides_at_event: dict[str, Any] = {}
        for fish_idx, fish_id in enumerate(identities):
            ypos = float(trajectories[event_frame, fish_idx, 1])
            fish_sides_at_event[fish_id] = assign_zone(
                ypos, y0_px=zone_y0_px, flow_is_higher_y=zone_flow_higher_y,
            )

        win_start = max(0, event_frame - window_frames)
        win_end = min(n_frames - 1, event_frame + window_frames)

        # Build all-valid fish (pass validity threshold) and side-fish (valid
        # AND on the annotated side).  When side_aware is False all valid fish
        # are candidates (original behaviour).
        all_valid_fish: list[tuple[int, str]] = []
        side_fish: list[tuple[int, str]] = []
        for fish_idx, fish_id in enumerate(identities):
            positions_window = trajectories[win_start:win_end + 1, fish_idx, :]
            valid = ~np.isnan(positions_window).any(axis=1)
            if int(valid.sum()) < min_valid_frames_in_window:
                continue
            all_valid_fish.append((fish_idx, fish_id))
            if not side_aware or fish_sides_at_event[fish_id] == annotated_side:
                side_fish.append((fish_idx, fish_id))

        valid_fish_count = len(all_valid_fish)
        n_side_fish = len(side_fish)  # captured for output column

        # Determine effective candidate set (hybrid rule).
        # Use side-filtered set when ≥ min_valid_fish same-side candidates
        # exist; otherwise fall back to all valid fish (legacy kinematic rule).
        used_legacy_fallback = False
        if side_aware:
            if n_side_fish >= min_valid_fish:
                candidate_fish = side_fish
            else:
                # Hybrid fallback: insufficient same-side candidates.
                used_legacy_fallback = True
                candidate_fish = all_valid_fish
                annotated_side = None  # annotation unused; emit as NA
                if len(candidate_fish) < min_valid_fish:
                    skipped["too_few_valid_fish"] += 1
                    continue
        else:
            candidate_fish = all_valid_fish
            if valid_fish_count < min_valid_fish:
                skipped["too_few_valid_fish"] += 1
                continue

        # Attacker = candidate with highest acceleration in window
        max_accel = -np.inf
        attacker_id: str | None = None
        attacker_idx: int | None = None
        attacker_max_accel_frame = event_frame
        for fish_idx, fish_id in candidate_fish:
            accel_window = fish_acceleration[fish_id][win_start:win_end + 1]
            finite = accel_window[np.isfinite(accel_window)]
            if finite.size == 0:
                continue
            local_max = float(np.nanmax(accel_window))
            if local_max > max_accel:
                max_accel = local_max
                attacker_id = fish_id
                attacker_idx = fish_idx
                arg = int(np.nanargmax(accel_window))
                attacker_max_accel_frame = win_start + arg

        if attacker_id is None:
            skipped["no_attacker_found"] += 1
            continue

        # Verify attacker's side at event frame (hard-filter invariant).
        # Skipped for legacy-fallback events where no side constraint is active.
        if side_aware and not used_legacy_fallback:
            atk_side_at_event = fish_sides_at_event[attacker_id]
            if not isinstance(atk_side_at_event, str):
                skipped["attacker_position_nan_at_event"] += 1
                continue

        attacker_positions = trajectories[:, attacker_idx, :]
        endpoint = get_trajectory_endpoint(
            attacker_positions, attacker_max_accel_frame, endpoint_frames,
        )
        if endpoint is None or np.isnan(endpoint).any():
            skipped["invalid_endpoint"] += 1
            continue

        # Attacked = closest other CANDIDATE fish to the endpoint
        min_distance = np.inf
        attacked_id: str | None = None
        attacked_idx: int | None = None
        for fish_idx, fish_id in candidate_fish:
            if fish_id == attacker_id:
                continue
            fish_pos = trajectories[attacker_max_accel_frame, fish_idx, :]
            if np.isnan(fish_pos).any():
                continue
            d = float(euclidean(fish_pos, endpoint))
            if d < min_distance:
                min_distance = d
                attacked_id = fish_id
                attacked_idx = fish_idx
        if attacked_id is None:
            skipped["no_attacked_found"] += 1
            continue

        # Attacker / attacked velocities at the event frame (px/s)
        atk_vel_px_s = float(fish_velocity[attacker_id][event_frame])
        vct_vel_px_s = float(fish_velocity[attacked_id][event_frame])

        intensity = classify_intensity(max_accel, interaction_thresholds)

        # Effective side where the aggression occurred:
        #   - side-aware events  -> human-annotated side
        #   - legacy-fallback    -> actual zone of the selected attacker
        effective_side = (
            annotated_side if not used_legacy_fallback
            else fish_sides_at_event.get(attacker_id)
        )
        n_calm = sum(1 for s in fish_sides_at_event.values() if s == CALM_LABEL)
        n_flow = sum(1 for s in fish_sides_at_event.values() if s == FLOW_LABEL)

        row: dict[str, Any] = {
            "trial": trial,
            "segment": segment_idx,
            "treatment": treatment,
            "boris_observation": str(ev["Observation id"]),
            "time_absolute_sec": t_abs,
            "time_within_segment_sec": t_within,
            "frame": event_frame,
            "fps": fps,
            "body_length_px": body_length_px,
            "attacker_id": str(attacker_id),
            "attacked_id": str(attacked_id),
            "attacker_max_acceleration_px_per_frame2": max_accel,
            "attacker_max_acceleration_BL_per_s2": (max_accel * bl_per_s2),
            "interaction_type": intensity,
            "attacker_velocity_at_event_px_per_s": atk_vel_px_s,
            "attacker_velocity_at_event_BL_per_s": atk_vel_px_s / body_length_px,
            "attacked_velocity_at_event_px_per_s": vct_vel_px_s,
            "attacked_velocity_at_event_BL_per_s": vct_vel_px_s / body_length_px,
            "distance_at_endpoint_px": min_distance,
            "distance_at_endpoint_BL": min_distance / body_length_px,
            "acceleration_frame": attacker_max_accel_frame,
            "valid_fish_in_window": valid_fish_count,
            # Side-aware diagnostics
            "annotated_side": annotated_side,  # NA when legacy fallback was used
            "side_status": "fallback_legacy" if used_legacy_fallback else side_status,
            "side_filter_applied": not used_legacy_fallback,
            "attacker_side_at_event": fish_sides_at_event.get(attacker_id),
            "attacked_side_at_event": fish_sides_at_event.get(attacked_id),
            # n_side_fish = same-side candidates (0 or 1 for fallback events,
            # showing why the fallback was triggered).
            "n_candidates_on_side": n_side_fish if side_aware else valid_fish_count,
            # Effective side and per-zone fish counts at the event frame
            "side": effective_side,
            "individuals_in_calm": n_calm,
            "individuals_in_flow": n_flow,
        }
        # all-fish positions + sides at the event frame
        for fish_idx, fish_id in enumerate(identities):
            pos = trajectories[event_frame, fish_idx, :]
            row[f"salmon_{fish_id}_x"] = (
                float(pos[0]) if not np.isnan(pos[0]) else None
            )
            row[f"salmon_{fish_id}_y"] = (
                float(pos[1]) if not np.isnan(pos[1]) else None
            )
            side_lbl = fish_sides_at_event.get(fish_id)
            row[f"salmon_{fish_id}_side"] = (
                side_lbl if isinstance(side_lbl, str) else None
            )
        results.append(row)

    df = pd.DataFrame(results)
    return df, skipped


# ---------------------------------------------------------------------------
# idmatcher attachment
# ---------------------------------------------------------------------------

def attach_idmatcher_columns(
    events: pd.DataFrame,
    idmatcher_map: pd.DataFrame | None,
    covered_trial_set: set[int],
    bl_diag: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Left-merge idmatcher map onto attacker_id and attacked_id.

    Adds 6 columns plus optional body-length-cm enrichment. NA matched IDs are
    explained by ``attacker_match_quality`` sentinel values:
      - 'segment1_no_cross_match'
      - 'no_idmatcher_data'
      - 'fish_not_in_idmatcher_map'
      - or the original quality string (incl. 'rejected_conflict').
    """
    out = events.copy()

    new_cols = [
        "attacker_fish_id_matched_trial",
        "attacker_match_quality",
        "attacker_match_accepted",
        "attacked_fish_id_matched_trial",
        "attacked_match_quality",
        "attacked_match_accepted",
    ]
    if idmatcher_map is None:
        for c in new_cols:
            out[c] = pd.NA
        out["attacker_match_quality"] = "idmatcher_disabled"
        out["attacked_match_quality"] = "idmatcher_disabled"
        return out

    def _merge_side(side: str) -> pd.DataFrame:
        suffix = f"{side}_"
        map_renamed = idmatcher_map.rename(columns={
            "fish_id_segment": f"{side}_id",
            "fish_id_matched_trial": f"{suffix}fish_id_matched_trial",
            "match_quality": f"{suffix}match_quality",
            "match_accepted": f"{suffix}match_accepted",
        })
        return out.merge(
            map_renamed,
            on=["trial", "segment", f"{side}_id"],
            how="left",
        )

    out = _merge_side("attacker")
    out = _merge_side("attacked")

    # Fill sentinel quality for rows where merge produced no match
    for side in ("attacker", "attacked"):
        q_col = f"{side}_match_quality"
        m_col = f"{side}_fish_id_matched_trial"
        a_col = f"{side}_match_accepted"

        is_seg1 = out["segment"] == 1
        no_data = ~out["trial"].isin(covered_trial_set)
        was_no_merge = out[q_col].isna()

        # Segment 1 -> always 'segment1_no_cross_match'
        out.loc[is_seg1, q_col] = "segment1_no_cross_match"
        out.loc[is_seg1, m_col] = pd.NA
        out.loc[is_seg1, a_col] = False

        # Trial absent from idmatcher data
        sel = ~is_seg1 & no_data
        out.loc[sel, q_col] = "no_idmatcher_data"
        out.loc[sel, m_col] = pd.NA
        out.loc[sel, a_col] = False

        # Trial present but specific (trial, segment, fish) wasn't mapped
        sel = ~is_seg1 & ~no_data & was_no_merge
        out.loc[sel, q_col] = "fish_not_in_idmatcher_map"
        out.loc[sel, m_col] = pd.NA
        out.loc[sel, a_col] = False

    out = attach_body_length_cm(out, bl_diag)
    return out


# ---------------------------------------------------------------------------
# Scope filtering (applied to event dataframe just after extract_trial_info)
# ---------------------------------------------------------------------------

def apply_scope_filters(
    events: pd.DataFrame,
    *,
    tp1_switch: bool,
    drop_pharma_groups: bool,
    segment_schedule: list[tuple[float, float]],
    obs_start_lookup: dict[str, float],
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Drop rows according to the scope-filter switches.

    Adds two telemetry counts to the returned dict:
      - ``skipped_tp1_switch``: events in seg1 window when TP1_SWITCH is on.
      - ``skipped_pharma_group``: events whose trial parsed to a pharma label.
    """
    info = {"skipped_tp1_switch": 0, "skipped_pharma_group": 0}
    df = events.copy()

    if drop_pharma_groups:
        mask = df.apply(
            lambda r: is_pharma(r.get("treatment")) or is_pharma(r.get("condition")),
            axis=1,
        )
        info["skipped_pharma_group"] = int(mask.sum())
        df = df.loc[~mask].copy()

    if tp1_switch:
        def _is_seg1(r: pd.Series) -> bool:
            obs = r.get("boris_observation") or r.get("Observation id")
            obs_start = obs_start_lookup.get(obs)
            if obs_start is None:
                return False
            return in_segment1(float(r["Start (s)"]), obs_start, segment_schedule)

        mask = df.apply(_is_seg1, axis=1) if not df.empty else pd.Series([], dtype=bool)
        info["skipped_tp1_switch"] = int(mask.sum())
        df = df.loc[~mask].copy()

    return df, info
