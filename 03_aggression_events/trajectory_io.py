"""Trajectory file I/O: H5/NPY reading, NPY writing, and identity-label cleanup.

Reads idtracker.ai outputs (either ``trajectories.h5`` or ``trajectories.npy``)
and writes the canonical "fat NPY" used by the rest of the pipeline:
a zero-d object array containing a dict with ``trajectories`` (frames, fish, 2),
``frames_per_second``, ``body_length``, ``identities_labels``, and any other
metadata available.

Identity-label cleanup
----------------------
idtracker.ai sometimes emits 5 labels (e.g. ``['0','1','2','3','4']``) when only
4 fish are in the tank. The extra label is almost always ``'0'`` but can be any
identity that was barely tracked. :func:`clean_identity_labels` reduces the set
to the canonical four by:

1. Preferring to drop ``'0'`` if it's present alongside 4 other labels.
2. Otherwise dropping the identity whose valid-position fraction across all
   frames is lowest.

A short log line records the original label set, the dropped label, and the
reason. Any file that cannot be reduced to exactly four identities raises a
``RuntimeError`` rather than guessing.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import h5py
import numpy as np

logger = logging.getLogger(__name__)


# Keys the script will look for inside an H5 file when extracting metadata.
H5_TRAJECTORY_KEYS = ("trajectories", "tracking", "tracks", "positions", "data")
H5_IDENTITY_KEYS = ("identities_labels", "identities", "ids", "labels", "animal_ids")
H5_FPS_KEYS = ("frames_per_second", "fps", "frame_rate", "framerate")
H5_BODY_LEN_KEYS = ("body_length", "body_length_px", "bodylength")

# Trial / segment parsing.
# Matches either "session_trial5_Segment2" (raw session folder) or
# "trajectories_trial5_Segment2" (staged file), case-insensitive.
_SESSION_RE = re.compile(r"trial(\d+)_Segment(\d+)", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def parse_session_name(path_or_name: str | Path) -> tuple[int, int] | None:
    """Return ``(trial_number, segment_number)`` parsed from a session name."""
    name = Path(path_or_name).name if isinstance(path_or_name, Path) else str(path_or_name)
    m = _SESSION_RE.search(name)
    if m is None:
        return None
    return int(m.group(1)), int(m.group(2))


def canonical_npy_name(trial: int, segment: int) -> str:
    return f"trajectories_trial{trial}_Segment{segment}.npy"


# ---------------------------------------------------------------------------
# Identity cleanup
# ---------------------------------------------------------------------------

CANONICAL_LABELS = ("1", "2", "3", "4")


def clean_identity_labels(
    labels: list[str],
    trajectories: np.ndarray,
    file_label: str = "<unknown>",
    canonical: tuple[str, ...] = CANONICAL_LABELS,
) -> tuple[list[str], np.ndarray, dict[str, Any]]:
    """Reduce the identity-label set so it is a subset of ``canonical``.

    Rules:

    1. Already a subset and matches order -> noop.
    2. Already a subset but wrong order -> reorder.
    3. Contains ``'0'`` plus a 4-fish canonical superset -> drop ``'0'``.
    4. Has more labels than canonical -> drop the lowest-valid-fraction extras
       until what remains is a subset of canonical.
    5. Resulting set must have at least two identities; otherwise raises.

    Returns the labels in their original numerical order (a subset of
    ``canonical``).
    """
    labels = [str(x) for x in labels]
    original = list(labels)
    canon_set = set(canonical)

    valid_frac = _valid_fractions(trajectories, labels)

    # Case A: already a subset of canonical (1..4)
    if set(labels) <= canon_set:
        # Reorder by canonical sequence (drop missing entries naturally)
        new_labels = [c for c in canonical if c in labels]
        order = [labels.index(c) for c in new_labels]
        new_traj = trajectories[:, order, :]
        reason = (
            "already_canonical" if new_labels == labels
            else "reordered_to_canonical_subset"
        )
        info = {
            "original_labels": original,
            "dropped": None,
            "reason": reason,
            "valid_fraction": valid_frac,
        }
        if len(new_labels) < 2:
            raise RuntimeError(
                f"[{file_label}] only {len(new_labels)} canonical identities "
                f"present ({new_labels!r}); need at least 2"
            )
        if reason != "already_canonical":
            logger.info(
                "[%s] identity-label cleanup: %r -> %r (%s)",
                file_label, original, new_labels, reason,
            )
        return new_labels, new_traj, info

    # Case B: extras present (e.g. '0' or '5'). Drop until remainder is a
    # subset of canonical.
    extras = [lbl for lbl in labels if lbl not in canon_set]

    # Preference: '0' first if present
    drop_order: list[str] = []
    if "0" in extras:
        drop_order.append("0")
    # then other extras by lowest valid fraction
    other_extras = sorted(
        (lbl for lbl in extras if lbl != "0"),
        key=lambda lbl: valid_frac.get(lbl, 0.0),
    )
    drop_order.extend(other_extras)

    # If we still have more than len(canonical) labels after dropping extras,
    # drop further by lowest valid fraction.
    after = [lbl for lbl in labels if lbl not in set(drop_order)]
    if len(after) > len(canonical):
        # rank remaining by valid fraction and drop the worst
        extra_ranked = sorted(after, key=lambda lbl: valid_frac.get(lbl, 0.0))
        n_extra = len(after) - len(canonical)
        drop_order.extend(extra_ranked[:n_extra])
        after = [lbl for lbl in labels if lbl not in set(drop_order)]

    if not set(after) <= canon_set:
        raise RuntimeError(
            f"[{file_label}] after dropping {drop_order!r}, surviving labels "
            f"{after!r} are not a subset of canonical {list(canonical)!r}. "
            f"Original: {original!r}, valid fractions: {valid_frac!r}"
        )

    new_labels = [c for c in canonical if c in after]
    order = [labels.index(c) for c in new_labels]
    new_traj = trajectories[:, order, :]
    reason = "drop_zero" if drop_order == ["0"] else "drop_extras_by_valid_fraction"
    info = {
        "original_labels": original,
        "dropped": ",".join(drop_order),
        "reason": reason,
        "valid_fraction": valid_frac,
    }
    logger.info(
        "[%s] identity-label cleanup: %r -> %r (dropped %r, %s, "
        "valid_fraction=%s)",
        file_label, original, new_labels, drop_order, reason,
        {k: round(v, 3) for k, v in valid_frac.items()},
    )
    if len(new_labels) < 2:
        raise RuntimeError(
            f"[{file_label}] cleanup left {len(new_labels)} identities; "
            f"need at least 2"
        )
    return new_labels, new_traj, info


def _valid_fractions(trajectories: np.ndarray, labels: list[str]) -> dict[str, float]:
    """Return the fraction of frames with non-NaN positions, per identity."""
    out: dict[str, float] = {}
    n_frames = trajectories.shape[0]
    for i, lbl in enumerate(labels):
        positions = trajectories[:, i, :]
        valid = ~np.isnan(positions).any(axis=1)
        out[lbl] = float(valid.sum() / max(n_frames, 1))
    return out


# ---------------------------------------------------------------------------
# NPY <-> H5
# ---------------------------------------------------------------------------

def load_npy_payload(npy_path: Path) -> dict[str, Any]:
    """Load the canonical fat-NPY format used by this pipeline.

    Returns a dict with at minimum ``trajectories``, ``frames_per_second``,
    ``body_length``, ``identities_labels``.
    """
    data = np.load(npy_path, allow_pickle=True)
    if data.dtype == object and data.ndim == 0:
        return data.item()
    raise RuntimeError(
        f"Unexpected NPY format at {npy_path}: dtype={data.dtype}, ndim={data.ndim}"
    )


def save_npy_payload(npy_path: Path, payload: dict[str, Any]) -> None:
    npy_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(npy_path, np.array(payload, dtype=object), allow_pickle=True)


def convert_h5_to_npy_payload(
    h5_path: Path,
    fps_default: float,
    body_len_default_px: float,
) -> dict[str, Any]:
    """Read an idtracker.ai H5 file into the canonical fat-NPY dict."""
    with h5py.File(h5_path, "r") as f:
        traj = _find_3d_array(f)
        payload: dict[str, Any] = {"trajectories": traj}

        for key in H5_FPS_KEYS:
            if key in f.attrs:
                payload["frames_per_second"] = float(f.attrs[key])
                break
            if key in f:
                payload["frames_per_second"] = float(np.array(f[key]).item())
                break
        if "frames_per_second" not in payload:
            payload["frames_per_second"] = fps_default
            logger.warning("[%s] fps not found; using default %s",
                           h5_path.name, fps_default)

        for key in H5_BODY_LEN_KEYS:
            if key in f.attrs:
                payload["body_length"] = float(f.attrs[key])
                break
            if key in f:
                payload["body_length"] = float(np.array(f[key]).item())
                break
        if "body_length" not in payload:
            payload["body_length"] = body_len_default_px
            logger.warning("[%s] body_length not found; using default %s px",
                           h5_path.name, body_len_default_px)

        for key in H5_IDENTITY_KEYS:
            if key in f.attrs:
                payload["identities_labels"] = [
                    str(x) for x in np.atleast_1d(f.attrs[key]).tolist()
                ]
                break
            if key in f:
                payload["identities_labels"] = [
                    str(x) for x in np.array(f[key]).tolist()
                ]
                break
        if "identities_labels" not in payload:
            n_animals = traj.shape[1]
            payload["identities_labels"] = [str(i + 1) for i in range(n_animals)]
            logger.warning("[%s] identities_labels not found; defaulting to %r",
                           h5_path.name, payload["identities_labels"])

    return payload


def _find_3d_array(h5_file: h5py.File) -> np.ndarray:
    """Locate the (frames, animals, 2) array inside an H5 file."""
    for key in H5_TRAJECTORY_KEYS:
        if key in h5_file and h5_file[key].ndim == 3 and h5_file[key].shape[-1] == 2:
            return np.array(h5_file[key], dtype=np.float64)

    candidates: list[tuple[str, tuple[int, ...]]] = []

    def visit(name: str, obj: h5py.HLObject) -> None:
        if isinstance(obj, h5py.Dataset) and obj.ndim == 3 and obj.shape[-1] == 2:
            candidates.append((name, obj.shape))

    h5_file.visititems(visit)
    if not candidates:
        raise RuntimeError(
            f"No (frames, animals, 2) dataset found in {h5_file.filename}"
        )
    # Prefer the largest by frame count
    name, _ = max(candidates, key=lambda c: c[1][0])
    logger.info("[%s] using dataset %s for trajectories",
                Path(h5_file.filename).name, name)
    return np.array(h5_file[name], dtype=np.float64)


# ---------------------------------------------------------------------------
# Session-folder discovery / staging
# ---------------------------------------------------------------------------

def discover_session_files(
    sessions_dir: Path, excluded_trials: list[int] | None = None
) -> list[dict[str, Any]]:
    """Walk ``sessions_dir`` and return a list of session entries.

    Each entry is a dict with:
        - ``trial`` (int)
        - ``segment`` (int)
        - ``session_dir`` (Path)
        - ``npy_path`` (Path | None) -- prefer this if present
        - ``h5_path`` (Path | None) -- fallback if NPY missing

    Parameters
    ----------
    sessions_dir
        Root folder containing session subfolders.
    excluded_trials
        Trial numbers to skip entirely. Defaults to None (no exclusions).
    """
    excluded_trials = excluded_trials or []
    entries: list[dict[str, Any]] = []
    if not sessions_dir.exists():
        raise FileNotFoundError(f"SESSIONS_DIR does not exist: {sessions_dir}")
    for sub in sorted(sessions_dir.iterdir()):
        if not sub.is_dir():
            continue
        parsed = parse_session_name(sub.name)
        if parsed is None:
            continue
        trial, segment = parsed
        if trial in excluded_trials:
            logger.info("[%s] trial %d in excluded list; skipping", sub.name, trial)
            continue
        traj_dir = sub / "trajectories"
        npy_path = traj_dir / "trajectories.npy"
        h5_path = traj_dir / "trajectories.h5"
        if not npy_path.exists() and not h5_path.exists():
            logger.warning("[%s] no trajectories.{npy,h5}; skipping", sub.name)
            continue
        entries.append({
            "trial": trial,
            "segment": segment,
            "session_dir": sub,
            "npy_path": npy_path if npy_path.exists() else None,
            "h5_path": h5_path if h5_path.exists() else None,
        })
    return entries


def stage_session(
    entry: dict[str, Any],
    out_dir: Path,
    fps_default: float,
    body_len_default_px: float,
    overwrite: bool = False,
) -> tuple[Path, dict[str, Any]]:
    """Copy or convert the NPY for one session into ``out_dir``.

    Returns
    -------
    tuple
        ``(out_path, cleanup_info)``.
    """
    out_path = out_dir / canonical_npy_name(entry["trial"], entry["segment"])
    if out_path.exists() and not overwrite:
        payload = load_npy_payload(out_path)
        labels = [str(x) for x in payload.get("identities_labels", [])]
        info = {
            "original_labels": labels,
            "dropped": None,
            "reason": "already_staged",
            "valid_fraction": {},
        }
        return out_path, info

    if entry["npy_path"] is not None:
        payload = load_npy_payload(entry["npy_path"])
        source = entry["npy_path"]
    else:
        payload = convert_h5_to_npy_payload(
            entry["h5_path"],
            fps_default=fps_default,
            body_len_default_px=body_len_default_px,
        )
        source = entry["h5_path"]

    labels = [str(x) for x in payload.get("identities_labels", [])]
    traj = payload["trajectories"]
    new_labels, new_traj, info = clean_identity_labels(
        labels, traj, file_label=out_path.name,
    )
    payload["trajectories"] = new_traj
    payload["identities_labels"] = new_labels
    payload["identity_cleanup_info"] = info
    payload["source_file"] = str(source)
    save_npy_payload(out_path, payload)
    return out_path, info
