import os
import sys
import io
import math
import pickle
import json
import csv
import argparse
from pathlib import Path

import numpy as np

# Try to import idtrackerai if available. When it is not, list_of_blobs.pickle
# is loaded through _SafeUnpickler below instead of a bare pickle.load — see
# docs/EXTRACT_BBOXES_FIX.md section 5 for why a bare pickle.load is not a
# real fallback (it requires idtrackerai to be importable exactly like the
# primary path does).
try:
    from idtrackerai.list_of_blobs import ListOfBlobs
except ImportError:
    ListOfBlobs = None

try:
    import cv2
except ImportError:
    cv2 = None


# -----------------------------------------------------------------------------
# Safe pickle loading (docs/EXTRACT_BBOXES_FIX.md section 5)
# -----------------------------------------------------------------------------

class _Stub:
    """Placeholder for any idtrackerai.* class encountered while unpickling.
    __setstate__ deposits the real pickled attributes onto this instance;
    no idtrackerai code ever runs."""

    def __setstate__(self, state):
        if isinstance(state, dict):
            self.__dict__.update(state)
        else:
            self.__dict__["_state"] = state


class _SafeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.startswith("idtrackerai"):
            return _Stub
        if module.startswith("numpy"):
            return super().find_class(module, name)
        raise pickle.UnpicklingError(
            f"Refusing to unpickle {module}.{name} (not in the allowlist)."
        )


def load_list_of_blobs_safely(path):
    with open(path, "rb") as f:
        return _SafeUnpickler(f).load()


def resolve_pickle_path(session_path: Path):
    """Pick the base list_of_blobs.pickle, skipping macOS resource forks
    (._list_of_blobs.pickle) that can sit next to the real file on a
    OneDrive/iCloud-synced corpus. Deliberately does NOT consider
    list_of_blobs_validated.pickle -- see docs/EXTRACT_BBOXES_FIX.md."""
    p = session_path / "preprocessing" / "list_of_blobs.pickle"
    if p.exists() and p.stat().st_size > 0 and not p.name.startswith("._"):
        return p
    return None


def load_list_of_blobs(path: Path):
    """Returns (list_of_blobs_obj, loader_name) where loader_name is
    'idtrackerai' (real classes, methods available) or 'stub' (state-only,
    no idtrackerai import). Both paths expose the same pickled attributes."""
    path = str(path)
    if ListOfBlobs is not None and hasattr(ListOfBlobs, "load"):
        try:
            return ListOfBlobs.load(path), "idtrackerai"
        except Exception:
            pass

    return load_list_of_blobs_safely(path), "stub"


# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------

def bbox_diag_from_contour(contour):
    """Diagonal of the inclusive bounding box, matching idtracker.ai's own
    `extension` property exactly (verified: hypot(np.ptp(contour, axis=0))
    equals blob.extension to the ULP on a real corpus). Used unconditionally
    -- for both the live idtrackerai loader and the stub loader -- so a run
    with idtrackerai installed and a run without it produce identical output.

    NOTE: np.ndarray.ptp was removed in NumPy 2.0; must call np.ptp(arr, axis).
    """
    c = np.asarray(contour)
    if c.ndim != 2 or c.shape[0] == 0:
        return None
    d = np.ptp(c, axis=0)
    return float(math.hypot(float(d[0]), float(d[1])))


def area_from_contour(contour):
    """Returns (area_px2, method). cv2.contourArea when available, otherwise
    the shoelace formula (same result for a simple closed polygon)."""
    c = np.asarray(contour, dtype=np.float64)
    if c.ndim != 2 or c.shape[0] < 3:
        return None, None
    if cv2 is not None:
        try:
            return float(cv2.contourArea(c.astype(np.float32))), "cv2.contourArea"
        except Exception:
            pass
    x = c[:, 0]
    y = c[:, 1]
    shoelace = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return float(shoelace), "shoelace"


# -----------------------------------------------------------------------------
# Identity resolution (docs/EXTRACT_BBOXES_FIX.md N2)
# -----------------------------------------------------------------------------

def _clean_ident_list(raw, n_animals):
    """Coerce a raw identity-like value (scalar, list, generator, 0-d numpy
    array) to a list of valid 1..n_animals ints. Drops None, non-integers,
    <=0 (the -1 "not assigned" sentinel and 0), and > n_animals."""
    if raw is None:
        return []
    if isinstance(raw, (int, np.integer)):
        items = [raw]
    elif isinstance(raw, (str, bytes)):
        items = [raw]
    else:
        try:
            items = list(raw)
        except TypeError:
            items = [raw]

    out = []
    for v in items:
        if v is None:
            continue
        try:
            iv = int(v)
        except (TypeError, ValueError):
            continue
        if iv <= 0:
            continue
        if n_animals is not None and iv > n_animals:
            continue
        out.append(iv)
    return out


def resolve_identities(blob, n_animals):
    """Returns (identities, source) where source is one of
    'all_final_identities', 'user_generated_identities', 'identity', 'none'.
    Deliberately never reads blob.final_identities: on idtrackerai 6.0.14 it
    is a generator that is always truthy, so `if blob.final_identities` is
    dead code, not a real fallback gate -- see docs/EXTRACT_BBOXES_FIX.md."""
    aid = getattr(blob, "all_final_identities", None)
    ids = _clean_ident_list(aid, n_animals)
    if ids:
        return ids, "all_final_identities"

    ugi = getattr(blob, "user_generated_identities", None)
    ids = _clean_ident_list(ugi, n_animals)
    if ids:
        return ids, "user_generated_identities"

    ident = getattr(blob, "identity", None)
    ids = _clean_ident_list(ident, n_animals)
    if ids:
        return ids, "identity"

    return [], "none"


def get_n_animals(session_path: Path, attrs: dict):
    """n_animals from session.json (present in every real session sampled),
    falling back to attributes.json's identities_labels (attributes.json has
    no number_of_animals key)."""
    session_json = session_path / "session.json"
    if session_json.exists():
        try:
            with open(session_json, "r", encoding="utf-8") as f:
                sess = json.load(f)
            for key in ("n_animals", "number_of_animals"):
                v = sess.get(key)
                if isinstance(v, (int, float)) and v > 0:
                    return int(v), float(sess.get("median_body_length")) if sess.get("median_body_length") else None
        except Exception:
            pass

    labels = attrs.get("identities_labels")
    if isinstance(labels, list) and labels:
        return len(labels), None

    return None, None


# -----------------------------------------------------------------------------
# Main extraction
# -----------------------------------------------------------------------------

def process_session(session_path, output_dir, min_certainty=0.0):
    session_path = Path(session_path)
    session_id = session_path.name
    output_dir = Path(output_dir)

    lob_path = resolve_pickle_path(session_path)
    attrs_path = session_path / "trajectories" / "trajectories_csv" / "attributes.json"

    summary = {
        "session_id": session_id,
        "status": "skipped",
        "reason": "",
        "n_frames_in_blobs": 0,
        "n_rows_written": 0,
        "metadata": {},
        "pickle_file": None,
        "loader": None,
    }

    # 0. Check requirements
    if lob_path is None:
        summary["reason"] = "Missing list_of_blobs.pickle"
        return summary

    if not attrs_path.exists():
        summary["reason"] = "Missing attributes.json"
        return summary

    # 1. Load metadata
    try:
        with open(attrs_path, "r", encoding="utf-8") as f:
            attrs = json.load(f)

        length_unit = attrs.get("length_unit", 0)
        fps = attrs.get("frames_per_second", None)

        try:
            length_unit = float(length_unit)
            if not (length_unit > 0):
                length_unit = None
        except Exception:
            length_unit = None

        n_animals, tracker_median_bl = get_n_animals(session_path, attrs)
        if n_animals is None:
            summary["reason"] = "Could not determine n_animals from session.json or attributes.json"
            return summary

        summary["metadata"] = {
            "length_unit": length_unit,
            "fps": fps,
            "n_animals": n_animals,
            "min_certainty": min_certainty,
        }

    except Exception as e:
        summary["reason"] = f"Error loading attributes: {e}"
        return summary

    # 2. Load blobs
    try:
        lob, loader = load_list_of_blobs(lob_path)
        blobs_in_video = lob.blobs_in_video
        n_frames = len(blobs_in_video)  # NOT lob.number_of_frames -- a property, absent under stub-loading
        summary["n_frames_in_blobs"] = int(n_frames)
        summary["pickle_file"] = lob_path.name
        summary["loader"] = loader
    except Exception as e:
        summary["reason"] = f"Error loading blobs: {e}"
        return summary

    # 3. Iterate frames
    candidate_rows = []

    exclusion_counts = {
        "no_identity": 0,
        "merged_multi_identity": 0,
        "not_unicity_frame": 0,
        "not_seems_like_individual": 0,
        "is_a_crossing": 0,
        "missing_bbox": 0,
        "below_min_certainty": 0,
    }
    identity_source_counts = {
        "all_final_identities": 0,
        "user_generated_identities": 0,
        "identity": 0,
        "none": 0,
    }
    certainty_lt_05 = 0
    certainty_lt_005 = 0
    n_unicity_frames = 0

    for f in range(n_frames):
        blobs_in_frame = blobs_in_video[f]
        is_unicity_frame = (len(blobs_in_frame) == n_animals)
        if is_unicity_frame:
            n_unicity_frames += 1

        for blob in blobs_in_frame:
            identities, id_source = resolve_identities(blob, n_animals)
            identity_source_counts[id_source] += 1

            if not identities:
                exclusion_counts["no_identity"] += 1
                continue

            if len(identities) > 1:
                exclusion_counts["merged_multi_identity"] += 1
                continue

            ident = identities[0]

            if not is_unicity_frame:
                exclusion_counts["not_unicity_frame"] += 1
                continue

            seems_like_individual = getattr(blob, "seems_like_individual", None)
            if seems_like_individual is None:
                exclusion_counts["not_seems_like_individual"] += 1
                continue
            if not bool(seems_like_individual):
                exclusion_counts["not_seems_like_individual"] += 1
                continue

            is_crossing = getattr(blob, "is_a_crossing", None)
            if is_crossing is not None and bool(is_crossing):
                exclusion_counts["is_a_crossing"] += 1
                continue

            contour = getattr(blob, "contour", None)
            body_length_px = bbox_diag_from_contour(contour) if contour is not None else None
            if body_length_px is None:
                exclusion_counts["missing_bbox"] += 1
                continue

            certainty = getattr(blob, "identity_certainty", None)
            try:
                certainty = float(certainty) if certainty is not None else None
                if certainty is not None and math.isnan(certainty):
                    certainty = None
            except (TypeError, ValueError):
                certainty = None

            if certainty is not None:
                if certainty < 0.5:
                    certainty_lt_05 += 1
                if certainty < 0.05:
                    certainty_lt_005 += 1

            if min_certainty > 0.0 and certainty is not None and certainty < min_certainty:
                exclusion_counts["below_min_certainty"] += 1
                continue

            area_px2, area_method = area_from_contour(contour)

            body_length_cm = None
            if body_length_px is not None and length_unit is not None:
                body_length_cm = body_length_px / length_unit

            candidate_rows.append({
                "session_id": session_id,
                "frame": int(f),
                "identity": int(ident),
                "body_length_px": body_length_px,
                "body_length_cm": body_length_cm,
                "is_individual": True,
                "is_crossing": False,
                "identity_certainty": certainty,
                "area_px2": area_px2,
                "identity_source": id_source,
            })

    # 4. Write CSV
    output_dir.mkdir(parents=True, exist_ok=True)

    if not candidate_rows:
        summary["status"] = "completed_empty"
        summary["reason"] = "No valid rows extracted"
    else:
        candidate_rows.sort(key=lambda x: (x["frame"], x["identity"]))

        out_name = f"{session_id}_bboxes.csv"
        out_path = output_dir / out_name

        fieldnames = [
            "session_id", "frame", "identity",
            "body_length_px", "body_length_cm",
            "is_individual", "is_crossing",
            "identity_certainty", "area_px2", "identity_source",
        ]

        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(candidate_rows)

        summary["n_rows_written"] = len(candidate_rows)
        summary["status"] = "completed"

        # Per-identity stats -- aggregate body_length_px (always present when
        # a row was written at all), convert to cm only if length_unit holds.
        # docs/EXTRACT_BBOXES_FIX.md N6: aggregating the already-converted
        # (and usually absent) cm column made identity_stats report
        # {"n_frames": 0} on every real session in this corpus.
        rows_by_id = {}
        for r in candidate_rows:
            rows_by_id.setdefault(r["identity"], []).append(r)

        ident_stats = {}
        perid_medians_px = {}
        for ident, rows in rows_by_id.items():
            px_vals = np.array([r["body_length_px"] for r in rows], dtype=float)
            cm_vals = np.array(
                [r["body_length_cm"] for r in rows if r["body_length_cm"] is not None],
                dtype=float,
            )
            q25, q50, q75 = np.percentile(px_vals, [25, 50, 75])
            perid_medians_px[ident] = float(q50)
            stats = {
                "n_frames": int(len(px_vals)),
                "pct_valid": round(len(px_vals) / float(n_frames) * 100, 2),
                "median_body_length_px": round(float(q50), 4),
                "iqr_body_length_px": round(float(q75 - q25), 4),
            }
            if len(cm_vals) > 0:
                cq25, cq50, cq75 = np.percentile(cm_vals, [25, 50, 75])
                stats["median_body_length_cm"] = round(float(cq50), 4)
                stats["iqr_body_length_cm"] = round(float(cq75 - cq25), 4)
            ident_stats[int(ident)] = stats

        summary["identity_stats"] = ident_stats

        all_px = np.array([r["body_length_px"] for r in candidate_rows], dtype=float)
        median_all_px = float(np.median(all_px))
        summary["median_body_length_px_all"] = round(median_all_px, 4)
        summary["tracker_median_body_length"] = tracker_median_bl
        if tracker_median_bl:
            summary["reconciliation_pct"] = round(
                100.0 * (median_all_px / tracker_median_bl - 1.0), 3
            )
        else:
            summary["reconciliation_pct"] = None

        if len(perid_medians_px) >= 2:
            mx = max(perid_medians_px.values())
            mn = min(perid_medians_px.values())
            summary["perid_median_ratio"] = round(mx / mn, 3) if mn > 0 else None
        else:
            summary["perid_median_ratio"] = None

    summary["n_unicity_frames"] = n_unicity_frames
    summary["pct_unicity_frames"] = round(100.0 * n_unicity_frames / n_frames, 2) if n_frames else 0.0
    summary["exclusion_counts"] = exclusion_counts
    summary["identity_source_counts"] = identity_source_counts
    summary["pct_certainty_lt_0.5"] = round(
        100.0 * certainty_lt_05 / len(candidate_rows), 2
    ) if candidate_rows else None
    summary["pct_certainty_lt_0.05"] = round(
        100.0 * certainty_lt_005 / len(candidate_rows), 2
    ) if candidate_rows else None

    # 5. Write summary JSON
    sum_name = f"{session_id}_bboxes_summary.json"
    with open(output_dir / sum_name, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    ap = argparse.ArgumentParser(description="Extract per-frame bbox proxy lengths from idtracker.ai list_of_blobs.pickle")
    ap.add_argument("--session_dir", required=True, help="Path to one session_* folder")
    ap.add_argument("--out_dir", required=True, help="Folder where to write <session_id>_bboxes.csv and summary")
    ap.add_argument("--min_certainty", type=float, default=0.0,
                     help="Drop rows with identity_certainty below this value (default 0.0 = keep all; certainty is always recorded and reported regardless)")
    args = ap.parse_args()

    session = Path(args.session_dir)
    out_dir = Path(args.out_dir)

    if not session.exists() or not session.is_dir():
        print(f"Session not found or not a folder: {session}", file=sys.stderr)
        sys.exit(2)

    print(f"Processing {session.name}...")
    res = process_session(session, out_dir, min_certainty=args.min_certainty)
    print(f"  -> {res['status']} ({res['n_rows_written']} rows)")
    if res.get("reason"):
        print(f"     Reason: {res['reason']}")
    if res.get("reconciliation_pct") is not None:
        print(f"     Reconciliation vs tracker median_body_length: {res['reconciliation_pct']}%")

    # non-zero only for real failures
    if res["status"] not in ("completed", "completed_empty"):
        sys.exit(3)


if __name__ == "__main__":
    main()
