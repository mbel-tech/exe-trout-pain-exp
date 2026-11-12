"""Trajectory jump-protection and PCHIP gap interpolation.

Adapted from block C STEP 3 of the original Colab notebook
(``interpolate_trajectory`` at line 1064 of the source script). Operates on
arrays of shape ``(frames, fish, 2)``.

The flow per fish coordinate:

1. **Jump protection** -- compute a rolling median over ``jump_window`` frames;
   flag points whose deviation from the rolling median exceeds
   ``jump_threshold`` as NaN.
2. **PCHIP fill** -- for gaps of length <= ``max_gap`` frames, fill with
   :class:`scipy.interpolate.PchipInterpolator` on float64 (the float32 path
   in older scipy fitpack can raise on highly-clustered knots).
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator

logger = logging.getLogger(__name__)


def interpolate_trajectory(
    traj: np.ndarray,
    max_gap: int,
    remove_jumps: bool = True,
    jump_threshold: float = 80.0,
    jump_window: int = 11,
) -> np.ndarray:
    """Interpolate per-fish trajectory gaps with optional jump removal.

    Parameters
    ----------
    traj:
        Array of shape ``(frames, fish, 2)``.
    max_gap:
        Largest gap (in frames) that will be filled. Longer gaps are left as
        NaN.
    remove_jumps:
        When True, points that deviate from the rolling median by more than
        ``jump_threshold`` are first replaced with NaN.
    jump_threshold:
        Pixel distance threshold used by the rolling-median jump detector.
    jump_window:
        Window length (frames) for the rolling-median jump detector.

    Returns
    -------
    np.ndarray
        Interpolated copy of ``traj`` (float64).
    """
    if traj.ndim != 3 or traj.shape[-1] != 2:
        raise ValueError(f"expected (frames, fish, 2), got {traj.shape}")
    out = traj.astype(np.float64, copy=True)
    n_frames, n_fish, _ = out.shape

    for fish_idx in range(n_fish):
        for coord in (0, 1):
            series = out[:, fish_idx, coord]
            if remove_jumps:
                series = _flag_jumps(series, threshold=jump_threshold,
                                     window=jump_window)
            series = _pchip_fill(series, max_gap=max_gap)
            out[:, fish_idx, coord] = series
    return out


def _rolling_median(x: np.ndarray, window: int) -> np.ndarray:
    """Centred rolling median that propagates NaN through gaps."""
    if window < 3:
        return x.copy()
    half = window // 2
    n = x.shape[0]
    pad = np.concatenate([
        np.full(half, np.nan),
        x,
        np.full(half, np.nan),
    ])
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        seg = pad[i:i + window]
        out[i] = np.nanmedian(seg) if np.any(~np.isnan(seg)) else np.nan
    return out


def _flag_jumps(series: np.ndarray, threshold: float, window: int) -> np.ndarray:
    """Mark points >threshold from the rolling median as NaN."""
    rolling = _rolling_median(series, window)
    diff = np.abs(series - rolling)
    bad = diff > threshold
    if bad.any():
        series = series.copy()
        series[bad] = np.nan
    return series


def _pchip_fill(series: np.ndarray, max_gap: int) -> np.ndarray:
    """Fill NaN gaps up to ``max_gap`` frames using PCHIP interpolation."""
    n = series.shape[0]
    valid = ~np.isnan(series)
    if valid.sum() < 2:
        return series
    valid_idx = np.nonzero(valid)[0]
    valid_vals = series[valid_idx].astype(np.float64)
    try:
        interp = PchipInterpolator(valid_idx, valid_vals, extrapolate=False)
    except ValueError as exc:
        logger.warning("PCHIP failed (%s); leaving series unchanged", exc)
        return series

    out = series.copy()
    # Iterate over NaN runs; only fill if run length <= max_gap and
    # bracketed by valid samples (no extrapolation).
    i = 0
    while i < n:
        if not np.isnan(out[i]):
            i += 1
            continue
        j = i
        while j < n and np.isnan(out[j]):
            j += 1
        gap_len = j - i
        if gap_len <= max_gap and i > 0 and j < n:
            idx = np.arange(i, j)
            filled = interp(idx)
            mask = ~np.isnan(filled)
            out[idx[mask]] = filled[mask]
        i = j
    return out


def summarize_interpolation(before: np.ndarray, after: np.ndarray) -> dict[str, Any]:
    """Compact summary for logging."""
    return {
        "n_frames": int(before.shape[0]),
        "n_fish": int(before.shape[1]),
        "nan_before": int(np.isnan(before).any(axis=-1).sum()),
        "nan_after": int(np.isnan(after).any(axis=-1).sum()),
    }
