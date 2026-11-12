"""Calm / Flow zone assignment from raw pixel y-coordinate.

Ported from the R main pipeline
``01_core_pipeline/01_idtracker_pipeline_GOT.R::zone_occupancy_y_boundary``
(lines 1408-1443). The rule is deterministic: a single horizontal pixel-y
threshold (``y0_px``, default 480) splits the arena into a low-y half and a
high-y half. ``flow_is_higher_y`` selects which half is labelled "flow".

Defaults match the R-side defaults
(``DEFAULT_Y0_PX = mean(c(483, 477)) = 480`` and
``DEFAULT_FLOW_HIGHER_Y = TRUE``). No per-trial calibration is required;
no margin buffer is applied (matching the R function's actual behaviour,
which also does not enforce the documented `DEFAULT_MARGIN_CM = 10`).
"""

from __future__ import annotations

from typing import Any

import numpy as np

FLOW_LABEL = "flow"
CALM_LABEL = "calm"
Y_THRESHOLD_DEFAULT = 480.0
FLOW_IS_HIGHER_Y_DEFAULT = True


def assign_zone(
    y: Any,
    y0_px: float = Y_THRESHOLD_DEFAULT,
    flow_is_higher_y: bool = FLOW_IS_HIGHER_Y_DEFAULT,
) -> Any:
    """Classify a y-coordinate (or array of them) as ``"flow"`` / ``"calm"``.

    Parameters
    ----------
    y
        Scalar float or numpy array of pixel y-coordinates. ``NaN`` propagates.
    y0_px
        Threshold separating the two zones (default 480).
    flow_is_higher_y
        If ``True`` (default), ``y >= y0_px`` is ``"flow"`` and ``y < y0_px``
        is ``"calm"``. If ``False``, the labels are swapped.

    Returns
    -------
    str or numpy.ndarray of dtype object
        The label(s). ``NaN`` inputs return ``np.nan``.
    """
    arr = np.asarray(y, dtype=float)
    scalar = arr.ndim == 0

    if flow_is_higher_y:
        high_label, low_label = FLOW_LABEL, CALM_LABEL
    else:
        high_label, low_label = CALM_LABEL, FLOW_LABEL

    if scalar:
        v = float(arr)
        if np.isnan(v):
            return np.nan
        return high_label if v >= y0_px else low_label

    out = np.empty(arr.shape, dtype=object)
    nan_mask = np.isnan(arr)
    high_mask = (~nan_mask) & (arr >= y0_px)
    low_mask = (~nan_mask) & (arr < y0_px)
    out[nan_mask] = np.nan
    out[high_mask] = high_label
    out[low_mask] = low_label
    return out


def assign_zone_for_all_fish(
    positions_xy: np.ndarray,
    y0_px: float = Y_THRESHOLD_DEFAULT,
    flow_is_higher_y: bool = FLOW_IS_HIGHER_Y_DEFAULT,
) -> list[Any]:
    """Convenience: classify every fish's position at a single frame.

    Parameters
    ----------
    positions_xy
        Array of shape ``(n_fish, 2)`` with ``(x, y)`` for each fish at one
        frame. Only the second column (y) is used.
    y0_px, flow_is_higher_y
        Same as :func:`assign_zone`.

    Returns
    -------
    list
        One entry per fish: either ``"flow"``, ``"calm"``, or ``np.nan``.
    """
    arr = np.asarray(positions_xy, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(
            f"positions_xy must have shape (n_fish, 2); got {arr.shape}"
        )
    ys = arr[:, 1]
    labels = assign_zone(ys, y0_px=y0_px, flow_is_higher_y=flow_is_higher_y)
    return list(labels)
