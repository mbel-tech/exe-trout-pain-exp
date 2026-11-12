"""Shared visual conventions of the Python-drawn manuscript figures, so that every figure follows
the same rules as the R-drawn locomotion figure (figure_3_locomotion.R):

 - four combinations, always in the order Control/No Choice, Exercise/No Choice,
   Control/Exercise Choice, Exercise/Exercise Choice;
 - colour = combination (dark/light blue = No Choice, yellow/orange = Exercise Choice),
   shape = rearing treatment (triangle = Control, circle = Exercise);
 - points = observational units, black crossbar = mean, black error bar = SEM, both centred on
   the points of the cell;
 - panels are identified by a bold letter only (no panel titles);
 - one legend type everywhere, placed to the right of the panels (add_combined_legend): the
   "Treatment" block (two shapes). The "Group" block (two colour dots per row) is added only
   where the groups are not named on the x axis, i.e. in the interval plots (groups=True);
   in the strip plots every cell is named under the axis, so the Group block is left out;
 - bold axis titles; "No Choice" / "Exercise Choice" capitalised as in the text;
 - letters for pairwise differences are bold and come from cld_letters.py."""
import numpy as np
import matplotlib.pyplot as plt

COL = {("No choice", "Control"): "#0072B2", ("No choice", "Exercise"): "#56B4E9",
       ("Exercise choice", "Control"): "#E69F00", ("Exercise choice", "Exercise"): "#D55E00"}
MK = {"Control": "^", "Exercise": "o"}
CELLS = [("No choice", "Control"), ("No choice", "Exercise"),
         ("Exercise choice", "Control"), ("Exercise choice", "Exercise")]
LAB = ["Control\nNo Choice", "Exercise\nNo Choice", "Control\nExercise\nChoice", "Exercise\nExercise\nChoice"]
NAME = {c: l.replace("\n", " ") for c, l in zip(CELLS, LAB)}
LEGEND_W_IN = 1.55          # width reserved for the legend at the right of a figure (inches)


def apply_style(font_size=10):
    plt.rcParams.update({"font.size": font_size, "font.family": "Arial", "axes.labelweight": "bold",
                         "axes.spines.top": False, "axes.spines.right": False})


def panel_letter(ax, letter):
    """Bold panel letter at the top left; panels carry no title."""
    ax.set_title(letter, loc="left", fontweight="bold")


def strip(ax, vals_by_cell, rng, size=42, tick_size=9):
    """One cell per x position: points (jittered), black mean crossbar and SEM bar centred on them.
    Returns the top of each cell (highest point or upper SEM tip), for placing letters."""
    tops = []
    for i, (c, v) in enumerate(zip(CELLS, vals_by_cell)):
        v = np.asarray(v, float); v = v[~np.isnan(v)]
        ax.scatter(i + rng.uniform(-.12, .12, len(v)), v, s=size, marker=MK[c[1]], color=COL[c],
                   alpha=.85, zorder=3, linewidths=0)
        m, se = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
        ax.errorbar(i, m, yerr=se, color="black", capsize=6, lw=1.4, zorder=4)
        ax.plot([i - .2, i + .2], [m, m], color="black", lw=2.6, zorder=4)
        tops.append(max(v.max(), m + se))
    ax.set_xticks(range(4)); ax.set_xticklabels(LAB, fontsize=tick_size); ax.set_xlim(-.6, 3.6)
    return tops


def legend_rect(fig):
    """rect for fig.tight_layout() that leaves room for the legend on the right."""
    return (0, 0, 1 - LEGEND_W_IN / fig.get_figwidth(), 1)


def add_combined_legend(fig, y_center=0.5, groups=False):
    """Legend to the right of the panels, in the style of the locomotion figure: Treatment shapes,
    plus the Group colour dots when groups=True (needed only where the x axis does not name the
    groups). Call after fig.tight_layout(rect=legend_rect(fig))."""
    w, h = fig.get_figwidth(), fig.get_figheight()
    if not groups:
        height_in = 0.87
        ax = fig.add_axes([1 - LEGEND_W_IN / w, y_center - height_in / h / 2, LEGEND_W_IN / w, height_in / h])
        ax.axis("off"); ax.set_xlim(0, LEGEND_W_IN); ax.set_ylim(0, height_in)
        y = height_in - 0.14 - 0.29 * np.arange(3)
        kw = dict(va="center", ha="left", fontsize=9)
        ax.text(0.05, y[0], "Treatment", fontweight="bold", **kw)
        for yy, trt in [(y[1], "Control"), (y[2], "Exercise")]:
            ax.scatter([0.14], [yy], s=70, marker=MK[trt], color="black", linewidths=0)
            ax.text(0.30, yy, trt, **kw)
        return ax
    height_in = 1.75
    ax = fig.add_axes([1 - LEGEND_W_IN / w, y_center - height_in / h / 2, LEGEND_W_IN / w, height_in / h])
    ax.axis("off"); ax.set_xlim(0, LEGEND_W_IN); ax.set_ylim(0, height_in)
    y = height_in - 0.14 - 0.29 * np.arange(6)       # Group, No Choice, Exercise Choice, Treatment, Control, Exercise
    kw = dict(va="center", ha="left", fontsize=9)
    ax.text(0.05, y[0], "Group", fontweight="bold", **kw); ax.text(0.05, y[3], "Treatment", fontweight="bold", **kw)
    for yy, grp, lab in [(y[1], "No choice", "No Choice"), (y[2], "Exercise choice", "Exercise Choice")]:
        ax.scatter([0.14, 0.30], [yy, yy], s=70, color=[COL[(grp, "Control")], COL[(grp, "Exercise")]], linewidths=0)
        ax.text(0.44, yy, lab, **kw)
    for yy, trt in [(y[4], "Control"), (y[5], "Exercise")]:
        ax.scatter([0.14], [yy], s=70, marker=MK[trt], color="black", linewidths=0)
        ax.text(0.30, yy, trt, **kw)
    return ax
