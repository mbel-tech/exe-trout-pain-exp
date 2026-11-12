"""Figure 4 of the article: school-level indicators that changed between the two observation
intervals (A polarisation, B school area, C school speed); means +- SEM of the trials per
combination and interval.
Usage: python 07_figures/figure_4_school_structure.py <behaviour_models.R run dir> <out png>"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
beh, png = sys.argv[1:3]
COL = {("No choice", "Control"): "#0072B2", ("No choice", "Exercise"): "#56B4E9",
       ("Exercise choice", "Control"): "#E69F00", ("Exercise choice", "Exercise"): "#D55E00"}
MK = {"Control": "^", "Exercise": "o"}
CELLS = list(COL)
NAME = dict(zip(CELLS, ["Control No Choice", "Exercise No Choice", "Control Exercise Choice", "Exercise Exercise Choice"]))
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "font.family": "Arial"})
gd = pd.read_csv(f"{beh}/datasets/gd_interval.csv")
gd["Group"] = gd.Group.astype(str).str.strip().str.lower().map(lambda x: "No choice" if x.startswith("no") else "Exercise choice")
panels = [("mean_polarization", "Polarisation"), ("mean_hull_area_cm2", "School area (cm$^2$)"), ("mean_centroid_speed_cm_s", "School speed (cm s$^{-1}$)")]
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
for ax, (v, yl), L in zip(axes, panels, "ABC"):
    for j, c in enumerate(CELLS):
        d = gd[(gd.Group == c[0]) & (gd.Treatment == c[1])].groupby("timepoint")[v].agg(["mean", "sem"])
        ax.errorbar(np.array([0, 1]) + (j - 1.5) * 0.07, d["mean"].values, yerr=d["sem"].values, color=COL[c], marker=MK[c[1]], ms=8, lw=1.6, capsize=3, label=NAME[c])
    ax.set_xticks([0, 1]); ax.set_xticklabels(["20–30 min", "40–50 min"]); ax.set_xlim(-0.4, 1.4); ax.set_ylabel(yl); ax.set_title(L, loc="left", fontweight="bold")
h, l = axes[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, frameon=False)
fig.tight_layout(rect=(0, 0.09, 1, 1)); fig.savefig(png, dpi=300); plt.close(fig)
