"""Figures 4 (school-level indicators that changed between intervals; file
Figure_6_school_changed.png), 6 (aggression in time and space; file Figure_8_aggression_time_space.png)
and 7 (plasma cortisol and lactate; file Figure_9_plasma.png) of the article. The script also draws
two working figures that are not in the article (Figure_4_effect_summary.png, a summary of
standardised trial-level effects, and Figure_6_school.png, all six school indicators).
Usage: python 07_figures/figures_4_6_7_school_aggression_plasma.py <behaviour run dir>
       <plasma/monoamine run dir> <aggr grid csv> <aggression events csv> <out dir>
Letters (compact letter display) are drawn in the school and plasma figures only for panels whose cells differ in
the pipeline Tukey CLD (<run dir>/cld/); panels in which all cells share a letter carry none."""
import sys, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cld_letters import letters_from_cld_csv, draw_letters, statement_p, significant
beh, endo, gridf, evf, out = sys.argv[1:6]
from figure_style import COL, MK, CELLS, LAB, NAME, apply_style, panel_letter, legend_rect, add_combined_legend
from figure_style import strip as style_strip
apply_style(10)
rng = np.random.default_rng(1)
canon = lambda s: s.astype(str).str.strip().str.lower().map(lambda x: "No choice" if x.startswith("no") else "Exercise choice")

# ------------------------------------------------------------------ trial-level table
ind = pd.read_csv(f"{beh}/datasets/ind_interval.csv"); gd = pd.read_csv(f"{beh}/datasets/gd_interval.csv")
for d in (ind, gd): d["Group"] = canon(d.Group)
drop = {}
for line in open(f"{endo}/reports/outlier_report.md", encoding="utf-8"):
    p = [x.strip() for x in line.split("|")]
    if len(p) > 4 and p[3].isdigit(): drop.setdefault((p[1], p[2]), set()).add(int(p[3]))
def endo_trial(fname, module, col, label, idcol):
    d = pd.read_csv(f"{endo}/datasets/{fname}"); d["Group"] = canon(d.Group)
    d = d[~d[idcol].isin(drop.get((module, label), set()))]
    return d.groupby(["trial_uid", "Group", "Treatment"])[col].mean().reset_index().rename(columns={col: "v"})
def beh_trial(d, v):
    return d.groupby(["trial_id", "Group", "Treatment"])[v].mean().reset_index().rename(columns={v: "v"})
def beh_change(d, v):
    w = d.groupby(["trial_id", "Group", "Treatment", "timepoint"])[v].mean().unstack("timepoint")
    return (w[3] - w[2]).rename("v").reset_index()
def effects(t):
    cm = t.groupby(["Group", "Treatment"]).v.agg(["mean", "var", "count"])
    s = np.sqrt(((cm["count"] - 1) * cm["var"]).sum() / (cm["count"] - 1).sum()); df = int((cm["count"] - 1).sum())
    m = cm["mean"]; n = cm["count"]; tq = stats.t.ppf(0.975, df)
    def con(w):
        est = sum(w[c] * m[c] for c in w) / s; se = np.sqrt(sum(w[c] ** 2 / n[c] for c in w)); return est, est - tq * se, est + tq * se
    NC, EC, C, E = "No choice", "Exercise choice", "Control", "Exercise"
    return dict(treatment=con({(NC, E): .5, (EC, E): .5, (NC, C): -.5, (EC, C): -.5}),
                group=con({(NC, E): .5, (NC, C): .5, (EC, E): -.5, (EC, C): -.5}),
                interaction=con({(NC, E): .5, (NC, C): -.5, (EC, E): -.5, (EC, C): .5}))
rowsF = []
def add(family, label, t, change=None):
    e = effects(t); r = dict(family=family, label=label, **{k: v for k, v in e.items()})
    if change is not None: r["gxi"] = effects(change)["group"]
    rowsF.append(r)
for v, l in [("distance_m", "Distance travelled"), ("speed_cm_s_mean", "Swimming speed"), ("accel_cm_s2_mean", "Acceleration"),
             ("time_moving_TS1_prop", "Time swimming > 1 BL s$^{-1}$"), ("n_zone_transitions", "Zone transitions"), ("time_spent_zone_flow_prop", "Time in flow zone")]:
    add("Locomotion", l, beh_trial(ind, v), beh_change(ind, v))
for v, l in [("mean_iid_cm", "Inter-individual distance"), ("mean_nnd_cm", "Nearest-neighbour distance"), ("mean_polarization", "Polarisation"),
             ("mean_hull_area_cm2", "School area"), ("mean_centroid_speed_cm_s", "School speed"), ("zone_synchrony", "Zone synchrony")]:
    add("School", l, beh_trial(gd, v), beh_change(gd, v))
g = pd.read_csv(gridf); g["Group"] = canon(g.group); g["Treatment"] = g.treatment
ta = g.groupby(["trial", "Group", "Treatment"]).n.sum().reset_index(); ta["v"] = ta.n.astype(float)   # untransformed (single Box-Cox rule: per-trial counts meet both assumptions)
w = g.groupby(["trial", "Group", "Treatment", "segment"]).n.sum().unstack("segment"); ca = (w[3] - w[2]).rename("v").reset_index()
add("Aggression", "Aggressive events", ta, ca)
pl = pd.read_csv(f"{endo}/datasets/plasma.csv"); idc = "Sample" if "Sample" in pl.columns else [c for c in pl.columns if c.lower().startswith("sample")][0]
add("Plasma", "Cortisol", endo_trial("plasma.csv", "Plasma", "plasma_cortisol", "Plasma cortisol", idc))
add("Plasma", "Lactate", endo_trial("plasma.csv", "Plasma", "Lactate_mmol_L", "Lactate (mmol/L)", idc))
for reg, f in [("DM", "mono_DM.csv"), ("RN", "mono_RN.csv")]:
    for col, l in [("HT5", "5-HT"), ("HIAA5", "5-HIAA"), ("ratio_5HIAA_5HT", "5-HIAA/5-HT"), ("DA", "DA"), ("DOPAC", "DOPAC"), ("ratio_DOPAC_DA", "DOPAC/DA")]:
        add("Dm" if reg == "DM" else "Raphe nuclei", l, endo_trial(f, f"Monoamine {reg}", col, l, "Sample_ID"))
# status of each effect in the main analysis: R = significant and robust, S = significant, not robust
STATUS = {("Aggressive events", "interaction"): "R", ("Zone transitions", "gxi"): "R",
          ("Distance travelled", "gxi"): "S", ("Swimming speed", "gxi"): "S", ("Time swimming > 1 BL s$^{-1}$", "gxi"): "S",
          ("Polarisation", "group"): "S", ("5-HIAA/5-HT|Dm", "interaction"): "S"}
pd.DataFrame([{**{"family": r["family"], "label": r["label"]}, **{f"{k}_{s}": r[k][i] for k in ["treatment", "group", "interaction", "gxi"] if k in r for i, s in enumerate(["d", "lo", "hi"])}} for r in rowsF]).to_csv(f"{out}/standardised_effects_trial_level.csv", index=False)

fig, axes = plt.subplots(1, 4, figsize=(13.5, 9.2), sharey=True)
titles = [("treatment", "Rearing treatment\n(Exercise − Control)"), ("group", "Arena condition\n(No Choice − Exercise Choice)"),
          ("interaction", "Treatment × group\n(half difference of differences)"), ("gxi", "Group × interval\n(change, No Choice − Exercise Choice)")]
ys, y, ylab, fam_pos = [], 0, [], []
last = None
for r in rowsF:
    if r["family"] != last: y += 0.8; fam_pos.append((y - 0.75, r["family"])); last = r["family"]
    ys.append(y); ylab.append(r["label"]); y += 1
for ax, (k, ttl) in zip(axes, titles):
    ax.axvline(0, color="0.6", lw=0.8)
    for r, yy in zip(rowsF, ys):
        if k not in r: continue
        d, lo, hi = r[k]; key = r["label"] + ("|Dm" if r["family"] == "Dm" and r["label"] == "5-HIAA/5-HT" else "")
        st = STATUS.get((key, k)) or STATUS.get((r["label"], k) if r["family"] not in ("Dm", "Raphe nuclei") else None)
        c = "#C0392B" if st == "R" else ("#C0392B" if st == "S" else "0.35")
        ax.plot([lo, hi], [yy, yy], color=c, lw=1.4)
        ax.plot(d, yy, "o", ms=7, mfc=(c if st == "R" else "white" if st == "S" else c), mec=c, mew=1.5)
    ax.set_title(ttl, fontsize=10); ax.set_xlabel("Standardised effect (trial-level SD)"); ax.set_xlim(-4.2, 4.2)
axes[0].set_yticks(ys); axes[0].set_yticklabels(ylab); axes[0].invert_yaxis()
for yy, f in fam_pos: axes[0].text(-0.02, yy, f, transform=axes[0].get_yaxis_transform(), ha="right", va="center", fontweight="bold", fontsize=10)
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([], [], marker="o", color="#C0392B", mfc="#C0392B", lw=1.4, label="significant and robust"),
                    Line2D([], [], marker="o", color="#C0392B", mfc="white", lw=1.4, label="significant, not robust"),
                    Line2D([], [], marker="o", color="0.35", mfc="0.35", lw=1.4, label="not significant")], loc="lower center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0.035, 1, 1)); fig.savefig(f"{out}/Figure_4_effect_summary.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Figure 6: school structure
SCHOOL = {"mean_polarization": "Polarisation", "mean_nnd_cm": "Nearest-neighbour distance (cm)", "mean_iid_cm": "Inter-individual distance (cm)",
          "mean_hull_area_cm2": "School area (cm$^2$)", "mean_centroid_speed_cm_s": "School speed (cm s$^{-1}$)", "zone_synchrony": "Zone synchrony"}
def school_figure(variables, nrows, ncols, figsize, fname):
    """Interval plot (means +- SEM per combination and interval) of school-level indicators."""
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    for ax, v, L in zip(np.ravel(axes), variables, "ABCDEF"):
        cld = letters_from_cld_csv(f"{beh}/cld/interval_GD_{v}.csv", by=("Group", "Treatment", "interval"), order=[(g_, t_, tp) for tp in (2, 3) for g_, t_ in CELLS]); lx, lt, ll = [], [], []
        if not significant(statement_p(f"{beh}/anovas/anova_stmt_interval_TreatmentxGroupxInterval.csv", f"GD_{v}")): cld = None   # letters only if treatment x group x interval is significant
        for j, c in enumerate(CELLS):
            d = gd[(gd.Group == c[0]) & (gd.Treatment == c[1])].groupby("timepoint")[v].agg(["mean", "sem"])
            x = np.array([0, 1]) + (j - 1.5) * 0.15   # series spread apart as in the locomotion figure
            ax.errorbar(x, d["mean"].values, yerr=d["sem"].values, color=COL[c], marker=MK[c[1]], ms=8, lw=1.6, capsize=3)
            if cld: lx += list(x); lt += list((d["mean"] + d["sem"]).values); ll += [cld.get((c[0], c[1], int(tp)), "") for tp in d.index]
        ax.set_xticks([0, 1]); ax.set_xticklabels(["20–30 min", "40–50 min"], fontweight="bold"); ax.set_xlim(-0.45, 1.45); ax.set_ylabel(SCHOOL[v]); panel_letter(ax, L)
        if cld: draw_letters(ax, lx, lt, ll, frac=0.015, fontsize=9.5)
    fig.tight_layout(rect=legend_rect(fig)); add_combined_legend(fig, groups=True)   # groups are not named on the x axis here, so the Group block stays
    fig.savefig(f"{out}/{fname}", dpi=300); plt.close(fig)
school_figure(list(SCHOOL), 2, 3, (12.5, 7), "Figure_6_school.png")
# the three indicators that changed between intervals (figure of the Methods-plus-Results draft)
school_figure(["mean_polarization", "mean_hull_area_cm2", "mean_centroid_speed_cm_s"], 1, 3, (12.5, 3.9), "Figure_6_school_changed.png")

# ------------------------------------------------------------------ Figure 9: plasma
def strip(ax, vals, letters=None):
    tops = style_strip(ax, vals, rng, size=40)
    if letters: draw_letters(ax, range(4), tops, [letters.get(c, "") for c in CELLS])
pl["Group"] = canon(pl.Group)
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))
for ax, (col, lab, yl, L) in zip(axes, [("plasma_cortisol", "Plasma cortisol", "Plasma cortisol (ng mL$^{-1}$)", "A"), ("Lactate_mmol_L", "Lactate (mmol/L)", "Plasma lactate (mmol L$^{-1}$)", "B")]):
    d = pl[~pl[idc].isin(drop.get(("Plasma", lab), set()))]
    strip(ax, [d[(d.Group == g_) & (d.Treatment == t_)][col] for g_, t_ in CELLS], letters_from_cld_csv(f"{endo}/cld/{col}.csv", order=CELLS) if significant(statement_p(f"{endo}/anovas/anova_stmt_TreatmentxGroup.csv", col, "key")) else None); ax.set_ylabel(yl); panel_letter(ax, L)
fig.tight_layout(rect=legend_rect(fig)); add_combined_legend(fig)
fig.savefig(f"{out}/Figure_9_plasma.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Figure 8: aggression in time and space
ev = pd.read_csv(evf); ev["Group"] = canon(ev.group); ev["Treatment"] = ev.treatment_clean
order = ev.drop_duplicates("trial").assign(k=lambda d: d.apply(lambda r: CELLS.index((r.Group, r.Treatment)), axis=1)).sort_values(["k", "trial"])
fig = plt.figure(figsize=(13.5, 6.2)); gs = fig.add_gridspec(2, 4, width_ratios=[2.2, 2.2, 1.5, 1.5])
ax = fig.add_subplot(gs[:, :2])
for i, (_, r) in enumerate(order.iterrows()):
    e = ev[ev.trial == r.trial]; c = COL[(r.Group, r.Treatment)]
    x = e.time_within_segment_sec + np.where(e.segment == 3, 660, 0)
    ax.vlines(x, i - 0.38, i + 0.38, color=c, lw=0.9)
ax.axvspan(600, 660, color="0.92", zorder=0); ax.set_xlim(0, 1260); ax.set_ylim(-0.7, 11.7); ax.invert_yaxis()
ax.set_xticks([0, 300, 600, 660, 960, 1260]); ax.set_xticklabels(["20", "25", "30", "40", "45", "50"]); ax.set_xlabel("Time after transfer to the arena (min)")
ax.set_yticks(range(12)); ax.set_yticklabels([f"{NAME[(r.Group, r.Treatment)]}  (trial {r.trial})" for _, r in order.iterrows()], fontsize=8.5)
ax.set_title("A", loc="left", fontweight="bold")
ev["ax_x"] = [r.get(f"salmon_{int(r.attacker_id)}_x", np.nan) if pd.notna(r.attacker_id) else np.nan for _, r in ev.iterrows()]
ev["ax_y"] = [r.get(f"salmon_{int(r.attacker_id)}_y", np.nan) if pd.notna(r.attacker_id) else np.nan for _, r in ev.iterrows()]
xl = (np.nanmin(ev.ax_x), np.nanmax(ev.ax_x)); yl = (np.nanmin(ev.ax_y), np.nanmax(ev.ax_y))
bound = np.nanmean([ev[ev.attacker_side_at_event == "calm"].ax_y.max(), ev[ev.attacker_side_at_event == "flow"].ax_y.min()])
for j, c in enumerate(CELLS):
    a = fig.add_subplot(gs[j // 2, 2 + j % 2]); e = ev[(ev.Group == c[0]) & (ev.Treatment == c[1])]
    a.scatter(e.ax_x, e.ax_y, s=9, color=COL[c], alpha=0.6, linewidths=0)
    a.axhline(bound, color="0.3", lw=0.8, ls="--"); a.set_xlim(*xl); a.set_ylim(*yl); a.invert_yaxis(); a.set_xticks([]); a.set_yticks([])
    for sp in a.spines.values(): sp.set_visible(True)
    a.set_title(("B   " if j == 0 else "") + NAME[c], fontsize=9, loc="left", fontweight="bold" if j == 0 else "normal")
    box = dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.0)   # keeps the zone labels readable over the points
    a.text(0.02, 0.04, "flow zone", transform=a.transAxes, fontsize=7.5, color="0.3", bbox=box, zorder=5); a.text(0.02, 0.92, "calm zone", transform=a.transAxes, fontsize=7.5, color="0.3", bbox=box, zorder=5)
fig.tight_layout(); fig.savefig(f"{out}/Figure_8_aggression_time_space.png", dpi=300); plt.close(fig)
print("boundary y (px):", round(bound, 1), "| events with position:", ev.ax_x.notna().sum())
