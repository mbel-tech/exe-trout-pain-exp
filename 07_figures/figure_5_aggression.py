"""Figure 5 of the article: aggressive events, 2 x 2 panels (file Figure_7_aggression.png).
A all events per trial, B events in the calm zone, C events in the flow zone, D proportion of the
events of a trial that occurred in the flow zone. In every panel points are trials and the black
crossbar and error bar are the mean +- SEM of those points, centred on them (figure_style.strip).
Panels carry a bold letter only; the Group / Treatment legend is the one of the locomotion figure.

Letters appear in a panel only if the full interaction it displays is significant (cld_letters.py):
A treatment x group (count model on totals), B and C treatment x group x zone (count model by
zone), D treatment x group (binomial model R2b of the aggression run). Letters come from the
adjusted pairwise contrasts of the design random-effects refit
(output/aggression_design_re/contrasts_refit.csv, from aggression_models_design_random_effects.R).

Usage: python 07_figures/figure_5_aggression.py <aggression_models.R run dir> <out png> [<contrasts_refit.csv>]"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from figure_style import CELLS, apply_style, panel_letter, strip, legend_rect, add_combined_legend
from cld_letters import letters_from_pvalues, relabel_in_plot_order, draw_letters, anova_table_p, significant
run, out = sys.argv[1:3]
contr = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "..", "output", "aggression_design_re", "contrasts_refit.csv")
refit_anova = os.path.join(os.path.dirname(contr), "anova_refit.csv")
apply_style(11); rng = np.random.default_rng(1)


def agg_letters(model, side=None):
    """Letters from the refit contrasts, only if the panel's full interaction is significant."""
    term = "group:treatment:side" if side else "group:treatment"
    if not os.path.exists(contr) or not significant(anova_table_p(refit_anova, term, model)):
        return None
    c = pd.read_csv(contr); c = c[c.model == model]
    c = c[c.side == side] if side else c[c.side.isna()]
    name = {f"{g} {t}": (g, t) for g, t in CELLS}
    p = {tuple(name[x] for x in r.contrast.split(" / ")): r["p.value"] for _, r in c.iterrows()}
    return letters_from_pvalues(CELLS, p) if len(p) == 6 else None


def add_letters(ax, tops, letters):
    if letters:
        draw_letters(ax, range(4), tops, [letters.get(c, "") for c in CELLS])


d = pd.read_csv(f"{run}/datasets/counts_trial_segment_side.csv").rename(columns={"treatment_clean": "treatment"})
tr = d.groupby(["trial", "group", "treatment"], as_index=False).n_events.sum()
tz = d.groupby(["trial", "group", "treatment", "side"], as_index=False).n_events.sum()
assert tr.n_events.sum() == 442, tr.n_events.sum()
sel = lambda df, g, t: df[(df.group == g) & (df.treatment == t)].n_events

fig, axs = plt.subplots(2, 2, figsize=(11, 9)); ax = axs.ravel()
ax[0].set_ylim(-2, 85)
add_letters(ax[0], strip(ax[0], [sel(tr, g, t) for g, t in CELLS], rng), agg_letters("B_trial_x_interval_total"))
for k, side in enumerate(["calm", "flow"]):
    s = tz[tz.side == side]; a = ax[k + 1]; a.set_ylim(-2, 52)
    add_letters(a, strip(a, [sel(s, g, t) for g, t in CELLS], rng), agg_letters("A_trial_x_interval_x_zone", side))
ax[0].set_ylabel("Aggressive events per trial (20 min)")
ax[1].set_ylabel("Aggressive events in the calm zone\nper trial (20 min)")
ax[2].set_ylabel("Aggressive events in the flow zone\nper trial (20 min)")

# D: per-trial proportion of events in the flow zone, same crossbar mean +- SEM as A-C
p = pd.read_csv(f"{run}/datasets/probs_trial.csv").rename(columns={"treatment_clean": "treatment"})
p["prop"] = p.n_flow / p.n_total
a = ax[3]; a.set_ylim(0, 1.0)
tops = strip(a, [p[(p.group == g) & (p.treatment == t)].prop for g, t in CELLS], rng)
a.axhline(.5, ls="--", lw=.8, color="grey", zorder=1)
a.set_ylabel("Proportion of events in the flow zone")
txg_D = anova_table_p(f"{run}/anovas/ANOVA_R2b.csv", "group:treatment_clean")
if significant(txg_D):                       # letters of the binomial model, named in plot order
    e = pd.read_csv(f"{run}/cld/R2b_cld.csv")
    let = {(r.group, r.treatment_clean): str(r[".group"]).strip() for _, r in e.iterrows()}
    if len(set(let.values())) > 1:
        add_letters(a, tops, relabel_in_plot_order(let, CELLS))
for a_, L in zip(ax, "ABCD"):
    panel_letter(a_, L)
fig.tight_layout(rect=legend_rect(fig)); add_combined_legend(fig)
fig.savefig(out, dpi=300); plt.close(fig)
print("treatment x group p: all events", anova_table_p(refit_anova, "group:treatment", "B_trial_x_interval_total"),
      "| by zone (three-way)", anova_table_p(refit_anova, "group:treatment:side", "A_trial_x_interval_x_zone"),
      "| D", txg_D)
