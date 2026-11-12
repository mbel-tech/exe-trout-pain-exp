"""Figure 8 of the article: serotonergic system in the dorsomedial pallium (Dm): A 5-HT,
B 5-HIAA, C 5-HIAA/5-HT (written as Figure_10_Dm_serotonin.png, its working file name).
Usage: python 07_figures/figure_8_dm_serotonin.py <aggression grid csv> <mono_DM.csv> <outlier_report.md> <out dir>
       [<aggression contrasts_refit.csv>] [<mono/plasma cld dir>]
<aggression grid csv> and contrasts_refit.csv come from aggression_models_design_random_effects.R
(default: output/aggression_design_re/); mono_DM.csv and outlier_report.md from
plasma_monoamine_models.R. Letters (compact letter display) are drawn only when the cells of a
panel differ, from the Tukey CLD of the plasma/monoamine run (default: ../cld next to mono_DM.csv)."""
import sys, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cld_letters import letters_from_pvalues, letters_from_cld_csv, draw_letters, statement_p, anova_table_p, significant
grid, mono, outl, out = sys.argv[1:5]
contr = sys.argv[5] if len(sys.argv) > 5 else os.path.join(HERE, "..", "output", "aggression_design_re", "contrasts_refit.csv")
cld_dir = sys.argv[6] if len(sys.argv) > 6 else os.path.join(os.path.dirname(os.path.abspath(mono)), "..", "cld")
from figure_style import CELLS, apply_style, panel_letter, legend_rect, add_combined_legend
from figure_style import strip as style_strip
apply_style(11); rng = np.random.default_rng(1)
def strip(ax, vals_by_cell): return style_strip(ax, vals_by_cell, rng, tick_size=9.5)
def add_letters(ax, tops, letters):
    """letters: {cell: label} or None (None = cells do not differ, nothing is drawn)."""
    if letters: draw_letters(ax, range(4), tops, [letters.get(c, "") for c in CELLS])
# The aggression figure is drawn by make_figure_aggression_2x2.py (it used to be drawn here as three panels).
m = pd.read_csv(mono); m["group"] = m.Group.str.replace("Choice","choice"); 
drop = {}
for line in open(outl, encoding="utf-8"):
    p = [x.strip() for x in line.split("|")]
    if len(p) > 4 and p[1] == "Monoamine DM": drop.setdefault(p[2], set()).add(int(p[3]))
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3))
for a,(col,name,yl,ttl) in zip(ax, [("HT5","5-HT","5-HT (ng mg$^{-1}$ protein)","A"),("HIAA5","5-HIAA","5-HIAA (ng mg$^{-1}$ protein)","B"),("ratio_5HIAA_5HT","5-HIAA/5-HT","5-HIAA/5-HT","C")]):
    x = m[~m.Sample_ID.isin(drop.get(name, set()))]
    tops = strip(a, [x[(x.group==g)&(x.Treatment==t)][col] for g,t in CELLS]); a.set_ylabel(yl); panel_letter(a, ttl)
    txg = statement_p(os.path.join(cld_dir, "..", "anovas", "anova_stmt_TreatmentxGroup.csv"), f"DM_{col}", "key")   # letters only if treatment x group is significant
    add_letters(a, tops, letters_from_cld_csv(os.path.join(cld_dir, f"DM_{col}.csv"), order=CELLS) if significant(txg) else None)
fig.tight_layout(rect=legend_rect(fig)); add_combined_legend(fig)
fig.savefig(f"{out}/Figure_10_Dm_serotonin.png", dpi=300); plt.close(fig)
print({k:sorted(v) for k,v in drop.items()})
