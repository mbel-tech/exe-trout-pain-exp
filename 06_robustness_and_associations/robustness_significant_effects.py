"""Robustness triad for every significant effect (same logic as Bellio et al. 2026):
(1) exact permutation test on trial-level values, (2) exact Wilcoxon test on trial-level
values, (3) leave-one-trial-out. Usage: python 06_robustness_and_associations/robustness_significant_effects.py <behaviour run dir> <aggr grid csv> <mono_DM.csv> <outlier_report.md> <out csv>"""
import sys, itertools, numpy as np, pandas as pd
from scipy import stats
beh, grid, mono, outl, out = sys.argv[1:6]
rows = []
def perm2(x, y):               # exact two-sample permutation on difference of means
    z = np.r_[x, y]; n = len(x); obs = abs(np.mean(x) - np.mean(y)); c = t = 0
    for idx in itertools.combinations(range(len(z)), n):
        m = np.zeros(len(z), bool); m[list(idx)] = True
        c += abs(z[m].mean() - z[~m].mean()) >= obs - 1e-12; t += 1
    return c / t
def signflip(d):               # exact paired permutation
    d = np.asarray(d, float); obs = abs(d.mean()); c = t = 0
    for s in itertools.product([1, -1], repeat=len(d)):
        c += abs((d * s).mean()) >= obs - 1e-12; t += 1
    return c / t
def loo(vals, f):              # f(values without trial i) -> p
    ps = [f(vals.drop(i)) for i in vals.index]
    return min(ps), max(ps), sum(p < 0.05 for p in ps), len(ps)
def inter_perm(df, col):       # within-treatment (within-tank) permutation of group labels
    def stat(d): m = d.groupby(["Treatment", "Group"])[col].mean(); return (m["Exercise", "No choice"] - m["Exercise", "Exercise choice"]) - (m["Control", "No choice"] - m["Control", "Exercise choice"])
    obs = abs(stat(df)); c = t = 0
    blocks = [df[df.Treatment == tr] for tr in ["Control", "Exercise"]]
    def labelings(b):
        n = len(b); k = int((b.Group == "No choice").sum())
        for idx in itertools.combinations(range(n), k):
            g = np.array(["Exercise choice"] * n, dtype=object); g[list(idx)] = "No choice"; yield g
    for g0 in labelings(blocks[0]):
        for g1 in labelings(blocks[1]):
            d = pd.concat(blocks).copy(); d["Group"] = np.r_[g0, g1]
            c += abs(stat(d)) >= obs - 1e-12; t += 1
    return c / t
def inter_day(df, col):        # day-level: (Exercise - Control) per day, 3 vs 3 days
    w = df.pivot_table(index=["date", "Group"], columns="Treatment", values=col).reset_index()
    w["D"] = w["Exercise"] - w["Control"]
    a, b = w[w.Group == "No choice"].D.values, w[w.Group == "Exercise choice"].D.values
    return perm2(a, b), stats.mannwhitneyu(a, b, method="exact").pvalue

# ---- behaviour: interval and group x interval effects (within-trial contrasts) ----
ind = pd.read_csv(f"{beh}/datasets/ind_interval.csv"); gd = pd.read_csv(f"{beh}/datasets/gd_interval.csv")
def trial_diff(d, v):
    w = d.groupby(["trial_id", "Group", "Treatment", "timepoint"])[v].mean().unstack("timepoint")
    return (w[3] - w[2]).rename("d").reset_index().set_index("trial_id")
for d, v, eff in [(ind, "distance_m", "both"), (ind, "speed_cm_s_mean", "both"), (ind, "time_moving_TS1_prop", "both"),
                  (ind, "n_zone_transitions", "both"), (ind, "accel_cm_s2_mean", "interval"), (ind, "time_spent_zone_flow_prop", "interval"),
                  (gd, "mean_centroid_speed_cm_s", "interval"), (gd, "mean_hull_area_cm2", "interval"), (gd, "mean_polarization", "interval")]:
    td = trial_diff(d, v)
    f = lambda x: stats.wilcoxon(x.d.values, method="exact").pvalue
    rows.append(dict(outcome=v, effect="interval", n_trials=len(td), mean_change=td.d.mean(), perm_p=signflip(td.d.values), wilcoxon_p=f(td), loo=loo(td, f)))
    if eff == "both":
        a = lambda x: x[x.Group == "No choice"].d.values; b = lambda x: x[x.Group == "Exercise choice"].d.values
        f = lambda x: stats.mannwhitneyu(a(x), b(x), method="exact").pvalue
        rows.append(dict(outcome=v, effect="group x interval", n_trials=len(td), mean_change=a(td).mean() - b(td).mean(), perm_p=perm2(a(td), b(td)), wilcoxon_p=f(td), loo=loo(td, f)))
# ---- aggression: treatment x group on events per trial ----
g = pd.read_csv(grid); tr = g.groupby(["trial", "group", "treatment", "date"], as_index=False).n.sum()
# Single Box-Cox rule (2026-11-12): events per trial pass Shapiro-Wilk and Levene untransformed
# (analysis/boxcox_secondary), so the robustness tests use the raw counts, not log counts.
tr = tr.rename(columns={"group": "Group", "treatment": "Treatment"}).set_index("trial"); tr["n"] = tr.n.astype(float)
pd_, wd_ = inter_day(tr.reset_index(), "n")
rows.append(dict(outcome="aggressive events per trial", effect="treatment x group", n_trials=12, mean_change=np.nan, perm_p=inter_perm(tr, "n"),
                 wilcoxon_p=wd_, loo=loo(tr, lambda x: inter_perm(x, "n")), day_perm_p=pd_))
a = tr[(tr.Group == "No choice") & (tr.Treatment == "Exercise")].n.values; b = tr[(tr.Group == "No choice") & (tr.Treatment == "Control")].n.values
rows.append(dict(outcome="aggressive events per trial", effect="Exercise vs Control within No Choice", n_trials=6, mean_change=a.mean() - b.mean(), perm_p=perm2(a, b), wilcoxon_p=stats.mannwhitneyu(a, b, method="exact").pvalue, loo=None))
# ---- Dm 5-HIAA/5-HT: treatment x group on trial means ----
m = pd.read_csv(mono); drop = set()
for line in open(outl, encoding="utf-8"):
    p = [x.strip() for x in line.split("|")]
    if len(p) > 4 and p[1] == "Monoamine DM" and p[2] == "5-HIAA/5-HT": drop.add(int(p[3]))
m = m[~m.Sample_ID.isin(drop)].dropna(subset=["ratio_5HIAA_5HT"]); m["Group"] = m.Group.str.replace("Choice", "choice")
tm = m.groupby(["trial_uid", "date_uid", "Treatment", "Group"], as_index=False).ratio_5HIAA_5HT.mean().rename(columns={"date_uid": "date"}).set_index("trial_uid")
pd_, wd_ = inter_day(tm.reset_index(), "ratio_5HIAA_5HT")
rows.append(dict(outcome="Dm 5-HIAA/5-HT (trial means)", effect="treatment x group", n_trials=len(tm), mean_change=np.nan, perm_p=inter_perm(tm, "ratio_5HIAA_5HT"),
                 wilcoxon_p=wd_, loo=loo(tm, lambda x: inter_perm(x, "ratio_5HIAA_5HT")), day_perm_p=pd_))
r = pd.DataFrame(rows)
r[["loo_min_p", "loo_max_p", "loo_n_sig", "loo_n"]] = r.loo.apply(lambda x: pd.Series(x if x else [np.nan] * 4)); r = r.drop(columns="loo")
r.to_csv(out, index=False); pd.set_option("display.width", 250); print(r.round(4).to_string(index=False))
