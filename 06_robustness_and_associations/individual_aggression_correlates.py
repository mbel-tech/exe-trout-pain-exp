"""Body size and individual correlates of aggression (GOT experiment).
Usage: python 06_robustness_and_associations/individual_aggression_correlates.py <GOT datasets dir> <c23_behaviour_gated.csv> <mono_DM.csv>
       <outlier_report.md> <out dir>
Links aggression events (attacker/target identities) to sampled fish through the body-length
match table, then tests, with permutation WITHIN trial (so that differences between trials,
treatments and groups cannot generate the association):
  attacks given / received / net  vs  relative body mass, Dm 5-HIAA/5-HT, plasma cortisol,
  proportion of time in the flow zone.
Also tests trial-level aggression against mean body mass and body-mass variation."""
import sys, re, numpy as np, pandas as pd
from scipy import stats
data, c23f, monof, outl, out = sys.argv[1:6]
import os; os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(20261111)

ev = pd.read_csv(f"{data}/aggression_annotated_june2026.csv")
c23 = pd.read_csv(c23f)
c23["seg3_id"] = c23.tracking_fish_id.str.extract(r"_F(\d+)$").astype(int)
# segment-3 local id -> matched master id (from events), per trial
m3 = pd.concat([ev[ev.segment == 3][["trial", "attacker_id", "attacker_fish_id_matched_trial"]].set_axis(["trial", "lid", "mf"], axis=1),
                ev[ev.segment == 3][["trial", "attacked_id", "attacked_fish_id_matched_trial"]].set_axis(["trial", "lid", "mf"], axis=1)]).dropna().drop_duplicates()
amb = m3.groupby(["trial", "lid"]).mf.nunique().max()
fish = c23.merge(m3, left_on=["trial_id", "seg3_id"], right_on=["trial", "lid"], how="left")
print("matched fish:", len(fish), "| with master id:", fish.mf.notna().sum(), "| max ids per local id:", amb)

acc = ev[(ev.attacker_match_accepted == True) & (ev.attacked_match_accepted == True)]
given = acc.groupby("attacker_fish_id_matched_trial").size(); recv = acc.groupby("attacked_fish_id_matched_trial").size()
fish["given"] = fish.mf.map(given).fillna(0); fish["received"] = fish.mf.map(recv).fillna(0)
fish.loc[fish.mf.isna(), ["given", "received"]] = np.nan
fish["net"] = fish.given - fish.received

mono = pd.read_csv(monof); drop = set()
for line in open(outl, encoding="utf-8"):
    p = [x.strip() for x in line.split("|")]
    if len(p) > 4 and p[1] == "Monoamine DM" and p[2] == "5-HIAA/5-HT": drop.add(int(p[3]))
mono = mono[~mono.Sample_ID.isin(drop)].set_index("Sample_ID")
fish["Dm_ratio"] = fish.Sample_ID.map(mono.ratio_5HIAA_5HT)
cort_out = set()
for line in open(outl, encoding="utf-8"):
    p = [x.strip() for x in line.split("|")]
    if len(p) > 4 and p[1] == "Plasma" and "cortisol" in p[2].lower(): cort_out.add(int(p[3]))
fish["cortisol"] = fish.plasma_cortisol.where(~fish.Sample_ID.isin(cort_out))
fish["rel_weight"] = fish.weight_g - fish.groupby("trial_id").weight_g.transform("mean")
fish["weight_rank"] = fish.groupby("trial_id").weight_g.rank()
fish.to_csv(f"{out}/fish_level_linked.csv", index=False)

def within_trial(d, x, y, nperm=20000):
    d = d.dropna(subset=[x, y]); d = d[d.groupby("trial_id")[x].transform("count") >= 3]
    cx = d[x] - d.groupby("trial_id")[x].transform("mean"); cy = d[y] - d.groupby("trial_id")[y].transform("mean")
    # rank within trial for a Spearman-type statistic
    rx = d.groupby("trial_id")[x].rank(); ry = d.groupby("trial_id")[y].rank()
    rx = rx - rx.groupby(d.trial_id).transform("mean"); ry = ry - ry.groupby(d.trial_id).transform("mean")
    obs = np.corrcoef(rx, ry)[0, 1]; g = d.trial_id.values; ryv = ry.values.copy(); rxv = rx.values; c = 0
    idx = [np.where(g == t)[0] for t in np.unique(g)]
    for _ in range(nperm):
        p = ryv.copy()
        for i in idx: p[i] = rng.permutation(ryv[i])
        c += abs(np.corrcoef(rxv, p)[0, 1]) >= abs(obs) - 1e-12
    rho_raw, p_raw = stats.spearmanr(d[x], d[y])
    return dict(x=x, y=y, n_fish=len(d), n_trials=d.trial_id.nunique(), rho_within_trial=obs, p_perm_within_trial=(c + 1) / (nperm + 1), rho_pooled=rho_raw, p_pooled=p_raw)
rows = [within_trial(fish, x, y) for y in ["given", "received", "net"] for x in ["rel_weight", "Dm_ratio", "cortisol", "time_spent_zone_flow_prop"]]
r = pd.DataFrame(rows); r.to_csv(f"{out}/fish_level_associations.csv", index=False)
pd.set_option("display.width", 250); print(r.round(3).to_string(index=False))
# heaviest fish: share of attacks given
h = fish.dropna(subset=["given"]).copy(); tot = h.groupby("trial_id").given.transform("sum")
h["share"] = h.given / tot
big = h[h.weight_rank == h.groupby("trial_id").weight_rank.transform("max")]
print("share of attacks given by the heaviest matched fish per trial: mean %.2f (chance %.2f), n=%d trials; Wilcoxon vs chance p=%.3f" % (
    big.share.mean(), (1 / h.groupby("trial_id").given.transform("count")).mean(), len(big),
    stats.wilcoxon(big.share - 1 / h.groupby("trial_id").given.transform("count")[big.index]).pvalue))

# ---- trial level: aggression vs body mass ----
samp = pd.read_csv(f"{data}/GOT_monoamine_data_AUG2026.csv").drop_duplicates("sample_id")
samp["trial_key"] = samp.sampling_date.astype(str) + "_" + samp.treatment.str.lower()
key2trial = fish.merge(samp[["sample_id", "trial_key"]], left_on="Sample_ID", right_on="sample_id").groupby("trial_key").trial_id.agg(lambda s: s.mode()[0])
samp["trial_id"] = samp.trial_key.map(key2trial)
tw = samp.groupby("trial_id").weight_g.agg(mean_w="mean", sd_w="std", min_w="min", max_w="max", n_fish="count"); tw["cv_w"] = tw.sd_w / tw.mean_w
ne = ev.groupby(["trial", "group", "treatment_clean"]).size().rename("n").reset_index().set_index("trial")
t = ne.join(tw); t["logn"] = np.log(t.n); t.to_csv(f"{out}/trial_level_mass_aggression.csv")
print(t.round(2).to_string())
for v in ["mean_w", "cv_w"]:
    rho, p = stats.spearmanr(t[v], t.n); print(f"trial aggression vs {v}: Spearman rho={rho:.2f}, p={p:.3f}")
import statsmodels.formula.api as smf, statsmodels.api as sm
t["G"] = (t.group == "No choice") - 0.5; t["T"] = (t.treatment_clean == "Exercise") - 0.5
# untransformed counts: per-trial totals meet both assumptions (single Box-Cox rule, Supplementary Table S10)
for f in ["n ~ G*T", "n ~ G*T + mean_w", "n ~ G*T + cv_w"]:
    m = smf.ols(f, t).fit(); print(f, "| interaction F p=%.3f" % m.pvalues["G:T"], "| covariate p=%s" % (("%.3f" % m.pvalues.iloc[-1]) if "+" in f else "-"), "| df resid", int(m.df_resid))
print(samp.groupby(["treatment", "group"]).weight_g.agg(["mean", "sem", "count"]).round(1).to_string())
wm = smf.ols("weight_g ~ C(treatment)*C(group)", samp).fit(); print(sm.stats.anova_lm(wm, typ=2).round(3).to_string())
