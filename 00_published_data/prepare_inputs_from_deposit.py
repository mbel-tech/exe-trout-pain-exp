"""Rebuild the input files of the statistics and figure scripts from the published data.

The data behind the article are deposited on Zenodo (doi:10.5281/zenodo.24136844) as tidy
tables (got_*.csv, lower-case labels, empty cells for missing values, a data dictionary).
The scripts in 05_statistics/ to 07_figures/ were written against the working files of the
experiment. This script converts the deposit tables into those working files, with the same
file names, column names and labels, so that the scripts run unchanged.

Usage (from the repository root):
    python 00_published_data/prepare_inputs_from_deposit.py <unzipped deposit folder>

It writes:
    data/aggression_annotated_june2026.csv          aggression events (2.8.3)
    data/GOT_monoamine_data_AUG2026.csv             one row per fish x brain region (Dm, Rn)
    data/GOT glu_lac_cort updated march 2025.xlsx   plasma and body measures, sheet "results"
    data/monoamine concentration ALL GOT.xlsx       monoamines, sheet "sheet1"
    02_trajectory_pipeline/STEP6_output/from_deposit/01_data/
        session_gd_summary_TP2_TP3_ExeChoice_NoChoice.csv    school-level behaviour
        NND_fish_by_timepoint.csv                            nearest-neighbour distance per fish
        individual_analysis_dataset_TP2_TP3_ExeChoice_NoChoice.csv   fish-level behaviour
    02_trajectory_pipeline/STEP1_5/STEP1_5_deposit/
        idmatcher_id_map_deposit.csv, reliability_tier_deposit.csv   identity across intervals
    02_trajectory_pipeline/STEP1_6_output/STEP1_6_output_20260819_181051/
        match_table_export_20260819_181445.xlsx      tracked fish -> sampled fish (body length)

The deposit holds no raw trajectories, video or BORIS files, so steps 01-04 (tracking,
trajectory pipeline, aggression-event extraction, sample-dataset build) cannot be run from
it; their outputs are what the deposit contains. The Dl, Vv and POA monoamines are not in the
deposit, so plasma_monoamine_models.R analyses Dm and Rn only when run from these files.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if len(sys.argv) < 2:
    sys.exit(__doc__)
DEP = Path(sys.argv[1])
ROOT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
STEP6 = ROOT / "02_trajectory_pipeline" / "STEP6_output" / "from_deposit" / "01_data"
STEP15 = ROOT / "02_trajectory_pipeline" / "STEP1_5" / "STEP1_5_deposit"
STEP16 = ROOT / "02_trajectory_pipeline" / "STEP1_6_output" / "STEP1_6_output_20260819_181051"
for d in (DATA, STEP6, STEP15, STEP16):
    d.mkdir(parents=True, exist_ok=True)


def read(name):
    return pd.read_csv(DEP / f"{name}.csv")


def cap(s):
    """'exercise choice' -> 'Exercise choice', 'control' -> 'Control' (labels of the working files)."""
    return s.str[0].str.upper() + s.str[1:]


def title(s):
    """'exercise choice' -> 'Exercise Choice' (labels of the laboratory workbooks)."""
    return s.str.title()


def write_csv(df, path):
    df.to_csv(path, index=False)
    print(f"  {path.relative_to(ROOT)}  ({len(df)} rows)")


trials = read("got_trials")
fish = read("got_fish")
plasma = read("got_plasma")
mono_long = read("got_monoamines_long")
mono_samples = read("got_monoamine_samples")
school = read("got_behaviour_school_interval")
ind = read("got_behaviour_fish_interval")
cross = read("got_tracking_to_sample_crosswalk")
events = read("got_aggression_events")

print("Writing input files:")

# ---- behaviour (05_statistics/behaviour_models.R) ---------------------------------------
lab = lambda d: d.assign(Group=cap(d.group), Treatment=cap(d.treatment)).drop(columns=["group", "treatment"])

gd = lab(school)[["trial_id", "timepoint", "Group", "Treatment", "trial_date", "mean_polarization",
                  "mean_iid_cm", "zone_synchrony", "mean_centroid_speed_cm_s", "mean_hull_area_cm2"]]
write_csv(gd, STEP6 / "session_gd_summary_TP2_TP3_ExeChoice_NoChoice.csv")

# per-fish NND, keyed by the fish number within the recording (T5_S2_F3 -> 3)
nnd = lab(ind)[["trial_id", "timepoint", "tracking_fish_id", "mean_nnd_cm", "Group", "Treatment", "trial_date"]]
nnd = nnd.assign(fish_id=nnd.tracking_fish_id.str.extract(r"_F(\d+)$", expand=False).astype(int))
nnd = nnd[nnd.mean_nnd_cm.notna()][["trial_id", "timepoint", "fish_id", "mean_nnd_cm", "Group", "Treatment", "trial_date"]]
write_csv(nnd, STEP6 / "NND_fish_by_timepoint.csv")

# fish-level indicators; identity, NND and proportions are re-derived by behaviour_models.R
derived = ["mean_nnd_cm", "master_fish_id", "id_match_quality", "id_reliability_tier",
           "time_spent_zone_flow_prop", "time_moving_TS1_prop"]
ind_out = lab(ind).rename(columns={"tracking_fish_id": "fish_id"})
ind_out = ind_out.drop(columns=derived + [c for c in ind_out if c.endswith("_missing_reason")])
first = ["fish_id", "trial_id", "timepoint", "Group", "Treatment", "trial_date"]
ind_out = ind_out[first + [c for c in ind_out if c not in first]]
write_csv(ind_out, STEP6 / "individual_analysis_dataset_TP2_TP3_ExeChoice_NoChoice.csv")

# fish the matcher never scored (no match quality) were not in the identity map
idmap = ind[ind.id_match_quality.notna()][["trial_id", "timepoint", "tracking_fish_id", "master_fish_id", "id_match_quality"]].rename(columns={
    "timepoint": "segment", "tracking_fish_id": "tracker_fish_id",
    "master_fish_id": "fish_id_matched_trial", "id_match_quality": "match_quality"})
write_csv(idmap, STEP15 / "idmatcher_id_map_deposit.csv")

tier = ind[["trial_id", "id_reliability_tier"]].drop_duplicates().rename(columns={"id_reliability_tier": "idmatcher_tier"})
assert tier.trial_id.is_unique, "one identity-reliability tier per trial expected"
write_csv(tier, STEP15 / "reliability_tier_deposit.csv")

# ---- aggression events (05_statistics/aggression_*.R, 06_, 07_) -------------------------
ev = events.rename(columns={"trial_id": "trial", "timepoint": "segment", "trial_date": "date"})
ev["treatment_clean"] = cap(ev.treatment)
ev["treatment"] = ev.boris_treatment_label          # treatment label as typed in BORIS
ev["group"] = cap(ev.group)
ev = ev.rename(columns={c: c.replace("fish_", "salmon_", 1) for c in ev if c.startswith("fish_") and c[5].isdigit()})
ev = ev.drop(columns=["event_id", "boris_treatment_label"] + [c for c in ev if c.endswith("_missing_reason")])
first = ["trial", "segment", "group", "treatment_clean", "date", "treatment"]
ev = ev[first + [c for c in ev if c not in first]]
write_csv(ev, DATA / "aggression_annotated_june2026.csv")

# ---- monoamines, one row per fish x region ----------------------------------------------
ANALYTE = {"5-HT": "HT5", "5-HIAA": "HIAA5", "5-HIAA/5-HT": "HIAA5_HT5", "DA": "DA", "DOPAC": "DOPAC", "DOPAC/DA": "DOPAC_DA"}
mw = (mono_long.assign(analyte=mono_long.analyte.map(ANALYTE))
      .pivot(index=["sample_id", "region"], columns="analyte", values="value")  # keeps all-missing rows
      .reindex(columns=list(ANALYTE.values())).reset_index())
mw = (mw.merge(mono_samples[["sample_id", "region", "operator", "punches", "protein_ug", "tissue_status"]],
               on=["sample_id", "region"], how="left")
        .merge(fish[["sample_id", "trial_id", "group", "treatment", "replicate", "sampling_date", "sex",
                     "weight_g", "length_cm"]], on="sample_id", how="left"))
mono_csv = mw.rename(columns={"region": "brain_region", "replicate": "trial"})
mono_csv = mono_csv[["sample_id", "brain_region", "operator", "punches", "protein_ug", "HT5", "HIAA5", "HIAA5_HT5",
                     "DA", "DOPAC", "DOPAC_DA", "tissue_status", "group", "treatment", "trial", "trial_id",
                     "sampling_date", "sex", "weight_g", "length_cm"]].sort_values(["sample_id", "brain_region"])
write_csv(mono_csv, DATA / "GOT_monoamine_data_AUG2026.csv")

# ---- laboratory-workbook layouts (plasma_monoamine_models.R, link_behaviour_physiology.R) ---
pf = plasma.merge(fish, on=["sample_id", "trial_id", "group", "treatment"], how="left")
hct = pf[["hematocrit_1_pct", "hematocrit_2_pct"]].mean(axis=1)
results = pd.DataFrame({
    "Sample": pf.sample_id, "Group": title(pf.group), "tank": cap(pf.treatment),
    "plasma_cortisol": pf.cortisol, "Glucose_mmol_L": pf.glucose_mmol_L, "glu_plate": pf.glucose_plate,
    "Lactate_mmol_L": pf.lactate_mmol_L, "lac_plate": pf.lactate_plate, "plate": pf.cortisol_plate,
    "plate_date": pd.to_datetime(pf.cortisol_plate_date), "sampling_date": pd.to_datetime(pf.sampling_date),
    "weight": pf.weight_g, "length": pf.length_cm, "hematocrit avg.": hct, "haemoglobin": pf.haemoglobin,
    "sex": pf.sex,
})
# The original sheet repeats the "Sample" column further right; link_behaviour_physiology.R
# relies on the resulting name of the first one (sample_1 after readxl + janitor).
results.insert(len(results.columns), "Sample ", pf.sample_id)
results = results.rename(columns={"Sample ": "Sample"})
with pd.ExcelWriter(DATA / "GOT glu_lac_cort updated march 2025.xlsx") as xw:
    results.to_excel(xw, sheet_name="results", index=False)
print(f"  data/GOT glu_lac_cort updated march 2025.xlsx  ({len(results)} rows)")

sheet1 = pd.DataFrame({
    "Sample_ID": mono_csv.sample_id, "Brain_area": mono_csv.brain_region, "user": mono_csv.operator,
    "Group": title(mono_csv.group), "Treatment": cap(mono_csv.treatment), "Sampling_date": pd.to_datetime(mono_csv.sampling_date).dt.strftime("%m/%d/%Y"),
    "Trial": mono_csv.trial, "5-HT": mono_csv.HT5, "5-HIAA": mono_csv.HIAA5, "5-HIAA/5-HT": mono_csv.HIAA5_HT5,
    "DA": mono_csv.DA, "DOPAC": mono_csv.DOPAC, "DOPAC/DA": mono_csv.DOPAC_DA,
})
with pd.ExcelWriter(DATA / "monoamine concentration ALL GOT.xlsx") as xw:
    sheet1.to_excel(xw, sheet_name="sheet1", index=False)
print(f"  data/monoamine concentration ALL GOT.xlsx  ({len(sheet1)} rows)")

# ---- tracked fish -> sampled fish (body-length match, 2.4.3) ----------------------------
# Reliability tier as in 02_trajectory_pipeline/01_6_match_body_length.R: relative margin =
# ambiguity margin / smallest difference between the sampled lengths of the trial;
# < 0.4 unreliable, < 0.8 borderline, otherwise reliable, and conservative when the
# tracked-sampled length difference is also <= 0.5 cm.
spacing = fish.groupby("trial_id").length_cm.apply(lambda s: np.diff(np.sort(s.values)).min())
rel = cross.ambiguity_margin_cm / cross.trial_id.map(spacing)
tier_bl = np.select([rel < 0.4, rel < 0.8, cross.abs_length_difference_cm <= 0.5],
                    ["unreliable", "borderline", "conservative"], "reliable")
match = pd.DataFrame({
    "trial_id": cross.trial_id, "timepoint": cross.timepoint,
    "sampling FISH ID": cross.sample_id, "tracking FISH ID": cross.tracking_fish_id,
    "sampling lenght (cm)": cross.sampling_length_cm, "tracking length (cm)": cross.tracking_length_cm,
    "absolute length difference (cm)": cross.abs_length_difference_cm,
    "ambiguity margin (cm)": cross.ambiguity_margin_cm, "minimum spacing (cm)": cross.trial_id.map(spacing),
    "relative margin": rel, "reliability_tier": tier_bl,
    "IRS identity recovery probability (%)": cross.irs_recovery_probability_pct,
    "IRS reliability bin": cross.irs_reliability_bin,
})
with pd.ExcelWriter(STEP16 / "match_table_export_20260819_181445.xlsx") as xw:
    match.to_excel(xw, index=False)
print(f"  {(STEP16 / 'match_table_export_20260819_181445.xlsx').relative_to(ROOT)}  ({len(match)} rows)")

print("""
Then run, from the repository root (bash syntax):
  GOT_STEP6_DIR=02_trajectory_pipeline/STEP6_output/from_deposit/01_data \\
  GOT_STEP1_5_DIR=02_trajectory_pipeline/STEP1_5/STEP1_5_deposit GOT_STEP1_5_STAMP=deposit \\
    Rscript 05_statistics/behaviour_models.R
  Rscript 05_statistics/plasma_monoamine_models.R
and the remaining steps 5-7 as listed in README.md.""")
