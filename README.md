# Analysis code: exercise, exercise choice and coping with a nociceptive procedure in juvenile rainbow trout (GOT experiment)

This repository holds the code used to analyse the data of the GOT experiment. Juvenile rainbow trout
(*Oncorhynchus mykiss*) were reared for 40 days with (Exercise) or without (Control) water flow. They then
received an acetic-acid injection in the lips and spent 1 h in groups of four in an arena that either offered
a flow zone (Exercise Choice) or did not (No Choice). That gives a 2 × 2 design with 12 trials and 48 fish.

The code is organised in the order the article describes the analysis. Section numbers below
(2.8.1, 2.9.2, …) refer to the Materials and Methods of the article. "Supplementary Material" refers to its
supplementary file.

> **Article:** *reference and DOI to be added on publication.*
> **Data:** not included in this repository. *Location of the data deposit to be added.*

## Contents

| Folder | Article section | What it does |
|---|---|---|
| `01_tracking/` | 2.8.1; Supp. "Automated identity tracking workflow", "Identity matching across intervals", "Body-length estimation" | Helpers around idtracker.ai: bounding-box extraction used for body-length estimation, and the Google Colab notebook that runs idmatcher.ai on the two intervals of each trial |
| `02_trajectory_pipeline/` | 2.8.1, 2.8.2; Supp. "Trajectory processing in R", "Identity matching across intervals" | R pipeline from idtracker.ai trajectories to individual and school-level indicators: preprocessing (jump rejection, gap interpolation ≤ 1 s, cm conversion, body length), identity-match annotation (Rules 1–3, Supp. Table SM1), body-length match to the sampled fish, individual kinematic indicators, school-level indicators (swaRm/swaRmverse) |
| `03_aggression_events/` | 2.8.3; Supp. "Quantification of aggression events from BORIS annotations and idtracker.ai trajectories" | Python workflow that assigns each BORIS-scored aggression event to an attacker and a target (peak acceleration within ±1 s, projected end point, hybrid side filter) and to a zone |
| `04_datasets/` | 2.5–2.7 | Builds the per-fish sample dataset (plasma, body mass, monoamines per brain region) from the laboratory workbooks, verifies it independently, and adds the Dl, Vv and POA monoamine concentrations |
| `05_statistics/` | 2.9.1, 2.9.2 | Mixed models of behaviour, aggression, plasma and monoamines (random-effect reduction, Box–Cox rule, Kenward–Roger F-tests, Wald χ², parametric-bootstrap LRT, Tukey contrasts, Rosner outlier screening), trial-level ANOVAs, hierarchy and timing analyses |
| `06_robustness_and_associations/` | 2.9.2, 2.9.3 | Robustness tests of every significant effect (exact permutation, exact Wilcoxon, leave-one-trial-out; Supp. Table S12) and the within-trial permutation tests of individual correlates of aggression (Supp. Table S13) |
| `07_figures/` | Figures 3–8 | Scripts that draw the data figures of the article |

## Article → script map

| Article item | Script | Notes |
|---|---|---|
| 2.8.1 Video tracking: per-clip idtracker.ai sessions | idtracker.ai 6.0.14 (GUI, one `.toml` per clip) | The short Python script that cut the full videos into 10-min clips is **not included** (see "Not included") |
| 2.8.1 Trajectory preprocessing (jumps, interpolation ≤ 1 s, cm, body length) | `02_trajectory_pipeline/01_preprocess_trajectories.R` (calls `01_tracking/extract_bounding_boxes.py`) | Step 1 of the pipeline |
| 2.8.1 Identity matching between intervals (idmatcher.ai) | `01_tracking/idmatcher_by_trial_colab.ipynb` → `02_trajectory_pipeline/01_5_annotate_identity_matches.R` (+ `helpers/annotate_idmatcher_scores.R`) | Annotation cascade, Rules 1–3, Supp. Table SM1 |
| Body length from tracking reconciled with body length at sampling | `02_trajectory_pipeline/01_6_match_body_length.R` | Hungarian matching; reliability tiers |
| 2.8.2 Individual indicators | `02_trajectory_pipeline/02_individual_indicators.R` | Step 2 |
| 2.8.2 School-level indicators | `02_trajectory_pipeline/06_school_indicators.R` | Step 6 (IID, NND, polarisation, school area, school speed, zone synchrony) |
| 2.8.3 Attacker, target and zone of each aggression event | `03_aggression_events/aggression_pipeline.py` (settings in `config.py`) | Output `master_events_v2_noTP1_noPharma_side_hybrid.csv` is the event table used by the statistics (`data/aggression_annotated_june2026.csv`) |
| 2.5–2.7 Sample dataset (plasma, body mass, monoamines) | `04_datasets/build_sample_dataset.py`, `verify_sample_dataset.py`, `merge_dl_vv_poa_concentrations.py` | |
| 2.9.1–2.9.2 Individual and school-level indicators (Supp. Tables S4, S5, S10) | `05_statistics/behaviour_models.R` | Reads the step 6 output of the pipeline |
| 2.9.2 Aggression counts with exposure offset, binomial flow-zone models (Supp. Table S6) | `05_statistics/aggression_models.R`, `05_statistics/aggression_models_design_random_effects.R` | The second script refits the count models with the design random effects (day, trial) |
| 2.9.2 Events per trial (two-way ANOVA) and body mass (two-way ANOVA), Box–Cox rule | `05_statistics/trial_level_anovas.R` | |
| 2.9.2 Hierarchy indices (permutation tests) and temporal clustering (Supp. Table S7) | `05_statistics/aggression_hierarchy_and_timing.R` | Exploratory |
| 2.9.2 Plasma and monoamines, Rosner outlier screening (Supp. Tables S8, S9, S11) | `05_statistics/plasma_monoamine_models.R` | |
| 2.9.2 Linking behaviour to the sampled fish | `05_statistics/link_behaviour_physiology.R` | Writes `c23_behaviour_gated.csv` used below |
| 2.9.2 Individual correlates of aggression (Supp. Table S13) | `06_robustness_and_associations/individual_aggression_correlates.py` | Within-trial permutation, 20 000 permutations |
| 2.9.3 Robustness of significant effects (Supp. Table S12) | `06_robustness_and_associations/robustness_significant_effects.py` | |
| Figure 3 (locomotion) | `07_figures/figure_3_locomotion.R` | Written as `Figure_5_*` (working number) |
| Figure 4 (school structure) | `07_figures/figure_4_school_structure.py` | |
| Figure 5 (aggressive events) | `07_figures/figure_5_aggression.py` | Written as `Figure_7_aggression.png` |
| Figures 6 (aggression in time and space) and 7 (plasma) | `07_figures/figures_6_7_aggression_time_space_and_plasma.py` | Written as `Figure_8_*` and `Figure_9_*` |
| Figure 8 (Dm serotonergic system) | `07_figures/figure_8_dm_serotonin.py` | Written as `Figure_10_Dm_serotonin.png` |

Shared code: `05_statistics/shared/house_style.R` (figure style, number formatting, effect sizes,
random-effect reduction, parametric bootstrap) and `05_statistics/shared/sensitivity_checks.R` (the
sensitivity checks written with every run) are sourced by the statistics scripts. `07_figures/figure_style.py`
and `07_figures/cld_letters.py` are shared by the Python figures. `05_statistics/validate_stat_types.R` checks
a finished run folder, and `05_statistics/tests/` holds unit tests for the random-effect reduction, the
1-df parametric bootstrap and the endocrine model selection (`Rscript 05_statistics/tests/test_*.R`).

Steps 4 and 5 of the trajectory pipeline (`04_exploratory_graphs.R`, `05_exploratory_individual_stats.R`) are
exploratory and kept so that the pipeline runs as it did. The statistics reported in the article come from
`05_statistics/`.

## Running the code

All paths are relative to the repository root. Run every command from there.

**Input data.** Put the data in a folder `data/` at the repository root:

```
data/
  tracking_sessions/                  idtracker.ai session folders, one per trial × interval
                                      (session_trial{N}_Segment{M}/trajectories/trajectories.npy|h5)
  idmatcher/idmatcherai outputs from colab/<timestamp>/   output of the idmatcher.ai notebook
  trials summary GOT.xlsx             trial, group, treatment and date per trial
  sampling data GOT.xlsx              body mass and length of the sampled fish
  GOT glu_lac_cort updated march 2025.xlsx   plasma cortisol, glucose, lactate
  monoamine concentration ALL GOT.xlsx
  raw/                                laboratory workbooks read by 04_datasets/ (HPLC concentrations,
                                      Bradford protein, punches) and Dl_Vv_POA_concentrations.csv
  aggression/                         BORIS export (new_aggression_15.6.2026.xlsx), identity-match and
                                      body-length tables read by 03_aggression_events/config.py
  aggression_annotated_june2026.csv   event table written by 03_aggression_events/ (see above)
```

**Order of execution.**

```
# 1. Tracking (idtracker.ai GUI), then identity matching: run 01_tracking/idmatcher_by_trial_colab.ipynb in Colab
# 2. Trajectory pipeline (steps 1, 1.5, 1.6, 2, 4, 5, 6)
Rscript -e 'source("02_trajectory_pipeline/00_run_trajectory_pipeline.R")'
# 3. Aggression events
cd 03_aggression_events && python aggression_pipeline.py all && cd ..
# 4. Sample dataset
python 04_datasets/build_sample_dataset.py
python 04_datasets/verify_sample_dataset.py
python 04_datasets/merge_dl_vv_poa_concentrations.py
# 5. Statistics (each run writes a time-stamped folder under output/)
Rscript 05_statistics/behaviour_models.R
GOT_MONO_CSV=data/GOT_monoamine_data_with_Dl_Vv_POA.csv Rscript 05_statistics/plasma_monoamine_models.R
Rscript 05_statistics/aggression_models.R
Rscript 05_statistics/aggression_models_design_random_effects.R data output/aggression_design_re
Rscript 05_statistics/aggression_hierarchy_and_timing.R
Rscript 05_statistics/trial_level_anovas.R output/aggression_design_re data/GOT_monoamine_data_AUG2026.csv
Rscript 05_statistics/link_behaviour_physiology.R
# 6. Robustness and individual associations; 7. figures
#    (each script prints its arguments: run folders from step 5)
```

`behaviour_models.R` reads, by default, the trajectory-pipeline runs used for the article
(`STEP6_output_20260821_124225`, `STEP1_5_20260821_122105`). After a new pipeline run, set `GOT_STEP6_DIR`,
`GOT_STEP1_5_DIR` and `GOT_STEP1_5_STAMP` to the new folders. `GOT_ROOT` (default: the working directory) and
`GOT_OUT_BASE` (default: `output/`) can be set the same way.

**Note.** In this repository the scripts were reorganised and their hard-coded paths were replaced by
paths relative to the repository root. The analysis code itself is unchanged. All R files parse, all
Python files compile and the unit tests pass. Several scripts were re-run on the original data in the new
layout and gave identical results (robustness tests, trial-level ANOVAs, aggression design random-effects
refit). The whole chain has not been re-run end to end in this layout.

Code comments sometimes cite numbered rules of the internal analysis plan (for example "HANDOFF S7",
"CONVENTIONS C12", "LOGS.md"). Those working documents are not part of this archive. The rules they refer
to are described in the Materials and Methods.

## Software

R 4.6.0 with lme4 2.0.1, lmerTest 3.2.1, glmmTMB 1.1.14, car 3.1.5, emmeans 2.0.3, multcomp 1.4.30,
DHARMa 0.4.7, pbkrtest 0.5.5, MuMIn 1.48.19, EnvStats 3.1.0, performance 0.16.0, survival 3.8.6,
effectsize 1.0.2, swaRm 0.6.0, swaRmverse 0.1.1, clue 0.3.68, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0,
purrr 1.2.2, data.table 1.18.2.1, jsonlite 2.0.0, ggplot2 4.0.3 and patchwork 1.3.2 (Supplementary Table S3).
The scripts also use readxl, writexl, officer, flextable, multcompView, cowplot, ggtext, scales, stringr,
janitor, future, future.apply and e1071 for input/output, reports and plots.

Python 3.13 with NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1, statsmodels 0.15.0, Matplotlib 3.11.2, openpyxl,
h5py and tabulate. The bounding-box helper and the idmatcher.ai notebook run in the idtracker.ai
environment (idtracker.ai 6.0.14). Aggression was scored in BORIS 8.21.10.

## Not included

- The Python script that cut the full trial videos into 10-min clips (Supplementary Material, "Automated
  identity tracking workflow").
- Scripts that only assembled the manuscript (Word documents, supplementary tables, Supplementary Figure S1)
  and earlier versions of the analysis.
- The data.

## Licence

MIT (see `LICENSE`).
