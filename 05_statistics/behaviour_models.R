################################################################################
# behaviour_models.R
# ─────────────────────────────────────────────────────────────────────────────
# Replication of FINAL_easy_script_JUN2026_AICc.R (easy_script_reference.R)
# applied to GOT fish-behaviour data.
#
# Analysis units:
#   A) Group-dynamics (GD): one row per trial × interval (session-level school)
#   B) Individual fish:     one row per fish × interval (TP2 + TP3)
#
# Two analysis levels per unit:
#   1) Interval-level:    fixed = Treatment * Group * interval
#   2) Trial-aggregated:  fixed = Treatment * Group
#
# Model selection: AICc over random-effect candidate structures (lmerTest/glmmTMB)
# Inference: Type III ANOVA (car::Anova for GLMM, anova(type=3) for LMM)
# Post-hoc:  emmeans + Tukey pairwise contrasts
#
# Outputs (all written to output/<timestamp>/):
#   - 4 analysis-ready CSVs (gd_interval, gd_trial, ind_interval, ind_trial)
#   - 2 markdown model reports (interval-level, trial-aggregated)
#   - 3 ANOVA-statement CSVs (UTF-8 BOM)
#   - 1 Word document with model-decision tables
#
# ── Deviations from reference ──────────────────────────────────────────────
# 1. GOT 2×2 factorial (Treatment × Group) replaces single treatment predictor.
# 2. No 'tank' RE; phys_trial → trial_id; RE sets have no tank candidate.
# 3. GD and individual intervals are {2, 3} (TP2/TP3 subset — pharma excluded).
# 4. Angular momentum dropped (all NA).
# 5. Fish identity: cross-segment master IDs from idmatcher (STEP1_5) link the same
#    fish across TP2/TP3 within a trial. The individual-interval RE set therefore
#    includes (1|trial_id)+(1|fish_id_matched_trial) as repeated-measures candidates
#    (AICc-selected). Unmatched fish (rejected_conflict / untracked, e.g. trial 8)
#    fall back to a unique singleton ID. Fish are never reused across trials.
# 6. Data subsetted to Exercise choice and No choice groups (Naloxone/Morphine
#    excluded, consistent with existing pipeline practice).
# 7. p-value floor for "< X" notation: 0.001 (handoff/CONVENTIONS.md C7/A-9;
#    superseded the prior 0.0001 project convention 2026-08-24).
# 8. Bounded-duration individual indicators (time in flow zone, time moving >1 BL/s)
#    are modelled as proportions of the recording window via Beta GLMM (not Gaussian):
#    interval bound = 600 s; trial bound = n_intervals × 600 s (= 1200 s here, since
#    only TP2+TP3 are present — 2 intervals). Raw seconds retained for reporting.
# 9. Time-in-calm-zone removed as an indicator (per request). Angular momentum
#    likewise dropped (all-NA). All other indicators retained and populated —
#    individual mean_nnd_cm is populated from per-fish NND (NND_fish_by_timepoint).
# 10. Transformations (Gaussian LME, to reduce residual non-normality):
#     - Box-Cox (car::bcPower) for Distance and Mean speed at the interval level.
#     - Yeo-Johnson (car::yjPower, lambda by ML via car::powerTransform family
#       "yjPower") for Mean acceleration at both interval and trial levels (handles
#       its negative values natively; analogous to the bcPower Box-Cox path).
#     These are applied as the selected transform for the named indicators.
# 11. Output binned into subfolders: datasets/ , anovas/ , reports/. Each markdown
#     report has an analogous Word (.docx) report with identical content.
################################################################################


################################################################################
###  LIBRARIES
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(car)        # Levene's test, Type-III ANOVA
library(tibble)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MuMIn)      # AICc()
library(emmeans)
library(MASS)
library(officer)    # Word document creation
library(flextable)  # Pretty tables in Word
library(pbkrtest)   # Kenward-Roger ddf
library(performance) # icc(), r2()

options(contrasts = c("contr.sum", "contr.poly"))
# Re-assert dplyr::select after MASS masks it
select <- dplyr::select

# -- House style / sensitivity modules (single source of truth; see the
#    approved implementation plan for why these replace what used to be six
#    independent local copies of the same formatters) --------------------
.env_or <- function(var, default) { v <- Sys.getenv(var, ""); if (nzchar(v)) v else default }
GOT_ROOT <- .env_or("GOT_ROOT", ".")
source(file.path(GOT_ROOT, "05_statistics/shared/house_style.R"))
source(file.path(GOT_ROOT, "05_statistics/shared/sensitivity_checks.R"))


################################################################################
###  PATHS
################################################################################

# ── Input data ────────────────────────────────────────────────────────────────
# The article used the trajectory-pipeline runs STEP6_output_20260821_124225 and
# STEP1_5_20260821_122105 (the default paths below). After a new run of the
# pipeline, point GOT_STEP6_DIR and GOT_STEP1_5_DIR (and GOT_STEP1_5_STAMP, the
# timestamp in the STEP 1.5 file names) at the new output folders.
STEP6_DIR <- .env_or("GOT_STEP6_DIR", file.path(
  GOT_ROOT, "02_trajectory_pipeline/STEP6_output",
  "STEP6_output_20260821_124225/01_data"
))

GD_SESSION_CSV <- file.path(
  STEP6_DIR, "session_gd_summary_TP2_TP3_ExeChoice_NoChoice.csv"
)
GD_NND_CSV <- file.path(STEP6_DIR, "NND_fish_by_timepoint.csv")

# Individual: pre-joined dataset with NND, IID per fish (its fish_id_matched_trial
# column is a bogus trial+fish concatenation and is rebuilt below from idmatcher).
IND_CSV <- file.path(
  STEP6_DIR, "individual_analysis_dataset_TP2_TP3_ExeChoice_NoChoice.csv"
)

# idmatcher (STEP1_5): cross-segment fish identity. tracker_fish_id (e.g. T5_S2_F1)
# → master fish_id_matched_trial (e.g. T5_MF1) stable across TP2/TP3 within a trial.
STEP1_5_STAMP <- .env_or("GOT_STEP1_5_STAMP", "20260821_122105")
STEP1_5_DIR <- .env_or("GOT_STEP1_5_DIR", file.path(
  GOT_ROOT, "02_trajectory_pipeline/STEP1_5", paste0("STEP1_5_", STEP1_5_STAMP)
))
IDMAP_CSV <- file.path(STEP1_5_DIR, paste0("idmatcher_id_map_", STEP1_5_STAMP, ".csv"))
TIER_CSV  <- file.path(STEP1_5_DIR, paste0("reliability_tier_", STEP1_5_STAMP, ".csv"))

# ── Output ────────────────────────────────────────────────────────────────────
# Results are binned into three subfolders: datasets / anovas / reports.
TIMESTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(
  .env_or("GOT_OUT_BASE", file.path(GOT_ROOT, "output")), TIMESTAMP
)
DIR_DATASETS <- file.path(OUTPUT_DIR, "datasets")
DIR_ANOVAS   <- file.path(OUTPUT_DIR, "anovas")
DIR_REPORTS  <- file.path(OUTPUT_DIR, "reports")
for (.d in c(OUTPUT_DIR, DIR_DATASETS, DIR_ANOVAS, DIR_REPORTS))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# Accumulator for aicc_selection.csv (CONVENTIONS.md C12): one row set per
# indicator, flushed to a single CSV in DIR_ANOVAS once all indicators have
# run (see the write-out beside the anova_stmt_*.csv exports below).
.AICC_LOG <- list()
.aicc_log_add <- function(label, aicc_tbl) {
  aicc_tbl$label <- label
  # HANDOFF S2/S5: every candidate level in this design is a DESIGN level, so
  # the candidate set is a single forced specification and AICc selects
  # nothing. Recorded explicitly so a reader cannot mistake this file for the
  # record of a search that never happened.
  aicc_tbl$nothing_selectable <- nrow(aicc_tbl) == 1L
  .AICC_LOG[[length(.AICC_LOG) + 1]] <<- aicc_tbl
}

# Accumulators for anova_satterthwaite.csv (HANDOFF S7) and the per-indicator
# anova.csv (effect sizes -- HANDOFF S6 [INV], CONVENTIONS.md C12).
.SATT_LOG <- list()
.satt_log_add <- function(label, satt_tbl) {
  satt_tbl$label <- label
  .SATT_LOG[[length(.SATT_LOG) + 1]] <<- satt_tbl
}
.ANOVA_CSV_LOG <- list()
.anova_csv_log_add <- function(label, tbl) {
  tbl$label <- label
  .ANOVA_CSV_LOG[[length(.ANOVA_CSV_LOG) + 1]] <<- tbl
}

# finalize_gaussian_inference(): HANDOFF S5/S7 -- Kenward-Roger is the
# PRIMARY denominator-df approximation for every Gaussian LMM; Satterthwaite
# is refit alongside it as the named sensitivity check and written to
# anova_satterthwaite.csv. Disagreement across the alpha=0.05 boundary is
# flagged (not buried -- HANDOFF S7 [INV]).
finalize_gaussian_inference <- function(model_final, label) {
  # A model whose random structure was reduced away entirely (S3.2b) is an
  # ordinary linear model. Kenward-Roger has nothing to correct there -- the
  # denominator df ARE the residual df -- but the indicator must still be
  # reported. Before this branch existed, such a model produced no ANOVA at
  # all and the indicator vanished from anova.csv without a word, which is a
  # worse failure than a wrong number: mean acceleration disappeared from the
  # results set entirely and nothing said so.
  if (!inherits(model_final, "merMod") && !inherits(model_final, "lmerModLmerTest")) {
    at <- tryCatch(as.data.frame(car::Anova(model_final, type = 3)),
                   error = function(e) NULL)
    if (is.null(at)) return(NULL)
    at$term <- rownames(at)
    df_res <- tryCatch(stats::df.residual(model_final), error = function(e) NA_real_)
    at <- at[!at$term %in% c("Residuals", "(Intercept)"), , drop = FALSE]
    out <- data.frame(
      NumDF = at$Df, DenDF = df_res,
      `F value` = at[["F value"]], `Pr(>F)` = at[["Pr(>F)"]],
      check.names = FALSE, row.names = at$term)
    cat("  (no random effect retained -- residual df =", df_res,
        "; Kenward-Roger not applicable)\n")
    return(out)
  }
  anova_kr <- tryCatch(anova(model_final, type = 3, ddf = "Kenward-Roger"),
                       error = function(e) NULL)
  anova_sat <- tryCatch(anova(model_final, type = 3, ddf = "Satterthwaite"),
                        error = function(e) NULL)
  if (!is.null(anova_sat)) {
    satt_tbl <- as.data.frame(anova_sat)
    satt_tbl$term <- rownames(satt_tbl)
    if (!is.null(anova_kr)) {
      kr_tbl <- as.data.frame(anova_kr); kr_tbl$term <- rownames(kr_tbl)
      m <- merge(satt_tbl, kr_tbl[, c("term", "Pr(>F)")], by = "term",
                suffixes = c("_satterthwaite", "_kr"))
      m$p_crosses_alpha <- (m$`Pr(>F)_satterthwaite` < 0.05) != (m$`Pr(>F)_kr` < 0.05)
      if (any(m$p_crosses_alpha, na.rm = TRUE))
        cat("  *** KR/Satterthwaite disagree across p=0.05 for:",
            paste(m$term[which(m$p_crosses_alpha)], collapse = ", "), "***\n")
      satt_tbl <- m
    }
    .satt_log_add(label, satt_tbl)
  }
  if (is.null(anova_kr)) anova_kr <- anova_sat   # KR unavailable -> fall back, still report
  anova_kr
}

# anova_effect_size_table(): one row per term -> anova.csv (CONVENTIONS.md
# C12). Dispatches on the table shape the same way format_anova_statement()
# does (NumDF/DenDF = KR/Satterthwaite F row; Df = car::Anova Wald-chisq row).
anova_effect_size_table <- function(anova_obj, model_obj, model_type, label, r2m = NA_real_,
                                    r2c = NA_real_, icc = NA_real_) {
  if (is.null(anova_obj)) return(invisible(NULL))
  # Name what was actually fitted, derived from the model rather than from the
  # caller's branch. A Gaussian model whose random structure was reduced away
  # (S3.2b) is an ordinary linear model, and anova.csv must say so: its
  # denominator counts rows, not trials. Deriving it here rather than in each
  # fitting branch is deliberate -- there are four Gaussian branches and an
  # earlier per-branch version of this relabel silently covered only one of
  # them, leaving mean acceleration described as a mixed model when it had no
  # random effect left.
  if (!is.null(model_obj) && grepl("^Gaussian", model_type) &&
      !inherits(model_obj, "merMod") && !inherits(model_obj, "lmerModLmerTest"))
    model_type <- "Gaussian LM (no RE retained)"
  tab <- as.data.frame(anova_obj); tab$term <- rownames(tab)
  rows <- lapply(seq_len(nrow(tab)), function(i) {
    row <- tab[i, , drop = FALSE]
    if (all(c("NumDF", "DenDF", "F value", "Pr(>F)") %in% names(row))) {
      es <- effect_size("F", F = row$`F value`, df1 = row$NumDF, df2 = row$DenDF)
      data.frame(term = row$term, stat_type = "F", statistic = row$`F value`,
                df1 = row$NumDF, df2 = row$DenDF, p = row$`Pr(>F)`,
                eta2_p = es$es, omega2_p = es$es_secondary)
    } else if (all(c("Chisq", "Df", "Pr(>Chisq)") %in% names(row))) {
      data.frame(term = row$term, stat_type = "chisq", statistic = row$Chisq,
                df1 = row$Df, df2 = NA_real_, p = row$`Pr(>Chisq)`,
                eta2_p = NA_real_, omega2_p = NA_real_)
    } else NULL
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df)) return(invisible(NULL))
  df$model_type <- model_type; df$R2m <- r2m; df$R2c <- r2c; df$ICC <- icc
  .anova_csv_log_add(label, df)
  invisible(df)
}
cat("\n=== behaviour_models.R ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("  subfolders: datasets / anovas / reports\n")


################################################################################
###  DATA IMPORT AND DATASET CONSTRUCTION
################################################################################

cat("\n--- Loading data ---\n")

# ── A. Group-dynamics session-level (interval-level) ──────────────────────────

gd_raw <- read_csv(GD_SESSION_CSV, show_col_types = FALSE)

# Per-fish NND (keyed trial_id × timepoint × numeric fish_id 1..4) — used both for
# the session-level GD mean and to populate the individual-level NND indicator.
nnd_fish_raw <- read_csv(GD_NND_CSV, show_col_types = FALSE) %>%
  filter(Group %in% c("Exercise choice", "No choice"))

# Session-level mean NND from per-fish NND data
nnd_session <- nnd_fish_raw %>%
  group_by(trial_id, timepoint) %>%
  summarise(mean_nnd_cm = mean(mean_nnd_cm, na.rm = TRUE), .groups = "drop")

# Per-fish NND keyed for the individual dataset join (fish number 1..4)
nnd_fish <- nnd_fish_raw %>%
  transmute(trial_id, timepoint, fish_num = as.integer(fish_id), mean_nnd_cm)

# idmatcher cross-segment identity map: fish_id (= tracker_fish_id) → master fish ID.
# master_fish_id is stable across TP2/TP3 within a trial; NA where idmatcher could
# not match a fish (rejected_conflict, or fish absent from a segment).
idmap <- read_csv(IDMAP_CSV, show_col_types = FALSE) %>%
  transmute(fish_id          = tracker_fish_id,
            master_fish_id   = fish_id_matched_trial,
            id_match_quality = match_quality)

# Per-trial idmatcher reliability tier (conservative > reliable > borderline > unreliable)
id_tier <- read_csv(TIER_CSV, show_col_types = FALSE) %>%
  transmute(trial_id, id_reliability_tier = idmatcher_tier)

first_val <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA else x[1] }

gd_interval <- gd_raw %>%
  # Drop angular momentum (all NA) and time_num (internal label)
  select(-any_of(c("mean_ang_momentum", "time_num"))) %>%
  # Join session-mean NND
  left_join(nnd_session, by = c("trial_id", "timepoint")) %>%
  mutate(
    interval   = as.factor(timepoint),
    trial_id   = as.factor(trial_id),
    trial_date = as.factor(as.Date(trial_date)),
    Treatment  = as.factor(Treatment),
    Group      = as.factor(Group)
  )

cat("gd_interval: N =", nrow(gd_interval), "rows |",
    "trials:", nlevels(gd_interval$trial_id),
    "intervals:", nlevels(gd_interval$interval),
    "groups:", nlevels(gd_interval$Group),
    "treatments:", nlevels(gd_interval$Treatment), "\n")

# ── A-agg. GD trial-aggregated ───────────────────────────────────────────────
gd_trial <- gd_interval %>%
  group_by(trial_id) %>%
  summarise(
    trial_date               = first_val(as.character(trial_date)),
    Treatment                = first_val(as.character(Treatment)),
    Group                    = first_val(as.character(Group)),
    intervals_included       = paste(sort(unique(as.character(interval))), collapse = ", "),
    n_intervals              = n_distinct(interval),
    mean_iid_cm              = mean(mean_iid_cm,              na.rm = TRUE),
    mean_nnd_cm              = mean(mean_nnd_cm,              na.rm = TRUE),
    mean_centroid_speed_cm_s = mean(mean_centroid_speed_cm_s, na.rm = TRUE),
    mean_hull_area_cm2       = mean(mean_hull_area_cm2,       na.rm = TRUE),
    mean_polarization        = mean(mean_polarization,        na.rm = TRUE),
    zone_synchrony           = mean(zone_synchrony,           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trial_date = as.factor(trial_date),
    Treatment  = as.factor(Treatment),
    Group      = as.factor(Group),
    trial_id   = as.factor(trial_id)
  )

cat("gd_trial: N =", nrow(gd_trial), "rows\n")

# ── B. Individual interval-level ──────────────────────────────────────────────
ind_raw <- read_csv(IND_CSV, show_col_types = FALSE)

# Each interval (timepoint) is a fixed-duration recording window.
# time_spent_zone_flow_s and time_moving_TS1_s are durations bounded [0, INTERVAL_SECONDS]
# → modelled as proportions of the window via Beta GLMM (see notes below).
INTERVAL_SECONDS <- 600   # 10-min recording window per interval

.clip01 <- function(p) pmin(pmax(p, 0), 1)

# Filter out any remaining exclusion flags if present in this file.
# - The source mean_nnd_cm column is all-NA; drop it and re-populate from per-fish NND.
# - The source fish_id_matched_trial is a bogus trial+fish concatenation; drop it and
#   rebuild a TRUE cross-segment master identity from idmatcher (stable across TP2/TP3).
ind_interval <- ind_raw %>%
  { if ("exclude_row" %in% names(.)) filter(., !exclude_row) else . } %>%
  { if ("flag_dead_immobile" %in% names(.)) filter(., !flag_dead_immobile) else . } %>%
  select(-any_of(c("mean_nnd_cm", "fish_id_matched_trial"))) %>%
  mutate(fish_num = as.integer(sub(".*_F(\\d+)$", "\\1", fish_id))) %>%
  left_join(nnd_fish, by = c("trial_id", "timepoint", "fish_num")) %>%
  left_join(idmap,   by = "fish_id") %>%
  left_join(id_tier, by = "trial_id") %>%
  mutate(
    # Cross-segment master fish identity. Use the idmatcher master ID where available;
    # otherwise fall back to the segment-specific fish_id (a unique singleton fish).
    fish_matched_flag     = !is.na(master_fish_id),
    fish_id_matched_trial = ifelse(fish_matched_flag, master_fish_id, fish_id),
    interval   = as.factor(timepoint),
    trial_id   = as.factor(trial_id),
    trial_date = as.factor(as.Date(trial_date)),
    Treatment  = as.factor(Treatment),
    Group      = as.factor(Group),
    fish_id    = as.factor(fish_id),
    fish_id_matched_trial = as.factor(fish_id_matched_trial),
    # Bounded-duration → proportion of the interval window (raw seconds retained above)
    time_spent_zone_flow_prop = .clip01(time_spent_zone_flow_s / INTERVAL_SECONDS),
    time_moving_TS1_prop      = .clip01(time_moving_TS1_s      / INTERVAL_SECONDS)
  )

# Repeated-measures coverage: how many master fish are observed at both TP2 and TP3
.fish_obs <- ind_interval %>%
  dplyr::count(fish_id_matched_trial, name = "n_obs")
cat("ind_interval: N =", nrow(ind_interval), "rows |",
    "trials:", nlevels(ind_interval$trial_id),
    "intervals:", nlevels(ind_interval$interval),
    "\n  fish identity: matched-to-master rows =", sum(ind_interval$fish_matched_flag),
    "| distinct master fish =", nlevels(ind_interval$fish_id_matched_trial),
    "| paired (2 obs) =", sum(.fish_obs$n_obs == 2),
    "| singletons (1 obs) =", sum(.fish_obs$n_obs == 1),
    "\n  NND non-NA:", sum(!is.na(ind_interval$mean_nnd_cm)), "of", nrow(ind_interval), "\n")

# ── B-agg. Individual trial-aggregated ───────────────────────────────────────
# Aggregate over all fish × intervals within each trial → one row per trial.
# Counts summed; continuous variables averaged.
#
# Bounded-duration indicators (flow time, time moving > 1 BL/s) are handled
# separately below: per the experimental definition, the trial-level value is the
# SUM across intervals.  Because fish are NOT matched across intervals in this
# dataset (each fish_id appears in a single interval), we sum the per-interval
# fish-mean across intervals → seconds bounded [0, n_intervals × INTERVAL_SECONDS].
# Data has 2 intervals (TP2, TP3) ⇒ bound = 1200 s.  The proportion (÷ that bound)
# is the response actually modelled (Beta GLMM); the summed seconds are retained
# for reporting.
ind_bounded_trial <- ind_interval %>%
  group_by(trial_id, interval) %>%
  summarise(
    flow_int = mean(time_spent_zone_flow_s, na.rm = TRUE),
    move_int = mean(time_moving_TS1_s,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(trial_id) %>%
  summarise(
    n_intervals_bounded    = n(),
    time_spent_zone_flow_s = sum(flow_int, na.rm = TRUE),
    time_moving_TS1_s      = sum(move_int, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    .bound = n_intervals_bounded * INTERVAL_SECONDS,
    time_spent_zone_flow_prop = .clip01(time_spent_zone_flow_s / .bound),
    time_moving_TS1_prop      = .clip01(time_moving_TS1_s      / .bound)
  ) %>%
  select(-.bound)

ind_trial <- ind_interval %>%
  group_by(trial_id) %>%
  summarise(
    trial_date               = first_val(as.character(trial_date)),
    Treatment                = first_val(as.character(Treatment)),
    Group                    = first_val(as.character(Group)),
    n_fish_intervals         = n(),
    distance_m               = mean(distance_m,              na.rm = TRUE),
    speed_cm_s_mean          = mean(speed_cm_s_mean,         na.rm = TRUE),
    accel_cm_s2_mean         = mean(accel_cm_s2_mean,        na.rm = TRUE),
    n_zone_transitions       = sum(n_zone_transitions,       na.rm = TRUE),
    mean_nnd_cm              = mean(mean_nnd_cm,             na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(ind_bounded_trial, by = "trial_id") %>%
  mutate(
    trial_date = as.factor(trial_date),
    Treatment  = as.factor(Treatment),
    Group      = as.factor(Group),
    trial_id   = as.factor(trial_id)
  )

cat("ind_trial: N =", nrow(ind_trial), "rows\n")

# Save analysis-ready datasets → datasets/ subfolder
write_csv(gd_interval,  file.path(DIR_DATASETS, "gd_interval.csv"))
write_csv(gd_trial,     file.path(DIR_DATASETS, "gd_trial.csv"))
write_csv(ind_interval, file.path(DIR_DATASETS, "ind_interval.csv"))
write_csv(ind_trial,    file.path(DIR_DATASETS, "ind_trial.csv"))
cat("Datasets saved to", DIR_DATASETS, "\n")

cat("\nGroup counts (gd_interval):\n")
print(table(gd_interval$Group, gd_interval$Treatment))
cat("\nGroup counts (ind_interval):\n")
print(table(ind_interval$Group, ind_interval$Treatment))


################################################################################
###  RANDOM-EFFECT CANDIDATE SETS
################################################################################
#
# RANDOM-EFFECT STRUCTURE -- design-determined, not AICc-selected.
#
# The GOT design, established from datasets/gd_trial.csv and design_balance.csv
# rather than assumed:
#
#   housing tank (2)  -> Treatment   n = 1 unit per level   [confounded]
#        |
#   date (6)          -> Group       3 dates per level      [Group varies BETWEEN dates]
#        |
#   trial (12)        -> 2 per date, one Control + one Exercise
#        |                           [Treatment varies WITHIN date]
#   fish (48)         -> 4 per trial
#        |
#   interval (2)      -> within fish
#
# trial_id is NESTED in trial_date (2 trials per date) and trial_id labels are
# globally unique, so (1|trial_date) + (1|trial_id) already expresses
# (1|trial_date/trial_id) -- no implicit-nesting trap.
#
# HANDOFF.md S2 [INV]: "the unit of replication is the level at which treatment
# was applied ... forced into every candidate, never selected by model
# comparison." Applying that to BOTH factors gives two design random effects,
# not one:
#
#   (1|trial_date)  -- unit of replication for GROUP (Group varies only between
#                      dates; the 6 Exercise-choice trials come from just 3
#                      dates, so testing Group against between-TRIAL variance
#                      is pseudoreplication and yields an anticonservative
#                      F(1,8) where the honest denominator is ~4 df)
#   (1|trial_id)    -- level at which TREATMENT varies (within date)
#
# An earlier revision of this script treated trial_date as an AICc-selectable
# nuisance term on the grounds that it is 1:1 with Group and therefore
# "competes with" the Group fixed effect. That reasoning was backwards:
# competing with Group is precisely what a whole-plot error term does, and
# Group remains estimable because dates are NESTED in Group, not crossed with
# it.
#
# SINGULARITY. (1|trial_date) is singular in every candidate it appeared in
# (8/8 interval-level, 8/8 trial-level -- between-date variance estimated at
# zero throughout). This does NOT license dropping it: with 6 dates and 3 per
# arm, a REML boundary estimate of zero is what the design produces whether or
# not the true variance is zero. S2 forces the term; S5's guard retains a
# singular structure when no eligible alternative exists. The boundary estimate
# is reported (aicc_selection.csv) and the design-based split-plot degrees of
# freedom ship as a declared comparator in sensitivity_forced_re.csv.
#
# REPEATED MEASURES. Verified against datasets/ind_interval.csv:
#   - fish_id is SEGMENT-SCOPED: 96 unique values for 96 rows. As a grouping
#     variable it gives 96 singleton clusters, perfectly confounded with
#     residual. It is never used as an RE (asserted below).
#   - fish_id_matched_trial is the cross-interval identity and is INCOMPLETE:
#     54 clusters -- 42 fish observed at both intervals, 12 singletons where
#     the matcher failed (fish_matched_flag 84 TRUE / 12 FALSE). The
#     repeated-measures term is informed by 42 of 48 fish.
#   - With k = 2 intervals a random intercept, compound symmetry, AR(1) and an
#     unstructured correlation are MATHEMATICALLY IDENTICAL. There is no
#     correlation structure to choose, and none is chosen.
#   - A random slope for interval over fish is not identifiable (2 parameters
#     against 2 observations per fish) and is not fitted.
#   - A fish intercept spanning intervals is itself an IDENTITY-DEPENDENT
#     modelling choice under HANDOFF S4.1 [INV]. Identity is validated (the
#     matcher writes a per-fish reliability tier), so it is licensed -- and the
#     licence is checked, not asserted, by the identity-tier refit in the
#     sensitivity block.
#
# CONSEQUENCE: every candidate level in this design is a design level, so the
# AICc candidate set collapses to ONE specification per tier and S5's selection
# step is a NO-OP for this dataset. That is a finding, not an omission:
# aicc_selection.csv records the forced structure and its singularity flags
# with an explicit nothing_selectable column, and the Methods says AICc
# selected nothing because nothing was selectable.
#
# Trial-level tiers (one row per trial) keep trial_id as a FIXED effect
# (implicit in Treatment*Group) since a random intercept over singleton
# clusters is unidentifiable -- HANDOFF S2 [DR]. trial_date remains the design
# RE for Group at that tier.

# GD interval: one row per trial x interval (24 rows). No fish level exists.
RE_GD_INTERVAL <- c(
  forced_design = "(1 | trial_date) + (1 | trial_id)"
)

# GD trial: one row per trial (12 rows). trial_id unidentifiable -> fixed.
RE_GD_TRIAL <- c(
  forced_design = "(1 | trial_date)"
)

# Individual interval: one row per fish x interval (96 rows).
RE_IND_INTERVAL <- c(
  forced_design = "(1 | trial_date) + (1 | trial_id) + (1 | fish_id_matched_trial)"
)

# Individual trial: one row per trial (12 rows). trial_id unidentifiable ->
# fixed; fish identity not meaningful post-aggregation.
RE_IND_TRIAL <- c(
  forced_design = "(1 | trial_date)"
)

# Declared comparators for sensitivity_forced_re.csv. The primary structure
# cannot be varied by selection any more, so the robustness question moves
# here: what would the conclusion have been under each structure the design
# argument rules out? (HANDOFF S7 forced-RE check.)
RE_COMPARATORS_GD_INTERVAL <- c(
  date_dropped = "(1 | trial_id)"
)
RE_COMPARATORS_IND_INTERVAL <- c(
  date_dropped = "(1 | trial_id) + (1 | fish_id_matched_trial)",
  fish_dropped = "(1 | trial_date) + (1 | trial_id)"
)

# Guard: no RE may ever be built on the segment-scoped fish_id (96 singleton
# clusters). Fails the run rather than fitting a degenerate model.
.assert_no_raw_fish_re <- function(re_strings) {
  bad <- grep("fish_id", re_strings, value = TRUE)
  bad <- bad[!grepl("fish_id_matched_trial", bad)]
  if (length(bad))
    stop("Random effect on segment-scoped fish_id (96 singleton clusters): ",
         paste(bad, collapse = "; "), call. = FALSE)
  invisible(TRUE)
}
.assert_no_raw_fish_re(c(RE_GD_INTERVAL, RE_GD_TRIAL, RE_IND_INTERVAL, RE_IND_TRIAL,
                         RE_COMPARATORS_GD_INTERVAL, RE_COMPARATORS_IND_INTERVAL))


################################################################################
###  SHARED HELPER FUNCTIONS
###  Ported from easy_script_reference.R; adapted for GOT Treatment/Group 2×2.
################################################################################


# Smithson-Verkuilen beta adjustment
adjust_beta_01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

# Levene's test — returns p-value (NA if < 2 groups or error)
levene_p <- function(y, group) {
  d <- data.frame(y = y, g = as.factor(group))[is.finite(y), ]
  if (nlevels(d$g) < 2) return(NA_real_)
  tryCatch(
    car::leveneTest(y ~ g, data = d)$`Pr(>F)`[1],
    error = function(e) NA_real_
  )
}

# Shapiro-Wilk on model residuals
shapiro_residual_p <- function(model) {
  res <- residuals(model)
  res <- res[is.finite(res)]
  if (length(res) < 3) return(NA_real_)
  if (length(res) > 5000) res <- sample(res, 5000)
  tryCatch(shapiro.test(res)$p.value, error = function(e) NA_real_)
}

# Box-Cox transform estimated by maximum likelihood to maximise the normality of
# OLS residuals under the given fixed-effects model (car::powerTransform).
#   - strictly positive y  → classic Box-Cox          (car::bcPower)
#   - y with zeros/negatives → Box-Cox with negatives  (car::bcnPower; estimates a
#                              location parameter gamma so the transform is defined)
# Returns list(y_bc, name) on the finite subset of y, or NULL on failure.
estimate_boxcox <- function(y_valid, fixed_str, data_fit) {
  # lambda chosen as in Bellio et al. (2026): grid -2 to 2 in steps of 0.1, keeping the value that
  # maximises the Shapiro-Wilk p of the residuals of the fixed-effects model. Responses with zeros or
  # negatives use the Box-Cox-with-negatives family (Hawkins & Weisberg 2017), with the location
  # parameter gamma estimated by maximum likelihood (car::powerTransform) and lambda from the grid.
  d <- data_fit
  fml <- as.formula(paste(".y ~", fixed_str))
  lam_grid <- seq(-2, 2, by = 0.1)
  grid_best <- function(tf) {
    sw <- vapply(lam_grid, function(l) {
      yt <- tryCatch(tf(l), error = function(e) NA)
      if (!all(is.finite(yt))) return(NA_real_)
      d$.y <- yt
      tryCatch(shapiro.test(residuals(lm(fml, data = d)))$p.value, error = function(e) NA_real_)
    }, numeric(1))
    if (all(is.na(sw))) return(NA_real_)
    round(lam_grid[which.max(sw)], 1)
  }
  tryCatch({
    if (min(y_valid, na.rm = TRUE) > 0) {
      lam <- grid_best(function(l) car::bcPower(y_valid, l))
      if (is.na(lam)) return(NULL)
      y_bc <- car::bcPower(y_valid, lam)
      if (!all(is.finite(y_bc))) return(NULL)
      list(y_bc = y_bc, lambda = lam, name = sprintf("box-cox (lambda=%.1f)", lam))
    } else {
      d$.y <- y_valid
      pt  <- car::powerTransform(lm(fml, data = d), family = "bcnPower")
      cf  <- coef(pt); cfv <- if (is.matrix(cf)) cf[1, ] else cf
      gam <- as.numeric(cfv[["gamma"]])
      lam <- grid_best(function(l) car::bcnPower(y_valid, lambda = l, gamma = gam))
      if (is.na(lam)) return(NULL)
      y_bc <- car::bcnPower(y_valid, lambda = lam, gamma = gam)
      if (!all(is.finite(y_bc))) return(NULL)
      list(y_bc = y_bc, lambda = lam,
           name = sprintf("box-cox-neg (lambda=%.1f, gamma=%.3f)", lam, gam))
    }
  }, error = function(e) NULL)
}

# Yeo-Johnson transform with ML-estimated lambda via car::powerTransform
# (family = "yjPower"), estimated on the OLS residuals under the fixed-effects model
# — exactly analogous to the Box-Cox path (bcPower). Natively handles the negative
# values in acceleration. Returns list(y_t, name) on the finite subset, or NULL.
estimate_yeojohnson <- function(y_valid, fixed_str, data_fit) {
  d <- data_fit
  d$.y <- y_valid
  fml <- as.formula(paste(".y ~", fixed_str))
  tryCatch({
    m   <- lm(fml, data = d)
    pt  <- car::powerTransform(m, family = "yjPower")
    lam <- as.numeric(coef(pt))[1]
    y_t <- car::yjPower(y_valid, lambda = lam)
    if (!all(is.finite(y_t))) return(NULL)
    list(y_t = y_t, name = sprintf("yeo-johnson (lambda=%.3f)", lam))
  }, error = function(e) NULL)
}

# Normality check + optional log1p/sqrt/box-cox/yeo-johnson transform
# Adapted: GOT uses Treatment × Group as Levene grouping (2×2 factorial)
# try_boxcox / try_yeojohnson: when TRUE the corresponding transform is computed and,
# if successful, becomes the selected transform (explicit per-indicator request).
check_and_transform <- function(y, fixed_str, data,
                                try_boxcox = FALSE, try_yeojohnson = FALSE) {
  y_valid <- y[is.finite(y)]
  if (length(y_valid) < 8)
    return(list(y_trans = y, transform_name = "none",
                sw_p_raw = NA, sw_p_trans = NA,
                levene_p_raw = NA, levene_p_trans = NA,
                normality_flag = FALSE))

  .fit_resid <- function(y_v, d) {
    d$.y <- y_v
    tryCatch({
      mod <- lm(as.formula(paste(".y ~", fixed_str)), data = d)
      r <- residuals(mod)[is.finite(residuals(mod))]
      if (length(r) > 5000) r <- sample(r, 5000)
      r
    }, error = function(e) y_v - mean(y_v, na.rm = TRUE))
  }

  .sw_p <- function(r) tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_)

  data_fit <- data[is.finite(y), ]
  # GOT: Treatment × Group grouping (replaces reference's treatment × zone)
  .grp <- tryCatch(
    interaction(data_fit[, intersect(c("Treatment", "Group"), names(data_fit))],
                drop = TRUE),
    error = function(e) factor(rep("all", nrow(data_fit)))
  )

  r_raw   <- .fit_resid(y_valid, data_fit)
  sw_raw  <- .sw_p(r_raw)
  lev_raw <- levene_p(y_valid, .grp)

  # SINGLE-TRANSFORMATION RULE (2026-11-12): if the raw response fails either
  # assumption (Shapiro-Wilk on OLS residuals p < 0.05, or Levene across the
  # Treatment x Group cells p < 0.05), it is Box-Cox transformed with lambda
  # estimated by maximum likelihood (car::powerTransform; bcnPower with an
  # estimated location gamma when the response has zero or negative values).
  # No other transformation is considered; the try_boxcox / try_yeojohnson
  # per-indicator flags are no longer used.
  best_name <- "none"
  best_sw   <- sw_raw
  best_lev  <- lev_raw
  best_lambda <- NA_real_
  y_trans   <- y
  fails <- (!is.na(sw_raw) && sw_raw < 0.05) || (!is.na(lev_raw) && lev_raw < 0.05)
  if (fails) {
    bc <- estimate_boxcox(y_valid, fixed_str, data_fit)
    if (!is.null(bc)) {
      best_sw   <- .sw_p(.fit_resid(bc$y_bc, data_fit))
      best_lev  <- levene_p(bc$y_bc, .grp)
      best_name <- bc$name
      best_lambda <- bc$lambda
      y_trans[is.finite(y)] <- bc$y_bc
    }
  }

  norm_flag <- (!is.na(best_sw) && best_sw < 0.05) ||
               (!is.na(best_lev) && best_lev < 0.05)

  list(y_trans = y_trans, transform_name = best_name, lambda = best_lambda,
       sw_p_raw = sw_raw, sw_p_trans = best_sw,
       levene_p_raw = lev_raw, levene_p_trans = best_lev,
       normality_flag = norm_flag)
}

# Gaussian LMM: AICc-based RE selection + normality check
# Adapted: checks !is.na(Treatment) instead of !is.na(treatment)
fit_gaussian_aicc <- function(data, response_var, fixed_str, re_cands,
                              label = "", skip_transform = FALSE,
                              try_boxcox = FALSE, try_yeojohnson = FALSE) {
  y_raw <- data[[response_var]]
  valid <- is.finite(y_raw) & !is.na(data$Treatment)
  d_fit <- data[valid, , drop = FALSE]
  y_v   <- y_raw[valid]

  if (skip_transform) {
    sw_raw <- tryCatch(shapiro.test(y_v)$p.value, error = function(e) NA_real_)
    trans <- list(y_trans = y_v, transform_name = "none",
                  sw_p_raw = sw_raw, sw_p_trans = sw_raw,
                  levene_p_raw = NA_real_, levene_p_trans = NA_real_,
                  normality_flag = (!is.na(sw_raw) && sw_raw < 0.05))
  } else {
    trans <- check_and_transform(y_v, fixed_str, d_fit,
                                 try_boxcox = try_boxcox,
                                 try_yeojohnson = try_yeojohnson)
  }
  d_fit$.y <- trans$y_trans

  cat("\n  [", label, "] transform:", trans$transform_name,
      " | SW(raw) p=", round(trans$sw_p_raw, 3),
      " | SW(trans) p=", round(trans$sw_p_trans, 3),
      " | Levene(trans) p=", round(trans$levene_p_trans, 3),
      if (!is.na(trans$normality_flag) && trans$normality_flag)
        "  *** NORMALITY FLAG ***" else "", "\n")

  # HANDOFF.md S5 step 3: select among RE candidates by AICc, fit with ML --
  # REML likelihoods are not comparable across differing RE structures, which
  # is what the previous REML-throughout version of this fitter did.
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(".y ~", fixed_str, "+", re_cands[[rn]]))
    m    <- tryCatch(lmerTest::lmer(fmla, data = d_fit, REML = FALSE),
                     error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)

  if (length(cand_results) == 0) {
    warning("No LMM converged for ", label); return(NULL)
  }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")

  # HANDOFF.md S5 step 4 -- the singularity guard (house_style.R A11):
  # some-but-not-all singular -> singular candidates ineligible; all singular
  # -> keep the full set rather than discarding everything.
  guard <- apply_singularity_guard(cand_results, lme4::isSingular)
  eligible_idx <- which(guard$eligible)
  best_idx     <- eligible_idx[which.min(aicc_vals[eligible_idx])]
  best         <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE            = sapply(cand_results, `[[`, "name"),
    re_formula    = sapply(cand_results, `[[`, "re"),
    fixed_formula = fixed_str,
    AICc          = round(aicc_vals, 2),
    delta_AICc    = round(aicc_vals - min(aicc_vals[eligible_idx]), 2),
    singular      = guard$singular,
    eligible      = guard$eligible,
    selected      = seq_along(cand_results) == best_idx
  )
  cat("  AICc table (ML selection):\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name, "(", best$re, ")",
      if (guard$all_singular) " -- ALL CANDIDATES SINGULAR, guard kept full set" else "", "\n")
  .aicc_log_add(label, aicc_tbl)

  # HANDOFF.md S5 step 5: refit the SELECTED structure with REML before
  # inference -- then REDUCE it if the fit is singular.
  #
  # A singular fit means the data do not identify that random-effect structure.
  # Kenward-Roger, which supplies the denominator df for every Gaussian model
  # here, is derived assuming the variance parameters lie in the INTERIOR of
  # the parameter space; at a boundary estimate the correction is undefined, so
  # keeping a zero-variance term corrupts the very machinery it was kept to
  # protect. reduce_singular_re() drops the offending component and refits,
  # iteratively, and logs every drop. See Methods 2.9.1.
  .fit_reml <- function(re_str) {
    f <- if (is.null(re_str)) as.formula(paste(".y ~", fixed_str))
         else as.formula(paste(".y ~", fixed_str, "+", re_str))
    if (is.null(re_str)) stats::lm(f, data = d_fit)
    else lmerTest::lmer(f, data = d_fit, REML = TRUE)
  }
  .red <- reduce_singular_re(.fit_reml, split_re_terms(best$re), label = label)
  model_final <- .red$model
  if (is.null(model_final)) model_final <- best$model   # ML fit as a last-resort fallback
  re_final <- if (!is.null(.red$re) && !is.na(.red$re)) .red$re else best$re
  singular_final <- tryCatch(
    if (inherits(model_final, "merMod")) lme4::isSingular(model_final) else FALSE,
    error = function(e) NA)
  if (length(.red$dropped)) {
    cat("  RE reduced:", best$re, "->", re_final,
        "(dropped", length(.red$dropped), "singular component(s))\n")
    re_reduction_log_add(label, "Gaussian", best$re, re_final, .red$dropped, .red$reasons)
  }

  list(model          = model_final,
       re_chosen      = if (length(.red$dropped)) paste0(best$name, "_reduced") else best$name,
       re_formula     = re_final,
       re_dropped     = if (length(.red$dropped)) paste(.red$dropped, collapse = " + ") else NA_character_,
       transform_name = trans$transform_name,
       sw_p_raw       = trans$sw_p_raw,
       sw_p_trans     = trans$sw_p_trans,
       levene_p_raw   = trans$levene_p_raw,
       bc_lambda      = if (is.null(trans$lambda)) NA_real_ else trans$lambda,
       levene_p_trans = trans$levene_p_trans,
       normality_flag = trans$normality_flag,
       singular       = singular_final,
       aicc_table     = aicc_tbl)
}

# Beta GLMM: AICc-based RE selection with the same forced-design-RE candidate
# sets and singularity guard as the Gaussian fitter. glmmTMB already fits by
# ML (no REML concept for these families), so no separate refit step is
# needed -- the AICc-selected model IS the inference model.
fit_beta_aicc <- function(data, response_var, fixed_str, re_cands, label = "") {
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(response_var, "~", fixed_str, "+", re_cands[[rn]]))
    m    <- tryCatch(
      glmmTMB(fmla, data = data, family = beta_family(link = "logit")),
      error = function(e) NULL
    )
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)
  if (length(cand_results) == 0) {
    warning("No Beta GLMM converged for ", label); return(NULL)
  }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  guard <- apply_singularity_guard(cand_results, function(m)
    tryCatch(glmmTMB::isSingular(m), error = function(e) FALSE))
  eligible_idx <- which(guard$eligible)
  best_idx     <- eligible_idx[which.min(aicc_vals[eligible_idx])]
  best         <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE         = sapply(cand_results, `[[`, "name"),
    re_formula = sapply(cand_results, `[[`, "re"),
    AICc       = round(aicc_vals, 2),
    delta_AICc = round(aicc_vals - min(aicc_vals[eligible_idx]), 2),
    singular   = guard$singular,
    eligible   = guard$eligible,
    selected   = seq_along(cand_results) == best_idx
  )
  cat("  AICc table (beta):\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name,
      if (guard$all_singular) " -- ALL CANDIDATES SINGULAR, guard kept full set" else "", "\n")
  .aicc_log_add(label, aicc_tbl)

  # Reduce a singular structure -- same argument as the Gaussian fitter. A
  # boundary variance estimate is a statement that the data do not identify the
  # term, and glmmTMB's standard errors (hence every Wald chi-square reported
  # from these models) come from a Hessian that is not positive definite when a
  # component sits on the boundary.
  .fit_beta <- function(re_str) {
    f <- if (is.null(re_str)) as.formula(paste(response_var, "~", fixed_str))
         else as.formula(paste(response_var, "~", fixed_str, "+", re_str))
    glmmTMB(f, data = data, family = beta_family(link = "logit"))
  }
  .red <- reduce_singular_re(.fit_beta, split_re_terms(best$re), label = label)
  m_final  <- if (!is.null(.red$model)) .red$model else best$model
  re_final <- if (!is.null(.red$re) && !is.na(.red$re)) .red$re else best$re
  if (length(.red$dropped)) {
    cat("  RE reduced:", best$re, "->", re_final,
        "(dropped", length(.red$dropped), "singular component(s))\n")
    re_reduction_log_add(label, "Beta GLMM", best$re, re_final, .red$dropped, .red$reasons)
  }
  list(model = m_final,
       re_chosen = if (length(.red$dropped)) paste0(best$name, "_reduced") else best$name,
       re_formula = re_final,
       re_dropped = if (length(.red$dropped)) paste(.red$dropped, collapse = " + ") else NA_character_,
       singular = isTRUE(.red$singular), aicc_table = aicc_tbl)
}

# Negative-binomial GLMM: AICc RE selection, same forced-RE candidates and
# singularity guard as the other two fitters.
fit_negbin_aicc <- function(data, response_var, fixed_str, re_cands, label = "") {
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(response_var, "~", fixed_str, "+", re_cands[[rn]]))
    m    <- tryCatch(
      glmmTMB(fmla, data = data, family = nbinom2(link = "log")),
      error = function(e) NULL
    )
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)
  if (length(cand_results) == 0) {
    warning("No NB GLMM converged for ", label); return(NULL)
  }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  guard <- apply_singularity_guard(cand_results, function(m)
    tryCatch(glmmTMB::isSingular(m), error = function(e) FALSE))
  eligible_idx <- which(guard$eligible)
  best_idx     <- eligible_idx[which.min(aicc_vals[eligible_idx])]
  best         <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE         = sapply(cand_results, `[[`, "name"),
    re_formula = sapply(cand_results, `[[`, "re"),
    AICc       = round(aicc_vals, 2),
    delta_AICc = round(aicc_vals - min(aicc_vals[eligible_idx]), 2),
    singular   = guard$singular,
    eligible   = guard$eligible,
    selected   = seq_along(cand_results) == best_idx
  )
  cat("  AICc table (negbin):\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name,
      if (guard$all_singular) " -- ALL CANDIDATES SINGULAR, guard kept full set" else "", "\n")
  .aicc_log_add(label, aicc_tbl)

  # Reduce a singular structure -- same argument as the other two fitters.
  .fit_nb <- function(re_str) {
    f <- if (is.null(re_str)) as.formula(paste(response_var, "~", fixed_str))
         else as.formula(paste(response_var, "~", fixed_str, "+", re_str))
    glmmTMB(f, data = data, family = glmmTMB::nbinom2())
  }
  .red <- reduce_singular_re(.fit_nb, split_re_terms(best$re), label = label)
  m_final  <- if (!is.null(.red$model)) .red$model else best$model
  re_final <- if (!is.null(.red$re) && !is.na(.red$re)) .red$re else best$re
  if (length(.red$dropped)) {
    cat("  RE reduced:", best$re, "->", re_final,
        "(dropped", length(.red$dropped), "singular component(s))\n")
    re_reduction_log_add(label, "Negative binomial GLMM", best$re, re_final,
                          .red$dropped, .red$reasons)
  }
  list(model = m_final,
       re_chosen = if (length(.red$dropped)) paste0(best$name, "_reduced") else best$name,
       re_formula = re_final,
       re_dropped = if (length(.red$dropped)) paste(.red$dropped, collapse = " + ") else NA_character_,
       singular = isTRUE(.red$singular), aicc_table = aicc_tbl)
}


################################################################################
###  GENERIC ANALYSIS LOOP FUNCTION
###  Dispatches to the correct model family based on variable type.
###  Returns a named list: model, anova, emm, posthoc, aicc_tbl, model_row
################################################################################

run_indicator <- function(v, data, fixed_str, re_cands,
                          count_vars, force_beta_vars, force_gaussian_vars,
                          emm_spec, boxcox_vars = character(0),
                          yeojohnson_vars = character(0)) {
  cat("\n\n============================================================\n")
  cat("Variable:", v, "\n")
  cat("============================================================\n")

  use_boxcox <- v %in% boxcox_vars
  use_yj     <- v %in% yeojohnson_vars

  # Select only needed columns for the model data frame
  keep_cols <- intersect(
    c("trial_id", "trial_date", "Treatment", "Group", "interval", "fish_id",
      "fish_id_matched_trial"),
    names(data)
  )
  df_var <- data %>%
    dplyr::select(all_of(keep_cols), all_of(v)) %>%
    filter(!is.na(.data[[v]]))

  var_min <- min(df_var[[v]], na.rm = TRUE)
  var_max <- max(df_var[[v]], na.rm = TRUE)
  cat("N =", nrow(df_var),
      "| Range:", round(var_min, 4), "to", round(var_max, 4), "\n")

  .ass <- list()
  model_final <- NULL; model_type <- ""; transform <- ""; shapiro_p <- NA_real_
  norm_flag <- NA; re_chosen <- ""; re_formula_str <- ""; aicc_tbl <- NULL
  anova_final <- NULL; emm_final <- NULL; posthoc_final <- NULL

  # ── Count variable ─────────────────────────────────────────────────────────
  if (v %in% count_vars) {
    cat("Decision: count — negative binomial GLMM\n")
    fit_res <- fit_negbin_aicc(df_var, v, fixed_str, re_cands, label = v)
    if (is.null(fit_res)) return(NULL)
    model_final   <- fit_res$model; model_type <- "Negative binomial GLMM"
    transform     <- "None"; shapiro_p <- NA_real_; norm_flag <- NA
    re_chosen     <- fit_res$re_chosen; re_formula_str <- fit_res$re_formula
    anova_final   <- tryCatch(car::Anova(model_final, type = 3), error = function(e) NULL)
    emm_final     <- tryCatch(emmeans(model_final, emm_spec, type = "response"),
                              error = function(e) NULL)
    posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"),
                              error = function(e) NULL)
    aicc_tbl      <- fit_res$aicc_table
    r2_nb <- tryCatch(MuMIn::r.squaredGLMM(model_final), error = function(e) NULL)
    anova_effect_size_table(anova_final, model_final, model_type, label = v,
                            r2m = if (!is.null(r2_nb)) unname(r2_nb[1, "R2m"]) else NA_real_,
                            r2c = if (!is.null(r2_nb)) unname(r2_nb[1, "R2c"]) else NA_real_,
                            icc = tryCatch(performance::icc(model_final)$ICC_adjusted, error = function(e) NA_real_))

  # ── Force beta (proportion 0-1) ────────────────────────────────────────────
  } else if (v %in% force_beta_vars) {
    cat("Decision: proportion — beta GLMM\n")
    response_beta <- paste0(v, "_beta")
    df_var <- df_var %>% mutate("{response_beta}" := adjust_beta_01(.data[[v]]))
    fit_res <- fit_beta_aicc(df_var, response_beta, fixed_str, re_cands, label = v)
    if (is.null(fit_res)) return(NULL)
    model_final   <- fit_res$model; model_type <- "Beta GLMM"
    transform     <- "Adjusted 0-1 beta"; shapiro_p <- NA_real_; norm_flag <- NA
    re_chosen     <- fit_res$re_chosen; re_formula_str <- fit_res$re_formula
    anova_final   <- tryCatch(car::Anova(model_final, type = 3), error = function(e) NULL)
    emm_final     <- tryCatch(emmeans(model_final, emm_spec, type = "response"),
                              error = function(e) NULL)
    posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"),
                              error = function(e) NULL)
    aicc_tbl      <- fit_res$aicc_table
    r2_beta <- tryCatch(MuMIn::r.squaredGLMM(model_final), error = function(e) NULL)
    anova_effect_size_table(anova_final, model_final, model_type, label = v,
                            r2m = if (!is.null(r2_beta)) unname(r2_beta[1, "R2m"]) else NA_real_,
                            r2c = if (!is.null(r2_beta)) unname(r2_beta[1, "R2c"]) else NA_real_,
                            icc = tryCatch(performance::icc(model_final)$ICC_adjusted, error = function(e) NA_real_))

  # ── Force Gaussian (pre-transformed; skip normality check) ─────────────────
  } else if (v %in% force_gaussian_vars) {
    cat("Decision: forced Gaussian LME (pre-transformed)\n")
    fit_res <- fit_gaussian_aicc(df_var, v, fixed_str, re_cands,
                                 label = v, skip_transform = TRUE)
    if (is.null(fit_res)) return(NULL)
    model_final   <- fit_res$model; model_type <- "Gaussian LME"
    .ass <- fit_res
    transform     <- fit_res$transform_name; shapiro_p <- shapiro_residual_p(model_final)
    norm_flag     <- fit_res$normality_flag
    re_chosen     <- fit_res$re_chosen; re_formula_str <- fit_res$re_formula
    anova_final   <- finalize_gaussian_inference(model_final, label = v)
    emm_final     <- tryCatch(emmeans(model_final, emm_spec), error = function(e) NULL)
    posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"), error = function(e) NULL)
    aicc_tbl      <- fit_res$aicc_table
    r2_g <- tryCatch(MuMIn::r.squaredGLMM(model_final), error = function(e) NULL)
    anova_effect_size_table(anova_final, model_final, model_type, label = v,
                            r2m = if (!is.null(r2_g)) unname(r2_g[1, "R2m"]) else NA_real_,
                            r2c = if (!is.null(r2_g)) unname(r2_g[1, "R2c"]) else NA_real_,
                            icc = tryCatch(performance::icc(model_final)$ICC_adjusted, error = function(e) NA_real_))

  # ── Mixed-sign: Gaussian (Yeo-Johnson / Box-Cox-neg if requested) ──────────
  } else if (var_min < 0 && var_max > 0) {
    cat("Decision: mixed-sign — Gaussian LME",
        if (use_yj) "+ Yeo-Johnson\n" else if (use_boxcox) "+ Box-Cox (bcnPower)\n"
        else ", no transform\n")
    fit_res <- fit_gaussian_aicc(df_var, v, fixed_str, re_cands, label = v,
                                 try_boxcox = use_boxcox, try_yeojohnson = use_yj)
    if (is.null(fit_res)) return(NULL)
    model_final   <- fit_res$model; model_type <- "Gaussian LME"
    .ass <- fit_res
    transform     <- fit_res$transform_name; shapiro_p <- shapiro_residual_p(model_final)
    norm_flag     <- fit_res$normality_flag
    re_chosen     <- fit_res$re_chosen; re_formula_str <- fit_res$re_formula
    anova_final   <- finalize_gaussian_inference(model_final, label = v)
    emm_final     <- tryCatch(emmeans(model_final, emm_spec), error = function(e) NULL)
    posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"), error = function(e) NULL)
    aicc_tbl      <- fit_res$aicc_table
    r2_g <- tryCatch(MuMIn::r.squaredGLMM(model_final), error = function(e) NULL)
    anova_effect_size_table(anova_final, model_final, model_type, label = v,
                            r2m = if (!is.null(r2_g)) unname(r2_g[1, "R2m"]) else NA_real_,
                            r2c = if (!is.null(r2_g)) unname(r2_g[1, "R2c"]) else NA_real_,
                            icc = tryCatch(performance::icc(model_final)$ICC_adjusted, error = function(e) NA_real_))

  # ── Positive continuous: Gaussian + transform check ────────────────────────
  } else if (var_min >= 0 && var_max > 1) {
    cat("Decision: positive continuous — Gaussian LME + transform check",
        if (use_yj) "(incl. Yeo-Johnson)\n" else if (use_boxcox) "(incl. Box-Cox)\n" else "\n")
    fit_res <- fit_gaussian_aicc(df_var, v, fixed_str, re_cands, label = v,
                                 try_boxcox = use_boxcox, try_yeojohnson = use_yj)
    if (is.null(fit_res)) return(NULL)
    model_final   <- fit_res$model; model_type <- "Gaussian LME"
    .ass <- fit_res
    transform     <- fit_res$transform_name; shapiro_p <- shapiro_residual_p(model_final)
    norm_flag     <- fit_res$normality_flag
    re_chosen     <- fit_res$re_chosen; re_formula_str <- fit_res$re_formula
    anova_final   <- finalize_gaussian_inference(model_final, label = v)
    emm_final     <- tryCatch(emmeans(model_final, emm_spec), error = function(e) NULL)
    posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"), error = function(e) NULL)
    aicc_tbl      <- fit_res$aicc_table
    r2_g <- tryCatch(MuMIn::r.squaredGLMM(model_final), error = function(e) NULL)
    anova_effect_size_table(anova_final, model_final, model_type, label = v,
                            r2m = if (!is.null(r2_g)) unname(r2_g[1, "R2m"]) else NA_real_,
                            r2c = if (!is.null(r2_g)) unname(r2_g[1, "R2c"]) else NA_real_,
                            icc = tryCatch(performance::icc(model_final)$ICC_adjusted, error = function(e) NA_real_))

  # ── Fallback Gaussian ──────────────────────────────────────────────────────
  } else {
    cat("Decision: fallback — Gaussian LME\n")
    fit_res <- fit_gaussian_aicc(df_var, v, fixed_str, re_cands, label = v,
                                 try_boxcox = use_boxcox, try_yeojohnson = use_yj)
    if (is.null(fit_res)) return(NULL)
    model_final   <- fit_res$model; model_type <- "Gaussian LME"
    .ass <- fit_res
    transform     <- fit_res$transform_name; shapiro_p <- shapiro_residual_p(model_final)
    norm_flag     <- fit_res$normality_flag
    re_chosen     <- fit_res$re_chosen; re_formula_str <- fit_res$re_formula
    anova_final   <- finalize_gaussian_inference(model_final, label = v)
    emm_final     <- tryCatch(emmeans(model_final, emm_spec), error = function(e) NULL)
    posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"), error = function(e) NULL)
    aicc_tbl      <- fit_res$aicc_table
    r2_g <- tryCatch(MuMIn::r.squaredGLMM(model_final), error = function(e) NULL)
    anova_effect_size_table(anova_final, model_final, model_type, label = v,
                            r2m = if (!is.null(r2_g)) unname(r2_g[1, "R2m"]) else NA_real_,
                            r2c = if (!is.null(r2_g)) unname(r2_g[1, "R2c"]) else NA_real_,
                            icc = tryCatch(performance::icc(model_final)$ICC_adjusted, error = function(e) NA_real_))
  }

  if (!is.null(model_final) && grepl("^Gaussian", model_type) &&
      !inherits(model_final, "merMod") && !inherits(model_final, "lmerModLmerTest"))
    model_type <- "Gaussian LM (no RE retained)"

  model_row <- tibble(
    variable           = v,
    min_value          = var_min,
    max_value          = var_max,
    model_type         = model_type,
    transformation     = transform,
    bc_lambda          = if (is.null(.ass$bc_lambda)) NA_real_ else .ass$bc_lambda,
    sw_p_raw           = if (is.null(.ass$sw_p_raw)) NA_real_ else .ass$sw_p_raw,
    levene_p_raw       = if (is.null(.ass$levene_p_raw)) NA_real_ else .ass$levene_p_raw,
    sw_p_trans         = if (is.null(.ass$sw_p_trans)) NA_real_ else .ass$sw_p_trans,
    levene_p_trans     = if (is.null(.ass$levene_p_trans)) NA_real_ else .ass$levene_p_trans,
    re_selected        = re_chosen,
    re_formula         = re_formula_str,
    residual_shapiro_p = shapiro_p,
    # assump_violated: TRUE = normality or homoscedasticity assumption VIOLATED
    # (SW p < 0.05 or Levene p < 0.05 on OLS residuals after best transform)
    # FALSE = assumptions satisfied; NA = not applicable (Beta/NB GLMM)
    assump_violated    = if (is.na(norm_flag)) NA_character_ else as.character(norm_flag)
  )

  cat("\nModel row:\n"); print(model_row)
  if (!is.null(anova_final)) { cat("\nANOVA:\n"); print(anova_final) }
  if (!is.null(posthoc_final)) { cat("\nTukey post hoc:\n"); print(posthoc_final) }

  list(
    model     = model_final,
    anova     = anova_final,
    emm       = emm_final,
    posthoc   = posthoc_final,
    aicc_tbl  = aicc_tbl,
    model_row = model_row
  )
}


################################################################################
###  INDICATOR DEFINITIONS
################################################################################

# Group-dynamics (school-level) indicators
GD_VARS      <- c("mean_iid_cm", "mean_nnd_cm", "mean_centroid_speed_cm_s",
                  "mean_hull_area_cm2", "mean_polarization", "zone_synchrony")
GD_FORCE_BETA  <- c("mean_polarization", "zone_synchrony")  # proportions 0-1
GD_COUNT_VARS  <- character(0)
GD_FORCE_GAUSS <- character(0)
# Box-Cox for school-level Mean NND at the interval level (maximise residual SW p).
GD_BOXCOX_INTERVAL <- c("mean_nnd_cm")

# Individual fish-level indicators
# time_spent_zone_flow_prop and time_moving_TS1_prop are bounded continuous
# variables (proportion of the recording window, [0,1]) → Beta GLMM, NOT Gaussian.
# Raw-second columns are retained in the dataset; the proportions are modelled.
# mean_nnd_cm (individual) is populated by joining per-fish NND (NND_fish_by_timepoint).
IND_VARS     <- c("distance_m", "speed_cm_s_mean", "accel_cm_s2_mean",
                  "n_zone_transitions",
                  "time_spent_zone_flow_prop", "time_moving_TS1_prop", "mean_nnd_cm")
IND_COUNT_VARS  <- c("n_zone_transitions")   # count → NB GLMM
IND_FORCE_BETA  <- c("time_spent_zone_flow_prop", "time_moving_TS1_prop")  # bounded [0,1]
IND_FORCE_GAUSS <- character(0)

# Transformations to minimise residual non-normality (Gaussian LME only):
#   Box-Cox (car::bcPower) for strictly positive responses (distance, speed, per-fish NND).
#   Yeo-Johnson (car::yjPower, ML lambda via car::powerTransform) for acceleration.
IND_BOXCOX_INTERVAL <- c("distance_m", "speed_cm_s_mean", "mean_nnd_cm")
IND_BOXCOX_TRIAL    <- character(0)
IND_YJ_INTERVAL     <- c("accel_cm_s2_mean")
IND_YJ_TRIAL        <- c("accel_cm_s2_mean")


################################################################################
###  ANALYSIS LOOPS
################################################################################

run_all <- function(vars, data, fixed_str, re_cands,
                    count_vars, force_beta, force_gauss, emm_spec,
                    boxcox_vars = character(0), yeojohnson_vars = character(0)) {
  results <- list()
  for (v in vars) {
    res <- tryCatch(
      run_indicator(v, data, fixed_str, re_cands,
                    count_vars, force_beta, force_gauss, emm_spec,
                    boxcox_vars = boxcox_vars,
                    yeojohnson_vars = yeojohnson_vars),
      error = function(e) { message("ERROR for ", v, ": ", conditionMessage(e)); NULL }
    )
    results[[v]] <- res
  }
  results
}

# ── 1. GD interval-level ─────────────────────────────────────────────────────
cat("\n\n################################################################\n")
cat("GD INTERVAL-LEVEL  |  fixed = Treatment * Group * interval\n")
cat("################################################################\n")
gd_int_res <- run_all(GD_VARS, gd_interval,
                      "Treatment * Group * interval", RE_GD_INTERVAL,
                      GD_COUNT_VARS, GD_FORCE_BETA, GD_FORCE_GAUSS,
                      ~ Treatment * Group * interval,
                      boxcox_vars = GD_BOXCOX_INTERVAL)

# ── 2. GD trial-aggregated ────────────────────────────────────────────────────
cat("\n\n################################################################\n")
cat("GD TRIAL-AGGREGATED  |  fixed = Treatment * Group\n")
cat("################################################################\n")
gd_tri_res <- run_all(GD_VARS, gd_trial,
                      "Treatment * Group", RE_GD_TRIAL,
                      GD_COUNT_VARS, GD_FORCE_BETA, GD_FORCE_GAUSS,
                      ~ Treatment * Group)

# ── 3. Individual interval-level ──────────────────────────────────────────────
cat("\n\n################################################################\n")
cat("INDIVIDUAL INTERVAL-LEVEL  |  fixed = Treatment * Group * interval\n")
cat("################################################################\n")
ind_int_res <- run_all(IND_VARS, ind_interval,
                       "Treatment * Group * interval", RE_IND_INTERVAL,
                       IND_COUNT_VARS, IND_FORCE_BETA, IND_FORCE_GAUSS,
                       ~ Treatment * Group * interval,
                       boxcox_vars = IND_BOXCOX_INTERVAL,
                       yeojohnson_vars = IND_YJ_INTERVAL)

# ── 4. Individual trial-aggregated ────────────────────────────────────────────
cat("\n\n################################################################\n")
cat("INDIVIDUAL TRIAL-AGGREGATED  |  fixed = Treatment * Group\n")
cat("################################################################\n")
ind_tri_res <- run_all(IND_VARS, ind_trial,
                       "Treatment * Group", RE_IND_TRIAL,
                       IND_COUNT_VARS, IND_FORCE_BETA, IND_FORCE_GAUSS,
                       ~ Treatment * Group,
                       boxcox_vars = IND_BOXCOX_TRIAL,
                       yeojohnson_vars = IND_YJ_TRIAL)


################################################################################
###  EXTRACT RESULT LISTS
################################################################################

extract_list <- function(results_list, element) {
  lapply(results_list, function(r) if (!is.null(r)) r[[element]] else NULL)
}

gd_int_model_dt  <- bind_rows(Filter(Negate(is.null), extract_list(gd_int_res,  "model_row")))
gd_tri_model_dt  <- bind_rows(Filter(Negate(is.null), extract_list(gd_tri_res,  "model_row")))
ind_int_model_dt <- bind_rows(Filter(Negate(is.null), extract_list(ind_int_res, "model_row")))
ind_tri_model_dt <- bind_rows(Filter(Negate(is.null), extract_list(ind_tri_res, "model_row")))
# Assumption-check table (single Box-Cox rule): one CSV with all four analysis levels.
write.csv(bind_rows(
  gd_int_model_dt  %>% mutate(unit = "Group dynamics",  level = "interval"),
  gd_tri_model_dt  %>% mutate(unit = "Group dynamics",  level = "trial"),
  ind_int_model_dt %>% mutate(unit = "Individual fish", level = "interval"),
  ind_tri_model_dt %>% mutate(unit = "Individual fish", level = "trial")),
  file.path(DIR_ANOVAS, "model_decisions_assumptions.csv"), row.names = FALSE)


################################################################################
###  OUTPUT HELPERS  (p floor: p < 0.001, handoff/CONVENTIONS.md C7/A-9)
################################################################################

pretty_indicator_labels <- c(
  mean_iid_cm              = "Mean IID (cm)",
  mean_nnd_cm              = "Mean NND (cm)",
  mean_centroid_speed_cm_s = "Mean centroid speed (cm/s)",
  mean_hull_area_cm2       = "Mean hull area (cm²)",
  mean_polarization        = "Polarisation",
  zone_synchrony           = "Zone synchrony",
  distance_m               = "Distance (m)",
  speed_cm_s_mean          = "Mean speed (cm/s)",
  accel_cm_s2_mean         = "Mean acceleration (cm/s²)",
  n_zone_transitions       = "Zone transitions (count)",
  time_spent_zone_flow_s   = "Time in flow zone (s)",
  time_moving_TS1_s        = "Time moving >1 BL/s (s)",
  time_spent_zone_flow_prop = "Flow-zone time (proportion)",
  time_moving_TS1_prop      = "Time moving >1 BL/s (proportion)"
)

pretty_variable_name <- function(x)
  ifelse(x %in% names(pretty_indicator_labels), pretty_indicator_labels[x], x)

# p-value floor is p < 0.001 (fmt_p_house(), house_style.R A5;
# CONVENTIONS.md C7/A-9 -- the old "< 0.0001" floor is superseded). This
# wrapper keeps the "p " prefix + fmt_p_house's "p = "/"p < 0.001" body
# consistent with the old call sites below, which prepend "p " themselves.
format_p_statement <- function(p)
  ifelse(is.na(p), "NA", sub("^p ", "", vapply(p, fmt_p_house, character(1))))
# to_subscript_statement() is sourced from house_style.R (A7).

# add_pct_change_column(): RESULTS_STYLE.md R5 [INV] -- "% change against the
# named reference group... mandatory in text whenever means are given". The
# GOT reference cell is Control:No-choice (TG_BREAKS_PLOT[1], the CLD anchor
# cell). % change is computed against the reference row holding every other
# grouping column (e.g. interval) fixed, using pct_change() from
# house_style.R.
add_pct_change_column <- function(emm_df) {
  val_col <- intersect(c("emmean", "response", "prob", "rate"), names(emm_df))
  if (!length(val_col) || !all(c("Treatment", "Group") %in% names(emm_df))) return(emm_df)
  val_col <- val_col[1]
  drop_cols  <- c("Treatment", "Group", val_col, "SE", "df", "asymp.LCL",
                 "asymp.UCL", "lower.CL", "upper.CL", "null")
  match_cols <- setdiff(names(emm_df), drop_cols)
  ref_rows <- emm_df[emm_df$Treatment == "Control" & emm_df$Group == "No choice", , drop = FALSE]
  emm_df$pct_change_vs_ControlNoChoice <- NA_real_
  if (!nrow(ref_rows)) return(emm_df)
  for (i in seq_len(nrow(emm_df))) {
    ref_i <- if (length(match_cols)) {
      which(apply(ref_rows[, match_cols, drop = FALSE], 1,
                  function(r) isTRUE(all(r == emm_df[i, match_cols]))))
    } else seq_len(nrow(ref_rows))
    if (length(ref_i) == 1)
      emm_df$pct_change_vs_ControlNoChoice[i] <- pct_change(emm_df[[val_col]][i], ref_rows[[val_col]][ref_i])
  }
  emm_df
}

format_numeric_columns <- function(tab, digits = 4)
  tab %>% mutate(across(where(is.numeric), ~ round(.x, digits)))

format_p_columns <- function(tab) {
  p_cols <- names(tab)[grepl("p|Pr\\(", names(tab), ignore.case = TRUE)]
  for (pc in p_cols) {
    if (is.numeric(tab[[pc]])) {
      tab[[pc]] <- ifelse(is.na(tab[[pc]]), NA, sub("^p = ", "", vapply(tab[[pc]], fmt_p_house, character(1))))
    }
  }
  tab
}

clean_table_for_markdown <- function(tab) {
  tab <- as.data.frame(tab)
  if ("variable" %in% names(tab)) tab$variable <- pretty_variable_name(tab$variable)
  # Rename columns for clarity in markdown output
  if ("normality_flag"  %in% names(tab))
    names(tab)[names(tab) == "normality_flag"]  <- "norm_violated"
  if ("assump_violated" %in% names(tab))
    names(tab)[names(tab) == "assump_violated"] <- "norm_violated"
  if ("residual_shapiro_p" %in% names(tab))
    names(tab)[names(tab) == "residual_shapiro_p"] <- "sw_p_LMM_resid"
  tab <- format_numeric_columns(tab)
  tab <- format_p_columns(tab)
  tab
}

make_markdown_table <- function(tab) {
  tab <- clean_table_for_markdown(tab)
  if (nrow(tab) == 0) return("_No rows._")
  tab[] <- lapply(tab, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  hdr  <- paste0("| ", paste(names(tab), collapse = " | "), " |")
  sep  <- paste0("| ", paste(rep("---", ncol(tab)), collapse = " | "), " |")
  body <- apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

anova_to_table <- function(anova_object) {
  tab <- as.data.frame(anova_object)
  tab$term <- rownames(tab); rownames(tab) <- NULL
  tab %>% relocate(term)
}

posthoc_to_table <- function(posthoc_object, significant_only = FALSE) {
  tab <- as.data.frame(posthoc_object)
  if (significant_only && "p.value" %in% names(tab))
    tab <- tab %>% filter(p.value < 0.05)
  tab
}

write_model_markdown_report <- function(
    file_path, title, analysis_description,
    model_decision_table_object, anova_list_object,
    emmeans_list_object, posthoc_list_object,
    aicc_list_object = NULL,
    companion_table = NULL, companion_label = NULL) {

  md <- c(
    paste0("# ", title), "",
    analysis_description, "",
    paste0("Generated on: ", Sys.time()), "",
    "---", "",
    "## Model-selection summary", "",
    make_markdown_table(model_decision_table_object), ""
  )

  if (!is.null(companion_table)) {
    lbl <- if (!is.null(companion_label)) companion_label else "Companion analysis"
    md <- c(md, paste0("### ", lbl), "", make_markdown_table(companion_table), "")
  }

  for (v in names(anova_list_object)) {
    if (is.null(anova_list_object[[v]])) next
    pretty_v <- pretty_variable_name(sub("^(GD_|Ind_)", "", v))
    unit_tag  <- if (grepl("^GD_", v)) " [Group dynamics]" else " [Individual]"
    md <- c(md, "", "---", "",
            paste0("## ", pretty_v, unit_tag), "",
            paste0("**Variable key:** `", v, "`"), "")

    if (!is.null(aicc_list_object) && !is.null(aicc_list_object[[v]]))
      md <- c(md, "### AICc RE selection", "",
              make_markdown_table(aicc_list_object[[v]]), "")

    row_v <- sub("^(GD_|Ind_)", "", v)
    if ("variable" %in% names(model_decision_table_object)) {
      dr <- model_decision_table_object %>% filter(variable == row_v)
      if (nrow(dr) > 0)
        md <- c(md, "### Selected model", "", make_markdown_table(dr), "")
    }

    md <- c(md, "### Type III ANOVA", "",
            make_markdown_table(anova_to_table(anova_list_object[[v]])), "")

    if (!is.null(emmeans_list_object[[v]]))
      md <- c(md, "### Estimated marginal means", "",
              make_markdown_table(add_pct_change_column(as.data.frame(emmeans_list_object[[v]]))), "")

    if (!is.null(posthoc_list_object[[v]])) {
      sig_ph <- posthoc_to_table(posthoc_list_object[[v]], significant_only = TRUE)
      md <- c(md, "### Significant Tukey contrasts (p < 0.05)", "",
              make_markdown_table(sig_ph), "")
      all_ph <- posthoc_to_table(posthoc_list_object[[v]], significant_only = FALSE)
      md <- c(md, "### All Tukey contrasts", "",
              make_markdown_table(all_ph), "")
    }
  }

  writeLines(md, con = file_path, useBytes = TRUE)
  normalizePath(file_path)
}

# ── Word report mirroring the markdown report ────────────────────────────────
# Renders the SAME content as write_model_markdown_report (model-selection summary,
# optional companion table, and per-indicator AICc / selected model / Type III
# ANOVA / EMMs / Tukey contrasts) into a .docx via officer + flextable.
.df_to_flex <- function(df, font_size = 8L) {
  df <- clean_table_for_markdown(df)          # same prettifying as the markdown tables
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  if (nrow(df) == 0) df <- data.frame(` ` = "No rows.", check.names = FALSE)
  ft <- flextable(df)
  ft <- bold(ft, part = "header")
  ft <- bg(ft, part = "header", bg = "#2E75B6")
  ft <- color(ft, part = "header", color = "white")
  ft <- fontsize(ft, part = "all", size = font_size)
  ft <- font(ft, part = "all", fontname = "Calibri")
  ft <- align(ft, part = "header", align = "center")
  ft <- align(ft, part = "body",   align = "left")
  ft <- set_table_properties(ft, width = 1, layout = "autofit")
  ft
}

write_model_word_report <- function(
    file_path, title, analysis_description,
    model_decision_table_object, anova_list_object,
    emmeans_list_object, posthoc_list_object,
    aicc_list_object = NULL,
    companion_table = NULL, companion_label = NULL) {

  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, analysis_description, style = "Normal")
  doc <- body_add_par(doc, paste0("Generated on: ", Sys.time()), style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")

  doc <- body_add_par(doc, "Model-selection summary", style = "heading 2")
  doc <- body_add_flextable(doc, .df_to_flex(model_decision_table_object))

  if (!is.null(companion_table)) {
    lbl <- if (!is.null(companion_label)) companion_label else "Companion analysis"
    doc <- body_add_par(doc, lbl, style = "heading 3")
    doc <- body_add_flextable(doc, .df_to_flex(companion_table))
  }

  for (v in names(anova_list_object)) {
    if (is.null(anova_list_object[[v]])) next
    pretty_v <- pretty_variable_name(sub("^(GD_|Ind_)", "", v))
    unit_tag <- if (grepl("^GD_", v)) " [Group dynamics]" else " [Individual]"
    doc <- body_add_par(doc, paste0(pretty_v, unit_tag), style = "heading 2")
    doc <- body_add_par(doc, paste0("Variable key: ", v), style = "Normal")

    if (!is.null(aicc_list_object) && !is.null(aicc_list_object[[v]])) {
      doc <- body_add_par(doc, "AICc RE selection", style = "heading 3")
      doc <- body_add_flextable(doc, .df_to_flex(aicc_list_object[[v]]))
    }

    row_v <- sub("^(GD_|Ind_)", "", v)
    if ("variable" %in% names(model_decision_table_object)) {
      dr <- model_decision_table_object %>% filter(variable == row_v)
      if (nrow(dr) > 0) {
        doc <- body_add_par(doc, "Selected model", style = "heading 3")
        doc <- body_add_flextable(doc, .df_to_flex(dr))
      }
    }

    doc <- body_add_par(doc, "Type III ANOVA", style = "heading 3")
    doc <- body_add_flextable(doc, .df_to_flex(anova_to_table(anova_list_object[[v]])))

    if (!is.null(emmeans_list_object[[v]])) {
      doc <- body_add_par(doc, "Estimated marginal means", style = "heading 3")
      doc <- body_add_flextable(doc, .df_to_flex(add_pct_change_column(as.data.frame(emmeans_list_object[[v]]))))
    }

    if (!is.null(posthoc_list_object[[v]])) {
      doc <- body_add_par(doc, "Significant Tukey contrasts (p < 0.05)", style = "heading 3")
      doc <- body_add_flextable(doc,
        .df_to_flex(posthoc_to_table(posthoc_list_object[[v]], significant_only = TRUE)))
      doc <- body_add_par(doc, "All Tukey contrasts", style = "heading 3")
      doc <- body_add_flextable(doc,
        .df_to_flex(posthoc_to_table(posthoc_list_object[[v]], significant_only = FALSE)))
    }
  }

  print(doc, target = file_path)
  normalizePath(file_path)
}

format_anova_statement <- function(anova_obj, model_obj = NULL,
                                   indicator, term, term_label) {
  if (is.null(anova_obj)) return(NA_character_)
  tab <- as.data.frame(anova_obj); tab$term <- rownames(tab)
  row <- tab[tab$term == term, , drop = FALSE]
  if (nrow(row) == 0) return(NA_character_)

  # Delegates the statistic/p/effect-size string to anova_statement()
  # (house_style.R A7) -- this is what makes "every F carries eta2p,
  # every chisq carries none" structural rather than a rule to remember
  # (HANDOFF.md S6 [INV], CONVENTIONS.md C6).
  if (all(c("NumDF", "DenDF", "F value", "Pr(>F)") %in% names(row))) {
    anova_statement(term_label, "F", stat = row$`F value`,
                    df1 = row$NumDF, df2 = row$DenDF, p = row$`Pr(>F)`)
  } else if (all(c("Df", "F value", "Pr(>F)") %in% names(row))) {
    df2 <- if (!is.null(model_obj)) df.residual(model_obj) else NA_real_
    anova_statement(term_label, "F", stat = row$`F value`,
                    df1 = row$Df, df2 = df2, p = row$`Pr(>F)`)
  } else if (all(c("Chisq", "Df", "Pr(>Chisq)") %in% names(row))) {
    anova_statement(term_label, "chisq", stat = row$Chisq,
                    df1 = row$Df, p = row$`Pr(>Chisq)`)
  } else NA_character_
}

make_statement_table <- function(anova_list_obj, model_list_obj, term, term_label) {
  indicators <- names(anova_list_obj)
  tibble(
    variable  = indicators,
    indicator = pretty_variable_name(sub("^(GD_|Ind_)", "", indicators)),
    unit      = ifelse(grepl("^GD_",  indicators), "Group dynamics", "Individual"),
    statement = map_chr(indicators, function(v)
      tryCatch(
        format_anova_statement(anova_list_obj[[v]], model_list_obj[[v]],
                               v, term, term_label),
        error = function(e) NA_character_
      )
    )
  )
}

write_csv_utf8bom <- function(df, path) {
  tmp <- tempfile(fileext = ".csv")
  write.csv(df, tmp, row.names = FALSE, fileEncoding = "UTF-8")
  raw_content <- readBin(tmp, raw(), file.info(tmp)$size)
  writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), raw_content), path)
  unlink(tmp)
}

pretty_decision_cols <- c(
  variable           = "Variable",
  min_value          = "Min",
  max_value          = "Max",
  model_type         = "Model type",
  transformation     = "Transformation",
  re_selected        = "RE selected",
  re_formula         = "RE formula",
  residual_shapiro_p = "SW p (LMM resid.)",    # p ≥ 0.05 → assumption satisfied
  assump_violated    = "Norm. violated?"        # TRUE = SW or Levene p < 0.05 (OLS screen)
)

format_decision_table_for_word <- function(dt) {
  dt2 <- as.data.frame(dt)
  dt2$variable <- pretty_variable_name(dt2$variable)
  if ("min_value" %in% names(dt2)) dt2$min_value <- round(dt2$min_value, 3)
  if ("max_value" %in% names(dt2)) dt2$max_value <- round(dt2$max_value, 3)
  if ("residual_shapiro_p" %in% names(dt2)) {
    pv <- suppressWarnings(as.numeric(dt2$residual_shapiro_p))
    dt2$residual_shapiro_p <- ifelse(is.na(pv), "", vapply(pv, fmt_p_house, character(1)))
  }
  for (old_nm in names(pretty_decision_cols)) {
    new_nm <- pretty_decision_cols[[old_nm]]
    idx    <- which(names(dt2) == old_nm)
    if (length(idx) == 1L) names(dt2)[idx] <- new_nm
  }
  dt2
}

make_decision_flextable <- function(dt) {
  formatted <- format_decision_table_for_word(dt)
  ft <- flextable(formatted)
  ft <- bold(ft, part = "header")
  ft <- bg(ft, part = "header", bg = "#2E75B6")
  ft <- color(ft, part = "header", color = "white")
  ft <- fontsize(ft, part = "all", size = 10)
  ft <- font(ft, part = "all", fontname = "Calibri")
  ft <- align(ft, part = "header", align = "center")
  ft <- align(ft, part = "body",   align = "left")
  ft <- set_table_properties(ft, width = 1, layout = "autofit")
  ft
}


################################################################################
###  BUILD MERGED RESULT OBJECTS FOR REPORTING
################################################################################

# Add unit prefix to names so GD and individual results don't collide
prefix_results <- function(results, prefix) {
  stats::setNames(results, paste0(prefix, names(results)))
}

gd_int_res_p  <- prefix_results(gd_int_res,  "GD_")
gd_tri_res_p  <- prefix_results(gd_tri_res,  "GD_")
ind_int_res_p <- prefix_results(ind_int_res, "Ind_")
ind_tri_res_p <- prefix_results(ind_tri_res, "Ind_")

# Combined decision tables (with unit column)
gd_int_dt_lab  <- gd_int_model_dt  %>% mutate(unit = "Group dynamics")
ind_int_dt_lab <- ind_int_model_dt %>% mutate(unit = "Individual fish")
interval_combined_dt <- bind_rows(gd_int_dt_lab, ind_int_dt_lab) %>% relocate(unit)

gd_tri_dt_lab  <- gd_tri_model_dt  %>% mutate(unit = "Group dynamics")
ind_tri_dt_lab <- ind_tri_model_dt %>% mutate(unit = "Individual fish")
trial_combined_dt <- bind_rows(gd_tri_dt_lab, ind_tri_dt_lab) %>% relocate(unit)

# Merged ANOVA / EMM / posthoc / AICc lists (both units)
all_int <- c(gd_int_res_p, ind_int_res_p)
all_tri <- c(gd_tri_res_p, ind_tri_res_p)


################################################################################
###  EXPORT MARKDOWN + WORD REPORTS  (→ reports/ subfolder)
################################################################################

int_desc <- paste(
  "Fixed effects: Treatment × Group × interval.",
  "RE: AICc-selected from trial_id / trial_date / trial_id+trial_date.",
  "GD indicators: session-level school metrics (TP2+TP3, ExeChoice+NoChoice).",
  "Individual indicators: per-fish metrics (TP2+TP3, ExeChoice+NoChoice).",
  "Transform screening: SW + Levene on OLS residuals; log1p/sqrt/Box-Cox/Yeo-Johnson as configured.",
  "Column key — sw_p_LMM_resid: SW p on final LMM residuals (p >= 0.05 = normality satisfied;",
  "H0: residuals are normally distributed).",
  "norm_violated: OLS-screening flag (TRUE = assumption issue detected before/after transform;",
  "FALSE = screening passed; NA = non-Gaussian family)."
)
tri_desc <- paste(
  "Fixed effects: Treatment × Group.",
  "RE: AICc-selected (trial_date only at trial level).",
  "GD: trial-mean of session-level school metrics (N ≈ 12 trials).",
  "Individual: trial-mean over fish × intervals (N ≈ 12 trials).",
  "Transform screening: SW + Levene on OLS residuals; log1p/sqrt/Box-Cox/Yeo-Johnson as configured.",
  "Column key — sw_p_LMM_resid: SW p on final LMM residuals (p >= 0.05 = normality satisfied).",
  "norm_violated: OLS-screening flag (TRUE = assumption violated; FALSE = passed; NA = GLMM)."
)

int_md_path  <- file.path(DIR_REPORTS, "GOT_interval_Treatment_x_Group_x_interval.md")
int_docx_path <- file.path(DIR_REPORTS, "GOT_interval_Treatment_x_Group_x_interval.docx")
tri_md_path  <- file.path(DIR_REPORTS, "GOT_trial_aggregated_Treatment_x_Group.md")
tri_docx_path <- file.path(DIR_REPORTS, "GOT_trial_aggregated_Treatment_x_Group.docx")

# Interval-level — markdown then analogous Word document
write_model_markdown_report(
  file_path = int_md_path,
  title     = "GOT interval-level model results (AICc RE selection)",
  analysis_description = int_desc,
  model_decision_table_object = interval_combined_dt,
  anova_list_object   = extract_list(all_int, "anova"),
  emmeans_list_object = extract_list(all_int, "emm"),
  posthoc_list_object = extract_list(all_int, "posthoc"),
  aicc_list_object    = extract_list(all_int, "aicc_tbl")
)
write_model_word_report(
  file_path = int_docx_path,
  title     = "GOT interval-level model results (AICc RE selection)",
  analysis_description = int_desc,
  model_decision_table_object = interval_combined_dt,
  anova_list_object   = extract_list(all_int, "anova"),
  emmeans_list_object = extract_list(all_int, "emm"),
  posthoc_list_object = extract_list(all_int, "posthoc"),
  aicc_list_object    = extract_list(all_int, "aicc_tbl")
)
cat("\nInterval-level reports:\n  ", int_md_path, "\n  ", int_docx_path, "\n")

# Trial-aggregated — markdown then analogous Word document
write_model_markdown_report(
  file_path = tri_md_path,
  title     = "GOT trial-aggregated model results (AICc RE selection)",
  analysis_description = tri_desc,
  model_decision_table_object = trial_combined_dt,
  anova_list_object   = extract_list(all_tri, "anova"),
  emmeans_list_object = extract_list(all_tri, "emm"),
  posthoc_list_object = extract_list(all_tri, "posthoc"),
  aicc_list_object    = extract_list(all_tri, "aicc_tbl"),
  companion_table     = interval_combined_dt,
  companion_label     = "Interval-level model-selection summary (for reference)"
)
write_model_word_report(
  file_path = tri_docx_path,
  title     = "GOT trial-aggregated model results (AICc RE selection)",
  analysis_description = tri_desc,
  model_decision_table_object = trial_combined_dt,
  anova_list_object   = extract_list(all_tri, "anova"),
  emmeans_list_object = extract_list(all_tri, "emm"),
  posthoc_list_object = extract_list(all_tri, "posthoc"),
  aicc_list_object    = extract_list(all_tri, "aicc_tbl"),
  companion_table     = interval_combined_dt,
  companion_label     = "Interval-level model-selection summary (for reference)"
)
cat("Trial-aggregated reports:\n  ", tri_md_path, "\n  ", tri_docx_path, "\n")


################################################################################
###  EXPORT ANOVA STATEMENT CSVs (UTF-8 BOM)
################################################################################

all_int_anova <- extract_list(all_int, "anova")
all_tri_anova <- extract_list(all_tri, "anova")
all_int_model <- extract_list(all_int, "model")
all_tri_model <- extract_list(all_tri, "model")

# Treatment × Group interaction (main effect of interest at both levels)
stmt_int_TxG   <- make_statement_table(all_int_anova, all_int_model,
                                        "Treatment:Group", "Treatment × Group")
stmt_tri_TxG   <- make_statement_table(all_tri_anova, all_tri_model,
                                        "Treatment:Group", "Treatment × Group")
# Three-way interaction at interval level
stmt_int_TxGxI <- make_statement_table(all_int_anova, all_int_model,
                                        "Treatment:Group:interval",
                                        "Treatment × Group × interval")
# Main effects
stmt_int_T     <- make_statement_table(all_int_anova, all_int_model,
                                        "Treatment", "Treatment")
stmt_int_G     <- make_statement_table(all_int_anova, all_int_model,
                                        "Group", "Group")
# Trial-level main effects. Added 2026-08-21: these were previously produced
# only inside the markdown report, so figure_3_locomotion.R's trial-level caption
# statistics had to be hardcoded (and silently went stale). Exporting them here
# lets those captions be read from disk like every other statistic.
stmt_tri_T     <- make_statement_table(all_tri_anova, all_tri_model,
                                        "Treatment", "Treatment")
stmt_tri_G     <- make_statement_table(all_tri_anova, all_tri_model,
                                        "Group", "Group")

################################################################################
###  C12 SENSITIVITY ARTIFACTS  (HANDOFF.md S7, CONVENTIONS.md C12)
###  Written to DIR_ANOVAS alongside the anova_stmt_*.csv exports below.
###  Primary endpoint (plan Open question 1): time_spent_zone_flow_prop at the
###  trial level -- GOT's structural analogue of CHOICE's flow-preference
###  endpoint (a bounded proportion of trial time in the flow zone).
################################################################################

PRIMARY_ENDPOINT <- "time_spent_zone_flow_prop"

# -- aicc_selection.csv, anova_satterthwaite.csv, anova.csv: flush the
#    accumulators populated during the indicator loops above. --------------
if (length(.AICC_LOG))
  write.csv(bind_rows(.AICC_LOG), file.path(DIR_ANOVAS, "aicc_selection.csv"), row.names = FALSE)

# re_reduction.csv -- which design-motivated random effects the data did not
# identify, and were therefore dropped. Written whether or not anything was
# dropped, so an empty file means "nothing was singular", not "nobody looked".
re_reduction_log_write(DIR_ANOVAS)
if (length(.SATT_LOG))
  write.csv(bind_rows(.SATT_LOG), file.path(DIR_ANOVAS, "anova_satterthwaite.csv"), row.names = FALSE)
if (length(.ANOVA_CSV_LOG))
  write.csv(bind_rows(.ANOVA_CSV_LOG), file.path(DIR_ANOVAS, "anova.csv"), row.names = FALSE)
cat("  wrote aicc_selection.csv / anova_satterthwaite.csv / anova.csv\n")

# -- design_balance.csv: run over the factor set actually used as RE
#    candidates, on both analysis units. This is where GOT's trial_date <->
#    tank aliasing (LOGS.md 2026-05) must surface -- see the plan's Open
#    question 4; it is reported here, not silently assumed away. -----------
db_gd  <- sensitivity_design_balance(gd_interval,  c("trial_id", "trial_date", "Treatment", "Group"),
                                     label = "GD_interval", run_dir = DIR_ANOVAS)
db_ind <- sensitivity_design_balance(ind_interval, c("trial_id", "trial_date", "Treatment", "Group"),
                                     label = "Ind_interval", run_dir = DIR_ANOVAS)
write.csv(bind_rows(db_gd, db_ind), file.path(DIR_ANOVAS, "design_balance.csv"), row.names = FALSE)
cat("  wrote design_balance.csv\n")

# -- metric_collinearity.csv: every other individual indicator screened
#    (Spearman) against the primary endpoint, at the tier that is REPORTED
#    (fish x interval), not one level up. ------------------------------------
.collin_candidates <- setdiff(IND_VARS, PRIMARY_ENDPOINT)
if (PRIMARY_ENDPOINT %in% names(ind_interval))
  metric_collinearity(ind_interval, PRIMARY_ENDPOINT, .collin_candidates, run_dir = DIR_ANOVAS)

# -- sensitivity_loo.csv: leave-one-TRIAL-out stability of the primary
#    endpoint's fixed effects. The model is the fish x interval fit; the
#    leave-one-out unit is the trial (12 of them), not the row. --------------
.primary_model     <- ind_int_res[[PRIMARY_ENDPOINT]]$model      # fish x interval (primary)
.primary_model_agg <- ind_tri_res[[PRIMARY_ENDPOINT]]$model      # trial (comparator)
if (!is.null(.primary_model))
  sensitivity_loo(.primary_model, ind_interval, trial_col = "trial_id",
                  label = PRIMARY_ENDPOINT, run_dir = DIR_ANOVAS)

################################################################################
# sensitivity_forced_re.csv -- HANDOFF S7 forced-RE check.
#
# The primary RE structure is now design-determined and cannot be varied by
# selection (see the RANDOM-EFFECT STRUCTURE block above), so the robustness
# question moves here: what would the conclusion have been under each
# structure the design argument rules out? Three declared comparators per
# outcome, fixed effects untouched throughout:
#
#   (i)   date_dropped  -- (1|trial_date) removed. This is what the previous
#         vintage reported: Group at df2 ~ 8 rather than ~ 4.
#   (ii)  fish_dropped  -- (1|fish_id_matched_trial) removed. Interval terms
#         at df2 ~ 80 rather than ~ 40.
#   (iii) splitplot_df  -- the DESIGN-BASED denominator df (Group and Treatment
#         on 4 df: 6 dates - 2 groups; 12 trials - 6 dates - 2 terms),
#         recomputed from the primary F regardless of the boundary variance
#         estimate. With 6 dates and 3 per arm a REML estimate of zero is what
#         the design produces whether or not the true variance is zero, so the
#         design-based df ship alongside whatever Kenward-Roger returns.
################################################################################
cat("\n-- forced-RE sensitivity (three declared comparators) --\n")

.SPLITPLOT_DF <- c(Treatment = 4, Group = 4)
.SPLITPLOT_DF["Treatment:Group"] <- 4

.pval_from <- function(an, trm) {
  if (is.null(an)) return(NA_real_)
  if (!trm %in% rownames(an)) return(NA_real_)
  pcol <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p"), names(an))
  if (!length(pcol)) return(NA_real_)
  as.numeric(an[trm, pcol[1]])
}
.icc_of <- function(m) tryCatch(performance::icc(m)$ICC_adjusted, error = function(e) NA_real_)
.sing_of <- function(m) tryCatch(
  if (inherits(m, "glmmTMB")) NA else lme4::isSingular(m), error = function(e) NA)

# The comparator is refitted from the PRIMARY MODEL'S OWN frame and its own
# fixed formula, so the response scale, the transform already applied, and the
# row exclusions are identical by construction -- only the random part moves.
# HANDOFF S7: "never replace, never touch the fixed effects."
.model_frame_of <- function(m) if (inherits(m, "glmmTMB")) m$frame else m@frame
.fixed_formula_of <- function(m)
  paste(deparse(lme4::nobars(stats::formula(m))), collapse = " ")

.forced_re_rows <- list()
.forced_re_one <- function(label, res_primary, comparators, is_beta = FALSE) {
  m0 <- res_primary$model
  if (is.null(m0)) return(invisible(NULL))
  d_fit <- tryCatch(.model_frame_of(m0), error = function(e) NULL)
  fixed_lhs_rhs <- tryCatch(.fixed_formula_of(m0), error = function(e) NULL)
  if (is.null(d_fit) || is.null(fixed_lhs_rhs)) return(invisible(NULL))
  re_primary_str <- paste(
    vapply(lme4::findbars(stats::formula(m0)),
           function(b) paste0("(", paste(deparse(b), collapse = " "), ")"), character(1)),
    collapse = " + ")
  an0 <- tryCatch(as.data.frame(car::Anova(m0, type = 3)), error = function(e) NULL)
  for (cn in names(comparators)) {
    fmla <- tryCatch(as.formula(paste(fixed_lhs_rhs, "+", comparators[[cn]])),
                     error = function(e) NULL)
    if (is.null(fmla)) next
    m1 <- tryCatch(
      if (is_beta) glmmTMB(fmla, data = d_fit, family = beta_family(link = "logit"))
      else lmerTest::lmer(fmla, data = d_fit, REML = TRUE),
      error = function(e) NULL)
    if (is.null(m1)) next
    an1 <- tryCatch(as.data.frame(car::Anova(m1, type = 3)), error = function(e) NULL)
    if (is.null(an0) || is.null(an1)) next
    for (trm in intersect(rownames(an0), rownames(an1))) {
      if (identical(trm, "(Intercept)")) next
      p0 <- .pval_from(an0, trm); p1 <- .pval_from(an1, trm)
      .forced_re_rows[[length(.forced_re_rows) + 1]] <<- data.frame(
        label = label, term = trm, comparator = cn,
        re_primary    = re_primary_str,
        re_comparator = comparators[[cn]],
        p_primary = p0, p_comparator = p1,
        icc_primary = .icc_of(m0), icc_comparator = .icc_of(m1),
        singular_primary = .sing_of(m0), singular_comparator = .sing_of(m1),
        p_crosses_alpha = isTRUE((p0 < 0.05) != (p1 < 0.05)),
        stringsAsFactors = FALSE)
    }
  }
  invisible(NULL)
}

for (.lab in names(ind_int_res)) {
  .r <- ind_int_res[[.lab]]
  if (is.null(.r)) next
  tryCatch(.forced_re_one(.lab, .r, RE_COMPARATORS_IND_INTERVAL,
                          is_beta = .lab %in% IND_FORCE_BETA),
           error = function(e) message("forced-RE failed for ", .lab, ": ", conditionMessage(e)))
}
for (.lab in names(gd_int_res)) {
  .r <- gd_int_res[[.lab]]
  if (is.null(.r)) next
  tryCatch(.forced_re_one(.lab, .r, RE_COMPARATORS_GD_INTERVAL,
                          is_beta = .lab %in% GD_FORCE_BETA),
           error = function(e) message("forced-RE failed for ", .lab, ": ", conditionMessage(e)))
}

# Comparator (iii): design-based split-plot df, recomputed from the primary F.
.splitplot_rows <- list()
.anova_all <- if (length(.ANOVA_CSV_LOG)) bind_rows(.ANOVA_CSV_LOG) else NULL
if (!is.null(.anova_all)) {
  .sp <- .anova_all[.anova_all$stat_type == "F" &
                    .anova_all$term %in% names(.SPLITPLOT_DF), , drop = FALSE]
  for (.i in seq_len(nrow(.sp))) {
    .Fv <- .sp$statistic[.i]; .d1 <- .sp$df1[.i]
    .d2sp <- unname(.SPLITPLOT_DF[.sp$term[.i]])
    .p_sp <- stats::pf(.Fv, .d1, .d2sp, lower.tail = FALSE)
    .splitplot_rows[[length(.splitplot_rows) + 1]] <- data.frame(
      label = .sp$label[.i], term = .sp$term[.i], comparator = "splitplot_df",
      re_primary = "forced design RE",
      re_comparator = "design-based df (6 dates - 2 groups; 12 trials - 6 dates - 2 terms)",
      p_primary = .sp$p[.i], p_comparator = .p_sp,
      icc_primary = .sp$ICC[.i], icc_comparator = NA_real_,
      singular_primary = NA, singular_comparator = NA,
      p_crosses_alpha = isTRUE((.sp$p[.i] < 0.05) != (.p_sp < 0.05)),
      stringsAsFactors = FALSE)
  }
}

.forced_re_all <- bind_rows(c(.forced_re_rows, .splitplot_rows))
if (!length(.forced_re_all) || !nrow(.forced_re_all))
  .forced_re_all <- data.frame(label = NA_character_, term = NA_character_,
                               comparator = "none computable",
                               re_primary = NA_character_, re_comparator = NA_character_,
                               p_primary = NA_real_, p_comparator = NA_real_,
                               icc_primary = NA_real_, icc_comparator = NA_real_,
                               singular_primary = NA, singular_comparator = NA,
                               p_crosses_alpha = NA)
write.csv(.forced_re_all, file.path(DIR_ANOVAS, "sensitivity_forced_re.csv"), row.names = FALSE)
cat("  wrote sensitivity_forced_re.csv (", nrow(.forced_re_all), " rows, ",
    sum(.forced_re_all$p_crosses_alpha, na.rm = TRUE), " crossing alpha)\n", sep = "")

################################################################################
# pseudoreplication_check.csv -- HANDOFF S7 aggregation-level check.
#
# "Refit at the unit level AND the level above it; report both beta, both p,
# the beta ratio, and a direction-agree flag. A conclusion that flips with
# aggregation level is NOT reportable as a treatment effect."
#
# Two comparisons per individual indicator, since GOT has three tiers:
#   fish x interval (primary)  vs  trial x interval
#   fish x interval (primary)  vs  trial-aggregated
################################################################################
cat("\n-- aggregation-level (pseudoreplication) check --\n")

# The trial x interval tier is built HERE, for the comparison only. It is NOT a
# reported model: aggregating to it would collapse the interval terms from
# ~80 df to ~8 and discard genuine within-trial replication (see the
# RANDOM-EFFECT STRUCTURE block). It exists solely as the level-above
# comparator this check requires.
ind_trial_interval <- ind_interval %>%
  group_by(trial_id, trial_date, Treatment, Group, interval) %>%
  summarise(across(all_of(intersect(IND_VARS, names(ind_interval))),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

.beta_p_for <- function(model, term) {
  an <- tryCatch(as.data.frame(car::Anova(model, type = 3)), error = function(e) NULL)
  p <- .pval_from(an, term)
  b <- tryCatch({
    fx <- if (inherits(model, "glmmTMB")) glmmTMB::fixef(model)$cond else lme4::fixef(model)
    unname(fx[grep(paste0("^", term), names(fx))][1])
  }, error = function(e) NA_real_)
  c(as.numeric(b), p)
}

.pseudo_rows <- list()
for (.v in intersect(IND_VARS, names(ind_trial_interval))) {
  .m_unit <- ind_int_res[[.v]]$model
  .m_ti <- tryCatch({
    d <- as.data.frame(ind_trial_interval)
    d <- d[is.finite(d[[.v]]), , drop = FALSE]
    d$.y <- if (.v %in% IND_FORCE_BETA) adjust_beta_01(d[[.v]]) else d[[.v]]
    if (.v %in% IND_FORCE_BETA)
      glmmTMB(.y ~ Treatment * Group * interval + (1 | trial_date) + (1 | trial_id),
              data = d, family = beta_family(link = "logit"))
    else
      lmerTest::lmer(.y ~ Treatment * Group * interval + (1 | trial_date) + (1 | trial_id),
                     data = d, REML = TRUE)
  }, error = function(e) NULL)
  .m_agg <- ind_tri_res[[.v]]$model
  .pairs <- list(list(nm = "fish_interval_vs_trial_interval", m = .m_ti),
                 list(nm = "fish_interval_vs_trial_aggregated", m = .m_agg))
  for (.pp in .pairs) {
    if (is.null(.m_unit) || is.null(.pp$m)) next
    u <- .beta_p_for(.m_unit, "Treatment"); a <- .beta_p_for(.pp$m, "Treatment")
    .pseudo_rows[[length(.pseudo_rows) + 1]] <- data.frame(
      label = .v, term = "Treatment", comparison = .pp$nm,
      beta_unit_level = u[1], p_unit_level = u[2],
      beta_agg_level = a[1], p_agg_level = a[2],
      beta_ratio = if (is.finite(a[1]) && a[1] != 0) u[1] / a[1] else NA_real_,
      direction_agree = isTRUE(sign(u[1]) == sign(a[1])),
      reportable = isTRUE(sign(u[1]) == sign(a[1])) &&
                   isTRUE((u[2] < 0.05) == (a[2] < 0.05)),
      stringsAsFactors = FALSE)
  }
}
.pseudo_all <- if (length(.pseudo_rows)) bind_rows(.pseudo_rows) else
  data.frame(label = NA_character_, term = NA_character_, comparison = "not computable",
             beta_unit_level = NA_real_, p_unit_level = NA_real_,
             beta_agg_level = NA_real_, p_agg_level = NA_real_,
             beta_ratio = NA_real_, direction_agree = NA, reportable = NA)
write.csv(.pseudo_all, file.path(DIR_ANOVAS, "pseudoreplication_check.csv"), row.names = FALSE)
cat("  wrote pseudoreplication_check.csv (", nrow(.pseudo_all), " rows, ",
    sum(!.pseudo_all$reportable, na.rm = TRUE), " non-reportable)\n", sep = "")

################################################################################
# sensitivity_identity_tier.csv -- HANDOFF S4.1 [INV].
#
# A fish random intercept spanning intervals is an IDENTITY-DEPENDENT modelling
# choice: it asserts that fish_id_matched_trial refers to the same physical
# animal at TP2 and TP3. Identity is validated (the matcher writes a per-fish
# reliability tier), so the choice is licensed -- but the licence is CHECKED
# here, not asserted, by refitting the primary without the rows the matcher
# itself flagged 'unreliable'.
################################################################################
.id_rows <- list()
if ("id_reliability_tier" %in% names(ind_interval)) {
  .ind_rel <- ind_interval[ind_interval$id_reliability_tier != "unreliable", , drop = FALSE]
  cat("\n-- identity-tier sensitivity: dropping ",
      nrow(ind_interval) - nrow(.ind_rel), " 'unreliable' rows of ",
      nrow(ind_interval), " --\n", sep = "")
  for (.v in intersect(IND_VARS, names(.ind_rel))) {
    .m_full <- ind_int_res[[.v]]$model
    .m_rel <- tryCatch({
      d <- as.data.frame(.ind_rel)
      d <- d[is.finite(d[[.v]]), , drop = FALSE]
      d$.y <- if (.v %in% IND_FORCE_BETA) adjust_beta_01(d[[.v]]) else d[[.v]]
      f <- as.formula(paste(".y ~ Treatment * Group * interval +",
                            RE_IND_INTERVAL[["forced_design"]]))
      if (.v %in% IND_FORCE_BETA) glmmTMB(f, data = d, family = beta_family(link = "logit"))
      else lmerTest::lmer(f, data = d, REML = TRUE)
    }, error = function(e) NULL)
    if (is.null(.m_full) || is.null(.m_rel)) next
    u <- .beta_p_for(.m_full, "Treatment"); r <- .beta_p_for(.m_rel, "Treatment")
    .id_rows[[length(.id_rows) + 1]] <- data.frame(
      label = .v, term = "Treatment",
      n_rows_full = nrow(ind_interval), n_rows_reliable = nrow(.ind_rel),
      p_full = u[2], p_reliable = r[2],
      p_crosses_alpha = isTRUE((u[2] < 0.05) != (r[2] < 0.05)),
      stringsAsFactors = FALSE)
  }
}
.id_all <- if (length(.id_rows)) bind_rows(.id_rows) else
  data.frame(label = NA_character_, term = NA_character_, n_rows_full = NA_integer_,
             n_rows_reliable = NA_integer_, p_full = NA_real_, p_reliable = NA_real_,
             p_crosses_alpha = NA)
write.csv(.id_all, file.path(DIR_ANOVAS, "sensitivity_identity_tier.csv"), row.names = FALSE)
cat("  wrote sensitivity_identity_tier.csv (", nrow(.id_all), " rows)\n", sep = "")

################################################################################
# pb_lrt_<var>_<term>.csv -- HANDOFF S7.
#
# Kenward-Roger's protection DOES NOT EXIST for a Wald chi-square: it carries
# no denominator df, so a Wald test on a 96-row fish-level beta/NB GLMM gets no
# downweighting whatever for clustering and is ANTICONSERVATIVE on the
# between-trial terms. The parametric-bootstrap LRT is the honest counterpart,
# and it must be fitted on the tier that is REPORTED -- previously it ran on
# ind_trial while the manuscript quoted interval-level Wald statistics.
#
# Coverage extends to every Wald-chi-square outcome carrying a between-trial
# claim, not just the two individual beta indicators: zone_synchrony and
# mean_polarization (GD, beta) and n_zone_transitions (NB) get no df
# protection either.
#
# N is a REPORTING CONSTRAINT (A-14): the smallest reportable p is 1/(N+1), so
# a reduced-N pass locates borderline cases and only those are re-run at full N
# (S7 [DR]). Convergence below 50% is flagged as an abort, not reported as a p.
################################################################################
PB_LRT_NSIM_REDUCED <- 200L
PB_LRT_NSIM_FULL    <- 1000L
PB_LRT_SEED         <- 20260824L

.run_pb_lrt_for <- function(v, data, re_str, run_dir, fam = "beta") {
  d <- as.data.frame(data)
  d <- d[is.finite(d[[v]]), , drop = FALSE]
  # pb_lrt() (house_style.R A12) hardcodes the response column name
  # as ".y" in its internal formula construction.
  d$.y <- if (fam == "beta") adjust_beta_01(d[[v]]) else round(d[[v]])
  # +/-1 sum-coded numeric copies of the two-level factors, for the 1-df marginal
  # tests below (house_style.R A12b).
  d <- add_sum_coded_cols(d)
  family_obj <- if (fam == "beta") beta_family(link = "logit") else glmmTMB::nbinom2()
  # WHAT THESE TESTS ACTUALLY COMPARE.
  #
  # A likelihood-ratio test compares two nested models, so dropping a main
  # effect from a fully-crossed fixed structure necessarily drops its
  # interactions too -- a reduced model retaining Group:interval but not Group
  # violates marginality and is not a model worth fitting. Each row below is
  # therefore a test of EVERY TERM INVOLVING the named factor, not of that
  # factor's marginal effect, and `terms_dropped` / `df_lrt` record which.
  #
  # This is a real difference from the Kenward-Roger F-tests in anova.csv,
  # which ARE marginal (Type III under sum-to-zero coding). Where the two
  # disagree, they are answering different questions and both are reported.
  #
  # The marginal question is answered by the two extra rows added below
  # ("Group_1df", "Treatment_1df"): the factors are coded +/-1 and ONLY the
  # main-effect column is removed from the reduced model, so the interactions
  # stay in both models and the test has exactly 1 df. See pb_1df_specs().
  terms <- list(
    Treatment = list(full = "Treatment * Group * interval",
                     reduced = "Group * interval",
                     dropped = "Treatment, Treatment:Group, Treatment:interval, Treatment:Group:interval",
                     df = 4L),
    Group     = list(full = "Treatment * Group * interval",
                     reduced = "Treatment * interval",
                     dropped = "Group, Treatment:Group, Group:interval, Treatment:Group:interval",
                     df = 4L))
  terms[["Treatment:Group"]] <- list(
    full = "Treatment * Group * interval",
    reduced = "Treatment + Group + interval + Treatment:interval + Group:interval",
    dropped = "Treatment:Group, Treatment:Group:interval",
    df = 2L)
  for (.tt in c("Group", "Treatment")) terms[[paste0(.tt, "_1df")]] <- pb_1df_specs(.tt)
  for (term_name in names(terms)) {
    spec <- terms[[term_name]]
    res <- pb_lrt(spec$full, spec$reduced, re_str, d, family = family_obj,
                  nsim = PB_LRT_NSIM_REDUCED, seed = PB_LRT_SEED)
    if (is.finite(res$p_pb) && res$p_pb < 0.10) {
      cat("  [", v, "/", term_name, "] p_pb=", round(res$p_pb, 4),
          " near threshold at N=", PB_LRT_NSIM_REDUCED,
          " -- re-running at N=", PB_LRT_NSIM_FULL, "\n", sep = "")
      res <- pb_lrt(spec$full, spec$reduced, re_str, d, family = family_obj,
                    nsim = PB_LRT_NSIM_FULL, seed = PB_LRT_SEED)
    }
    conv_frac <- if (is.finite(res$N) && res$N > 0) res$nsim_ok / res$N else NA_real_
    out <- data.frame(variable = v, term = term_name,
                      terms_dropped = spec$dropped, df_lrt = spec$df,
                      is_marginal_test = spec$df == 1L,
                      lrt_obs = res$lrt_obs,
                      p_pb = res$p_pb, nsim_ok = res$nsim_ok, N = res$N,
                      conv_frac = conv_frac,
                      aborted_low_convergence = isTRUE(conv_frac < 0.5),
                      p_floor = 1 / (res$N + 1),
                      seed = res$seed, method = res$method)
    write.csv(out, file.path(run_dir, sprintf("pb_lrt_%s_%s.csv", v, gsub("[: ]", "", term_name))),
              row.names = FALSE)
    cat("  wrote pb_lrt_", v, "_", gsub("[: ]", "", term_name), ".csv (p_pb=",
        round(res$p_pb, 4), ", N=", res$N, ", conv=", round(conv_frac, 3), ")\n", sep = "")
  }
}

# Individual-tier Wald-chi-square outcomes, fitted on the REPORTED tier.
for (.v in IND_FORCE_BETA)
  tryCatch(.run_pb_lrt_for(.v, ind_interval, RE_IND_INTERVAL[["forced_design"]],
                           DIR_ANOVAS, fam = "beta"),
           error = function(e) message("PB-LRT failed for ", .v, ": ", conditionMessage(e)))
for (.v in IND_COUNT_VARS)
  tryCatch(.run_pb_lrt_for(.v, ind_interval, RE_IND_INTERVAL[["forced_design"]],
                           DIR_ANOVAS, fam = "nb"),
           error = function(e) message("PB-LRT failed for ", .v, ": ", conditionMessage(e)))
# Group-dynamics beta outcomes (trial x interval -- no fish level exists).
for (.v in GD_FORCE_BETA)
  tryCatch(.run_pb_lrt_for(.v, gd_interval, RE_GD_INTERVAL[["forced_design"]],
                           DIR_ANOVAS, fam = "beta"),
           error = function(e) message("PB-LRT failed for ", .v, ": ", conditionMessage(e)))

# -- multiverse_<endpoint>.csv: response-scale x aggregation x RE-structure
#    sweep for the primary endpoint (HANDOFF S7 -- ships as a supplementary
#    table). Spans all three aggregation tiers, since aggregation level is one
#    of the choices the data did not force. -----------------------------------
.mv_grid <- list(
  list(tier = "fish_interval",  which = "ind_interval",
       re = RE_IND_INTERVAL[["forced_design"]], fixed = "Treatment * Group * interval"),
  list(tier = "fish_interval",  which = "ind_interval",
       re = "(1 | trial_id)", fixed = "Treatment * Group * interval"),
  list(tier = "trial_interval", which = "ind_trial_interval",
       re = "(1 | trial_date) + (1 | trial_id)", fixed = "Treatment * Group * interval"),
  list(tier = "trial",          which = "ind_trial",
       re = "(1 | trial_date)", fixed = "Treatment * Group"))
.mv_rows <- lapply(seq_along(.mv_grid), function(i) {
  g <- .mv_grid[[i]]
  dd <- switch(g$which, ind_interval = ind_interval,
               ind_trial_interval = ind_trial_interval, ind_trial = ind_trial)
  d <- as.data.frame(dd)
  d <- d[is.finite(d[[PRIMARY_ENDPOINT]]), , drop = FALSE]
  d$.y <- adjust_beta_01(d[[PRIMARY_ENDPOINT]])
  fmla <- as.formula(paste(".y ~", g$fixed, "+", g$re))
  m <- tryCatch(glmmTMB(fmla, data = d, family = beta_family(link = "logit")),
                error = function(e) NULL)
  if (is.null(m)) return(data.frame(endpoint = PRIMARY_ENDPOINT, tier = g$tier,
                                    re_spec = g$re, response_scale = "beta(0,1)",
                                    n_rows = nrow(d), p_treatment = NA_real_))
  an <- tryCatch(as.data.frame(car::Anova(m, type = 3)), error = function(e) NULL)
  p_t <- .pval_from(an, "Treatment")
  data.frame(endpoint = PRIMARY_ENDPOINT, tier = g$tier, re_spec = g$re,
             response_scale = "beta(0,1)", n_rows = nrow(d), p_treatment = p_t)
})
write.csv(bind_rows(.mv_rows), file.path(DIR_ANOVAS, sprintf("multiverse_%s.csv", PRIMARY_ENDPOINT)),
          row.names = FALSE)
cat("  wrote multiverse_", PRIMARY_ENDPOINT, ".csv\n", sep = "")

# -- equivalence_tost.csv: TOST for the primary endpoint's Treatment marginal
#    contrast, on the REPORTED tier. SESOI is BORROWED from the CHOICE bundle
#    (d = 0.70, ICC = 0.20; handoff/HANDOFF.md S6 EX) because no GOT-specific
#    power analysis exists in this repo -- labelled as borrowed so a GOT null
#    tested against it is never mis-reported as a GOT-native pre-registered
#    SESOI. -------------------------------------------------------------------
if (!is.null(.primary_model)) {
  .marg_contrast <- tryCatch(
    as.data.frame(pairs(emmeans(.primary_model, ~ Treatment, type = "link"))),
    error = function(e) NULL)
  if (!is.null(.marg_contrast) && nrow(.marg_contrast) >= 1) {
    sensitivity_tost(
      label = paste0(PRIMARY_ENDPOINT, "_Treatment_marginal_fish_interval"),
      est = .marg_contrast$estimate[1], se = .marg_contrast$SE[1],
      df = .marg_contrast$df[1], sesoi_d = 0.70,
      sesoi_source = "borrowed_from_CHOICE_bundle_not_GOT_native",
      run_dir = DIR_ANOVAS)
    cat("  wrote equivalence_tost.csv (SESOI borrowed from CHOICE, d=0.70)\n")
  }
}

write_csv_utf8bom(stmt_int_TxG,
  file.path(DIR_ANOVAS, "anova_stmt_interval_TreatmentxGroup.csv"))
write_csv_utf8bom(stmt_tri_TxG,
  file.path(DIR_ANOVAS, "anova_stmt_trial_TreatmentxGroup.csv"))
write_csv_utf8bom(stmt_int_TxGxI,
  file.path(DIR_ANOVAS, "anova_stmt_interval_TreatmentxGroupxInterval.csv"))
write_csv_utf8bom(
  bind_cols(stmt_int_T %>% rename(stmt_Treatment = statement),
            stmt_int_G %>% select(stmt_Group = statement)),
  file.path(DIR_ANOVAS, "anova_stmt_interval_main_effects.csv")
)
write_csv_utf8bom(
  bind_cols(stmt_tri_T %>% rename(stmt_Treatment = statement),
            stmt_tri_G %>% select(stmt_Group = statement)),
  file.path(DIR_ANOVAS, "anova_stmt_trial_main_effects.csv")
)

cat("\nANOVA statement CSVs written.\n")

################################################################################
###  POST-HOC (Tukey) + COMPACT-LETTER-DISPLAY  (→ post-hoc/ and cld/ subfolders)
################################################################################
# For every indicator, the full-interaction emmeans (interval level =
# Treatment×Group×interval; trial level = Treatment×Group) are subjected to
# Tukey pairwise post-hoc; compact-letter-display (CLD) summarises which cells
# differ. One CSV per indicator per level is written:
#   post-hoc/<level>_<key>.csv : full Tukey pairwise table
#   cld/<level>_<key>.csv      : one row per emmeans cell with its CLD letter
DIR_POSTHOC <- file.path(OUTPUT_DIR, "post-hoc")
DIR_CLD     <- file.path(OUTPUT_DIR, "cld")
for (.d in c(DIR_POSTHOC, DIR_CLD)) dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# Reconstruct a symmetric p-value matrix from a pairs()-on-flat-grid table.
# pairs()/cld() enumerate the SAME reference grid in the SAME row order, and
# pairs() (no `by=`) emits choose(n,2) contrasts in combn(n,2) order -- so
# this mapping is exact, not inferred from parsing contrast-name text.
.build_pmat_from_pairs <- function(n, p_values) {
  if (length(p_values) != choose(n, 2)) return(NULL)
  m <- matrix(1, n, n)
  idx <- utils::combn(n, 2)
  for (k in seq_len(ncol(idx))) {
    i <- idx[1, k]; j <- idx[2, k]
    m[i, j] <- p_values[k]; m[j, i] <- p_values[k]
  }
  m
}

export_posthoc_cld <- function(res_list, level) {
  emms <- extract_list(res_list, "emm")
  phos <- extract_list(res_list, "posthoc")
  np <- 0L; nc <- 0L
  for (v in names(res_list)) {
    ind_pretty <- pretty_variable_name(sub("^(GD_|Ind_)", "", v))
    unit_lab   <- if (grepl("^GD_", v)) "Group dynamics" else "Individual"

    if (!is.null(phos[[v]])) {
      ph <- tryCatch(as.data.frame(phos[[v]]), error = function(e) NULL)
      if (!is.null(ph) && nrow(ph) > 0) {
        ph <- cbind(key = v, level = level, unit = unit_lab, indicator = ind_pretty, ph)
        write_csv_utf8bom(ph, file.path(DIR_POSTHOC, paste0(level, "_", v, ".csv")))
        np <- np + 1L
      }
    }
    if (!is.null(emms[[v]])) {
      cl <- tryCatch(as.data.frame(multcomp::cld(emms[[v]], adjust = "tukey",
                                                 Letters = letters, alpha = 0.05)),
                     error = function(e) NULL)
      if (!is.null(cl) && ".group" %in% names(cl)) {
        cl$.group <- trimws(cl$.group)

        # CONVENTIONS.md C5 four-step pipeline: steps 2-4 (drop redundant
        # letters -> remap control-first -> guard), gated on being able to
        # build a Treatment:Group[:interval] cell key and a matching pmat.
        # MUST run before the .group -> cld rename below, since cld_pipeline()
        # reads the `.group` column by that name.
        if (all(c("Treatment", "Group") %in% names(cl))) {
          cl$cell <- paste0(as.character(cl$Treatment), ":", as.character(cl$Group))
          cell_order <- TG_BREAKS_PLOT
          if ("interval" %in% names(cl)) {
            ivs <- sort(unique(as.character(cl$interval)))
            cl$cell <- paste0(cl$cell, ":", as.character(cl$interval))
            cell_order <- as.vector(sapply(TG_BREAKS_PLOT, function(cc) paste0(cc, ":", ivs)))
          }
          ph_v <- tryCatch(as.data.frame(phos[[v]]), error = function(e) NULL)
          pmat <- if (!is.null(ph_v)) .build_pmat_from_pairs(nrow(cl), ph_v$p.value) else NULL
          if (!is.null(pmat)) {
            res <- cld_pipeline(cl, cell_ids = seq_len(nrow(cl)), pmat = pmat,
                                label = paste0(level, "_", v), cell_col = "cell",
                                cell_order = cell_order)
            if (!res$omitted) cl <- res$cld_df
          }
        }
        names(cl)[names(cl) == ".group"] <- "cld"
        cl <- cbind(key = v, level = level, unit = unit_lab, indicator = ind_pretty, cl)
        write_csv_utf8bom(cl, file.path(DIR_CLD, paste0(level, "_", v, ".csv")))
        nc <- nc + 1L
      }
    }
  }
  c(posthoc = np, cld = nc)
}
ph_int <- export_posthoc_cld(all_int, "interval")
ph_tri <- export_posthoc_cld(all_tri, "trial")
cat(sprintf("post-hoc/ : %d Tukey tables | cld/ : %d CLD tables\n",
            ph_int[["posthoc"]] + ph_tri[["posthoc"]],
            ph_int[["cld"]] + ph_tri[["cld"]]))


################################################################################
###  MARGINAL EMMEANS FOR SIGNIFICANT ANOVA TERMS  (→ marginal-emm/ subfolder)
################################################################################
# For every interval-level indicator with ANOVA terms significant at p < 0.05,
# compute marginal emmeans (averaging over non-focal factors), Tukey pairwise
# post-hoc, and CLD. One CSV per indicator × term saved to:
#   marginal-emm/posthoc/<level>_<key>_<term>.csv
#   marginal-emm/cld/<level>_<key>_<term>.csv
# Column `term` carries the ANOVA term string (e.g. "Group:interval").

DIR_MARGINAL     <- file.path(OUTPUT_DIR, "marginal-emm")
DIR_MARGINAL_PH  <- file.path(DIR_MARGINAL, "posthoc")
DIR_MARGINAL_CLD <- file.path(DIR_MARGINAL, "cld")
for (.d in c(DIR_MARGINAL_PH, DIR_MARGINAL_CLD))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# Convert "Group:interval" → formula ~ Group * interval
term_to_emm_spec <- function(term_str) {
  vars <- trimws(strsplit(term_str, ":")[[1]])
  as.formula(paste("~", paste(vars, collapse = " * ")))
}

# Extract p-value from a single ANOVA row data frame
anova_row_p <- function(row) {
  if ("Pr(>F)"    %in% names(row)) return(as.numeric(row[["Pr(>F)"]][1]))
  if ("Pr(>Chisq)"%in% names(row)) return(as.numeric(row[["Pr(>Chisq)"]][1]))
  NA_real_
}

export_marginal_emm <- function(res_list, level, p_thresh = 0.05) {
  n_ph <- 0L; n_cld <- 0L
  for (v in names(res_list)) {
    r <- res_list[[v]]
    if (is.null(r) || is.null(r$model) || is.null(r$anova)) next

    anv        <- as.data.frame(r$anova)
    anv$term   <- rownames(anv)
    skip_pat   <- "^\\(Intercept\\)$|^Intercept$"
    sig_terms  <- anv$term[!grepl(skip_pat, anv$term)]
    sig_terms  <- sig_terms[vapply(sig_terms, function(t) {
      p <- anova_row_p(anv[anv$term == t, , drop = FALSE])
      !is.na(p) && p < p_thresh
    }, logical(1))]

    if (length(sig_terms) == 0) next

    ind_pretty <- pretty_variable_name(sub("^(GD_|Ind_)", "", v))
    unit_lab   <- if (grepl("^GD_", v)) "Group dynamics" else "Individual"

    for (tm in sig_terms) {
      safe_nm <- gsub(":", "x", tm)          # safe for file names
      spec    <- tryCatch(term_to_emm_spec(tm), error = function(e) NULL)
      if (is.null(spec)) next

      emm_obj <- tryCatch(
        emmeans(r$model, spec, type = "response"),
        error = function(e) NULL
      )
      if (is.null(emm_obj)) next

      # Post-hoc pairwise
      ph <- tryCatch(as.data.frame(pairs(emm_obj, adjust = "tukey")),
                     error = function(e) NULL)
      if (!is.null(ph) && nrow(ph) > 0) {
        ph <- cbind(key = v, level = level, unit = unit_lab,
                    indicator = ind_pretty, term = tm, ph)
        write_csv_utf8bom(ph,
          file.path(DIR_MARGINAL_PH,
                    paste0(level, "_", v, "_", safe_nm, ".csv")))
        n_ph <- n_ph + 1L
      }

      # CLD -- CONVENTIONS.md C5 steps 2-4 applied wherever the marginal grid
      # still carries both Treatment and Group (e.g. the Treatment:Group and
      # Treatment:Group:interval terms); a grid marginalised down to a single
      # factor (e.g. Group:interval, no Treatment) has no "control-first"
      # anchor to remap onto, so only step 1 (cld() itself, unchanged) runs
      # for those -- stated here rather than silently forcing a remap that
      # would not mean what C5's worked example means.
      cl <- tryCatch(
        as.data.frame(multcomp::cld(emm_obj, adjust = "tukey",
                                    Letters = letters, alpha = 0.05)),
        error = function(e) NULL
      )
      if (!is.null(cl) && ".group" %in% names(cl)) {
        cl$.group <- trimws(cl$.group)

        # See export_posthoc_cld() above: this MUST run before the .group ->
        # cld rename, since cld_pipeline() reads the `.group` column by name.
        if (all(c("Treatment", "Group") %in% names(cl))) {
          cl$cell <- paste0(as.character(cl$Treatment), ":", as.character(cl$Group))
          cell_order <- TG_BREAKS_PLOT
          if ("interval" %in% names(cl)) {
            ivs <- sort(unique(as.character(cl$interval)))
            cl$cell <- paste0(cl$cell, ":", as.character(cl$interval))
            cell_order <- as.vector(sapply(TG_BREAKS_PLOT, function(cc) paste0(cc, ":", ivs)))
          }
          pmat <- if (!is.null(ph)) .build_pmat_from_pairs(nrow(cl), ph$p.value) else NULL
          if (!is.null(pmat)) {
            res <- cld_pipeline(cl, cell_ids = seq_len(nrow(cl)), pmat = pmat,
                                label = paste0(level, "_", v, "_", safe_nm),
                                cell_col = "cell", cell_order = cell_order)
            if (!res$omitted) cl <- res$cld_df
          }
        }
        names(cl)[names(cl) == ".group"] <- "cld"

        cl <- cbind(key = v, level = level, unit = unit_lab,
                    indicator = ind_pretty, term = tm, cl)
        write_csv_utf8bom(cl,
          file.path(DIR_MARGINAL_CLD,
                    paste0(level, "_", v, "_", safe_nm, ".csv")))
        n_cld <- n_cld + 1L
      }
    }
  }
  c(posthoc = n_ph, cld = n_cld)
}

mg_int <- export_marginal_emm(all_int, "interval")
mg_tri <- export_marginal_emm(all_tri, "trial")
cat(sprintf("marginal-emm/ : %d post-hoc | %d CLD tables\n",
            mg_int[["posthoc"]] + mg_tri[["posthoc"]],
            mg_int[["cld"]] + mg_tri[["cld"]]))


################################################################################
###  ADDITIONAL ANOVA STATEMENT CSVs (Group×interval + interval main effect)
################################################################################
# These two terms drive the row-2 emmeans panels in Figure 5 (interval-individual).
stmt_int_GxI <- make_statement_table(all_int_anova, all_int_model,
                                      "Group:interval", "Group × interval")
stmt_int_I   <- make_statement_table(all_int_anova, all_int_model,
                                      "interval", "Interval")
write_csv_utf8bom(stmt_int_GxI,
  file.path(DIR_ANOVAS, "anova_stmt_interval_GroupxInterval.csv"))
write_csv_utf8bom(stmt_int_I,
  file.path(DIR_ANOVAS, "anova_stmt_interval_Interval.csv"))
cat("Extra ANOVA statement CSVs (Group×interval, Interval) written.\n")


################################################################################
###  EXPORT WORD DECISION TABLES  (→ reports/ subfolder)
################################################################################

docx_path <- file.path(DIR_REPORTS, "GOT_easy_script_model_decisions.docx")

doc <- read_docx()
doc <- body_add_par(doc, "GOT Model Decision Tables", style = "heading 1")
doc <- body_add_par(doc,
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")), style = "Normal")
doc <- body_add_par(doc,
  paste("Fixed: Treatment × Group × interval (interval-level);",
        "Treatment × Group (trial-aggregated). RE selected by AICc.",
        "Data: ExeChoice + NoChoice groups only (pharma excluded). TP2 + TP3."),
  style = "Normal")
doc <- body_add_par(doc,
  paste("Column key — SW p (LMM resid.): Shapiro-Wilk p on final LMM residuals;",
        "p >= 0.05 = normality satisfied (fail to reject H0), p < 0.05 = violated.",
        "Norm. violated?: TRUE = OLS-residual screening flagged a problem",
        "(SW p < 0.05 or Levene p < 0.05 on best-transformed OLS residuals, used",
        "for transform selection before LMM fitting); FALSE = screening passed;",
        "NA = not applicable for Beta / NB GLMM families."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc, "Interval-level — Group dynamics",   style = "heading 2")
doc <- body_add_flextable(doc, make_decision_flextable(gd_int_model_dt))
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc, "Interval-level — Individual fish",  style = "heading 2")
doc <- body_add_flextable(doc, make_decision_flextable(ind_int_model_dt))
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc, "Trial-aggregated — Group dynamics", style = "heading 2")
doc <- body_add_flextable(doc, make_decision_flextable(gd_tri_model_dt))
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc, "Trial-aggregated — Individual fish",style = "heading 2")
doc <- body_add_flextable(doc, make_decision_flextable(ind_tri_model_dt))

print(doc, target = docx_path)
cat("\nWord document:", docx_path, "\n")


################################################################################
###  FINAL SUMMARY
################################################################################

cat("\n\n=== FINAL MODEL DECISION SUMMARIES ===\n")
cat("\n--- GD interval-level ---\n");    print(gd_int_model_dt)
cat("\n--- GD trial-aggregated ---\n");  print(gd_tri_model_dt)
cat("\n--- Individual interval-level ---\n"); print(ind_int_model_dt)
cat("\n--- Individual trial-aggregated ---\n"); print(ind_tri_model_dt)

cat("\n\n=== FILES WRITTEN TO", OUTPUT_DIR, "===\n")
cat("datasets/\n")
cat("  gd_interval.csv  gd_trial.csv  ind_interval.csv  ind_trial.csv\n")
cat("anovas/\n")
cat("  anova_stmt_interval_TreatmentxGroup.csv\n")
cat("  anova_stmt_interval_TreatmentxGroupxInterval.csv\n")
cat("  anova_stmt_interval_main_effects.csv\n")
cat("  anova_stmt_trial_TreatmentxGroup.csv\n")
cat("  anova_stmt_trial_main_effects.csv\n")
cat("reports/\n")
cat("  GOT_interval_Treatment_x_Group_x_interval.md\n")
cat("  GOT_interval_Treatment_x_Group_x_interval.docx\n")
cat("  GOT_trial_aggregated_Treatment_x_Group.md\n")
cat("  GOT_trial_aggregated_Treatment_x_Group.docx\n")
cat("  GOT_easy_script_model_decisions.docx\n")

# HANDOFF.md S10 [INV]: a completed-run sentinel, checked by downstream
# auto-discovery before a run directory is trusted (see figure_3_locomotion.R).
writeLines(as.character(Sys.time()), file.path(OUTPUT_DIR, "_RUN_COMPLETE"))


cat("\n=== behaviour_models.R complete ===\n")
