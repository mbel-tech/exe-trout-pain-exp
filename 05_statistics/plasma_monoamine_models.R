################################################################################
# plasma_monoamine_models.R
# ─────────────────────────────────────────────────────────────────────────────
# Variant of an earlier plasma/monoamine script (not included) with Rosner's generalised ESD
# outlier test (EnvStats::rosnerTest, k = 3 max outliers) applied to EACH indicator
# BEFORE model fitting. True outliers are logged (MD + Word table) and dropped from
# that indicator's analysis dataset. Indicators are processed concentrations-first,
# ratios last: before testing a ratio, the true outliers of its component
# concentrations are dropped first (e.g. 5-HT and 5-HIAA outliers are removed before
# Rosner's test on 5-HIAA/5-HT), then the ratio is tested for its own outliers.
# Rosner's test is run on the raw indicator scale.
#
# Primary analysis of GOT plasma metabolites/cortisol and brain monoamines,
# adapted from behaviour_models.R (methodology + reporting/output style).
#
# Design (confirmed):
#   - Fixed effects: indicator ~ Treatment * Group   (no interval/repeated measures)
#   - Group scope: Exercise Choice + No Choice only (pharma & plasma baselines dropped)
#   - Monoamines analysed per brain region (DM, RN) separately (no cross-area model)
#   - Monoamine concentrations analysed as provided (no protein normalization);
#     the two turnover ratios (5-HIAA/5-HT, DOPAC/DA) analysed directly
#
# Modelling (per indicator):
#   1. Fit the design MIXED model (lmer); AICc ranks the candidate RE structures
#      (ML). The chosen structure is refit by REML and, if SINGULAR, REDUCED: the
#      offending variance component is dropped and the model refitted until it is
#      non-singular (a singular fit is never retained -- LOGS.md 2026-08-26
#      (later)). If every random effect is dropped the model is an ordinary lm,
#      labelled 'Gaussian LM (no RE retained)'. Every drop is written to
#      anovas/re_reduction.csv.
#   2. Shapiro-Wilk on selected-model residuals.
#   3. If non-normal (p < 0.05): optimise a Box-Cox lambda that MAXIMISES the
#      residual Shapiro p (grid search; lambda ~ 0 -> natural log); refit.
#   4. Adopt the transformed model only if residual Shapiro p improves.
#   Family: Gaussian LMM/LM after Box-Cox (all indicators positive continuous).
#           Indicators that remain non-normal are flagged (Gamma GLMM = documented
#           future fallback).
#
# Inference: Type III ANOVA (lmer -> lmerTest Satterthwaite; lm -> car::Anova);
#            emmeans ~ Treatment*Group; Tukey pairwise (focal = Treatment×Group).
#
# Outputs (output/<timestamp>/):
#   datasets/  cleaned analysis-ready CSVs (plasma, mono_DM, mono_RN)
#   anovas/    ANOVA-statement CSVs (UTF-8 BOM): Treatment×Group, Treatment, Group
#   reports/   per-module markdown + Word reports; combined model-decision table
#              (markdown + Word); raw_means_SE.md (+ CSV, figure-ready)
################################################################################

suppressMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble)
  library(lme4); library(lmerTest); library(glmmTMB); library(MuMIn); library(emmeans)
  library(car); library(officer); library(flextable)
  library(pbkrtest); library(performance)
})
# NOTE: EnvStats is NOT attached with library() because it masks base::print,
# which breaks officer's print.rdocx (Word output). It is used namespace-qualified
# as EnvStats::rosnerTest(). Ensure it is installed.
if (!requireNamespace("EnvStats", quietly = TRUE))
  stop("Package 'EnvStats' is required (install.packages('EnvStats')).")

options(contrasts = c("contr.sum", "contr.poly"), stringsAsFactors = FALSE)
select <- dplyr::select   # guard against MASS/other masking

set.seed(20260608)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Paths can be overridden from the environment (GOT_ROOT, GOT_DATA_DIR,
# GOT_OUT_BASE) so the script runs from a worktree or a re-lettered drive.
.env_or <- function(var, default) { v <- Sys.getenv(var, ""); if (nzchar(v)) v else default }
GOT_ROOT  <- .env_or("GOT_ROOT", ".")
DATA_DIR  <- .env_or("GOT_DATA_DIR", file.path(GOT_ROOT, "data"))
PLASMA_XLSX <- file.path(DATA_DIR, "GOT glu_lac_cort updated march 2025.xlsx")
# Monoamines are read from the repaired tidy dataset, NOT from
# "monoamine concentration ALL GOT.xlsx": that workbook's DM 5-HT column is a
# fill-drag artefact in 88 of 92 rows (see README_GOT_monoamine_data_AUG2026.txt).
# Reading it made every true DM 5-HT value look like an outlier.
MONO_CSV    <- .env_or("GOT_MONO_CSV", file.path(DATA_DIR, "GOT_monoamine_data_AUG2026.csv"))
# Further regions (Dl, Vv, POA; HPLC July 2026) are analysed with the default
# strategy when the dataset carries analyte values for them. Each indicator needs
# at least MIN_CELL_EXTRA fish in every Treatment x Group cell, before and after
# removal of outliers; the ones that do not are listed in
# anovas/extra_regions_not_analysed.csv and are not modelled.
EXTRA_REGIONS  <- c(DL = "Dl", VV = "Vv", POA = "POA")
MIN_CELL_EXTRA <- 3L

TIMESTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(.env_or("GOT_OUT_BASE",
                                file.path(GOT_ROOT, "output")),
                        paste0("rosner_", TIMESTAMP))
DIR_DATASETS <- file.path(OUTPUT_DIR, "datasets")
DIR_ANOVAS   <- file.path(OUTPUT_DIR, "anovas")
DIR_REPORTS  <- file.path(OUTPUT_DIR, "reports")
DIR_DIAG     <- file.path(OUTPUT_DIR, "diagnostics")
for (.d in c(OUTPUT_DIR, DIR_DATASETS, DIR_ANOVAS, DIR_REPORTS, DIR_DIAG))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# -- House style / sensitivity modules (single source of truth) -------------
source(file.path(GOT_ROOT, "05_statistics/shared/house_style.R"))
source(file.path(GOT_ROOT, "05_statistics/shared/sensitivity_checks.R"))

# Accumulators for aicc_selection.csv / anova_satterthwaite.csv / anova.csv
# (CONVENTIONS.md C12), flushed once all indicators have run.
.AICC_LOG <- list()
.aicc_log_add <- function(label, aicc_tbl) {
  aicc_tbl$label <- label
  .AICC_LOG[[length(.AICC_LOG) + 1]] <<- aicc_tbl
}
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

# HANDOFF S5/S7: Kenward-Roger primary, Satterthwaite the named sensitivity.
finalize_gaussian_inference <- function(model_final, label) {
  if (!inherits(model_final, "lmerModLmerTest") && !inherits(model_final, "merMod"))
    return(tryCatch(car::Anova(model_final, type = 3), error = function(e) NULL))  # lm fallback
  anova_kr  <- tryCatch(anova(model_final, type = 3, ddf = "Kenward-Roger"), error = function(e) NULL)
  anova_sat <- tryCatch(anova(model_final, type = 3, ddf = "Satterthwaite"), error = function(e) NULL)
  if (!is.null(anova_sat)) {
    satt_tbl <- as.data.frame(anova_sat); satt_tbl$term <- rownames(satt_tbl)
    if (!is.null(anova_kr)) {
      kr_tbl <- as.data.frame(anova_kr); kr_tbl$term <- rownames(kr_tbl)
      m <- merge(satt_tbl, kr_tbl[, c("term", "Pr(>F)")], by = "term",
                suffixes = c("_satterthwaite", "_kr"))
      m$p_crosses_alpha <- (m$`Pr(>F)_satterthwaite` < 0.05) != (m$`Pr(>F)_kr` < 0.05)
      if (any(m$p_crosses_alpha, na.rm = TRUE))
        cat("  *** KR/Satterthwaite disagree across p=0.05 for:",
            paste(m$term[which(m$p_crosses_alpha)], collapse = ", "), "***
")
      satt_tbl <- m
    }
    .satt_log_add(label, satt_tbl)
  }
  if (is.null(anova_kr)) anova_kr <- anova_sat
  anova_kr
}

anova_effect_size_table <- function(anova_obj, model_obj, model_type, label, r2m = NA_real_,
                                    r2c = NA_real_, icc = NA_real_) {
  if (is.null(anova_obj)) return(invisible(NULL))
  tab <- as.data.frame(anova_obj); tab$term <- rownames(tab)
  rows <- lapply(seq_len(nrow(tab)), function(i) {
    row <- tab[i, , drop = FALSE]
    if (all(c("NumDF", "DenDF", "F value", "Pr(>F)") %in% names(row))) {
      es <- effect_size("F", F = row$`F value`, df1 = row$NumDF, df2 = row$DenDF)
      data.frame(term = row$term, stat_type = "F", statistic = row$`F value`,
                df1 = row$NumDF, df2 = row$DenDF, p = row$`Pr(>F)`,
                eta2_p = es$es, omega2_p = es$es_secondary)
    } else if (all(c("Df", "F value", "Pr(>F)") %in% names(row))) {
      df2 <- if (!is.null(model_obj)) df.residual(model_obj) else NA_real_
      es <- effect_size("F", F = row$`F value`, df1 = row$Df, df2 = df2)
      data.frame(term = row$term, stat_type = "F", statistic = row$`F value`,
                df1 = row$Df, df2 = df2, p = row$`Pr(>F)`,
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

cat("\n=== plasma_monoamine_models.R (Rosner ESD outlier removal) ===\n")
cat("Output:", OUTPUT_DIR, "\n")
ROSNER_K <- 3   # maximum number of outliers tested per indicator

GROUP_KEEP <- c("Exercise Choice", "No Choice")

################################################################################
###  DATA INGESTION & CLEANING
################################################################################

cat("\n--- Loading & cleaning data ---\n")

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Read an Excel sheet and keep the FIRST occurrence of each named column,
# dropping blank/"None"/duplicate-named columns (sheets have junk trailing cols).
read_clean_xlsx <- function(path, sheet) {
  df <- suppressMessages(read_excel(path, sheet = sheet, .name_repair = "minimal"))
  nm <- trimws(names(df))
  keep <- nm != "" & !is.na(nm) & nm != "None" & !duplicated(nm)
  df <- df[, keep, drop = FALSE]
  names(df) <- nm[keep]
  df
}

# ── Plasma: one row per fish; `tank` col encodes Treatment; `Group ` (trailing
#    space) holds condition. Indicators: cortisol, glucose, lactate. Per-indicator
#    assay plate = batch RE candidate. ────────────────────────────────────────
plasma_raw <- read_clean_xlsx(PLASMA_XLSX, "results")

plasma <- plasma_raw %>%
  transmute(
    Sample          = as.character(Sample),
    Treatment       = str_to_title(trimws(as.character(tank))),   # tank col = treatment
    Group           = trimws(as.character(Group)),
    plasma_cortisol = num(plasma_cortisol),
    Glucose_mmol_L  = num(Glucose_mmol_L),
    Lactate_mmol_L  = num(Lactate_mmol_L),
    cort_plate      = as.character(plate),       # cort assay plate
    glu_plate       = as.character(glu_plate),
    lac_plate       = as.character(lac_plate),
    sampling_date   = as.character(sampling_date)
  ) %>%
  filter(Group %in% GROUP_KEEP) %>%
  mutate(
    Treatment     = factor(Treatment, levels = c("Control", "Exercise")),
    Group         = factor(Group,     levels = GROUP_KEEP),
    sampling_date = factor(sampling_date)
  )

cat("plasma: N =", nrow(plasma), "| cells:\n"); print(table(plasma$Treatment, plasma$Group))

# Per-indicator batch (assay plate) column map for the (1|plate) RE candidate
PLASMA_BATCH <- c(plasma_cortisol = "cort_plate",
                  Glucose_mmol_L  = "glu_plate",
                  Lactate_mmol_L  = "lac_plate")

# ── Monoamine: one row per fish per Brain_area; analyse DM/RN separately. Two
#    turnover ratios analysed directly. Normalize case; make R-safe names. ─────
mono_raw <- suppressMessages(readr::read_csv(MONO_CSV, show_col_types = FALSE, guess_max = 5000))

.canon_group <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(grepl("^exe", x) ~ "Exercise Choice",
                   grepl("^no", x)  ~ "No Choice",
                   TRUE ~ str_to_title(x))
}
mono <- mono_raw %>%
  transmute(
    Sample_ID      = as.character(sample_id),
    Brain_area     = trimws(as.character(brain_region)),
    user           = trimws(as.character(operator)),
    Treatment      = str_to_title(trimws(as.character(treatment))),
    Group          = .canon_group(group),
    Trial          = as.character(trial),
    HT5            = num(HT5),
    HIAA5          = num(HIAA5),
    ratio_5HIAA_5HT = num(HIAA5_HT5),
    DA             = num(DA),
    DOPAC          = num(DOPAC),
    ratio_DOPAC_DA = num(DOPAC_DA)
  ) %>%
  filter(Group %in% GROUP_KEEP, Brain_area %in% c("DM", "RN", unname(EXTRA_REGIONS))) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Exercise")),
    Group     = factor(Group,     levels = GROUP_KEEP),
    Trial     = factor(Trial),
    user      = factor(user)
  )

mono_DM <- mono %>% filter(Brain_area == "DM") %>% droplevels()
mono_RN <- mono %>% filter(Brain_area == "RN") %>% droplevels()

cat("\nmono_DM: N =", nrow(mono_DM), "| cells:\n"); print(table(mono_DM$Treatment, mono_DM$Group))
cat("\nmono_RN: N =", nrow(mono_RN), "| cells:\n"); print(table(mono_RN$Treatment, mono_RN$Group))

################################################################################
###  DESIGN KEYS -- recovering the physical trial and testing date
################################################################################
# HANDOFF.md S4.3 [INV]: "A behavioural trial and an endocrine sample are linked
# by (tank, testing date, treatment), NOT by the raw numeric `trial` field alone
# -- the same integer label is independently assigned within each data stream
# and does not denote the same physical trial across streams. Joining on the
# bare integer silently mis-pairs a fraction of the dataset."
#
# That is exactly what was happening here. The monoamine sheet's `Trial` column
# runs 1-3, and every one of those three levels POOLS FISH FROM ALL FOUR
# Treatment x Group cells (14-16 fish per level). A random intercept on bare
# `Trial` therefore groups together animals from different treatments and
# different testing days -- it is not the design random effect, it accounts for
# no trial clustering whatever, and it was the RE in force for every monoamine
# model in the previous vintage.
#
# The physical trial is the COMPOSITE key. Verified against the data:
#   mono   (Trial, Treatment, Group) -> 12 groups of ~4 fish  = the 12 trials
#   mono   (Trial, Group)            ->  6 groups of ~8 fish  = the 6 dates
#   plasma (sampling_date, Treatment)-> 12 groups of  4 fish  = the 12 trials
#   plasma  sampling_date            ->  6 levels             = the 6 dates
#
# And the two streams agree: mono's (Trial, Group) maps 1:1 onto the same six
# testing dates the plasma sheet and the behaviour side carry
# (2024-11-14/15/20/22/26/27), cross-checked through the shared sample
# identifier (plasma$Sample == mono$Sample_ID, 48/48 overlap).
#
# This also RESOLVES the standing "plasma has no design RE" limitation: plasma
# does carry the trial, encoded jointly across two columns rather than in one.
.add_design_keys_mono <- function(d) {
  d$trial_uid <- factor(paste(d$Trial, d$Treatment, d$Group, sep = "_"))
  d$date_uid  <- factor(paste(d$Trial, d$Group, sep = "_"))
  d
}
.add_design_keys_plasma <- function(d) {
  d$trial_uid <- factor(paste(d$sampling_date, d$Treatment, sep = "_"))
  d$date_uid  <- factor(d$sampling_date)
  d
}
mono_DM <- .add_design_keys_mono(mono_DM)
mono_RN <- .add_design_keys_mono(mono_RN)
plasma  <- .add_design_keys_plasma(plasma)

# Fail the run rather than fit a model on a mis-specified design key.
.assert_design_keys <- function(d, label, expect_trials = 12L, expect_dates = 6L) {
  nt <- nlevels(droplevels(d$trial_uid)); nd <- nlevels(droplevels(d$date_uid))
  if (nt != expect_trials || nd != expect_dates)
    stop(sprintf("%s: design keys resolve to %d trials / %d dates, expected %d / %d",
                 label, nt, nd, expect_trials, expect_dates), call. = FALSE)
  cat(sprintf("  %s design keys OK: %d trials, %d dates (mean %.1f fish per trial)\n",
              label, nt, nd, nrow(d) / nt))
  invisible(TRUE)
}
cat("\nDesign-key resolution:\n")
.assert_design_keys(mono_DM, "mono_DM")
.assert_design_keys(mono_RN, "mono_RN")
.assert_design_keys(plasma,  "plasma")

# Extra regions: one data frame per region that has any analyte value.
mono_X <- list()
for (.k in names(EXTRA_REGIONS)) {
  .dx <- mono %>% filter(Brain_area == EXTRA_REGIONS[[.k]]) %>% droplevels()
  if (!nrow(.dx) || !any(is.finite(unlist(.dx[c("HT5", "HIAA5", "DA", "DOPAC")])))) next
  mono_X[[.k]] <- .add_design_keys_mono(.dx)
  .assert_design_keys(mono_X[[.k]], paste0("mono_", .k))
}

write_csv(plasma,  file.path(DIR_DATASETS, "plasma.csv"))
write_csv(mono_DM, file.path(DIR_DATASETS, "mono_DM.csv"))
write_csv(mono_RN, file.path(DIR_DATASETS, "mono_RN.csv"))
for (.k in names(mono_X)) write_csv(mono_X[[.k]], file.path(DIR_DATASETS, paste0("mono_", .k, ".csv")))

################################################################################
###  RANDOM-EFFECT CANDIDATE SETS (mixed-only; lm = singularity fallback)
################################################################################
# Both endocrine streams share the behaviour side's design, so they take the
# same forced structure (see the behaviour script's RANDOM-EFFECT STRUCTURE
# block for the full argument):
#
#   (1|date_uid)   -- unit of replication for GROUP (Group varies only between
#                     testing dates; 3 dates per level)
#   (1|trial_uid)  -- level at which TREATMENT varies (within date)
#
# Both are DESIGN random effects under HANDOFF S2 [INV] and are forced, never
# AICc-selected. The unit here is one fish, so no third level is needed: unlike
# the behavioural data there is one measurement per animal, so no
# repeated-measures term exists.
#
# `user` (HPLC operator) is RETIRED as a candidate, not out-competed. The
# design_balance.csv screen found it has exactly ONE level within each of the
# DM and RN subsets -- it was never an independent random effect there, only
# nominally a candidate, and offering AICc a degenerate alternative to choose
# against misrepresents what was checked.
#
# `batch` (assay plate) IS a genuine nuisance grouping and remains the one
# thing AICc is allowed to decide, per S5's "AICc is only ever allowed to
# decide nuisance terms."
RE_PLASMA <- c(
  forced_design       = "(1 | date_uid) + (1 | trial_uid)",
  forced_design_plate = "(1 | date_uid) + (1 | trial_uid) + (1 | batch)"
)
RE_MONO <- c(
  forced_design = "(1 | date_uid) + (1 | trial_uid)"
)

# Declared comparators for sensitivity_forced_re.csv (HANDOFF S7). The primary
# structure is design-determined, so the robustness question is what the
# conclusion would have been under the structures the design argument rules
# out -- including, pointedly, the bare-`Trial` RE the previous vintage used.
RE_COMPARATORS_MONO <- c(
  date_dropped     = "(1 | trial_uid)",
  bare_trial_index = "(1 | Trial)"
)
RE_COMPARATORS_PLASMA <- c(
  date_dropped = "(1 | trial_uid)",
  no_design_re = "(1 | batch)"
)

################################################################################
###  HELPER FUNCTIONS
################################################################################

shapiro_residual_p <- function(model) {
  res <- residuals(model); res <- res[is.finite(res)]
  if (length(res) < 3) return(NA_real_)
  if (length(res) > 5000) res <- sample(res, 5000)
  tryCatch(shapiro.test(res)$p.value, error = function(e) NA_real_)
}

# Box-Cox lambda chosen to MAXIMISE the Shapiro-Wilk p of OLS residuals under the
# fixed-effects model. Strictly positive y only (classic Box-Cox). lambda~0 -> log.
# Returns list(y_bc, lambda, name, sw) or NULL.
estimate_boxcox_shapiro <- function(y_valid, fixed_str, data_fit) {
  if (any(!is.finite(y_valid)) || min(y_valid) <= 0) return(NULL)
  d <- data_fit
  bct <- function(y, lam) if (abs(lam) < 1e-8) log(y) else (y^lam - 1) / lam
  .sw <- function(yt) {
    d$.y <- yt
    r <- tryCatch(residuals(lm(as.formula(paste(".y ~", fixed_str)), data = d)),
                  error = function(e) NULL)
    if (is.null(r)) return(NA_real_)
    r <- r[is.finite(r)]
    if (length(r) < 3) return(NA_real_)
    tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_)
  }
  coarse <- seq(-2, 2, by = 0.1)
  swc <- vapply(coarse, function(l) .sw(bct(y_valid, l)), numeric(1))
  if (all(is.na(swc))) return(NULL)
  l0 <- coarse[which.max(replace(swc, is.na(swc), -Inf))]
  fine <- seq(l0 - 0.1, l0 + 0.1, by = 0.02)
  swf <- vapply(fine, function(l) .sw(bct(y_valid, l)), numeric(1))
  lam <- fine[which.max(replace(swf, is.na(swf), -Inf))]
  if (abs(lam) < 0.05) {
    lam <- 0; y_bc <- log(y_valid); nm <- "box-cox (log, lambda=0)"
  } else {
    y_bc <- (y_valid^lam - 1) / lam; nm <- sprintf("box-cox (lambda=%.2f)", lam)
  }
  list(y_bc = y_bc, lambda = lam, name = nm, sw = max(swf, na.rm = TRUE))
}

# ---- Shared pieces of the three fit_select* functions ------------------------
# Rank the converged candidates by AICc (ML) and pick the structure the
# reduction starts from: the lowest-AICc NON-SINGULAR one, else the lowest-AICc
# one, else (nothing converged) the first candidate's structure.
.pick_start <- function(cand, re_cands) {
  aicc_tbl <- NULL; start <- NULL
  if (length(cand) > 0) {
    aiccv <- vapply(cand, `[[`, numeric(1), "aicc")
    sing  <- vapply(cand, `[[`, logical(1), "singular")
    pool  <- if (any(!sing)) which(!sing) else seq_along(cand)
    # a non-finite AICc (degenerate fit) must not crash the pick: rank it last.
    rk <- aiccv[pool]; rk[!is.finite(rk)] <- Inf
    start <- cand[[pool[which.min(rk)]]]
    aicc_tbl <- tibble(
      RE = vapply(cand, `[[`, character(1), "name"),
      re_formula = vapply(cand, `[[`, character(1), "re"),
      AICc = round(aiccv, 2),
      delta_AICc = round(aiccv - min(aiccv), 2),
      singular = sing
    )
  }
  list(name = if (!is.null(start)) start$name else names(re_cands)[1],
       re   = if (!is.null(start)) start$re   else re_cands[[1]],
       aicc_tbl = aicc_tbl)
}

# Package the outcome of reduce_singular_re(): name the result, add the no-RE
# row to the AICc table when every random effect was dropped, mark the selected
# row, and return the reduction trail for re_reduction.csv.
.reduction_result <- function(red, pick, no_re_name, no_re_aicc) {
  none <- is.na(red$re) || identical(red$re, "none")
  aicc_tbl <- pick$aicc_tbl
  if (none) {
    row <- tibble(RE = no_re_name, re_formula = "(none)", AICc = round(no_re_aicc, 2),
                  delta_AICc = NA_real_, singular = FALSE)
    aicc_tbl <- if (is.null(aicc_tbl)) row else bind_rows(aicc_tbl, row)
  }
  if (!is.null(aicc_tbl))
    aicc_tbl$selected <- aicc_tbl$RE == (if (none) no_re_name else pick$name)
  list(re_chosen    = if (none) no_re_name
                      else if (length(red$dropped)) paste0(pick$name, "_reduced") else pick$name,
       re_formula   = if (none) "(none)" else red$re,
       re_attempted = pick$re,
       re_retained  = if (none) "none" else red$re,
       dropped = red$dropped, reasons = red$reasons,
       aicc_table = aicc_tbl)
}

# Mixed-first AICc model selection, then singular-RE REDUCTION.
#   1. Fit each mixed RE candidate by ML; rank by AICc (.pick_start).
#   2. Refit the chosen structure by REML for inference through
#      reduce_singular_re(): a singular fit is NOT retained. The component with
#      the smallest estimated SD is dropped and the model refitted, until the fit
#      is non-singular. If every random effect is dropped the result is an
#      ordinary lm.
# Why (LOGS.md 2026-08-26, "Singular random effects are now dropped, not
# retained"): a singular fit has a variance estimate on the boundary of the
# parameter space, which Kenward-Roger's derivation excludes, and a structure
# the data cannot identify costs power without buying Type I error protection.
# The design structure is still ATTEMPTED first, and every drop is returned in
# `dropped` / `reasons` so run_indicator() can write re_reduction.csv.
# `respcol` is the response column name already present in `d`.
fit_select <- function(d, respcol, fixed_str, re_cands, label = "") {
  d$.resp <- d[[respcol]]
  d <- d[is.finite(d$.resp), , drop = FALSE]

  # HANDOFF S5 step 3: candidates are compared by AICc fit with ML -- REML
  # likelihoods are not comparable across differing RE structures.
  cand <- lapply(names(re_cands), function(rn) {
    fml <- as.formula(paste(".resp ~", fixed_str, "+", re_cands[[rn]]))
    m <- tryCatch(suppressMessages(suppressWarnings(
            lmerTest::lmer(fml, data = d, REML = FALSE))), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf),
         singular = lme4::isSingular(m))
  })
  pick <- .pick_start(Filter(Negate(is.null), cand), re_cands)

  # HANDOFF S5 step 5: REML refit for inference; fit_fn(NULL) is the no-RE terminal model.
  fit_fn <- function(re_str) {
    if (is.null(re_str))
      return(stats::lm(as.formula(paste(".resp ~", fixed_str)), data = d))
    suppressMessages(suppressWarnings(lmerTest::lmer(
      as.formula(paste(".resp ~", fixed_str, "+", re_str)), data = d, REML = TRUE)))
  }
  red <- reduce_singular_re(fit_fn, split_re_terms(pick$re), label = label)
  if (is.null(red$model)) stop("fit_select: no model could be fitted for ", respcol)
  c(list(model = red$model, is_lm = inherits(red$model, "lm"), singular = FALSE,
         r2m = NA_real_, r2c = NA_real_, re_r2_delta = NA_real_),
    .reduction_result(red, pick, "lm (no RE retained)",
                      tryCatch(MuMIn::AICc(red$model), error = function(e) NA_real_)))
}
# Gamma GLMM (log link) RE selection and singular-RE reduction -- same rule as
# fit_select, using glmmTMB::glmmTMB + Gamma(link="log"); the no-RE terminal
# model is a plain Gamma glm. Applied ONLY when Gaussian+transform still fails
# Shapiro-Wilk. `respcol` must be a column that is strictly positive (min > 0).
# One singularity definition for every family: the reduction's own check.
.is_singular_gamma <- function(m) .is_singular_any(m)

fit_select_gamma <- function(d, respcol, fixed_str, re_cands, label = "") {
  d$.resp <- d[[respcol]]
  d <- d[is.finite(d$.resp) & d$.resp > 0, , drop = FALSE]

  cand <- lapply(names(re_cands), function(rn) {
    fml <- as.formula(paste(".resp ~", fixed_str, "+", re_cands[[rn]]))
    m <- tryCatch(suppressMessages(suppressWarnings(
            glmmTMB::glmmTMB(fml, data = d, family = Gamma(link = "log")))),
            error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf),
         singular = .is_singular_gamma(m))
  })
  pick <- .pick_start(Filter(Negate(is.null), cand), re_cands)

  fit_fn <- function(re_str) {
    if (is.null(re_str))
      return(stats::glm(as.formula(paste(".resp ~", fixed_str)), data = d, family = Gamma(link = "log")))
    suppressMessages(suppressWarnings(glmmTMB::glmmTMB(
      as.formula(paste(".resp ~", fixed_str, "+", re_str)), data = d, family = Gamma(link = "log"))))
  }
  red <- reduce_singular_re(fit_fn, split_re_terms(pick$re), label = label)
  if (is.null(red$model)) return(NULL)   # not even the no-RE glm converges
  c(list(model = red$model, is_lm = FALSE, is_glm = inherits(red$model, "glm"), singular = FALSE),
    .reduction_result(red, pick, "glm (no RE retained)",
                      tryCatch(MuMIn::AICc(red$model), error = function(e) NA_real_)))
}
# Gamma-GLMM assumption checks (run only when a Gamma GLMM/GLM is fitted).
# Returns a semicolon-delimited string listing any violated assumptions, or
# "all assumptions met". Assumptions checked:
#   (1) values must be strictly positive;
#   (2) zeros are not allowed (incompatible with the log link / Gamma support);
#   (3) residual variance scales ~ square of the mean (Gamma: Var ∝ μ²) — assessed
#       as the slope of log(cell variance) ~ log(cell mean) across Treatment×Group
#       cells, expected ≈ 2;
#   (4) sensitivity to very small near-zero values (smallest value ≪ median);
#   (5) random-effect structure must remain identifiable (not singular / not
#       collapsed to a no-RE GLM).
gamma_assumption_checks <- function(d, v, is_glm, singular) {
  y <- d[[v]]; y <- y[is.finite(y)]
  flags <- character(0)

  # (1) strictly positive
  if (length(y) == 0 || min(y) < 0) flags <- c(flags, "non-positive values present")
  # (2) exact zeros
  if (any(y == 0)) flags <- c(flags, "exact zeros present (invalid for log link)")
  # (3) mean-variance power ≈ 2 across Treatment×Group cells
  mv <- tryCatch(
    d %>% dplyr::filter(is.finite(.data[[v]])) %>%
      dplyr::group_by(Treatment, Group) %>%
      dplyr::summarise(m = mean(.data[[v]]), vv = var(.data[[v]]), .groups = "drop") %>%
      dplyr::filter(is.finite(m), is.finite(vv), m > 0, vv > 0),
    error = function(e) NULL)
  if (!is.null(mv) && nrow(mv) >= 3) {
    slope <- tryCatch(unname(coef(lm(log(vv) ~ log(m), data = mv))[2]),
                      error = function(e) NA_real_)
    if (!is.na(slope) && (slope < 1 || slope > 3))
      flags <- c(flags, sprintf("mean-variance power %.2f (Gamma=2; Tweedie 1-2)", slope))
  } else {
    flags <- c(flags, "mean-variance scaling not assessable (<3 usable cells)")
  }
  # (4) near-zero sensitivity
  med <- stats::median(y)
  if (length(y) > 0 && med > 0 && (min(y) / med) < 0.01)
    flags <- c(flags, sprintf("near-zero values (min/median = %.4f)", min(y) / med))
  # (5) random-effect identifiability
  if (isTRUE(is_glm))        flags <- c(flags, "all random effects dropped (singular) -> GLM")
  else if (isTRUE(singular)) flags <- c(flags, "RE reduced (singular component(s) dropped)")

  if (length(flags) == 0) "all assumptions met" else paste(flags, collapse = "; ")
}

# Tweedie GLMM (log link, estimated power 1<p<2) RE selection and singular-RE
# reduction -- same rule; used as the terminal step of the 5-HT cascade.
# Strictly positive responses; the no-RE terminal model is a Tweedie glmmTMB
# with no random part.
fit_select_tweedie <- function(d, respcol, fixed_str, re_cands, label = "") {
  d$.resp <- d[[respcol]]
  d <- d[is.finite(d$.resp) & d$.resp > 0, , drop = FALSE]
  fit_tw <- function(fml) suppressMessages(suppressWarnings(
    glmmTMB::glmmTMB(fml, data = d, family = glmmTMB::tweedie(link = "log"))))

  cand <- lapply(names(re_cands), function(rn) {
    m <- tryCatch(fit_tw(as.formula(paste(".resp ~", fixed_str, "+", re_cands[[rn]]))),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf),
         singular = .is_singular_gamma(m))
  })
  pick <- .pick_start(Filter(Negate(is.null), cand), re_cands)

  fit_fn <- function(re_str) {
    fit_tw(as.formula(if (is.null(re_str)) paste(".resp ~", fixed_str)
                      else paste(".resp ~", fixed_str, "+", re_str)))
  }
  red <- reduce_singular_re(fit_fn, split_re_terms(pick$re), label = label)
  if (is.null(red$model)) return(NULL)
  rr <- .reduction_result(red, pick, "tweedie (no RE)",
                          tryCatch(MuMIn::AICc(red$model), error = function(e) NA_real_))
  c(list(model = red$model, is_lm = FALSE, is_glm = identical(rr$re_retained, "none"),
         singular = FALSE), rr)
}
# Diagnostic figure: Normal Q-Q of residuals + residuals-vs-fitted (PNG).
# Uses Pearson residuals for GLM-family models, raw residuals otherwise.
save_diag_plot <- function(model, title, path, glmm = FALSE) {
  r <- tryCatch(if (glmm) residuals(model, type = "pearson") else residuals(model),
                error = function(e) NULL)
  f <- tryCatch(fitted(model), error = function(e) NULL)
  if (is.null(r) || is.null(f) || length(r) < 3) return(invisible(FALSE))
  tryCatch({
    grDevices::png(path, width = 1200, height = 520, res = 110)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.4, 1))
    stats::qqnorm(r, main = paste0(title, "\nNormal Q-Q (residuals)"), pch = 19, cex = 0.6)
    stats::qqline(r, col = "red", lwd = 2)
    plot(f, r, xlab = "Fitted", ylab = "Residuals", main = "Residuals vs fitted", pch = 19, cex = 0.6)
    graphics::abline(h = 0, col = "red", lty = 2)
    tryCatch(graphics::lines(stats::lowess(f, r), col = "blue", lwd = 2), error = function(e) NULL)
    TRUE
  }, error = function(e) FALSE)
}

################################################################################
###  ROSNER'S GENERALISED ESD OUTLIER DETECTION
################################################################################
# Rosner's test on the raw indicator scale; up to k=ROSNER_K outliers. Returns the
# sample IDs flagged as TRUE outliers (those whose test statistic exceeds the
# critical value), with their values and per-step statistics.
rosner_flag <- function(df, valcol, idcol, k = ROSNER_K) {
  sub <- df[is.finite(df[[valcol]]), c(idcol, valcol)]
  x <- sub[[valcol]]; n <- length(x)
  empty <- data.frame(sample_id = character(0), value = numeric(0),
                      R = numeric(0), lambda = numeric(0), stringsAsFactors = FALSE)
  if (n < 10) return(empty)                         # too few obs for a reliable ESD test
  kk <- max(1L, min(as.integer(k), n - 2L))
  rt <- tryCatch(EnvStats::rosnerTest(x, k = kk, warn = FALSE), error = function(e) NULL)
  if (is.null(rt) || is.null(rt$all.stats)) return(empty)
  st <- rt$all.stats
  flagged <- st$Obs.Num[st$Outlier %in% TRUE]
  if (length(flagged) == 0) return(empty)
  data.frame(
    sample_id = as.character(sub[[idcol]][flagged]),
    value     = as.numeric(sub[[valcol]][flagged]),
    R         = round(as.numeric(st[st$Outlier %in% TRUE, grep("^R\\.", names(st))[1]]), 3),
    lambda    = round(as.numeric(st[st$Outlier %in% TRUE, grep("^lambda", names(st))[1]]), 3),
    stringsAsFactors = FALSE
  )
}

# Per-module outlier detection. Concentrations first; ratios last with the
# component-outlier samples removed before the ratio is tested. Returns:
#   $log     : data.frame(module, indicator, sample_id, value, R, lambda, reason)
#   $exclude : named list  indicator -> character vector of sample IDs to drop
compute_module_outliers <- function(data, id_col, conc_vars, ratio_map, module, k = ROSNER_K) {
  data[[id_col]] <- as.character(data[[id_col]])
  logs <- list(); exclude <- list(); conc_ids <- list()

  for (v in conc_vars) {
    fl <- rosner_flag(data, v, id_col, k)
    conc_ids[[v]] <- fl$sample_id
    exclude[[v]]  <- fl$sample_id
    if (nrow(fl) > 0)
      logs[[length(logs) + 1]] <- data.frame(module = module, indicator = v,
        sample_id = fl$sample_id, value = fl$value, R = fl$R, lambda = fl$lambda,
        reason = "Rosner ESD true outlier", stringsAsFactors = FALSE)
  }

  for (r in names(ratio_map)) {
    comps    <- ratio_map[[r]]
    comp_ids <- unique(unlist(conc_ids[comps]))
    rat_ids  <- as.character(data[[id_col]][is.finite(data[[r]])])
    comp_drop <- intersect(comp_ids, rat_ids)
    # log component-derived exclusions
    for (id in comp_drop) {
      wc  <- comps[vapply(comps, function(cc) id %in% conc_ids[[cc]], logical(1))]
      val <- data[[r]][data[[id_col]] == id][1]
      logs[[length(logs) + 1]] <- data.frame(module = module, indicator = r,
        sample_id = id, value = as.numeric(val), R = NA_real_, lambda = NA_real_,
        reason = paste0("dropped: outlier in component ", paste(wc, collapse = "/")),
        stringsAsFactors = FALSE)
    }
    # Rosner on the ratio AFTER removing component outliers
    sub <- data[!(as.character(data[[id_col]]) %in% comp_ids), ]
    fl  <- rosner_flag(sub, r, id_col, k)
    if (nrow(fl) > 0)
      logs[[length(logs) + 1]] <- data.frame(module = module, indicator = r,
        sample_id = fl$sample_id, value = fl$value, R = fl$R, lambda = fl$lambda,
        reason = "Rosner ESD true outlier", stringsAsFactors = FALSE)
    exclude[[r]] <- union(comp_ids, fl$sample_id)
  }

  log_df <- if (length(logs)) do.call(rbind, logs) else
    data.frame(module = character(0), indicator = character(0), sample_id = character(0),
               value = numeric(0), R = numeric(0), lambda = numeric(0),
               reason = character(0), stringsAsFactors = FALSE)
  list(log = log_df, exclude = exclude)
}

################################################################################
###  PER-INDICATOR ANALYSIS
################################################################################

run_indicator <- function(v, data, re_cands, batch_col = NULL,
                          is_ratio = FALSE, unit_label = "",
                          selection = "best", fallback = "gamma", diag_dir = NULL,
                          id_col = NULL, exclude_ids = character(0)) {
  # selection: "best"     = adopt the transform with the highest residual Shapiro p
  #            "priority" = adopt the FIRST transform reaching Shapiro p >= 0.05 in
  #                         the order raw -> log -> sqrt -> Box-Cox; else best.
  # fallback : "gamma"   = Gamma GLMM (log link) if still non-normal & positive
  #            "tweedie" = Tweedie GLMM (log link) if still non-normal & positive
  #            "none"    = force Gaussian LME with the best transform (no GLM fallback)
  # exclude_ids: Rosner ESD true-outlier sample IDs (in id_col) dropped before fitting.
  cat("\n  [", unit_label, "]", v, "")
  d <- data
  if (!is.null(batch_col)) d$batch <- factor(d[[batch_col]])
  keep <- is.finite(d[[v]]) & !is.na(d$Treatment) & !is.na(d$Group)
  if (!is.null(id_col) && length(exclude_ids) > 0)
    keep <- keep & !(as.character(d[[id_col]]) %in% as.character(exclude_ids))
  d <- d[keep, , drop = FALSE]
  d$Treatment <- droplevels(d$Treatment); d$Group <- droplevels(d$Group)
  n_excluded <- if (!is.null(id_col))
    length(intersect(as.character(data[[id_col]]), as.character(exclude_ids))) else 0L
  n <- nrow(d); vmin <- min(d[[v]]); vmax <- max(d[[v]])
  fixed_str <- "Treatment * Group"

  # ── STEP 1: RE selection on original scale (Gaussian LMM) ──────────────────
  d$.orig <- d[[v]]
  selA <- fit_select(d, ".orig", fixed_str, re_cands)
  sw_raw <- shapiro_residual_p(selA$model)

  # ── STEP 2: SINGLE-TRANSFORMATION RULE (2026-11-12) ─────────────────────────
  # Assumptions are checked on the raw response: Shapiro-Wilk on the residuals of
  # the fixed-effects (OLS) model and Levene's test across the Treatment x Group
  # cells. If either p < 0.05, the response is Box-Cox transformed with the lambda
  # (grid -2 to 2, step 0.1) that maximises the Shapiro-Wilk p of the OLS residuals,
  # as in Bellio et al. (2026); bcnPower with a location gamma estimated by maximum
  # likelihood when the response contains zeros. No other
  # transformation and no GLM-family fallback is used; the same checks are
  # repeated on the transformed response and reported (assumption table).
  sel <- selA; transform_name <- "none"; lambda <- NA_real_
  is_gamma <- FALSE; is_tweedie <- FALSE; sw_trans <- sw_raw
  .grp4 <- interaction(d$Treatment, d$Group, drop = TRUE)
  .ols_sw <- function(yv) { d$.z <- yv
    r <- tryCatch(residuals(lm(as.formula(paste(".z ~", fixed_str)), data = d)), error = function(e) NULL)
    if (is.null(r) || length(r) < 3) NA_real_ else tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_) }
  .lev <- function(yv) tryCatch(car::leveneTest(yv ~ .grp4)$`Pr(>F)`[1], error = function(e) NA_real_)
  ols_sw_raw <- .ols_sw(d[[v]]); lev_raw <- .lev(d[[v]])
  ols_sw_trans <- NA_real_; lev_trans <- NA_real_
  if ((!is.na(ols_sw_raw) && ols_sw_raw < 0.05) || (!is.na(lev_raw) && lev_raw < 0.05)) {
    y <- d[[v]]
    bc <- tryCatch({
      m0 <- lm(as.formula(paste(".orig ~", fixed_str)), data = d)
      lam_grid <- seq(-2, 2, by = 0.1)   # grid as in Bellio et al. (2026): lambda maximising the Shapiro-Wilk p
      grid_best <- function(tf) {
        sw <- vapply(lam_grid, function(l) { yt <- tryCatch(tf(l), error = function(e) NA)
          if (!all(is.finite(yt))) NA_real_ else .ols_sw(yt) }, numeric(1))
        if (all(is.na(sw))) NA_real_ else round(lam_grid[which.max(sw)], 1)
      }
      if (min(y) > 0) {
        lam <- grid_best(function(l) car::bcPower(y, l))
        list(y = car::bcPower(y, lam), lambda = lam, name = sprintf("box-cox (lambda=%.1f)", lam))
      } else {
        cf <- coef(car::powerTransform(m0, family = "bcnPower")); cfv <- if (is.matrix(cf)) cf[1, ] else cf
        gam <- as.numeric(cfv[["gamma"]])   # location by ML (Hawkins & Weisberg 2017); lambda from the grid
        lam <- grid_best(function(l) car::bcnPower(y, lambda = l, gamma = gam))
        list(y = car::bcnPower(y, lambda = lam, gamma = gam), lambda = lam,
             name = sprintf("box-cox (lambda=%.1f, gamma=%.3g)", lam, gam))
      }
    }, error = function(e) NULL)
    if (!is.null(bc) && all(is.finite(bc$y))) {
      d$.try <- bc$y
      sel <- fit_select(d, ".try", fixed_str, re_cands)
      transform_name <- bc$name; lambda <- bc$lambda
      ols_sw_trans <- .ols_sw(bc$y); lev_trans <- .lev(bc$y)
    }
  }
  sw_trans <- shapiro_residual_p(sel$model)
  fallback <- "none"

  # ── STEP 3: GLM-family fallback (only if residuals STILL non-normal & y > 0) ─
  #   fallback = "gamma"   -> Gamma GLMM (log link)
  #   fallback = "tweedie" -> Tweedie GLMM (log link, estimated power)
  #   fallback = "none"    -> keep the best Gaussian LME (forced Gaussian)
  if (!is.na(sw_trans) && sw_trans < 0.05 && min(d[[v]]) > 0 && fallback != "none") {
    if (fallback == "tweedie") {
      cat(" [->Tweedie GLMM]")
      sg <- tryCatch(fit_select_tweedie(d, v, fixed_str, re_cands),
                     error = function(e) { message(" [Tweedie error: ", conditionMessage(e), "]"); NULL })
      if (!is.null(sg)) { sel <- sg; transform_name <- "none"; lambda <- NA_real_; is_tweedie <- TRUE }
    } else if (fallback == "gamma") {
      cat(" [->Gamma GLMM]")
      sg <- tryCatch(fit_select_gamma(d, v, fixed_str, re_cands),
                     error = function(e) { message(" [Gamma error: ", conditionMessage(e), "]"); NULL })
      if (!is.null(sg)) { sel <- sg; transform_name <- "none"; lambda <- NA_real_; is_gamma <- TRUE }
    }
  }

  # ── STEP 4: Inference (family-aware) ──────────────────────────────────────
  model <- sel$model; is_glmfam <- is_gamma || is_tweedie
  # Derived from the fitted object, not from which branch produced it: a Gaussian
  # model with every random effect reduced away is a plain lm.
  is_lm <- !is_glmfam && inherits(model, "lm")
  .indicator_label <- paste0(unit_label, "_", v)
  if (!is.null(sel$aicc_table)) .aicc_log_add(.indicator_label, sel$aicc_table)
  if (length(sel$dropped))   # logged once, for the SELECTED model only
    re_reduction_log_add(.indicator_label,
                         if (is_tweedie) "Tweedie" else if (is_gamma) "Gamma" else "Gaussian",
                         sel$re_attempted, sel$re_retained, sel$dropped, sel$reasons)
  if (is_glmfam) {
    anova_final <- tryCatch(car::Anova(model, type = 3), error = function(e) NULL)
    emm_final   <- tryCatch(emmeans(model, ~ Treatment * Group, type = "response"), error = function(e) NULL)
  } else if (is_lm) {
    anova_final <- tryCatch(car::Anova(model, type = 3), error = function(e) NULL)
    emm_final   <- tryCatch(emmeans(model, ~ Treatment * Group), error = function(e) NULL)
  } else {
    # HANDOFF S5/S7: Kenward-Roger primary, Satterthwaite the named
    # sensitivity (written to anova_satterthwaite.csv).
    anova_final <- finalize_gaussian_inference(model, label = .indicator_label)
    emm_final   <- tryCatch(emmeans(model, ~ Treatment * Group), error = function(e) NULL)
  }
  posthoc_final <- tryCatch(pairs(emm_final, adjust = "tukey"), error = function(e) NULL)
  icc_final <- tryCatch(performance::icc(model)$ICC_adjusted, error = function(e) NA_real_)
  r2_final  <- tryCatch({rv <- MuMIn::r.squaredGLMM(model)
                         c(R2m = as.numeric(rv[1, "R2m"]), R2c = as.numeric(rv[1, "R2c"]))},
                        error = function(e) c(R2m = NA_real_, R2c = NA_real_))

  sw_final <- if (is_glmfam) NA_real_ else shapiro_residual_p(model)

  mtype <- if (is_tweedie && isTRUE(sel$is_glm))   "Tweedie GLM (log link, no RE retained)"
           else if (is_tweedie)                     "Tweedie GLMM (log link)"
           else if (is_gamma && isTRUE(sel$is_glm)) "Gamma GLM (log link, no RE retained)"
           else if (is_gamma)                       "Gamma GLMM (log link)"
           else if (is_lm)                          "Gaussian LM (no RE retained)"
           else                                     "Gaussian LMM"

  # R2m / R2c / delta no longer decide anything (the singular-RE rule replaced the
  # R2-based retention), so these columns are NA by construction.
  re_r2m <- NA_real_; re_r2c <- NA_real_; re_r2_dlt <- NA_real_
  if (length(sel$dropped))
    cat(sprintf(" [RE reduced: dropped %s%s]", paste(sel$dropped, collapse = " + "),
                if (is_lm || isTRUE(sel$is_glm)) " -> no RE" else ""))
  # GLM-family assumption checks (Gamma/Tweedie); flagged in the modelling report.
  glm_flags <- if (is_glmfam)
    gamma_assumption_checks(d, v, isTRUE(sel$is_glm),
                            length(sel$dropped) > 0 && !isTRUE(sel$is_glm)) else ""

  # Diagnostic figure (Normal Q-Q + residuals vs fitted) for the final model.
  if (!is.null(diag_dir)) {
    fn <- file.path(diag_dir,
                    paste0(gsub("[^A-Za-z0-9]+", "_", paste(unit_label, v)), ".png"))
    save_diag_plot(model, paste(unit_label, "-", v, "[", mtype, "]"), fn, glmm = is_glmfam)
  }

  model_row <- tibble(
    variable       = v,
    unit           = unit_label,
    n              = n,
    n_outliers_dropped = n_excluded,
    min_value      = vmin,
    max_value      = vmax,
    is_ratio       = is_ratio,
    model_type     = mtype,
    transformation = transform_name,
    lambda         = if (is.na(lambda)) NA_real_ else round(lambda, 2),
    re_selected    = sel$re_chosen,
    re_formula     = sel$re_formula,
    r2m            = if (is.na(re_r2m))  NA_real_ else round(re_r2m,  3),
    r2c            = if (is.na(re_r2c))  NA_real_ else round(re_r2c,  3),
    re_r2_delta    = if (is.na(re_r2_dlt)) NA_real_ else round(re_r2_dlt, 3),
    sw_raw         = round(sw_raw, 3),
    ols_sw_raw     = ols_sw_raw,
    levene_raw     = lev_raw,
    ols_sw_trans   = ols_sw_trans,
    levene_trans   = lev_trans,
    sw_trans       = if (is_glmfam || transform_name == "none") NA_real_ else round(sw_trans, 3),
    sw_final       = if (is.na(sw_final)) NA_real_ else round(sw_final, 3),
    norm_ok        = is_glmfam || (!is.na(sw_final) && sw_final >= 0.05),
    singular       = sel$singular,
    glm_assumption_flags = glm_flags
  )

  cat("->", mtype, "|", transform_name,
      "| RE:", sel$re_chosen,
      "| SW", round(sw_raw, 3), "->", round(if (is.na(sw_final)) sw_trans else sw_final, 3))
  if (is_glmfam) cat(" | GLM checks:", glm_flags)
  anova_effect_size_table(anova_final, model, mtype, label = .indicator_label,
                          r2m = unname(r2_final["R2m"]), r2c = unname(r2_final["R2c"]),
                          icc = icc_final)
  list(model = model, anova = anova_final, emm = emm_final,
       posthoc = posthoc_final, aicc_tbl = sel$aicc_table, model_row = model_row)
}

run_all <- function(vars, data, re_cands, batch_map = NULL, ratio_vars = character(0),
                    unit_label = "", strategy_map = list(), diag_dir = NULL,
                    id_col = NULL, exclude_map = list()) {
  res <- list()
  for (v in vars) {
    bc <- if (!is.null(batch_map) && v %in% names(batch_map)) batch_map[[v]] else NULL
    strat <- if (!is.null(strategy_map[[v]])) strategy_map[[v]]
             else list(selection = "best", fallback = "gamma")
    ex <- if (!is.null(exclude_map[[v]])) exclude_map[[v]] else character(0)
    res[[v]] <- tryCatch(
      run_indicator(v, data, re_cands, batch_col = bc,
                    is_ratio = v %in% ratio_vars, unit_label = unit_label,
                    selection = strat$selection, fallback = strat$fallback,
                    diag_dir = diag_dir, id_col = id_col, exclude_ids = ex),
      error = function(e) { message("ERROR ", v, ": ", conditionMessage(e)); NULL })
  }
  res
}

################################################################################
###  INDICATOR DEFINITIONS & RUN
################################################################################

PLASMA_VARS <- c("plasma_cortisol", "Glucose_mmol_L", "Lactate_mmol_L")
MONO_VARS   <- c("HT5", "HIAA5", "DA", "DOPAC", "ratio_5HIAA_5HT", "ratio_DOPAC_DA")
MONO_RATIOS <- c("ratio_5HIAA_5HT", "ratio_DOPAC_DA")

# Per-indicator modelling overrides for Mono DM (user-specified):
#   5-HT  : priority cascade raw -> log -> sqrt -> Box-Cox -> Tweedie GLMM
#   ratios: force Gaussian LME with the best (Shapiro-maximising) transform, no GLM
DM_STRATEGY <- list(
  HT5             = list(selection = "priority", fallback = "tweedie"),
  ratio_5HIAA_5HT = list(selection = "best",     fallback = "none"),
  ratio_DOPAC_DA  = list(selection = "best",     fallback = "none")
)

# Concentrations (tested first) vs ratios (tested last, after dropping component
# outliers). Component map: ratio -> its (numerator, denominator) concentrations.
MONO_CONC  <- c("HT5", "HIAA5", "DA", "DOPAC")
MONO_RATIO_MAP <- list(ratio_5HIAA_5HT = c("HIAA5", "HT5"),
                       ratio_DOPAC_DA  = c("DOPAC", "DA"))

# ── Rosner's ESD outlier detection per module (concentrations first, ratios last) ──
cat("\n\n### ROSNER ESD OUTLIER DETECTION (k =", ROSNER_K, ") ###\n")
plasma_out <- compute_module_outliers(plasma,  "Sample",    PLASMA_VARS, list(),         "Plasma",       ROSNER_K)
monoDM_out <- compute_module_outliers(mono_DM, "Sample_ID", MONO_CONC,   MONO_RATIO_MAP, "Monoamine DM", ROSNER_K)
monoRN_out <- compute_module_outliers(mono_RN, "Sample_ID", MONO_CONC,   MONO_RATIO_MAP, "Monoamine RN", ROSNER_K)
# Extra regions: indicators with too few values are set aside before Rosner's test.
.extra_cells <- function(d, v, drop_ids = character(0)) {
  dd <- d[is.finite(d[[v]]) & !(as.character(d$Sample_ID) %in% drop_ids), , drop = FALSE]
  table(dd$Treatment, dd$Group)
}
.extra_ok <- function(d, v, drop_ids = character(0)) all(.extra_cells(d, v, drop_ids) >= MIN_CELL_EXTRA)
.extra_skip_row <- function(k, v, stage, drop_ids = character(0)) {
  tab <- .extra_cells(mono_X[[k]], v, drop_ids)
  data.frame(region = k, variable = v, stage = stage, n = sum(tab), min_cell_n = min(tab),
             stringsAsFactors = FALSE)
}
monoX_vars <- list(); monoX_out <- list(); .extra_skipped <- list()
for (.k in names(mono_X)) {
  .ok <- vapply(MONO_VARS, function(v) .extra_ok(mono_X[[.k]], v), logical(1))
  monoX_vars[[.k]] <- MONO_VARS[.ok]
  for (v in MONO_VARS[!.ok])
    .extra_skipped[[length(.extra_skipped) + 1]] <- .extra_skip_row(.k, v, "before outlier removal")
  .conc  <- intersect(MONO_CONC, monoX_vars[[.k]])
  .rmap  <- MONO_RATIO_MAP[intersect(names(MONO_RATIO_MAP), monoX_vars[[.k]])]
  monoX_out[[.k]] <- compute_module_outliers(mono_X[[.k]], "Sample_ID", .conc, .rmap,
                                             paste("Monoamine", .k), ROSNER_K)
  # the rule is applied again to what is left after outlier removal
  .ok2 <- vapply(monoX_vars[[.k]], function(v)
    .extra_ok(mono_X[[.k]], v, as.character(monoX_out[[.k]]$exclude[[v]])), logical(1))
  for (v in monoX_vars[[.k]][!.ok2])
    .extra_skipped[[length(.extra_skipped) + 1]] <-
      .extra_skip_row(.k, v, "after outlier removal", as.character(monoX_out[[.k]]$exclude[[v]]))
  monoX_vars[[.k]] <- monoX_vars[[.k]][.ok2]
}
write.csv(if (length(.extra_skipped)) bind_rows(.extra_skipped) else
            data.frame(region = character(0), variable = character(0), stage = character(0), n = integer(0), min_cell_n = integer(0)),
          file.path(DIR_ANOVAS, "extra_regions_not_analysed.csv"), row.names = FALSE)
outlier_log <- bind_rows(c(list(plasma_out$log, monoDM_out$log, monoRN_out$log),
                           lapply(monoX_out, `[[`, "log")))
cat("Total logged outlier/exclusion records:", nrow(outlier_log), "\n")
if (nrow(outlier_log) > 0) print(as.data.frame(outlier_log), row.names = FALSE)

cat("\n\n### PLASMA ###")
plasma_res <- run_all(PLASMA_VARS, plasma, RE_PLASMA, batch_map = PLASMA_BATCH,
                      unit_label = "Plasma", diag_dir = DIR_DIAG,
                      id_col = "Sample", exclude_map = plasma_out$exclude)
cat("\n\n### MONOAMINE DM ###")
monoDM_res <- run_all(MONO_VARS, mono_DM, RE_MONO, ratio_vars = MONO_RATIOS,
                      unit_label = "Monoamine DM", strategy_map = DM_STRATEGY,
                      diag_dir = DIR_DIAG, id_col = "Sample_ID",
                      exclude_map = monoDM_out$exclude)
cat("\n\n### MONOAMINE RN ###")
monoRN_res <- run_all(MONO_VARS, mono_RN, RE_MONO, ratio_vars = MONO_RATIOS,
                      unit_label = "Monoamine RN", diag_dir = DIR_DIAG,
                      id_col = "Sample_ID", exclude_map = monoRN_out$exclude)
monoX_res <- list()
for (.k in names(mono_X)) {
  cat("\n\n### MONOAMINE", .k, "###")
  monoX_res[[.k]] <- run_all(monoX_vars[[.k]], mono_X[[.k]], RE_MONO, ratio_vars = MONO_RATIOS,
                             unit_label = paste("Monoamine", .k), diag_dir = DIR_DIAG,
                             id_col = "Sample_ID", exclude_map = monoX_out[[.k]]$exclude)
}
cat("\n")

extract_list <- function(rs, el) lapply(rs, function(r) if (!is.null(r)) r[[el]] else NULL)
mrow <- function(rs) bind_rows(Filter(Negate(is.null), extract_list(rs, "model_row")))

################################################################################
###  OUTPUT HELPERS (style ported from behaviour_models.R)
################################################################################

pretty_indicator_labels <- c(
  plasma_cortisol = "Plasma cortisol", Glucose_mmol_L = "Glucose (mmol/L)",
  Lactate_mmol_L = "Lactate (mmol/L)",
  HT5 = "5-HT", HIAA5 = "5-HIAA", DA = "DA", DOPAC = "DOPAC",
  ratio_5HIAA_5HT = "5-HIAA/5-HT", ratio_DOPAC_DA = "DOPAC/DA"
)
pretty_variable_name <- function(x)
  ifelse(x %in% names(pretty_indicator_labels), pretty_indicator_labels[x], x)

# p floor is p < 0.001 (fmt_p_house(), house_style.R A5; CONVENTIONS.md
# C7/A-9 -- supersedes the old "< 0.0001" floor). to_subscript_statement() is
# sourced from house_style.R (A7); not redefined here.
format_p_statement <- function(p)
  ifelse(is.na(p), "NA", sub("^p ", "", vapply(p, fmt_p_house, character(1))))
format_numeric_columns <- function(tab, digits = 3) tab %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
format_p_columns <- function(tab) {
  pc <- names(tab)[grepl("p|Pr\\(", names(tab), ignore.case = TRUE)]
  for (c in pc) if (is.numeric(tab[[c]]))
    tab[[c]] <- ifelse(is.na(tab[[c]]), NA, sub("^p = ", "", vapply(tab[[c]], fmt_p_house, character(1))))
  tab
}
clean_table_for_markdown <- function(tab) {
  tab <- as.data.frame(tab)
  if ("variable" %in% names(tab)) tab$variable <- pretty_variable_name(tab$variable)
  format_p_columns(format_numeric_columns(tab))
}
make_markdown_table <- function(tab) {
  tab <- clean_table_for_markdown(tab)
  if (nrow(tab) == 0) return("_No rows._")
  # escape literal pipes (e.g. in RE formulas "(1 | batch)") so columns stay aligned
  tab[] <- lapply(tab, function(x) {
    x <- as.character(x); x[is.na(x)] <- ""; gsub("|", "\\|", x, fixed = TRUE) })
  hdr <- paste0("| ", paste(names(tab), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(tab)), collapse = " | "), " |")
  body <- apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}
anova_to_table <- function(a) { tab <- as.data.frame(a); tab$term <- rownames(tab); rownames(tab) <- NULL; relocate(tab, term) }
posthoc_sig <- function(ph) { tab <- as.data.frame(ph); if ("p.value" %in% names(tab)) tab <- tab %>% filter(p.value < 0.05); tab }

.df_to_flex <- function(df, font_size = 8L) {
  df <- clean_table_for_markdown(df)
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  if (nrow(df) == 0) df <- data.frame(` ` = "No rows.", check.names = FALSE)
  ft <- flextable(df); ft <- bold(ft, part = "header"); ft <- bg(ft, part = "header", bg = "#2E75B6")
  ft <- color(ft, part = "header", color = "white"); ft <- fontsize(ft, part = "all", size = font_size)
  ft <- font(ft, part = "all", fontname = "Calibri"); ft <- align(ft, part = "header", align = "center")
  ft <- align(ft, part = "body", align = "left"); set_table_properties(ft, width = 1, layout = "autofit")
}

# Per-module report: model-selection summary + per-indicator AICc / selected model /
# Type III ANOVA / EMMs / significant Tukey contrasts (focal = Treatment×Group).
write_module_reports <- function(md_path, docx_path, title, desc, model_dt, results) {
  anv <- extract_list(results, "anova"); emm <- extract_list(results, "emm")
  pho <- extract_list(results, "posthoc"); aic <- extract_list(results, "aicc_tbl")

  md <- c(paste0("# ", title), "", desc, "", paste0("Generated on: ", Sys.time()), "",
          "---", "", "## Model-selection summary", "", make_markdown_table(model_dt), "")
  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, desc, style = "Normal")
  doc <- body_add_par(doc, "Model-selection summary", style = "heading 2")
  doc <- body_add_flextable(doc, .df_to_flex(model_dt))

  for (v in names(anv)) {
    if (is.null(anv[[v]])) next
    pv <- pretty_variable_name(v)
    md  <- c(md, "", "---", "", paste0("## ", pv), "", paste0("**Variable key:** `", v, "`"), "")
    doc <- body_add_par(doc, pv, style = "heading 2")
    if (!is.null(aic[[v]])) {
      md  <- c(md, "### AICc RE selection", "", make_markdown_table(aic[[v]]), "")
      doc <- body_add_par(doc, "AICc RE selection", style = "heading 3")
      doc <- body_add_flextable(doc, .df_to_flex(aic[[v]]))
    }
    dr <- model_dt %>% filter(variable == v)
    if (nrow(dr) > 0) {
      md  <- c(md, "### Selected model", "", make_markdown_table(dr), "")
      doc <- body_add_par(doc, "Selected model", style = "heading 3")
      doc <- body_add_flextable(doc, .df_to_flex(dr))
    }
    md  <- c(md, "### Type III ANOVA", "", make_markdown_table(anova_to_table(anv[[v]])), "")
    doc <- body_add_par(doc, "Type III ANOVA", style = "heading 3")
    doc <- body_add_flextable(doc, .df_to_flex(anova_to_table(anv[[v]])))
    if (!is.null(emm[[v]])) {
      md  <- c(md, "### Estimated marginal means", "", make_markdown_table(as.data.frame(emm[[v]])), "")
      doc <- body_add_par(doc, "Estimated marginal means", style = "heading 3")
      doc <- body_add_flextable(doc, .df_to_flex(as.data.frame(emm[[v]])))
    }
    if (!is.null(pho[[v]])) {
      md  <- c(md, "### Significant Tukey contrasts (p < 0.05)", "", make_markdown_table(posthoc_sig(pho[[v]])), "")
      doc <- body_add_par(doc, "Significant Tukey contrasts (p < 0.05)", style = "heading 3")
      doc <- body_add_flextable(doc, .df_to_flex(posthoc_sig(pho[[v]])))
    }
  }
  writeLines(md, con = md_path, useBytes = TRUE)
  print(doc, target = docx_path)
}

# ANOVA statement helpers (manuscript-ready phrases) ---------------------------
# Delegates the statistic/p/effect-size string to anova_statement()
# (house_style.R A7) -- every F-row gets eta2p, every chisq-row gets
# none, structurally (HANDOFF.md S6 [INV], CONVENTIONS.md C6).
format_anova_statement <- function(a, model_obj, term, term_label) {
  if (is.null(a)) return(NA_character_)
  tab <- as.data.frame(a); tab$term <- rownames(tab)
  row <- tab[tab$term == term, , drop = FALSE]; if (nrow(row) == 0) return(NA_character_)
  if (all(c("NumDF","DenDF","F value","Pr(>F)") %in% names(row))) {
    # lmerTest Kenward-Roger / Satterthwaite
    anova_statement(term_label, "F", stat = row$`F value`,
                    df1 = row$NumDF, df2 = row$DenDF, p = row$`Pr(>F)`)
  } else if (all(c("Df","F value","Pr(>F)") %in% names(row)) || all(c("Sum Sq","F value","Pr(>F)") %in% names(row))) {
    # car::Anova lm (F test)
    df1 <- if ("Df" %in% names(row)) row$Df else NA_real_
    df2 <- if (!is.null(model_obj)) tryCatch(df.residual(model_obj), error=function(e) NA_real_) else NA_real_
    anova_statement(term_label, "F", stat = row$`F value`, df1 = df1, df2 = df2, p = row$`Pr(>F)`)
  } else if (all(c("Chisq","Df","Pr(>Chisq)") %in% names(row))) {
    # car::Anova Gamma GLMM / glmmTMB (Wald chi-square)
    anova_statement(term_label, "chisq", stat = row$Chisq, df1 = row$Df, p = row$`Pr(>Chisq)`)
  } else NA_character_
}
make_statement_table <- function(anv, mdl, term, term_label) {
  ind <- names(anv)
  tibble(variable = ind, indicator = pretty_variable_name(ind),
         statement = map_chr(ind, function(v)
           tryCatch(format_anova_statement(anv[[v]], mdl[[v]], term, term_label),
                    error = function(e) NA_character_)))
}
write_csv_utf8bom <- function(df, path) {
  tmp <- tempfile(fileext = ".csv"); write.csv(df, tmp, row.names = FALSE, fileEncoding = "UTF-8")
  rc <- readBin(tmp, raw(), file.info(tmp)$size)
  writeBin(c(as.raw(c(0xEF,0xBB,0xBF)), rc), path); unlink(tmp)
}

# Raw means ± SE per Treatment×Group cell (original scale), AFTER Rosner outlier
# removal (per-indicator exclusions so means match the fitted, cleaned data).
se_fun <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
raw_means_block <- function(data, vars, id_col = NULL, exclude_map = list()) {
  out <- list()
  for (v in vars) {
    dd <- data
    ex <- if (!is.null(exclude_map[[v]])) exclude_map[[v]] else character(0)
    if (!is.null(id_col) && length(ex) > 0)
      dd <- dd[!(as.character(dd[[id_col]]) %in% as.character(ex)), , drop = FALSE]
    tab <- dd %>% filter(is.finite(.data[[v]])) %>%
      group_by(Treatment, Group) %>%
      summarise(n = n(), mean = mean(.data[[v]], na.rm = TRUE), SE = se_fun(.data[[v]]), .groups = "drop") %>%
      arrange(Treatment, Group)
    out[[v]] <- tab
  }
  out
}
# fmt3() is sourced from house_style.R (A5) -- local duplicate removed.
raw_means_md <- function(blocks_named) {
  md <- c("# Raw descriptive means and standard errors", "",
          "Means and SE (= sd/√n) on the original measurement scale, by Treatment × Group cell.",
          "", paste0("Generated: ", Sys.time()), "", "---", "")
  for (mod in names(blocks_named)) {
    md <- c(md, paste0("## ", mod), "")
    for (v in names(blocks_named[[mod]])) {
      tab <- blocks_named[[mod]][[v]]
      md <- c(md, paste0("### ", pretty_variable_name(v), "  (`", v, "`)"), "",
              "| Treatment | Group | n | Raw mean | Raw SE |", "| --- | --- | --- | --- | --- |")
      for (i in seq_len(nrow(tab)))
        md <- c(md, sprintf("| %s | %s | %s | %s | %s |", tab$Treatment[i], tab$Group[i],
                            tab$n[i], fmt3(tab$mean[i]), fmt3(tab$SE[i])))
      md <- c(md, "")
    }
    md <- c(md, "---", "")
  }
  md
}

################################################################################
###  WRITE OUTPUTS
################################################################################

cat("\n--- Writing outputs ---\n")

# ── Outlier report (Rosner ESD true outliers + ratio component exclusions) ────
outlier_tbl <- outlier_log
if (nrow(outlier_tbl) > 0)
  outlier_tbl <- outlier_tbl %>%
    mutate(indicator_label = pretty_variable_name(indicator),
           value = round(value, 4)) %>%
    select(module, indicator = indicator_label, sample_id, value, R, lambda, reason)
outlier_md <- c(
  "# Rosner's generalised ESD outlier report",
  "",
  paste0("Test: EnvStats::rosnerTest, k = ", ROSNER_K,
         " (max outliers per indicator), alpha = 0.05, on the raw indicator scale."),
  "Concentrations were tested first; ratios last after dropping their component",
  "concentration outliers. \"Rosner ESD true outlier\" = test statistic R exceeds the",
  "critical value lambda; such samples are logged and dropped before model fitting.",
  "\"dropped: outlier in component …\" = sample removed from a ratio because it was a",
  "true outlier in the ratio's numerator/denominator concentration.",
  "", paste0("Generated: ", Sys.time()), "", "---", "",
  if (nrow(outlier_tbl) == 0) "_No outliers detected._" else make_markdown_table(outlier_tbl)
)
writeLines(outlier_md, file.path(DIR_REPORTS, "outlier_report.md"), useBytes = TRUE)
odoc <- read_docx()
odoc <- body_add_par(odoc, "Rosner's generalised ESD outlier report", style = "heading 1")
odoc <- body_add_par(odoc, paste0(
  "EnvStats::rosnerTest, k = ", ROSNER_K, " (max outliers/indicator), alpha = 0.05, raw scale. ",
  "Concentrations tested first; ratios last after dropping component-concentration outliers. ",
  "True outliers (R > lambda) are dropped before model fitting."), style = "Normal")
odoc <- body_add_flextable(odoc,
  if (nrow(outlier_tbl) == 0) flextable(data.frame(Note = "No outliers detected."))
  else .df_to_flex(outlier_tbl, font_size = 9L))
print(odoc, target = file.path(DIR_REPORTS, "outlier_report.docx"))
cat("Outlier report written (", nrow(outlier_tbl), "records).\n")

plasma_dt <- mrow(plasma_res); monoDM_dt <- mrow(monoDM_res); monoRN_dt <- mrow(monoRN_res)
monoX_dt <- lapply(monoX_res, mrow)
write.csv(bind_rows(c(list(plasma_dt, monoDM_dt, monoRN_dt), monoX_dt)), file.path(DIR_ANOVAS, "model_decisions_assumptions.csv"), row.names = FALSE)

DESC <- paste(
  "Primary analysis: indicator ~ Treatment * Group (Exercise Choice + No Choice).",
  "Rosner's generalised ESD outlier test (k = 3) applied per indicator before fitting;",
  "true outliers dropped (see outlier_report); ratios tested after dropping their",
  "component-concentration outliers. Mixed model (lmer) always fitted; AICc selects the",
  "random-effect structure; a singular structure is reduced (offending component dropped and refit; lm if none remains).",
  "Residual normality by Shapiro-Wilk; if non-normal, Box-Cox lambda chosen to",
  "maximise residual Shapiro p (lambda~0 -> log), adopted only if it improves.",
  "Means reported as raw means ± SE on the original scale (post-outlier-removal)."
)

write_module_reports(file.path(DIR_REPORTS, "GOT_plasma_Treatment_x_Group.md"),
                     file.path(DIR_REPORTS, "GOT_plasma_Treatment_x_Group.docx"),
                     "GOT plasma — Treatment × Group", DESC, plasma_dt, plasma_res)
write_module_reports(file.path(DIR_REPORTS, "GOT_monoamine_DM_Treatment_x_Group.md"),
                     file.path(DIR_REPORTS, "GOT_monoamine_DM_Treatment_x_Group.docx"),
                     "GOT monoamine (DM) — Treatment × Group", DESC, monoDM_dt, monoDM_res)
write_module_reports(file.path(DIR_REPORTS, "GOT_monoamine_RN_Treatment_x_Group.md"),
                     file.path(DIR_REPORTS, "GOT_monoamine_RN_Treatment_x_Group.docx"),
                     "GOT monoamine (RN) — Treatment × Group", DESC, monoRN_dt, monoRN_res)
for (.k in names(monoX_res))
  write_module_reports(file.path(DIR_REPORTS, paste0("GOT_monoamine_", .k, "_Treatment_x_Group.md")),
                       file.path(DIR_REPORTS, paste0("GOT_monoamine_", .k, "_Treatment_x_Group.docx")),
                       paste0("GOT monoamine (", .k, ") — Treatment × Group"), DESC,
                       monoX_dt[[.k]], monoX_res[[.k]])

# Combined model-decision table (markdown + Word) ------------------------------
decision_dt <- bind_rows(c(list(plasma_dt, monoDM_dt, monoRN_dt), monoX_dt))
writeLines(c("# Model-decision summary (all modules)", "",
             paste0("Generated: ", Sys.time()), "", make_markdown_table(decision_dt)),
           file.path(DIR_REPORTS, "model_decision_summary.md"), useBytes = TRUE)
ddoc <- read_docx()
ddoc <- body_add_par(ddoc, "GOT mono/plasma — Model-decision summary", style = "heading 1")
ddoc <- body_add_par(ddoc, DESC, style = "Normal")
ddoc <- body_add_flextable(ddoc, .df_to_flex(decision_dt, font_size = 8L))
print(ddoc, target = file.path(DIR_REPORTS, "model_decision_summary.docx"))

# ANOVA-statement CSVs (UTF-8 BOM): Treatment×Group, Treatment, Group ---------
.x_el <- function(el) do.call(c, c(list(list()), lapply(names(monoX_res), function(k)
  unname(extract_list(monoX_res[[k]], el)))))
.x_keys  <- unlist(lapply(names(monoX_res), function(k) paste0(k, "_", names(monoX_res[[k]]))))
.x_units <- unlist(lapply(names(monoX_res), function(k) rep(paste0("Mono", k), length(monoX_res[[k]]))))
all_anv <- c(extract_list(plasma_res,"anova"), extract_list(monoDM_res,"anova"), extract_list(monoRN_res,"anova"), .x_el("anova"))
all_mdl <- c(extract_list(plasma_res,"model"), extract_list(monoDM_res,"model"), extract_list(monoRN_res,"model"), .x_el("model"))
# disambiguate duplicate indicator names across modules (mono DM vs RN)
unit_of <- c(setNames(rep("Plasma", length(plasma_res)), names(plasma_res)),
             setNames(rep("MonoDM", length(monoDM_res)), paste0("DM_", names(monoDM_res))),
             setNames(rep("MonoRN", length(monoRN_res)), paste0("RN_", names(monoRN_res))),
             setNames(.x_units, .x_keys))
names(all_anv) <- c(names(plasma_res), paste0("DM_", names(monoDM_res)), paste0("RN_", names(monoRN_res)), .x_keys)
names(all_mdl) <- names(all_anv)
mk_stmt <- function(term, lab) {
  ind <- names(all_anv)
  tibble(key = ind, unit = unit_of[ind],
         indicator = pretty_variable_name(sub("^(DM_|RN_|DL_|VV_|POA_)", "", ind)),
         statement = map_chr(ind, function(v)
           tryCatch(format_anova_statement(all_anv[[v]], all_mdl[[v]], term, lab),
                    error = function(e) NA_character_)))
}
################################################################################
###  C12 SENSITIVITY ARTIFACTS  (HANDOFF.md S7, CONVENTIONS.md C12)
################################################################################

if (length(.AICC_LOG))
  write.csv(bind_rows(.AICC_LOG), file.path(DIR_ANOVAS, "aicc_selection.csv"), row.names = FALSE)
# re_reduction.csv -- which design random effects the data did not identify and were
# therefore dropped. Written whether or not anything was dropped, so an empty file
# means "nothing was singular", not "nobody looked".
re_reduction_log_write(DIR_ANOVAS)
if (length(.SATT_LOG))
  write.csv(bind_rows(.SATT_LOG), file.path(DIR_ANOVAS, "anova_satterthwaite.csv"), row.names = FALSE)
if (length(.ANOVA_CSV_LOG))
  write.csv(bind_rows(.ANOVA_CSV_LOG), file.path(DIR_ANOVAS, "anova.csv"), row.names = FALSE)
cat("  wrote aicc_selection.csv / anova_satterthwaite.csv / anova.csv\n")

# design_balance.csv: run over the RECOVERED design keys (trial_uid, date_uid)
# alongside the raw columns they were derived from, so the bare-`Trial`
# pooling defect this run fixes stays visible in the artifact rather than
# disappearing once corrected.
db_rows <- list(
  sensitivity_design_balance(mono_DM, c("trial_uid", "date_uid", "Trial", "user",
                                        "Treatment", "Group"),
                             label = "MonoDM", run_dir = DIR_ANOVAS),
  sensitivity_design_balance(mono_RN, c("trial_uid", "date_uid", "Trial", "user",
                                        "Treatment", "Group"),
                             label = "MonoRN", run_dir = DIR_ANOVAS),
  sensitivity_design_balance(plasma, c("trial_uid", "date_uid", "sampling_date",
                                       "Treatment", "Group"),
                             label = "Plasma", run_dir = DIR_ANOVAS)
)
write.csv(bind_rows(db_rows), file.path(DIR_ANOVAS, "design_balance.csv"), row.names = FALSE)
cat("  wrote design_balance.csv\n")

################################################################################
# sensitivity_forced_re.csv + pseudoreplication_check.csv -- HANDOFF S7.
#
# forced-RE: three declared comparators against the forced design structure --
#   date_dropped     : (1|date_uid) removed; Group loses its whole-plot error
#   bare_trial_index : the (1|Trial) RE the previous vintage used, which pools
#                      fish across all four Treatment x Group cells (mono only)
#   splitplot_df     : design-based denominator df (Group and Treatment on 4)
#                      recomputed from the primary F, so a boundary variance
#                      estimate cannot silently hand out more df than the
#                      design supports
#
# pseudoreplication: the unit here is one fish and treatment was applied at the
# tank, so S2 [DR]'s aggregation question applies -- refit at the trial level
# (mean per trial_uid) and report both betas, both p, the ratio and a
# direction-agree flag.
################################################################################
cat("\n-- forced-RE and aggregation-level sensitivity (endocrine) --\n")

.SPLITPLOT_DF_ENDO <- c(Treatment = 4, Group = 4)
.SPLITPLOT_DF_ENDO["Treatment:Group"] <- 4

.pval_from <- function(an, trm) {
  if (is.null(an) || !trm %in% rownames(an)) return(NA_real_)
  pcol <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p"), names(an))
  if (!length(pcol)) return(NA_real_)
  as.numeric(an[trm, pcol[1]])
}
.model_frame_of <- function(m) if (inherits(m, "glmmTMB")) m$frame else
  tryCatch(m@frame, error = function(e) stats::model.frame(m))
.fixed_formula_of <- function(m)
  paste(deparse(lme4::nobars(stats::formula(m))), collapse = " ")
.beta_p_for <- function(model, term) {
  an <- tryCatch(as.data.frame(car::Anova(model, type = 3)), error = function(e) NULL)
  p <- .pval_from(an, term)
  b <- tryCatch({
    fx <- if (inherits(model, "glmmTMB")) glmmTMB::fixef(model)$cond else
          if (inherits(model, "lm") && !inherits(model, "merMod")) stats::coef(model)
          else lme4::fixef(model)
    unname(fx[grep(paste0("^", term), names(fx))][1])
  }, error = function(e) NA_real_)
  c(as.numeric(b), p)
}

.endo_forced_rows <- list()
.endo_forced_one <- function(label, res_primary, comparators) {
  m0 <- res_primary$model
  # A plain lm fallback has no random part to vary -- recorded as such.
  if (is.null(m0) || !inherits(m0, c("merMod", "lmerModLmerTest", "glmmTMB"))) {
    .endo_forced_rows[[length(.endo_forced_rows) + 1]] <<- data.frame(
      label = label, term = NA_character_, comparator = "not applicable",
      re_primary = "lm (no RE retained -- no random part)",
      re_comparator = NA_character_, p_primary = NA_real_, p_comparator = NA_real_,
      p_crosses_alpha = NA, stringsAsFactors = FALSE)
    return(invisible(NULL))
  }
  d_fit <- tryCatch(.model_frame_of(m0), error = function(e) NULL)
  fx    <- tryCatch(.fixed_formula_of(m0), error = function(e) NULL)
  if (is.null(d_fit) || is.null(fx)) return(invisible(NULL))
  an0 <- tryCatch(as.data.frame(car::Anova(m0, type = 3)), error = function(e) NULL)
  for (cn in names(comparators)) {
    m1 <- tryCatch(lmerTest::lmer(as.formula(paste(fx, "+", comparators[[cn]])),
                                  data = d_fit, REML = TRUE), error = function(e) NULL)
    if (is.null(m1)) next
    an1 <- tryCatch(as.data.frame(car::Anova(m1, type = 3)), error = function(e) NULL)
    for (trm in intersect(rownames(an0), rownames(an1))) {
      if (identical(trm, "(Intercept)")) next
      p0 <- .pval_from(an0, trm); p1 <- .pval_from(an1, trm)
      .endo_forced_rows[[length(.endo_forced_rows) + 1]] <<- data.frame(
        label = label, term = trm, comparator = cn,
        re_primary = paste(vapply(lme4::findbars(stats::formula(m0)),
                                  function(b) paste0("(", paste(deparse(b), collapse = " "), ")"),
                                  character(1)), collapse = " + "),
        re_comparator = comparators[[cn]],
        p_primary = p0, p_comparator = p1,
        p_crosses_alpha = isTRUE((p0 < 0.05) != (p1 < 0.05)),
        stringsAsFactors = FALSE)
    }
  }
  invisible(NULL)
}

for (.nm in names(monoDM_res))
  tryCatch(.endo_forced_one(paste0("Monoamine DM_", .nm), monoDM_res[[.nm]], RE_COMPARATORS_MONO),
           error = function(e) message("forced-RE DM ", .nm, ": ", conditionMessage(e)))
for (.nm in names(monoRN_res))
  tryCatch(.endo_forced_one(paste0("Monoamine RN_", .nm), monoRN_res[[.nm]], RE_COMPARATORS_MONO),
           error = function(e) message("forced-RE RN ", .nm, ": ", conditionMessage(e)))
for (.k in names(monoX_res)) for (.nm in names(monoX_res[[.k]]))
  tryCatch(.endo_forced_one(paste0("Monoamine ", .k, "_", .nm), monoX_res[[.k]][[.nm]], RE_COMPARATORS_MONO),
           error = function(e) message("forced-RE ", .k, " ", .nm, ": ", conditionMessage(e)))
for (.nm in names(plasma_res))
  tryCatch(.endo_forced_one(paste0("Plasma_", .nm), plasma_res[[.nm]], RE_COMPARATORS_PLASMA),
           error = function(e) message("forced-RE plasma ", .nm, ": ", conditionMessage(e)))

# splitplot_df comparator, recomputed from the primary F.
.sp_rows <- list()
if (length(.ANOVA_CSV_LOG)) {
  .aa <- bind_rows(.ANOVA_CSV_LOG)
  .aa <- .aa[.aa$stat_type == "F" & .aa$term %in% names(.SPLITPLOT_DF_ENDO), , drop = FALSE]
  for (.i in seq_len(nrow(.aa))) {
    .d2 <- unname(.SPLITPLOT_DF_ENDO[.aa$term[.i]])
    .psp <- stats::pf(.aa$statistic[.i], .aa$df1[.i], .d2, lower.tail = FALSE)
    .sp_rows[[length(.sp_rows) + 1]] <- data.frame(
      label = .aa$label[.i], term = .aa$term[.i], comparator = "splitplot_df",
      re_primary = "(1 | date_uid) + (1 | trial_uid)",
      re_comparator = "design-based df (6 dates - 2 groups; 12 trials - 6 dates - 2 terms)",
      p_primary = .aa$p[.i], p_comparator = .psp,
      p_crosses_alpha = isTRUE((.aa$p[.i] < 0.05) != (.psp < 0.05)),
      stringsAsFactors = FALSE)
  }
}
.endo_forced_all <- bind_rows(c(.endo_forced_rows, .sp_rows))
if (!length(.endo_forced_all) || !nrow(.endo_forced_all))
  .endo_forced_all <- data.frame(label = NA_character_, term = NA_character_,
                                 comparator = "none computable", re_primary = NA_character_,
                                 re_comparator = NA_character_, p_primary = NA_real_,
                                 p_comparator = NA_real_, p_crosses_alpha = NA)
write.csv(.endo_forced_all, file.path(DIR_ANOVAS, "sensitivity_forced_re.csv"), row.names = FALSE)
cat("  wrote sensitivity_forced_re.csv (", nrow(.endo_forced_all), " rows, ",
    sum(.endo_forced_all$p_crosses_alpha, na.rm = TRUE), " crossing alpha)\n", sep = "")

# pseudoreplication_check.csv: fish level vs trial-aggregated.
.endo_pseudo <- list()
.endo_pseudo_one <- function(stream_label, d, vars, res_list) {
  for (v in intersect(vars, names(d))) {
    m_unit <- res_list[[v]]$model
    if (is.null(m_unit)) next
    d_agg <- d %>%
      group_by(trial_uid, date_uid, Treatment, Group) %>%
      summarise(.y = mean(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
      as.data.frame()
    m_agg <- tryCatch(lmerTest::lmer(.y ~ Treatment * Group + (1 | date_uid),
                                     data = d_agg, REML = TRUE), error = function(e) NULL)
    if (is.null(m_agg)) next
    u <- .beta_p_for(m_unit, "Treatment"); a <- .beta_p_for(m_agg, "Treatment")
    .endo_pseudo[[length(.endo_pseudo) + 1]] <<- data.frame(
      label = paste0(stream_label, "_", v), term = "Treatment",
      comparison = "fish_vs_trial_aggregated",
      beta_unit_level = u[1], p_unit_level = u[2],
      beta_agg_level = a[1], p_agg_level = a[2],
      beta_ratio = if (is.finite(a[1]) && a[1] != 0) u[1] / a[1] else NA_real_,
      direction_agree = isTRUE(sign(u[1]) == sign(a[1])),
      reportable = isTRUE(sign(u[1]) == sign(a[1])) &&
                   isTRUE((u[2] < 0.05) == (a[2] < 0.05)),
      stringsAsFactors = FALSE)
  }
}
tryCatch(.endo_pseudo_one("Monoamine DM", mono_DM, MONO_VARS, monoDM_res),
         error = function(e) message("pseudorep DM: ", conditionMessage(e)))
tryCatch(.endo_pseudo_one("Monoamine RN", mono_RN, MONO_VARS, monoRN_res),
         error = function(e) message("pseudorep RN: ", conditionMessage(e)))
for (.k in names(monoX_res))
  tryCatch(.endo_pseudo_one(paste("Monoamine", .k), mono_X[[.k]], monoX_vars[[.k]], monoX_res[[.k]]),
           error = function(e) message("pseudorep ", .k, ": ", conditionMessage(e)))
tryCatch(.endo_pseudo_one("Plasma", plasma, PLASMA_VARS, plasma_res),
         error = function(e) message("pseudorep plasma: ", conditionMessage(e)))
.endo_pseudo_all <- if (length(.endo_pseudo)) bind_rows(.endo_pseudo) else
  data.frame(label = NA_character_, term = NA_character_, comparison = "not computable",
             beta_unit_level = NA_real_, p_unit_level = NA_real_, beta_agg_level = NA_real_,
             p_agg_level = NA_real_, beta_ratio = NA_real_, direction_agree = NA,
             reportable = NA)
write.csv(.endo_pseudo_all, file.path(DIR_ANOVAS, "pseudoreplication_check.csv"),
          row.names = FALSE)
cat("  wrote pseudoreplication_check.csv (", nrow(.endo_pseudo_all), " rows, ",
    sum(!.endo_pseudo_all$reportable, na.rm = TRUE), " non-reportable)\n", sep = "")

# metric_collinearity.csv: within each region, every other analyte/ratio
# screened (Spearman) against that region's primary turnover ratio -- the
# construct each region's own manuscript panel foregrounds.
.mono_dm_vars <- c("HT5", "HIAA5", "DA", "DOPAC", "ratio_5HIAA_5HT", "ratio_DOPAC_DA")
if ("ratio_DOPAC_DA" %in% names(mono_DM))
  metric_collinearity(mono_DM, "ratio_DOPAC_DA", setdiff(.mono_dm_vars, "ratio_DOPAC_DA"),
                      run_dir = DIR_ANOVAS)

writeLines(as.character(Sys.time()), file.path(OUTPUT_DIR, "_RUN_COMPLETE"))

write_csv_utf8bom(mk_stmt("Treatment:Group", "Treatment × Group"),
                  file.path(DIR_ANOVAS, "anova_stmt_TreatmentxGroup.csv"))
write_csv_utf8bom(mk_stmt("Treatment", "Treatment"),
                  file.path(DIR_ANOVAS, "anova_stmt_Treatment.csv"))
write_csv_utf8bom(mk_stmt("Group", "Group"),
                  file.path(DIR_ANOVAS, "anova_stmt_Group.csv"))

# ── Post-hoc (Tukey) + compact-letter-display per indicator -------------------
# For every indicator the Treatment×Group emmeans are subjected to Tukey pairwise
# post-hoc comparisons; compact-letter-display (CLD) summarises which of the four
# Treatment:Group cells differ. One CSV per indicator is written to:
#   post-hoc/<key>.csv : full Tukey pairwise table (contrast, estimate, t, p, ...)
#   cld/<key>.csv      : one row per Treatment:Group cell with its CLD letter
# Keys match the ANOVA-statement CSVs (plasma as-is, mono DM_/RN_ prefixed).
# Reconstruct a symmetric p-value matrix from a pairs()-on-flat-grid table --
# see behaviour_models.R's export_posthoc_cld() for the identical
# rationale (pairs()/cld() enumerate the same reference grid in the same row
# order; pairs() with no `by=` emits choose(n,2) contrasts in combn(n,2) order).
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
DIR_POSTHOC <- file.path(OUTPUT_DIR, "post-hoc")
DIR_CLD     <- file.path(OUTPUT_DIR, "cld")
for (.d in c(DIR_POSTHOC, DIR_CLD)) dir.create(.d, showWarnings = FALSE, recursive = TRUE)

all_emm <- c(extract_list(plasma_res, "emm"),
             extract_list(monoDM_res, "emm"),
             extract_list(monoRN_res, "emm"), .x_el("emm"))
all_pho <- c(extract_list(plasma_res, "posthoc"),
             extract_list(monoDM_res, "posthoc"),
             extract_list(monoRN_res, "posthoc"), .x_el("posthoc"))
names(all_emm) <- names(all_anv)
names(all_pho) <- names(all_anv)

n_cld <- 0L; n_pho <- 0L
for (v in names(all_anv)) {
  ind_pretty <- pretty_variable_name(sub("^(DM_|RN_|DL_|VV_|POA_)", "", v))
  unit_lab   <- unit_of[[v]]

  # Tukey pairwise post-hoc table
  if (!is.null(all_pho[[v]])) {
    ph_df <- tryCatch(as.data.frame(all_pho[[v]]), error = function(e) NULL)
    if (!is.null(ph_df) && nrow(ph_df) > 0) {
      ph_df <- cbind(key = v, unit = unit_lab, indicator = ind_pretty, ph_df)
      write_csv_utf8bom(ph_df, file.path(DIR_POSTHOC, paste0(v, ".csv")))
      n_pho <- n_pho + 1L
    }
  }

  # CLD letters per Treatment:Group cell (Tukey-adjusted)
  if (!is.null(all_emm[[v]])) {
    cld_df <- tryCatch(
      as.data.frame(multcomp::cld(all_emm[[v]], adjust = "tukey",
                                  Letters = letters, alpha = 0.05)),
      error = function(e) NULL)
    if (!is.null(cld_df) && ".group" %in% names(cld_df)) {
      cld_df$.group <- trimws(cld_df$.group)

      # CONVENTIONS.md C5 four-step pipeline (steps 2-4): MUST run before the
      # .group -> cld rename below, since cld_pipeline() reads `.group` by name.
      if (all(c("Treatment", "Group") %in% names(cld_df))) {
        cld_df$cell <- paste0(as.character(cld_df$Treatment), ":", as.character(cld_df$Group))
        ph_v <- tryCatch(as.data.frame(all_pho[[v]]), error = function(e) NULL)
        pmat <- if (!is.null(ph_v)) .build_pmat_from_pairs(nrow(cld_df), ph_v$p.value) else NULL
        if (!is.null(pmat)) {
          res <- cld_pipeline(cld_df, cell_ids = seq_len(nrow(cld_df)), pmat = pmat,
                              label = v, cell_col = "cell", cell_order = TG_BREAKS_PLOT)
          if (!res$omitted) cld_df <- res$cld_df
        }
      }
      names(cld_df)[names(cld_df) == ".group"] <- "cld"
      cld_df <- cbind(key = v, unit = unit_lab, indicator = ind_pretty, cld_df)
      write_csv_utf8bom(cld_df, file.path(DIR_CLD, paste0(v, ".csv")))
      n_cld <- n_cld + 1L
    }
  }
}
cat(sprintf("post-hoc/ : %d Tukey tables | cld/ : %d CLD tables\n", n_pho, n_cld))

# Raw means ± SE (figure-ready + manuscript means) ----------------------------
rm_blocks <- list(
  "Plasma"        = raw_means_block(plasma,  PLASMA_VARS, "Sample",    plasma_out$exclude),
  "Monoamine DM"  = raw_means_block(mono_DM, MONO_VARS,   "Sample_ID", monoDM_out$exclude),
  "Monoamine RN"  = raw_means_block(mono_RN, MONO_VARS,   "Sample_ID", monoRN_out$exclude)
)
for (.k in names(mono_X))
  rm_blocks[[paste("Monoamine", .k)]] <- raw_means_block(mono_X[[.k]], monoX_vars[[.k]], "Sample_ID",
                                                         monoX_out[[.k]]$exclude)
writeLines(raw_means_md(rm_blocks), file.path(DIR_REPORTS, "raw_means_SE.md"), useBytes = TRUE)
rm_long <- bind_rows(lapply(names(rm_blocks), function(mod)
  bind_rows(lapply(names(rm_blocks[[mod]]), function(v) {
    t <- rm_blocks[[mod]][[v]]; t$module <- mod; t$variable <- v
    t$indicator <- pretty_variable_name(v); t }))))
write_csv(rm_long, file.path(DIR_REPORTS, "raw_means_SE.csv"))

################################################################################
###  SUMMARY
################################################################################
cat("\n\n=== DECISION SUMMARY ===\n"); print(as.data.frame(decision_dt))
cat("\n=== FILES WRITTEN TO", OUTPUT_DIR, "===\n")
cat("datasets/ : plasma.csv, mono_DM.csv, mono_RN.csv\n")
cat("anovas/   : anova_stmt_TreatmentxGroup.csv, anova_stmt_Treatment.csv, anova_stmt_Group.csv\n")
cat("reports/  : outlier_report.(md|docx), GOT_plasma_*, GOT_monoamine_DM_*,\n")
cat("            GOT_monoamine_RN_* (.md+.docx), model_decision_summary.(md|docx),\n")
cat("            raw_means_SE.(md|csv)\n")
cat("diagnostics/ : per-indicator Normal Q-Q + residuals-vs-fitted PNGs\n")
cat("\n=== plasma_monoamine_models.R complete ===\n")
