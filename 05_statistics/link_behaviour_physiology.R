################################################################################
# link_behaviour_physiology.R
# ─────────────────────────────────────────────────────────────────────────────
# Morphometry <-> neurochemistry <-> behaviour correlation module.
#
# Modernises 99_archive_old_versions/correlation_analysis_19.03.R (dead F:/ and
# OneDrive paths, top-level install.packages()) and, for the first time in this
# tree, actually produces the fish-level joins that script only sketched — see
# LOGS.md 2026-08-19 for the audit that found the body-length matcher had never
# run (fixed the same day in 01_6_match_body_length.R).
#
# THREE LAYERS, each reported separately:
#   C1  morphometry <-> neurochemistry   (n up to 48 after dropping the
#       Naloxone/Morphine pharma arm — see GROUP_KEEP; no BL match needed,
#       both sides key off the fish's own Sample_ID from the sampling event)
#   C2  behaviour   <-> neurochemistry   (n up to 45; requires the BL match,
#       gated to reliability_tier %in% c("reliable","conservative"))
#       — reported BOTH as an exploratory grid AND as a small confirmatory
#       subset with its own FDR pool (section 4b), since correcting a
#       hypothesis-driven comparison inside a 200+ test pool would bury it
#   C3  morphometry <-> behaviour        (same gate as C2)
#
# JOIN CHAIN (assembled once, exported to datasets/, then split by layer):
#   ind_interval.csv (latest behaviour_models.R run), filtered to
#   timepoint == 3 (segment 3 — temporally closest to sampling/euthanasia) ->
#   its `fish_id` is already in tracker form (e.g. T5_S3_F1), which is exactly
#   the "tracking FISH ID" key in match_table_export ->
#   match_table_export (latest STEP1_6_output run) -> "sampling FISH ID"
#   (= Sample_ID) + reliability_tier ->
#   GOT glu_lac_cort updated march 2025.xlsx (Sample) — morphometry (length,
#   weight, hematocrit, haemoglobin, sex) AND plasma (cortisol, glucose,
#   lactate); this file was checked against datasets/sampling data GOT.xlsx
#   and carries identical length/weight values for the same Sample_ID, so it
#   is used as the single morphometry+plasma source rather than joining both ->
#   monoamine concentration ALL GOT.xlsx (Sample_ID, pivoted wide per
#   Brain_area DM/RN) — neurochemistry.
#
# STATISTICS — both descriptive and model-based, per layer:
#   Descriptive: Spearman correlation (robust to the skew already documented
#     in these indicators elsewhere in this pipeline), BH/FDR-corrected across
#     each layer's own comparisons (add_bh_column(), sourced from
#     stats_helpers.R — the one genuinely shared helper file in
#     this codebase; everything else here is deliberately self-contained,
#     matching this project's easy-script convention).
#   Model-based: per pair, glmmTMB LMM  y ~ scale(x) + Group*Treatment +
#     (1|trial_id), so a correlation isn't just re-expressing the known design
#     effect. This is a SINGLE random-effect candidate with a singularity
#     fallback to lm(), not the multi-candidate AICc search the primary
#     outcome scripts use — a deliberate, documented scope choice: this module
#     screens ~170 pairs across three layers, and full AICc RE selection per
#     pair would be disproportionate for an exploratory correlation screen.
#     Where raw and design-adjusted association disagree, both are reported.
#
# GUARDS:
#   - Fails loudly (not silently) if no STEP1_6_output run exists.
#   - Asserts and reports row counts at every join; never silently drops rows.
#   - Refuses C2/C3 if fewer than MIN_RELIABLE_FISH pass the reliability gate.
#   - Built-in canary: tracked body_length_cm (from the match table) vs sampled
#     `length` (glu_lac_cort) is the SAME quantity measured two ways — if this
#     pair does not correlate strongly, the match itself is wrong and every
#     other C2/C3 result is void. Checked first, reported first, and the run
#     is flagged prominently in the report if it fails.
#
# Outputs (output/correlations_output/<timestamp>/):
#   datasets/  assembled per-fish table (C1 and the gated C2/C3 table)
#   anovas/    per-pair Type III ANOVA + AICc-style model summary CSVs
#   post-hoc/  per-layer Spearman correlation matrices (raw + BH-adjusted)
#   reports/   markdown + Word report; figures/ heatmaps + scatter panels
################################################################################

suppressMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble)
  library(glmmTMB); library(car); library(emmeans)
  library(ggplot2); library(scales); library(patchwork)
  library(officer); library(flextable)
})

options(contrasts = c("contr.sum", "contr.poly"), stringsAsFactors = FALSE)
select <- dplyr::select
set.seed(20260819)

MIN_RELIABLE_FISH <- 10   # minimum fish surviving the BL-match gate before C2/C3 run at all
BL_GATE_TIERS <- c("reliable", "conservative")

################################################################################
###  PATHS
################################################################################

ROOT <- "."
DATA_DIR <- file.path(ROOT, "datasets")
EASY_OUT <- file.path(ROOT, "output")
EFF_DIR  <- file.path(ROOT, "02_trajectory_pipeline")

source(file.path(EFF_DIR, "helpers/stats_helpers.R"))  # add_bh_column() etc.

# ---- latest behaviour_models.R run (per-fish behaviour) ------------
.behav_runs <- list.dirs(EASY_OUT, recursive = FALSE, full.names = TRUE)
.behav_runs <- .behav_runs[grepl("/[0-9]{8}_[0-9]{6}$", .behav_runs)]
if (length(.behav_runs) == 0) stop("No behaviour run (YYYYMMDD_HHMMSS) found under ", EASY_OUT)
BEHAV_RUN_DIR <- sort(.behav_runs)[length(.behav_runs)]
IND_INTERVAL_CSV <- file.path(BEHAV_RUN_DIR, "datasets", "ind_interval.csv")
if (!file.exists(IND_INTERVAL_CSV)) stop("ind_interval.csv not found in latest behaviour run: ", BEHAV_RUN_DIR)

# ---- latest STEP1_6_output run (BL match table) -----------------------------
STEP1_6_ROOT <- file.path(EFF_DIR, "STEP1_6_output")
if (!dir.exists(STEP1_6_ROOT)) {
  stop(
    "No STEP1_6_output/ found under ", EFF_DIR, " — the body-length matcher has ",
    "never been run. This module requires it for the C2/C3 layers (C1 alone does ",
    "not need it, but this script does not degrade to C1-only silently — re-run ",
    "01_6_match_body_length.R first, or ask for a C1-only variant explicitly)."
  )
}
.step16_runs <- list.dirs(STEP1_6_ROOT, recursive = FALSE, full.names = TRUE)
.step16_runs <- .step16_runs[grepl("STEP1_6_output_[0-9]{8}_[0-9]{6}$", .step16_runs)]
if (length(.step16_runs) == 0) stop("No STEP1_6_output_* run found under ", STEP1_6_ROOT)
STEP1_6_RUN_DIR <- sort(.step16_runs)[length(.step16_runs)]
.match_xlsx <- list.files(STEP1_6_RUN_DIR, pattern = "^match_table_export_.*\\.xlsx$", full.names = TRUE)
if (length(.match_xlsx) == 0) stop("No match_table_export_*.xlsx found in ", STEP1_6_RUN_DIR)
MATCH_XLSX <- sort(.match_xlsx)[length(.match_xlsx)]

MORPH_PLASMA_XLSX <- file.path(DATA_DIR, "GOT glu_lac_cort updated march 2025.xlsx")
MONO_XLSX <- file.path(DATA_DIR, "monoamine concentration ALL GOT.xlsx")

TIMESTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(EASY_OUT, "correlations_output", TIMESTAMP)
DIR_DATASETS <- file.path(OUTPUT_DIR, "datasets")
DIR_ANOVAS   <- file.path(OUTPUT_DIR, "anovas")
DIR_POSTHOC  <- file.path(OUTPUT_DIR, "post-hoc")
DIR_REPORTS  <- file.path(OUTPUT_DIR, "reports")
DIR_FIGURES  <- file.path(OUTPUT_DIR, "figures")
for (.d in c(OUTPUT_DIR, DIR_DATASETS, DIR_ANOVAS, DIR_POSTHOC, DIR_REPORTS, DIR_FIGURES))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

cat("\n=== link_behaviour_physiology.R ===\n")
cat("Behaviour run :", BEHAV_RUN_DIR, "\n")
cat("BL match run  :", STEP1_6_RUN_DIR, "\n")
cat("Output        :", OUTPUT_DIR, "\n\n")

################################################################################
###  SMALL HELPERS (mirrors the report-writing style of the sibling easy-scripts;
###  see plasma_monoamine_models.R for the originals — copied
###  and lightly adapted here rather than sourced, per this codebase's
###  established self-contained-script convention)
################################################################################

write_df <- function(df, name, dir = DIR_DATASETS) {
  path <- file.path(dir, name)
  tryCatch(readr::write_csv(df, path), error = function(e) message("  write_df failed for ", name, ": ", e$message))
  invisible(path)
}

format_p_statement <- function(p)
  ifelse(is.na(p), "NA", ifelse(p < 0.0001, "< 0.0001", paste0("= ", sprintf("%.3f", p))))

format_numeric_columns <- function(tab, digits = 3) tab %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
format_p_columns <- function(tab) {
  pc <- names(tab)[grepl("^p$|p_value|p\\.value|p_raw|p_bh|p_adj", names(tab), ignore.case = TRUE)]
  for (c in pc) if (is.numeric(tab[[c]]))
    tab[[c]] <- ifelse(is.na(tab[[c]]), NA, ifelse(tab[[c]] < 0.0001, "< 0.0001", sprintf("%.3f", tab[[c]])))
  tab
}
clean_table_for_markdown <- function(tab) format_p_columns(format_numeric_columns(as.data.frame(tab)))
make_markdown_table <- function(tab) {
  tab <- clean_table_for_markdown(tab)
  if (nrow(tab) == 0) return("_No rows._")
  tab[] <- lapply(tab, function(x) { x <- as.character(x); x[is.na(x)] <- ""; gsub("|", "\\|", x, fixed = TRUE) })
  hdr <- paste0("| ", paste(names(tab), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(tab)), collapse = " | "), " |")
  body <- apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}
.df_to_flex <- function(df, font_size = 8L) {
  df <- clean_table_for_markdown(df)
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  if (nrow(df) == 0) df <- data.frame(` ` = "No rows.", check.names = FALSE)
  ft <- flextable(df); ft <- bold(ft, part = "header"); ft <- bg(ft, part = "header", bg = "#2E75B6")
  ft <- color(ft, part = "header", color = "white"); ft <- fontsize(ft, part = "all", size = font_size)
  ft <- font(ft, part = "all", fontname = "Calibri"); ft <- align(ft, part = "header", align = "center")
  ft <- align(ft, part = "body", align = "left"); set_table_properties(ft, width = 1, layout = "autofit")
}

################################################################################
###  1) LOAD + ASSEMBLE
################################################################################

cat("--- Loading source data ---\n")

behav_raw <- read_csv(IND_INTERVAL_CSV, show_col_types = FALSE) %>%
  filter(timepoint == 3) %>%   # segment 3: temporally closest to sampling
  transmute(
    tracking_fish_id = as.character(fish_id),
    trial_id = as.integer(trial_id),
    # str_to_title normalises "Exercise choice"/"No choice" (this file's own
    # casing) to "Exercise Choice"/"No Choice", matching GROUP_KEEP below and
    # the monoamine sheet's labels — needed for the label cross-check further
    # down to compare like with like instead of flagging a casing difference
    # as a data mismatch.
    Group = str_to_title(as.character(Group)), Treatment = str_to_title(as.character(Treatment)),
    distance_m, speed_cm_s_mean, accel_cm_s2_mean,
    n_zone_transitions, time_spent_zone_flow_prop, time_moving_TS1_prop,
    mean_nnd_cm, mean_iid_cm_fish
  )
cat("Behaviour (S3 only):", nrow(behav_raw), "fish rows\n")
stopifnot("ind_interval.csv is missing expected columns" = all(
  c("tracking_fish_id", "trial_id", "distance_m") %in% names(behav_raw)))

match_raw <- read_excel(MATCH_XLSX) %>%
  transmute(
    tracking_fish_id = as.character(`tracking FISH ID`),
    Sample_ID = as.character(as.integer(`sampling FISH ID`)),
    reliability_tier = as.character(reliability_tier),
    tracked_body_length_cm = as.numeric(`tracking length (cm)`),
    sampled_length_cm_at_match = as.numeric(`sampling lenght (cm)`)
  )
cat("BL match table:", nrow(match_raw), "rows, tier distribution:\n")
print(table(match_raw$reliability_tier))

morph_plasma_raw <- read_excel(MORPH_PLASMA_XLSX, sheet = "results") %>%
  janitor::clean_names() %>%
  # NOTE: this sheet has a duplicate "Sample" header (once near the start, once
  # near the end next to a second glucose reading); readxl/janitor disambiguate
  # to sample_1 / sample_22. Both hold identical values (verified) — sample_1
  # (the primary Sample_ID, paired with plasma_cortisol/glucose/lactate) is used.
  filter(!is.na(sample_1)) %>%
  transmute(
    Sample_ID = as.character(as.integer(sample_1)),
    length_cm = as.numeric(length),
    weight_g = as.numeric(weight),
    hematocrit_pct = as.numeric(hematocrit_avg),
    haemoglobin = as.numeric(haemoglobin),
    sex = as.character(sex),
    plasma_cortisol = as.numeric(plasma_cortisol),
    Glucose_mmol_L = as.numeric(glucose_mmol_l),
    Lactate_mmol_L = as.numeric(lactate_mmol_l)
  ) %>%
  distinct(Sample_ID, .keep_all = TRUE) %>%
  mutate(condition_factor_K = 100 * weight_g / (length_cm^3))  # Fulton's K
cat("Morphometry/plasma (glu_lac_cort):", nrow(morph_plasma_raw), "unique samples\n")

GROUP_KEEP <- c("Exercise Choice", "No Choice")  # matches plasma_monoamine_models.R's
                                                   # GROUP_KEEP — drops the Naloxone/Morphine pharma
                                                   # arm, which isn't crossed with Treatment the same
                                                   # way and would break the Group*Treatment model term

mono_src <- read_excel(MONO_XLSX, sheet = "sheet1") %>%
  transmute(
    Sample_ID = as.character(as.integer(Sample_ID)),
    Brain_area = trimws(as.character(Brain_area)),
    Group = str_to_title(trimws(as.character(Group))),
    Treatment = str_to_title(trimws(as.character(Treatment))),
    HT5 = as.numeric(`5-HT`), HIAA5 = as.numeric(`5-HIAA`),
    ratio_5HIAA_5HT = as.numeric(`5-HIAA/5-HT`),
    DA = as.numeric(DA), DOPAC = as.numeric(DOPAC),
    # DOPAC/DA is the dopaminergic turnover ratio, the direct analogue of
    # 5-HIAA/5-HT. Included so the dopamine hypotheses in the confirmatory
    # subset below can be framed on turnover (comparable to the serotonin
    # ones) rather than on raw concentration alone.
    ratio_DOPAC_DA = as.numeric(`DOPAC/DA`)
  ) %>%
  filter(Brain_area %in% c("DM", "RN"), Group %in% GROUP_KEEP) %>%
  mutate(Group = factor(Group, levels = GROUP_KEEP), Treatment = factor(Treatment, levels = c("Control", "Exercise")))

mono_raw <- mono_src %>%
  pivot_wider(
    id_cols = c(Sample_ID, Group, Treatment), names_from = Brain_area,
    values_from = c(HT5, HIAA5, ratio_5HIAA_5HT, DA, DOPAC, ratio_DOPAC_DA),
    names_glue = "{Brain_area}_{.value}"
  )
cat("Neurochemistry (pivoted DM/RN, pharma arm dropped):", nrow(mono_raw), "unique samples\n")

################################################################################
###  2) C1 — morphometry <-> neurochemistry (no BL match needed)
################################################################################

c1_data <- morph_plasma_raw %>%
  inner_join(mono_raw, by = "Sample_ID")
cat("\nC1 (morphometry<->neurochemistry) assembled n =", nrow(c1_data),
    " (", nrow(morph_plasma_raw), "morphometry samples inner-joined to",
    nrow(mono_raw), "neurochemistry samples)\n")
write_df(c1_data, "c1_morphometry_neurochemistry.csv")

################################################################################
###  3) C2/C3 — behaviour <-> neurochemistry, morphometry <-> behaviour
###     (both gated on BL-match reliability tier)
################################################################################

bl_gated <- match_raw %>% filter(reliability_tier %in% BL_GATE_TIERS)
cat("\nBL-match gate: keeping tiers {", paste(BL_GATE_TIERS, collapse = ", "),
    "} ->", nrow(bl_gated), "of", nrow(match_raw), "tracked fish\n")

if (nrow(bl_gated) < MIN_RELIABLE_FISH) {
  stop(
    "Only ", nrow(bl_gated), " fish passed the BL-match reliability gate ",
    "(< MIN_RELIABLE_FISH = ", MIN_RELIABLE_FISH, "). Refusing to run C2/C3 — ",
    "the behaviour<->neurochemistry and morphometry<->behaviour layers would ",
    "not be trustworthy at this n. C1 has already been written; re-run after ",
    "investigating why so few fish matched reliably."
  )
}

c23_data <- bl_gated %>%
  inner_join(behav_raw, by = "tracking_fish_id") %>%   # Group/Treatment from HERE (tracking side) are authoritative
  left_join(morph_plasma_raw, by = "Sample_ID") %>%
  left_join(mono_raw %>% select(-Group, -Treatment), by = "Sample_ID")

# Cross-check: the monoamine sheet records its own independent Group/Treatment
# per Sample_ID. If it disagrees with the tracking-side labels for the same
# fish, that's a real data-integrity issue worth surfacing, not silently
# picking one source.
.label_check <- bl_gated %>% inner_join(behav_raw, by = "tracking_fish_id") %>%
  select(Sample_ID, Group_tracking = Group, Treatment_tracking = Treatment) %>%
  inner_join(mono_src %>% distinct(Sample_ID, Group, Treatment), by = "Sample_ID") %>%
  filter(as.character(Group_tracking) != as.character(Group) | as.character(Treatment_tracking) != as.character(Treatment))
if (nrow(.label_check) > 0) {
  cat("\n*** WARNING: ", nrow(.label_check), " fish have mismatched Group/Treatment between the",
      "tracking pipeline and the monoamine sheet for the same Sample_ID: ***\n")
  print(.label_check)
} else {
  cat("Group/Treatment cross-check: tracking-side and monoamine-sheet labels agree for all matched fish.\n")
}

cat("C2/C3 assembled n =", nrow(c23_data), "fish (",
    sum(!is.na(c23_data$DM_HT5) | !is.na(c23_data$RN_HT5)), "with neurochemistry,",
    sum(!is.na(c23_data$length_cm)), "with morphometry)\n")
write_df(c23_data, "c23_behaviour_gated.csv")

# ---- CANARY: tracked body_length_cm vs sampled length — same quantity, two
# measurement methods. If this does not correlate strongly, the match itself
# is wrong and every C2/C3 result below is void. Checked and reported FIRST.
canary <- c23_data %>% filter(!is.na(tracked_body_length_cm), !is.na(length_cm))
canary_cor <- if (nrow(canary) >= 4) cor.test(canary$tracked_body_length_cm, canary$length_cm, method = "spearman") else NULL
CANARY_OK <- !is.null(canary_cor) && !is.na(canary_cor$estimate) && canary_cor$estimate > 0.5 && canary_cor$p.value < 0.05
cat("\n*** CANARY CHECK (tracked vs sampled length, n =", nrow(canary), ") ***\n")
if (!is.null(canary_cor)) {
  cat("    Spearman rho =", round(canary_cor$estimate, 3), ", p", format_p_statement(canary_cor$p.value), "\n")
}
cat("    Match quality:", if (CANARY_OK) "OK — BL match is trustworthy" else "FAILED — DO NOT TRUST C2/C3 RESULTS", "\n\n")

################################################################################
###  4) STATISTICS ENGINE — descriptive (Spearman+BH) and model-based (glmmTMB)
################################################################################

spearman_matrix <- function(data, x_vars, y_vars, layer_label) {
  purrr::map_dfr(x_vars, function(xv) {
    purrr::map_dfr(y_vars, function(yv) {
      d <- data %>% select(x = all_of(xv), y = all_of(yv)) %>% filter(is.finite(x), is.finite(y))
      if (nrow(d) < 4) return(tibble(layer = layer_label, x = xv, y = yv, n = nrow(d), rho = NA_real_, p_raw = NA_real_))
      ct <- suppressWarnings(cor.test(d$x, d$y, method = "spearman"))
      tibble(layer = layer_label, x = xv, y = yv, n = nrow(d), rho = unname(ct$estimate), p_raw = ct$p.value)
    })
  }) %>% add_bh_column(p_col = "p_raw", out_col = "p_BH")
}

# Per-pair LMM: y ~ scale(x) + Group*Treatment + (1|trial_id) when trial_id is
# available, with a singularity/too-few-groups fallback to lm(). Single RE
# candidate by design — see header.
#
# trial_id availability differs by layer: C2/C3 carry the tracking pipeline's
# trial_id (5,6,7,...28) via the behaviour join. C1 (pure sampling-event data,
# morphometry+neurochemistry, no tracking involved) has no equivalent — the
# monoamine sheet's own "Trial" column (values 1-3) was checked and rejected:
# each value spans many different sampling dates, so it is a processing
# batch/round label, not a tank-cohort trial identifier, and using it as
# (1|trial_id) would misrepresent the data. C1 therefore always fits a plain
# lm() with Group*Treatment only — labeled as such, not silently degraded.
model_pair <- function(data, xv, yv, layer_label) {
  has_trial <- "trial_id" %in% names(data)
  d <- data %>%
    transmute(
      x = .data[[xv]], y = .data[[yv]],
      trial_id = if (has_trial) factor(trial_id) else NA_character_,
      Group, Treatment
    ) %>%
    filter(is.finite(x), is.finite(y))
  n_trials <- if (has_trial) dplyr::n_distinct(d$trial_id) else NA_integer_
  if (nrow(d) < 8 || (has_trial && n_trials < 2)) {
    return(tibble(layer = layer_label, x = xv, y = yv, n = nrow(d), n_trials = n_trials,
                  model = "skipped (n<8 or <2 trials)", x_estimate = NA_real_, x_p = NA_real_, singular = NA))
  }
  d$xz <- as.numeric(scale(d$x))

  m <- NULL; used_lm <- !has_trial
  if (has_trial) {
    m <- tryCatch(
      suppressWarnings(glmmTMB(y ~ xz + Group * Treatment + (1 | trial_id), data = d, REML = TRUE)),
      error = function(e) NULL
    )
    is_sing <- if (!is.null(m)) tryCatch(any(abs(unlist(VarCorr(m)$cond)) < 1e-6), error = function(e) TRUE) else TRUE
    used_lm <- is.null(m) || isTRUE(is_sing)
  }
  if (used_lm) {
    m <- tryCatch(lm(y ~ xz + Group * Treatment, data = d), error = function(e) NULL)
  }
  if (is.null(m)) {
    return(tibble(layer = layer_label, x = xv, y = yv, n = nrow(d), n_trials = n_trials,
                  model = "failed to fit", x_estimate = NA_real_, x_p = NA_real_, singular = NA))
  }
  co <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  xrow <- if (!is.null(co) && "xz" %in% rownames(co)) co["xz", ] else NULL
  tibble(
    layer = layer_label, x = xv, y = yv, n = nrow(d), n_trials = n_trials,
    model = if (!has_trial) "lm (no trial identifier available for this layer)"
            else if (used_lm) "lm (RE singular/unfit — trial_id dropped)"
            else "glmmTMB (1|trial_id)",
    x_estimate = if (!is.null(xrow)) unname(xrow[1]) else NA_real_,
    x_p = if (!is.null(xrow)) unname(xrow[ncol(co)]) else NA_real_,
    singular = used_lm
  )
}

run_layer <- function(data, x_vars, y_vars, layer_label) {
  desc <- spearman_matrix(data, x_vars, y_vars, layer_label)
  mdl  <- purrr::map_dfr(seq_len(nrow(desc)), function(i)
    model_pair(data, desc$x[i], desc$y[i], layer_label))
  list(descriptive = desc, model = mdl %>% add_bh_column(p_col = "x_p", out_col = "x_p_BH"))
}

MORPH_VARS <- c("length_cm", "weight_g", "condition_factor_K", "hematocrit_pct", "haemoglobin")
NEURO_VARS <- c("DM_HT5", "DM_HIAA5", "DM_ratio_5HIAA_5HT", "DM_DA", "DM_DOPAC", "DM_ratio_DOPAC_DA",
                "RN_HT5", "RN_HIAA5", "RN_ratio_5HIAA_5HT", "RN_DA", "RN_DOPAC", "RN_ratio_DOPAC_DA",
                "plasma_cortisol", "Glucose_mmol_L", "Lactate_mmol_L")

# BEHAV_VARS: selected 2026-08-19 by a collinearity screen over all 28 primary
# behavioural metrics in ind_interval.csv (S3, n=48; hierarchical clustering on
# 1-|r|).
#
# IMPORTANT — the screen is run on SPEARMAN, not Pearson. A first pass used
# Pearson and produced a materially wrong answer, because the descriptive
# statistics in this module are rank-based: two metrics can look only moderately
# related in Pearson while being perfectly rank-equivalent, in which case they
# are literally the same test here. Concretely:
#   distance_m vs speed_cm_s_mean : Pearson 0.694  BUT  Spearman 1.000
# Under Pearson these looked like two distinct axes; under Spearman they are one.
# Any future revision of this list must re-screen on Spearman for the same reason.
#
# Rank-identical clusters (Spearman r >= 0.999) — exactly ONE member kept from each:
#   {distance_m, speed_cm_s_mean}                                    -> distance_m
#   {time_moving_TS1_s, time_moving_TS1_prop, ...centered}           -> time_moving_TS1_prop
#   {time_spent_zone_calm_s, ..._flow_s, ..._flow_prop, ...centered} -> time_spent_zone_flow_prop
#
# Genuine Spearman singletons (independent of everything at r < 0.80), kept:
#   accel_cm_s2_mean    burst//acceleration intensity. NOTE: this is a Spearman
#                       singleton even though Pearson put it in the speed cluster
#                       at r=1.00 — it carries distinct RANK information, so the
#                       earlier Pearson-based decision to drop it was wrong.
#   n_zone_transitions  zone-switching frequency
#   mean_iid_cm_fish    mean distance to ALL other fish (note r=0.02 vs
#                       mean_nnd_cm despite the similar name — NND is distance to
#                       the single nearest fish, a different construct)
#
# `mean_nnd_cm` is also a Spearman singleton but is deliberately EXCLUDED on data-
# quality grounds, not redundancy: it correlates at Pearson r=0.836 with
# time_spent_zone_missing_s (tracking-occlusion time), which suggests it may be
# measuring tracking quality rather than social spacing. Flagged in LOGS.md
# 2026-08-19 as worth checking in the NND computation itself before it is trusted.
BEHAV_VARS <- c("distance_m", "accel_cm_s2_mean", "n_zone_transitions",
                "time_spent_zone_flow_prop", "time_moving_TS1_prop", "mean_iid_cm_fish")

cat("--- Running C1: morphometry <-> neurochemistry (n =", nrow(c1_data), ") ---\n")
C1 <- run_layer(c1_data, MORPH_VARS, NEURO_VARS, "C1_morphometry_neurochemistry")

cat("--- Running C2: behaviour <-> neurochemistry (n =", nrow(c23_data), ", gated) ---\n")
C2 <- run_layer(c23_data, BEHAV_VARS, NEURO_VARS, "C2_behaviour_neurochemistry")

cat("--- Running C3: morphometry <-> behaviour (n =", nrow(c23_data), ", gated) ---\n")
C3 <- run_layer(c23_data, MORPH_VARS, BEHAV_VARS, "C3_morphometry_behaviour")

all_descriptive <- bind_rows(C1$descriptive, C2$descriptive, C3$descriptive)
all_model       <- bind_rows(C1$model, C2$model, C3$model)
write_df(all_descriptive, "spearman_correlations_all_layers.csv", DIR_POSTHOC)
write_df(all_model, "lmm_trial_adjusted_all_layers.csv", DIR_ANOVAS)

stopifnot(
  "FDR-adjusted p must be >= raw p everywhere" =
    all(is.na(all_descriptive$p_raw) | all_descriptive$p_BH >= all_descriptive$p_raw - 1e-9)
)

################################################################################
###  4b) C2 CONFIRMATORY SUBSET — pre-specified pairs, own FDR pool
################################################################################
# The C2 grid above is exploratory: 6 behaviours x 15 neurochemical markers = 90
# tests, where FDR correction is so punishing that only a very large effect could
# survive. This subset instead pre-specifies a small set of hypotheses drawn from
# the fish neuro-behaviour literature, and corrects them in their OWN (much
# smaller) pool — so a real moderate effect has a fair chance of being detected.
#
# This is a genuinely different inferential claim from the grid, not a re-slice of
# it: these pairs are chosen a priori by hypothesis, so their p-values mean
# something the grid's cannot. Reported as its own clearly-labelled section, and
# NOT merged into the exploratory FDR pool (doing so would destroy the very
# advantage that makes it worth running).
#
# Behaviour variables here are the collinearity-screened survivors (see
# BEHAV_VARS note above); neurochemistry is framed on TURNOVER ratios for both
# serotonin and dopamine, so the two transmitter systems are treated symmetrically.
C2_CONFIRMATORY <- tribble(
  ~hypothesis,                                                    ~neuro,                 ~behav,
  "H1 serotonergic turnover (raphe) ~ social spacing",            "RN_ratio_5HIAA_5HT",   "mean_iid_cm_fish",
  "H1 serotonergic turnover (raphe) ~ flow-zone preference",      "RN_ratio_5HIAA_5HT",   "time_spent_zone_flow_prop",
  "H2 serotonergic turnover (forebrain DM) ~ social spacing",     "DM_ratio_5HIAA_5HT",   "mean_iid_cm_fish",
  "H2 serotonergic turnover (forebrain DM) ~ flow-zone pref.",    "DM_ratio_5HIAA_5HT",   "time_spent_zone_flow_prop",
  "H3 dopaminergic turnover (raphe) ~ locomotor activity",        "RN_ratio_DOPAC_DA",    "distance_m",
  "H4 dopaminergic turnover (forebrain DM) ~ locomotor activity", "DM_ratio_DOPAC_DA",    "distance_m",
  "H5 HPI stress axis ~ flow-zone preference",                    "plasma_cortisol",      "time_spent_zone_flow_prop",
  "H5 HPI stress axis ~ zone-switching restlessness",             "plasma_cortisol",      "n_zone_transitions"
)
# Two revisions vs. the first 12-pair draft, both made after checking the data
# rather than on style grounds:
#
# 1. The dopamine hypotheses (H3/H4) originally tested BOTH distance_m and
#    speed_cm_s_mean. Those are Spearman r = 1.000 (see BEHAV_VARS note above),
#    so the second pair was not an independent test — it was the same test run
#    twice, inflating the correction denominator for nothing. Keeping distance_m
#    only.
#
# 2. An "H6 metabolic validity check" (plasma lactate ~ activity) was dropped
#    entirely. It was included as a mechanistic sanity check — more active fish
#    should show higher lactate — and it FAILS: lactate vs every activity measure
#    is essentially zero (Spearman +0.004 to -0.076, n=42). Most likely the
#    plasma sample reflects handling/euthanasia stress rather than the preceding
#    hour of swimming. Reporting a "null" for a test that cannot detect its own
#    positive control would be misleading, so it is removed rather than reported.
#    Recorded in LOGS.md 2026-08-19 as a data-interpretation caveat for plasma
#    lactate generally, not just for this module.

cat("--- Running C2 CONFIRMATORY subset:", nrow(C2_CONFIRMATORY), "pre-specified pairs ---\n")

.missing_conf <- setdiff(c(C2_CONFIRMATORY$neuro, C2_CONFIRMATORY$behav), names(c23_data))
if (length(.missing_conf) > 0)
  stop("Confirmatory subset references columns absent from c23_data: ", paste(.missing_conf, collapse = ", "))

c2_conf_desc <- purrr::pmap_dfr(C2_CONFIRMATORY, function(hypothesis, neuro, behav) {
  d <- c23_data %>% select(x = all_of(behav), y = all_of(neuro)) %>% filter(is.finite(x), is.finite(y))
  if (nrow(d) < 4) return(tibble(hypothesis = hypothesis, behaviour = behav, neurochemical = neuro,
                                 n = nrow(d), rho = NA_real_, p_raw = NA_real_))
  ct <- suppressWarnings(cor.test(d$x, d$y, method = "spearman"))
  tibble(hypothesis = hypothesis, behaviour = behav, neurochemical = neuro,
         n = nrow(d), rho = unname(ct$estimate), p_raw = ct$p.value)
}) %>% add_bh_column(p_col = "p_raw", out_col = "p_BH")   # own small pool, NOT the grid's ~200

c2_conf_model <- purrr::pmap_dfr(C2_CONFIRMATORY, function(hypothesis, neuro, behav) {
  model_pair(c23_data, behav, neuro, "C2_confirmatory") %>% mutate(hypothesis = hypothesis, .before = 1)
}) %>% add_bh_column(p_col = "x_p", out_col = "x_p_BH")

write_df(c2_conf_desc, "c2_confirmatory_spearman.csv", DIR_POSTHOC)
write_df(c2_conf_model, "c2_confirmatory_lmm.csv", DIR_ANOVAS)

stopifnot(
  "Confirmatory FDR-adjusted p must be >= raw p" =
    all(is.na(c2_conf_desc$p_raw) | c2_conf_desc$p_BH >= c2_conf_desc$p_raw - 1e-9)
)

.n_conf_sig <- sum(!is.na(c2_conf_desc$p_BH) & c2_conf_desc$p_BH < 0.05)
cat("    Confirmatory subset:", .n_conf_sig, "of", nrow(c2_conf_desc),
    "pairs significant after BH correction within the confirmatory pool\n")

################################################################################
###  5) FIGURES
################################################################################

save_heatmap <- function(desc, title, fname) {
  d <- desc %>% mutate(sig = ifelse(!is.na(p_BH) & p_BH < 0.05, "*", ""))
  p <- ggplot(d, aes(x = y, y = x, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(sprintf("%.2f", rho), sig)), size = 2.6) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
    labs(title = title, x = NULL, y = NULL, fill = "Spearman\nrho") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(DIR_FIGURES, fname), p, width = max(6, length(unique(d$y)) * 0.55),
         height = max(4, length(unique(d$x)) * 0.4), dpi = 300, limitsize = FALSE)
}
save_heatmap(C1$descriptive, paste0("C1: Morphometry <-> Neurochemistry (n=", nrow(c1_data), ")"), "C1_heatmap.png")
save_heatmap(C2$descriptive, paste0("C2: Behaviour <-> Neurochemistry (BL-match gated, n=", nrow(c23_data), ")"), "C2_heatmap.png")
save_heatmap(C3$descriptive, paste0("C3: Morphometry <-> Behaviour (BL-match gated, n=", nrow(c23_data), ")"), "C3_heatmap.png")

# Canary scatter — always produced, regardless of pass/fail, since it's the
# single most important sanity check in this module.
if (nrow(canary) >= 4) {
  p_canary <- ggplot(canary, aes(x = length_cm, y = tracked_body_length_cm)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#0072B2", size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00") +
    labs(title = paste0("Canary: sampled vs tracked length (rho=", round(canary_cor$estimate, 3),
                        ", p", format_p_statement(canary_cor$p.value), ")"),
         x = "Sampled length, calipers (cm)", y = "Tracked length, idtracker.ai (cm)") +
    theme_classic(base_size = 10)
  ggsave(file.path(DIR_FIGURES, "canary_length_agreement.png"), p_canary, width = 5, height = 5, dpi = 300)
}

################################################################################
###  6) REPORTS (markdown + Word)
################################################################################

md <- c(
  "# Correlations — Morphometry, Neurochemistry, Behaviour",
  "",
  paste0("_Generated ", Sys.time(), ". Behaviour run: ", basename(BEHAV_RUN_DIR),
        ". BL match run: ", basename(STEP1_6_RUN_DIR), "._"),
  "",
  "## Canary check (read this first)",
  "",
  if (nrow(canary) >= 4) paste0(
    "Tracked body length (idtracker.ai) vs sampled length (calipers), same fish, n = ", nrow(canary),
    ": Spearman rho = ", round(canary_cor$estimate, 3), ", p ", format_p_statement(canary_cor$p.value), ". ",
    "**", if (CANARY_OK) "PASSED — the BL match is trustworthy; C2/C3 below can be relied on."
          else "FAILED — the BL match does not reproduce a known-true correlation. Do NOT trust C2/C3 below.", "**"
  ) else "Canary could not be computed (insufficient overlapping data).",
  "",
  "## Sample sizes",
  "",
  paste0("- C1 (morphometry<->neurochemistry): n = ", nrow(c1_data), " (no BL match required)"),
  paste0("- C2/C3 assembled pool (BL-match gated, tiers {", paste(BL_GATE_TIERS, collapse=", "), "}): n = ", nrow(c23_data),
        " of ", nrow(match_raw), " tracked fish"),
  "",
  paste0("## C2 CONFIRMATORY — ", nrow(c2_conf_desc), " pre-specified pairs (read before the exploratory grids)"),
  "",
  paste0("Six a-priori hypotheses, two pairs each, corrected within their **own** pool of ",
        nrow(c2_conf_desc), " tests — not pooled with the ", nrow(all_descriptive),
        "-test exploratory grid below. These pairs were chosen by hypothesis before looking at ",
        "their results, so unlike the grid their p-values support confirmatory inference."),
  "",
  paste0("Behaviour variables are the survivors of a full collinearity screen over all 28 primary ",
        "behavioural metrics (hierarchical clustering on 1-|r|, n=48): `accel_cm_s2_mean` and ",
        "`mean_nnd_cm` were dropped as near-duplicates of `speed_cm_s_mean` (r = 1.00 and 0.97). ",
        "Both serotonin and dopamine are framed on turnover ratios for symmetry."),
  "",
  paste0("**Result: ", .n_conf_sig, " of ", nrow(c2_conf_desc),
        " pairs significant after BH correction within this pool.**"),
  "",
  make_markdown_table(c2_conf_desc %>% arrange(p_raw)), "",
  "### C2 confirmatory — trial-adjusted (glmmTMB, Group*Treatment controlled)", "",
  make_markdown_table(c2_conf_model %>% arrange(x_p)), "",
  "---", "",
  "## C1 — Morphometry <-> Neurochemistry (exploratory grid)", "",
  make_markdown_table(C1$descriptive %>% arrange(p_BH)), "",
  "### C1 trial-adjusted (glmmTMB, Group*Treatment controlled)", "",
  make_markdown_table(C1$model %>% arrange(x_p_BH)), "",
  "## C2 — Behaviour <-> Neurochemistry (BL-match gated, exploratory grid)", "",
  make_markdown_table(C2$descriptive %>% arrange(p_BH)), "",
  "### C2 trial-adjusted", "",
  make_markdown_table(C2$model %>% arrange(x_p_BH)), "",
  "## C3 — Morphometry <-> Behaviour (BL-match gated)", "",
  make_markdown_table(C3$descriptive %>% arrange(p_BH)), "",
  "### C3 trial-adjusted", "",
  make_markdown_table(C3$model %>% arrange(x_p_BH))
)
writeLines(md, file.path(DIR_REPORTS, "correlations_report.md"))

doc <- read_docx()
doc <- body_add_par(doc, "Correlations — Morphometry, Neurochemistry, Behaviour", style = "heading 1")
doc <- body_add_par(doc, paste0("Generated ", Sys.time()), style = "Normal")
doc <- body_add_par(doc, "Canary check", style = "heading 2")
doc <- body_add_par(doc, if (nrow(canary) >= 4) paste0(
    "Tracked vs sampled length, n=", nrow(canary), ": rho=", round(canary_cor$estimate, 3),
    ", p", format_p_statement(canary_cor$p.value), " — ",
    if (CANARY_OK) "PASSED" else "FAILED, DO NOT TRUST C2/C3"
  ) else "Not computable.", style = "Normal")
doc <- body_add_par(doc, paste0("C2 CONFIRMATORY — ", nrow(c2_conf_desc), " pre-specified pairs"), style = "heading 2")
doc <- body_add_par(doc, paste0(
  "Six a-priori hypotheses, two pairs each, corrected within their own pool of ", nrow(c2_conf_desc),
  " tests — not pooled with the ", nrow(all_descriptive), "-test exploratory grid. Behaviour variables ",
  "are the survivors of a full collinearity screen (accel_cm_s2_mean and mean_nnd_cm dropped as ",
  "near-duplicates of speed_cm_s_mean, r = 1.00 and 0.97). Result: ", .n_conf_sig, " of ",
  nrow(c2_conf_desc), " significant after BH correction within this pool."), style = "Normal")
doc <- body_add_flextable(doc, .df_to_flex(c2_conf_desc %>% arrange(p_raw)))
doc <- body_add_par(doc, "C2 confirmatory — trial-adjusted (glmmTMB)", style = "heading 2")
doc <- body_add_flextable(doc, .df_to_flex(c2_conf_model %>% arrange(x_p)))

for (lbl in c("C1", "C2", "C3")) {
  obj <- get(lbl)
  doc <- body_add_par(doc, paste0(lbl, " — descriptive, exploratory grid (Spearman, BH-adjusted)"), style = "heading 2")
  doc <- body_add_flextable(doc, .df_to_flex(obj$descriptive %>% arrange(p_BH)))
  doc <- body_add_par(doc, paste0(lbl, " — trial-adjusted (glmmTMB)"), style = "heading 2")
  doc <- body_add_flextable(doc, .df_to_flex(obj$model %>% arrange(x_p_BH)))
}
print(doc, target = file.path(DIR_REPORTS, "correlations_report.docx"))

cat("\n=== FILES WRITTEN TO", OUTPUT_DIR, "===\n")
cat("datasets/: c1_morphometry_neurochemistry.csv, c23_behaviour_gated.csv\n")
cat("post-hoc/: c2_confirmatory_spearman.csv (pre-specified pairs, own FDR pool)\n")
cat("anovas/  : c2_confirmatory_lmm.csv\n")
cat("post-hoc/: spearman_correlations_all_layers.csv\n")
cat("anovas/  : lmm_trial_adjusted_all_layers.csv\n")
cat("figures/ : C1/C2/C3 heatmaps, canary_length_agreement.png\n")
cat("reports/ : correlations_report.md, correlations_report.docx\n")
cat("\n=== link_behaviour_physiology.R complete ===\n")
