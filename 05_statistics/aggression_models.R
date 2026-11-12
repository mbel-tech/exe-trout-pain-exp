################################################################################
# aggression_models.R
# ─────────────────────────────────────────────────────────────────────────────
# Single-file port of the GOT aggression pipeline
#   (an earlier set of R scripts, 00_load ... 06_sensitivity, not included)
# rebuilt in the style of behaviour_models.R (timestamped output root,
# datasets/anovas/reports/post-hoc/cld/figures subfolders, Word + CSV reports).
#
# Input (event-level master CSV):
#   data/aggression_annotated_june2026.csv
#
# Analyses ("everything"):
#   R1  count aggression rate  (n_events on factorial grids + exposure offset)
#         R1a trial×segment×side, R1b trial×side, R1c fish×segment×side
#         (attacker+attacked), R1d cross-segment fish×side (attacker+attacked)
#   R2  P(flow side)           binomial logit  (trial/fish × segment/aggregated)
#   R3  crowding → side        binomial logit  (event/fish level)
#   + distribution audit, AICc random-effect selection, descriptive raw figures,
#     sensitivity S1–S5 and leave-one-trial-out stability.
#
# ── Deviations from the source pipeline ───────────────────────────────────────
#  1. COUNT FAMILY SELECTION uses a CORRECTED, RE-ORDERED DECISION TREE
#     (fit_aggression_count(), applied to the R1 count responses only). Each test
#     is applied to the model it actually pertains to, so every branch is live:
#       (1) marginal mean/var/(var/mean) is LOGGED as a descriptor only (it is the
#           UNconditional ratio, not a family decision — see note below).
#       (2) start from a POISSON GLMM; AICc random-effect selection within it.
#       (3) performance::check_overdispersion() on the POISSON fit: if overdispersed
#           (p < .05) → refit Negative-Binomial GLMM (nbinom2) + AICc RE selection.
#       (4) DHARMa::testZeroInflation() on the current best: if zero-inflated
#           (p < .05) → refit Zero-Inflated NB (ziformula = ~1) + AICc RE selection.
#       (5) fit assessment with DHARMa's actual targets — testUniformity (KS vs
#           uniform), dispersion, zero-inflation, outliers — LOGGED, not used to
#           trigger a transform. Final family ∈ {Poisson, NB, ZINB}.
#       (6) inference: car::Anova(type 3) + emmeans pairwise + CLD on the COUNT
#           (response / IRR) scale, identical scale across every R1 response.
#     EVERY decision is recorded in the modelling-choices report (.csv + .docx).
#     WHY this ordering: the marginal var/mean is inflated by the offset (log 600 vs
#     1200), the factorial mean structure and structural zeros, so it always exceeds
#     1 and cannot decide the family; overdispersion is a CONDITIONAL property tested
#     on the Poisson fit; an NB fit absorbs overdispersion so ZINB must be reached via
#     a dedicated zero-inflation test; and DHARMa scaled residuals are UNIFORM (not
#     normal) under a correct model, so Shapiro-normality is the wrong test and the
#     sqrt→Gaussian fallback it forced is removed (it caused singular (1|trial) fits
#     and put emmeans on an incomparable sqrt scale).
#  2. Binomial responses (R2, R3) are kept as binomial (logit) GLMMs exactly as in
#     the source; the count decision tree is logged as "N/A (binomial)" for them.
#  3. One monolithic script + in-script Word/CSV reports replace the modular
#     numbered scripts, the cache/*.rds intermediates, and the Quarto report.
################################################################################


################################################################################
###  LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(data.table)   # fast load + grid construction (source uses data.table)
  library(dplyr)
  library(glmmTMB)      # Poisson / nbinom2 / ZINB / binomial GLMMs
  library(lme4)         # isSingular()
  library(car)          # Type II/III ANOVA
  library(MuMIn)        # AICc()
  library(emmeans)      # estimated marginal means + Tukey
  library(multcomp)     # compact-letter displays (cld)
  library(performance)  # check_overdispersion(), r2()
  library(DHARMa)       # simulateResiduals() + dispersion / zero-infl tests
  library(e1071)        # skewness (distribution audit)
  library(ggplot2)
  library(scales)
  library(officer)      # Word reports
  library(flextable)    # pretty Word tables
})

options(contrasts = c("contr.sum", "contr.poly"))  # sensible Type-III ANOVA


################################################################################
###  CONFIG
################################################################################

# Tolerance for "variance ≈ mean" in decision-tree step (1). var/mean ratios at or
# below this pick Poisson; above it pick Negative Binomial.
DISP_TOL <- 1.05

# Significance threshold used throughout the decision tree.
ALPHA <- 0.05

# DHARMa simulations for residual diagnostics.
N_SIM <- 500


################################################################################
###  PATHS
################################################################################

MASTER_CSV <- file.path(
  "data",
  "aggression_annotated_june2026.csv"
)

TIMESTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(
  "output/aggression_output", TIMESTAMP
)
DIR_DATASETS <- file.path(OUTPUT_DIR, "datasets")
DIR_ANOVAS   <- file.path(OUTPUT_DIR, "anovas")
DIR_REPORTS  <- file.path(OUTPUT_DIR, "reports")
DIR_POSTHOC  <- file.path(OUTPUT_DIR, "post-hoc")
DIR_CLD      <- file.path(OUTPUT_DIR, "cld")
DIR_FIGURES  <- file.path(OUTPUT_DIR, "figures")
DIR_DHARMA   <- file.path(DIR_FIGURES, "dharma")
for (.d in c(OUTPUT_DIR, DIR_DATASETS, DIR_ANOVAS, DIR_REPORTS,
             DIR_POSTHOC, DIR_CLD, DIR_FIGURES, DIR_DHARMA))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# -- House style / sensitivity modules (single source of truth) -------------
GOT_ROOT <- "."
source(file.path(GOT_ROOT, "05_statistics/shared/house_style.R"))
source(file.path(GOT_ROOT, "05_statistics/shared/sensitivity_checks.R"))

cat("\n=== aggression_models.R ===\n")
cat("Input :", MASTER_CSV, "\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("  subfolders: datasets / anovas / reports / post-hoc / cld / figures\n\n")


################################################################################
###  SMALL I/O HELPERS
################################################################################

# Write a data frame as a CSV into one of the output subfolders (UTF-8).
write_df <- function(df, name, dir = OUTPUT_DIR) {
  path <- file.path(dir, name)
  tryCatch(
    data.table::fwrite(as.data.table(df), path),
    error = function(e) message("  write_df failed for ", name, ": ", e$message)
  )
  invisible(path)
}

safe_write_lines <- function(lines, path) {
  tryCatch(writeLines(lines, path),
           error = function(e) message("  write failed ", path, ": ", e$message))
}

# Vectorised number formatter. p < 0.001 (fmt_p_house(), house_style.R
# A5; CONVENTIONS.md C7/A-9) -- supersedes the old "<0.0001" floor. Applied
# per-column below: p-named columns get the p floor, other numeric columns
# (estimates, SEs, counts) keep max-3dp with no floor collapse -- C7's
# descriptives-out-of-scope carve-out.
fmt_num <- function(x, digits = 3L, na_str = "NA") {
  x <- as.numeric(x)
  result <- format(round(x, digits), trim = TRUE, scientific = FALSE)
  result[is.na(x)] <- na_str
  result
}
# fmt_num_col(): dispatches to fmt_p_house() for a column whose NAME marks it
# as a p-value (matches the same "p|Pr(" pattern the behaviour/rosner scripts
# use), fmt_num() otherwise.
fmt_num_col <- function(col_name, x, digits = 3L, na_str = "NA") {
  if (grepl("p|Pr\\(", col_name, ignore.case = TRUE))
    ifelse(is.na(x), na_str, sub("^p = ", "", vapply(x, fmt_p_house, character(1))))
  else fmt_num(x, digits = digits, na_str = na_str)
}

# Convert a data frame to a GitHub-flavored Markdown table (AI-readable).
# Numerics formatted to max 3 dp; p-named columns floored at p < 0.001 via fmt_num_col(); pipes escaped.
df_to_md <- function(df, digits = 3, na = "NA") {
  df <- as.data.frame(df)
  if (is.null(df) || nrow(df) == 0L) return("_(no rows)_\n")
  cells <- as.data.frame(mapply(function(col, nm) {
    if (is.numeric(col)) {
      fmt_num_col(nm, col, digits = digits, na_str = na)
    } else {
      out <- as.character(col); out[is.na(out)] <- na
      gsub("\\|", "\\\\|", out)
    }
  }, df, names(df), SIMPLIFY = FALSE), stringsAsFactors = FALSE, optional = TRUE)
  names(cells) <- names(df)
  header <- paste0("| ", paste(names(cells), collapse = " | "), " |")
  sep    <- paste0("| ", paste(rep("---", ncol(cells)), collapse = " | "), " |")
  body   <- vapply(seq_len(nrow(cells)),
                   function(i) paste0("| ", paste(unlist(cells[i, ]), collapse = " | "), " |"),
                   character(1))
  paste0(paste(c(header, sep, body), collapse = "\n"), "\n")
}

# Flextable styled like behaviour_models.R (blue header, Calibri).
.df_to_flex <- function(df, font_size = 8L) {
  df <- as.data.frame(df)
  df[] <- mapply(function(col, nm) if (is.numeric(col)) fmt_num_col(nm, col) else as.character(col),
                df, names(df), SIMPLIFY = FALSE)
  ft <- flextable::flextable(df)
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::bg(ft, part = "header", bg = "#2E75B6")
  ft <- flextable::color(ft, part = "header", color = "white")
  ft <- flextable::fontsize(ft, part = "all", size = font_size)
  ft <- flextable::font(ft, part = "all", fontname = "Calibri")
  ft <- flextable::align(ft, part = "header", align = "center")
  ft <- flextable::set_table_properties(ft, width = 1, layout = "autofit")
  ft
}


################################################################################
###  DATA IMPORT AND VALIDATION  (ported from 00_load.R)
################################################################################

cat("--- Loading and validating master CSV ---\n")
if (!file.exists(MASTER_CSV)) stop("Master CSV not found: ", MASTER_CSV)

dt <- fread(MASTER_CSV, encoding = "UTF-8")
cat(sprintf("  Loaded %d rows × %d columns\n", nrow(dt), ncol(dt)))

stopifnot(
  "Expected 442 rows"   = nrow(dt) == 442,
  "Expected 54 columns" = ncol(dt) == 54
)

# Type coercions + reference levels (No choice / Control / calm).
dt[, trial   := as.integer(trial)]
dt[, segment := factor(segment, levels = c(2L, 3L))]
dt[, group   := factor(group, levels = c("No choice", "Exercise choice"))]
dt[, treatment_clean := factor(treatment_clean, levels = c("Control", "Exercise"))]
dt[, side    := factor(side, levels = c("calm", "flow"))]
dt[, side_filter_applied     := as.logical(side_filter_applied)]
dt[, attacker_match_accepted := as.logical(attacker_match_accepted)]
dt[, attacked_match_accepted := as.logical(attacked_match_accepted)]
dt[, date := as.Date(date)]
dt[, individuals_in_calm := as.integer(individuals_in_calm)]
dt[, individuals_in_flow := as.integer(individuals_in_flow)]

stopifnot(
  "segment values"       = all(dt$segment %in% c(2L, 3L)),
  "side values"          = all(dt$side %in% c("calm", "flow")),
  "no missing trial"     = all(!is.na(dt$trial)),
  "no missing group"     = all(!is.na(dt$group)),
  "no missing treatment" = all(!is.na(dt$treatment_clean)),
  "valid_fish always 4"  = all(dt$valid_fish_in_window == 4L),
  "all high_intensity"   = all(dt$interaction_type == "high_intensity")
)
cat(sprintf("  Validation passed. Trials: %s\n",
            paste(sort(unique(dt$trial)), collapse = ", ")))
cat(sprintf("  Side distribution: %d calm / %d flow\n",
            sum(dt$side == "calm"), sum(dt$side == "flow")))


################################################################################
###  DATA DERIVATION  (ported from 01_derive.R + crowding z from 01b_transform.R)
################################################################################

cat("\n--- Building derived datasets ---\n")

trial_design <- unique(dt[, .(trial, group, treatment_clean, date)])
stopifnot(nrow(trial_design) == uniqueN(dt$trial))

# 1. events_full — primary inference dataset (all events)
events_full <- copy(dt)
events_full[, fish_attacker := fifelse(
  !is.na(attacker_fish_id_matched_trial),
  attacker_fish_id_matched_trial,
  paste(trial, segment, attacker_id, sep = "_seg"))]
events_full[, fish_attacked := fifelse(
  !is.na(attacked_fish_id_matched_trial),
  attacked_fish_id_matched_trial,
  paste(trial, segment, attacked_id, sep = "_seg"))]
events_full[, match_accepted_both := attacker_match_accepted & attacked_match_accepted]
events_full[, crowding_diff := individuals_in_flow - individuals_in_calm]
events_full[, side_flow := as.integer(side == "flow")]

# 2. events_strict — side filter applied & both IDs accepted
events_strict <- events_full[side_filter_applied == TRUE & match_accepted_both == TRUE]

# 3. events_sensitivity_idmatcher — drop rejected_conflict rows
events_sensitivity_idmatcher <- events_full[
  !attacker_match_quality %in% c("rejected_conflict")]

# Helper: expand a count grid to include structural zeros
expand_count_grid <- function(events, key_cols, count_col = "n_events",
                              design = trial_design) {
  key_vals  <- lapply(key_cols, function(k) sort(unique(events[[k]])))
  names(key_vals) <- key_cols
  grid <- do.call(CJ, c(key_vals, list(sorted = FALSE)))
  if ("trial" %in% key_cols) grid <- merge(grid, design, by = "trial", all.x = TRUE)
  counts <- events[, .(n_events = .N), by = key_cols]
  grid   <- merge(grid, counts, by = key_cols, all.x = TRUE)
  grid[is.na(n_events), n_events := 0L]
  setnames(grid, "n_events", count_col)
  grid
}

# 4. R1a — trial × segment × side  (full grid = n_trials × 2 seg × 2 side, offset log 600)
counts_trial_segment_side <- expand_count_grid(events_full, c("trial", "segment", "side"))
counts_trial_segment_side[, offset_log := log(600)]
stopifnot(nrow(counts_trial_segment_side) ==
          uniqueN(events_full$trial) * uniqueN(events_full$segment) * uniqueN(events_full$side))

# 5. R1b — trial × side  (full grid = n_trials × 2 side, offset log 1200)
counts_trial_side <- expand_count_grid(events_full, c("trial", "side"))
counts_trial_side[, offset_log := log(1200)]
stopifnot(nrow(counts_trial_side) == uniqueN(events_full$trial) * uniqueN(events_full$side))

# 6a/6b. R1c — fish × segment × side (attacker / attacked)
build_fish_seg_grid <- function(id_col, fish_name) {
  fts <- unique(events_full[!is.na(get(id_col)), .(trial, segment, fish = get(id_col))])
  grid <- rbindlist(lapply(c("calm", "flow"), function(s)
    copy(fts)[, side := factor(s, levels = c("calm", "flow"))]))
  grid <- merge(grid, trial_design, by = "trial")
  cnt  <- events_full[!is.na(get(id_col)),
                      .(n_events = .N), by = .(trial, segment, side, fish = get(id_col))]
  out  <- merge(grid, cnt, by = c("trial", "segment", "side", "fish"), all.x = TRUE)
  out[is.na(n_events), n_events := 0L]
  out[, offset_log := log(600)]
  setnames(out, "fish", fish_name)
  out
}
counts_fish_segment_side_attacker <- build_fish_seg_grid(
  "attacker_fish_id_matched_trial", "fish_attacker")
counts_fish_segment_side_attacked <- build_fish_seg_grid(
  "attacked_fish_id_matched_trial", "fish_attacked")

# 7a/7b. R1d — cross-segment fish × side (attacker / attacked), aggregated
build_fish_side_grid <- function(id_col, fish_name) {
  ftm  <- unique(events_full[!is.na(get(id_col)), .(trial, fish = get(id_col))])
  grid <- rbindlist(lapply(c("calm", "flow"), function(s)
    copy(ftm)[, side := factor(s, levels = c("calm", "flow"))]))
  grid <- merge(grid, trial_design, by = "trial")
  cnt  <- events_full[!is.na(get(id_col)),
                      .(n_events = .N), by = .(fish = get(id_col), side, trial)]
  out  <- merge(grid, cnt, by = c("trial", "fish", "side"), all.x = TRUE)
  out[is.na(n_events), n_events := 0L]
  out[, offset_log := log(1200)]
  setnames(out, "fish", fish_name)
  out
}
counts_fish_side_attacker <- build_fish_side_grid(
  "attacker_fish_id_matched_trial", "fish_attacker")
counts_fish_side_attacked <- build_fish_side_grid(
  "attacked_fish_id_matched_trial", "fish_attacked")

# 8. R2 probability cells (binomial)
make_prob_cell <- function(events, key_cols) {
  agg <- events[, .(n_flow = sum(side == "flow"), n_total = .N), by = key_cols]
  if ("trial" %in% key_cols) agg <- merge(agg, trial_design, by = "trial", all.x = TRUE)
  agg[, n_calm := n_total - n_flow]
  agg
}
probs_trial_segment <- make_prob_cell(events_full, c("trial", "segment"))
probs_trial         <- make_prob_cell(events_full, "trial")
probs_fish_segment  <- events_full[!is.na(attacker_fish_id_matched_trial),
  .(n_flow = sum(side == "flow"), n_total = .N, trial = first(trial),
    group = first(group), treatment_clean = first(treatment_clean)),
  by = .(fish_attacker = attacker_fish_id_matched_trial, segment)][n_total > 0]
probs_fish <- events_full[!is.na(attacker_fish_id_matched_trial),
  .(n_flow = sum(side == "flow"), n_total = .N, trial = first(trial),
    group = first(group), treatment_clean = first(treatment_clean)),
  by = .(fish_attacker = attacker_fish_id_matched_trial)][n_total > 0]

# 9. R3 crowding datasets (z-standardised crowding_diff, from 01b)
cw_mean <- mean(events_full$crowding_diff)
cw_sd   <- sd(events_full$crowding_diff)

crowding_events <- events_full[, .(
  trial, segment, group, treatment_clean, fish_attacker, side_flow,
  crowding_diff, individuals_in_flow, individuals_in_calm, side_filter_applied)]
crowding_events[, crowding_diff_z := (crowding_diff - cw_mean) / cw_sd]

crowding_fish_segment <- events_full[!is.na(attacker_fish_id_matched_trial), .(
  n_flow_attacks = sum(side == "flow"), n_total = .N,
  mean_crowding_diff = mean(crowding_diff), trial = first(trial),
  group = first(group), treatment_clean = first(treatment_clean)),
  by = .(fish_attacker = attacker_fish_id_matched_trial, segment)][n_total > 0]
crowding_fish_segment[, mean_crowding_diff_z := (mean_crowding_diff - cw_mean) / cw_sd]

crowding_fish <- events_full[!is.na(attacker_fish_id_matched_trial), .(
  n_flow_attacks = sum(side == "flow"), n_total = .N,
  mean_crowding_diff = mean(crowding_diff), trial = first(trial),
  group = first(group), treatment_clean = first(treatment_clean)),
  by = .(fish_attacker = attacker_fish_id_matched_trial)][n_total > 0]
crowding_fish[, mean_crowding_diff_z := (mean_crowding_diff - cw_mean) / cw_sd]

# Persist the analysis-ready datasets.
write_df(events_full,               "events_full.csv",               DIR_DATASETS)
write_df(events_strict,             "events_strict.csv",             DIR_DATASETS)
write_df(counts_trial_segment_side, "counts_trial_segment_side.csv", DIR_DATASETS)
write_df(counts_trial_side,         "counts_trial_side.csv",         DIR_DATASETS)
write_df(counts_fish_side_attacker, "counts_fish_side_attacker.csv", DIR_DATASETS)
write_df(probs_trial,               "probs_trial.csv",               DIR_DATASETS)
write_df(crowding_events,           "crowding_events.csv",           DIR_DATASETS)
cat(sprintf("  events_full=%d  R1a=%d cells  R1b=%d cells  (zeros R1a=%d)\n",
            nrow(events_full), nrow(counts_trial_segment_side), nrow(counts_trial_side),
            sum(counts_trial_segment_side$n_events == 0)))


################################################################################
###  MODEL-FITTING HELPERS  (ported / adapted from helpers/model_utils.R)
################################################################################

# Family fitting wrappers (all glmmTMB; signature shared for aicc_re_select()).
fit_pois_model  <- function(formula, data, reml = FALSE)
  glmmTMB(formula, data = data, family = poisson, REML = reml)
fit_count_model <- function(formula, data, zi = ~0, reml = FALSE)
  glmmTMB(formula, data = data, family = nbinom2, ziformula = zi, REML = reml)
fit_zinb_model  <- function(formula, data, reml = FALSE)
  glmmTMB(formula, data = data, family = nbinom2, ziformula = ~1, REML = reml)
fit_binom_model <- function(formula, data, reml = FALSE)
  glmmTMB(formula, data = data, family = binomial("logit"), REML = reml)
# Singularity check for glmmTMB (near-zero RE variance) or lme4 models.
is_singular_model <- function(m) {
  tryCatch({
    if (inherits(m, "glmmTMB")) {
      vc <- unlist(VarCorr(m))
      any(vc < 1e-6, na.rm = TRUE) || isTRUE(m$fit$convergence != 0)
    } else {
      lme4::isSingular(m)
    }
  }, error = function(e) NA)
}

# DHARMa diagnostics: returns tests + saves the 4-panel plot.
# uniformity_p is the key fit diagnostic (KS test vs uniform distribution).
test_assumptions <- function(model, n_sim = N_SIM, save_path = NULL, label = "") {
  res <- tryCatch(simulateResiduals(model, n = n_sim, plot = FALSE),
                  error = function(e) NULL)
  if (is.null(res)) return(list(residuals = NULL, dispersion = NULL, zero_infl = NULL,
                                uniformity = NULL, uniformity_p = NA_real_, any_sig = NA))
  disp <- tryCatch(testDispersion(res, plot = FALSE),    error = function(e) NULL)
  zinf <- tryCatch(testZeroInflation(res, plot = FALSE), error = function(e) NULL)
  unif <- tryCatch(testUniformity(res, plot = FALSE),    error = function(e) NULL)
  gp <- function(x) if (is.null(x)) NA_real_ else x$p.value
  if (!is.null(save_path)) {
    grDevices::png(save_path, width = 1200, height = 600)
    tryCatch(plot(res, main = label), error = function(e) plot.new())
    grDevices::dev.off()
  }
  list(residuals = res, dispersion = disp, zero_infl = zinf, uniformity = unif,
       uniformity_p = gp(unif),
       any_sig = any(c(gp(disp), gp(zinf), gp(unif)) < ALPHA, na.rm = TRUE))
}

# AICc random-effect selection. Fits every RE candidate with `family_fn`, ranks by
# AICc, and returns the chosen model, preferring:
#   (a) a converged, non-singular candidate within ΔAICc ≤ 2 of the AICc-best;
#   (b) if none is non-singular within that window, the simplest candidate overall.
# Returns list(model, re_chosen, re_formula, delta_best, aicc_table, singular_dropped).
aicc_re_select <- function(fixed_str, re_candidates, data, family_fn, label = "") {
  fits <- lapply(names(re_candidates), function(nm) {
    fml <- stats::as.formula(paste(fixed_str, "+", re_candidates[[nm]]))
    m <- tryCatch(family_fn(fml, data, reml = FALSE),
                  error = function(e) { message("  RE '", nm, "' failed: ",
                                                 conditionMessage(e)); NULL })
    if (is.null(m)) return(NULL)
    conv <- tryCatch(isTRUE(m$fit$convergence == 0), error = function(e) NA)
    sing <- isTRUE(is_singular_model(m))
    list(name = nm, re = re_candidates[[nm]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf),
         converged = conv, singular = sing)
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) { warning("No RE candidate converged for ", label); return(NULL) }

  aicc     <- vapply(fits, `[[`, numeric(1), "aicc")
  conv     <- vapply(fits, function(x) isTRUE(x$converged), logical(1))
  singular <- vapply(fits, function(x) isTRUE(x$singular),  logical(1))
  ok       <- conv & !singular   # fully healthy fits

  best <- which.min(aicc)
  singular_dropped <- FALSE

  # Prefer a healthy candidate within ΔAICc ≤ 2.
  close_ok <- which(ok & (aicc - min(aicc) <= 2))
  if (length(close_ok) > 0) {
    best <- close_ok[which.min(aicc[close_ok])]
  } else if (singular[best] || !conv[best]) {
    # No healthy option within ΔAICc ≤ 2: fall back to simplest candidate
    # (first element of re_candidates, assumed simplest by convention).
    best <- 1L
    singular_dropped <- TRUE
    message("  [!] ", label, ": no non-singular RE within ΔAICc≤2; fell back to simplest RE (",
            fits[[1]]$name, ")")
  }

  tbl <- data.frame(
    RE             = vapply(fits, `[[`, character(1), "name"),
    re_formula     = vapply(fits, `[[`, character(1), "re"),
    AICc           = round(aicc, 2),
    delta_AICc     = round(aicc - min(aicc), 2),
    converged      = conv,
    singular       = singular,
    selected       = seq_along(fits) == best,
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(tbl$AICc), ]
  list(model = fits[[best]]$model, re_chosen = fits[[best]]$name,
       re_formula = fits[[best]]$re, delta_best = round(min(aicc) - sort(aicc)[1], 2),
       aicc_table = tbl, singular_dropped = singular_dropped)
}

# emmeans contrasts as a tidy data frame, response scale.
# Uses "sidak" throughout (emmeans auto-downgrades "tukey" to "sidak" for
# conditional specs like ~ group*treatment | side anyway; be explicit).
extract_emmeans_df <- function(model, specs, type = "response", adjust = "sidak") {
  emm  <- tryCatch(emmeans::emmeans(model, specs = specs, type = type),
                   error = function(e) NULL)
  if (is.null(emm)) return(list(emmeans = NULL, contrasts = NULL, emm_obj = NULL))
  cont <- tryCatch(emmeans::contrast(emm, method = "pairwise", adjust = adjust),
                   error = function(e) NULL)
  cont_df <- if (is.null(cont)) NULL else as.data.frame(summary(cont, infer = TRUE))
  if (!is.null(cont_df)) {
    est_col <- intersect(c("estimate", "ratio", "odds.ratio"), names(cont_df))[1]
    if (!is.na(est_col) && est_col != "estimate")
      names(cont_df)[names(cont_df) == est_col] <- "estimate"
    if ("asymp.LCL" %in% names(cont_df)) {
      names(cont_df)[names(cont_df) == "asymp.LCL"] <- "lower.CL"
      names(cont_df)[names(cont_df) == "asymp.UCL"] <- "upper.CL"
    }
  }
  list(emmeans = as.data.frame(emm), contrasts = cont_df, emm_obj = emm)
}

# Compact-letter display from an emmeans object (sidak-adjusted, consistent with contrasts).
extract_cld <- function(emm_obj) {
  if (is.null(emm_obj)) return(NULL)
  tryCatch(
    as.data.frame(multcomp::cld(emm_obj, adjust = "sidak", Letters = letters,
                                 sort = FALSE)),
    error = function(e) NULL
  )
}

# Leave-one-trial-out stability of fixed effects.
loo_trial_stability <- function(model, data, trial_col = "trial") {
  full_est <- tryCatch(fixef(model)$cond, error = function(e) return(NULL))
  if (is.null(full_est)) return(NULL)
  trials <- unique(data[[trial_col]])
  loo <- lapply(trials, function(t) {
    d <- data[data[[trial_col]] != t, ]
    m <- tryCatch(update(model, data = d), error = function(e) NULL)
    if (is.null(m)) return(rep(NA_real_, length(full_est)))
    tryCatch(fixef(m)$cond, error = function(e) rep(NA_real_, length(full_est)))
  })
  mat <- do.call(rbind, loo)
  data.frame(coef = names(full_est), full_est = full_est,
             min_loo = apply(mat, 2, min, na.rm = TRUE),
             max_loo = apply(mat, 2, max, na.rm = TRUE),
             range_loo = apply(mat, 2, function(x) diff(range(x, na.rm = TRUE))),
             row.names = NULL)
}


################################################################################
###  PLOTTING HELPERS  (ported from helpers/plotting.R)
################################################################################

# theme_aggression() is sourced from house_style.R (A2) -- the local
# theme_bw()-based definition previously here is retired (C8 A-5: one theme
# dialect for the whole chain).
GROUP_COLS   <- c("No choice" = "#56B4E9", "Exercise choice" = "#D55E00")
TREAT_SHAPES <- c("Control" = 17L, "Exercise" = 19L)
SIDE_TREAT_SHAPES <- c("Control.flow" = 17L, "Exercise.flow" = 19L,
                       "Control.calm" = 24L, "Exercise.calm" = 21L)

plot_rate_strip <- function(rates_df, title = "Aggression rate") {
  rates_df <- as.data.table(rates_df)
  ggplot(rates_df, aes(x = group, y = rate, colour = group, shape = treatment_clean)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.42),
               size = 2.5, alpha = 0.85) +
    stat_summary(aes(group = interaction(group, treatment_clean)), fun = mean,
                 geom = "crossbar", width = 0.28, linewidth = 1.6,
                 position = position_dodge(0.42)) +
    stat_summary(aes(group = interaction(group, treatment_clean)), fun.data = mean_se,
                 geom = "errorbar", width = 0.14, linewidth = 0.7,
                 position = position_dodge(0.42)) +
    scale_colour_manual(values = GROUP_COLS, name = "Group") +
    scale_shape_manual(values = TREAT_SHAPES, name = "Treatment") +
    labs(title = title, x = NULL, y = "Aggression acts per min") + theme_aggression()
}

plot_rate_strip_side <- function(rates_df, title = "Aggression rate (by side)") {
  rates_df <- copy(as.data.table(rates_df))
  rates_df[, treat_side := paste(treatment_clean, side, sep = ".")]
  rates_df[, pt_fill := ifelse(side == "calm", "white", GROUP_COLS[as.character(group)])]
  ggplot(rates_df, aes(x = group, y = rate)) +
    geom_point(aes(shape = treat_side, colour = group, fill = pt_fill),
               size = 2.8, stroke = 1.8,
               position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.52)) +
    stat_summary(aes(group = interaction(group, treatment_clean, side), colour = group),
                 fun = mean, geom = "crossbar", width = 0.18, linewidth = 1.3,
                 position = position_dodge(0.52)) +
    stat_summary(aes(group = interaction(group, treatment_clean, side), colour = group),
                 fun.data = mean_se, geom = "errorbar", width = 0.09, linewidth = 0.65,
                 position = position_dodge(0.52)) +
    scale_colour_manual(values = GROUP_COLS, name = "Group") +
    scale_shape_manual(values = SIDE_TREAT_SHAPES, name = "Treatment × Side") +
    scale_fill_identity() +
    labs(title = title, x = NULL, y = "Aggression acts per min") + theme_aggression()
}

plot_side_mosaic <- function(events_df) {
  summ <- as.data.table(events_df)[, .N, by = .(group, treatment_clean, side)]
  summ[, total := sum(N), by = .(group, treatment_clean)]
  summ[, pct := N / total]
  ggplot(summ, aes(x = interaction(group, treatment_clean, sep = "\n"),
                   y = pct, fill = side)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(N, "\n(", round(pct * 100), "%)")),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(calm = "#4E84C4", flow = "#D16103"), name = "Side") +
    scale_y_continuous(labels = percent_format()) +
    labs(title = "Distribution of aggression events by arena side",
         x = "Group × Treatment", y = "Proportion of events (raw)") + theme_aggression()
}

plot_coverage_bars <- function(events_df) {
  summ <- as.data.table(events_df)[, .(n_total = .N, n_side_ok = sum(side_filter_applied),
    n_id_ok = sum(attacker_match_accepted)), by = .(trial, segment)]
  summ[, pct_side := n_side_ok / n_total]; summ[, pct_id := n_id_ok / n_total]
  long <- rbind(
    data.frame(trial = summ$trial, interval = summ$segment,
               metric = "Side annotation OK", pct = summ$pct_side),
    data.frame(trial = summ$trial, interval = summ$segment,
               metric = "Attacker ID matched", pct = summ$pct_id))
  ggplot(long, aes(x = factor(trial), y = pct, fill = metric)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    facet_wrap(~ paste("Interval", interval)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(title = "Data-quality coverage by trial and interval",
         x = "Trial", y = "Proportion of events") + theme_aggression()
}

.parse_contrast <- function(lbl) {
  first <- trimws(sub("(\\s*/\\s*|\\s*-\\s*).*", "", as.character(lbl)))
  first <- gsub("[()]", "", first)
  group <- ifelse(grepl("Exercise choice", first), "Exercise choice", "No choice")
  rest  <- sub("Exercise choice", "", first)
  treatment <- ifelse(grepl("Exercise", rest), "Exercise", "Control")
  data.frame(group = group, treatment_clean = treatment, stringsAsFactors = FALSE)
}

plot_forest <- function(contrasts_df, xlab = "Incidence Rate Ratio", vline = 1,
                        title = NULL) {
  df <- as.data.frame(contrasts_df)
  parsed <- .parse_contrast(df$contrast)
  df$group <- parsed$group; df$treatment_clean <- parsed$treatment_clean
  ggplot(df, aes(x = estimate, xmin = lower.CL, xmax = upper.CL,
                 y = contrast, colour = group)) +
    geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50") +
    geom_errorbar(width = 0.25, linewidth = 0.6, orientation = "y") +
    geom_point(aes(shape = treatment_clean), size = 3.0) +
    scale_x_log10() +
    scale_colour_manual(values = GROUP_COLS, name = "Group") +
    scale_shape_manual(values = TREAT_SHAPES, name = "Treatment") +
    labs(title = title, x = xlab, y = NULL) + theme_aggression()
}


################################################################################
###  MODELLING-CHOICES DECISION LOG
################################################################################

DECISION_LOG <- list()

# Append one fully-specified decision row. Missing fields default to NA.
# Column guide:
#   dispersion_p    — performance::check_overdispersion() on the Poisson fit
#   zero_infl_p     — DHARMa::testZeroInflation() on the Poisson/NB fit
#   zinb_triggered  — was ZINB refitted?
#   dharma_uniformity_p — KS (testUniformity) on final model; final quality check, not a decision
#   final_model     — Poisson | NB | ZINB  (count); Binomial GLMM  (R2/R3)
log_decision <- function(family, level, response, n = NA, mean_y = NA, var_y = NA,
                         ratio = NA, base_family = NA, family_reason = NA,
                         re_selected = NA, re_delta = NA, singular_dropped = NA,
                         dispersion_ratio = NA, dispersion_p = NA,
                         zero_infl_p = NA, zinb_triggered = NA,
                         dharma_uniformity_p = NA,
                         final_model = NA, link = NA, notes = NA) {
  DECISION_LOG[[length(DECISION_LOG) + 1]] <<- data.frame(
    family = family, level = level, response = response, n = n,
    mean = round_num(mean_y), var = round_num(var_y), var_over_mean = round_num(ratio),
    base_family = base_family, family_reason = family_reason,
    re_selected = re_selected, re_delta_AICc = re_delta,
    singular_dropped = singular_dropped,
    dispersion_ratio = round_num(dispersion_ratio),
    dispersion_p = round_num(dispersion_p),
    zero_infl_p = round_num(zero_infl_p), zinb_triggered = zinb_triggered,
    dharma_uniformity_p = round_num(dharma_uniformity_p),
    final_model = final_model, link = link, notes = notes,
    stringsAsFactors = FALSE)
  invisible(NULL)
}
round_num <- function(x) {
  if (!is.numeric(x) || is.na(x)) return(x)
  round(x, 3)   # tiny values are floored to "p < 0.001" by fmt_num_col()/fmt_p_house() at render time
}


################################################################################
###  DECISION-TREE COUNT CORE  (the prescribed modelling logic)
################################################################################

# Fit one R1 count response by the CORRECTED decision tree and run inference.
#   label         short id, e.g. "R1b"
#   level         human-readable scope, e.g. "trial × side (aggregated)"
#   data          the count grid (must contain n_events, offset_log, factors)
#   fixed_str     fixed-effects formula string INCLUDING + offset(offset_log)
#   re_candidates named list of RE strings for AICc selection
#   emm_spec      emmeans specification (formula)
# Returns list(model, family_final, anova, emm, cld, aicc_table).
fit_aggression_count <- function(label, level, data, fixed_str, re_candidates,
                                 emm_spec) {
  cat("\n============================================================\n")
  cat("R1 count response:", label, "—", level, "\n")
  cat("============================================================\n")
  data <- as.data.frame(data)
  y <- data$n_events

  # ── Step 1: marginal descriptors (logged, not used for family decision) ──────
  mu <- mean(y); v <- var(y); ratio <- v / mu
  cat(sprintf("  [1] mean=%.3f var=%.3f ratio=%.3f (marginal descriptor)\n", mu, v, ratio))

  # ── Step 2: start from Poisson; AICc RE selection ─────────────────────────
  sel_pois <- aicc_re_select(fixed_str, re_candidates, data, fit_pois_model,
                             paste0(label, "_POIS"))
  if (is.null(sel_pois)) { warning("R1 ", label, ": Poisson fit failed; skipping."); return(NULL) }
  best <- sel_pois$model; re_sel <- sel_pois$re_chosen
  re_delta <- if (nrow(sel_pois$aicc_table) > 1) sort(sel_pois$aicc_table$delta_AICc)[2] else 0
  sing_dropped <- isTRUE(sel_pois$singular_dropped)
  cat(sprintf("  [2] Poisson AICc RE: %s (next ΔAICc=%.2f%s)\n",
              re_sel, re_delta, if (sing_dropped) "; singular fallback" else ""))
  write_df(sel_pois$aicc_table, sprintf("AICc_RE_%s_POIS.csv", label), DIR_ANOVAS)
  current_family <- "Poisson"

  # ── Step 3: conditional overdispersion on the Poisson fit → NB ─────────────
  od <- tryCatch(performance::check_overdispersion(best), error = function(e) NULL)
  od_ratio <- if (!is.null(od)) as.numeric(od$dispersion_ratio) else NA_real_
  od_p     <- if (!is.null(od)) as.numeric(od$p_value)          else NA_real_
  if (!is.na(od_p) && od_p < ALPHA) {
    cat(sprintf("  [3] overdispersion p=%.4f (ratio=%.2f) → NB refit\n", od_p, od_ratio))
    sel_nb <- aicc_re_select(fixed_str, re_candidates, data, fit_count_model,
                             paste0(label, "_NB"))
    if (!is.null(sel_nb)) {
      best <- sel_nb$model; re_sel <- sel_nb$re_chosen
      sing_dropped <- isTRUE(sel_nb$singular_dropped)
      re_delta <- if (nrow(sel_nb$aicc_table) > 1) sort(sel_nb$aicc_table$delta_AICc)[2] else 0
      current_family <- "NB"
      write_df(sel_nb$aicc_table, sprintf("AICc_RE_%s_NB.csv", label), DIR_ANOVAS)
    }
  } else {
    cat(sprintf("  [3] overdispersion p=%s (ratio=%s) → Poisson adequate\n",
                format(round(od_p, 3)), format(round(od_ratio, 2))))
  }

  # ── Step 4: zero-inflation test on current best → ZINB ───────────────────
  res4 <- tryCatch(simulateResiduals(best, n = N_SIM, plot = FALSE), error = function(e) NULL)
  zi_test <- if (!is.null(res4))
    tryCatch(testZeroInflation(res4, plot = FALSE), error = function(e) NULL) else NULL
  zi_p <- if (!is.null(zi_test)) zi_test$p.value else NA_real_
  zinb_triggered <- FALSE
  if (!is.na(zi_p) && zi_p < ALPHA) {
    zinb_triggered <- TRUE
    cat(sprintf("  [4] zero-inflation p=%.4f → ZINB refit\n", zi_p))
    sel_zi <- aicc_re_select(fixed_str, re_candidates, data, fit_zinb_model,
                             paste0(label, "_ZINB"))
    if (!is.null(sel_zi)) {
      best <- sel_zi$model; re_sel <- sel_zi$re_chosen
      sing_dropped <- isTRUE(sel_zi$singular_dropped)
      re_delta <- if (nrow(sel_zi$aicc_table) > 1) sort(sel_zi$aicc_table$delta_AICc)[2] else 0
      current_family <- "ZINB"
      write_df(sel_zi$aicc_table, sprintf("AICc_RE_%s_ZINB.csv", label), DIR_ANOVAS)
    }
  } else {
    cat(sprintf("  [4] zero-inflation p=%s → no ZINB needed\n", format(round(zi_p, 3))))
  }

  # ── Step 5: final DHARMa fit assessment (logged, no transform) ───────────
  diag <- test_assumptions(best, save_path = file.path(DIR_DHARMA,
                           sprintf("dharma_%s.png", label)), label = label)
  unif_p <- diag$uniformity_p
  cat(sprintf("  [5] final family: %s | DHARMa uniformity p=%s%s\n",
              current_family, format(round(unif_p, 3)),
              if (isTRUE(diag$any_sig)) " [!]" else ""))

  # ── Step 6: inference always on count/IRR scale ──────────────────────────
  an  <- tryCatch(car::Anova(best, type = 3), error = function(e) {
           message("  ANOVA failed: ", e$message); NULL })
  emm <- extract_emmeans_df(best, emm_spec, type = "response")
  cld <- extract_cld(emm$emm_obj)
  if (!is.null(an))  write_df(cbind(term = rownames(an), as.data.frame(an)),
                              sprintf("ANOVA_%s.csv", label), DIR_ANOVAS)
  if (!is.null(emm$contrasts)) write_df(emm$contrasts,
                              sprintf("%s_tukey_pairs.csv", label), DIR_POSTHOC)
  if (!is.null(cld)) write_df(cld, sprintf("%s_cld.csv", label), DIR_CLD)
  # Group × Treatment MARGINAL CLD (averaged over side/segment), used by the
  # manuscript trial-level panels which pool over side. Exported alongside the
  # side-conditioned CLD above.
  gt_emm <- tryCatch(emmeans::emmeans(best, ~ group * treatment_clean, type = "response"),
                     error = function(e) NULL)
  gt_cld <- extract_cld(gt_emm)
  if (!is.null(gt_cld)) write_df(gt_cld, sprintf("%s_grp_trt_cld.csv", label), DIR_CLD)
  # Group × Treatment × INTERVAL (segment) CLD, for models that carry a segment
  # term (e.g. R1a, trial × segment × side). Used by the manuscript Panel B so
  # every group×treatment point per interval gets its own compact letter.
  if ("segment" %in% names(data)) {
    gti_emm <- tryCatch(
      emmeans::emmeans(best, ~ group * treatment_clean * segment, type = "response"),
      error = function(e) NULL)
    gti_cld <- extract_cld(gti_emm)
    if (!is.null(gti_cld))
      write_df(gti_cld, sprintf("%s_grp_trt_int_cld.csv", label), DIR_CLD)
    # Pairwise Tukey contrasts for group × treatment × interval (all 28 pairs
    # across the 8 cells, Tukey-adjusted). Lets reviewers see exact p-values
    # for each cross-interval and within-interval comparison.
    gti_pairs <- tryCatch({
      ct <- emmeans::contrast(gti_emm, method = "pairwise", adjust = "tukey")
      as.data.frame(ct)
    }, error = function(e) NULL)
    if (!is.null(gti_pairs))
      write_df(gti_pairs, sprintf("%s_grp_trt_int_pairs.csv", label), DIR_POSTHOC)
  }
  safe_write_lines(capture.output(summary(best)),
                   file.path(DIR_REPORTS, sprintf("best_model_summary_%s.txt", label)))

  # ── Log the full decision trail ────────────────────────────────────────────
  log_decision(
    family = "R1 (count)", level = level, response = label,
    n = nrow(data), mean_y = mu, var_y = v, ratio = ratio,
    base_family = current_family,
    family_reason = sprintf(
      "Poisson→dispersion p=%s→%s%s",
      format(round(od_p, 3)),
      if (current_family %in% c("NB", "ZINB")) "NB" else "Poisson",
      if (zinb_triggered) sprintf("→zero-infl p=%s→ZINB", format(round(zi_p, 3))) else ""
    ),
    re_selected = re_sel, re_delta = re_delta, singular_dropped = sing_dropped,
    dispersion_ratio = od_ratio, dispersion_p = od_p,
    zero_infl_p = zi_p, zinb_triggered = zinb_triggered,
    dharma_uniformity_p = unif_p,
    final_model = current_family, link = "log",
    notes = if (isTRUE(diag$any_sig)) "DHARMa flagged ≥1 diagnostic" else NA
  )

  list(model = best, family_final = current_family, anova = an, emm = emm,
       cld = cld, aicc_table = sel_pois$aicc_table, label = label)
}

# Fit a binomial R2/R3 response (count tree N/A — logged as such).
fit_aggression_binom <- function(label, family, level, data, formula, emm_spec,
                                  notes = NA) {
  cat("\n--- ", family, " ", label, " — ", level, " ---\n", sep = "")
  m <- tryCatch(fit_binom_model(stats::as.formula(formula), data),
                error = function(e) { message("  fit failed: ", e$message); NULL })
  if (is.null(m)) return(NULL)
  diag <- test_assumptions(m, save_path = file.path(DIR_DHARMA,
                           sprintf("dharma_%s.png", label)), label = label)
  an  <- tryCatch(car::Anova(m, type = 3), error = function(e) NULL)
  emm <- extract_emmeans_df(m, emm_spec, type = "response")
  cld <- extract_cld(emm$emm_obj)
  if (!is.null(an))  write_df(cbind(term = rownames(an), as.data.frame(an)),
                              sprintf("ANOVA_%s.csv", label), DIR_ANOVAS)
  if (!is.null(emm$contrasts)) write_df(emm$contrasts,
                              sprintf("%s_tukey_pairs.csv", label), DIR_POSTHOC)
  if (!is.null(cld)) write_df(cld, sprintf("%s_cld.csv", label), DIR_CLD)

  log_decision(family = family, level = level, response = label, n = nrow(data),
               base_family = "binomial", family_reason = "binary/proportion response",
               dispersion_ratio = if (!is.null(diag$dispersion)) diag$dispersion$statistic else NA,
               dispersion_p     = if (!is.null(diag$dispersion)) diag$dispersion$p.value else NA,
               dharma_uniformity_p = diag$uniformity_p,
               final_model = "Binomial GLMM", link = "logit",
               notes = paste(na.omit(c("count decision tree N/A (binomial)", notes)),
                             collapse = "; "))
  list(model = m, anova = an, emm = emm, cld = cld, label = label)
}


################################################################################
###  ANALYSIS 1 — DISTRIBUTION / FAMILY AUDIT  (ported from 01b_transform.R)
################################################################################

cat("\n\n########## DISTRIBUTION AUDIT ##########\n")

dist_summary <- function(x, label) {
  x <- x[!is.na(x)]
  data.frame(response = label,
             type = ifelse(all(x %% 1 == 0), "count_integer", "continuous"),
             min = min(x), max = max(x), n = length(x), n_zeros = sum(x == 0),
             pct_zeros = 100 * mean(x == 0), mean = mean(x), sd = sd(x),
             skewness = e1071::skewness(x), stringsAsFactors = FALSE)
}

audit_rows <- list()
for (nm in c("R1a_trial_seg_side", "R1b_trial_side")) {
  dat <- if (nm == "R1a_trial_seg_side") counts_trial_segment_side else counts_trial_side
  ds  <- dist_summary(dat$n_events, nm)
  mu <- mean(dat$n_events); v <- var(dat$n_events)
  ds$var_over_mean <- v / mu
  # NOTE: marginal var/mean is a descriptor only (inflated by offset/mean-structure).
  # Family decision is made from conditional overdispersion in fit_aggression_count().
  ds$marginal_family_descriptor <- ifelse(v / mu <= DISP_TOL, "poisson (marginal)", "NB (marginal)")
  audit_rows[[nm]] <- ds
}
for (nm in c("R2a_trial_seg", "R2b_trial")) {
  dat <- if (nm == "R2a_trial_seg") probs_trial_segment else probs_trial
  ds <- dist_summary(dat$n_flow / dat$n_total, nm)
  ds$type <- "proportion_binomial"; ds$var_over_mean <- NA
  ds$marginal_family_descriptor <- "binomial"
  audit_rows[[nm]] <- ds
}
ds_cw <- dist_summary(crowding_events$crowding_diff, "crowding_diff_predictor")
ds_cw$var_over_mean <- NA; ds_cw$marginal_family_descriptor <- "predictor (z-scored)"
audit_rows[["crowding_diff_predictor"]] <- ds_cw

audit_table <- rbindlist(lapply(audit_rows, as.data.table), fill = TRUE)
write_df(audit_table, "transformations_audit.csv", DIR_ANOVAS)
cat("  wrote transformations_audit.csv\n"); print(audit_table)


################################################################################
###  ANALYSIS 2 — RANDOM-EFFECT SELECTION TABLE  (ported from 05_random_effects.R)
################################################################################

cat("\n\n########## RANDOM-EFFECT SELECTION (AICc) ##########\n")

re_report <- list()
run_re_report <- function(fixed_str, re_candidates, data, family_fn, label) {
  sel <- aicc_re_select(fixed_str, re_candidates, data, family_fn, label)
  if (is.null(sel)) return(NULL)
  tbl <- sel$aicc_table; tbl$response <- label; tbl
}
re_report$R1a <- run_re_report(
  "n_events ~ group * treatment_clean * side + segment + offset(offset_log)",
  list(trial_only = "(1 | trial)", trial_nested = "(1 | trial/segment)"),
  counts_trial_segment_side, fit_count_model, "R1a_trial_seg_side")
re_report$R1d_atk <- run_re_report(
  "n_events ~ group * treatment_clean * side + offset(offset_log)",
  list(trial_only = "(1 | trial)",
       trial_fish = "(1 | trial) + (1 | fish_attacker)",
       trial_nested = "(1 | trial/fish_attacker)"),
  counts_fish_side_attacker, fit_count_model, "R1d_attacker")
re_report$R3a <- run_re_report(
  "side_flow ~ crowding_diff_z + group * treatment_clean + segment",
  list(trial_only = "(1 | trial)",
       trial_fish = "(1 | trial) + (1 | fish_attacker)",
       trial_nested = "(1 | trial/fish_attacker)"),
  crowding_events, fit_binom_model, "R3a_event")
re_selection_table <- rbindlist(lapply(Filter(Negate(is.null), re_report),
                                       as.data.table), fill = TRUE)
write_df(re_selection_table, "re_selection_table.csv", DIR_ANOVAS)
cat("  wrote re_selection_table.csv\n")


################################################################################
###  ANALYSIS 3 — R1 COUNT MODELS via the DECISION TREE
################################################################################

cat("\n\n########## R1 COUNT MODELS (decision tree) ##########\n")

emm_side <- ~ group * treatment_clean | side
R1 <- list()

R1$R1a <- fit_aggression_count("R1a", "trial × segment × side",
  counts_trial_segment_side,
  "n_events ~ group * treatment_clean * side + segment + offset(offset_log)",
  list(trial_only = "(1 | trial)", trial_nested = "(1 | trial/segment)"), emm_side)

R1$R1b <- fit_aggression_count("R1b", "trial × side (aggregated)",
  counts_trial_side,
  "n_events ~ group * treatment_clean * side + offset(offset_log)",
  list(trial_only = "(1 | trial)"), emm_side)

R1$R1c_attacker <- fit_aggression_count("R1c_attacker", "fish × segment × side (attacker)",
  counts_fish_segment_side_attacker,
  "n_events ~ group * treatment_clean * side + segment + offset(offset_log)",
  list(trial_only = "(1 | trial)",
       trial_fish = "(1 | trial) + (1 | fish_attacker)"), emm_side)

R1$R1c_attacked <- fit_aggression_count("R1c_attacked", "fish × segment × side (attacked)",
  counts_fish_segment_side_attacked,
  "n_events ~ group * treatment_clean * side + segment + offset(offset_log)",
  list(trial_only = "(1 | trial)",
       trial_fish = "(1 | trial) + (1 | fish_attacked)"), emm_side)

R1$R1d_attacker <- fit_aggression_count("R1d_attacker", "cross-segment fish × side (attacker)",
  counts_fish_side_attacker,
  "n_events ~ group * treatment_clean * side + offset(offset_log)",
  list(trial_only = "(1 | trial)",
       trial_fish = "(1 | trial) + (1 | fish_attacker)",
       trial_nested = "(1 | trial/fish_attacker)"), emm_side)

R1$R1d_attacked <- fit_aggression_count("R1d_attacked", "cross-segment fish × side (attacked)",
  counts_fish_side_attacked,
  "n_events ~ group * treatment_clean * side + offset(offset_log)",
  list(trial_only = "(1 | trial)",
       trial_fish = "(1 | trial) + (1 | fish_attacked)",
       trial_nested = "(1 | trial/fish_attacked)"), emm_side)

# Leave-one-trial-out stability on the most interpretable fish-level count model.
loo_r1d <- if (!is.null(R1$R1d_attacker))
  tryCatch(loo_trial_stability(R1$R1d_attacker$model,
           as.data.frame(counts_fish_side_attacker)), error = function(e) NULL) else NULL
if (!is.null(loo_r1d)) write_df(loo_r1d, "loo_R1d_attacker.csv", DIR_ANOVAS)


################################################################################
###  ANALYSIS 3b — FULL FOUR-WAY INTERACTION TEST
###  group × treatment × side × interval(segment)
###
###  The primary R1 count models carry `segment` only as an ADDITIVE covariate,
###  so the group×treatment×side×interval interaction is never estimated there
###  (and the question "does the group×treatment×side structure differ between
###  intervals?" is left untested). This module fits the SATURATED fixed-effects
###  structure
###     n_events ~ group * treatment_clean * side * segment + offset(offset_log)
###  to obtain a formal Type-III test of every interaction term — in particular
###  the 4-way group:treatment_clean:side:segment. Family is selected by the same
###  Poisson→dispersion→NB rule as the decision tree; RE selection is reused. The
###  ANOVA table and a plain-language verdict are written for auditability.
################################################################################

# Fit the saturated group×treatment×side×segment count model + Type-III ANOVA.
#   label         short id for the output file, e.g. "R1a_4way"
#   level         human-readable scope
#   data          count grid carrying group, treatment_clean, side, segment
#   re_candidates named list of RE strings for AICc selection
# Returns list(model, family, anova, fourway_term, fourway_p, fourway_sig).
fit_fourway_interaction <- function(label, level, data, re_candidates) {
  cat("\n============================================================\n")
  cat("4-way interaction test:", label, "—", level, "\n")
  cat("============================================================\n")
  data <- as.data.frame(data)
  fixed_str <- "n_events ~ group * treatment_clean * side * segment + offset(offset_log)"

  # Family: start Poisson + AICc RE selection, then conditional dispersion → NB.
  sel <- aicc_re_select(fixed_str, re_candidates, data, fit_pois_model,
                        paste0(label, "_POIS"))
  if (is.null(sel)) { warning("4-way ", label, ": Poisson fit failed; skipping."); return(NULL) }
  best <- sel$model; fam <- "Poisson"
  od <- tryCatch(performance::check_overdispersion(best), error = function(e) NULL)
  od_p <- if (!is.null(od)) as.numeric(od$p_value) else NA_real_
  if (!is.na(od_p) && od_p < ALPHA) {
    sel_nb <- aicc_re_select(fixed_str, re_candidates, data, fit_count_model,
                             paste0(label, "_NB"))
    if (!is.null(sel_nb)) { best <- sel_nb$model; fam <- "NB" }
  }
  cat(sprintf("  family: %s (conditional dispersion p=%s)\n",
              fam, format(round(od_p, 3))))

  an <- tryCatch(car::Anova(best, type = 3), error = function(e) {
          message("  ANOVA failed: ", e$message); NULL })
  if (is.null(an)) return(NULL)
  an_df <- cbind(term = rownames(an), as.data.frame(an), row.names = NULL)
  write_df(an_df, sprintf("ANOVA_%s.csv", label), DIR_ANOVAS)

  fw_term <- "group:treatment_clean:side:segment"
  fw_p <- suppressWarnings(as.numeric(an_df[["Pr(>Chisq)"]][an_df$term == fw_term]))
  fw_p <- if (length(fw_p) == 1) fw_p else NA_real_
  fw_sig <- isTRUE(fw_p < ALPHA)
  cat(sprintf("  4-way %s: p=%s → %s\n", fw_term, format(round(fw_p, 4)),
              if (fw_sig) "SIGNIFICANT" else "not significant"))

  # Plain-language verdict line (appended to a shared report file).
  verdict <- sprintf(
    "%s (%s) | family=%s | 4-way %s p=%s → %s",
    label, level, fam, fw_term, format(round(fw_p, 4)),
    if (fw_sig) "interaction SUPPORTED" else "interaction NOT supported")
  safe_write_lines(verdict,
                   file.path(DIR_REPORTS, sprintf("fourway_interaction_%s.txt", label)))

  list(model = best, family = fam, anova = an_df,
       fourway_term = fw_term, fourway_p = fw_p, fourway_sig = fw_sig)
}

cat("\n\n########## 4-WAY INTERACTION TEST (group×treatment×side×interval) ##########\n")
R1_4way <- list()
R1_4way$R1a <- fit_fourway_interaction(
  "R1a_4way", "trial × segment × side (saturated)",
  counts_trial_segment_side,
  list(trial_only = "(1 | trial)", trial_nested = "(1 | trial/segment)"))


################################################################################
###  ANALYSIS 4 — R2 / R3 BINOMIAL MODELS  (logit; count tree N/A)
################################################################################

cat("\n\n########## R2 P(flow) + R3 crowding (binomial) ##########\n")

emm_gt <- ~ group * treatment_clean
R2 <- list(); R3 <- list()

R2$R2a <- fit_aggression_binom("R2a", "R2 (binomial)", "P(flow) trial × segment",
  probs_trial_segment, "cbind(n_flow, n_calm) ~ group * treatment_clean + segment + (1 | trial)",
  emm_gt)
R2$R2b <- fit_aggression_binom("R2b", "R2 (binomial)", "P(flow) trial (aggregated)",
  probs_trial, "cbind(n_flow, n_calm) ~ group * treatment_clean + (1 | trial)", emm_gt)
R2$R2c <- fit_aggression_binom("R2c", "R2 (binomial)", "P(flow) fish × segment",
  probs_fish_segment[n_total >= 2],
  "cbind(n_flow, n_total - n_flow) ~ group * treatment_clean + segment + (1 | trial) + (1 | fish_attacker)",
  emm_gt, notes = "cells with n_total >= 2")
R2$R2d <- fit_aggression_binom("R2d", "R2 (binomial)", "P(flow) cross-segment fish",
  probs_fish[n_total >= 2],
  "cbind(n_flow, n_total - n_flow) ~ group * treatment_clean + (1 | trial)",
  emm_gt, notes = "cells with n_total >= 2")

R3$R3a <- fit_aggression_binom("R3a", "R3 (binomial)", "crowding → side, event seg-resolved",
  crowding_events, "side_flow ~ crowding_diff_z + group * treatment_clean + segment + (1 | trial)",
  emm_gt)
R3$R3b <- fit_aggression_binom("R3b", "R3 (binomial)", "crowding → side, event aggregated",
  crowding_events, "side_flow ~ crowding_diff_z + group * treatment_clean + (1 | trial)",
  emm_gt)
R3$R3c <- fit_aggression_binom("R3c", "R3 (binomial)", "crowding → side, fish × segment",
  crowding_fish_segment[n_total >= 2],
  "cbind(n_flow_attacks, n_total - n_flow_attacks) ~ mean_crowding_diff_z + group * treatment_clean + segment + (1 | trial) + (1 | fish_attacker)",
  emm_gt, notes = "cells with n_total >= 2")
R3$R3d <- fit_aggression_binom("R3d", "R3 (binomial)", "crowding → side, cross-segment fish",
  crowding_fish[n_total >= 2],
  "cbind(n_flow_attacks, n_total - n_flow_attacks) ~ mean_crowding_diff_z + group * treatment_clean + (1 | trial)",
  emm_gt, notes = "cells with n_total >= 2")


################################################################################
###  ANALYSIS 5 — DESCRIPTIVE RAW-DATA FIGURES  (ported from 02_describe.R)
################################################################################

cat("\n\n########## DESCRIPTIVE FIGURES ##########\n")

save_fig <- function(p, name, w = 7, h = 5) {
  tryCatch(ggsave(file.path(DIR_FIGURES, name), p, width = w, height = h, dpi = 150),
           error = function(e) message("  fig failed ", name, ": ", e$message))
}
rates_int <- counts_trial_segment_side[, .(rate = sum(n_events) / 10),
                                       by = .(trial, segment, group, treatment_clean)]
save_fig(plot_rate_strip(rates_int, "Aggression rate — all intervals, sides combined (raw)"),
         "fig01_rate_interval_resolved.png")
rates_int_side <- counts_trial_segment_side[, .(trial, segment, side, group,
                                                 treatment_clean, rate = n_events / 10)]
save_fig(plot_rate_strip_side(rates_int_side, "Aggression rate — by side, all intervals (raw)"),
         "fig01a_rate_interval_by_side.png")
rates_agg <- counts_trial_side[, .(rate = sum(n_events) / 20),
                               by = .(trial, group, treatment_clean)]
save_fig(plot_rate_strip(rates_agg, "Aggression rate — intervals combined, sides combined (raw)"),
         "fig02_rate_aggregated.png")
rates_agg_side <- counts_trial_side[, .(trial, side, group, treatment_clean,
                                        rate = n_events / 20)]
save_fig(plot_rate_strip_side(rates_agg_side, "Aggression rate — by side, intervals combined (raw)"),
         "fig02b_rate_aggregated_by_side.png")
save_fig(plot_side_mosaic(events_full), "fig03_side_distribution.png", w = 8)
save_fig(plot_coverage_bars(events_full), "fig07_coverage_bars.png", w = 10)

# Forest panels from emmeans contrasts (when available).
build_forest <- function(emm_obj, model_label) {
  if (is.null(emm_obj) || is.null(emm_obj$contrasts)) return(NULL)
  d <- emm_obj$contrasts; d$model <- model_label; d
}
forest_rows <- Filter(Negate(is.null), list(
  if (!is.null(R1$R1d_attacker)) build_forest(R1$R1d_attacker$emm, "R1d: aggressive rate"),
  if (!is.null(R2$R2b))         build_forest(R2$R2b$emm, "R2b: P(flow)")))
if (length(forest_rows) > 0) {
  fdat <- rbindlist(forest_rows, fill = TRUE)
  if (all(c("estimate", "lower.CL", "upper.CL", "contrast") %in% names(fdat)))
    save_fig(plot_forest(fdat, xlab = "IRR / OR (response scale)",
             title = "Model contrasts (R1d, R2b)") + facet_wrap(~ model, scales = "free_x"),
             "fig05_forest.png", w = 10)
}


################################################################################
###  ANALYSIS 6 — SENSITIVITY S1–S5  (ported from 06_sensitivity.R)
################################################################################

cat("\n\n########## SENSITIVITY ANALYSES ##########\n")

extract_key_coefs <- function(model, label) {
  cf <- tryCatch(as.data.frame(summary(model)$coefficients$cond),
                 error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  cf$coef <- rownames(cf); cf$model <- label
  cf[grepl("group|treatment|side", cf$coef, ignore.case = TRUE), ]
}
sens <- list()

# S1 — events_strict subset (R1b structure)
ct_t_strict <- events_strict[, .(n_events = .N, group = first(group),
  treatment_clean = first(treatment_clean)), by = .(trial, side)]
ct_t_strict[, offset_log := log(1200)]
m_s1 <- tryCatch(fit_count_model(
  n_events ~ group * treatment_clean * side + offset(offset_log) + (1 | trial),
  data = ct_t_strict), error = function(e) NULL)
if (!is.null(m_s1)) sens$S1_R1b_strict <- extract_key_coefs(m_s1, "S1 R1b (events_strict)")
if (!is.null(R1$R1b)) sens$S1_R1b_primary <- extract_key_coefs(R1$R1b$model, "Primary R1b")

# S2 — exclude rejected_conflict
ct_t_s2 <- events_sensitivity_idmatcher[, .(n_events = .N, group = first(group),
  treatment_clean = first(treatment_clean)), by = .(trial, side)]
ct_t_s2[, offset_log := log(1200)]
m_s2 <- tryCatch(fit_count_model(
  n_events ~ group * treatment_clean * side + offset(offset_log) + (1 | trial),
  data = ct_t_s2), error = function(e) NULL)
if (!is.null(m_s2)) sens$S2_R1b <- extract_key_coefs(m_s2, "S2 R1b (no rejected_conflict)")

# S3 — add side_filter_applied covariate to R3a
m_s3 <- tryCatch(fit_binom_model(
  side_flow ~ crowding_diff_z + group * treatment_clean + segment +
    side_filter_applied + (1 | trial), data = crowding_events),
  error = function(e) NULL)
if (!is.null(m_s3)) {
  cf <- as.data.frame(summary(m_s3)$coefficients$cond); cf$coef <- rownames(cf)
  cf$model <- "S3 R3a + side_filter_applied"
  sens$S3 <- cf[grepl("side_filter|crowding", cf$coef), ]
}

# S4 — leave-one-cell-out (low-N cells)
s4_excl <- list(
  exclude_NoChoice_Exercise = list(group = "No choice", treatment_clean = "Exercise"),
  exclude_ExChoice_Control  = list(group = "Exercise choice", treatment_clean = "Control"))
for (nm in names(s4_excl)) {
  e <- s4_excl[[nm]]
  ct_loo <- counts_trial_side[!(group == e$group & treatment_clean == e$treatment_clean)]
  m <- tryCatch(fit_count_model(
    n_events ~ group * treatment_clean * side + offset(offset_log) + (1 | trial),
    data = ct_loo), error = function(e) NULL)
  if (!is.null(m)) sens[[paste0("S4_", nm)]] <- extract_key_coefs(m, nm)
}

# S5 — valid_fish_in_window has no variance (always 4): documented, skipped.
sens$S5_note <- data.frame(model = "S5", coef = "valid_fish_in_window",
  note = "always 4 for all events; no variance; predictor excluded",
  stringsAsFactors = FALSE)

sensitivity_table <- rbindlist(lapply(Filter(Negate(is.null), sens), as.data.table),
                               fill = TRUE)
write_df(sensitivity_table, "sensitivity_table.csv", DIR_ANOVAS)
cat("  wrote sensitivity_table.csv\n")


################################################################################
###  REPORTS — modelling-choices log (CSV + Word) + summary
################################################################################

cat("\n\n########## REPORTS ##########\n")

decision_df  <- rbindlist(DECISION_LOG, fill = TRUE)
decision_dfr <- as.data.frame(decision_df)   # plain frame for report helpers
write_df(decision_df, "modelling_choices.csv", DIR_REPORTS)
cat("  wrote reports/modelling_choices.csv (", nrow(decision_df), " responses)\n", sep = "")

# Shared, up-to-date description of the modelling workflow (used in every report).
TREE_DESC <- paste0(
  "Count responses (R1): start from a Poisson GLMM; test CONDITIONAL overdispersion on the ",
  "Poisson fit (performance::check_overdispersion) -> refit Negative Binomial if overdispersed; ",
  "test zero-inflation (DHARMa::testZeroInflation) -> refit Zero-Inflated NB if zero-inflated; ",
  "final fit checked with DHARMa::testUniformity (KS). Random-effect structure chosen by AICc ",
  "(singular fits rejected). All inference on the count / IRR (response) scale. ",
  "Binomial responses (R2, R3): logit GLMM; the count decision tree is not applicable.")


###  Module A — modelling choices as Markdown  ###############################

safe_write_lines(c(
  "# Aggression Analysis — Modelling Choices Log",
  "",
  sprintf("Generated %s from `%s`.", TIMESTAMP, basename(MASTER_CSV)),
  "",
  TREE_DESC,
  "",
  "## Decision log (one row per response)",
  "",
  df_to_md(decision_dfr),
  "",
  "## Random-effect selection (AICc)",
  "",
  df_to_md(as.data.frame(re_selection_table)),
  "",
  "## Distribution audit (marginal descriptors only)",
  "",
  df_to_md(as.data.frame(audit_table))
), file.path(DIR_REPORTS, "modelling_choices.md"))
cat("  wrote reports/modelling_choices.md\n")

# Word modelling-choices report (kept; description now matches the corrected tree).
tryCatch({
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Aggression analysis — modelling choices report",
                               style = "heading 1")
  doc <- officer::body_add_par(doc, sprintf("Generated %s. Input: %s.",
                                            TIMESTAMP, basename(MASTER_CSV)), style = "Normal")
  doc <- officer::body_add_par(doc, TREE_DESC, style = "Normal")
  doc <- officer::body_add_par(doc, "Decision log", style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(decision_dfr, font_size = 7L))
  doc <- officer::body_add_par(doc, "Random-effect selection (AICc)", style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(re_selection_table, font_size = 7L))
  doc <- officer::body_add_par(doc, "Distribution audit", style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(audit_table, font_size = 7L))
  print(doc, target = file.path(DIR_REPORTS, "modelling_choices_report.docx"))
  cat("  wrote reports/modelling_choices_report.docx\n")
}, error = function(e) message("  Word modelling-choices report failed: ", e$message))


###  Module B — full analysis report (Word + Markdown)  ######################

# Ordered inventory of fitted responses with human-readable scope + meaning.
report_items <- Filter(function(x) !is.null(x$res), list(
  list(label = "R1a",          fam = "R1", kind = "count", scope = "Trial × segment × side",          res = R1$R1a),
  list(label = "R1b",          fam = "R1", kind = "count", scope = "Trial × side (intervals aggregated)", res = R1$R1b),
  list(label = "R1c_attacker", fam = "R1", kind = "count", scope = "Fish × segment × side — attacker identity", res = R1$R1c_attacker),
  list(label = "R1c_attacked", fam = "R1", kind = "count", scope = "Fish × segment × side — attacked identity", res = R1$R1c_attacked),
  list(label = "R1d_attacker", fam = "R1", kind = "count", scope = "Cross-segment fish × side — attacker", res = R1$R1d_attacker),
  list(label = "R1d_attacked", fam = "R1", kind = "count", scope = "Cross-segment fish × side — attacked", res = R1$R1d_attacked),
  list(label = "R2a", fam = "R2", kind = "binom", scope = "P(flow): trial × segment",        res = R2$R2a),
  list(label = "R2b", fam = "R2", kind = "binom", scope = "P(flow): trial (aggregated)",      res = R2$R2b),
  list(label = "R2c", fam = "R2", kind = "binom", scope = "P(flow): fish × segment",          res = R2$R2c),
  list(label = "R2d", fam = "R2", kind = "binom", scope = "P(flow): cross-segment fish",      res = R2$R2d),
  list(label = "R3a", fam = "R3", kind = "binom", scope = "Crowding → side: event, segment-resolved", res = R3$R3a),
  list(label = "R3b", fam = "R3", kind = "binom", scope = "Crowding → side: event, aggregated",       res = R3$R3b),
  list(label = "R3c", fam = "R3", kind = "binom", scope = "Crowding → side: fish × segment",          res = R3$R3c),
  list(label = "R3d", fam = "R3", kind = "binom", scope = "Crowding → side: cross-segment fish",      res = R3$R3d)
))

fam_title <- function(f) switch(f,
  R1 = "Aggression rate (count models)",
  R2 = "Side preference P(flow) (binomial)",
  R3 = "Crowding → side (binomial)", f)

# Family-level overview table (the "R" introduction table).
r_overview <- data.frame(
  R          = c("R1", "R2", "R3"),
  Meaning    = c(
    "Aggression rate: number of high-intensity aggression events per unit time (exposure offset).",
    "Side preference P(flow): probability an aggression event occurs on the FLOW side vs the calm side.",
    "Crowding → side: whether the local crowding difference (individuals in flow minus calm) predicts aggression on the flow side."),
  Response   = c("event counts + offset", "binomial: flow vs calm", "binomial: flow (0/1) / proportion"),
  `Model family` = c("Poisson / NB / ZINB GLMM (decision tree)", "Binomial (logit) GLMM", "Binomial (logit) GLMM"),
  Predictors = c("group * treatment * side (+ segment)",
                 "group * treatment (+ segment)",
                 "crowding_diff_z + group * treatment (+ segment)"),
  check.names = FALSE, stringsAsFactors = FALSE)

# Decision-log row lookup (as a named list).
decrow <- function(lbl) {
  r <- decision_dfr[decision_dfr$response == lbl, , drop = FALSE]
  if (nrow(r) == 0) NULL else as.list(r[1, ])
}

# Response inventory table (one row per fitted model).
scope_tbl <- do.call(rbind, lapply(report_items, function(it) {
  d <- decrow(it$label)
  data.frame(Response = it$label, `R-family` = it$fam, Scope = it$scope,
             `Final model` = if (is.null(d)) NA else d$final_model,
             `N (cells/obs)` = if (is.null(d)) NA else d$n,
             check.names = FALSE, stringsAsFactors = FALSE)
}))

# Per-response "modelling steps" tables.
count_steps <- function(d) data.frame(
  Step = c("1. Marginal dispersion (descriptor only)",
           "2. Base model + random effects (AICc)",
           "3. Overdispersion test on Poisson fit",
           "4. Zero-inflation test",
           "5. Residual fit (DHARMa KS uniformity)",
           "Final model"),
  Detail = c(
    sprintf("mean=%s, var=%s, var/mean=%s", d$mean, d$var, d$var_over_mean),
    sprintf("Poisson; RE = %s%s", d$re_selected,
            if (isTRUE(as.logical(d$singular_dropped))) " (singular RE dropped to simplest)" else ""),
    sprintf("dispersion ratio=%s, p=%s → %s", d$dispersion_ratio, d$dispersion_p,
            if (!is.na(d$dispersion_p) && d$dispersion_p < ALPHA) "refit Negative Binomial" else "Poisson adequate"),
    sprintf("p=%s → %s", d$zero_infl_p,
            if (isTRUE(as.logical(d$zinb_triggered))) "refit Zero-Inflated NB" else "no zero-inflation"),
    sprintf("KS p=%s (%s)", d$dharma_uniformity_p,
            if (!is.na(d$dharma_uniformity_p) && d$dharma_uniformity_p >= ALPHA) "adequate" else "flagged"),
    d$final_model),
  stringsAsFactors = FALSE)

binom_steps <- function(d) data.frame(
  Step = c("Response type", "Model family", "Link",
           "Overdispersion (DHARMa)", "Residual fit (DHARMa KS uniformity)", "Final model"),
  Detail = c("binary / binomial proportion", "Binomial GLMM", "logit",
             sprintf("ratio=%s, p=%s", d$dispersion_ratio, d$dispersion_p),
             sprintf("KS p=%s", d$dharma_uniformity_p),
             d$final_model),
  stringsAsFactors = FALSE)

# Tidy ANOVA / contrasts / CLD for display.
clean_anova <- function(an) {
  if (is.null(an)) return(NULL)
  d <- as.data.frame(an); d <- cbind(Term = rownames(d), d); rownames(d) <- NULL
  pc <- grep("^Pr", names(d)); if (length(pc)) names(d)[pc] <- "p_value"
  d
}
clean_contrasts <- function(ct) {
  if (is.null(ct)) return(NULL)
  keep <- intersect(c("contrast", "side", "estimate", "SE", "lower.CL", "upper.CL",
                      "z.ratio", "p.value"), names(ct))
  ct[, keep, drop = FALSE]
}
clean_cld <- function(cld) {
  if (is.null(cld)) return(NULL)
  d <- as.data.frame(cld); names(d) <- sub("^\\.group$", "group_letters", names(d))
  d[] <- lapply(d, function(c) if (is.character(c)) trimws(c) else c)
  d
}

# Collect significant (p <= 0.05) effects. DECISION (2026-08-24): the
# "marginal" (0.05 < p <= 0.095) verbal significance tier is retired --
# RESULTS_STYLE.md R11 forbids a verbal tier beyond the exact p-value
# ("no trend towards significance / marginally significant"); every p is
# reported exactly (or at its floor) regardless of which side of 0.05 it
# falls on, and the reader draws their own conclusion from the number.
collect_sig <- function(items, thr_sig = 0.05) {
  out <- list()
  for (it in items) {
    an <- clean_anova(it$res$anova)
    if (!is.null(an) && "p_value" %in% names(an)) {
      stat_col <- intersect(c("Chisq", "LR Chisq", "F"), names(an))[1]
      an <- an[an$Term != "(Intercept)", , drop = FALSE]   # intercept is not a result
      sub <- an[!is.na(an$p_value) & an$p_value <= thr_sig, , drop = FALSE]
      for (i in seq_len(nrow(sub)))
        out[[length(out) + 1]] <- data.frame(
          response = it$label, family = it$fam, source = "ANOVA (Type III)",
          effect = sub$Term[i],
          statistic = if (!is.na(stat_col)) round(sub[[stat_col]][i], 3) else NA,
          p_value = sub$p_value[i],
          stringsAsFactors = FALSE)
    }
    ct <- it$res$emm$contrasts
    if (!is.null(ct) && "p.value" %in% names(ct)) {
      sub <- ct[!is.na(ct$p.value) & ct$p.value <= thr_sig, , drop = FALSE]
      for (i in seq_len(nrow(sub))) {
        lbl <- as.character(sub$contrast[i])
        if ("side" %in% names(sub)) lbl <- paste0(lbl, "  | side=", sub$side[i])
        out[[length(out) + 1]] <- data.frame(
          response = it$label, family = it$fam, source = "post-hoc contrast",
          effect = lbl,
          statistic = if ("estimate" %in% names(sub)) round(sub$estimate[i], 3) else NA,
          p_value = sub$p.value[i],
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) return(data.frame())
  res <- do.call(rbind, out)
  res[order(res$p_value), ]
}

sig_tbl <- collect_sig(report_items)

# ---- Markdown render (catered to AI comprehension) ----
ml <- c(
  "# Aggression Analysis — Statistical Report",
  "",
  sprintf("> Machine-readable report. Generated %s from `%s`. Each table is preceded by a one-line description of what it contains and how to read it. All tables are GitHub-flavored Markdown.",
          TIMESTAMP, basename(MASTER_CSV)),
  "",
  "## How to read this report",
  "",
  "- **Design / reference levels:** `group` ∈ {No choice (ref), Exercise choice}; `treatment` ∈ {Control (ref), Exercise}; `side` ∈ {calm (ref), flow}; `segment` ∈ {2, 3}. Effects/contrasts are expressed relative to these reference levels.",
  "- **Count models (R1):** estimated marginal contrasts are **Incidence Rate Ratios (IRR)**. IRR < 1 ⇒ the first level of the contrast has a *lower* event rate; IRR > 1 ⇒ *higher*. Confidence intervals that exclude 1 indicate significance.",
  "- **Binomial models (R2, R3):** contrasts are **Odds Ratios (OR)**; OR is interpreted like IRR but for odds. CIs excluding 1 indicate significance.",
  "- **p-value convention used in the summary:** `significant` = p ≤ 0.05, reported to the p < 0.001 floor (no verbal tier beyond the exact p-value -- handoff/RESULTS_STYLE.md R11).",
  "- **Modelling-steps tables** record how each count response's family was chosen by the decision tree (Poisson → NB → ZINB) and the final residual diagnostic.",
  "- **Abbreviations:** GLMM = generalized linear mixed model; NB = negative binomial; ZINB = zero-inflated NB; RE = random effect; AICc = corrected Akaike information criterion; ANOVA = Type III Wald χ² test of fixed effects; CLD = compact letter display (groups sharing a letter are not significantly different); KS = Kolmogorov–Smirnov uniformity test of DHARMa residuals.",
  "",
  TREE_DESC,
  "",
  "## 1. Response overview",
  "",
  "Each row is one analysis family (R1–R3) with its biological question, response type, and model family.",
  "",
  df_to_md(r_overview),
  "",
  "### Response inventory",
  "",
  "Every fitted model, its scope, the family chosen by the decision tree, and its sample size.",
  "",
  df_to_md(scope_tbl),
  "",
  "## 2. Per-response results",
  ""
)

cur_fam <- ""
for (it in report_items) {
  if (it$fam != cur_fam) {
    cur_fam <- it$fam
    ml <- c(ml, sprintf("### %s — %s", cur_fam, fam_title(cur_fam)), "")
  }
  d <- decrow(it$label)
  steps <- if (it$kind == "count") count_steps(d) else binom_steps(d)
  an <- clean_anova(it$res$anova)
  ct <- clean_contrasts(it$res$emm$contrasts)
  cl <- clean_cld(it$res$cld)
  interp <- if (it$kind == "count")
    "Contrasts below are **IRR** (rate ratios); IRR<1 = first level lower." else
    "Contrasts below are **OR** (odds ratios); OR<1 = first level lower odds."
  ml <- c(ml,
    sprintf("#### %s — %s", it$label, it$scope), "",
    sprintf("Final model: **%s**. %s", if (is.null(d)) "?" else d$final_model, interp), "",
    "**Modelling steps** — how the model/family was decided.", "", df_to_md(steps), "",
    "**Type III ANOVA** — Wald χ² test of each fixed effect (and interactions).", "",
    if (is.null(an)) "_(ANOVA unavailable)_\n" else df_to_md(an), "",
    "**Post-hoc pairwise contrasts** (Sidak-adjusted) — group×treatment comparisons, by side where applicable.", "",
    if (is.null(ct)) "_(no contrasts)_\n" else df_to_md(ct), "",
    "**Compact letter display (CLD)** — marginal means; levels sharing a letter do not differ.", "",
    if (is.null(cl)) "_(no CLD)_\n" else df_to_md(cl), "")
}

ml <- c(ml,
  "## 3. Significant results (p ≤ 0.05)",
  "",
  "Pooled across all responses: every ANOVA term and post-hoc contrast with p ≤ 0.05. `statistic` = χ² for ANOVA terms, IRR/OR estimate for contrasts.",
  "",
  if (nrow(sig_tbl) == 0) "_No effects reached p ≤ 0.05._\n" else df_to_md(sig_tbl))

safe_write_lines(ml, file.path(DIR_REPORTS, "analysis_report.md"))
cat("  wrote reports/analysis_report.md\n")

# ---- Word render (mirrors the Markdown) ----
tryCatch({
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Aggression Analysis — Statistical Report", style = "heading 1")
  doc <- officer::body_add_par(doc, sprintf("Generated %s. Input: %s.",
                                            TIMESTAMP, basename(MASTER_CSV)), style = "Normal")
  doc <- officer::body_add_par(doc, TREE_DESC, style = "Normal")
  doc <- officer::body_add_par(doc,
    "Reference levels: group=No choice, treatment=Control, side=calm. Count contrasts are IRR; binomial contrasts are OR (values <1 = first level lower). CLD: levels sharing a letter do not differ.",
    style = "Normal")

  doc <- officer::body_add_par(doc, "1. Response overview", style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(r_overview, 8L))
  doc <- officer::body_add_par(doc, "Response inventory", style = "heading 3")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(scope_tbl, 8L))

  doc <- officer::body_add_par(doc, "2. Per-response results", style = "heading 2")
  cur_fam <- ""
  for (it in report_items) {
    if (it$fam != cur_fam) {
      cur_fam <- it$fam
      doc <- officer::body_add_par(doc, sprintf("%s — %s", cur_fam, fam_title(cur_fam)),
                                   style = "heading 3")
    }
    d <- decrow(it$label)
    steps <- if (it$kind == "count") count_steps(d) else binom_steps(d)
    an <- clean_anova(it$res$anova)
    ct <- clean_contrasts(it$res$emm$contrasts)
    cl <- clean_cld(it$res$cld)
    doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(
      sprintf("%s — %s", it$label, it$scope),
      officer::fp_text(bold = TRUE, font.size = 11, font.family = "Calibri"))))
    doc <- officer::body_add_par(doc, "Modelling steps", style = "Normal")
    doc <- flextable::body_add_flextable(doc, .df_to_flex(steps, 8L))
    if (!is.null(an)) {
      doc <- officer::body_add_par(doc, "Type III ANOVA", style = "Normal")
      doc <- flextable::body_add_flextable(doc, .df_to_flex(an, 8L))
    }
    if (!is.null(ct)) {
      doc <- officer::body_add_par(doc, "Post-hoc pairwise contrasts (Sidak)", style = "Normal")
      doc <- flextable::body_add_flextable(doc, .df_to_flex(ct, 7L))
    }
    if (!is.null(cl)) {
      doc <- officer::body_add_par(doc, "Compact letter display", style = "Normal")
      doc <- flextable::body_add_flextable(doc, .df_to_flex(cl, 8L))
    }
  }

  doc <- officer::body_add_par(doc, "3. Significant results (p <= 0.05)",
                               style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(
    if (nrow(sig_tbl) == 0) data.frame(note = "No effects reached p <= 0.05") else sig_tbl, 8L))

  print(doc, target = file.path(DIR_REPORTS, "analysis_report.docx"))
  cat("  wrote reports/analysis_report.docx\n")
}, error = function(e) message("  Word analysis report failed: ", e$message))

# sessionInfo for reproducibility.
safe_write_lines(capture.output(sessionInfo()),
                 file.path(DIR_REPORTS, "sessionInfo.txt"))

cat("\n=== DONE ===\n")
cat("All outputs under:", OUTPUT_DIR, "\n")

# design_balance.csv (CONVENTIONS.md C12): trial is already the forced design
# RE in every RE candidate set in this script (R1a/R1b/R1c/R1d all carry
# (1|trial)); this documents what it is aliased with.
db <- tryCatch(
  sensitivity_design_balance(events_full, c("trial", "group", "treatment_clean", "segment", "side"),
                             label = "aggression_events", run_dir = DIR_ANOVAS),
  error = function(e) NULL)
if (!is.null(db)) cat("  wrote design_balance.csv\n")

# HANDOFF.md S10 [INV]: completed-run sentinel for downstream auto-discovery.
writeLines(as.character(Sys.time()), file.path(OUTPUT_DIR, "_RUN_COMPLETE"))

