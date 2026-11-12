# =============================================================================
# transformation_diagnostics.R
# Pre-modelling assumption diagnostics for the GOT pipeline
# =============================================================================
# Sourced by 00_run_trajectory_pipeline.R (before STEP5 and STEP6).
# Standalone guard at bottom allows this file to be sourced directly too.
#
# Functions provided:
#   build_diag_table_row()     — assemble the 1-row assumption-diagnostics table
#   run_diag_for_indicator()   — run SW + Levene on indicator, pick transformation
#   diag_row_beta()            — pre-filled row for Beta GLMM indicators
#   diag_row_count()           — pre-filled row for count GLMM indicators
#   run_agg_diag()             — post-hoc residual diagnostics for aggregated lm
# =============================================================================

# ---- Dependencies (all must be sourced/loaded before this file) ----
# run_shapiro_safe()        — stats_helpers.R
# run_levene_on_residuals() — stats_helpers.R
# These are accessed via match.fun() with a fallback inline guard.

# Ensure lme4 is available (used for diagnostic intercept models)
if (!requireNamespace("lme4", quietly = TRUE))
  stop("Package 'lme4' is required by transformation_diagnostics.R", call. = FALSE)


# =============================================================================
# Internal helpers
# =============================================================================

#' Format a p-value for display in the diagnostics table.
#' Returns "NA" for missing/non-finite; "< 0.0001" when p < 0.0001; else 3 dp.
.fmt_diag_p <- function(p) {
  if (is.null(p) || length(p) != 1L || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.0001) return("< 0.0001")
  formatC(p, format = "f", digits = 3L)
}

#' Run Shapiro-Wilk if run_shapiro_safe() is available, else inline fallback.
.safe_sw <- function(x) {
  if (exists("run_shapiro_safe", mode = "function"))
    return(run_shapiro_safe(x))
  x <- x[is.finite(x)]
  if (length(x) < 3L || length(x) > 5000L)
    return(list(statistic = NA_real_, p.value = NA_real_,
                note = "SW not run: n outside [3, 5000]"))
  res <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
  if (is.null(res)) return(list(statistic = NA_real_, p.value = NA_real_, note = "SW error"))
  list(statistic = unname(res$statistic), p.value = res$p.value, note = "")
}

#' Run Levene's test if run_levene_on_residuals() is available, else fallback.
.safe_levene <- function(resid_vec, grouping_fac) {
  resid_df <- data.frame(resid       = resid_vec,
                          group_combo = as.factor(as.character(grouping_fac)))
  if (exists("run_levene_on_residuals", mode = "function"))
    return(run_levene_on_residuals(resid_df))
  # Fallback: car::leveneTest
  out <- tryCatch(car::leveneTest(resid ~ group_combo, data = resid_df),
                  error = function(e) e)
  if (inherits(out, "error"))
    return(list(statistic = NA_real_, df1 = NA_real_, df2 = NA_real_,
                p.value = NA_real_, note = out$message))
  od  <- as.data.frame(out)
  fco <- grep("F value|F", names(od), value = TRUE)[1L]
  pco <- grep("^Pr\\(", names(od), value = TRUE)[1L]
  list(statistic = if (!is.na(fco)) od[1L, fco] else NA_real_,
       df1 = if ("Df" %in% names(od)) od[1L, "Df"] else NA_real_,
       df2 = if ("Df" %in% names(od) && nrow(od) >= 2L) od[2L, "Df"] else NA_real_,
       p.value  = if (!is.na(pco)) od[1L, pco] else NA_real_,
       note = "")
}

#' Fit a minimal lmer(y ~ RE) for residual extraction.
.fit_diag_lmer <- function(y_vec, data, re_form_str) {
  data$.diag_y_td <- y_vec
  # Fallback to intercept-only if re_form parse fails or lmer fails
  fml <- tryCatch(
    stats::as.formula(paste(".diag_y_td ~", re_form_str)),
    error = function(e) stats::as.formula(".diag_y_td ~ 1")
  )
  tryCatch(lme4::lmer(fml, data = data, REML = TRUE),
           error = function(e) NULL)
}

#' Evaluate one Gaussian candidate: transform y, fit intercept lmer, run SW + Levene.
#' Returns list(sw_p, lev_p, aicc) or NULL if transform is invalid.
.eval_gaussian_cand <- function(y_raw, y_tf, data, grouping_fac, re_form_str) {
  if (!all(is.finite(y_tf))) return(NULL)
  fit   <- .fit_diag_lmer(y_tf, data, re_form_str)
  resid <- if (!is.null(fit)) residuals(fit) else (y_tf - mean(y_tf, na.rm = TRUE))
  sw    <- .safe_sw(resid)
  lev   <- .safe_levene(resid, grouping_fac)
  aicc_val <- if (!is.null(fit))
    tryCatch(MuMIn::AICc(fit), error = function(e) NA_real_)
  else NA_real_
  list(sw_p = sw$p.value, lev_p = lev$p.value, aicc = aicc_val)
}


# =============================================================================
# build_diag_table_row()
# =============================================================================
#' Assemble a single-row data frame for the "Assumption diagnostics" table.
#'
#' @param analysis_unit   "session", "fish", or "trial (aggregated)"
#' @param family          e.g. "Gaussian LMM", "Beta GLMM", "Gamma LMM", "Gaussian LM"
#' @param n_obs           integer — number of observations in the diagnostic model
#' @param raw_shapiro_p   already-formatted string
#' @param raw_levene_p    already-formatted string
#' @param transform_selected  "identity", "log", "sqrt", "log1p", "Gamma LMM (log link)", etc.
#' @param transform_rationale  short explanation (≤ 20 words)
#' @param tf_shapiro_p    post-transform Shapiro p (formatted), or "not applicable"
#' @param tf_levene_p     post-transform Levene p (formatted), or "not applicable"
#' @param model_response  column name passed to the fitted model
#' @param flags           "" or comma-separated flag names
build_diag_table_row <- function(
    analysis_unit,
    family,
    n_obs,
    raw_shapiro_p,
    raw_levene_p,
    transform_selected,
    transform_rationale,
    tf_shapiro_p,
    tf_levene_p,
    model_response,
    flags = ""
) {
  data.frame(
    "Unit"             = as.character(analysis_unit),
    "Family"           = as.character(family),
    "N"                = as.integer(n_obs),
    "Shapiro p (raw)"  = as.character(raw_shapiro_p),
    "Levene p (raw)"   = as.character(raw_levene_p),
    "Transform"        = as.character(transform_selected),
    "Rationale"        = as.character(transform_rationale),
    "Shapiro p (tf)"   = as.character(tf_shapiro_p),
    "Levene p (tf)"    = as.character(tf_levene_p),
    "Model response"   = as.character(model_response),
    "Flags"            = as.character(flags),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# run_diag_for_indicator()
# =============================================================================
#' Run assumption diagnostics for a continuous (Gaussian) indicator in STEP6.
#'
#' Algorithm:
#'   1. Fit intercept + RE diagnostic model; extract residuals.
#'   2. Run Shapiro-Wilk (SW) and Levene on residuals.
#'   3. If both pass (p >= P_CUTOFF): identity selected.
#'   4. Otherwise build candidate set: identity, sqrt (if y >= 0), log (if y > 0),
#'      log1p (if y >= 0 and has zeros).
#'   5. Evaluate each Gaussian candidate with SW + Levene.
#'   6. Among candidates that pass both: prefer identity; else pick by max(min(SW_p, Lev_p)).
#'   7. If try_gamma = TRUE and all_positive and no Gaussian candidate passed:
#'      fit glmmTMB Gamma intercept model; run DHARMa KS test.
#'      If KS passes: set use_gamma = TRUE.
#'   8. If still nothing: use identity + flag "transform_did_not_improve".
#'
#' @param y             numeric vector (finite values only; same length as nrow(data))
#' @param data          data frame containing y and RE columns
#' @param grouping      factor for Levene (e.g. interaction(data$Group, data$Treatment))
#' @param re_form_str   RE formula string (e.g. "(1 | trial_id)") for diagnostic model
#' @param ind_nm        column name of y in data (for response_col naming)
#' @param analysis_unit "session" or "fish" — shown in diagnostics table
#' @param try_gamma     logical — attempt Gamma LMM as fallback when all Gaussian fail
#' @param P_CUTOFF      p-value threshold for assumption tests (default 0.05)
#'
#' @return list with:
#'   $row          — 1-row data.frame (the diagnostics table row)
#'   $transform    — character: "identity", "log", "sqrt", "log1p", or "gamma"
#'   $response_col — column name to pass to the main fitting function
#'   $use_gamma    — logical: TRUE when Gamma LMM should be used
run_diag_for_indicator <- function(y,
                                   data,
                                   grouping,
                                   re_form_str    = "(1 | trial_id)",
                                   ind_nm         = "outcome",
                                   analysis_unit  = "session",
                                   try_gamma      = TRUE,
                                   P_CUTOFF       = 0.05) {

  y <- y[is.finite(y)]
  n <- length(y)

  # ---- Minimum-data guard ---------------------------------------------------
  if (n < 3L) {
    row <- build_diag_table_row(
      analysis_unit, "Gaussian LMM", n,
      "NA", "NA",
      "identity", "Insufficient data for diagnostics",
      "NA", "NA", ind_nm, "small_n"
    )
    return(list(row = row, transform = "identity",
                response_col = ind_nm, use_gamma = FALSE))
  }

  # ---- Data-range flags -----------------------------------------------------
  flags_vec  <- character(0L)
  if (n < 20L) flags_vec <- c(flags_vec, "small_n")
  all_pos    <- all(y > 0,  na.rm = TRUE)
  all_nonneg <- all(y >= 0, na.rm = TRUE)
  pct_zero   <- mean(y == 0, na.rm = TRUE)
  if (pct_zero > 0.2) flags_vec <- c(flags_vec, "zero_inflated")

  # ---- Evaluate raw (identity) candidate -----------------------------------
  eval_raw  <- .eval_gaussian_cand(y, y, data, grouping, re_form_str)
  raw_sw_p  <- eval_raw$sw_p  %||% NA_real_
  raw_lev_p <- eval_raw$lev_p %||% NA_real_

  passes_raw <- isTRUE(!is.na(raw_sw_p)  && raw_sw_p  >= P_CUTOFF &&
                         !is.na(raw_lev_p) && raw_lev_p >= P_CUTOFF)

  if (passes_raw) {
    row <- build_diag_table_row(
      analysis_unit, "Gaussian LMM", n,
      .fmt_diag_p(raw_sw_p), .fmt_diag_p(raw_lev_p),
      "identity", "Raw residuals normal and homoscedastic",
      .fmt_diag_p(raw_sw_p), .fmt_diag_p(raw_lev_p),
      ind_nm, paste(flags_vec, collapse = ", ")
    )
    message("[DIAG] ", ind_nm, ": identity (SW p=", .fmt_diag_p(raw_sw_p),
            ", Lev p=", .fmt_diag_p(raw_lev_p), ")")
    return(list(row = row, transform = "identity",
                response_col = ind_nm, use_gamma = FALSE))
  }

  # ---- Build transformation candidate set ----------------------------------
  cands        <- list()   # named list: label → list(y_tf, sw_p, lev_p, aicc)
  eval_id      <- list(y_tf = y, sw_p = raw_sw_p, lev_p = raw_lev_p, aicc = eval_raw$aicc)
  cands[["identity"]] <- eval_id

  if (all_nonneg && !all(y == 0)) {
    y_sqrt    <- sqrt(y)
    ev_sqrt   <- .eval_gaussian_cand(y, y_sqrt, data, grouping, re_form_str)
    if (!is.null(ev_sqrt))
      cands[["sqrt"]] <- list(y_tf  = y_sqrt,
                               sw_p  = ev_sqrt$sw_p, lev_p = ev_sqrt$lev_p,
                               aicc  = ev_sqrt$aicc)
  }

  if (all_pos) {
    y_log   <- log(y)
    ev_log  <- .eval_gaussian_cand(y, y_log, data, grouping, re_form_str)
    if (!is.null(ev_log))
      cands[["log"]] <- list(y_tf = y_log,
                              sw_p = ev_log$sw_p, lev_p = ev_log$lev_p,
                              aicc = ev_log$aicc)
  } else if (all_nonneg && pct_zero > 0) {
    y_lp    <- log1p(y)
    ev_lp   <- .eval_gaussian_cand(y, y_lp, data, grouping, re_form_str)
    if (!is.null(ev_lp))
      cands[["log1p"]] <- list(y_tf = y_lp,
                                sw_p = ev_lp$sw_p, lev_p = ev_lp$lev_p,
                                aicc = ev_lp$aicc)
  }

  # ---- Summary table for Gaussian candidates --------------------------------
  cand_df <- do.call(rbind, lapply(names(cands), function(nm) {
    cc <- cands[[nm]]
    sw  <- if (is.finite(cc$sw_p))  cc$sw_p  else 0
    lev <- if (is.finite(cc$lev_p)) cc$lev_p else 0
    data.frame(label  = nm,
               sw_p   = sw,
               lev_p  = lev,
               aicc   = if (is.finite(cc$aicc)) cc$aicc else Inf,
               passes = (sw >= P_CUTOFF) & (lev >= P_CUTOFF),
               joint  = pmin(sw, lev),
               stringsAsFactors = FALSE)
  }))

  passing_g <- cand_df[cand_df$passes %in% TRUE, , drop = FALSE]

  # ---- Optional Gamma LMM diagnostic --------------------------------------
  gamma_ok   <- FALSE
  gamma_ks_p <- NA_real_
  if (try_gamma && all_pos && nrow(passing_g) == 0L &&
      requireNamespace("glmmTMB", quietly = TRUE) &&
      requireNamespace("DHARMa",  quietly = TRUE)) {
    data$.diag_y_gamma <- y
    fml_g <- tryCatch(
      stats::as.formula(paste(".diag_y_gamma ~", re_form_str)),
      error = function(e) stats::as.formula(".diag_y_gamma ~ 1")
    )
    fit_g <- tryCatch(
      glmmTMB::glmmTMB(
        fml_g, data = data,
        family  = Gamma(link = "log"),
        control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500L))
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_g)) {
      sim_g <- tryCatch(
        DHARMa::simulateResiduals(fit_g, n = 200L, plot = FALSE),
        error = function(e) NULL
      )
      if (!is.null(sim_g)) {
        gamma_ks_p <- tryCatch(
          DHARMa::testUniformity(sim_g, plot = FALSE)$p.value,
          error = function(e) NA_real_
        )
        gamma_ok <- isTRUE(!is.na(gamma_ks_p) && gamma_ks_p >= P_CUTOFF)
      }
    }
  }

  # ---- Decision: prefer Gamma when no Gaussian passes -----------------------
  if (gamma_ok && nrow(passing_g) == 0L) {
    flags_str <- paste(flags_vec, collapse = ", ")
    row <- build_diag_table_row(
      analysis_unit, "Gamma LMM (log link)", n,
      .fmt_diag_p(raw_sw_p), .fmt_diag_p(raw_lev_p),
      "Gamma LMM (log link)",
      paste0("Gaussian SW/Levene failed; Gamma DHARMa KS p = ", .fmt_diag_p(gamma_ks_p)),
      "not applicable", "not applicable",
      ind_nm, flags_str
    )
    message("[DIAG] ", ind_nm, ": Gamma LMM (KS p=", .fmt_diag_p(gamma_ks_p), ")")
    return(list(row = row, transform = "gamma",
                response_col = ind_nm, use_gamma = TRUE))
  }

  # ---- Decision: prefer identity if it passes ------------------------------
  if (any(passing_g$label == "identity")) {
    chosen_nm <- "identity"
  } else if (nrow(passing_g) > 0L) {
    passing_g <- passing_g[order(-passing_g$joint, passing_g$aicc), , drop = FALSE]
    chosen_nm <- passing_g$label[1L]
  } else {
    # None pass — pick best by joint score then AICc
    flags_vec <- c(flags_vec, "transform_did_not_improve")
    cand_df   <- cand_df[order(-cand_df$joint, cand_df$aicc), , drop = FALSE]
    chosen_nm <- cand_df$label[1L]
  }

  chosen_c   <- cands[[chosen_nm]]
  tf_sw_p    <- chosen_c$sw_p
  tf_lev_p   <- chosen_c$lev_p

  # Build response column name for non-identity transforms
  resp_col <- if (chosen_nm == "identity") {
    ind_nm
  } else {
    paste0(chosen_nm, "_", ind_nm)
  }

  rationale <- switch(chosen_nm,
    identity = "Raw residuals normal and homoscedastic",
    log      = paste0("log improved SW p: ", .fmt_diag_p(raw_sw_p),
                      " → ", .fmt_diag_p(tf_sw_p)),
    log1p    = paste0("log1p improved SW p: ", .fmt_diag_p(raw_sw_p),
                      " → ", .fmt_diag_p(tf_sw_p)),
    sqrt     = paste0("sqrt improved SW p: ", .fmt_diag_p(raw_sw_p),
                      " → ", .fmt_diag_p(tf_sw_p)),
    paste0(chosen_nm, " selected (best available; SW p raw=", .fmt_diag_p(raw_sw_p), ")")
  )

  flags_str <- paste(flags_vec, collapse = ", ")
  row <- build_diag_table_row(
    analysis_unit, "Gaussian LMM", n,
    .fmt_diag_p(raw_sw_p), .fmt_diag_p(raw_lev_p),
    chosen_nm, rationale,
    .fmt_diag_p(tf_sw_p), .fmt_diag_p(tf_lev_p),
    resp_col, flags_str
  )

  message("[DIAG] ", ind_nm, ": ", chosen_nm,
          " (SW p=", .fmt_diag_p(tf_sw_p),
          ", Lev p=", .fmt_diag_p(tf_lev_p), ")")
  list(row = row, transform = chosen_nm, response_col = resp_col, use_gamma = FALSE)
}


# =============================================================================
# Convenience rows for non-Gaussian indicators
# =============================================================================

#' Pre-filled diagnostics row for Beta GLMM indicators (bounded [0, 1]).
#' @param n_obs         integer
#' @param analysis_unit "session" or "fish"
#' @param ind_nm        raw column name
#' @param clipped       logical — TRUE if any values were exactly 0 or 1 before clipping
diag_row_beta <- function(n_obs, analysis_unit, ind_nm, clipped = FALSE) {
  build_diag_table_row(
    analysis_unit, "Beta GLMM (logit link)", n_obs,
    "not applicable", "not applicable",
    "Beta GLMM (logit link)", "Bounded [0, 1] outcome — Beta family used",
    "not applicable", "not applicable",
    ind_nm,
    if (clipped) "boundary_clip_applied" else ""
  )
}

#' Pre-filled diagnostics row for count GLMM indicators.
#' @param n_obs         integer
#' @param analysis_unit "session" or "fish"
#' @param ind_nm        raw column name
#' @param overdisp      overdispersion ratio (numeric, optional)
diag_row_count <- function(n_obs, analysis_unit, ind_nm, overdisp = NA_real_) {
  flag <- if (!is.na(overdisp) && is.finite(overdisp) && overdisp > 3) "severe_overdispersion" else ""
  build_diag_table_row(
    analysis_unit, "Count GLMM (Poisson / NB)", n_obs,
    "not applicable", "not applicable",
    "Poisson / NB GLMM", "Count outcome — Poisson / NB family used",
    "not applicable", "not applicable",
    ind_nm, flag
  )
}


# =============================================================================
# run_agg_diag()
# =============================================================================
#' Post-hoc residual diagnostics for trial-aggregated lm models.
#'
#' Fits no new models — tests residuals of the already-fitted lm.
#' Results are informational only; no transformation is applied retroactively.
#'
#' @param fit       fitted lm object (or NULL)
#' @param grouping  factor vector aligned with model data rows (Group × Treatment)
#' @param ind_nm    indicator name (for model_response column)
#' @param n_obs     integer — number of trial-level observations
#' @return 1-row data.frame from build_diag_table_row()
run_agg_diag <- function(fit, grouping, ind_nm, n_obs = NA_integer_) {

  n_obs_int <- as.integer(n_obs)

  if (is.null(fit)) {
    return(build_diag_table_row(
      "trial (aggregated)", "Gaussian LM", n_obs_int,
      "NA", "NA", "identity", "Model unavailable",
      "NA", "NA", ind_nm, "manual_review"
    ))
  }

  resid_vec <- tryCatch(residuals(fit), error = function(e) NULL)
  if (is.null(resid_vec) || length(resid_vec) == 0L) {
    return(build_diag_table_row(
      "trial (aggregated)", "Gaussian LM", n_obs_int,
      "NA", "NA", "identity", "Residuals unavailable",
      "NA", "NA", ind_nm, "manual_review"
    ))
  }

  sw  <- .safe_sw(resid_vec)
  # For Levene, align grouping to model rows via names if possible
  grp_agg <- tryCatch({
    if (!is.null(names(resid_vec)) && !is.null(names(grouping))) {
      grouping[match(names(resid_vec), names(grouping))]
    } else grouping[seq_along(resid_vec)]
  }, error = function(e) grouping[seq_along(resid_vec)])
  lev <- .safe_levene(resid_vec, grp_agg)

  sw_p  <- sw$p.value
  lev_p <- lev$p.value

  flags_vec <- character(0L)
  if (length(resid_vec) < 20L) flags_vec <- c(flags_vec, "small_n")
  if (isTRUE(!is.na(sw_p)  && sw_p  < 0.05) ||
      isTRUE(!is.na(lev_p) && lev_p < 0.05))
    flags_vec <- c(flags_vec, "manual_review")

  build_diag_table_row(
    "trial (aggregated)", "Gaussian LM",
    length(resid_vec),
    .fmt_diag_p(sw_p), .fmt_diag_p(lev_p),
    "identity",
    "Post-hoc check only — no transformation applied to aggregated model",
    .fmt_diag_p(sw_p), .fmt_diag_p(lev_p),
    ind_nm,
    paste(flags_vec, collapse = ", ")
  )
}


# =============================================================================
# diag_row_from_choose_tf()
# =============================================================================
#' Build a build_diag_table_row() output from a choose_transformation() result.
#'
#' Used by STEP5's analyze_continuous_indicator() to extract the assumption-
#' diagnostics table row from the already-computed choose_transformation() result,
#' avoiding re-computation.
#'
#' @param tf           output of choose_transformation() — or NULL
#' @param indicator    raw column name (e.g. "distance_m")
#' @param n_obs        number of fish-level observations used in modelling
#' @param analysis_unit "fish" or "session"
#' @return 1-row data.frame from build_diag_table_row(), or NULL if tf is NULL
diag_row_from_choose_tf <- function(tf, indicator, n_obs,
                                     analysis_unit = "fish") {
  if (is.null(tf)) return(NULL)

  chosen    <- tf$chosen_transform %||% "identity"
  raw_sw_p  <- tryCatch(tf$raw_eval$shapiro$p.value,   error = function(e) NA_real_)
  raw_lev_p <- tryCatch(tf$raw_eval$levene$p.value,    error = function(e) NA_real_)

  # Non-Gaussian fallbacks: these are handled by their own return paths
  if (chosen %in% c("gamma_glmm_shifted", "gamma_diagnostic", "beta_glmm")) {
    fam <- if (grepl("gamma", chosen, ignore.case = TRUE))
      "Gamma GLMM (log link)" else "Beta GLMM (logit link)"
    return(build_diag_table_row(
      analysis_unit, fam, as.integer(n_obs),
      "not applicable", "not applicable",
      fam, "Gaussian LMMs failed — GLMM family selected",
      "not applicable", "not applicable",
      indicator, ""
    ))
  }

  # Gaussian path: read post-transform diagnostics from transformed_eval
  tf_eval   <- tf$transformed_eval %||% tf$raw_eval
  tf_sw_p   <- tryCatch(tf_eval$shapiro$p.value, error = function(e) NA_real_)
  tf_lev_p  <- tryCatch(tf_eval$levene$p.value,  error = function(e) NA_real_)

  resp_col  <- if (identical(chosen, "identity")) indicator
               else paste0(chosen, "_", indicator)

  # Truncate decision note to <= 15 words
  note_raw   <- tf$decision_note %||% ""
  note_words <- strsplit(trimws(note_raw), "\\s+")[[1L]]
  rationale  <- paste(note_words[seq_len(min(15L, length(note_words)))], collapse = " ")

  flags_str <- if (!is.na(n_obs) && as.integer(n_obs) < 20L) "small_n" else ""

  build_diag_table_row(
    analysis_unit, "Gaussian LMM", as.integer(n_obs),
    .fmt_diag_p(raw_sw_p), .fmt_diag_p(raw_lev_p),
    chosen, rationale,
    .fmt_diag_p(tf_sw_p), .fmt_diag_p(tf_lev_p),
    resp_col, flags_str
  )
}


# =============================================================================
# Standalone guard — source dependencies if run outside the master pipeline
# =============================================================================
if (!exists("run_shapiro_safe", mode = "function") ||
    !exists("run_levene_on_residuals", mode = "function")) {
  .td_dir <- tryCatch({
    frames    <- sys.frames()
    ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
    if (length(ofile_envs) > 0L) {
      normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                    winslash = "/", mustWork = FALSE)
    } else getwd()
  }, error = function(e) getwd())

  .sh_path <- file.path(.td_dir, "stats_helpers.R")
  if (file.exists(.sh_path)) {
    source(.sh_path, local = FALSE)
  } else {
    message("Note: stats_helpers.R not found — using inline SW/Levene fallbacks.")
  }
  rm(.td_dir, .sh_path)
}

# Internal null-coalescing operator (only define if not already in environment)
if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b
