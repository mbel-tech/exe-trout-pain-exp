# =============================================================================
# 00_stats_helpers.R — Shared statistical helpers for the GOT pipeline
# =============================================================================
# Source this file before running any step script that needs these functions.
# 00_master_pipeline_GOT.R sources it automatically.
#
# Functions provided:
#   apply_transform()            — apply raw / log / log10 / sqrt / box_cox
#   run_boxcox_grid()            — grid search for optimal Box-Cox lambda (SW-p)
#   screen_transformation()      — choose best transform per outcome vector
#   gate_fish_id_by_bl_tier()    — filter trials to those meeting MIN_BL_RELIABILITY_TIER
#   aggregate_trial_level()      — aggregate fish×TP data to one row per trial
#   add_bh_column()              — add Benjamini-Hochberg adjusted p-value column
#   fit_dharma_summary()         — run DHARMa diagnostics and return a summary df
#   add_r2_row()                 — append R²m / R²c (MuMIn) row to an anova table
# =============================================================================


# ---- Transformation utilities ------------------------------------------------

#' Apply a named transformation to a numeric vector.
#' @param y      numeric vector (must be > 0 for log / box_cox)
#' @param method one of "raw", "log", "log10", "sqrt", "box_cox"
#' @param lambda Box-Cox lambda (required when method == "box_cox")
apply_transform <- function(y, method, lambda = NULL) {
  if (method == "raw")   return(y)
  if (method == "log10") return(log10(y))
  if (method == "log")   return(log(y))
  if (method == "sqrt")  return(sqrt(y))
  if (method == "box_cox") {
    if (is.null(lambda) || !is.finite(lambda))
      stop("box_cox requires a finite lambda")
    return(if (abs(lambda) < 0.01) log(y) else (y^lambda - 1) / lambda)
  }
  stop("Unknown transformation: ", method)
}

#' Grid search for the Box-Cox lambda that maximises Shapiro-Wilk p.
#' Returns list(lambda, sw_p). Both NA when insufficient data or all y <= 0.
run_boxcox_grid <- function(y, lambdas = seq(-2, 2, by = 0.1)) {
  y <- y[is.finite(y) & y > 0]
  if (length(y) < 3L) return(list(lambda = NA_real_, sw_p = NA_real_))
  best_lambda <- NA_real_
  best_p      <- -Inf
  for (lam in lambdas) {
    yt <- if (abs(lam) < 0.01) log(y) else (y^lam - 1) / lam
    sw <- tryCatch(shapiro.test(yt), error = function(e) NULL)
    if (!is.null(sw) && is.finite(sw$p.value) && sw$p.value > best_p) {
      best_p      <- sw$p.value
      best_lambda <- lam
    }
  }
  list(lambda = best_lambda,
       sw_p   = if (is.finite(best_p)) best_p else NA_real_)
}

#' Choose the best transformation for a numeric outcome vector.
#'
#' Candidates: raw, sqrt (if all y >= 0), log / log10 (if all y > 0), box_cox.
#' Box-Cox wins over log only when SW-p(box_cox) > SW-p(log) * 1.5 AND |lambda| >= 0.1.
#'
#' @param y          numeric vector
#' @param outcome_nm name label used in the returned summary row
#' @return list with:
#'   selected  — character method name
#'   bc_lambda — numeric lambda (NA unless box_cox selected)
#'   sw_p      — Shapiro-Wilk p for the selected transform
#'   table     — data.frame summary of all candidates
screen_transformation <- function(y, outcome_nm = "outcome") {
  y_clean <- y[is.finite(y)]

  sw_p_of <- function(vec) {
    vec <- vec[is.finite(vec)]
    if (length(vec) < 3L) return(NA_real_)
    tryCatch(shapiro.test(vec)$p.value, error = function(e) NA_real_)
  }

  candidates <- list()
  candidates[["raw"]]  <- list(y = y_clean, lambda = NA_real_)
  if (all(y_clean >= 0))
    candidates[["sqrt"]] <- list(y = sqrt(y_clean), lambda = NA_real_)
  if (all(y_clean > 0)) {
    candidates[["log"]]   <- list(y = log(y_clean),   lambda = NA_real_)
    candidates[["log10"]] <- list(y = log10(y_clean), lambda = NA_real_)
    bc <- run_boxcox_grid(y_clean)
    if (!is.na(bc$lambda))
      candidates[["box_cox"]] <- list(y = apply_transform(y_clean, "box_cox", bc$lambda),
                                      lambda = bc$lambda)
  }

  sw_ps <- vapply(candidates, function(c) sw_p_of(c$y), numeric(1))

  tbl <- data.frame(
    outcome   = outcome_nm,
    method    = names(candidates),
    sw_p      = sw_ps,
    bc_lambda = vapply(candidates, function(c) c$lambda, numeric(1)),
    selected  = FALSE,
    stringsAsFactors = FALSE
  )

  # Preference rule: box_cox beats log only when SW-p ratio >= 1.5 AND |lambda| >= 0.1
  best_nm <- names(which.max(sw_ps))
  if (!is.na(sw_ps["box_cox"]) && !is.na(sw_ps["log"]) &&
      is.finite(sw_ps["box_cox"]) && is.finite(sw_ps["log"])) {
    bc_lambda_val <- candidates[["box_cox"]]$lambda
    if (sw_ps["box_cox"] > sw_ps["log"] * 1.5 && abs(bc_lambda_val) >= 0.1) {
      best_nm <- "box_cox"
    } else if (best_nm == "box_cox") {
      # box_cox didn't beat log by 1.5x — fall back to log
      best_nm <- "log"
    }
  }

  tbl$selected[tbl$method == best_nm] <- TRUE
  bc_lambda_sel <- if (best_nm == "box_cox") candidates[["box_cox"]]$lambda else NA_real_

  list(
    selected  = best_nm,
    bc_lambda = bc_lambda_sel,
    sw_p      = sw_ps[best_nm],
    table     = tbl
  )
}


# ---- BL-tier gating ----------------------------------------------------------

# Reliability tier ordering (ascending)
.BL_TIER_ORDER <- c("unreliable", "borderline", "reliable", "conservative")

#' R4: Resolve the BL reliability tier CSV path with stable canonical fallback.
#'
#' Lookup order:
#'   1. Explicit `path` argument (if provided and exists)
#'   2. BL_MATCH_PATH from .GlobalEnv (if set and exists)
#'   3. Canonical: PIPELINE_DATA_DIR/master_outputs/reliability/reliability_tier_LATEST.csv
#'      (written by STEP 1.5 and/or STEP 1.6 — most recent producer wins)
#'   4. NA_character_ if nothing is found.
#'
#' This allows STEP 5/6 + plasma/monoamine/aggression side scripts to read the
#' BL tier table in a fresh R session without having to first run STEP 1.5/1.6.
resolve_bl_match_path <- function(path = NULL) {
  if (!is.null(path) && length(path) == 1L && !is.na(path) && nzchar(path) &&
      file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  if (exists("BL_MATCH_PATH", envir = .GlobalEnv, inherits = FALSE)) {
    p <- get("BL_MATCH_PATH", envir = .GlobalEnv)
    if (!is.null(p) && length(p) == 1L && !is.na(p) && nzchar(p) &&
        file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    pdd <- get("PIPELINE_DATA_DIR", envir = .GlobalEnv)
    canonical <- file.path(pdd, "master_outputs", "reliability",
                           "reliability_tier_LATEST.csv")
    if (file.exists(canonical)) {
      return(normalizePath(canonical, winslash = "/", mustWork = FALSE))
    }
  }
  NA_character_
}

#' Return the set of trial_ids that meet or exceed the minimum BL match tier.
#'
#' @param bl_match_df  data.frame with at least columns trial_id and a tier column
#' @param min_tier     character: one of "unreliable","borderline","reliable","conservative"
#' @param tier_col     name of the column containing tier values (default "reliability_tier").
#'                     Use "hungarian_tier" for side-workflow gating,
#'                     "idmatcher_tier" for core STEP5/STEP6 RE gating.
#' @return integer vector of qualifying trial_ids
qualifying_trial_ids <- function(bl_match_df, min_tier = NULL,
                                 tier_col = "reliability_tier") {
  if (is.null(min_tier)) {
    min_tier <- if (exists("MIN_BL_RELIABILITY_TIER", envir = .GlobalEnv, inherits = FALSE))
      get("MIN_BL_RELIABILITY_TIER", envir = .GlobalEnv)
    else "reliable"
  }
  tier_idx <- match(tolower(trimws(min_tier)), .BL_TIER_ORDER)
  if (is.na(tier_idx)) stop("Unknown BL tier: ", min_tier,
                             ". Use one of: ", paste(.BL_TIER_ORDER, collapse = ", "))

  if (!(tier_col %in% names(bl_match_df)))
    stop("bl_match_df must contain a '", tier_col, "' column.")
  if (!("trial_id" %in% names(bl_match_df)))
    stop("bl_match_df must contain a 'trial_id' column.")

  tier_vals <- match(tolower(trimws(bl_match_df[[tier_col]])), .BL_TIER_ORDER)
  bl_match_df$trial_id[!is.na(tier_vals) & tier_vals >= tier_idx]
}


#' Gate idmatcher RE candidates: replace (1|fish_id) with (1|fish_id_matched_trial)
#' for trials that pass the idmatcher tier, keep (1|trial_id) only for the rest.
#'
#' @param candidate_list  named list of lme4 formula strings (the full RE candidate set)
#' @param bl_match_df     merged tier data.frame (from STEP 1.5); must have
#'                        trial_id + idmatcher_tier columns
#' @param qualifying_ids  integer vector of trial_ids that pass the idmatcher tier
#'                        (pre-computed by qualifying_trial_ids(..., tier_col="idmatcher_tier"))
#' @return named list — for qualifying trials: fish_id candidates replaced by
#'         fish_id_matched_trial candidates; for non-qualifying: (1|trial_id) only
gate_random_candidates_idmatcher <- function(candidate_list, bl_match_df,
                                             qualifying_ids) {
  out <- lapply(candidate_list, function(form_str) {
    # Replace every occurrence of fish_id (as an RE grouping factor) with
    # fish_id_matched_trial when the trial qualifies; otherwise remove entirely.
    if (grepl("\\bfish_id\\b", form_str)) {
      if (length(qualifying_ids) > 0) {
        gsub("\\bfish_id\\b", "fish_id_matched_trial", form_str)
      } else {
        NULL   # drop fish_id candidates when no trial qualifies
      }
    } else {
      form_str
    }
  })
  # Remove NULLs
  Filter(Negate(is.null), out)
}

#' TRUE/FALSE indicator: does this data.frame contain fish_id RE-eligible trials?
#' Returns a logical vector aligned with df's rows.
#'
#' @param df           data.frame that contains a trial_id column
#' @param bl_match_df  data.frame with trial_id + reliability_tier columns
#' @param min_tier     minimum tier (defaults to MIN_BL_RELIABILITY_TIER global)
rows_eligible_for_fish_re <- function(df, bl_match_df, min_tier = NULL) {
  ok_ids <- qualifying_trial_ids(bl_match_df, min_tier)
  df$trial_id %in% ok_ids
}


# ---- Trial-level aggregation -------------------------------------------------

#' Aggregate fish×TP data to one row per trial.
#'
#' Steps:
#'   1. Keep only TP2 and TP3 (or whatever tps_to_avg specifies).
#'   2. Average across fish within each (trial_id, timepoint).
#'   3. Average across timepoints within each trial_id.
#'
#' Returns a data.frame with one row per trial, columns:
#'   trial_id, Treatment, Group, trial_date, <outcome_cols>
#'
#' @param df           fish×TP data.frame
#' @param outcome_cols character vector of numeric outcome column names to aggregate
#' @param tps_to_avg   integer vector of timepoints to include (default c(2L, 3L))
aggregate_trial_level <- function(df, outcome_cols, tps_to_avg = c(2L, 3L)) {
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("dplyr is required for aggregate_trial_level()")

  tp_col <- if ("timepoint" %in% names(df)) "timepoint" else
    stop("df must contain a 'timepoint' column")

  # Normalise timepoint to integer for comparison
  df[[tp_col]] <- suppressWarnings(as.integer(df[[tp_col]]))

  # Step 1: filter to requested timepoints
  df_tp <- df[!is.na(df[[tp_col]]) & df[[tp_col]] %in% tps_to_avg, , drop = FALSE]

  if (nrow(df_tp) == 0)
    stop("No rows remain after filtering to timepoints: ", paste(tps_to_avg, collapse = ", "))

  id_cols <- intersect(c("trial_id", "Treatment", "Group", "trial_date"), names(df_tp))

  # Step 2: mean across fish within (trial_id, timepoint)
  agg_cols <- c(id_cols, tp_col, outcome_cols)
  agg_cols <- intersect(agg_cols, names(df_tp))
  df_tp    <- df_tp[, agg_cols, drop = FALSE]

  by_tp <- dplyr::group_by(df_tp, dplyr::across(dplyr::all_of(c("trial_id", tp_col))))
  by_tp <- dplyr::summarise(
    by_tp,
    dplyr::across(dplyr::all_of(
      intersect(c("Treatment", "Group", "trial_date"), names(df_tp))
    ), dplyr::first),
    dplyr::across(dplyr::all_of(outcome_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

  # Step 3: mean across timepoints within trial_id
  by_trial <- dplyr::group_by(by_tp, trial_id)
  by_trial <- dplyr::summarise(
    by_trial,
    dplyr::across(dplyr::all_of(
      intersect(c("Treatment", "Group", "trial_date"), names(by_tp))
    ), dplyr::first),
    dplyr::across(dplyr::all_of(outcome_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

  as.data.frame(by_trial)
}


# ---- Multiple-comparison correction -----------------------------------------

#' Add a Benjamini-Hochberg adjusted p-value column to an ANOVA summary table.
#'
#' @param anova_df   data.frame with a column named p_value (or p.value)
#' @param p_col      name of the raw p-value column (default "p_value")
#' @param out_col    name of the new BH column (default "p_BH")
#' @return anova_df with an additional column
add_bh_column <- function(anova_df, p_col = "p_value", out_col = "p_BH") {
  # Accept either p_value or p.value
  if (!(p_col %in% names(anova_df))) {
    alt <- if (p_col == "p_value") "p.value" else "p_value"
    if (alt %in% names(anova_df)) p_col <- alt else
      stop("Column '", p_col, "' not found in anova_df. Available: ",
           paste(names(anova_df), collapse = ", "))
  }
  ps <- anova_df[[p_col]]
  anova_df[[out_col]] <- p.adjust(ps, method = "BH")
  anova_df
}


# ---- DHARMa diagnostics ------------------------------------------------------

#' Run DHARMa simulated-residual diagnostics and return a one-row summary df.
#'
#' @param model      a fitted lmer / glmer / glmmTMB object
#' @param outcome_nm label for the 'outcome' column in the returned df
#' @param n_sim      number of simulations (default 500)
#' @param png_path   optional file path; if given, saves a 2-panel PNG
#' @return data.frame with columns: outcome, KS_stat, KS_p, dispersion_stat,
#'         dispersion_p, outlier_p
fit_dharma_summary <- function(model, outcome_nm = "outcome", n_sim = 500L,
                                png_path = NULL) {
  if (!requireNamespace("DHARMa", quietly = TRUE))
    stop("DHARMa is required for fit_dharma_summary()")

  sim <- tryCatch(
    DHARMa::simulateResiduals(fittedModel = model, n = n_sim, plot = FALSE),
    error = function(e) {
      warning("DHARMa simulation failed for '", outcome_nm, "': ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(sim)) {
    return(data.frame(
      outcome = outcome_nm, KS_stat = NA_real_, KS_p = NA_real_,
      dispersion_stat = NA_real_, dispersion_p = NA_real_, outlier_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  ks  <- DHARMa::testUniformity(sim, plot = FALSE)
  dsp <- DHARMa::testDispersion(sim, plot = FALSE)
  out <- DHARMa::testOutliers(sim, plot = FALSE, type = "binomial")

  if (!is.null(png_path)) {
    grDevices::png(png_path, width = 1600, height = 800, res = 150)
    DHARMa::plot(sim, title = outcome_nm)
    grDevices::dev.off()
  }

  data.frame(
    outcome         = outcome_nm,
    KS_stat         = as.numeric(ks$statistic),
    KS_p            = ks$p.value,
    dispersion_stat = as.numeric(dsp$statistic),
    dispersion_p    = dsp$p.value,
    outlier_p       = out$p.value,
    stringsAsFactors = FALSE
  )
}


# ---- R-squared (MuMIn) -------------------------------------------------------

#' Append a row with R²m and R²c (MuMIn::r.squaredGLMM) to an ANOVA table.
#'
#' @param anova_df  data.frame (ANOVA results for one model)
#' @param model     the fitted model object
#' @param term_col  name of the column that holds term names (default "term")
#' @return anova_df with an extra row for R²
add_r2_row <- function(anova_df, model, term_col = "term") {
  if (!requireNamespace("MuMIn", quietly = TRUE)) {
    warning("MuMIn not available — skipping R² row.")
    return(anova_df)
  }
  r2 <- tryCatch(
    MuMIn::r.squaredGLMM(model),
    error = function(e) NULL
  )
  if (is.null(r2)) return(anova_df)

  r2m <- if (!is.null(r2) && "R2m" %in% colnames(r2)) r2[1, "R2m"] else NA_real_
  r2c <- if (!is.null(r2) && "R2c" %in% colnames(r2)) r2[1, "R2c"] else NA_real_

  r2_row <- anova_df[NA_integer_, , drop = FALSE]  # same structure, all NA
  r2_row[1, ] <- NA
  if (term_col %in% names(r2_row)) r2_row[[term_col]] <- "R2m / R2c"
  if ("estimate" %in% names(r2_row)) {
    r2_row$estimate <- NA_real_
  }
  # Store R2 values in any available numeric columns
  extra_cols <- setdiff(names(r2_row), term_col)
  num_cols   <- extra_cols[vapply(r2_row[extra_cols], is.numeric, logical(1))]
  if (length(num_cols) >= 1) r2_row[[num_cols[1]]] <- r2m
  if (length(num_cols) >= 2) r2_row[[num_cols[2]]] <- r2c

  rbind(anova_df, r2_row)
}


# ---- Overdispersion check (Poisson) -----------------------------------------

#' Compute Pearson chi-squared / df for a fitted Poisson (or NB) model.
#' Returns a one-row data.frame with columns: outcome, pearson_ratio, n_obs, df_resid.
#' A ratio > 1.5 suggests overdispersion; consider NB or quasi-Poisson.
pearson_overdispersion <- function(model, outcome_nm = "outcome") {
  resid_pearson <- tryCatch(
    residuals(model, type = "pearson"),
    error = function(e) NULL
  )
  df_r <- tryCatch(df.residual(model), error = function(e) NA_real_)

  if (is.null(resid_pearson) || is.na(df_r) || df_r <= 0) {
    return(data.frame(outcome = outcome_nm, pearson_ratio = NA_real_,
                      n_obs = NA_integer_, df_resid = NA_real_,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    outcome      = outcome_nm,
    pearson_ratio = sum(resid_pearson^2, na.rm = TRUE) / df_r,
    n_obs        = length(resid_pearson),
    df_resid     = df_r,
    stringsAsFactors = FALSE
  )
}
