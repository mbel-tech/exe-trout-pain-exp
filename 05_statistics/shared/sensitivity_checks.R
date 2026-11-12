# =============================================================================
# sensitivity_checks.R -- the C12 robustness-artifact battery for the GOT
# manuscript chain.
# =============================================================================
# One function per named check in handoff/HANDOFF.md S7 / CONVENTIONS.md C12.
# Each writes its FIXED filename into the caller's run directory, whether or
# not it changed anything -- "a check with no artifact did not happen" (C12
# [INV]). Source house_style.R first (fmt_p_boot, effect_size, pb_lrt
# all live there).
#
# Two functions here are PROMOTIONS of code that already existed in the GOT
# chain, not new inventions:
#   - sensitivity_loo()      <- loo_trial_stability() in aggression_models.R
#   - metric_collinearity()  <- the Spearman screen documented in
#                                link_behaviour_physiology.R (a generic re-
#                                implementation; that script itself is out of
#                                scope for this refactor and is untouched)
# =============================================================================

if (!exists("GOT_ROOT")) GOT_ROOT <- "."

.write_sens_csv <- function(df, filename, run_dir) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(run_dir, filename), row.names = FALSE)
  invisible(df)
}

# =============================================================================
# aicc_selection.csv -- RE-candidate selection with singular/eligible columns
# (A11 apply_singularity_guard() in house_style.R does the actual
# guard logic; this just standardises the CSV shape all six fitters emit).
# =============================================================================
sensitivity_aicc_selection <- function(cand_results, aicc_vals, guard, selected_name,
                                       label, run_dir) {
  df <- data.frame(
    label      = label,
    RE         = vapply(cand_results, `[[`, character(1), "name"),
    re_formula = vapply(cand_results, `[[`, character(1), "re"),
    AICc       = round(aicc_vals, 2),
    delta_AICc = round(aicc_vals - min(aicc_vals[guard$eligible], na.rm = TRUE), 2),
    singular   = guard$singular,
    eligible   = guard$eligible,
    selected   = vapply(cand_results, `[[`, character(1), "name") == selected_name,
    all_singular_fallback = guard$all_singular
  )
  .write_sens_csv(df, "aicc_selection.csv", run_dir)
}

# =============================================================================
# sensitivity_forced_re.csv -- add the design RE back to the AICc winner (if
# AICc dropped it); compare F/p/ICC/singular side by side. Never touches
# fixed effects.
# =============================================================================
sensitivity_forced_re <- function(label, term, winning_re, forced_re,
                                  winning_p, forced_p, winning_icc, forced_icc,
                                  winning_singular, forced_singular, run_dir) {
  df <- data.frame(
    label = label, term = term,
    winning_re = winning_re, forced_re = forced_re,
    p_winning = winning_p, p_forced = forced_p,
    icc_winning = winning_icc, icc_forced = forced_icc,
    singular_winning = winning_singular, singular_forced = forced_singular,
    p_crosses_alpha = (winning_p < 0.05) != (forced_p < 0.05)
  )
  .write_sens_csv(df, "sensitivity_forced_re.csv", run_dir)
}

# =============================================================================
# sensitivity_loo.csv -- leave-one-trial-out stability of fixed effects.
# Promoted verbatim from aggression_models.R's loo_trial_stability(),
# generalised to any model with a fixef()-like accessor (glmmTMB or
# lme4/lmerTest).
# =============================================================================
sensitivity_loo <- function(model, data = NULL, trial_col = "trial_id", label = "", run_dir) {
  get_fixef <- function(m) {
    if (inherits(m, "glmmTMB")) return(glmmTMB::fixef(m)$cond)
    stats::coef(summary(m))[, 1]
  }
  full_est <- tryCatch(get_fixef(model), error = function(e) NULL)
  if (is.null(full_est)) return(invisible(NULL))

  # Refit from the MODEL'S OWN frame, not from a caller-supplied data frame.
  # The fitters build their response into a ".y" column on a per-variable
  # working copy, so update(model, data = <the original dataset>) fails on
  # every replicate -- silently, because the failure is caught and turned into
  # NA. That produced a sensitivity_loo.csv full of Inf/-Inf which LOOKED like
  # a completed check: worse than an absent artifact, since C12's rule is that
  # a check with no artifact did not happen, and this one had an artifact and
  # still did not happen. Refitting from the frame keeps the response,
  # transform, and row exclusions identical by construction.
  d_all <- tryCatch(
    if (inherits(model, "glmmTMB")) model$frame else model@frame,
    error = function(e) NULL)
  if (is.null(d_all) && !is.null(data)) d_all <- data
  if (is.null(d_all) || !trial_col %in% names(d_all)) {
    # The grouping column is not in the frame (it was not a model term), so a
    # leave-one-trial-out refit is not constructible. Say so in the artifact.
    df <- data.frame(label = label, coef = names(full_est),
                     full_est = as.numeric(full_est),
                     min_loo = NA_real_, max_loo = NA_real_, range_loo = NA_real_,
                     n_refits_ok = 0L, n_trials = NA_integer_,
                     result = sprintf("not computable: '%s' absent from the model frame",
                                      trial_col),
                     row.names = NULL)
    return(.write_sens_csv(df, "sensitivity_loo.csv", run_dir))
  }

  trials <- unique(d_all[[trial_col]])
  loo <- lapply(trials, function(t) {
    d <- d_all[d_all[[trial_col]] != t, , drop = FALSE]
    m <- tryCatch(stats::update(model, data = d), error = function(e) NULL)
    if (is.null(m)) return(rep(NA_real_, length(full_est)))
    est <- tryCatch(get_fixef(m), error = function(e) NULL)
    if (is.null(est) || length(est) != length(full_est))
      return(rep(NA_real_, length(full_est)))
    est
  })
  mat <- do.call(rbind, loo)
  n_ok <- sum(stats::complete.cases(mat))
  safe <- function(x, f) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else f(x) }
  df <- data.frame(
    label = label, coef = names(full_est), full_est = as.numeric(full_est),
    min_loo   = apply(mat, 2, safe, min),
    max_loo   = apply(mat, 2, safe, max),
    range_loo = apply(mat, 2, function(x) { x <- x[is.finite(x)]
                                            if (length(x) < 2) NA_real_ else diff(range(x)) }),
    n_refits_ok = n_ok, n_trials = length(trials),
    # A sign change across leave-one-out refits means the direction of the
    # estimate depends on a single trial.
    sign_change = apply(mat, 2, function(x) { x <- x[is.finite(x)]
                                              length(x) > 0 && any(x > 0) && any(x < 0) }),
    result = if (n_ok == 0) "not computable: every refit failed"
             else sprintf("%d of %d refits converged", n_ok, length(trials)),
    row.names = NULL
  )
  if (n_ok == 0)
    warning("sensitivity_loo: every leave-one-out refit failed for ", label,
            " -- artifact records the failure rather than an empty range.",
            call. = FALSE)
  .write_sens_csv(df, "sensitivity_loo.csv", run_dir)
}

# =============================================================================
# multiverse_<endpoint>.csv -- specification-curve sweep for the PRIMARY
# endpoint only. response scale x aggregation level x RE structure, all
# combinations that are actually fittable.
# =============================================================================
sensitivity_multiverse <- function(endpoint, data, response_var, fixed_str,
                                   re_candidates, fit_fn, run_dir) {
  agg_levels <- if (!is.null(data$.agg_level)) unique(data$.agg_level) else "single"
  rows <- list()
  for (rl in names(re_candidates)) {
    m <- tryCatch(fit_fn(data, response_var, fixed_str, re_candidates[[rl]]),
                  error = function(e) NULL)
    if (is.null(m)) next
    an <- tryCatch(as.data.frame(car::Anova(m, type = 3)), error = function(e) NULL)
    p_treat <- if (!is.null(an) && "Pr(>Chisq)" %in% names(an) && "treatment_clean" %in% rownames(an))
      an["treatment_clean", "Pr(>Chisq)"] else NA_real_
    rows[[rl]] <- data.frame(endpoint = endpoint, re_spec = rl, p_treatment = p_treat)
  }
  df <- if (length(rows)) do.call(rbind, rows) else
    data.frame(endpoint = endpoint, re_spec = NA_character_, p_treatment = NA_real_)
  .write_sens_csv(df, sprintf("multiverse_%s.csv", endpoint), run_dir)
}

# =============================================================================
# equivalence_tost.csv -- TOST against a declared SESOI. SESOI must be
# supplied by the caller (declared BEFORE looking at the data -- HANDOFF S7
# [INV]); this function does not invent one. Uses a Welch-style TOST on the
# standardised mean difference via effectsize/emmeans contrast SEs.
# =============================================================================
sensitivity_tost <- function(label, est, se, df, sesoi_d, sesoi_source, run_dir) {
  if (!is.finite(est) || !is.finite(se) || is.na(df)) {  # df = Inf is a valid asymptotic (z) test
    df_out <- data.frame(label = label, sesoi_d = sesoi_d, sesoi_source = sesoi_source,
                         result = "not computable", reason = "missing est/se/df")
    return(.write_sens_csv(df_out, "equivalence_tost.csv", run_dir))
  }
  t_lo <- (est - (-sesoi_d)) / se
  t_hi <- (est - sesoi_d) / se
  p_lo <- stats::pt(t_lo, df, lower.tail = FALSE)
  p_hi <- stats::pt(t_hi, df, lower.tail = TRUE)
  p_tost <- max(p_lo, p_hi)
  df_out <- data.frame(
    label = label, estimate = est, se = se, df = df,
    sesoi_d = sesoi_d, sesoi_source = sesoi_source,
    p_tost = p_tost,
    result = ifelse(p_tost < 0.05, "equivalent (confirmed absence)", "inconclusive")
  )
  .write_sens_csv(df_out, "equivalence_tost.csv", run_dir)
}

# =============================================================================
# pseudoreplication_check.csv -- refit at the unit level AND the level above
# it; report both beta, both p, the beta ratio, and a direction-agree flag.
# =============================================================================
sensitivity_pseudoreplication <- function(label, term, beta_unit, p_unit,
                                          beta_agg, p_agg, run_dir) {
  df <- data.frame(
    label = label, term = term,
    beta_unit_level = beta_unit, p_unit_level = p_unit,
    beta_agg_level  = beta_agg,  p_agg_level  = p_agg,
    beta_ratio = if (is.finite(beta_agg) && beta_agg != 0) beta_unit / beta_agg else NA_real_,
    direction_agree = sign(beta_unit) == sign(beta_agg),
    reportable = (sign(beta_unit) == sign(beta_agg)) &&
                 ((p_unit < 0.05) == (p_agg < 0.05))
  )
  .write_sens_csv(df, "pseudoreplication_check.csv", run_dir)
}

# =============================================================================
# metric_collinearity.csv -- every secondary metric screened against the
# PRIMARY endpoint by Spearman |r|. Generic re-implementation of the manual
# screen documented in link_behaviour_physiology.R:415-449 (Pearson gave a
# materially wrong answer there -- rank-identical metrics looked distinct
# under Pearson; this function defaults to Spearman for the same reason).
# =============================================================================
metric_collinearity <- function(data, primary_endpoint, candidate_vars,
                                threshold = 0.70, method = "spearman", run_dir) {
  rows <- lapply(candidate_vars, function(v) {
    if (!v %in% names(data) || !primary_endpoint %in% names(data)) return(NULL)
    ok <- is.finite(data[[v]]) & is.finite(data[[primary_endpoint]])
    if (sum(ok) < 3) return(data.frame(metric = v, primary = primary_endpoint,
                                       r = NA_real_, method = method,
                                       flagged = NA, n = sum(ok)))
    r <- suppressWarnings(stats::cor(data[[v]][ok], data[[primary_endpoint]][ok], method = method))
    data.frame(metric = v, primary = primary_endpoint, r = r, method = method,
              flagged = is.finite(r) && abs(r) >= threshold, n = sum(ok))
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  .write_sens_csv(df, "metric_collinearity.csv", run_dir)
}

# =============================================================================
# design_balance.csv -- design-confound screen, run BEFORE any model is
# fitted. Flags factor<->factor 1:1 mapping, |cor| > 0.9, unintended nesting.
# GOT-specific: this is where trial_date <-> tank aliasing (LOGS.md 2026-05)
# must surface -- see plan Open question 4.
# =============================================================================
sensitivity_design_balance <- function(data, factor_cols, label = "", run_dir) {
  factor_cols <- intersect(factor_cols, names(data))
  combos <- utils::combn(factor_cols, 2, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    a <- as.character(data[[pair[1]]]); b <- as.character(data[[pair[2]]])
    tab <- table(a, b)
    # 1:1 mapping: every row of the table has exactly one nonzero cell
    one_to_one <- all(rowSums(tab > 0) <= 1) || all(colSums(tab > 0) <= 1)
    n_a <- length(unique(a)); n_b <- length(unique(b))
    data.frame(label = label, factor_1 = pair[1], factor_2 = pair[2],
              levels_1 = n_a, levels_2 = n_b,
              one_to_one_mapping = one_to_one,
              flag = one_to_one && n_a > 1 && n_b > 1)
  })
  df <- if (length(rows)) do.call(rbind, rows) else
    data.frame(label = label, factor_1 = NA, factor_2 = NA, levels_1 = NA,
              levels_2 = NA, one_to_one_mapping = NA, flag = NA)
  .write_sens_csv(df, "design_balance.csv", run_dir)
}

# =============================================================================
# stat_type_validation.csv -- confirms the string in anova_stmt_*.csv matches
# the test actually run in anova.csv. See validate_stat_types.R, which is the
# script-level driver that calls this per run directory.
# =============================================================================
sensitivity_stat_type_validation <- function(anova_df, run_dir) {
  # anova_df must carry: variable, term, stat_type ("F"|"chisq"|"pb_lrt"),
  # statement (the printed caption string).
  df <- anova_df
  df$check_F     <- ifelse(df$stat_type == "F",     grepl("^[^\\u03c7]*F", df$statement), NA)
  df$check_chisq <- ifelse(df$stat_type == "chisq", grepl("\u03c7", df$statement) & !grepl("\u03b7", df$statement), NA)
  df$check_boot  <- ifelse(df$stat_type == "pb_lrt", grepl("\u2264|N =", df$statement), NA)
  df$mismatch    <- !(coalesce3(df$check_F, df$check_chisq, df$check_boot))
  .write_sens_csv(df, "stat_type_validation.csv", run_dir)
}
coalesce3 <- function(a, b, c) ifelse(!is.na(a), a, ifelse(!is.na(b), b, ifelse(!is.na(c), c, TRUE)))
