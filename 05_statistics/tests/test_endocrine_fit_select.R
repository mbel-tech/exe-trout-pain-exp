# Tests for fit_select() in the endocrine engine: a singular random-effect
# structure must be REDUCED (dropped and refitted), never retained.
# Run: Rscript "05_statistics/tests/test_endocrine_fit_select.R"
here <- (function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)]))) })()
source(file.path(here, "helpers.R"))
suppressMessages({ library(lme4); library(lmerTest); library(MuMIn); library(dplyr); library(tibble); library(car) })
options(contrasts = c("contr.sum", "contr.poly"))
suppressMessages(source(file.path(here, "..", "shared", "house_style.R")))
ENDO <- file.path(here, "..", "plasma_monoamine_models.R")
for (fn in c(".pick_start", ".reduction_result", ".is_singular_gamma")) invisible(load_fn(ENDO, fn))
fit_select <- load_fn(file.path(here, "..", "plasma_monoamine_models.R"), "fit_select")

RE <- list(forced_design = "(1 | date_uid) + (1 | trial_uid)")
fs <- function(d) suppressMessages(suppressWarnings(fit_select(d, "y", "Treatment * Group", RE)))
is_singular_lmer <- function(m) inherits(m, "merMod") && lme4::isSingular(m)

# 1. no between-trial / between-date variance: both design REs are dropped -> lm,
#    and the indicator must still produce Type III inference (it must not vanish).
res <- fs(make_design_data(trial_sd = 0))
stopifnot(isTRUE(res$is_lm), inherits(res$model, "lm"), !inherits(res$model, "merMod"),
          length(res$dropped) == 2L, length(res$reasons) == 2L,
          identical(res$re_formula, "(none)"), isFALSE(res$singular))
an <- car::Anova(res$model, type = 3)
stopifnot("Treatment:Group" %in% rownames(an), is.finite(an["Treatment:Group", "Pr(>F)"]))

# 2. a real trial variance: a mixed model is kept, it is NOT singular, trial_uid survives.
res <- fs(make_design_data(trial_sd = 3, seed = 2))
stopifnot(!isTRUE(res$is_lm), inherits(res$model, "merMod"), !is_singular_lmer(res$model),
          grepl("trial_uid", res$re_formula, fixed = TRUE))

# 3. whatever the data, a singular lmer is never returned.
for (s in 1:20) {
  res <- fs(make_design_data(trial_sd = 0.3, seed = s))
  stopifnot(!is_singular_lmer(res$model), !is.null(res$model))
}

# ---- Gamma and Tweedie selectors: same rule -----------------------------------
suppressMessages(library(glmmTMB))
fit_select_gamma   <- load_fn(file.path(here, "..", "plasma_monoamine_models.R"), "fit_select_gamma")
fit_select_tweedie <- load_fn(file.path(here, "..", "plasma_monoamine_models.R"), "fit_select_tweedie")
fg <- function(d) suppressMessages(suppressWarnings(fit_select_gamma(d, "y", "Treatment * Group", RE)))
ft <- function(d) suppressMessages(suppressWarnings(fit_select_tweedie(d, "y", "Treatment * Group", RE)))

# 4. Gamma, no between-trial/date variance: both REs dropped -> a plain Gamma glm
#    (never NULL, never a singular GLMM), and the drops are recorded.
res <- fg(make_positive_design_data(trial_sd = 0))
stopifnot(!is.null(res), isTRUE(res$is_glm), inherits(res$model, "glm"),
          !inherits(res$model, "glmmTMB"), length(res$dropped) == 2L,
          identical(res$re_formula, "(none)"), isFALSE(res$singular),
          !is.null(car::Anova(res$model, type = 3)))
# 5. Gamma, real trial variance: a GLMM is kept, it is not singular, trial_uid survives.
res <- fg(make_positive_design_data(trial_sd = 0.6, seed = 3))
stopifnot(!isTRUE(res$is_glm), inherits(res$model, "glmmTMB"),
          !.is_singular_any(res$model), grepl("trial_uid", res$re_formula, fixed = TRUE))
# 6. Tweedie: same two cases; the no-RE terminal model is a glmmTMB with no random part.
res <- ft(make_positive_design_data(trial_sd = 0))
stopifnot(!is.null(res), isTRUE(res$is_glm), inherits(res$model, "glmmTMB"),
          length(re_sd_components(res$model)) == 0L, length(res$dropped) == 2L,
          identical(res$re_formula, "(none)"))
res <- ft(make_positive_design_data(trial_sd = 0.6, seed = 3))
stopifnot(!isTRUE(res$is_glm), inherits(res$model, "glmmTMB"), !.is_singular_any(res$model))

cat("PASS: test_endocrine_fit_select.R\n")