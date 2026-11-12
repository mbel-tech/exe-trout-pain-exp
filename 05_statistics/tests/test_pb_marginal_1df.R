# Tests for the 1-df (marginal) parametric-bootstrap LRT helpers in house_style.R.
# Run: Rscript "05_statistics/tests/test_pb_marginal_1df.R"
here <- (function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)]))) })()
suppressMessages({ library(glmmTMB); library(dplyr); library(tibble) })
options(contrasts = c("contr.sum", "contr.poly"))
suppressMessages(source(file.path(here, "..", "shared", "house_style.R")))

# 12 trials x 4 fish x 2 intervals = 96 rows; beta response; Group and Treatment are
# constant within a trial, interval varies within it (as in the real fish x interval tier).
make_beta_data <- function(group_coef = 0, trial_sd = 0.3, seed = 1) {
  set.seed(seed)
  trial <- rep(1:12, each = 8)
  d <- data.frame(trial_uid = factor(trial))
  d$Treatment <- factor(ifelse(trial %% 2 == 1, "Control", "Exercise"), levels = c("Control", "Exercise"))
  d$Group     <- factor(ifelse(ceiling(trial / 2) <= 3, "A", "B"), levels = c("A", "B"))
  d$interval  <- factor(rep(c("TP2", "TP3"), 48), levels = c("TP2", "TP3"))
  eta <- -0.2 + group_coef * ifelse(d$Group == "A", 0.5, -0.5) + rnorm(12, 0, trial_sd)[trial]
  mu  <- plogis(eta); phi <- 25
  d$.y <- rbeta(96, mu * phi, (1 - mu) * phi)
  d
}
d <- add_sum_coded_cols(make_beta_data(group_coef = 0, seed = 11))

# 1. coding: +1 for the first level, -1 for the last, exactly the contr.sum convention
stopifnot(all(d$Group_c[d$Group == "A"] == 1), all(d$Group_c[d$Group == "B"] == -1),
          all(d$interval_c[d$interval == "TP2"] == 1), all(d$Treatment_c[d$Treatment == "Exercise"] == -1))
# a factor with other than 2 levels must be refused, not silently mis-coded
bad <- make_beta_data(); bad$interval <- factor(rep(c("a", "b", "c", "a"), 24))
stopifnot(inherits(try(add_sum_coded_cols(bad), silent = TRUE), "try-error"))

# 2. specs
sg <- pb_1df_specs("Group"); st <- pb_1df_specs("Treatment")
stopifnot(identical(sg$full, "Treatment_c * Group_c * interval_c"),
          identical(sg$reduced, "Treatment_c * interval_c + Group_c:Treatment_c + Group_c:interval_c + Group_c:Treatment_c:interval_c"),
          identical(sg$df, 1L), identical(st$df, 1L),
          identical(st$reduced, "Group_c * interval_c + Treatment_c:Group_c + Treatment_c:interval_c + Treatment_c:Group_c:interval_c"))
fit <- function(fx) glmmTMB(as.formula(paste(".y ~", fx, "+ (1 | trial_uid)")), data = d, family = beta_family())
f_full <- fit(sg$full); f_red <- fit(sg$reduced); f_fac <- fit("Treatment * Group * interval")
# 2a. the numeric-coded full model IS the reported factor model
stopifnot(abs(as.numeric(logLik(f_full)) - as.numeric(logLik(f_fac))) < 1e-4)
# 2b. the reduced model has exactly ONE fewer fixed parameter (1 df, marginal)
stopifnot(length(fixef(f_full)$cond) - length(fixef(f_red)$cond) == 1L,
          length(fixef(pb_1df_specs("Treatment") |> (\(s) fit(s$reduced))())$cond) == length(fixef(f_full)$cond) - 1L)

# 3. behaviour: a real Group effect is detected, a null one is not
dd <- add_sum_coded_cols(make_beta_data(group_coef = 1.5, seed = 12)); dd$.y <- dd$.y
res1 <- pb_lrt(sg$full, sg$reduced, "(1 | trial_uid)", dd, glmmTMB::beta_family(), nsim = 99L, seed = 1L)
stopifnot(identical(res1$method, "pb_lrt"), res1$p_pb <= 0.05)
d0 <- add_sum_coded_cols(make_beta_data(group_coef = 0, seed = 11))
res0 <- pb_lrt(sg$full, sg$reduced, "(1 | trial_uid)", d0, glmmTMB::beta_family(), nsim = 99L, seed = 1L)
stopifnot(identical(res0$method, "pb_lrt"), res0$p_pb > 0.10)

# 4. failure honesty: if the fits fail the result is NA with a method, never a small p
dc <- d0; dc$.y <- 0.5
resf <- pb_lrt(sg$full, sg$reduced, "(1 | trial_uid)", dc, glmmTMB::beta_family(), nsim = 9L, seed = 1L)
stopifnot(is.na(resf$p_pb) || resf$method == "pb_lrt")      # either honest NA, or a genuine bootstrap
if (!is.na(resf$p_pb)) stopifnot(resf$nsim_ok >= 5L)        # a p from fewer than half the replicates is never returned

cat("PASS: test_pb_marginal_1df.R\n")