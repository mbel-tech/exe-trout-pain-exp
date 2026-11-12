################################################################################
# aggression_models_design_random_effects.R
# Refit of the primary aggression count models under the design random effects
# used for behaviour and endocrine outcomes: (1 | date) + (1 | trial), reduced
# where singular. All 12 trials are present in aggression_annotated_june2026.csv
# (trial 8 = 49 events), so no zero-event trial needs re-injecting.
# Also fits a design-honest comparator: one row per trial, log rate, F on 8 df.
# Usage: Rscript 05_statistics/aggression_models_design_random_effects.R [<data dir> <output dir>]
################################################################################
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(glmmTMB)
  library(car); library(emmeans); library(DHARMa)})
options(contrasts = c("contr.sum", "contr.poly"))
set.seed(20261111)
args <- commandArgs(trailingOnly = TRUE)
DATA <- if (length(args) >= 1) args[1] else "data"
OUT  <- if (length(args) >= 2) args[2] else file.path("output", "aggression_design_re")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ev <- read.csv(file.path(DATA, "aggression_annotated_june2026.csv"))
roster <- ev %>% distinct(trial, group, treatment = treatment_clean, date)
stopifnot(nrow(roster) == 12)
grid <- roster %>% tidyr::crossing(segment = c(2, 3), side = c("calm", "flow"))
cnt <- ev %>% count(trial, segment, side, name = "n")
d <- grid %>% left_join(cnt, by = c("trial", "segment", "side")) %>%
  mutate(n = coalesce(n, 0L), off = log(600),
         group = factor(group, levels = c("No choice", "Exercise choice")),
         treatment = factor(treatment, levels = c("Control", "Exercise")),
         side = factor(side), segment = factor(segment),
         trial = factor(trial), date = factor(date))
write.csv(d, file.path(OUT, "grid_trial_interval_zone.csv"), row.names = FALSE)

fit_reduce <- function(fixed, data, family = nbinom2) {
  res <- list(c("(1 | date)", "(1 | trial)"), "(1 | trial)", "(1 | date)", character(0))
  for (re in res) {
    f <- as.formula(paste(fixed, if (length(re)) paste("+", paste(re, collapse = " + ")) else ""))
    m <- try(suppressWarnings(glmmTMB(f, data = data, family = family)), silent = TRUE)
    if (inherits(m, "try-error")) next
    ok <- isTRUE(m$sdr$pdHess) && m$fit$convergence == 0
    sds <- if (length(re)) unlist(lapply(VarCorr(m)$cond, function(v) attr(v, "stddev"))) else numeric(0)
    if (ok && all(sds > 1e-4)) return(list(m = m, re = if (length(re)) paste(re, collapse = " + ") else "(none)"))
  }
  stop("no structure converged")
}
report <- function(tag, fixed, data, specs) {
  # family: Poisson first, NB2 if overdispersed on the fitted Poisson
  fp <- fit_reduce(fixed, data, poisson)
  disp <- DHARMa::testDispersion(simulateResiduals(fp$m, n = 1000, plot = FALSE), plot = FALSE)
  fam <- "poisson"; f <- fp
  if (disp$p.value < 0.05) { f <- fit_reduce(fixed, data, nbinom2); fam <- "nbinom2" }
  sim <- simulateResiduals(f$m, n = 1000, plot = FALSE)
  a <- as.data.frame(car::Anova(f$m, type = 3)); a$term <- rownames(a)
  a$model <- tag; a$family <- fam; a$re <- f$re
  a$dharma_unif_p <- testUniformity(sim, plot = FALSE)$p.value
  a$dharma_disp_p <- testDispersion(sim, plot = FALSE)$p.value
  a$dharma_zi_p   <- testZeroInflation(sim, plot = FALSE)$p.value
  em <- emmeans(f$m, specs, type = "response", offset = log(600))
  e <- as.data.frame(em); e$model <- tag
  p <- as.data.frame(pairs(em, adjust = "sidak")); p$model <- tag
  list(a = a, e = e, p = p)
}
r1 <- report("A_trial_x_interval_x_zone", "n ~ group * treatment * side + segment + offset(off)", d,
             ~ group * treatment | side)
tot <- d %>% group_by(trial, group, treatment, date, segment) %>%
  summarise(n = sum(n), off = log(600), .groups = "drop")
r2 <- report("B_trial_x_interval_total", "n ~ group * treatment + segment + offset(off)", tot,
             ~ group * treatment)
write.csv(bind_rows(r1$a, r2$a), file.path(OUT, "anova_refit.csv"), row.names = FALSE)
write.csv(bind_rows(r1$e, r2$e), file.path(OUT, "emmeans_refit.csv"), row.names = FALSE)
write.csv(bind_rows(r1$p, r2$p), file.path(OUT, "contrasts_refit.csv"), row.names = FALSE)

# Design-honest comparator: one row per trial (events per 20 min), log scale.
tr <- d %>% group_by(trial, group, treatment, date) %>% summarise(n = sum(n), .groups = "drop")
lm1 <- lm(log(n) ~ group * treatment, data = tr)
a3 <- as.data.frame(car::Anova(lm1, type = 3)); a3$term <- rownames(a3)
write.csv(a3, file.path(OUT, "anova_trial_level_lm.csv"), row.names = FALSE)
# zone-resolved comparator: trial x zone, trial as error stratum for between-trial terms
tz <- d %>% group_by(trial, group, treatment, side) %>% summarise(n = sum(n), .groups = "drop")
aovz <- summary(aov(log1p(n) ~ group * treatment * side + Error(trial), data = tz))
capture.output(aovz, file = file.path(OUT, "aov_trial_zone_splitplot.txt"))
cell <- tr %>% group_by(group, treatment) %>%
  summarise(trials = n(), mean = mean(n), se = sd(n) / sqrt(n()), min = min(n), max = max(n), .groups = "drop")
cellz <- tz %>% group_by(group, treatment, side) %>%
  summarise(mean = mean(n), se = sd(n) / sqrt(n()), .groups = "drop")
write.csv(cell, file.path(OUT, "raw_cell_means_events_per_trial.csv"), row.names = FALSE)
write.csv(cellz, file.path(OUT, "raw_cell_means_events_per_trial_zone.csv"), row.names = FALSE)
print(bind_rows(r1$a, r2$a), digits = 4); print(a3, digits = 4); print(aovz)
print(as.data.frame(cell)); print(as.data.frame(cellz))
print(r2$e, digits = 4); print(r2$p, digits = 4); print(r1$e, digits = 4); print(r1$p, digits = 4)
sink(file.path(OUT, "sessionInfo.txt")); print(sessionInfo()); sink()
