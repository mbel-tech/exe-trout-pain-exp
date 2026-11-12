# Single-transformation rule applied to the secondary Gaussian analyses that are not
# run by the two main pipelines (behaviour; plasma/monoamines):
#   (1) aggression events per trial (12 trials), two-way ANOVA
#   (2) aggression events per trial and zone (24), split-plot ANOVA (trial as error stratum)
#   (3) body mass of the sampled fish (48), two-way ANOVA
# Rule (as in the main pipelines): Shapiro-Wilk on the residuals of the fixed-effects
# model and Levene's test across the Treatment x Group cells on the raw response; if
# either p < 0.05, Box-Cox transform with the lambda (grid -2 to 2, step 0.1) that maximises
# the Shapiro-Wilk p, as in Bellio et al. (2026) (bcnPower with ML location gamma when the
# response contains zeros); otherwise untransformed. The same tests are repeated after transformation.
# Usage: Rscript 05_statistics/trial_level_anovas.R <aggression design-RE output dir> <sample dataset csv> [<out dir>]
#   (the first argument is the output folder of aggression_models_design_random_effects.R,
#    which contains grid_trial_interval_zone.csv)
suppressPackageStartupMessages({ library(car); library(dplyr) })
options(contrasts = c("contr.sum", "contr.poly"))
a <- commandArgs(TRUE); refit_dir <- a[1]; mono_csv <- a[2]
out <- if (length(a) >= 3) a[3] else file.path("output", "trial_level_anovas"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

bc_rule <- function(y, X, cells, label, unit) {
  fml <- as.formula(paste("y ~", X$rhs)); d <- X$data; d$y <- y
  sw  <- function(v) { d$y <- v; shapiro.test(residuals(lm(fml, d)))$p.value }
  lev <- function(v) leveneTest(v ~ cells)$`Pr(>F)`[1]
  sw0 <- sw(y); lv0 <- lev(y); yt <- y; tr <- "none"; lam <- NA; sw1 <- NA; lv1 <- NA
  if (sw0 < 0.05 || lv0 < 0.05) {
    m0 <- lm(fml, d)
    lam_grid <- seq(-2, 2, by = 0.1)   # lambda maximising the Shapiro-Wilk p (Bellio et al., 2026)
    best <- function(tf) { p <- vapply(lam_grid, function(l) { v <- tf(l); if (all(is.finite(v))) sw(v) else NA_real_ }, numeric(1))
                           round(lam_grid[which.max(p)], 1) }
    if (min(y) > 0) {
      lam <- best(function(l) bcPower(y, l))
      yt <- bcPower(y, lam); tr <- sprintf("Box-Cox (lambda = %.1f)", lam)
    } else {
      cf <- coef(powerTransform(m0, family = "bcnPower")); cfv <- if (is.matrix(cf)) cf[1, ] else cf
      g <- cfv[["gamma"]]; lam <- best(function(l) bcnPower(y, lambda = l, gamma = g))
      yt <- bcnPower(y, lambda = lam, gamma = g); tr <- sprintf("Box-Cox (lambda = %.1f, gamma = %.2f)", lam, g)
    }
    sw1 <- sw(yt); lv1 <- lev(yt)
  }
  list(y = yt, row = data.frame(indicator = label, unit = unit, n = length(y), sw_raw = sw0, levene_raw = lv0,
                                transformation = tr, lambda = lam, sw_trans = sw1, levene_trans = lv1))
}

g <- read.csv(file.path(refit_dir, "grid_trial_interval_zone.csv"))
tr <- g %>% group_by(trial, group, treatment) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(group = factor(group), treatment = factor(treatment))
tz <- g %>% group_by(trial, group, treatment, side) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(group = factor(group), treatment = factor(treatment), side = factor(side), trial = factor(trial))
s <- read.csv(mono_csv) %>% distinct(sample_id, .keep_all = TRUE) %>%
  filter(!is.na(weight_g)) %>% mutate(group = factor(group), treatment = factor(treatment))

rows <- list(); res <- list()
r1 <- bc_rule(tr$n, list(rhs = "group * treatment", data = tr), interaction(tr$group, tr$treatment),
              "Aggressive events (per trial)", "Trial (n = 12)")
tr$y <- r1$y; rows[[1]] <- r1$row
a1 <- as.data.frame(Anova(lm(y ~ group * treatment, tr), type = 3)); a1$term <- rownames(a1); a1$analysis <- "events_per_trial"
# split-plot: rule applied on the trial x zone totals (zone within trial)
r2 <- bc_rule(tz$n, list(rhs = "group * treatment * side", data = tz), interaction(tz$group, tz$treatment),
              "Aggressive events (per trial and zone)", "Trial x zone (n = 24)")
tz$y <- r2$y; rows[[2]] <- r2$row
a2 <- summary(aov(y ~ group * treatment * side + Error(trial), data = tz))
r3 <- bc_rule(s$weight_g, list(rhs = "group * treatment", data = s), interaction(s$group, s$treatment),
              "Body mass (g)", "Fish (n = 48)")
s$y <- r3$y; rows[[3]] <- r3$row
a3 <- as.data.frame(Anova(lm(y ~ group * treatment, s), type = 3)); a3$term <- rownames(a3); a3$analysis <- "body_mass"

write.csv(do.call(rbind, rows), file.path(out, "assumptions_secondary.csv"), row.names = FALSE)
write.csv(rbind(a1, a3), file.path(out, "anova_secondary.csv"), row.names = FALSE)
capture.output(a2, file = file.path(out, "aov_trial_zone_splitplot_boxcox.txt"))
print(do.call(rbind, rows)); print(rbind(a1, a3)); print(a2)
