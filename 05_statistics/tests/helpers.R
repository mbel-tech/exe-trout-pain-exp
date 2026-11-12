# Shared helpers for the tests in this folder.

# Evaluate ONE top-level `name <- function(...)` out of a script that cannot be
# sourced whole (sourcing it would run the entire analysis), and return it.
load_fn <- function(file, name, envir = globalenv()) {
  exprs <- parse(file, keep.source = FALSE, encoding = "UTF-8")
  for (e in exprs) {
    if (is.call(e) && length(e) == 3 &&
        (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
        identical(e[[2]], as.name(name))) {
      eval(e, envir)
      return(get(name, envir = envir))
    }
  }
  stop("function ", name, " not found in ", file)
}

# Folder this test file lives in (works under Rscript).
test_dir <- function() {
  args <- commandArgs(FALSE)
  dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
}

# The 2 x 2 GOT design at the level the endocrine engines see it: 12 trials of 4
# fish, 6 testing dates of 2 trials each (one Control, one Exercise), Group
# constant within a date (dates 1-3 = A, 4-6 = B).
#   trial_sd = 0 : y is additive in Treatment and Group plus a within-trial
#                  pattern that sums to zero, so the between-trial AND
#                  between-date variance are EXACTLY zero (a singular fit).
#   trial_sd > 0 : random trial offsets of that SD plus small noise.
make_design_data <- function(trial_sd = 0, seed = 1) {
  set.seed(seed)
  d <- data.frame(trial_uid = rep(1:12, each = 4))
  d$date_uid  <- ceiling(d$trial_uid / 2)
  d$Group     <- factor(ifelse(d$date_uid <= 3, "A", "B"))
  d$Treatment <- factor(ifelse(d$trial_uid %% 2 == 1, "Control", "Exercise"))
  tr <- d$trial_uid
  d$trial_uid <- factor(d$trial_uid); d$date_uid <- factor(d$date_uid)
  base <- 10 + 0.5 * (d$Treatment == "Exercise") + 1 * (d$Group == "B")
  pat  <- rep(c(-1.5, -0.5, 0.5, 1.5), 12)
  if (trial_sd > 0) {
    d$y <- base + rnorm(12, 0, trial_sd)[tr] + pat + rnorm(48, 0, 0.3)
  } else {
    d$y <- base + pat
  }
  d
}

# Strictly positive version of the same design, for the Gamma / Tweedie selectors.
#   trial_sd = 0 : every trial in a Treatment x Group cell carries IDENTICAL data
#                  (same fish pattern), so between-trial and between-date
#                  variance are exactly zero while the fish still differ.
#   trial_sd > 0 : random trial offsets (log scale) of that SD plus small noise.
make_positive_design_data <- function(trial_sd = 0, seed = 1) {
  d <- make_design_data(trial_sd = 0, seed = seed)
  set.seed(seed)
  tr <- as.integer(d$trial_uid)
  eta <- log(10) + 0.05 * (d$Treatment == "Exercise") + 0.10 * (d$Group == "B") +
         rep(c(-0.3, -0.1, 0.15, 0.25), 12)
  if (trial_sd > 0) eta <- eta + rnorm(12, 0, trial_sd)[tr] + rnorm(48, 0, 0.05)
  d$y <- exp(eta)
  d
}