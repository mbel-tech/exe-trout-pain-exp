# Tests for the shared random-effect helpers in house_style.R.
# Run: Rscript "05_statistics/tests/test_split_re_terms.R"
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
suppressMessages(source(file.path(here, "..", "shared", "house_style.R")))

# --- split_re_terms ---------------------------------------------------------
stopifnot(identical(
  split_re_terms("(1 | trial_date) + (1 | trial_id) + (1 | fish_id_matched_trial)"),
  c("(1 | trial_date)", "(1 | trial_id)", "(1 | fish_id_matched_trial)")))
stopifnot(identical(split_re_terms("(1 | a)"), "(1 | a)"))
stopifnot(identical(split_re_terms(NULL), character(0)))
stopifnot(identical(split_re_terms(NA_character_), character(0)))
stopifnot(identical(split_re_terms(""), character(0)))

# --- re_reduction_log_add / re_reduction_log_write ---------------------------
# Column order must match re_reduction.csv of run 20260826_201247.
want_cols <- c("label", "model_type", "re_attempted", "re_retained",
               "n_dropped", "re_dropped", "drop_reason")

d0 <- file.path(tempdir(), "red_empty"); dir.create(d0)
suppressMessages(re_reduction_log_write(d0))
e <- read.csv(file.path(d0, "re_reduction.csv"), stringsAsFactors = FALSE)
stopifnot(identical(names(e), want_cols), nrow(e) == 0L)    # empty = "nobody was singular"

re_reduction_log_add("ind/x", "Gaussian", "(1 | a) + (1 | b)", "(1 | a)",
                     "(1 | b)", "singular; SD(b) = 0")
re_reduction_log_add("ind/y", "Beta GLMM", "(1 | a)", "(1 | a)", character(0), character(0))
d1 <- file.path(tempdir(), "red_two"); dir.create(d1)
suppressMessages(re_reduction_log_write(d1))
f <- read.csv(file.path(d1, "re_reduction.csv"), stringsAsFactors = FALSE)
stopifnot(identical(names(f), want_cols), nrow(f) == 2L,
          f$n_dropped[1] == 1L, f$re_dropped[1] == "(1 | b)",
          f$n_dropped[2] == 0L, is.na(f$re_dropped[2]), is.na(f$drop_reason[2]),
          f$re_retained[2] == "(1 | a)")
# a model with every RE dropped logs retained = "none"
re_reduction_log_add("ind/z", "Gaussian", "(1 | a)", "", "(1 | a)", "singular; SD(a) = 0")
d2 <- file.path(tempdir(), "red_three"); dir.create(d2)
suppressMessages(re_reduction_log_write(d2))
g <- read.csv(file.path(d2, "re_reduction.csv"), stringsAsFactors = FALSE)
stopifnot(g$re_retained[3] == "none")

cat("PASS: test_split_re_terms.R\n")