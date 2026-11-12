# compare_runs.R -- compare two runs' anova.csv row by row.
# Usage: Rscript compare_runs.R <old_run_dir> <new_run_dir> <out_csv>
# A run dir may be the run directory itself or its anovas/ folder.
# Rows are keyed on label + term + stat_type (+ an occurrence counter, in case a
# key repeats). Output: one row per key with p_old, p_new, delta_p,
# crosses_alpha (opposite sides of 0.05), status (both / only_old / only_new).
# A row present in only one run is the failure to watch for: an indicator that
# silently vanishes from anova.csv is worse than a wrong number.
# Exits 0 always; prints the counts so the caller can decide.
args <- commandArgs(TRUE)
if (length(args) != 3) stop("usage: compare_runs.R <old_run_dir> <new_run_dir> <out_csv>")

find_anova <- function(d) {
  cand <- c(file.path(d, "anovas", "anova.csv"), file.path(d, "anova.csv"))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) hit <- list.files(d, pattern = "^anova\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(hit)) stop("no anova.csv under ", d)
  hit[1]
}
load_run <- function(d) {
  x <- read.csv(find_anova(d), stringsAsFactors = FALSE)
  keys <- intersect(c("label", "term", "stat_type"), names(x))
  if (!all(c("term", "p") %in% names(x))) stop("anova.csv in ", d, " lacks term/p columns")
  x$.key <- do.call(paste, c(x[keys], sep = "\r"))
  x$.occ <- ave(seq_len(nrow(x)), x$.key, FUN = seq_along)
  x$.id  <- paste(x$.key, x$.occ, sep = "\r#")
  x
}
old <- load_run(args[1]); new <- load_run(args[2])
pick <- function(x, nm) if (nm %in% names(x)) x[[nm]] else rep(NA, nrow(x))
ids <- union(old$.id, new$.id)
o <- old[match(ids, old$.id), ]; n <- new[match(ids, new$.id), ]
src <- ifelse(!is.na(o$.id), o$.id, n$.id)
parts <- strsplit(sub("\r#.*$", "", src), "\r")
res <- data.frame(
  label     = vapply(parts, function(p) if (length(p) >= 1) p[1] else NA_character_, ""),
  term      = vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, ""),
  stat_type = vapply(parts, function(p) if (length(p) >= 3) p[3] else NA_character_, ""),
  occ       = as.integer(sub(".*\r#", "", src)),
  p_old = pick(o, "p"), p_new = pick(n, "p"),
  statistic_old = pick(o, "statistic"), statistic_new = pick(n, "statistic"),
  model_type_old = pick(o, "model_type"), model_type_new = pick(n, "model_type"),
  stringsAsFactors = FALSE)
res$delta_p <- res$p_new - res$p_old
res$crosses_alpha <- !is.na(res$p_old) & !is.na(res$p_new) & ((res$p_old < 0.05) != (res$p_new < 0.05))
res$status <- ifelse(is.na(o$.id), "only_new", ifelse(is.na(n$.id), "only_old", "both"))
write.csv(res, args[3], row.names = FALSE)
cat(sprintf("rows: %d | both: %d | only_old: %d | only_new: %d | crosses_alpha: %d | model_type changed: %d | max|delta_p|: %s\n",
  nrow(res), sum(res$status == "both"), sum(res$status == "only_old"), sum(res$status == "only_new"),
  sum(res$crosses_alpha), sum(res$status == "both" & !is.na(res$model_type_old) & res$model_type_old != res$model_type_new),
  format(suppressWarnings(max(abs(res$delta_p), na.rm = TRUE)), digits = 3)))
invisible(NULL)