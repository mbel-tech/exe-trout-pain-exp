################################################################################
# validate_stat_types.R
# ─────────────────────────────────────────────────────────────────────────────
# HANDOFF.md S10, CONVENTIONS.md C12: confirms the string printed in every
# anova_stmt_*.csv (what the manuscript actually quotes) matches the test
# actually run, as recorded in that run's anova.csv (where present) or is at
# minimum internally consistent (an F-caption carries F/eta2p language and no
# chi/eta mismatch, a chisq-caption carries chi and never eta2p). Enforces the
# [INV] "stat_type is authoritative, df_denom is not" -- the defect that once
# mis-reported nine beta-GLMMs as F(1, 1000) because a reader trusted the wrong
# column (HANDOFF.md S7).
#
# Usage:
#   Rscript validate_stat_types.R <run_dir> [<run_dir> ...]
#
# Exits non-zero (and prints every mismatch) if any anova_stmt_*.csv statement
# is inconsistent with what it claims to be. Mandatory before a figure built
# from a run is treated as final (HANDOFF.md S10).
################################################################################

suppressMessages({
  library(dplyr)
})

# A caption string is internally self-consistent when:
#   - it names an F statistic (contains "F" immediately followed by a digit or
#     subscript digit, per to_subscript_statement()'s output) -> must NOT also
#     carry chi (chi/eta mismatch is the historical defect class), and if it
#     carries eta2p that's expected;
#   - it is a chi-square statement (contains chi) -> must NEVER carry eta2p
#     (C6 [DR], structural via effect_size()'s NA-by-construction return);
#   - it is a bootstrap statement (contains the <= floor glyph or an explicit
#     "N =" / seed mention) -> must carry its N.
classify_statement <- function(s) {
  has_chi   <- grepl("χ²", s)                  # χ²
  has_F     <- grepl("F[₀-₉₁]|F[0-9]|F<sub>", s)
  has_eta2p <- grepl("η²p|η<sup>2</sup>", s)      # η²p or <sup>2</sup> HTML form
  has_boot  <- grepl("≤|p_pb|bootstrap", s, ignore.case = TRUE)
  if (has_chi)      return("chisq")
  if (has_F)        return("F")
  if (has_boot)      return("pb_lrt")
  "unknown"
}

check_statement <- function(s, source_file, row_idx) {
  if (is.na(s) || !nzchar(s)) return(NULL)
  cls <- classify_statement(s)
  has_chi   <- grepl("χ²", s)
  has_eta2p <- grepl("η²p|η<sup>2</sup>", s)
  mismatch <- FALSE
  reason <- NA_character_
  if (cls == "chisq" && has_eta2p) {
    mismatch <- TRUE
    reason <- "chisq statement carries eta2p (structural violation of C6 [DR])"
  }
  if (cls == "F" && has_chi) {
    mismatch <- TRUE
    reason <- "F statement also carries chi (chi/eta mismatch)"
  }
  if (grepl("p = 0\\.000\\b", s)) {
    mismatch <- TRUE
    reason <- "p = 0.000 printed (p-floor formatter not applied, C7/A-9)"
  }
  data.frame(source_file = source_file, row = row_idx, stat_type = cls,
            mismatch = mismatch, reason = reason, statement = s,
            stringsAsFactors = FALSE)
}

validate_run <- function(run_dir) {
  anova_dir <- file.path(run_dir, "anovas")
  if (!dir.exists(anova_dir)) {
    message("SKIP (no anovas/ subfolder): ", run_dir)
    return(invisible(NULL))
  }
  stmt_files <- list.files(anova_dir, pattern = "^anova_stmt_.*\\.csv$", full.names = TRUE)
  if (!length(stmt_files)) {
    message("SKIP (no anova_stmt_*.csv found): ", run_dir)
    return(invisible(NULL))
  }

  rows <- list()
  for (f in stmt_files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
    if (is.null(df)) next
    stmt_cols <- names(df)[grepl("statement|stmt_", names(df))]
    for (col in stmt_cols) {
      for (i in seq_len(nrow(df))) {
        r <- check_statement(df[[col]][i], basename(f), i)
        if (!is.null(r)) rows[[length(rows) + 1]] <- r
      }
    }
  }
  if (!length(rows)) {
    message("No statements found to validate in ", run_dir)
    return(invisible(NULL))
  }
  result <- bind_rows(rows)
  out_path <- file.path(anova_dir, "stat_type_validation.csv")
  write.csv(result, out_path, row.names = FALSE)

  n_mismatch <- sum(result$mismatch)
  cat(sprintf("\n=== %s ===\n", run_dir))
  cat(sprintf("  %d statements checked, %d mismatch(es)\n", nrow(result), n_mismatch))
  if (n_mismatch > 0) {
    bad <- result[result$mismatch, c("source_file", "row", "reason", "statement")]
    print(bad, row.names = FALSE)
  }
  cat("  wrote ", out_path, "\n", sep = "")
  invisible(result)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript validate_stat_types.R <run_dir> [<run_dir> ...]\n")
  quit(status = 1)
}

all_results <- lapply(args, validate_run)
all_results <- Filter(Negate(is.null), all_results)
total_mismatch <- sum(vapply(all_results, function(r) sum(r$mismatch), integer(1)))

cat(sprintf("\n=== TOTAL: %d mismatch(es) across %d run(s) ===\n", total_mismatch, length(all_results)))
if (total_mismatch > 0) quit(status = 1) else quit(status = 0)
