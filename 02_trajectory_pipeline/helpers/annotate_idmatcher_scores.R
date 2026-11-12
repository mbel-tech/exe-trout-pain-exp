# =============================================================================
# annotate_idmatcher_scores.R
#
# Reads idmatcher_identity_correspondences.csv, labels every matched
# segment-pair row with a match_quality / match_accepted flag, and applies
# tiered ambiguity-resolution rules so that as many rows as possible are
# salvaged without compromising correctness.
#
# Quality labels (highest confidence first):
#   "strong"        – both direct & indirect scores >= 0.5
#   "acceptable"    – both scores >= MIN_MATCH_SCORE, at least one < 0.5
#   "by_exclusion"  – score(s) below threshold, but all other fish in the
#                     same (trial × segment) group are accepted, so this
#                     identity is the only remaining possibility (bijective
#                     assignment guaranteed by the Hungarian algorithm)
#   "by_consensus"  – multiple-rejection group, BOTH the body-length tiebreaker
#                     AND the matrix tiebreaker confirm idmatcher's assignment
#   "by_matrix"     – matrix tiebreaker confirms idmatcher's assignment
#                     (body-length unavailable, unstable, or unused)
#   "by_bodylength" – body-length tiebreaker confirms idmatcher's assignment
#                     (matrix unavailable or unused)
#   "rejected_conflict" – tiebreaker(s) DISAGREE with idmatcher's assignment;
#                         post-processing in 00_idmatcher_helpers.R may resolve
#                         these when reliable anchors exist in the same trial
#   "rejected"      – multi-rejection group with no tiebreaker available
#   NA              – master segment row; no scores are produced for it
#
# Pipeline use:
#   Called via annotate_idmatcher_scores_run() from step1_5_annotate_idmatcher_GOT.R.
#   Can also be sourced interactively (standalone mode uses hardcoded paths below).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(clue)        # solve_LSAP (Hungarian)
})

# ── Default configuration (overridable via function parameters) ───────────────
MIN_MATCH_SCORE       <- 0.0
USE_BL_TIEBREAKER     <- TRUE
USE_MATRIX_TIEBREAKER <- TRUE
# Body-length tiebreaker decision margin (cm, per fish).
# If the Hungarian-optimal assignment is cheaper than idmatcher's choice
# by less than BL_MARGIN_CM × k, the BL tiebreaker abstains (returns NA).
BL_MARGIN_CM <- 2.0

# ── Standalone defaults (used only in interactive mode at bottom) ─────────────
PROJECT_ROOT <- "."
INPUT_PATH   <- file.path(PROJECT_ROOT, "datasets", "idmatcher_identity_correspondences.csv")
OUTPUT_PATH  <- file.path(PROJECT_ROOT, "datasets", "idmatcher_annotated.csv")
BL_CSV_PATH  <- file.path(PROJECT_ROOT, "Checked sessions GOT", "master_outputs",
                          "morph", "master_bodylength_by_fish.csv")
MATRICES_DIR <- file.path(PROJECT_ROOT, "datasets", "matrices")

# ── Tiebreaker helpers (module-level; used inside annotate_idmatcher_scores_run) ──

bodylength_agrees <- function(rejected_subset, bl_data, margin_cm = BL_MARGIN_CM) {
  trial_   <- rejected_subset$trial[1]
  master_s <- rejected_subset$master_segment[1]
  match_s  <- rejected_subset$segment[1]

  master_bl <- bl_data |>
    filter(trial_id == trial_, segment_id == master_s,
           fish_in_session %in% rejected_subset$fish_id_master) |>
    arrange(match(fish_in_session, rejected_subset$fish_id_master))
  matching_bl <- bl_data |>
    filter(trial_id == trial_, segment_id == match_s,
           fish_in_session %in% rejected_subset$fish_id) |>
    arrange(match(fish_in_session, rejected_subset$fish_id))

  if (nrow(master_bl) != nrow(rejected_subset) ||
      nrow(matching_bl) != nrow(rejected_subset)) return(NA)
  if (any(is.na(master_bl$body_length_cm)) ||
      any(is.na(matching_bl$body_length_cm))) return(NA)

  k <- nrow(rejected_subset)
  if (k < 2) return(NA)

  cost <- outer(master_bl$body_length_cm, matching_bl$body_length_cm,
                FUN = function(a, b) abs(a - b))
  hungarian      <- as.integer(solve_LSAP(cost, maximum = FALSE))
  optimal_cost   <- sum(cost[cbind(seq_len(k), hungarian)])
  idmatcher_cost <- sum(diag(cost))
  margin         <- idmatcher_cost - optimal_cost

  if (identical(hungarian, seq_len(k))) return(TRUE)
  if (margin < margin_cm * k) return(NA)
  FALSE
}

matrix_agrees <- function(rejected_subset, matrices_dir) {
  trial_   <- rejected_subset$trial[1]
  master_s <- rejected_subset$master_segment[1]
  match_s  <- rejected_subset$segment[1]

  stem <- sprintf("trial%02d_seg%d_vs_seg%d", trial_, match_s, master_s)
  direct_files <- list.files(matrices_dir,
    pattern = paste0("^", stem, ".*direct.*\\.csv$"), full.names = TRUE)
  indirect_files <- list.files(matrices_dir,
    pattern = paste0("^", stem, ".*indirect.*\\.csv$"), full.names = TRUE)
  if (length(direct_files) == 0 || length(indirect_files) == 0) return(NA)

  read_matrix <- function(path) {
    m <- read_csv(path, col_names = FALSE, show_col_types = FALSE, progress = FALSE)
    if (!is.numeric(m[[1]])) m <- m[, -1]
    as.matrix(m)
  }
  D <- read_matrix(direct_files[1])
  I <- read_matrix(indirect_files[1])
  if (any(dim(D) != dim(I))) return(NA)
  joined <- t(D) + I
  rows <- rejected_subset$fish_id_master
  cols <- rejected_subset$fish_id
  if (any(rows > nrow(joined)) || any(cols > ncol(joined))) return(NA)
  sub <- joined[rows, cols, drop = FALSE]
  if (nrow(sub) < 2) return(NA)
  hungarian <- as.integer(solve_LSAP(sub, maximum = TRUE))
  identical(hungarian, seq_len(nrow(rejected_subset)))
}


# =============================================================================
# annotate_idmatcher_scores_run()
#
# Main engine: reads the wide idmatcher correspondences CSV, pivots to long
# format, and labels each row with a match_quality / match_accepted flag.
#
# Parameters:
#   input_path    Path to idmatcher_identity_correspondences.csv
#   bl_csv_path   Path to master_bodylength_by_fish.csv (NULL to skip BL tiebreaker)
#   matrices_dir  Path to the matrices/ folder from the Colab run (NULL to skip)
#   output_path   If non-NULL, write the annotated long table to this path
#   min_match_score, use_bl_tiebreaker, use_matrix_tiebreaker, bl_margin_cm
#                 Override the module-level defaults for this call only
#
# Returns: tibble (annotated long table) invisibly
# =============================================================================
annotate_idmatcher_scores_run <- function(
    input_path,
    bl_csv_path        = NULL,
    matrices_dir       = NULL,
    output_path        = NULL,
    min_match_score    = MIN_MATCH_SCORE,
    use_bl_tiebreaker  = USE_BL_TIEBREAKER,
    use_matrix_tiebreaker = USE_MATRIX_TIEBREAKER,
    bl_margin_cm       = BL_MARGIN_CM
) {

  # ── Schema validation ──────────────────────────────────────────────────────
  raw <- read_csv(input_path, show_col_types = FALSE)
  required_cols <- c("trial", "fish_id_master", "master_segment", "master_session")
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols))
    stop("idmatcher CSV missing required columns: ", paste(missing_cols, collapse = ", "))
  if (!any(grepl("^seg\\d+_", names(raw))))
    stop("idmatcher CSV has no seg<N>_* columns — check the Colab output format.")

  # ── Wide → long ────────────────────────────────────────────────────────────
  long <- raw |>
    pivot_longer(
      cols          = matches("^seg\\d+_"),
      names_to      = c("segment", ".value"),
      names_pattern = "seg(\\d+)_(.*)"
    ) |>
    mutate(segment = as.integer(segment))

  annotated <- long |>
    mutate(
      match_quality = case_when(
        is.na(direct_score) | is.na(indirect_score) ~ NA_character_,
        direct_score >= 0.5 & indirect_score >= 0.5  ~ "strong",
        direct_score >= min_match_score & indirect_score >= min_match_score ~ "acceptable",
        TRUE ~ "rejected"
      ),
      match_accepted = case_when(
        is.na(match_quality)        ~ NA,
        match_quality == "rejected" ~ FALSE,
        TRUE                        ~ TRUE
      )
    ) |>
    filter(segment != master_segment)

  # ── Rule 1: single-rejection exclusion ────────────────────────────────────
  annotated <- annotated |>
    group_by(trial, segment) |>
    mutate(
      n_rej = sum(match_accepted == FALSE, na.rm = TRUE),
      n_acc = sum(match_accepted == TRUE,  na.rm = TRUE),
      match_quality = if_else(
        match_quality == "rejected" & n_rej == 1 & n_acc == n() - 1,
        "by_exclusion", match_quality
      ),
      match_accepted = if_else(match_quality == "by_exclusion", TRUE, match_accepted)
    ) |>
    select(-n_rej, -n_acc) |>
    ungroup()

  # ── Load optional reference data ───────────────────────────────────────────
  bl_data        <- NULL
  matrices_avail <- FALSE
  if (use_bl_tiebreaker && !is.null(bl_csv_path) && file.exists(bl_csv_path)) {
    bl_data <- read_csv(bl_csv_path, show_col_types = FALSE) |>
      select(trial_id, segment_id, fish_in_session,
             body_length_cm, body_length_unstable)
  }
  if (use_matrix_tiebreaker && !is.null(matrices_dir) && dir.exists(matrices_dir) &&
      length(list.files(matrices_dir, pattern = "\\.csv$")) > 0) {
    matrices_avail <- TRUE
  }

  # ── Apply tiebreakers to multi-rejection groups ────────────────────────────
  multi_rej_keys <- annotated |>
    group_by(trial, segment) |>
    summarise(n_rej = sum(match_quality == "rejected", na.rm = TRUE),
              .groups = "drop") |>
    filter(n_rej >= 2)

  if (nrow(multi_rej_keys) > 0) {
    cat(sprintf("\nResolving %d multi-rejection group(s) with tiebreakers...\n",
                nrow(multi_rej_keys)))
    for (i in seq_len(nrow(multi_rej_keys))) {
      t_i <- multi_rej_keys$trial[i]
      s_i <- multi_rej_keys$segment[i]
      subset_rows <- annotated$trial == t_i & annotated$segment == s_i &
                     annotated$match_quality == "rejected"
      rejected_subset <- annotated[subset_rows, , drop = FALSE]

      bl_ok  <- if (!is.null(bl_data))  bodylength_agrees(rejected_subset, bl_data, bl_margin_cm) else NA
      mat_ok <- if (matrices_avail)     matrix_agrees(rejected_subset, matrices_dir) else NA

      new_label <- case_when(
        isTRUE(bl_ok)  & isTRUE(mat_ok)  ~ "by_consensus",
        isTRUE(mat_ok) & is.na(bl_ok)    ~ "by_matrix",
        isTRUE(bl_ok)  & is.na(mat_ok)   ~ "by_bodylength",
        isTRUE(bl_ok)  & isFALSE(mat_ok) ~ "rejected_conflict",
        isFALSE(bl_ok) & isTRUE(mat_ok)  ~ "rejected_conflict",
        isFALSE(bl_ok) | isFALSE(mat_ok) ~ "rejected_conflict",
        TRUE                              ~ "rejected"
      )
      annotated$match_quality[subset_rows]  <- new_label
      annotated$match_accepted[subset_rows] <- new_label %in%
        c("by_consensus", "by_matrix", "by_bodylength")

      cat(sprintf("  trial %2d seg %d  (k=%d)  BL=%s  matrix=%s  -> %s\n",
                  t_i, s_i, nrow(rejected_subset),
                  ifelse(is.na(bl_ok), "n/a",  ifelse(bl_ok, "agree", "disagree")),
                  ifelse(is.na(mat_ok), "n/a", ifelse(mat_ok, "agree", "disagree")),
                  new_label))
    }
  }

  # ── Summary ────────────────────────────────────────────────────────────────
  cat("\n── idmatcher score annotation ─────────────────────────────────────\n")
  cat(sprintf("Threshold (min_match_score) : %.2f\n", min_match_score))
  cat(sprintf("Body-length tiebreaker      : %s\n",
              if (!is.null(bl_data)) "ACTIVE" else "skipped (file missing or disabled)"))
  cat(sprintf("Matrix tiebreaker           : %s\n",
              if (matrices_avail) "ACTIVE" else
                "skipped (matrices/ missing, empty, or disabled)"))
  cat(sprintf("Input rows (wide)           : %d\n", nrow(raw)))
  cat(sprintf("Output rows (long, matched) : %d\n", nrow(annotated)))
  cat("\nMatch quality breakdown:\n")
  print(
    annotated |>
      count(segment, match_quality) |>
      mutate(match_quality = replace_na(match_quality, "(master/NA)")) |>
      arrange(segment, match_quality)
  )
  cat(sprintf("\nAccepted matches : %d / %d\n",
              sum(annotated$match_accepted, na.rm = TRUE),
              sum(!is.na(annotated$match_accepted))))
  cat("───────────────────────────────────────────────────────────────────\n")

  # ── Write (optional) ──────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    write_csv(annotated, output_path)
    cat(sprintf("\nSaved: %s\n", output_path))
  }

  invisible(annotated)
}


# =============================================================================
# Standalone / interactive mode
# Runs the annotation with hardcoded defaults when this script is sourced
# directly (not called from the pipeline orchestrator).
# =============================================================================
if (!exists("ANNOTATE_IDMATCHER_CALLED_FROM_PIPELINE", inherits = FALSE)) {
  annotate_idmatcher_scores_run(
    input_path   = INPUT_PATH,
    bl_csv_path  = BL_CSV_PATH,
    matrices_dir = MATRICES_DIR,
    output_path  = OUTPUT_PATH
  )
}
