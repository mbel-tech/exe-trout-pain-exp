# =============================================================================
# STEP 1.6 — HUNGARIAN BODY-LENGTH ID MATCHING (R5: promoted from physiology/)
#
# Run by 00_run_trajectory_pipeline.R between STEP 1.5 and STEP 2.
# Reads PIPELINE_DATA_DIR and SAMPLING_DATA_XLSX from the master config.
# Clipboard input has been removed — set PIPELINE_DATA_DIR in master instead.
#
# What this script does:
# 1) Reads PIPELINE_DATA_DIR -> master_output (no clipboard input)
# 2) Recursively finds all CSV files named "master_bodylength_by_fish.csv"
# 3) Imports and combines them into master_bodylength_by_fish
# 4) Keeps only rows whose fish_id is exactly of form TX_SY_FZ and has S3
#    e.g. T1_S3_F1, T12_S3_F10
# 5) Imports Excel and fills "trial" for SAMPLE_ID >= 31
# 6) Matches sampled fish IDs to tracked fish IDs by body length
#    within each trial x segment_id using Hungarian assignment
# 7) Creates diagnostics and QC flags
# 8) Exports match_table_export_<timestamp>.xlsx in the same directory as morph
#    AND a companion bl_reliability_tier_<timestamp>.csv for downstream use
# 9) (R4) Also copies the BL-tier CSV to a stable canonical path:
#    PIPELINE_DATA_DIR/master_outputs/reliability/reliability_tier_LATEST.csv
#    so that STEP 5/6/side-workflows can fall back to that path when
#    BL_MATCH_PATH is not set in .GlobalEnv.
#
# Threshold strategy used here:
#   unreliable  : relative_margin < 0.4
#   borderline  : 0.4 <= relative_margin < 0.8
#   reliable    : relative_margin >= 0.8
#   conservative: reliable == YES AND abs_diff_cm <= 0.5
#
# Downstream use (monoamine / plasma / aggression scripts):
#   Set BL_MATCH_PATH to the path of "bl_reliability_tier_<timestamp>.csv".
#   That CSV contains columns:
#     trial_id        <integer>   — trial number
#     reliability_tier <character> — one of "unreliable", "borderline",
#                                    "reliable", "conservative"
#   qualifying_trial_ids() in 00_stats_helpers.R reads those two columns to
#   determine which trials pass the MIN_BL_RELIABILITY_TIER gate before adding
#   a (1|Sample_ID) or (1|fish_id) random-effect candidate.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(janitor)
  library(stringr)
  library(clue)      # solve_LSAP
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(writexl)
  library(tibble)
})

`%||%` <- function(x, y) {
  if (length(x) == 0 || all(is.na(x))) y else x
}

# -----------------------------------------------------------------------------
# 0) Resolve master_output from PIPELINE_DATA_DIR (R5: clipboard removed)
# -----------------------------------------------------------------------------
# Master sets PIPELINE_DATA_DIR before sourcing this file. The previous
# clipboard fallback for standalone runs has been removed; if a user wants to
# run this script outside the pipeline they must set PIPELINE_DATA_DIR first.
if (!exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
  stop(
    "PIPELINE_DATA_DIR is not set. Source 00_run_trajectory_pipeline.R ",
    "(which sets PIPELINE_DATA_DIR) before running this script, or assign ",
    "PIPELINE_DATA_DIR <- \"<path to session folder>\" manually.",
    call. = FALSE
  )
}
master_output <- normalizePath(get("PIPELINE_DATA_DIR", envir = .GlobalEnv),
                               winslash = "/", mustWork = FALSE)
assign("master_output", master_output, envir = .GlobalEnv)
message("master_output = ", master_output)

# -----------------------------------------------------------------------------
# 1) Find + import all "master_bodylength_by_fish.csv" under master_output
# -----------------------------------------------------------------------------
target_name <- "master_bodylength_by_fish.csv"

csv_paths <- list.files(
  path       = master_output,
  pattern    = "^master_bodylength_by_fish\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(csv_paths) == 0) {
  # -----------------------------------------------------------------------------
  # 1-fallback (added 2026-08-19): master_bodylength_by_fish.csv is not produced
  # anywhere in this tree (the intended morph/ output folder is empty — STEP 1
  # has not been re-run since May). The same per-fish lengths DO survive in
  # STEP2's bodylength_by_fish.csv (fish_id, body_length_cm, trial_id, timepoint),
  # so fall back to the newest STEP2_output run rather than requiring a fresh
  # STEP1 pass over ~13 GB of raw idtracker.ai sessions. This branch only
  # substitutes the INPUT source — everything below (Hungarian assignment,
  # reliability thresholds, export format) is unchanged. If a real
  # master_bodylength_by_fish.csv ever appears, the search above finds it first
  # and this fallback is never reached. See LOGS.md 2026-08-19.
  # Prefer the location relative to this script's own folder (.pipeline_dir is
  # set near the top of 00_run_trajectory_pipeline.R and sourced before
  # this file runs); fall back to a hardcoded absolute path, then cwd, so this
  # still works if the script is ever run standalone with only PIPELINE_DATA_DIR set.
  step2_root <- if (exists(".pipeline_dir", envir = .GlobalEnv, inherits = FALSE)) {
    file.path(get(".pipeline_dir", envir = .GlobalEnv), "STEP2_output")
  } else {
    file.path(get0(".repo_root", ifnotfound = ".."), "02_trajectory_pipeline/STEP2_output")
  }
  if (!dir.exists(step2_root)) {
    step2_root <- file.path(getwd(), "STEP2_output")
  }
  step2_runs <- list.dirs(step2_root, recursive = FALSE, full.names = TRUE)
  step2_runs <- step2_runs[grepl("STEP2_output_[0-9]{8}_[0-9]{6}$", step2_runs)]
  if (length(step2_runs) == 0) {
    stop(
      "No file named '", target_name, "' found under: ", master_output,
      ", and no STEP2_output_* run found under '", step2_root,
      "' to fall back to. Re-run STEP1 to produce the morph file, or point ",
      "PIPELINE_DATA_DIR somewhere that has one."
    )
  }
  step2_latest <- sort(step2_runs)[length(step2_runs)]
  fallback_csv <- file.path(step2_latest, "bodylength_by_fish.csv")
  if (!file.exists(fallback_csv)) {
    stop(
      "No file named '", target_name, "' found under: ", master_output,
      ", and the latest STEP2_output run ('", step2_latest,
      "') has no bodylength_by_fish.csv either."
    )
  }
  message(
    "No '", target_name, "' found under ", master_output, ".\n",
    "Falling back to STEP2 tracked lengths: ", fallback_csv
  )
  master_bodylength_by_fish <- readr::read_csv(fallback_csv, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      source_file = fallback_csv,
      # STEP2 stores the segment as `timepoint` (verified 1:1 with the fish_id's
      # own S<n> substring, e.g. fish_id T15_S2_F1 <-> timepoint 2); the matcher
      # downstream (see master_obs, ~line 250) needs a `segment_id` column.
      segment_id  = as.character(timepoint),
      .before = 1
    )
  csv_paths <- character(0)  # deliberately empty: keeps the later output-location
                              # logic (~line 1155) falling through to master_output,
                              # not into the STEP2_output run folder.
  message(
    "Imported master_bodylength_by_fish from STEP2 fallback: ",
    nrow(master_bodylength_by_fish), " rows, ",
    ncol(master_bodylength_by_fish), " cols."
  )
} else {

message("Found ", length(csv_paths), " file(s):")
message(paste0(" - ", csv_paths, collapse = "\n"))

expected_cols <- c(
  "fish_id", "fish_in_session", "n_bbox_frames", "px_median", "px_mean_trim",
  "body_length_sd_px", "body_length_method", "body_length_px",
  "body_length_heavy_tail", "body_length_unstable",
  "body_length_unstable_reason", "body_length_cm", "body_length_sd_cm",
  "session_id", "trial_id", "segment_id"
)

read_one_bodylength <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  
  df <- df %>%
    mutate(source_file = path, .before = 1)
  
  missing_cols <- setdiff(expected_cols, names(df))
  if (length(missing_cols) > 0) {
    warning(
      "File missing expected columns: ", basename(path), " -> ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df
}

master_bodylength_by_fish <- dplyr::bind_rows(lapply(csv_paths, read_one_bodylength))

assign("master_bodylength_by_fish", master_bodylength_by_fish, envir = .GlobalEnv)
message(
  "Imported master_bodylength_by_fish: ",
  nrow(master_bodylength_by_fish), " rows, ",
  ncol(master_bodylength_by_fish), " cols."
)

}  # end of "found master_bodylength_by_fish.csv normally" branch (1-fallback above)

assign("master_bodylength_by_fish", master_bodylength_by_fish, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# 1b) Keep only rows whose fish_id is exactly of form TX_SY_FZ and has S3
# -----------------------------------------------------------------------------
master_bodylength_by_fish <- master_bodylength_by_fish %>%
  mutate(fish_id = as.character(fish_id)) %>%
  filter(
    !is.na(fish_id),
    grepl("^T\\d+_S3_F\\d+$", fish_id)
  )

assign("master_bodylength_by_fish", master_bodylength_by_fish, envir = .GlobalEnv)
message(
  "Filtered master_bodylength_by_fish to fish_id pattern ^T\\d+_S3_F\\d+$ : ",
  nrow(master_bodylength_by_fish), " rows remain."
)

# =============================================================================
# EXCEL PART
# =============================================================================

# -----------------------------------------------------------------------------
# 2) Import Excel
# -----------------------------------------------------------------------------
xlsx_path <- if (exists("SAMPLING_DATA_XLSX", envir = .GlobalEnv, inherits = FALSE)) {
  get("SAMPLING_DATA_XLSX", envir = .GlobalEnv)
} else {
  file.path(get0(".repo_root", ifnotfound = ".."), "data/sampling data GOT.xlsx")
}

df <- readxl::read_excel(xlsx_path) %>%
  janitor::clean_names()

# -----------------------------------------------------------------------------
# Remove individuals marked as "dead in arena, not euthanized"
# -----------------------------------------------------------------------------
if ("notes" %in% names(df)) {
  
  dead_ids <- df %>%
    filter(
      !is.na(notes),
      str_detect(
        str_to_lower(notes),
        "dead in arena, not euthanized"
      )
    ) %>%
    pull(sample_id)
  
  if (length(dead_ids) > 0) {
    
    message(
      "\nRemoving SAMPLE_ID(s) marked 'dead in arena, not euthanized': ",
      paste(dead_ids, collapse = ", ")
    )
    
    df <- df %>%
      filter(!(sample_id %in% dead_ids))
  }
}

trial_counts <- df %>%
  count(trial)

if (any(trial_counts$n < 4)) {
  warning(
    "Some trials now have <4 fish after removing dead individuals:\n",
    paste(
      trial_counts$trial[trial_counts$n < 4],
      collapse = ", "
    )
  )
}

# -----------------------------------------------------------------------------
# 3) Fill trial from SAMPLE_ID >= 31 onward
#    - each trial has 4 fish
#    - start at trial 5
# -----------------------------------------------------------------------------
start_sample_id <- 31L
fish_per_trial  <- 4L
start_trial     <- 5L

if (!("sample_id" %in% names(df))) {
  stop("Column 'SAMPLE_ID' not found (expected as 'sample_id' after clean_names()).")
}
if (!("trial" %in% names(df))) {
  stop("Column 'trial' not found (expected as 'trial' after clean_names()).")
}
if (!("lenght" %in% names(df))) {
  stop("Column 'lenght' not found after clean_names().")
}

df <- df %>%
  mutate(
    sample_id_num = suppressWarnings(as.integer(sample_id)),
    trial = if_else(
      !is.na(sample_id_num) & sample_id_num >= start_sample_id,
      start_trial + ((sample_id_num - start_sample_id) %/% fish_per_trial),
      trial
    )
  ) %>%
  select(-sample_id_num)

assign("df", df, envir = .GlobalEnv)
message("Excel imported as df: ", nrow(df), " rows, ", ncol(df), " cols.")

# =============================================================================
# LENGTH-BASED IDENTITY MATCHING WORKFLOW
# =============================================================================

# ---- Prep "true" lengths from df (Excel) ----
df_true <- df %>%
  mutate(
    trial = suppressWarnings(as.integer(trial)),
    df_lenght_cm = suppressWarnings(as.numeric(lenght)),
    df_sample_id = as.character(sample_id)
  ) %>%
  filter(!is.na(trial), !is.na(df_lenght_cm)) %>%
  select(trial, df_sample_id, df_lenght_cm)

# ---- Prep observed lengths from master_bodylength_by_fish ----
master_obs <- master_bodylength_by_fish %>%
  transmute(
    trial          = suppressWarnings(as.integer(trial_id)),
    segment_id     = as.character(segment_id),
    master_fish_id = as.character(fish_id),
    bl_cm          = suppressWarnings(as.numeric(body_length_cm))
  ) %>%
  filter(!is.na(trial), !is.na(segment_id), !is.na(bl_cm))

# ---- Substitute trial-mean BL for matched fish when idmatcher ran ----
# idmatcher_id_map is written to GlobalEnv by STEP 1.5.  For fish whose
# seg2<->seg3 identity was accepted, body_length_mean_cm is the mean BL across
# both segments — a better estimate than the per-segment value alone.
# Unmatched fish (NA body_length_mean_cm) fall back to their per-segment bl_cm.
if (exists("idmatcher_id_map", envir = .GlobalEnv, inherits = FALSE)) {
  .idmap_bl <- get("idmatcher_id_map", envir = .GlobalEnv) %>%
    dplyr::select(tracker_fish_id, body_length_mean_cm) %>%
    dplyr::distinct()
  master_obs <- master_obs %>%
    dplyr::left_join(.idmap_bl, by = c("master_fish_id" = "tracker_fish_id")) %>%
    dplyr::mutate(bl_cm = dplyr::coalesce(body_length_mean_cm, bl_cm)) %>%
    dplyr::select(-body_length_mean_cm)
  rm(.idmap_bl)
}

# -----------------------------------------------------------------------------
# Helper: match ONE group (.x = one trial x segment_id) to true_df for that trial
# -----------------------------------------------------------------------------
match_one_trial_segment <- function(obs_df, true_df) {
  obs_df  <- obs_df  %>% arrange(master_fish_id)
  true_df <- true_df %>% arrange(df_sample_id)
  
  n_obs  <- nrow(obs_df)
  n_true <- nrow(true_df)
  
  if (n_obs == 0 || n_true == 0) {
    return(tibble(
      df_sample_id = character(),
      master_fish_id = character(),
      df_lenght_cm = numeric(),
      bl_cm = numeric(),
      abs_diff_cm = numeric(),
      ambiguity_margin_cm = numeric()
    ))
  }
  
  cost0 <- abs(outer(obs_df$bl_cm, true_df$df_lenght_cm, FUN = "-"))
  
  n <- max(n_obs, n_true)
  cost <- matrix(1e6, nrow = n, ncol = n)
  cost[seq_len(n_obs), seq_len(n_true)] <- cost0
  
  assignment <- solve_LSAP(cost)
  assigned_true_idx <- as.integer(assignment)[seq_len(n_obs)]
  assigned_true_idx[assigned_true_idx > n_true] <- NA_integer_
  
  best <- apply(cost0, 1, min, na.rm = TRUE)
  second_best <- apply(cost0, 1, function(v) {
    v <- sort(v)
    if (length(v) >= 2) v[2] else NA_real_
  })
  ambiguity_margin <- second_best - best
  
  obs_df %>%
    transmute(
      df_sample_id = true_df$df_sample_id[assigned_true_idx],
      master_fish_id = master_fish_id,
      df_lenght_cm = true_df$df_lenght_cm[assigned_true_idx],
      bl_cm = bl_cm,
      abs_diff_cm = abs(bl_cm - true_df$df_lenght_cm[assigned_true_idx]),
      ambiguity_margin_cm = ambiguity_margin
    )
}

# -----------------------------------------------------------------------------
# Run per trial x segment_id
# -----------------------------------------------------------------------------
match_table <- master_obs %>%
  group_by(trial, segment_id) %>%
  group_modify(function(.x, .y) {
    true_df_this_trial <- df_true %>% filter(trial == .y$trial[[1]])
    match_one_trial_segment(.x, true_df_this_trial)
  }) %>%
  ungroup() %>%
  arrange(trial, segment_id, master_fish_id) %>%
  select(
    trial, segment_id, df_sample_id, master_fish_id,
    df_lenght_cm, bl_cm, abs_diff_cm, ambiguity_margin_cm
  )

assign("match_table", match_table, envir = .GlobalEnv)

print(match_table)
message("\nSaved in Global Env: match_table")

# =============================================================================
# MATCH TABLE DIAGNOSTICS + QC FLAGS + EXCEL EXPORT
# =============================================================================

if (!exists("match_table", envir = .GlobalEnv)) {
  stop("match_table not found in Global Environment.")
}
if (!exists("df", envir = .GlobalEnv)) {
  stop("df (Excel table) not found in Global Environment.")
}

mt <- get("match_table", envir = .GlobalEnv)

req_cols <- c("trial", "segment_id", "abs_diff_cm", "ambiguity_margin_cm")
missing_req <- setdiff(req_cols, names(mt))
if (length(missing_req) > 0) {
  stop("match_table is missing required columns: ", paste(missing_req, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 4) Trial spacing + relative margin
# -----------------------------------------------------------------------------
df_true <- df %>%
  mutate(
    trial = suppressWarnings(as.integer(trial)),
    df_sample_id = as.character(sample_id),
    df_lenght_cm = suppressWarnings(as.numeric(lenght))
  ) %>%
  filter(!is.na(trial), !is.na(df_lenght_cm)) %>%
  select(trial, df_sample_id, df_lenght_cm)

trial_spacing <- df_true %>%
  group_by(trial) %>%
  summarise(
    n_fish = n(),
    min_spacing_cm = {
      x <- sort(df_lenght_cm)
      if (length(x) >= 2) min(diff(x)) else NA_real_
    },
    .groups = "drop"
  )

mt <- mt %>%
  left_join(trial_spacing %>% select(trial, min_spacing_cm), by = "trial") %>%
  mutate(
    relative_margin = dplyr::if_else(
      !is.na(min_spacing_cm) & is.finite(min_spacing_cm) & min_spacing_cm > 0 &
        !is.na(ambiguity_margin_cm) & is.finite(ambiguity_margin_cm),
      ambiguity_margin_cm / min_spacing_cm,
      NA_real_
    )
  )

# -----------------------------------------------------------------------------
# 5) Helpers for axis breaks
# -----------------------------------------------------------------------------
dynamic_breaks_dense <- function(x) {
  rng  <- range(x, na.rm = TRUE)
  span <- diff(rng)
  
  if (!all(is.finite(rng)) || !is.finite(span) || span <= 0) {
    return(rng[1])
  }
  
  by <- if (span > 10) span / 15 else 0.5
  
  seq(
    from = floor(rng[1] / by) * by,
    to   = ceiling(rng[2] / by) * by,
    by   = by
  )
}

one_decimal <- label_number(accuracy = 0.1)

# -----------------------------------------------------------------------------
# 6) Diagnostics plots
# -----------------------------------------------------------------------------
while (!is.null(dev.list())) {
  dev.off()
}

# Global plotting constants
threshold_line_size <- 2.2
axis_title_size     <- 24
axis_text_size      <- 20
plot_title_size     <- 16
plot_subtitle_size  <- 12
point_alpha         <- 0.72
axis_line_size      <- 1.4

# New thresholds
thr_unreliable   <- 0.4
thr_reliable     <- 0.8
thr_conservative <- 0.5

col_unreliable   <- "#D55E00"
col_borderline   <- "#E6AB02"
col_reliable     <- "#009E73"
col_conservative <- "#0072B2"

threshold_colors <- c(
  "Unreliable threshold (< 0.4)"              = col_unreliable,
  "Reliable threshold (>= 0.8)"               = col_reliable,
  "Conservative error threshold (<= 0.5 cm)"  = col_conservative
)

diag_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title.x  = element_text(size = axis_title_size),
    axis.title.y  = element_text(size = axis_title_size),
    axis.text.x   = element_text(size = axis_text_size),
    axis.text.y   = element_text(size = axis_text_size),
    plot.title    = element_text(size = plot_title_size),
    plot.subtitle = element_text(size = plot_subtitle_size),
    legend.position      = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.title         = element_blank(),
    legend.background    = element_rect(fill = "white", colour = "grey70"),
    legend.text          = element_text(size = 11),
    panel.grid.minor     = element_line(linewidth = 0.25),
    panel.grid.major     = element_line(linewidth = 0.35)
  )

mt <- mt %>%
  mutate(
    match_zone = case_when(
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin < thr_unreliable ~ "unreliable",
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin >= thr_reliable ~ "reliable",
      !is.na(relative_margin) & is.finite(relative_margin) ~ "borderline",
      TRUE ~ NA_character_
    )
  )

# Correlation for Graph B
cor_df_B <- mt %>%
  filter(
    !is.na(abs_diff_cm), is.finite(abs_diff_cm),
    !is.na(ambiguity_margin_cm), is.finite(ambiguity_margin_cm)
  )

cor_test_B <- cor.test(
  cor_df_B$ambiguity_margin_cm,
  cor_df_B$abs_diff_cm,
  method = "pearson"
)

cor_label_B <- paste0(
  "Correlation (r) = ", format(round(unname(cor_test_B$estimate), 3), nsmall = 3),
  "\np = ", format(round(cor_test_B$p.value, 3), nsmall = 3)
)

legend_x_lines_margin <- tibble::tibble(
  x = c(thr_unreliable, thr_reliable),
  threshold = factor(
    c("Unreliable threshold (< 0.4)", "Reliable threshold (>= 0.8)"),
    levels = c(
      "Unreliable threshold (< 0.4)",
      "Reliable threshold (>= 0.8)",
      "Conservative error threshold (<= 0.5 cm)"
    )
  )
)

legend_x_lines_relative <- tibble::tibble(
  x = c(thr_unreliable, thr_reliable),
  threshold = factor(
    c("Unreliable threshold (< 0.4)", "Reliable threshold (>= 0.8)"),
    levels = c(
      "Unreliable threshold (< 0.4)",
      "Reliable threshold (>= 0.8)",
      "Conservative error threshold (<= 0.5 cm)"
    )
  )
)

legend_y_lines_error <- tibble::tibble(
  y = thr_conservative,
  threshold = factor(
    "Conservative error threshold (<= 0.5 cm)",
    levels = c(
      "Unreliable threshold (< 0.4)",
      "Reliable threshold (>= 0.8)",
      "Conservative error threshold (<= 0.5 cm)"
    )
  )
)

# -----------------------------------------------------------------------------
# Plot A
# -----------------------------------------------------------------------------
legend_x_lines_margin_A <- tibble::tibble(
  x = c(thr_unreliable, thr_reliable),
  threshold = factor(
    c("Unreliable threshold (< 0.4)", "Reliable threshold (>= 0.8)"),
    levels = c("Unreliable threshold (< 0.4)", "Reliable threshold (>= 0.8)")
  )
)

threshold_colors_A <- c(
  "Unreliable threshold (< 0.4)" = col_unreliable,
  "Reliable threshold (>= 0.8)"  = col_reliable
)

histA <- hist(
  mt$ambiguity_margin_cm,
  breaks = seq(
    floor(min(mt$ambiguity_margin_cm, na.rm = TRUE) / 0.25) * 0.25,
    ceiling(max(mt$ambiguity_margin_cm, na.rm = TRUE) / 0.25) * 0.25 + 0.25,
    by = 0.25
  ),
  plot = FALSE
)

max_count_A <- max(histA$counts, na.rm = TRUE)

pA <- ggplot(mt, aes(x = ambiguity_margin_cm)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "grey75") +
  geom_vline(
    data = legend_x_lines_margin_A,
    aes(xintercept = x, color = threshold),
    linetype = "dashed",
    linewidth = threshold_line_size
  ) +
  scale_color_manual(values = threshold_colors_A, drop = FALSE) +
  scale_x_continuous(
    breaks = dynamic_breaks_dense,
    labels = one_decimal
  ) +
  scale_y_continuous(
    breaks = seq(0, max_count_A, by = 1),
    limits = c(0, max_count_A)
  ) +
  labs(
    title = "A) Distribution of ambiguity margins",
    x = "Ambiguity margin (cm)",
    y = "Number of matches"
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = "dashed",
    linewidth = threshold_line_size
  ))) +
  diag_theme

# -----------------------------------------------------------------------------
# Plot B
# -----------------------------------------------------------------------------
xmax_B <- max(mt$ambiguity_margin_cm, na.rm = TRUE)
ymax_B <- max(mt$abs_diff_cm, na.rm = TRUE)

quad_B <- tibble::tribble(
  ~xmin,             ~xmax,           ~ymin,             ~ymax,             ~fill,       ~label,                                       ~xlab,                                 ~ylab,
  -Inf,              thr_unreliable,  -Inf,              thr_conservative,  "#F6D6D6",   "",                                           thr_unreliable * 0.50,                thr_conservative * 0.50,
  thr_reliable,      Inf,             -Inf,              thr_conservative,  "#D8F0D2",   "Best matches",                               thr_reliable + (xmax_B - thr_reliable) * 0.55, thr_conservative * 0.50,
  -Inf,              thr_unreliable,  thr_conservative,  Inf,               "#F8D9C4",   "",                                           thr_unreliable * 0.50,                thr_conservative + (ymax_B - thr_conservative) * 0.55,
  thr_reliable,      Inf,             thr_conservative,  Inf,               "#D8E6F8",   "Clear identity but\nlarge measurement error", thr_reliable + (xmax_B - thr_reliable) * 0.55, thr_conservative + (ymax_B - thr_conservative) * 0.55
)

pB <- ggplot(mt, aes(x = ambiguity_margin_cm, y = abs_diff_cm)) +
  annotate("rect", xmin = quad_B$xmin[1], xmax = quad_B$xmax[1], ymin = quad_B$ymin[1], ymax = quad_B$ymax[1], fill = quad_B$fill[1], alpha = 1) +
  annotate("rect", xmin = quad_B$xmin[2], xmax = quad_B$xmax[2], ymin = quad_B$ymin[2], ymax = quad_B$ymax[2], fill = quad_B$fill[2], alpha = 1) +
  annotate("rect", xmin = quad_B$xmin[3], xmax = quad_B$xmax[3], ymin = quad_B$ymin[3], ymax = quad_B$ymax[3], fill = quad_B$fill[3], alpha = 1) +
  annotate("rect", xmin = quad_B$xmin[4], xmax = quad_B$xmax[4], ymin = quad_B$ymin[4], ymax = quad_B$ymax[4], fill = quad_B$fill[4], alpha = 1) +
  geom_point(alpha = point_alpha) +
  geom_vline(
    data = legend_x_lines_margin,
    aes(xintercept = x, color = threshold),
    linetype = "dashed",
    linewidth = threshold_line_size
  ) +
  geom_hline(
    data = legend_y_lines_error,
    aes(yintercept = y, color = threshold),
    linetype = "dashed",
    linewidth = threshold_line_size
  ) +
  annotate("text", x = quad_B$xlab[2], y = quad_B$ylab[2], label = quad_B$label[2], size = 4.7, fontface = "bold", alpha = 0.92) +
  annotate("text", x = quad_B$xlab[4], y = quad_B$ylab[4], label = quad_B$label[4], size = 4.7, fontface = "bold", alpha = 0.92) +
  annotate("text",
           x = thr_unreliable,
           y = thr_conservative + (ymax_B - thr_conservative) * 0.50,
           label = "Bad matches",
           angle = 90,
           vjust = -0.3,
           hjust = 0.5,
           size = 4.8,
           fontface = "bold") +
  annotate("text",
           x = Inf, y = Inf,
           label = cor_label_B,
           hjust = 1.03, vjust = 3.8,
           size = 4.2) +
  scale_color_manual(values = threshold_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = dynamic_breaks_dense,
    labels = one_decimal
  ) +
  scale_y_continuous(labels = one_decimal) +
  labs(
    title = "B) Ambiguity margin vs absolute length error",
    x = "Ambiguity margin (cm)",
    y = "Absolute length error (cm)"
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = "dashed",
    linewidth = threshold_line_size
  ))) +
  diag_theme +
  theme(
    axis.line = element_line(linewidth = axis_line_size, colour = "black"),
    axis.ticks = element_line(linewidth = axis_line_size, colour = "black")
  )

# -----------------------------------------------------------------------------
# Plot C
# -----------------------------------------------------------------------------
x_ticks_C <- trial_spacing %>%
  filter(!is.na(min_spacing_cm), is.finite(min_spacing_cm)) %>%
  count(min_spacing_cm, name = "n_trials") %>%
  filter(n_trials > 0) %>%
  arrange(min_spacing_cm) %>%
  pull(min_spacing_cm)

x_ticks_C <- x_ticks_C[round(x_ticks_C, 1) != 1.0]

spacing_vals <- trial_spacing$min_spacing_cm[
  is.finite(trial_spacing$min_spacing_cm)
]

max_spacing <- max(spacing_vals, na.rm = TRUE)

hist_C <- hist(
  spacing_vals,
  breaks = seq(
    0,
    ceiling(max_spacing / 0.1) * 0.1,
    by = 0.1
  ),
  plot = FALSE
)

yrange_C <- c(0, max(hist_C$counts, na.rm = TRUE))

grad_C <- expand.grid(
  x = seq(0, 1.9, length.out = 400),
  y = seq(0, yrange_C[2], length.out = 200)
)
grad_C$ease <- grad_C$x

pC <- ggplot() +
  geom_raster(
    data = grad_C,
    aes(x = x, y = y, fill = ease),
    alpha = 0.40,
    interpolate = TRUE,
    inherit.aes = FALSE
  ) +
  geom_histogram(
    data = trial_spacing,
    aes(x = min_spacing_cm),
    binwidth = 0.1,
    color = "black",
    fill = "grey75"
  ) +
  annotate("text", x = 0.15, y = yrange_C[2] * 0.92, label = "hard match", color = "#8B0000", fontface = "bold", size = 5, hjust = 0) +
  annotate("text", x = 1.85, y = yrange_C[2] * 0.92, label = "easy match", color = "#006400", fontface = "bold", size = 5, hjust = 1) +
  geom_segment(
    aes(x = 0, xend = 2.0, y = -0.25, yend = -0.25),
    linewidth = axis_line_size,
    colour = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_gradient(
    low = "#D73027",
    high = "#1A9850"
  ) +
  scale_x_continuous(
    limits = c(0, 2.0),
    breaks = x_ticks_C,
    labels = one_decimal
  ) +
  labs(
    title = "C) Trial-level fish separability",
    subtitle = "Minimum difference between sampled fish lengths within a trial",
    x = "Minimum inter-fish length spacing (cm)",
    y = "Number of trials"
  ) +
  coord_cartesian(
    ylim = c(-0.25, yrange_C[2]),
    clip = "off"
  ) +
  diag_theme +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Plot D
# -----------------------------------------------------------------------------
legend_x_lines_relative_D <- tibble::tibble(
  x = c(thr_unreliable, thr_reliable),
  threshold = factor(
    c("Unreliable threshold (< 0.4)", "Reliable threshold (>= 0.8)"),
    levels = c("Unreliable threshold (< 0.4)", "Reliable threshold (>= 0.8)")
  )
)

threshold_colors_D <- c(
  "Unreliable threshold (< 0.4)" = col_unreliable,
  "Reliable threshold (>= 0.8)"  = col_reliable
)

mt_rel_plot <- mt %>%
  filter(!is.na(relative_margin), is.finite(relative_margin), relative_margin > 0)

xrange_D <- range(mt_rel_plot$relative_margin, na.rm = TRUE)

grad_breaks_D <- exp(seq(log(xrange_D[1]), log(xrange_D[2]), length.out = 300))

grad_rect_D <- tibble::tibble(
  xmin = head(grad_breaks_D, -1),
  xmax = tail(grad_breaks_D, -1)
) %>%
  mutate(
    ymin = 0,
    ymax = 35,
    xmid = sqrt(xmin * xmax)
  )

# initialize
grad_rect_D$fill_col <- "#FFFFFF"

# region 1: x < thr_unreliable  --> red gradient
idx_red <- which(grad_rect_D$xmid < thr_unreliable)
if (length(idx_red) > 0) {
  red_fun <- scales::col_numeric(
    palette = c("#FDE0DD", "#D73027"),
    domain = range(grad_rect_D$xmid[idx_red], na.rm = TRUE)
  )
  grad_rect_D$fill_col[idx_red] <- red_fun(grad_rect_D$xmid[idx_red])
}

# region 2: thr_unreliable <= x < thr_reliable  --> yellow gradient
idx_yellow <- which(grad_rect_D$xmid >= thr_unreliable & grad_rect_D$xmid < thr_reliable)
if (length(idx_yellow) > 0) {
  yellow_fun <- scales::col_numeric(
    palette = c("#FFF7BC", "#FFD92F"),
    domain = range(grad_rect_D$xmid[idx_yellow], na.rm = TRUE)
  )
  grad_rect_D$fill_col[idx_yellow] <- yellow_fun(grad_rect_D$xmid[idx_yellow])
}

# region 3: x >= thr_reliable  --> green gradient
idx_green <- which(grad_rect_D$xmid >= thr_reliable)
if (length(idx_green) > 0) {
  green_fun <- scales::col_numeric(
    palette = c("#D9F0D3", "#1A9850"),
    domain = range(grad_rect_D$xmid[idx_green], na.rm = TRUE)
  )
  grad_rect_D$fill_col[idx_green] <- green_fun(grad_rect_D$xmid[idx_green])
}

pD <- ggplot() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x  = element_text(size = axis_title_size),
    axis.title.y  = element_text(size = axis_title_size),
    axis.text.x   = element_text(size = axis_text_size),
    axis.text.y   = element_text(size = axis_text_size),
    plot.title    = element_text(size = plot_title_size),
    plot.subtitle = element_text(size = plot_subtitle_size),
    legend.position = c(0.98, 0.72),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "grey70"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(colour = "grey85"),
    panel.grid.minor = element_line(colour = "grey92")
  ) +
  geom_rect(
    data = grad_rect_D,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = grad_rect_D$fill_col,
    alpha = 0.40,
    inherit.aes = FALSE,
    color = NA
  ) +
  geom_histogram(
    data = mt_rel_plot,
    aes(x = relative_margin),
    binwidth = 0.25,
    color = "black",
    fill = "grey75"
  ) +
  geom_segment(
    data = legend_x_lines_relative_D,
    aes(x = x, xend = x, y = 0, yend = 35, color = threshold),
    linetype = "dashed",
    linewidth = threshold_line_size,
    inherit.aes = FALSE
  ) +
  annotate("text", x = min(xrange_D) * 1.05, y = 33, label = "unreliable match", color = "#8B0000", fontface = "bold", size = 5, hjust = 0) +
  annotate("text", x = max(xrange_D) / 1.05, y = 33, label = "reliable match", color = "#006400", fontface = "bold", size = 5, hjust = 1) +
  scale_color_manual(values = threshold_colors_D, drop = FALSE) +
  scale_x_log10(
    breaks = sort(unique(c(0.4, 0.8, 1, 2, 5, 10, 20, 50))),
    labels = one_decimal
  ) +
  scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, by = 5),
    expand = c(0, 0)
  ) +
  labs(
    title = "D) Relative ambiguity",
    subtitle = "Ambiguity margin scaled by within-trial fish-length separation",
    x = "Relative ambiguity (dimensionless, log scale)",
    y = "Number of matches",
    color = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = "dashed",
        linewidth = threshold_line_size
      )
    )
  ) +
  coord_cartesian(xlim = xrange_D, ylim = c(0, 35), clip = "on")

# -----------------------------------------------------------------------------
# Plot E
# -----------------------------------------------------------------------------
xmax_E <- max(mt_rel_plot$relative_margin, na.rm = TRUE)
ymax_E <- max(mt_rel_plot$abs_diff_cm, na.rm = TRUE)

quad_E <- tibble::tribble(
  ~xmin,            ~xmax,           ~ymin,             ~ymax,             ~fill,       ~label,                               ~xlab,                                  ~ylab,
  0,                thr_reliable,    0,                 thr_conservative,  "#F7E2C8",   "borderline matches",                 sqrt(0.08 * thr_reliable),              thr_conservative * 0.50,
  thr_reliable,     xmax_E * 1.05,   0,                 thr_conservative,  "#D8F0D2",   "Best matches",                       sqrt(thr_reliable * xmax_E),            thr_conservative * 0.50,
  0,                thr_unreliable,  thr_conservative,  ymax_E * 1.02,     "#F6D6D6",   "Unreliable identity",                sqrt(0.05 * thr_unreliable),            thr_conservative + (ymax_E - thr_conservative) * 0.55,
  thr_reliable,     xmax_E * 1.05,   thr_conservative,  ymax_E * 1.02,     "#D8E6F8",   "Clear identity but\nlarge measurement error", sqrt(thr_reliable * xmax_E),   thr_conservative + (ymax_E - thr_conservative) * 0.55
)

pE <- ggplot(mt_rel_plot, aes(x = relative_margin, y = abs_diff_cm)) +
  annotate("rect", xmin = quad_E$xmin[1], xmax = quad_E$xmax[1], ymin = quad_E$ymin[1], ymax = quad_E$ymax[1], fill = quad_E$fill[1], alpha = 1) +
  annotate("rect", xmin = quad_E$xmin[2], xmax = quad_E$xmax[2], ymin = quad_E$ymin[2], ymax = quad_E$ymax[2], fill = quad_E$fill[2], alpha = 1) +
  annotate("rect", xmin = quad_E$xmin[3], xmax = quad_E$xmax[3], ymin = quad_E$ymin[3], ymax = quad_E$ymax[3], fill = quad_E$fill[3], alpha = 1) +
  annotate("rect", xmin = quad_E$xmin[4], xmax = quad_E$xmax[4], ymin = quad_E$ymin[4], ymax = quad_E$ymax[4], fill = quad_E$fill[4], alpha = 1) +
  geom_point(alpha = point_alpha) +
  geom_vline(
    data = legend_x_lines_relative,
    aes(xintercept = x, color = threshold),
    linetype = "dashed",
    linewidth = threshold_line_size
  ) +
  geom_hline(
    data = legend_y_lines_error,
    aes(yintercept = y, color = threshold),
    linetype = "dashed",
    linewidth = threshold_line_size
  ) +
  annotate("text", x = quad_E$xlab, y = quad_E$ylab, label = quad_E$label, size = 4.7, fontface = "bold", alpha = 0.92) +
  scale_color_manual(values = threshold_colors, drop = FALSE) +
  scale_x_log10(
    breaks = sort(unique(c(0.4, 0.8, 1, 2, 5, 10, 20, 50))),
    labels = one_decimal
  ) +
  scale_y_continuous(labels = one_decimal) +
  labs(
    title = "E) Relative ambiguity vs absolute length error",
    subtitle = "Single best summary diagnostic for identity certainty and measurement accuracy",
    x = "Relative ambiguity (dimensionless, log scale)",
    y = "Absolute length error (cm)"
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = "dashed",
    linewidth = threshold_line_size
  ))) +
  diag_theme +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = axis_line_size, colour = "black")
  )

# -----------------------------------------------------------------------------
# Main single diagnostic plot
# -----------------------------------------------------------------------------
p_main <- pE +
  geom_segment(
    aes(x = min(mt_rel_plot$relative_margin, na.rm = TRUE),
        xend = max(mt_rel_plot$relative_margin, na.rm = TRUE),
        y = 0, yend = 0),
    inherit.aes = FALSE,
    linewidth = axis_line_size,
    colour = "black"
  ) +
  labs(
    title = "Main diagnostic: relative ambiguity vs absolute length error",
    subtitle = "Recommended single summary plot for match reliability"
  ) +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = axis_line_size, colour = "black"),
    axis.ticks = element_line(linewidth = axis_line_size, colour = "black")
  )

print(pA)
print(pB)
print(pC)
print(pD)
print(pE)
print(p_main)

p_grid <- (pA | pB) / (pC | pD)
print(p_grid)

# -----------------------------------------------------------------------------
# 7) Add QC flags
#    Internal logic:
#      unreliable  : relative_margin < 0.4
#      borderline  : 0.4 <= relative_margin < 0.8
#      reliable    : relative_margin >= 0.8
#      conservative: reliable == YES AND abs_diff_cm <= 0.5
#
#    Export logic:
#      Keep only one zone column called "adaptive":
#        unreliable / borderline / reliable
# -----------------------------------------------------------------------------
mt2 <- mt %>%
  mutate(
    unreliable = if_else(
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin < thr_unreliable,
      "YES", "NO"
    ),
    borderline = if_else(
      !is.na(relative_margin) & is.finite(relative_margin) &
        relative_margin >= thr_unreliable & relative_margin < thr_reliable,
      "YES", "NO"
    ),
    reliable = if_else(
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin >= thr_reliable,
      "YES", "NO"
    ),
    conservative = if_else(
      reliable == "YES" & !is.na(abs_diff_cm) & is.finite(abs_diff_cm) & abs_diff_cm <= thr_conservative,
      "YES", "NO"
    ),
    adaptive = case_when(
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin < thr_unreliable ~ "unreliable",
      !is.na(relative_margin) & is.finite(relative_margin) &
        relative_margin >= thr_unreliable & relative_margin < thr_reliable ~ "borderline",
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin >= thr_reliable ~ "reliable",
      TRUE ~ NA_character_
    ),
    adaptive_color = case_when(
      adaptive == "unreliable" ~ "red",
      adaptive == "borderline" ~ "yellow",
      adaptive == "reliable" ~ "green",
      TRUE ~ NA_character_
    )
  ) %>%
  relocate(
    min_spacing_cm, relative_margin, adaptive, adaptive_color, conservative,
    .after = ambiguity_margin_cm
  )

assign("match_table", mt2, envir = .GlobalEnv)

message("\nAdded/updated columns in match_table: min_spacing_cm, relative_margin, adaptive, adaptive_color, conservative")

message("\nCounts (adaptive zone):")
print(table(mt2$adaptive, useNA = "ifany"))

message("\nCounts (conservative = reliable + abs_diff <= 0.5):")
print(table(mt2$conservative, useNA = "ifany"))

message("\nBy segment_id (mean abs diff, mean ambiguity, mean rel margin, % by adaptive zone):")
print(
  mt2 %>%
    group_by(segment_id) %>%
    summarise(
      mean_abs_diff = mean(abs_diff_cm, na.rm = TRUE),
      mean_margin = mean(ambiguity_margin_cm, na.rm = TRUE),
      mean_rel_margin = mean(relative_margin, na.rm = TRUE),
      pct_reliable = mean(adaptive == "reliable", na.rm = TRUE) * 100,
      pct_borderline = mean(adaptive == "borderline", na.rm = TRUE) * 100,
      pct_unreliable = mean(adaptive == "unreliable", na.rm = TRUE) * 100,
      pct_conservative = mean(conservative == "YES", na.rm = TRUE) * 100,
      .groups = "drop"
    )
)

# -----------------------------------------------------------------------------
# 8) Add trial metadata and build export table
#    Treatment = tank from sampling data
# -----------------------------------------------------------------------------
trial_metadata_raw <- df %>%
  mutate(
    trial = suppressWarnings(as.integer(trial)),
    Treatment = if ("tank" %in% names(.)) as.character(tank) else NA_character_,
    Group = if ("group" %in% names(.)) as.character(group) else NA_character_
  ) %>%
  select(any_of(c("trial", "Group", "Treatment"))) %>%
  distinct()

trial_metadata <- trial_metadata_raw %>%
  group_by(trial) %>%
  summarise(
    Group = {
      vals <- na.omit(as.character(Group))
      if (length(vals) == 0) NA_character_ else vals[1]
    },
    Treatment = {
      vals <- na.omit(as.character(Treatment))
      if (length(vals) == 0) NA_character_ else vals[1]
    },
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Build pretty export table
# - keep only adaptive (zone), not unreliable/borderline/reliable flags
# - add pretty labels for match metrics
# -----------------------------------------------------------------------------
match_table_export <- mt2 %>%
  left_join(trial_metadata, by = "trial") %>%
  mutate(
    trial_id = as.integer(suppressWarnings(as.integer(trial))),
    reliability_tier = dplyr::case_when(
      !is.na(relative_margin) & is.finite(relative_margin) & conservative == "YES" ~ "conservative",
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin >= thr_reliable ~ "reliable",
      !is.na(relative_margin) & is.finite(relative_margin) & relative_margin >= thr_unreliable ~ "borderline",
      !is.na(relative_margin) & is.finite(relative_margin) ~ "unreliable",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    trial,
    trial_id,
    Group,
    Treatment,
    segment_id,
    df_sample_id,
    master_fish_id,
    df_lenght_cm,
    bl_cm,
    abs_diff_cm,
    ambiguity_margin_cm,
    min_spacing_cm,
    relative_margin,
    reliability_tier,
    adaptive,
    adaptive_color,
    conservative
  ) %>%
  rename(
    timepoint = segment_id,
    `sampling FISH ID` = df_sample_id,
    `tracking FISH ID` = master_fish_id,
    `sampling lenght (cm)` = df_lenght_cm,
    `tracking length (cm)` = bl_cm,
    `absolute length difference (cm)` = abs_diff_cm,
    `ambiguity margin (cm)` = ambiguity_margin_cm,
    `minimum spacing (cm)` = min_spacing_cm,
    `relative margin` = relative_margin,
    `adaptive zone` = adaptive,
    `adaptive color` = adaptive_color,
    `conservative flag` = conservative
  )

assign("match_table_export", match_table_export, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# 9) Console summary by Group * Treatment
# -----------------------------------------------------------------------------
if (all(c("Group", "Treatment", "adaptive zone") %in% names(match_table_export))) {
  message("\nNumber of individuals by adaptive zone and Group * Treatment:")
  print(
    match_table_export %>%
      group_by(Group, Treatment) %>%
      summarise(
        n_individuals_reliable = dplyr::n_distinct(`sampling FISH ID`[`adaptive zone` == "reliable"]),
        n_individuals_borderline = dplyr::n_distinct(`sampling FISH ID`[`adaptive zone` == "borderline"]),
        n_individuals_borderline_plus_reliable = dplyr::n_distinct(
          `sampling FISH ID`[`adaptive zone` %in% c("borderline", "reliable")]
        ),
        .groups = "drop"
      ) %>%
      arrange(Group, Treatment)
  )
} else {
  message("\nGroup/Treatment metadata not available; skipped Group * Treatment adaptive-zone summary.")
}

# -----------------------------------------------------------------------------
# 10) Export match_table_export in same directory as morph/tracking output
# -----------------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_name  <- paste0("match_table_export_", timestamp, ".xlsx")

out_dir <- NULL

if (length(csv_paths) > 0) {
  out_dir <- dirname(csv_paths[[1]])
} else if (dir.exists(master_output)) {
  out_dir <- master_output
} else {
  out_dir <- getwd()
}

out_path <- file.path(out_dir, out_name)

writexl::write_xlsx(
  list(match_table_export = match_table_export),
  path = out_path
)

message("\nSaved Excel: ", out_path)

# Export a compact CSV for downstream BL-tier gating (point BL_MATCH_PATH here)
bl_tier_csv_name <- paste0("bl_reliability_tier_", timestamp, ".csv")
bl_tier_csv_path <- file.path(out_dir, bl_tier_csv_name)

bl_tier_export <- match_table_export %>%
  dplyr::select(
    trial_id,
    reliability_tier,
    `sampling FISH ID`,
    `tracking FISH ID`,
    timepoint
  ) %>%
  dplyr::filter(!is.na(trial_id)) %>%
  dplyr::arrange(trial_id, timepoint)

readr::write_csv(bl_tier_export, bl_tier_csv_path)
message("Saved BL-tier CSV (set as BL_MATCH_PATH in monoamine/plasma/aggression): ", bl_tier_csv_path)

# ---- R4: stable canonical copy + auto-set BL_MATCH_PATH --------------------
# Side workflows (and STEP 5/6) can fall back to this fixed path when
# BL_MATCH_PATH is not explicitly set in the R session.
.bl_canonical_dir  <- file.path(master_output, "master_outputs", "reliability")
dir.create(.bl_canonical_dir, recursive = TRUE, showWarnings = FALSE)
.bl_canonical_path <- file.path(.bl_canonical_dir, "reliability_tier_LATEST.csv")
file.copy(bl_tier_csv_path, .bl_canonical_path, overwrite = TRUE)
message("R4: also wrote canonical copy: ", .bl_canonical_path)

# Expose BL_MATCH_PATH to GlobalEnv so STEP 5/6 + side workflows pick it up.
# Use the canonical (stable) path so re-runs in fresh R sessions still work.
assign("BL_MATCH_PATH", .bl_canonical_path, envir = .GlobalEnv)
message("R4: BL_MATCH_PATH set to canonical: ", .bl_canonical_path)
rm(.bl_canonical_dir, .bl_canonical_path)

# Print tier distribution
message("\nReliability tier distribution across all rows:")
print(table(bl_tier_export$reliability_tier, useNA = "ifany"))

################################################################################

# =============================================================================
# 11) IDENTITY RECOVERY SIMULATION (IRS)  [wrapped in tryCatch — diagnostic only]
#     - plug-and-play validation of the matching algorithm
#     - estimates probability of recovering each matched identity
#     - ALL critical outputs (BL tier CSV, reliability_tier_LATEST.csv) were
#       already written before this section; a failure here is non-fatal.
# =============================================================================
tryCatch({

# -----------------------------------------------------------------------------
# Assumptions used here
# -----------------------------------------------------------------------------
# 1) The currently assigned df_sample_id <-> master_fish_id pairs in mt2 are used
#    as the provisional "true" identities for simulation.
# 2) Tracking noise is approximated from the empirical signed residuals:
#       residual_cm = bl_cm - df_lenght_cm
# 3) Simulations are run within each trial x segment_id group.
# 4) If a group-specific noise SD cannot be estimated, a global fallback SD is used.

# -----------------------------------------------------------------------------
# Prepare simulation input
# -----------------------------------------------------------------------------
irs_input <- mt2 %>%
  transmute(
    trial          = suppressWarnings(as.integer(trial)),
    segment_id     = as.character(segment_id),
    df_sample_id   = as.character(df_sample_id),
    master_fish_id = as.character(master_fish_id),
    true_length_cm = suppressWarnings(as.numeric(df_lenght_cm)),
    tracked_length_cm = suppressWarnings(as.numeric(bl_cm)),
    residual_cm    = tracked_length_cm - true_length_cm,
    abs_diff_cm    = suppressWarnings(as.numeric(abs_diff_cm)),
    ambiguity_margin_cm = suppressWarnings(as.numeric(ambiguity_margin_cm)),
    min_spacing_cm = suppressWarnings(as.numeric(min_spacing_cm)),
    relative_margin = suppressWarnings(as.numeric(relative_margin))
  ) %>%
  filter(
    !is.na(trial),
    !is.na(segment_id),
    !is.na(df_sample_id),
    !is.na(master_fish_id),
    !is.na(true_length_cm),
    is.finite(true_length_cm)
  ) %>%
  arrange(trial, segment_id, master_fish_id)

if (nrow(irs_input) == 0) {
  stop("IRS could not run: no valid rows in mt2.")
}

global_sigma_cm <- stats::sd(irs_input$residual_cm, na.rm = TRUE)

if (!is.finite(global_sigma_cm) || is.na(global_sigma_cm) || global_sigma_cm <= 0) {
  # fallback if residual SD cannot be estimated
  global_sigma_cm <- median(irs_input$abs_diff_cm, na.rm = TRUE)
}

if (!is.finite(global_sigma_cm) || is.na(global_sigma_cm) || global_sigma_cm <= 0) {
  global_sigma_cm <- 0.25
}

message("\nIRS global fallback sigma (cm): ", round(global_sigma_cm, 4))

# -----------------------------------------------------------------------------
# Helper: Hungarian assignment for one simulated group
# -----------------------------------------------------------------------------
irs_assign_group <- function(master_ids, observed_lengths, true_ids, true_lengths) {
  observed_lengths <- as.numeric(observed_lengths)
  true_lengths     <- as.numeric(true_lengths)
  
  n_obs  <- length(observed_lengths)
  n_true <- length(true_lengths)
  
  if (n_obs == 0 || n_true == 0) {
    return(tibble(
      master_fish_id = character(),
      assigned_df_sample_id = character()
    ))
  }
  
  cost0 <- abs(outer(observed_lengths, true_lengths, FUN = "-"))
  
  n <- max(n_obs, n_true)
  cost <- matrix(1e6, nrow = n, ncol = n)
  cost[seq_len(n_obs), seq_len(n_true)] <- cost0
  
  assignment <- clue::solve_LSAP(cost)
  assigned_true_idx <- as.integer(assignment)[seq_len(n_obs)]
  assigned_true_idx[assigned_true_idx > n_true] <- NA_integer_
  
  tibble(
    master_fish_id = as.character(master_ids),
    assigned_df_sample_id = as.character(true_ids[assigned_true_idx])
  )
}

# -----------------------------------------------------------------------------
# Plug-and-play IRS function
# -----------------------------------------------------------------------------
run_identity_recovery_simulation <- function(match_df,
                                             n_sim = 5000L,
                                             seed = 123,
                                             sigma_mode = c("group", "global"),
                                             sigma_floor_cm = 0.05,
                                             global_sigma_cm = NULL) {
  sigma_mode <- match.arg(sigma_mode)
  
  sim_df <- match_df %>%
    transmute(
      trial          = suppressWarnings(as.integer(trial)),
      segment_id     = as.character(segment_id),
      df_sample_id   = as.character(df_sample_id),
      master_fish_id = as.character(master_fish_id),
      true_length_cm = as.numeric(true_length_cm),
      residual_cm    = as.numeric(residual_cm),
      relative_margin = as.numeric(relative_margin),
      adaptive = if ("adaptive" %in% names(match_df)) as.character(adaptive) else NA_character_
    ) %>%
    filter(
      !is.na(trial),
      !is.na(segment_id),
      !is.na(df_sample_id),
      !is.na(master_fish_id),
      !is.na(true_length_cm),
      is.finite(true_length_cm)
    ) %>%
    arrange(trial, segment_id, master_fish_id)
  
  if (is.null(global_sigma_cm) || !is.finite(global_sigma_cm) || is.na(global_sigma_cm) || global_sigma_cm <= 0) {
    global_sigma_cm <- stats::sd(sim_df$residual_cm, na.rm = TRUE)
  }
  if (!is.finite(global_sigma_cm) || is.na(global_sigma_cm) || global_sigma_cm <= 0) {
    global_sigma_cm <- sigma_floor_cm
  }
  
  group_noise <- sim_df %>%
    group_by(trial, segment_id) %>%
    summarise(
      sigma_group_cm = stats::sd(residual_cm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sigma_group_cm = if_else(
        !is.na(sigma_group_cm) & is.finite(sigma_group_cm) & sigma_group_cm > 0,
        pmax(sigma_group_cm, sigma_floor_cm),
        NA_real_
      )
    )
  
  sim_df <- sim_df %>%
    left_join(group_noise, by = c("trial", "segment_id")) %>%
    mutate(
      sigma_used_cm = dplyr::case_when(
        sigma_mode == "group" &
          !is.na(sigma_group_cm) & is.finite(sigma_group_cm) & sigma_group_cm > 0 ~ sigma_group_cm,
        TRUE ~ global_sigma_cm
      ),
      sigma_used_cm = pmax(sigma_used_cm, sigma_floor_cm)
    )
  
  group_keys <- sim_df %>%
    distinct(trial, segment_id) %>%
    arrange(trial, segment_id)
  
  set.seed(seed)
  
  correct_counts <- sim_df %>%
    transmute(
      trial, segment_id, master_fish_id,
      correct_n = 0L
    )
  
  group_recovery <- vector("list", length = nrow(group_keys))
  
  for (g in seq_len(nrow(group_keys))) {
    tr  <- group_keys$trial[[g]]
    seg <- group_keys$segment_id[[g]]
    
    gdat <- sim_df %>%
      filter(trial == tr, segment_id == seg) %>%
      arrange(master_fish_id)
    
    master_ids   <- gdat$master_fish_id
    true_ids     <- gdat$df_sample_id
    true_lengths <- gdat$true_length_cm
    sigma_used   <- gdat$sigma_used_cm
    
    sim_correct_matrix <- matrix(FALSE, nrow = nrow(gdat), ncol = n_sim)
    
    for (i in seq_len(n_sim)) {
      sim_obs_lengths <- stats::rnorm(
        n = length(true_lengths),
        mean = true_lengths,
        sd = sigma_used
      )
      
      sim_assign <- irs_assign_group(
        master_ids = master_ids,
        observed_lengths = sim_obs_lengths,
        true_ids = true_ids,
        true_lengths = true_lengths
      ) %>%
        arrange(master_fish_id)
      
      sim_correct_matrix[, i] <- sim_assign$assigned_df_sample_id == true_ids
    }
    
    row_prob <- rowMeans(sim_correct_matrix, na.rm = TRUE)
    
    group_recovery[[g]] <- gdat %>%
      transmute(
        trial,
        segment_id,
        master_fish_id,
        df_sample_id,
        sigma_used_cm,
        irs_prob_correct = row_prob,
        irs_prob_incorrect = 1 - row_prob,
        irs_expected_recovery_group = mean(row_prob, na.rm = TRUE)
      )
  }
  
  irs_table <- bind_rows(group_recovery) %>%
    left_join(
      sim_df %>%
        select(trial, segment_id, master_fish_id, relative_margin, adaptive) %>%
        distinct(),
      by = c("trial", "segment_id", "master_fish_id")
    ) %>%
    mutate(
      irs_prob_correct_pct = irs_prob_correct * 100,
      irs_prob_bin = cut(
        irs_prob_correct,
        breaks = c(-Inf, 0.80, 0.95, 0.99, Inf),
        labels = c("<80%", "80-95%", "95-99%", ">=99%"),
        right = FALSE
      )
    ) %>%
    arrange(trial, segment_id, master_fish_id)
  
  irs_summary <- irs_table %>%
    summarise(
      n_matches = n(),
      mean_prob_correct = mean(irs_prob_correct, na.rm = TRUE),
      median_prob_correct = median(irs_prob_correct, na.rm = TRUE),
      min_prob_correct = min(irs_prob_correct, na.rm = TRUE),
      max_prob_correct = max(irs_prob_correct, na.rm = TRUE)
    )
  
  list(
    irs_table = irs_table,
    irs_summary = irs_summary,
    n_sim = n_sim,
    seed = seed,
    sigma_mode = sigma_mode,
    global_sigma_cm = global_sigma_cm
  )
}

# -----------------------------------------------------------------------------
# Run IRS on current dataset
# -----------------------------------------------------------------------------
irs_res <- run_identity_recovery_simulation(
  match_df = mt2 %>%
    transmute(
      trial,
      segment_id,
      df_sample_id,
      master_fish_id,
      true_length_cm = df_lenght_cm,
      residual_cm = bl_cm - df_lenght_cm,
      relative_margin,
      adaptive
    ),
  n_sim = 5000L,
  seed = 123,
  sigma_mode = "group",
  sigma_floor_cm = 0.05,
  global_sigma_cm = global_sigma_cm
)

identity_recovery_table <- mt2 %>%
  left_join(
    irs_res$irs_table %>%
      select(
        trial, segment_id, master_fish_id, df_sample_id,
        sigma_used_cm, irs_prob_correct, irs_prob_correct_pct,
        irs_prob_incorrect, irs_prob_bin, irs_expected_recovery_group
      ),
    by = c("trial", "segment_id", "master_fish_id", "df_sample_id")
  ) %>%
  arrange(trial, segment_id, master_fish_id)

assign("irs_res", irs_res, envir = .GlobalEnv)
assign("identity_recovery_table", identity_recovery_table, envir = .GlobalEnv)

message("\nIRS complete.")
print(irs_res$irs_summary)

# -----------------------------------------------------------------------------
# Publication-quality validation figure
# Figure shows:
#   - point-level estimated probability of correct identity
#   - adaptive zones
#   - binned summary with 95% CI
# -----------------------------------------------------------------------------
irs_plot_df <- identity_recovery_table %>%
  filter(
    !is.na(relative_margin),
    is.finite(relative_margin),
    relative_margin > 0,
    !is.na(irs_prob_correct),
    is.finite(irs_prob_correct)
  ) %>%
  mutate(
    adaptive = factor(
      adaptive,
      levels = c("unreliable", "borderline", "reliable")
    )
  )

# bin summary for smooth interpretation
irs_bin_summary <- irs_plot_df %>%
  mutate(
    rm_bin = cut(
      relative_margin,
      breaks = c(0, thr_unreliable, thr_reliable, Inf),
      include.lowest = TRUE,
      right = FALSE,
      labels = c(
        paste0("<", thr_unreliable),
        paste0(thr_unreliable, "–", thr_reliable),
        paste0("≥", thr_reliable)
      )
    )
  ) %>%
  group_by(rm_bin) %>%
  summarise(
    n = n(),
    mean_prob = mean(irs_prob_correct, na.rm = TRUE),
    se = stats::sd(irs_prob_correct, na.rm = TRUE) / sqrt(n),
    lwr = pmax(mean_prob - 1.96 * se, 0),
    upr = pmin(mean_prob + 1.96 * se, 1),
    .groups = "drop"
  ) %>%
  mutate(
    x = dplyr::case_when(
      as.character(rm_bin) == paste0("<", thr_unreliable) ~ thr_unreliable / 2,
      as.character(rm_bin) == paste0(thr_unreliable, "–", thr_reliable) ~ (thr_unreliable + thr_reliable) / 2,
      as.character(rm_bin) == paste0("≥", thr_reliable) ~ max(irs_plot_df$relative_margin, na.rm = TRUE) * 0.7,
      TRUE ~ NA_real_
    )
  )

p_irs_validation <- ggplot(irs_plot_df, aes(x = relative_margin, y = irs_prob_correct)) +
  annotate("rect", xmin = 0, xmax = thr_unreliable, ymin = 0, ymax = 1, fill = "#F6D6D6", alpha = 0.9) +
  annotate("rect", xmin = thr_unreliable, xmax = thr_reliable, ymin = 0, ymax = 1, fill = "#F7E2C8", alpha = 0.9) +
  annotate("rect", xmin = thr_reliable, xmax = max(irs_plot_df$relative_margin, na.rm = TRUE) * 1.05, ymin = 0, ymax = 1, fill = "#D8F0D2", alpha = 0.9) +
  geom_point(alpha = 0.75, size = 2.4) +
  geom_vline(xintercept = c(thr_unreliable, thr_reliable), linetype = "dashed", linewidth = 1.2) +
  geom_hline(yintercept = c(0.80, 0.95, 0.99), linetype = "dotted", linewidth = 0.7, color = "grey40") +
  geom_errorbar(
    data = irs_bin_summary,
    aes(x = x, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    width = 0.03,
    linewidth = 0.9
  ) +
  geom_point(
    data = irs_bin_summary,
    aes(x = x, y = mean_prob),
    inherit.aes = FALSE,
    shape = 23,
    size = 4,
    fill = "white",
    color = "black",
    stroke = 1
  ) +
  scale_x_log10(
    breaks = sort(unique(c(0.4, 0.8, 1, 2, 5, 10, 20, 50))),
    labels = one_decimal
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      unreliable = col_unreliable,
      borderline = col_borderline,
      reliable = col_reliable
    ),
    drop = FALSE
  ) +
  labs(
    title = "Identity recovery simulation validation",
    subtitle = "Estimated probability of correct identity recovery by relative margin",
    x = "Relative margin (log scale)",
    y = "Estimated probability of correct identity",
    color = "Adaptive zone"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey70")
  )

print(p_irs_validation)
assign("p_irs_validation", p_irs_validation, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# Save publication-quality IRS figure
# -----------------------------------------------------------------------------
irs_fig_name_png <- paste0("IRS_validation_figure_", timestamp, ".png")
irs_fig_name_pdf <- paste0("IRS_validation_figure_", timestamp, ".pdf")

irs_fig_path_png <- file.path(out_dir, irs_fig_name_png)
irs_fig_path_pdf <- file.path(out_dir, irs_fig_name_pdf)

ggsave(
  filename = irs_fig_path_png,
  plot = p_irs_validation,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = irs_fig_path_pdf,
  plot = p_irs_validation,
  width = 10,
  height = 7
)

message("\nSaved IRS validation figure:")
message(" - ", irs_fig_path_png)
message(" - ", irs_fig_path_pdf)

# -----------------------------------------------------------------------------
# Publication-ready IRS summary table by adaptive zone
# -----------------------------------------------------------------------------
irs_zone_summary <- identity_recovery_table %>%
  mutate(
    adaptive = factor(adaptive, levels = c("unreliable", "borderline", "reliable"))
  ) %>%
  group_by(adaptive) %>%
  summarise(
    n_matches = n(),
    mean_prob_correct = mean(irs_prob_correct, na.rm = TRUE),
    median_prob_correct = median(irs_prob_correct, na.rm = TRUE),
    min_prob_correct = min(irs_prob_correct, na.rm = TRUE),
    max_prob_correct = max(irs_prob_correct, na.rm = TRUE),
    .groups = "drop"
  )

assign("irs_zone_summary", irs_zone_summary, envir = .GlobalEnv)

message("\nIRS summary by adaptive zone:")
print(irs_zone_summary)

# -----------------------------------------------------------------------------
# Add IRS columns to export table if match_table_export already exists
# -----------------------------------------------------------------------------
if (exists("match_table_export", envir = .GlobalEnv)) {
  match_table_export <- match_table_export %>%
    left_join(
      identity_recovery_table %>%
        transmute(
          trial,
          timepoint = segment_id,
          `sampling FISH ID` = df_sample_id,
          `tracking FISH ID` = master_fish_id,
          `estimated ID recovery probability` = irs_prob_correct,
          `estimated ID recovery probability (%)` = irs_prob_correct_pct,
          `IRS noise sigma used (cm)` = sigma_used_cm,
          `IRS probability bin` = as.character(irs_prob_bin)
        ),
      by = c("trial", "timepoint", "sampling FISH ID", "tracking FISH ID")
    )
  
  assign("match_table_export", match_table_export, envir = .GlobalEnv)
  message("\nIRS columns appended to match_table_export.")
}

# -----------------------------------------------------------------------------
# Optional export of the per-match IRS table
# -----------------------------------------------------------------------------
irs_table_name <- paste0("identity_recovery_table_", timestamp, ".xlsx")
irs_table_path <- file.path(out_dir, irs_table_name)

writexl::write_xlsx(
  list(
    identity_recovery_table = identity_recovery_table,
    irs_zone_summary = irs_zone_summary
  ),
  path = irs_table_path
)

message("\nSaved IRS table: ", irs_table_path)

# -----------------------------------------------------------------------------
# PATCH: add Group + treatment + IRS columns to export,
#        save a copy to morph dir,
#        and print console summary for IRS > 0.8
# -----------------------------------------------------------------------------

export_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# sampling metadata: Group + tank as treatment
sampling_meta_export <- df %>%
  mutate(
    trial = suppressWarnings(as.integer(trial)),
    Group = if ("group" %in% names(.)) as.character(group) else NA_character_,
    treatment = if ("tank" %in% names(.)) as.character(tank) else NA_character_
  ) %>%
  select(trial, Group, treatment) %>%
  distinct() %>%
  group_by(trial) %>%
  summarise(
    Group = {
      vals <- na.omit(as.character(Group))
      if (length(vals) == 0) NA_character_ else vals[1]
    },
    treatment = {
      vals <- na.omit(as.character(treatment))
      if (length(vals) == 0) NA_character_ else vals[1]
    },
    .groups = "drop"
  )

# add IRS info + metadata + timestamp
match_table_export <- match_table_export %>%
  select(-any_of(c("Group", "treatment", "timestamp",
                   "IRS identity recovery probability (%)",
                   "IRS reliability bin"))) %>%
  left_join(
    sampling_meta_export,
    by = "trial"
  ) %>%
  left_join(
    identity_recovery_table %>%
      transmute(
        trial,
        timepoint = segment_id,
        `sampling FISH ID` = df_sample_id,
        `tracking FISH ID` = master_fish_id,
        `IRS identity recovery probability (%)` = round(irs_prob_correct_pct, 1),
        `IRS reliability bin` = as.character(irs_prob_bin)
      ),
    by = c("trial", "timepoint", "sampling FISH ID", "tracking FISH ID")
  ) %>%
  mutate(
    timestamp = export_timestamp
  ) %>%
  relocate(
    timestamp,
    trial,
    Group,
    treatment,
    timepoint,
    `sampling FISH ID`,
    `tracking FISH ID`,
    `IRS identity recovery probability (%)`,
    `IRS reliability bin`
  )

assign("match_table_export", match_table_export, envir = .GlobalEnv)

message("\nAdded Group, treatment, IRS columns, and timestamp to match_table_export.")

# -----------------------------------------------------------------------------
# Save a copy to morph dir
# -----------------------------------------------------------------------------
timestamp_file <- format(Sys.time(), "%Y%m%d_%H%M%S")
morph_out_name <- paste0("match_table_export_", timestamp_file, ".xlsx")

morph_dir <- NULL
if (exists("csv_paths", envir = .GlobalEnv) && length(csv_paths) > 0) {
  morph_dir <- dirname(csv_paths[[1]])
} else if (exists("master_output", envir = .GlobalEnv) && dir.exists(master_output)) {
  morph_dir <- master_output
} else {
  morph_dir <- getwd()
}

morph_out_path <- file.path(morph_dir, morph_out_name)

writexl::write_xlsx(
  list(match_table_export = match_table_export),
  path = morph_out_path
)

message("\nSaved match_table_export to morph dir: ", morph_out_path)

# -----------------------------------------------------------------------------
# Console summary: number of IDs with IRS probability > 0.8 by treatment x Group
# -----------------------------------------------------------------------------
if (all(c("Group", "treatment", "IRS identity recovery probability (%)", "sampling FISH ID") %in% names(match_table_export))) {
  message("\nNumber of unique IDs with IRS identity recovery probability > 0.8 by treatment x Group:")
  print(
    match_table_export %>%
      mutate(
        irs_prob_prop = `IRS identity recovery probability (%)` / 100
      ) %>%
      group_by(treatment, Group) %>%
      summarise(
        n_ids_irs_gt_0_8 = dplyr::n_distinct(`sampling FISH ID`[irs_prob_prop > 0.8]),
        .groups = "drop"
      ) %>%
      arrange(treatment, Group)
  )
} else {
  message("\nSkipped IRS > 0.8 summary: required columns not found in match_table_export.")
}

# --- close IRS tryCatch -------------------------------------------------------
}, error = function(.irs_err) {
  warning(
    "IRS section failed (non-fatal — BL tier outputs were already written): ",
    conditionMessage(.irs_err),
    call. = FALSE
  )
  message("IRS skipped due to error: ", conditionMessage(.irs_err))
})
# (Removed: stray install.packages / remotes::install_github calls that were
# left in from interactive development. They crashed non-interactive Rscript
# runs. Pipeline does not require ggedit.)