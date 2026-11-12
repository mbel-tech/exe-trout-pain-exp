###############################################################################
# GOT activity statistics pipeline (faster cleaned version)
# Design:
# - unit of analysis: one fish observation at one timepoint within one trial
# - repeated factor: timepoint repeated within trial_id
# - fish_id nested under trial_id is included as a candidate RE
# - model selection: candidate random structures compared by AICc
# - inference: Type III Wald chi-square ANOVA via car::Anova()
# - outputs saved into: GOT_activity_stats_<timestamp>
###############################################################################

############################
# 0) Package checks
############################

# ---- output helpers (standalone guard) ----
if (!exists(".got_make_output_dir", mode = "function")) {
  tryCatch(
    source(file.path(getwd(), "02_trajectory_pipeline", "helpers", "output_helpers.R"), local = FALSE),
    error = function(e) tryCatch(
      source(file.path(getwd(), "helpers/output_helpers.R"), local = FALSE),
      error = function(e2) stop(
        "Cannot find 00_output_helpers.R. Run via master pipeline or set working directory appropriately.",
        call. = FALSE
      )
    )
  )
}

# ---- pretty labels (R3: centralized — sourced by master) ----
# Normally available in .GlobalEnv via 00_run_trajectory_pipeline.R.
# Fallback below only fires when this script is run standalone.
if (!exists("pretty_indicator_name", mode = "function")) {
  .this_script_dir <- tryCatch({
    frames <- sys.frames()
    ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
    if (length(ofile_envs) > 0L) {
      normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                    winslash = "/", mustWork = FALSE)
    } else getwd()
  }, error = function(e) getwd())

  .pl_path <- file.path(.this_script_dir, "helpers/pretty_labels.R")
  if (file.exists(.pl_path)) {
    source(.pl_path, local = FALSE)
  } else {
    message("Note: pretty_labels.R not found at ", .pl_path,
            " — using built-in labels only.")
  }
  rm(.this_script_dir, .pl_path)
}

required_pkgs <- c(
  "dplyr", "tidyr", "tibble",
  "lme4", "lmerTest", "MuMIn", "DHARMa",
  "car", "emmeans",
  "ggplot2",
  "officer",
  "flextable",
  "stringr",
  "multcomp",
  "multcompView",
  "readxl",
  "future", "future.apply"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Install required packages first:\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(car)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(officer)
  library(flextable)
  library(ggplot2)
  library(DHARMa)
  library(readxl)
  library(future)
  library(future.apply)
})

############################
# 1) Global constants
############################

options(contrasts = c("contr.sum", "contr.poly"))

PIPELINE_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- .got_make_output_dir("STEP5_output")

DIR_REPORTS <- file.path(OUTDIR, "reports")
DIR_TABLES  <- file.path(OUTDIR, "tables")
DIR_PLOTS   <- file.path(OUTDIR, "plots")
DIR_RDS     <- file.path(OUTDIR, "r_objects")
DIR_LOGS    <- file.path(OUTDIR, "logs")

dir.create(DIR_REPORTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLES,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PLOTS,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_RDS,     recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOGS,    recursive = TRUE, showWarnings = FALSE)

MAX_DECIMALS <- 3
DOC_PAGE_WIDTH_IN <- 8.27
DOC_PAGE_HEIGHT_IN <- 11.69
P_VALUE_CUTOFF <- 0.05

# Full candidate set (fish_id_matched_trial version — used when idmatcher qualifies).
# (1|trial_id) is mandatory in every candidate — it is the experimental unit.
# fish_id_matched_trial = trial-scoped stable identity across seg2+seg3.
RANDOM_CANDIDATES_FULL <- c(
  "(1 | trial_id)",
  "(0 + timepoint | trial_id)",
  "(1 + timepoint | trial_id)",
  "(1 | trial_id) + (1 | fish_id_matched_trial)",
  "(1 | trial_id/fish_id_matched_trial)"
)

# Fallback when no idmatcher tier qualifies — trial_id RE only.
RANDOM_CANDIDATES_NO_FISH <- c(
  "(1 | trial_id)",
  "(0 + timepoint | trial_id)",
  "(1 + timepoint | trial_id)"
)

# Active candidate set — overwritten by gate_random_candidates() at startup
RANDOM_CANDIDATES <- RANDOM_CANDIDATES_FULL

# Gate fish RE candidates based on idmatcher reliability tier (primary) or
# Hungarian BL-match tier (fallback when idmatcher data not available).
#
# Idmatcher path (preferred):
#   - Reads idmatcher_tier from idmatcher_tier_df (GlobalEnv, written by STEP 1.5).
#   - If any trial meets IDMATCHER_MIN_TIER, use RANDOM_CANDIDATES_FULL with
#     (1|fish_id_matched_trial) — correct repeated-measures identity.
#   - If no trial qualifies, use RANDOM_CANDIDATES_NO_FISH.
#
# Hungarian fallback (bl_match_df, used only if idmatcher not available):
#   - Reads reliability_tier / hungarian_tier from the merged BL_MATCH_PATH CSV.
#   - Includes legacy (1|fish_id) candidates only if any trial meets min_tier.
gate_random_candidates <- function(bl_match_df = NULL, min_tier = NULL) {
  if (is.null(min_tier)) {
    min_tier <- if (exists("MIN_BL_RELIABILITY_TIER", envir = .GlobalEnv, inherits = FALSE))
      get("MIN_BL_RELIABILITY_TIER", envir = .GlobalEnv)
    else "reliable"
  }

  # ── Idmatcher path (primary) ─────────────────────────────────────────────
  idmatcher_min_tier <- if (exists("IDMATCHER_MIN_TIER", envir = .GlobalEnv, inherits = FALSE))
    get("IDMATCHER_MIN_TIER", envir = .GlobalEnv) else min_tier

  if (exists("idmatcher_tier_df", envir = .GlobalEnv, inherits = FALSE) &&
      is.data.frame(idmatcher_tier_df) &&
      "idmatcher_tier" %in% names(idmatcher_tier_df)) {
    qualifying_ids <- qualifying_trial_ids(idmatcher_tier_df, idmatcher_min_tier,
                                           tier_col = "idmatcher_tier")
    n_qualifying <- length(qualifying_ids)
    if (n_qualifying == 0) {
      message("STEP5: No trials meet idmatcher tier '", idmatcher_min_tier,
              "' — excluding fish RE candidates.")
      return(RANDOM_CANDIDATES_NO_FISH)
    }
    message("STEP5: ", n_qualifying, " trial(s) meet idmatcher tier '",
            idmatcher_min_tier, "' — using (1|fish_id_matched_trial) RE candidates.")
    return(RANDOM_CANDIDATES_FULL)
  }

  # ── Hungarian fallback ───────────────────────────────────────────────────
  if (is.null(bl_match_df) || !is.data.frame(bl_match_df)) {
    message("STEP5: No idmatcher or BL match table available — using full RE candidate set.")
    return(RANDOM_CANDIDATES_FULL)
  }
  # Try hungarian_tier column first (merged CSV), then reliability_tier (legacy)
  tier_col <- if ("hungarian_tier" %in% names(bl_match_df)) "hungarian_tier"
              else if ("reliability_tier" %in% names(bl_match_df)) "reliability_tier"
              else NULL
  if (is.null(tier_col)) {
    message("STEP5: BL match table has no tier column — using full RE candidate set.")
    return(RANDOM_CANDIDATES_FULL)
  }
  tier_order <- c("unreliable", "borderline", "reliable", "conservative")
  min_idx    <- match(tolower(trimws(min_tier)), tier_order)
  if (is.na(min_idx)) min_idx <- 3L
  tier_vals    <- match(tolower(trimws(bl_match_df[[tier_col]])), tier_order)
  n_qualifying <- sum(!is.na(tier_vals) & tier_vals >= min_idx)
  if (n_qualifying == 0) {
    message("STEP5: No trials meet Hungarian tier '", min_tier,
            "' — excluding fish RE candidates.")
    return(RANDOM_CANDIDATES_NO_FISH)
  }
  message("STEP5 (Hungarian fallback): ", n_qualifying, " trial(s) meet tier '", min_tier,
          "' — including (1|fish_id_matched_trial) RE candidates.")
  RANDOM_CANDIDATES_FULL
}

build_fixed_formula <- function(data) {
  tp_vals <- unique(as.character(stats::na.omit(data$timepoint)))
  tp_vals <- tp_vals[nzchar(tp_vals)]
  
  if (length(tp_vals) <= 1) {
    return("Treatment * Group")
  }
  
  "timepoint * Treatment * Group"
}

DIR_REPORTS <- NULL
DIR_TABLES  <- NULL
DIR_PLOTS   <- NULL
DIR_RDS     <- NULL
DIR_LOGS    <- NULL

set_output_dirs <- function(base_dir) {
  assign("DIR_REPORTS", file.path(base_dir, "reports"), envir = .GlobalEnv)
  assign("DIR_TABLES",  file.path(base_dir, "tables"),  envir = .GlobalEnv)
  assign("DIR_PLOTS",   file.path(base_dir, "plots"),   envir = .GlobalEnv)
  assign("DIR_RDS",     file.path(base_dir, "r_objects"), envir = .GlobalEnv)
  assign("DIR_LOGS",    file.path(base_dir, "logs"),    envir = .GlobalEnv)
  
  dir.create(get("DIR_REPORTS", envir = .GlobalEnv), recursive = TRUE, showWarnings = FALSE)
  dir.create(get("DIR_TABLES",  envir = .GlobalEnv), recursive = TRUE, showWarnings = FALSE)
  dir.create(get("DIR_PLOTS",   envir = .GlobalEnv), recursive = TRUE, showWarnings = FALSE)
  dir.create(get("DIR_RDS",     envir = .GlobalEnv), recursive = TRUE, showWarnings = FALSE)
  dir.create(get("DIR_LOGS",    envir = .GlobalEnv), recursive = TRUE, showWarnings = FALSE)
}

# Simpler candidate set used only during assumption diagnostics
ASSUMPTION_RANDOM_CANDIDATES <- c(
  "(1 | trial_id)"
)

ALLOW_LOG_TRANSFORM  <- TRUE
ALLOW_SQRT_TRANSFORM <- TRUE

ANALYZE_ANGULAR_AS_LINEAR <- FALSE
# Always run post-hoc and CLD regardless of interaction p-value (manuscript requirement)
RUN_CLD_ONLY_IF_SIG <- FALSE
RUN_POSTHOC_ONLY_IF_SIG <- FALSE

# Set to an integer to restrict the deadlike-fish exclusion to a specific trial only.
# Set to NULL to apply the rule to ALL trials (use with caution — broadens exclusion scope).
TRIAL_WITH_DEAD_FISH <- 10L

# --- Analysis scope switches ------------------------------------------------
# OPT-0: derive scope switches from master pipeline globals when available;
# fall back to the local defaults if running STEP5 standalone.
#
# KEEP_TP1 (global) = TRUE  → USE_TP1 = TRUE  (include TP1 subbranch)
# KEEP_TP1 (global) = FALSE → USE_TP1 = FALSE (skip TP1 subbranch)
USE_TP1 <- if (exists("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE)) {
  isTRUE(get("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE))
} else {
  FALSE   # standalone default: exclude TP1
}

# DROP_PHARMA (global) = TRUE  → EXCLUDE_PHARMA_GROUPS = TRUE  (exe_choice only)
# DROP_PHARMA (global) = FALSE → EXCLUDE_PHARMA_GROUPS = TRUE  (same default)
# (Set to FALSE explicitly in the global env if you want the all_groups layout.)
EXCLUDE_PHARMA_GROUPS <- if (exists("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE)) {
  isTRUE(get("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE))
} else {
  TRUE    # standalone default: exclude pharma groups
}
IMMOBILE_THRESHOLD <- 0.95
IMMOBILE_VARS <- c(
  "time_moving_TS0.5_s",
  "time_moving_TS1_s",
  "time_moving_TS1.5_s",
  "time_moving_TS2_s"
)

RECORDING_DURATION_SECONDS <- 600

ALL_INDICATORS <- c(
  "distance_m",
  "distance_BL",
  "speed_cm_s_mean",
  "speed_cm_s_median",
  "speed_cm_s_max",
  "speed_BL_s_mean",
  "speed_BL_s_median",
  "speed_BL_s_max",
  "accel_cm_s2_mean",
  "accel_cm_s2_median",
  "accel_cm_s2_max",
  "n_zone_transitions",
  "time_moving_TS0.5_s",
  "time_moving_TS1_s",
  "time_spent_zone_flow_s"
)

CONTINUOUS_INDICATORS <- c(
  "distance_m",
  "distance_BL",
  "speed_cm_s_mean",
  "speed_cm_s_median",
  "speed_cm_s_max",
  "speed_BL_s_mean",
  "speed_BL_s_median",
  "speed_BL_s_max",
  "accel_cm_s2_mean",
  "accel_cm_s2_median",
  "accel_cm_s2_max",
  "time_moving_TS0.5_s",
  "time_moving_TS1_s"
)

# Bounded [0, 1] outcomes modelled with Beta GLMM (glmmTMB::beta_family).
# The raw seconds column is converted to a proportion internally; plots always
# show the original seconds scale.
BETA_INDICATORS <- c("time_spent_zone_flow_s")

COUNT_INDICATORS <- c("n_zone_transitions")

ANGULAR_INDICATORS <- c(
  "heading_rad_mean",
  "heading_rad_median",
  "heading_rad_max",
  "turn_angle_rad_mean",
  "turn_angle_rad_median",
  "turn_angle_rad_max"
)

REFERENCE_DOCX <- normalizePath(
  if (exists("REFERENCE_DOCX_PATH", envir = .GlobalEnv, inherits = FALSE))
    get("REFERENCE_DOCX_PATH", envir = .GlobalEnv)
  else
    "data/reference_a4_portrait.docx",
  winslash = "/",
  mustWork = FALSE
)
if (!file.exists(REFERENCE_DOCX)) {
  if (exists("resolve_dataset_path", mode = "function")) {
    REFERENCE_DOCX <- resolve_dataset_path(REFERENCE_DOCX)
  } else {
    warning("REFERENCE_DOCX not found; Word reports will use a blank template: ", REFERENCE_DOCX,
            call. = FALSE)
  }
}

############################
# Officer compatibility patch
############################

.officer_read_docx <- function() {
  if (file.exists(REFERENCE_DOCX)) {
    return(officer::read_docx(path = REFERENCE_DOCX))
  }
  officer::read_docx()
}

.officer_body_add_par <- officer::body_add_par
.officer_body_add_img <- officer::body_add_img

.officer_print <- function(doc, path) {
  print(doc, target = path)
}

.officer_body_add_flextable <- function(doc, ft) {
  if ("body_add_flextable" %in% getNamespaceExports("officer")) {
    return(officer::body_add_flextable(doc, value = ft))
  }
  if ("body_add_flextable" %in% getNamespaceExports("flextable")) {
    return(flextable::body_add_flextable(doc, value = ft))
  }
  stop("No usable body_add_flextable() found.")
}

# No-op wrapper: do nothing and return the document unchanged
.officer_body_end_section_continuous <- function(doc, section) {
  doc
}

############################
# 3) Small utilities
############################

fmt_num <- function(x, digits = MAX_DECIMALS) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_p <- function(p, digits = 3) {
  # Report exact value to 3 decimal places when p > 0.0001.
  # Use threshold notation only when p <= 0.0001.
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- NA_character_
  out[!is.na(p) & p <= 0.0001] <- "< 0.0001"
  idx <- !is.na(p) & p > 0.0001
  if (any(idx)) {
    out[idx] <- formatC(p[idx], format = "f", digits = digits)
  }
  out
}

sig_symbol <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*",
                         ifelse(p < 0.1, ".", "")
                  )
           )
    )
  )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

INDICATOR_LABELS <- c(
  distance_m = "Distance travelled (m)",
  distance_BL = "Distance travelled (body lengths)",
  speed_cm_s_mean = "Mean speed (cm/s)",
  speed_cm_s_median = "Median speed (cm/s)",
  speed_cm_s_max = "Maximum speed (cm/s)",
  speed_BL_s_mean = "Mean speed (body lengths/s)",
  speed_BL_s_median = "Median speed (body lengths/s)",
  speed_BL_s_max = "Maximum speed (body lengths/s)",
  heading_rad_mean = "Mean heading (rad)",
  heading_rad_median = "Median heading (rad)",
  heading_rad_max = "Maximum heading (rad)",
  turn_angle_rad_mean = "Mean turn angle (rad)",
  turn_angle_rad_median = "Median turn angle (rad)",
  turn_angle_rad_max = "Maximum turn angle (rad)",
  accel_cm_s2_mean = "Mean acceleration (cm/s²)",
  accel_cm_s2_median = "Median acceleration (cm/s²)",
  accel_cm_s2_max = "Maximum acceleration (cm/s²)",
  n_zone_transitions = "Number of zone transitions",
  time_moving_TS0.5_s = "Time moving (TS0.5; s)",
  time_moving_TS1_s = "Time moving (TS1; s)",
  time_spent_zone_flow_s = "Time spent in flow zone (s)"
)

pretty_indicator_name <- function(x) {
  x <- as.character(x)[1]
  lab <- unname(INDICATOR_LABELS[x])
  if (length(lab) == 1 && !is.na(lab) && nzchar(lab)) {
    return(lab)
  }
  gsub("_", " ", x)
}

pretty_colname <- function(x) {
  # Strip leading dot (internal helper columns)
  x <- gsub("^\\.(.*)$", "\\1", x)

  # ANOVA-specific exact-match substitutions (applied before generic transforms
  # so that special characters like "Pr(>F)" survive title-casing intact).
  fixed_map <- c(
    "Pr(>Chisq)"       = "p",
    "Pr(>F)"           = "p",
    "Chisq"            = "χ²",   # χ²
    "F value"          = "F",
    "NumDF"            = "df (num)",
    "DenDF"            = "df (den)",
    "p_value_numeric"  = "p",
    "p_value_BH"       = "p (BH-adj)",
    "significant"      = "p < 0.05",
    "significant_BH"   = "p (BH) < 0.05",
    "significance_label" = "Sig.",
    "stat_col"         = "Stat. col",
    "p_col"            = "p col"
  )
  if (x %in% names(fixed_map)) return(fixed_map[[x]])

  # Generic transforms
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(x)
  x <- gsub("\\bAicc\\b",      "AICc",      x)
  x <- gsub("\\bDf\\b",        "df",         x)   # lower-case df
  x <- gsub("\\bP Value\\b",   "p",          x)
  x <- gsub("\\bCld Letter\\b","CLD letter", x)
  x
}

drop_all_empty_cols <- function(df) {
  if (is.null(df) || !is.data.frame(df) || ncol(df) == 0) return(df)
  
  keep <- vapply(df, function(col) {
    if (all(is.na(col))) return(FALSE)
    vals <- trimws(as.character(col))
    any(!is.na(vals) & nzchar(vals))
  }, logical(1))
  
  df[, keep, drop = FALSE]
}

prettify_table_for_report <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  # Drop internal helper columns (names starting with ".")
  internal_cols <- grep("^\\.", names(df), value = TRUE)
  if (length(internal_cols) > 0)
    df <- df[, !names(df) %in% internal_cols, drop = FALSE]
  df <- drop_all_empty_cols(df)
  if (ncol(df) == 0) return(df)
  # Pretty-print ANOVA term names in the "term" column (": " → " × ", timepoint → Interval)
  if ("term" %in% names(df) && exists("pretty_term_name", mode = "function")) {
    df$term <- vapply(as.character(df$term), pretty_term_name, character(1L))
  }
  # Format p-value columns to exactly 3 decimal places (as character) BEFORE renaming,
  # so flextable cannot display them with floating-point rounding artefacts.
  p_val_cols <- intersect(c("p_value_numeric", "p_value_BH",
                             "Pr(>F)", "Pr(>Chisq)", "p.value", "p"),
                          names(df))
  for (pc in p_val_cols) {
    raw <- suppressWarnings(as.numeric(df[[pc]]))
    df[[pc]] <- ifelse(is.na(raw), NA_character_,
                       ifelse(raw <= 0.0001, "< 0.0001",
                              formatC(raw, format = "f", digits = 3L)))
  }
  # Format remaining numeric test-statistic columns to 3 decimal places
  stat_cols <- intersect(c("F value", "Chisq", "F", "Chi Sq"),
                         names(df))
  for (sc in stat_cols) {
    raw <- suppressWarnings(as.numeric(df[[sc]]))
    df[[sc]] <- ifelse(is.na(raw), NA_character_,
                       formatC(round(raw, 3L), format = "f", digits = 3L))
  }
  names(df) <- vapply(names(df), pretty_colname, character(1))
  df
}

summarize_cld_availability <- function(res) {
  if (is.null(res$cld_tables)) {
    return("No CLD tables object was stored in the result.")
  }
  
  if (!is.list(res$cld_tables)) {
    return("CLD tables object was present but was not a list.")
  }
  
  if (length(res$cld_tables) == 0) {
    if (!is.null(res$anova_table) && is.data.frame(res$anova_table)) {
      if (!has_non_intercept_sig_terms(res$anova_table)) {
        return("CLD tables not generated (model may have failed or had no eligible terms).")
      }
    }
    return("CLD section omitted because no CLD tables were generated.")
  }
  
  msgs <- character(0)
  
  for (nm in names(res$cld_tables)) {
    tbl <- res$cld_tables[[nm]]
    term_label <- gsub(":", " × ", nm, fixed = TRUE)
    
    if (is.null(tbl)) {
      msgs <- c(msgs, paste0(term_label, ": CLD object was NULL."))
      next
    }
    
    if (!is.data.frame(tbl)) {
      msgs <- c(msgs, paste0(term_label, ": CLD object was not a data.frame."))
      next
    }
    
    if (nrow(tbl) == 0) {
      msgs <- c(msgs, paste0(term_label, ": CLD table had 0 rows."))
      next
    }
    
    if ("note" %in% names(tbl)) {
      note_vals <- unique(trimws(as.character(tbl$note)))
      note_vals <- note_vals[nzchar(note_vals)]
      if (length(note_vals) > 0) {
        msgs <- c(msgs, paste0(term_label, ": ", paste(note_vals, collapse = " | ")))
      }
    }
    
    if ("cld_letter" %in% names(tbl) && all(is.na(tbl$cld_letter) | !nzchar(trimws(as.character(tbl$cld_letter))))) {
      msgs <- c(msgs, paste0(term_label, ": CLD letters were all missing/empty."))
      next
    }
  }
  
  if (length(msgs) == 0) {
    return("CLD section omitted because no printable CLD table was available.")
  }
  
  paste(msgs, collapse = "\n")
}

has_non_intercept_sig_terms <- function(anova_df, alpha = 0.05) {
  p_col <- attr(anova_df, "p_col")
  if (is.null(p_col) || is.na(p_col) || !p_col %in% names(anova_df)) return(FALSE)
  keep <- !is.na(anova_df[[p_col]]) & anova_df[[p_col]] < alpha & anova_df$term != "(Intercept)"
  any(keep)
}

trim_table_numeric <- function(df, digits = MAX_DECIMALS) {
  out <- df
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  out
}

save_csv <- function(df, filename) {
  utils::write.csv(df, file = filename, row.names = FALSE, na = "")
}

as_ft <- function(df, add_sig_col = FALSE, p_col = NULL) {
  x <- trim_table_numeric(df, digits = MAX_DECIMALS)
  
  if (add_sig_col && !is.null(p_col) && p_col %in% names(x) && !"signif" %in% names(x)) {
    x$signif <- sig_symbol(x[[p_col]])
  }
  
  ft <- flextable::flextable(x)
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)
  ft <- flextable::align(ft, align = "center", part = "all")
  ft <- flextable::theme_booktabs(ft)
  ft
}

append_text <- function(doc, text, style = "Normal") {
  available_styles <- tryCatch(
    officer::styles_info(doc)$style_name,
    error = function(e) c("Normal", "heading 1", "heading 2", "heading 3")
  )
  
  style_map <- c(
    "heading 3" = "heading 3",
    "heading 3" = "heading 3",
    "heading4"  = "heading 3"
  )
  
  if (style %in% names(style_map)) {
    style <- unname(style_map[[style]])
  }
  
  if (!style %in% available_styles) {
    style <- "Normal"
  }
  
  .officer_body_add_par(doc, value = text, style = style)
}

append_ft <- function(doc, df, add_sig_col = FALSE, p_col = NULL) {
  df <- prettify_table_for_report(df)
  if (is.null(df) || !is.data.frame(df) || ncol(df) == 0 || nrow(df) == 0) {
    return(append_text(doc, "No table to display.", style = "Normal"))
  }
  ft <- as_ft(df, add_sig_col = add_sig_col, p_col = p_col)
  .officer_body_add_flextable(doc, ft)
}

append_plot <- function(doc, img_path, width = 6.5, height = 4.5) {
  .officer_body_add_img(doc, src = img_path, width = width, height = height)
}

a4_section <- officer::prop_section(
  page_size = officer::page_size(
    orient = "portrait",
    width = DOC_PAGE_WIDTH_IN,
    height = DOC_PAGE_HEIGHT_IN
  ),
  page_margins = officer::page_mar(
    top = 0.6, bottom = 0.6, left = 0.6, right = 0.6
  ),
  type = "continuous"
)

append_page_break <- function(doc) {
  if ("run_pagebreak" %in% getNamespaceExports("officer")) {
    return(officer::body_add(doc, officer::run_pagebreak()))
  }
  
  # fallback for older officer versions
  if ("body_add_break" %in% getNamespaceExports("officer")) {
    return(officer::body_add_break(doc))
  }
  
  warning("No page-break function available in this officer version; continuing without explicit page break.")
  doc
}
############################
# 4) Input loader
############################

read_clipboard_path <- function() {
  clip_txt <- NULL
  
  if (.Platform$OS.type == "windows") {
    clip_txt <- tryCatch(utils::readClipboard(), error = function(e) NULL)
  }
  
  if ((is.null(clip_txt) || length(clip_txt) == 0 || !nzchar(clip_txt[1])) &&
      requireNamespace("clipr", quietly = TRUE)) {
    clip_txt <- tryCatch(as.character(clipr::read_clip()), error = function(e) NULL)
  }
  
  if (is.null(clip_txt) || length(clip_txt) == 0) return(NULL)
  
  path <- paste(clip_txt, collapse = "")
  path <- trimws(path)
  path <- gsub('^"|"$', "", path)
  path <- gsub("\\\\", "/", path)
  if (!nzchar(path)) return(NULL)
  path
}

load_dataset_from_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  
  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  
  if (ext %in% c("xlsx", "xls")) {
    sheets <- readxl::excel_sheets(path)
    target_sheet <- sheets[grepl("fish_activity_summary", sheets, ignore.case = TRUE)]
    if (length(target_sheet) == 0) target_sheet <- sheets[1]
    return(as.data.frame(
      readxl::read_excel(path, sheet = target_sheet[1]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!is.data.frame(obj)) stop("RDS object is not a data.frame: ", path)
    return(obj)
  }
  
  if (ext %in% c("rdata", "rda")) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (exists("fish_activity_summary", envir = env, inherits = FALSE)) {
      obj <- get("fish_activity_summary", envir = env)
      if (is.data.frame(obj)) return(obj)
    }
    df_names <- ls(env)[vapply(ls(env), function(nm) is.data.frame(get(nm, envir = env)), logical(1))]
    if (length(df_names) == 1) return(get(df_names[1], envir = env))
    stop("RData/Rda did not contain a unique data.frame or fish_activity_summary object: ", path)
  }
  
  stop("Unsupported file type: ", path)
}

locate_input_data <- function() {
  if (exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
    dat <- get("fish_activity_summary", envir = .GlobalEnv)
    if (!is.data.frame(dat)) stop("fish_activity_summary exists but is not a data.frame.")
    return(list(
      data = dat,
      source = "global_environment",
      source_path = NA_character_
    ))
  }
  
  # Build a list of candidate search directories in priority order:
  # 1) PIPELINE_DATA_DIR/dataset (R1: canonical fixed-name handoff — preferred)
  # 2) PIPELINE_DATA_DIR/master_outputs (legacy STEP 4 → STEP 5 handoff)
  # 3) STEP4_output/ folders (most-recent fallback for benchmark/standalone runs)
  script_root <- "."
  candidate_dirs <- character(0)

  pipeline_data_dir <- if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE))
    get("PIPELINE_DATA_DIR", envir = .GlobalEnv) else "data/tracking_sessions"
  candidate_dirs <- c(candidate_dirs,
                       file.path(pipeline_data_dir, "dataset"),           # R1 canonical
                       file.path(pipeline_data_dir, "master_outputs"))    # legacy

  step4_root <- file.path(script_root, "STEP4_output")
  if (dir.exists(step4_root)) {
    step4_subdirs <- list.dirs(step4_root, recursive = FALSE, full.names = TRUE)
    step4_subdirs <- step4_subdirs[order(file.mtime(step4_subdirs), decreasing = TRUE)]
    candidate_dirs <- c(candidate_dirs, step4_subdirs)
  }

  for (default_data_dir in candidate_dirs) {
    if (!dir.exists(default_data_dir)) next
    default_files <- list.files(
      default_data_dir,
      pattern = "fish_activity_summary.*\\.(csv|tsv|txt|xlsx|xls|rds|RDS|rda|RData)$",
      full.names = TRUE, recursive = FALSE, ignore.case = TRUE
    )
    if (length(default_files) >= 1) {
      default_files <- default_files[order(file.mtime(default_files), decreasing = TRUE)]
      dat <- tryCatch(load_dataset_from_file(default_files[1]), error = function(e) NULL)
      if (!is.null(dat) && is.data.frame(dat)) {
        message("locate_input_data: loaded fish_activity_summary from: ", default_files[1])
        return(list(data = dat, source = "default_folder", source_path = default_files[1]))
      }
    }
  }

  clip_path <- read_clipboard_path()
  if (is.null(clip_path)) {
    stop("Could not find fish_activity_summary in .GlobalEnv and clipboard fallback did not yield a readable path.")
  }
  
  if (!file.exists(clip_path) && !dir.exists(clip_path)) {
    stop("Clipboard path does not exist: ", clip_path)
  }
  
  candidate_files <- character(0)
  
  if (file.exists(clip_path) && !dir.exists(clip_path)) {
    candidate_files <- clip_path
  } else {
    files <- list.files(
      clip_path,
      pattern = "\\.(csv|tsv|txt|xlsx|xls|rds|RDS|rda|RData)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    files1 <- files[grepl("fish_activity_summary", basename(files), ignore.case = TRUE)]
    candidate_files <- c(files1, setdiff(files, files1))
  }
  
  if (length(candidate_files) == 0) {
    stop("No supported data files found in clipboard folder: ", clip_path)
  }
  
  last_error <- NULL
  for (f in candidate_files) {
    dat <- tryCatch(load_dataset_from_file(f), error = function(e) { last_error <<- e$message; NULL })
    if (!is.null(dat) && is.data.frame(dat)) {
      return(list(
        data = dat,
        source = "clipboard_path",
        source_path = f
      ))
    }
  }
  
  stop("Could not load a dataset from clipboard path.\nLast error: ", last_error)
}

############################
# 5) Data preparation and QC
############################

check_required_columns <- function(dat, cols) {
  missing_cols <- setdiff(cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns:\n", paste(missing_cols, collapse = ", "))
  }
}

coerce_design_columns <- function(dat) {
  dat$trial_id <- as.factor(dat$trial_id)
  dat$timepoint <- as.factor(dat$timepoint)
  dat$Group <- as.factor(dat$Group)
  dat$Treatment <- as.factor(dat$Treatment)
  dat
}

compute_immobile_flag <- function(dat) {
  if (!all(IMMOBILE_VARS %in% names(dat))) {
    dat$flag_dead_immobile <- FALSE
    dat$flag_dead_immobile_reason <- "Immobility variables not all present; rule not applied"
    return(dat)
  }
  
  total_duration <- RECORDING_DURATION_SECONDS
  if (is.na(total_duration)) {
    total_duration <- suppressWarnings(max(unlist(dat[IMMOBILE_VARS]), na.rm = TRUE))
    if (!is.finite(total_duration)) {
      dat$flag_dead_immobile <- FALSE
      dat$flag_dead_immobile_reason <- "Could not derive recording duration; rule not applied"
      return(dat)
    }
  }
  
  immobile_prop_ok <- function(x) {
    ifelse(is.na(x), FALSE, (1 - (x / total_duration)) > IMMOBILE_THRESHOLD)
  }
  
  all_immobile <- Reduce(`&`, lapply(IMMOBILE_VARS, function(v) immobile_prop_ok(dat[[v]])))

  # If TRIAL_WITH_DEAD_FISH is NULL, flag any fish meeting the immobility criterion.
  # Otherwise restrict to the specified trial only.
  in_scope <- if (is.null(TRIAL_WITH_DEAD_FISH)) {
    rep(TRUE, nrow(dat))
  } else {
    as.integer(as.character(dat$trial_id)) == TRIAL_WITH_DEAD_FISH
  }

  dat$flag_dead_immobile <- in_scope & all_immobile
  dat$flag_dead_immobile[is.na(dat$flag_dead_immobile)] <- FALSE

  trial_scope_label <- if (is.null(TRIAL_WITH_DEAD_FISH)) "any trial" else paste0("trial ", TRIAL_WITH_DEAD_FISH)
  dat$flag_dead_immobile_reason <- ifelse(
    dat$flag_dead_immobile,
    paste0("Excluded: ", trial_scope_label, ", >", round(IMMOBILE_THRESHOLD * 100), "% immobile across all four time_moving thresholds"),
    ""
  )
  
  dat
}

run_qc <- function(dat) {
  dat_qc <- dat
  
  dat_qc$structural_missing_trial13_tp3 <- with(
    dat_qc,
    as.character(trial_id) == "13" & as.character(timepoint) == "3"
  )
  dat_qc$structural_missing_trial13_tp3[is.na(dat_qc$structural_missing_trial13_tp3)] <- FALSE
  
  dat_qc <- compute_immobile_flag(dat_qc)
  dat_qc$exclude_row <- dat_qc$flag_dead_immobile
  
  qc_exclusion_table <- dat_qc[dat_qc$exclude_row, c(
    intersect(c("fish_id", "trial_id", "timepoint", "Group", "Treatment"), names(dat_qc)),
    "flag_dead_immobile_reason"
  ), drop = FALSE]
  
  if (nrow(qc_exclusion_table) == 0) {
    qc_exclusion_table <- data.frame(
      note = "No rows excluded by the trial 10 immobility rule",
      stringsAsFactors = FALSE
    )
  }
  
  qc_counts_by_trial_timepoint <- aggregate(
    rep(1, nrow(dat_qc[!dat_qc$exclude_row, ])) ~ trial_id + timepoint,
    data = dat_qc[!dat_qc$exclude_row, , drop = FALSE],
    FUN = sum
  )
  names(qc_counts_by_trial_timepoint)[names(qc_counts_by_trial_timepoint) == "rep(1, nrow(dat_qc[!dat_qc$exclude_row, ]))"] <- "n_rows_retained"
  
  missingness_table <- data.frame(
    trial_id = "13",
    timepoint = "3",
    note = "Structurally missing: recording not made",
    stringsAsFactors = FALSE
  )
  
  dat_clean <- dat_qc[!dat_qc$exclude_row, , drop = FALSE]
  
  list(
    data = dat_clean,
    qc_exclusion_table = qc_exclusion_table,
    qc_counts_by_trial_timepoint = qc_counts_by_trial_timepoint,
    missingness_table = missingness_table
  )
}

############################
# 6) Model helpers
############################

calc_aicc <- function(model) {
  ll <- logLik(model)
  k <- attr(ll, "df")
  n <- stats::nobs(model)
  aic <- AIC(model)
  if (is.na(n) || is.na(k) || (n - k - 1) <= 0) return(NA_real_)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

fit_model_safely <- function(formula, data, model_kind = c("gaussian", "poisson", "nb"), REML = FALSE) {
  model_kind <- match.arg(model_kind)
  warn_vec <- character(0)
  
  fit <- withCallingHandlers(
    tryCatch({
      if (model_kind == "gaussian") {
        lme4::lmer(stats::as.formula(formula), data = data, REML = REML)
      } else if (model_kind == "poisson") {
        lme4::glmer(stats::as.formula(formula), data = data, family = poisson())
      } else if (model_kind == "nb") {
        lme4::glmer.nb(stats::as.formula(formula), data = data)
      }
    }, error = function(e) e),
    warning = function(w) {
      warn_vec <<- c(warn_vec, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE,
      model = NULL,
      error = fit$message,
      warnings = unique(warn_vec),
      singular = NA,
      aicc = NA_real_,
      formula = formula
    ))
  }
  
  singular_flag <- tryCatch(lme4::isSingular(fit, tol = 1e-4), error = function(e) NA)
  
  list(
    ok = TRUE,
    model = fit,
    error = NA_character_,
    warnings = unique(warn_vec),
    singular = singular_flag,
    aicc = calc_aicc(fit),
    formula = formula
  )
}

drop_random_structure <- function(formula_string) {
  f <- formula_string
  if (grepl("\\(1 \\+ timepoint \\| trial_id\\)", f, fixed = FALSE)) {
    return(gsub("\\(1 \\+ timepoint \\| trial_id\\)", "(1 | trial_id)", f))
  }
  if (grepl("\\(0 \\+ timepoint \\| trial_id\\)", f, fixed = FALSE)) {
    return(gsub("\\(0 \\+ timepoint \\| trial_id\\)", "(1 | trial_id)", f))
  }
  if (grepl("\\(1 \\| trial_id\\)", f, fixed = FALSE)) {
    return(gsub("\\s*\\+\\s*\\(1 \\| trial_id\\)", "", f))
  }
  f
}

refit_if_singular <- function(fit_obj, data, model_kind, REML = FALSE, max_steps = 2) {
  current <- fit_obj
  history <- list()
  
  for (i in seq_len(max_steps)) {
    if (!isTRUE(current$ok) || is.null(current$model) || !isTRUE(current$singular)) break
    new_formula <- drop_random_structure(current$formula)
    if (identical(new_formula, current$formula)) break
    
    history[[length(history) + 1]] <- data.frame(
      step = i,
      old_formula = current$formula,
      new_formula = new_formula,
      stringsAsFactors = FALSE
    )
    
    current <- fit_model_safely(
      formula = new_formula,
      data = data,
      model_kind = model_kind,
      REML = REML
    )
  }
  
  current$singularity_refit_history <- if (length(history) > 0) do.call(rbind, history) else NULL
  current
}

rank_candidate_models <- function(data, response, model_kind = c("gaussian", "poisson", "nb"),
                                  transformed_name = NULL, REML = FALSE,
                                  random_candidates = RANDOM_CANDIDATES,
                                  auto_drop_singular = TRUE,
                                  fixed_formula = NULL) {
  model_kind <- match.arg(model_kind)
  yname <- if (is.null(transformed_name)) response else transformed_name
  
  if (is.null(fixed_formula)) {
    fixed_formula <- build_fixed_formula(data)
  }
  
  fits <- lapply(random_candidates, function(re_term) {
    frm <- paste(yname, "~", fixed_formula, "+", re_term)
    fit <- fit_model_safely(
      formula = frm,
      data = data,
      model_kind = model_kind,
      REML = REML
    )
    if (auto_drop_singular && isTRUE(fit$ok) && isTRUE(fit$singular)) {
      fit <- refit_if_singular(fit, data = data, model_kind = model_kind, REML = REML, max_steps = 2)
    }
    fit
  })
  
  tbl <- do.call(rbind, lapply(fits, function(x) {
    data.frame(
      formula = x$formula,
      ok = x$ok,
      singular = ifelse(is.na(x$singular), NA, x$singular),
      aicc = x$aicc,
      error = ifelse(is.na(x$error), "", x$error),
      warnings = paste(unique(x$warnings), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  
  tbl <- tbl[order(tbl$aicc, na.last = TRUE), , drop = FALSE]
  if (nrow(tbl) > 0 && any(is.finite(tbl$aicc))) {
    best_aicc <- min(tbl$aicc, na.rm = TRUE)
    tbl$delta_aicc <- tbl$aicc - best_aicc
  } else {
    tbl$delta_aicc <- NA_real_
  }
  
  list(table = tbl, fits = fits, fixed_formula = fixed_formula)
}

select_best_model <- function(rank_obj, allow_singular = FALSE, delta_prefer_nonsingular = 2) {
  ok_idx <- which(vapply(rank_obj$fits, function(x) isTRUE(x$ok), logical(1)))
  if (length(ok_idx) == 0) return(NULL)
  
  # Build a ranked summary aligned to fits
  fit_tbl <- data.frame(
    idx = seq_along(rank_obj$fits),
    ok = vapply(rank_obj$fits, function(x) isTRUE(x$ok), logical(1)),
    singular = vapply(rank_obj$fits, function(x) isTRUE(x$singular), logical(1)),
    aicc = vapply(rank_obj$fits, function(x) x$aicc, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  fit_tbl <- fit_tbl[fit_tbl$ok & is.finite(fit_tbl$aicc), , drop = FALSE]
  if (nrow(fit_tbl) == 0) return(NULL)
  
  fit_tbl <- fit_tbl[order(fit_tbl$aicc), , drop = FALSE]
  fit_tbl$delta_aicc <- fit_tbl$aicc - fit_tbl$aicc[1]
  
  # Absolute best by AICc
  best_row <- fit_tbl[1, , drop = FALSE]
  
  # If singular models are allowed, just return the best AICc model
  if (allow_singular) {
    return(rank_obj$fits[[best_row$idx]])
  }
  
  # If the best model is already non-singular, use it
  if (!isTRUE(best_row$singular)) {
    return(rank_obj$fits[[best_row$idx]])
  }
  
  # Otherwise, search among near-equivalent models (delta < threshold)
  near_tbl <- fit_tbl[fit_tbl$delta_aicc < delta_prefer_nonsingular, , drop = FALSE]
  near_nonsingular <- near_tbl[!near_tbl$singular, , drop = FALSE]
  
  if (nrow(near_nonsingular) > 0) {
    # Choose the best AICc among the near-equivalent non-singular models
    chosen_idx <- near_nonsingular$idx[1]
    return(rank_obj$fits[[chosen_idx]])
  }
  
  # If no near-equivalent non-singular model exists, fall back to the singular AICc-best
  return(rank_obj$fits[[best_row$idx]])
}

###############################################################################
# UPDATED DIAGNOSTIC HELPERS
# Replace your current versions of:
#   - get_residual_df()
#   - save_gaussian_diag_plots()
###############################################################################

get_residual_df <- function(model, data = NULL, response_name = NULL) {
  # Residuals
  resid_vec <- tryCatch(
    stats::residuals(model),
    error = function(e) NULL
  )
  
  if (is.null(resid_vec)) {
    stop("Could not extract residuals from model.")
  }
  
  resid_vec <- as.numeric(resid_vec)
  n_resid <- length(resid_vec)
  
  # Fitted values: prefer predict(), fallback to fitted()
  fitted_vec <- tryCatch(
    stats::predict(model),
    error = function(e) NULL
  )
  
  if (is.null(fitted_vec) || length(fitted_vec) == 0) {
    fitted_vec <- tryCatch(
      stats::fitted(model),
      error = function(e) NULL
    )
  }
  
  if (is.null(fitted_vec) || length(fitted_vec) == 0) {
    fitted_vec <- rep(NA_real_, n_resid)
  } else {
    fitted_vec <- as.numeric(fitted_vec)
  }
  
  # Use the model frame when available, because it reflects rows actually used
  mf <- tryCatch(
    stats::model.frame(model),
    error = function(e) NULL
  )
  
  if (!is.null(mf) && nrow(mf) > 0) {
    if (all(c("timepoint", "Treatment", "Group") %in% names(mf))) {
      group_combo <- interaction(
        mf$timepoint,
        mf$Treatment,
        mf$Group,
        drop = TRUE
      )
    } else {
      group_combo <- factor(rep("all", nrow(mf)))
    }
    
    if (!is.null(response_name) && response_name %in% names(mf)) {
      response_vec <- mf[[response_name]]
    } else {
      response_vec <- rep(NA_real_, nrow(mf))
    }
  } else {
    if (!is.null(data) && all(c("timepoint", "Treatment", "Group") %in% names(data))) {
      group_combo <- interaction(
        data$timepoint,
        data$Treatment,
        data$Group,
        drop = TRUE
      )
    } else {
      group_combo <- factor(rep("all", n_resid))
    }
    
    if (!is.null(data) && !is.null(response_name) && response_name %in% names(data)) {
      response_vec <- data[[response_name]]
    } else {
      response_vec <- rep(NA_real_, n_resid)
    }
  }
  
  # Harmonize lengths
  lens <- c(
    fitted = length(fitted_vec),
    resid = length(resid_vec),
    group_combo = length(group_combo),
    response = length(response_vec)
  )
  
  positive_lens <- lens[lens > 0]
  if (length(positive_lens) == 0) {
    stop("Could not build residual data frame; all extracted vectors had length 0.")
  }
  
  n <- min(positive_lens)
  
  data.frame(
    fitted = fitted_vec[seq_len(n)],
    resid = resid_vec[seq_len(n)],
    group_combo = group_combo[seq_len(n)],
    response = response_vec[seq_len(n)],
    stringsAsFactors = FALSE
  )
}

save_gaussian_diag_plots <- function(model, data, response_name, prefix) {
  resid_df <- tryCatch(
    get_residual_df(model, data = data, response_name = response_name),
    error = function(e) e
  )
  
  if (inherits(resid_df, "error")) {
    return(list(
      resid_vs_fitted = NULL,
      qq = NULL,
      note = resid_df$message
    ))
  }
  
  if (!is.data.frame(resid_df) || nrow(resid_df) == 0) {
    return(list(
      resid_vs_fitted = NULL,
      qq = NULL,
      note = "Residual data frame had zero rows."
    ))
  }
  
  # ----------------------------
  # Residuals vs fitted
  # ----------------------------
  p1 <- ggplot(resid_df, aes(x = fitted, y = resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste(pretty_indicator_name(response_name), "- Residuals vs fitted"),
      x = "Fitted values",
      y = "Residuals"
    )
  
  f1 <- file.path(DIR_PLOTS, paste0(prefix, "_resid_vs_fitted.png"))
  ggsave(f1, p1, width = 6.5, height = 4.5, dpi = 300)
  
  # ----------------------------
  # Q-Q plot with proper 95% envelope
  # ----------------------------
  qq_vals <- as.numeric(resid_df$resid)
  qq_vals <- qq_vals[is.finite(qq_vals)]
  
  if (length(qq_vals) < 3) {
    return(list(
      resid_vs_fitted = f1,
      qq = NULL,
      note = "Too few finite residuals available for diagnostic plots."
    ))
  }
  
  # Standardize residuals
  qq_vals <- as.numeric(scale(qq_vals))
  qq_vals <- sort(qq_vals)
  n <- length(qq_vals)
  
  # Theoretical quantiles
  theo <- stats::qnorm(stats::ppoints(n))
  
  # Simulated pointwise 95% envelope
  nsim <- 2000L
  sim_mat <- replicate(nsim, sort(stats::rnorm(n)))
  env_lo <- apply(sim_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
  env_hi <- apply(sim_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  
  f2 <- file.path(DIR_PLOTS, paste0(prefix, "_qq.png"))
  png(f2, width = 1800, height = 1400, res = 220)
  
  plot(
    theo, qq_vals,
    pch = 16,
    main = paste(pretty_indicator_name(response_name), "- Normal Q-Q with 95% envelope"),
    xlab = "Theoretical quantiles",
    ylab = "Standardized residual quantiles"
  )
  
  polygon(
    x = c(theo, rev(theo)),
    y = c(env_lo, rev(env_hi)),
    col = grDevices::adjustcolor("grey70", alpha.f = 0.35),
    border = NA
  )
  
  abline(0, 1, lty = 2, lwd = 2)
  points(theo, qq_vals, pch = 16)
  
  dev.off()
  
  list(
    resid_vs_fitted = f1,
    qq = f2,
    note = ""
  )
}

run_shapiro_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(x) > 5000) {
    return(list(statistic = NA_real_, p.value = NA_real_, note = "Shapiro not run: residual count outside 3..5000"))
  }
  out <- stats::shapiro.test(x)
  list(statistic = unname(out$statistic), p.value = out$p.value, note = "")
}

run_levene_on_residuals <- function(resid_df) {
  out <- tryCatch(car::leveneTest(resid ~ group_combo, data = resid_df), error = function(e) e)
  
  if (inherits(out, "error")) {
    return(list(
      statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p.value = NA_real_, note = out$message
    ))
  }
  
  out_df <- as.data.frame(out)
  stat_col <- grep("F value|F", names(out_df), value = TRUE)[1]
  p_col <- grep("^Pr\\(", names(out_df), value = TRUE)[1]
  
  list(
    statistic = if (!is.na(stat_col)) out_df[1, stat_col] else NA_real_,
    df1 = if ("Df" %in% names(out_df)) out_df[1, "Df"] else NA_real_,
    df2 = if ("Df" %in% names(out_df) && nrow(out_df) >= 2) out_df[2, "Df"] else NA_real_,
    p.value = if (!is.na(p_col)) out_df[1, p_col] else NA_real_,
    note = ""
  )
}

save_glmm_diag_plots <- function(model, prefix) {
  sim <- DHARMa::simulateResiduals(model, plot = FALSE)
  
  f1 <- file.path(DIR_PLOTS, paste0(prefix, "_DHARMa_overview.png"))
  png(f1, width = 1800, height = 1400, res = 220)
  plot(sim)
  dev.off()
  
  list(dharma = f1)
}

############################
# 7) Transformation branch
############################

allowed_transformations <- function(x) {
  out <- list(identity = function(z) z)
  if (ALLOW_LOG_TRANSFORM && all(is.finite(x)) && all(x > 0, na.rm = TRUE)) {
    out$log <- function(z) log(z)
  }
  if (ALLOW_SQRT_TRANSFORM && all(is.finite(x) | is.na(x)) && all(x >= 0, na.rm = TRUE)) {
    out$sqrt <- function(z) sqrt(z)
  }
  out
}

evaluate_gaussian_assumptions <- function(data, response_name, transformed_values = NULL,
                                          transformed_label = "raw", save_plots = TRUE) {
  dat <- data
  dat$.resp_eval <- if (is.null(transformed_values)) dat[[response_name]] else transformed_values
  
  rank_obj <- rank_candidate_models(
    data = dat,
    response = response_name,
    model_kind = "gaussian",
    transformed_name = ".resp_eval",
    REML = FALSE,
    random_candidates = ASSUMPTION_RANDOM_CANDIDATES,
    auto_drop_singular = TRUE
  )
  
  best <- select_best_model(rank_obj, allow_singular = FALSE)
  if (is.null(best)) {
    return(list(
      ok = FALSE,
      label = transformed_label,
      selected_fit = NULL,
      candidate_table = rank_obj$table,
      shapiro = list(statistic = NA_real_, p.value = NA_real_, note = "No Gaussian assumption model converged"),
      levene = list(statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p.value = NA_real_, note = "No Gaussian assumption model converged"),
      diag_paths = NULL
    ))
  }
  
  resid_df <- get_residual_df(best$model, dat, ".resp_eval")
  shp <- run_shapiro_safe(resid_df$resid)
  lev <- run_levene_on_residuals(resid_df)
  
  diag_paths <- if (save_plots) {
    save_gaussian_diag_plots(
      model = best$model,
      data = dat,
      response_name = ".resp_eval",
      prefix = paste0(safe_name(response_name), "_", safe_name(transformed_label))
    )
  } else {
    NULL
  }
  
  list(
    ok = TRUE,
    label = transformed_label,
    selected_fit = best,
    candidate_table = rank_obj$table,
    shapiro = shp,
    levene = lev,
    diag_paths = diag_paths
  )
}

choose_transformation <- function(data, response_name) {
  x <- data[[response_name]]
  tfuns <- allowed_transformations(x)
  
  raw_eval <- evaluate_gaussian_assumptions(
    data, response_name,
    transformed_values = NULL,
    transformed_label = "raw",
    save_plots = TRUE
  )
  
  if (isTRUE(raw_eval$ok) &&
      !is.na(raw_eval$shapiro$p.value) && raw_eval$shapiro$p.value >= P_VALUE_CUTOFF &&
      !is.na(raw_eval$levene$p.value)  && raw_eval$levene$p.value  >= P_VALUE_CUTOFF) {
    return(list(
      decision_path = "A",
      decision_label = "No transformation applied",
      chosen_transform = "identity",
      transformed_values = data[[response_name]],
      raw_eval = raw_eval,
      transformed_eval = raw_eval,
      decision_note = "Raw residual diagnostics acceptable"
    ))
  }
  
  tf_names <- setdiff(names(tfuns), "identity")
  tf_evals <- list()
  
  if (length(tf_names) > 0) {
    for (nm in tf_names) {
      vals <- tfuns[[nm]](data[[response_name]])
      tf_evals[[nm]] <- evaluate_gaussian_assumptions(
        data = data,
        response_name = response_name,
        transformed_values = vals,
        transformed_label = nm,
        save_plots = FALSE
      )
      tf_evals[[nm]]$transformed_values <- vals
    }
  }
  
  if (length(tf_evals) == 0) {
    return(list(
      decision_path = "A",
      decision_label = "No transformation applied",
      chosen_transform = "identity",
      transformed_values = data[[response_name]],
      raw_eval = raw_eval,
      transformed_eval = raw_eval,
      decision_note = "No valid candidate transformation for observed support"
    ))
  }
  
  tf_summary <- do.call(rbind, lapply(names(tf_evals), function(nm) {
    ev <- tf_evals[[nm]]
    data.frame(
      transform = nm,
      ok = ev$ok,
      shapiro_p = ev$shapiro$p.value,
      levene_p = ev$levene$p.value,
      aicc = if (!is.null(ev$selected_fit)) ev$selected_fit$aicc else NA_real_,
      pass_both = ifelse(
        !is.na(ev$shapiro$p.value) & !is.na(ev$levene$p.value),
        ev$shapiro$p.value >= P_VALUE_CUTOFF & ev$levene$p.value >= P_VALUE_CUTOFF,
        FALSE
      ),
      stringsAsFactors = FALSE
    )
  }))
  
  passing <- tf_summary[tf_summary$pass_both, , drop = FALSE]
  
  if (nrow(passing) > 0) {
    passing <- passing[order(passing$aicc, na.last = TRUE), , drop = FALSE]
    best_tf_name <- passing$transform[1]
    best_tf <- tf_evals[[best_tf_name]]
    
    # Save plots only for chosen transformed model
    best_tf$diag_paths <- save_gaussian_diag_plots(
      model = best_tf$selected_fit$model,
      data = transform(data, .resp_eval = best_tf$transformed_values),
      response_name = ".resp_eval",
      prefix = paste0(safe_name(response_name), "_", safe_name(best_tf_name))
    )
    
    return(list(
      decision_path = "B",
      decision_label = "Transformation applied",
      chosen_transform = best_tf_name,
      transformed_values = best_tf$transformed_values,
      raw_eval = raw_eval,
      transformed_eval = best_tf,
      decision_note = "Transformation chosen because transformed residual diagnostics passed; selected by lowest AICc among acceptable transformations"
    ))
  }
  
  tf_summary$joint_score <- pmin(tf_summary$shapiro_p, tf_summary$levene_p, na.rm = TRUE)
  tf_summary <- tf_summary[order(-tf_summary$joint_score, tf_summary$aicc), , drop = FALSE]
  best_tf_name <- tf_summary$transform[1]
  best_tf <- tf_evals[[best_tf_name]]
  
  best_tf$diag_paths <- save_gaussian_diag_plots(
    model = best_tf$selected_fit$model,
    data = transform(data, .resp_eval = best_tf$transformed_values),
    response_name = ".resp_eval",
    prefix = paste0(safe_name(response_name), "_", safe_name(best_tf_name))
  )
  
  list(
    decision_path = "B",
    decision_label = "Transformation applied",
    chosen_transform = best_tf_name,
    transformed_values = best_tf$transformed_values,
    raw_eval = raw_eval,
    transformed_eval = best_tf,
    decision_note = "No candidate passed all residual checks; chose the best available transformation by residual test score, then AICc"
  )
}

############################
# 8) ANOVA / post hoc / CLD helpers
############################

extract_car_anova_table <- function(model) {
  # Use Kenward-Roger F-tests for LMMs when lmerTest + pbkrtest are available;
  # fall back to Wald chi-square otherwise.
  use_kr <- inherits(model, "lmerMod") &&
    requireNamespace("lmerTest", quietly = TRUE) &&
    requireNamespace("pbkrtest", quietly = TRUE)

  tab <- if (use_kr) {
    tryCatch(
      as.data.frame(car::Anova(lmerTest::as.lmerModLmerTest(model), type = 3,
                               test.statistic = "F")),
      error = function(e) as.data.frame(car::Anova(model, type = 3, test.statistic = "Chisq"))
    )
  } else {
    as.data.frame(car::Anova(model, type = 3, test.statistic = "Chisq"))
  }

  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  stat_col <- grep("(^F$|^F value$|^Chisq$|Chisq|LR Chisq|Wald)", names(tab), value = TRUE)[1]
  p_col    <- grep("^Pr\\(", names(tab), value = TRUE)[1]
  
  df_cols <- grep("Df|df|Res\\.Df|Res Df|NumDF|DenDF|Den Df|Num Df", names(tab), value = TRUE)
  df_cols <- setdiff(df_cols, p_col)
  
  if (is.na(stat_col)) stop("Could not detect chi-square statistic column in car::Anova() output.")
  if (is.na(p_col)) stop("Could not detect p-value column in car::Anova() output.")
  
  out <- tab
  # Add BH-adjusted p-values across terms within this model (non-primary; for reference)
  raw_pv <- suppressWarnings(as.numeric(out[[p_col]]))
  out$p_value_BH <- p.adjust(raw_pv, method = "BH")
  out$.stat_col <- stat_col
  out$.p_col <- p_col
  attr(out, "stat_col") <- stat_col
  attr(out, "p_col") <- p_col
  attr(out, "df_cols") <- df_cols
  out
}

build_anova_string <- function(anova_df, alpha = 0.05) {
  stat_col <- attr(anova_df, "stat_col")
  p_col    <- attr(anova_df, "p_col")
  df_cols  <- attr(anova_df, "df_cols")
  
  keep <- !is.na(anova_df[[p_col]]) &
    anova_df[[p_col]] < alpha &
    !is.na(anova_df$term) &
    trimws(as.character(anova_df$term)) != "(Intercept)"
  
  x <- anova_df[keep, , drop = FALSE]
  
  if (nrow(x) == 0) return("")
  
  pieces <- apply(x, 1, function(row_i) {
    row_i <- as.list(row_i)
    term_i <- as.character(row_i$term)
    term_i <- gsub(":", " × ", term_i, fixed = TRUE)
    stat_i <- as.numeric(row_i[[stat_col]])
    
    df_txt <- ""
    if (length(df_cols) > 0) {
      df_vals <- unlist(row_i[df_cols])
      df_vals <- df_vals[!is.na(df_vals)]
      if (length(df_vals) > 0) df_txt <- paste(as.character(df_vals), collapse = ", ")
    }
    if (!nzchar(df_txt)) df_txt <- "NA"
    
    p_i <- as.numeric(row_i[[p_col]])
    stat_fmt <- trimws(sub("0+$", "", sub("\\.$", "",
      format(signif(stat_i, 4), scientific = FALSE, trim = TRUE))))
    p_fmt  <- fmt_p(p_i)
    p_part <- if (!is.na(p_fmt) && startsWith(p_fmt, "<")) paste0("p ", p_fmt) else paste0("p = ", p_fmt)
    paste0(term_i, ": \u03C7\u00B2(", df_txt, ") = ", stat_fmt, ", ", p_part)
  })
  
  paste(pieces, collapse = "\n")
}

get_posthoc_terms_layout_i <- function(anova_df, alpha = 0.05) {
  # Always include full 3-way interaction; also include any additional sig terms.
  # (Does NOT early-return on empty sig_terms — full interaction is mandatory.)
  p_col <- attr(anova_df, "p_col")
  sig_terms <- anova_df$term[!is.na(anova_df[[p_col]]) & anova_df[[p_col]] < alpha]
  sig_terms <- setdiff(sig_terms, "(Intercept)")
  unique(c("timepoint:Treatment:Group", sig_terms))
}

get_posthoc_terms_layout_ii <- function(anova_df, alpha = 0.05) {
  # Always include full 3-way interaction; also include any additional sig terms.
  # (Does NOT early-return on empty sig_terms — full interaction is mandatory.)
  p_col <- attr(anova_df, "p_col")
  sig_terms <- anova_df$term[!is.na(anova_df[[p_col]]) & anova_df[[p_col]] < alpha]
  sig_terms <- setdiff(sig_terms, "(Intercept)")
  unique(c("timepoint:Treatment:Group", sig_terms))
}

emm_spec_from_term <- function(term, model = NULL) {
  vars_present <- NULL
  if (!is.null(model)) {
    vars_present <- all.vars(stats::formula(model))
  }
  
  if (!is.null(model) && !("timepoint" %in% vars_present)) {
    switch(term,
           "Treatment" = ~ Treatment,
           "Group" = ~ Group,
           "Treatment:Group" = ~ Treatment:Group,
           "timepoint:Treatment:Group" = ~ Treatment:Group,
           "timepoint:Treatment" = ~ Treatment,
           "timepoint:Group" = ~ Group,
           NULL
    )
  } else {
    switch(term,
           "timepoint" = ~ timepoint,
           "Treatment" = ~ Treatment,
           "Group" = ~ Group,
           "timepoint:Treatment" = ~ timepoint:Treatment,
           "timepoint:Group" = ~ timepoint:Group,
           "Treatment:Group" = ~ Treatment:Group,
           "timepoint:Treatment:Group" = ~ timepoint:Treatment:Group,
           NULL
    )
  }
}

run_posthoc_for_terms <- function(model, terms, type = "response") {
  out <- list()
  
  for (term in unique(terms)) {
    spec <- emm_spec_from_term(term, model = model)
    if (is.null(spec)) next
    
    emm <- tryCatch(
      emmeans::emmeans(model, specs = spec, type = type),
      error = function(e) e
    )
    if (inherits(emm, "error")) {
      out[[term]] <- list(error = emm$message, table = NULL)
      next
    }
    
    cmp <- tryCatch(
      as.data.frame(summary(emmeans::contrast(emm, method = "pairwise", adjust = "tukey"))),
      error = function(e) e
    )
    if (inherits(cmp, "error")) {
      out[[term]] <- list(error = cmp$message, table = NULL)
      next
    }
    
    p_col <- grep("^p\\.value$|^p.value$|p.value", names(cmp), value = TRUE)[1]
    
    out[[term]] <- list(
      error = NULL,
      table = cmp,
      p_col = ifelse(is.na(p_col), NULL, p_col)
    )
  }
  
  out
}

run_cld_all_combinations <- function(model, type = "response", alpha = 0.05, max_letters = 3) {
  cld_fallback_df <- function(msg) {
    data.frame(
      Group = NA_character_,
      Treatment = NA_character_,
      timepoint = NA_character_,
      cld_letter = NA_character_,
      note = as.character(msg),
      stringsAsFactors = FALSE
    )
  }
  
  trim_letters <- function(x, max_letters = 3) {
    x <- gsub("\\s+", "", as.character(x))
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    chars <- chars[nzchar(chars)]
    if (length(chars) <= max_letters) {
      paste(chars, collapse = "")
    } else {
      paste(chars[seq_len(max_letters)], collapse = "")
    }
  }
  
  vars_present <- all.vars(stats::formula(model))
  has_timepoint <- "timepoint" %in% vars_present
  
  spec <- if (has_timepoint) {
    ~ timepoint * Group * Treatment
  } else {
    ~ Group * Treatment
  }
  
  term_name <- if (has_timepoint) {
    "timepoint:Treatment:Group"
  } else {
    "Treatment:Group"
  }
  
  emm <- tryCatch(
    emmeans::emmeans(model, specs = spec, type = type),
    error = function(e) e
  )
  
  if (inherits(emm, "error")) {
    return(setNames(
      list(cld_fallback_df(paste("CLD failed:", emm$message))),
      term_name
    ))
  }
  
  cld_obj <- tryCatch(
    multcomp::cld(
      emm,
      adjust = "tukey",
      alpha = alpha,
      Letters = c(letters, LETTERS, 0:9)
    ),
    error = function(e) e
  )
  
  if (inherits(cld_obj, "error")) {
    return(setNames(
      list(cld_fallback_df(paste("CLD failed:", cld_obj$message))),
      term_name
    ))
  }
  
  cld_df <- tryCatch(as.data.frame(cld_obj), error = function(e) e)
  
  if (inherits(cld_df, "error")) {
    return(setNames(
      list(cld_fallback_df(paste("CLD failed:", cld_df$message))),
      term_name
    ))
  }
  
  if (!".group" %in% names(cld_df)) {
    return(setNames(
      list(cld_fallback_df("CLD output did not contain .group column.")),
      term_name
    ))
  }
  
  if (!"Group" %in% names(cld_df)) cld_df$Group <- NA_character_
  if (!"Treatment" %in% names(cld_df)) cld_df$Treatment <- NA_character_
  if (!"timepoint" %in% names(cld_df)) cld_df$timepoint <- NA_character_
  
  out <- data.frame(
    timepoint = as.character(cld_df$timepoint),
    Group = as.character(cld_df$Group),
    Treatment = as.character(cld_df$Treatment),
    cld_letter = vapply(cld_df$.group, trim_letters, character(1), max_letters = max_letters),
    note = "",
    stringsAsFactors = FALSE
  )
  
  out$term <- term_name
  attr(out, "warning_note") <- "CLD is based on Tukey-adjusted estimated marginal means; displayed letters are truncated."
  
  setNames(list(out), term_name)
}

run_cld_for_terms <- function(model, terms, type = "response", alpha = 0.05, max_letters = 3) {
  out <- list()
  
  trim_letters <- function(x, max_letters = 3) {
    x <- gsub("\\s+", "", as.character(x))
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    chars <- chars[nzchar(chars)]
    if (length(chars) <= max_letters) paste(chars, collapse = "") else paste(chars[seq_len(max_letters)], collapse = "")
  }
  
  vars_present <- all.vars(stats::formula(model))
  has_timepoint <- "timepoint" %in% vars_present
  
  for (term in unique(terms)) {
    spec <- emm_spec_from_term(term, model = model)
    if (is.null(spec)) next
    
    emm <- tryCatch(
      emmeans::emmeans(model, specs = spec, type = type),
      error = function(e) e
    )
    if (inherits(emm, "error")) {
      out[[term]] <- data.frame(
        note = paste("CLD failed:", emm$message),
        stringsAsFactors = FALSE
      )
      next
    }
    
    cld_obj <- tryCatch(
      multcomp::cld(
        emm,
        adjust = "tukey",
        alpha = alpha,
        Letters = c(letters, LETTERS, 0:9)
      ),
      error = function(e) e
    )
    if (inherits(cld_obj, "error")) {
      out[[term]] <- data.frame(
        note = paste("CLD failed:", cld_obj$message),
        stringsAsFactors = FALSE
      )
      next
    }
    
    cld_df <- tryCatch(as.data.frame(cld_obj), error = function(e) e)
    if (inherits(cld_df, "error")) {
      out[[term]] <- data.frame(
        note = paste("CLD failed:", cld_df$message),
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (!".group" %in% names(cld_df)) {
      out[[term]] <- data.frame(
        note = "CLD output did not contain .group column.",
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (!"Group" %in% names(cld_df)) cld_df$Group <- NA_character_
    if (!"Treatment" %in% names(cld_df)) cld_df$Treatment <- NA_character_
    if (!"timepoint" %in% names(cld_df)) cld_df$timepoint <- NA_character_
    
    keep_cols <- c("timepoint", "Group", "Treatment")
    if (!has_timepoint) keep_cols <- setdiff(keep_cols, "timepoint")
    
    out_df <- cld_df[, intersect(keep_cols, names(cld_df)), drop = FALSE]
    out_df$cld_letter <- vapply(cld_df$.group, trim_letters, character(1), max_letters = max_letters)
    out_df$note <- ""
    
    out[[term]] <- out_df
  }
  
  out
}

############################
# 9) Per-indicator analysis branches
############################

placeholder_result <- function(indicator, note) {
  res <- list(
    indicator = indicator,
    status = "skipped",
    note = note,
    candidate_table = data.frame(note = note, stringsAsFactors = FALSE),
    selected_model = NULL,
    selected_formula = NA_character_,
    selected_aicc = NA_real_,
    anova_table = data.frame(note = note, stringsAsFactors = FALSE),
    anova_string = "",
    posthoc = list(),
    cld = empty_cld_df(),
    cld_all_combinations = list(),
    cld_tables = list(),
    diagnostics = list()
  )
  attr(res, "is_placeholder") <- TRUE
  res
}

empty_cld_df <- function() {
  data.frame(
    timepoint = character(),
    Group = character(),
    Treatment = character(),
    cld_letter = character(),
    term = character(),
    stringsAsFactors = FALSE
  )
}

analyze_beta_indicator <- function(data, indicator, layout_label) {
  if (!requireNamespace("glmmTMB", quietly = TRUE))
    return(placeholder_result(indicator, "glmmTMB not available for Beta GLMM"))

  # Convert seconds to proportion; clamp to (0, 1) to avoid Beta boundary issues.
  raw_s    <- data[[indicator]]
  prop_col <- ".prop_beta"
  data[[prop_col]] <- pmax(pmin(raw_s / RECORDING_DURATION_SECONDS, 1 - 1e-4), 1e-4)
  data     <- data[is.finite(data[[prop_col]]), , drop = FALSE]
  if (nrow(data) < 5) return(placeholder_result(indicator, "Too few rows for Beta GLMM"))

  fixed_formula <- build_fixed_formula(data)
  fits <- lapply(RANDOM_CANDIDATES, function(re_term) {
    frm_str <- paste(prop_col, "~", fixed_formula, "+", re_term)
    fit <- tryCatch(
      glmmTMB::glmmTMB(
        stats::as.formula(frm_str), data = data,
        family = glmmTMB::beta_family(link = "logit"),
        control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000))
      ),
      error = function(e) NULL
    )
    list(formula = frm_str, model = fit,
         ok = !is.null(fit),
         aicc = if (!is.null(fit)) tryCatch(AICcmodavg::AICc(fit), error = function(e) {
           k <- length(stats::coef(fit)); n <- nrow(data)
           aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
           if (is.finite(aic) && n > k + 1) aic + 2*k*(k+1)/(n-k-1) else NA_real_
         }) else NA_real_)
  })

  ok_fits <- Filter(function(x) isTRUE(x$ok) && is.finite(x$aicc), fits)
  if (length(ok_fits) == 0) {
    # Fallback 1: zero-inflated Beta GLMM (handles boundary-heavy distributions)
    zi_fits <- lapply(RANDOM_CANDIDATES, function(re_term) {
      frm_str <- paste(prop_col, "~", fixed_formula, "+", re_term)
      fit <- tryCatch(
        glmmTMB::glmmTMB(
          stats::as.formula(frm_str), data = data,
          family = glmmTMB::beta_family(link = "logit"),
          ziformula = ~1,
          control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 2000))
        ),
        error = function(e) NULL
      )
      list(formula = frm_str, model = fit,
           ok = !is.null(fit),
           aicc = if (!is.null(fit)) tryCatch(AICcmodavg::AICc(fit), error = function(e) {
             k <- length(stats::coef(fit)); n <- nrow(data)
             aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
             if (is.finite(aic) && n > k + 1L) aic + 2*k*(k+1)/(n-k-1) else NA_real_
           }) else NA_real_)
    })
    zi_ok <- Filter(function(x) isTRUE(x$ok) && is.finite(x$aicc), zi_fits)
    if (length(zi_ok) > 0L) {
      fits    <- zi_fits
      ok_fits <- zi_ok
    } else {
      # Fallback 2: logit-transform + Gaussian LMM
      logit_y <- log(data[[prop_col]] / (1 - data[[prop_col]]))
      data[[".logit_y"]] <- logit_y
      data_logit <- data[is.finite(data[[".logit_y"]]), , drop = FALSE]
      gauss_fits <- if (nrow(data_logit) >= 5L) lapply(RANDOM_CANDIDATES, function(re_term) {
        frm_str <- paste(".logit_y ~", fixed_formula, "+", re_term)
        fit <- tryCatch(
          lme4::lmer(stats::as.formula(frm_str), data = data_logit, REML = FALSE),
          error = function(e) NULL
        )
        aicc_val <- if (!is.null(fit))
          tryCatch(AICcmodavg::AICc(fit), error = function(e) NA_real_) else NA_real_
        list(formula = frm_str, model = fit,
             ok = !is.null(fit) && is.finite(aicc_val), aicc = aicc_val)
      }) else list()
      gauss_ok <- Filter(function(x) isTRUE(x$ok) && is.finite(x$aicc), gauss_fits)
      if (length(gauss_ok) > 0L) {
        fits     <- gauss_fits
        ok_fits  <- gauss_ok
        data     <- data_logit
        prop_col <- ".logit_y"
      } else {
        return(placeholder_result(indicator,
          "No converged Beta GLMM, ZI-Beta GLMM, or logit-Gaussian model"))
      }
    }
  }

  best_fit <- ok_fits[[which.min(vapply(ok_fits, `[[`, numeric(1), "aicc"))]]
  final_model <- best_fit$model

  anova_df <- tryCatch(
    extract_car_anova_table(final_model),
    error = function(e) data.frame()
  )
  anova_string  <- if (nrow(anova_df) > 0) build_anova_caption_multiline(anova_df) else ""
  has_sig       <- nrow(anova_df) > 0 && has_non_intercept_sig_terms(anova_df)
  posthoc_terms <- if (layout_label == "all_groups")
    get_posthoc_terms_layout_i(anova_df) else get_posthoc_terms_layout_ii(anova_df)

  posthoc_list <- if (RUN_POSTHOC_ONLY_IF_SIG && !has_sig) list() else
    run_posthoc_for_terms(final_model, posthoc_terms, type = "response")
  cld_all <- if (RUN_CLD_ONLY_IF_SIG && !has_sig) list() else
    run_cld_all_combinations(final_model, type = "response")
  cld_tables <- if (RUN_CLD_ONLY_IF_SIG && !has_sig) list() else
    run_cld_for_terms(final_model, posthoc_terms, type = "response")
  cld_df <- if (length(cld_all) > 0) cld_all[[1]] else empty_cld_df()

  # DHARMa diagnostics
  dharma_summ <- NULL
  if (requireNamespace("DHARMa", quietly = TRUE)) {
    sim <- tryCatch(DHARMa::simulateResiduals(final_model, n = 500L, plot = FALSE),
                   error = function(e) NULL)
    if (!is.null(sim)) {
      png_path <- file.path(DIR_PLOTS,
                            paste0("dharma_", safe_name(indicator), "_", layout_label, ".png"))
      grDevices::png(png_path, width = 1600, height = 800, res = 150)
      try(DHARMa::plot(sim, title = indicator), silent = TRUE)
      grDevices::dev.off()
      dharma_summ <- tryCatch(
        data.frame(
          outcome         = indicator,
          KS_p            = DHARMa::testUniformity(sim, plot = FALSE)$p.value,
          dispersion_p    = DHARMa::testDispersion(sim, plot = FALSE)$p.value,
          stringsAsFactors = FALSE
        ),
        error = function(e) NULL
      )
    }
  }

  list(
    indicator = indicator, status = "ok",
    note = paste0("Beta GLMM on prop_in_flow = ", indicator, " / ", RECORDING_DURATION_SECONDS),
    transformation = list(chosen_transform = "beta_glmm", decision_path = "Beta GLMM for bounded [0,1]"),
    candidate_table = do.call(rbind, lapply(fits, function(x)
      data.frame(formula = x$formula, ok = isTRUE(x$ok), aicc = x$aicc, stringsAsFactors = FALSE))),
    selected_model = final_model,
    selected_formula = best_fit$formula,
    selected_aicc = best_fit$aicc,
    anova_table = anova_df, anova_string = anova_string,
    posthoc = posthoc_list, cld = cld_df,
    cld_all_combinations = cld_all, cld_tables = cld_tables,
    dharma = dharma_summ,
    diagnostics = list()
  )
}

analyze_continuous_indicator <- function(data, indicator, layout_label, tf = NULL) {

  # ---- Pre-modelling Gamma diagnostic check (Task 1) -------------------------
  # run_diag_for_indicator() tests a simple intercept + RE model for SW/Levene.
  # Only run when: tf not pre-supplied, module is loaded, and all y > 0.
  .s5_diag_row  <- NULL   # assumption-diagnostics table row (stored in $diag_row)
  .s5_use_gamma <- FALSE  # TRUE ⟹ skip Gaussian path and use Gamma LMM directly

  if (is.null(tf) && exists("run_diag_for_indicator", mode = "function")) {
    .s5_fin     <- is.finite(data[[indicator]])
    .s5_all_pos <- sum(.s5_fin) > 0L && all(data[[indicator]][.s5_fin] > 0)
    if (.s5_all_pos) {
      .s5_dat_fin <- data[.s5_fin, , drop = FALSE]
      .s5_grp <- tryCatch(
        interaction(.s5_dat_fin$Group, .s5_dat_fin$Treatment, drop = TRUE),
        error = function(e) rep("all", nrow(.s5_dat_fin))
      )
      .s5_diag_res <- tryCatch(
        run_diag_for_indicator(
          y             = data[[indicator]][.s5_fin],
          data          = .s5_dat_fin,
          grouping      = .s5_grp,
          re_form_str   = "(1 | trial_id)",
          ind_nm        = indicator,
          analysis_unit = "fish",
          try_gamma     = TRUE
        ),
        error = function(e) {
          warning("[DIAG S5] Gamma check failed for ", indicator, ": ",
                  e$message, call. = FALSE)
          NULL
        }
      )
      if (!is.null(.s5_diag_res) && isTRUE(.s5_diag_res$use_gamma)) {
        .s5_use_gamma <- TRUE
        .s5_diag_row  <- .s5_diag_res$row
      }
    }
  }

  # ---- Gamma path (diagnostic-driven, short-circuits Gaussian fitting) -------
  if (.s5_use_gamma && requireNamespace("glmmTMB", quietly = TRUE)) {
    raw_vals <- data[[indicator]]
    .shift   <- max(0, -min(raw_vals, na.rm = TRUE) + 1e-3)
    data[[".gamma_y"]] <- raw_vals + .shift
    data[[".gamma_y"]][!is.finite(data[[".gamma_y"]]) | data[[".gamma_y"]] <= 0] <- NA_real_
    data_g <- data[is.finite(data[[".gamma_y"]]), , drop = FALSE]
    if (nrow(data_g) >= 5L) {
      fml_g_fixed  <- build_fixed_formula(data_g)
      shift_note_d <- if (.shift > 0) paste0("; shifted +", round(.shift, 4)) else ""
      gamma_fits_d <- lapply(RANDOM_CANDIDATES, function(re_term) {
        frm_str <- paste(".gamma_y ~", fml_g_fixed, "+", re_term)
        fit <- tryCatch(
          glmmTMB::glmmTMB(
            stats::as.formula(frm_str), data = data_g,
            family  = Gamma(link = "log"),
            control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000))
          ),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        aicc_val <- tryCatch(AICcmodavg::AICc(fit), error = function(e) {
          k   <- length(stats::coef(fit)); n <- nrow(data_g)
          aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
          if (is.finite(aic) && n > k + 1L) aic + 2*k*(k+1)/(n-k-1) else NA_real_
        })
        list(formula = frm_str, model = fit, aicc = aicc_val)
      })
      gamma_ok_d <- Filter(function(x) !is.null(x) && is.finite(x$aicc), gamma_fits_d)
      if (length(gamma_ok_d) > 0L) {
        best_gd      <- gamma_ok_d[[which.min(vapply(gamma_ok_d, `[[`, numeric(1L), "aicc"))]]
        final_model  <- best_gd$model
        anova_df     <- tryCatch(extract_car_anova_table(final_model),
                                 error = function(e) data.frame())
        anova_string <- if (nrow(anova_df) > 0) build_anova_caption_multiline(anova_df) else ""
        has_sig      <- has_non_intercept_sig_terms(anova_df)
        ph_terms     <- if (layout_label == "all_groups")
          get_posthoc_terms_layout_i(anova_df) else get_posthoc_terms_layout_ii(anova_df)
        ph_list      <- if (RUN_POSTHOC_ONLY_IF_SIG && !has_sig) list() else
          run_posthoc_for_terms(final_model, ph_terms, type = "response")
        cld_all      <- if (RUN_CLD_ONLY_IF_SIG && !has_sig) list() else
          run_cld_all_combinations(final_model, type = "response")
        cld_tabs     <- if (RUN_CLD_ONLY_IF_SIG && !has_sig) list() else
          run_cld_for_terms(final_model, ph_terms, type = "response")
        cld_df_gd    <- if (length(cld_all) > 0L) cld_all[[1L]] else empty_cld_df()
        dharma_gd    <- NULL
        if (requireNamespace("DHARMa", quietly = TRUE)) {
          sim_gd <- tryCatch(
            DHARMa::simulateResiduals(final_model, n = 500L, plot = FALSE),
            error = function(e) NULL)
          if (!is.null(sim_gd)) {
            png_gd <- file.path(DIR_PLOTS, paste0("dharma_", safe_name(indicator),
                                                   "_diag_gamma_", layout_label, ".png"))
            grDevices::png(png_gd, width = 1600, height = 800, res = 150)
            try(DHARMa::plot(sim_gd, title = indicator), silent = TRUE)
            grDevices::dev.off()
            dharma_gd <- tryCatch(data.frame(
              outcome      = indicator,
              KS_p         = DHARMa::testUniformity(sim_gd, plot = FALSE)$p.value,
              dispersion_p = DHARMa::testDispersion(sim_gd, plot = FALSE)$p.value,
              stringsAsFactors = FALSE), error = function(e) NULL)
          }
        }
        r2_gd   <- tryCatch(MuMIn::r.squaredGLMM(final_model), error = function(e) NULL)
        cand_gd <- do.call(rbind, lapply(gamma_fits_d, function(x)
          if (is.null(x))
            data.frame(formula = NA_character_, aicc = NA_real_, stringsAsFactors = FALSE)
          else
            data.frame(formula = x$formula, aicc = x$aicc, stringsAsFactors = FALSE)))
        return(list(
          indicator    = indicator, status = "ok",
          note         = paste0("Gamma LMM (diagnostic-driven", shift_note_d,
                                "): pre-modelling SW/Levene failed for all Gaussian candidates"),
          transformation = list(chosen_transform = "gamma_diagnostic",
                                decision_path    = "Gamma LMM (diagnostic)",
                                decision_label   = "Gamma LMM (log link)",
                                decision_note    = "Pre-modelling diagnostics: Gamma preferred over Gaussian"),
          candidate_table      = cand_gd,
          selected_model       = final_model,
          selected_formula     = best_gd$formula,
          selected_aicc        = best_gd$aicc,
          anova_table          = anova_df,
          anova_string         = anova_string,
          posthoc              = ph_list,
          cld                  = cld_df_gd,
          cld_all_combinations = cld_all,
          cld_tables           = cld_tabs,
          diagnostics          = list(),
          diag_row             = .s5_diag_row,
          dharma               = dharma_gd,
          r2m = if (!is.null(r2_gd) && "R2m" %in% colnames(r2_gd)) r2_gd[1, "R2m"] else NA_real_,
          r2c = if (!is.null(r2_gd) && "R2c" %in% colnames(r2_gd)) r2_gd[1, "R2c"] else NA_real_
        ))
      }
    }
    message("[DIAG S5] Diagnostic-driven Gamma fitting failed for ", indicator,
            " — falling back to Gaussian transformation selection.")
    .s5_use_gamma <- FALSE
  }

  # ---- Gaussian path (original) ----------------------------------------------
  if (is.null(tf)) {
    tf <- choose_transformation(data, indicator)
  }

  # Task 2: Extract assumption diagnostics from choose_transformation() result
  if (is.null(.s5_diag_row) && exists("diag_row_from_choose_tf", mode = "function")) {
    .s5_diag_row <- tryCatch(
      diag_row_from_choose_tf(tf, indicator, nrow(data), "fish"),
      error = function(e) NULL
    )
  }

  dat <- data
  dat$.analysis_y <- tf$transformed_values

  final_rank <- rank_candidate_models(
    data = dat,
    response = indicator,
    model_kind = "gaussian",
    transformed_name = ".analysis_y",
    REML = FALSE,
    random_candidates = RANDOM_CANDIDATES,
    auto_drop_singular = TRUE
  )

  best_ml <- select_best_model(final_rank, allow_singular = FALSE)

  if (is.null(best_ml)) {
    # Gamma GLMM fallback: shift data to strictly-positive range, refit with glmmTMB
    if (requireNamespace("glmmTMB", quietly = TRUE)) {
      raw_vals <- data[[indicator]]
      .shift   <- max(0, -min(raw_vals, na.rm = TRUE) + 1e-3)
      data[[".gamma_y"]] <- raw_vals + .shift
      data[[".gamma_y"]][!is.finite(data[[".gamma_y"]]) | data[[".gamma_y"]] <= 0] <- NA_real_
      data_g <- data[is.finite(data[[".gamma_y"]]), , drop = FALSE]
      if (nrow(data_g) >= 5L) {
        fml_g_fixed <- build_fixed_formula(data_g)
        gamma_fits  <- lapply(RANDOM_CANDIDATES, function(re_term) {
          frm_str <- paste(".gamma_y ~", fml_g_fixed, "+", re_term)
          fit <- tryCatch(
            glmmTMB::glmmTMB(
              stats::as.formula(frm_str), data = data_g,
              family = Gamma(link = "log"),
              control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000))
            ),
            error = function(e) NULL
          )
          if (is.null(fit)) return(NULL)
          aicc_val <- tryCatch(AICcmodavg::AICc(fit), error = function(e) {
            k   <- length(stats::coef(fit)); n <- nrow(data_g)
            aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
            if (is.finite(aic) && n > k + 1L) aic + 2*k*(k+1)/(n-k-1) else NA_real_
          })
          list(formula = frm_str, model = fit, aicc = aicc_val)
        })
        gamma_ok <- Filter(function(x) !is.null(x) && is.finite(x$aicc), gamma_fits)
        if (length(gamma_ok) > 0L) {
          best_g      <- gamma_ok[[which.min(vapply(gamma_ok, `[[`, numeric(1), "aicc"))]]
          final_model <- best_g$model
          shift_note  <- if (.shift > 0) paste0("; shifted +", round(.shift, 4)) else ""
          anova_df    <- tryCatch(extract_car_anova_table(final_model),
                                  error = function(e) data.frame())
          anova_string <- if (nrow(anova_df) > 0) build_anova_caption_multiline(anova_df) else ""
          has_sig     <- has_non_intercept_sig_terms(anova_df)
          ph_terms    <- if (layout_label == "all_groups")
            get_posthoc_terms_layout_i(anova_df) else get_posthoc_terms_layout_ii(anova_df)
          ph_list     <- if (RUN_POSTHOC_ONLY_IF_SIG && !has_sig) list() else
            run_posthoc_for_terms(final_model, ph_terms, type = "response")
          cld_all     <- if (RUN_CLD_ONLY_IF_SIG && !has_sig) list() else
            run_cld_all_combinations(final_model, type = "response")
          cld_tabs    <- if (RUN_CLD_ONLY_IF_SIG && !has_sig) list() else
            run_cld_for_terms(final_model, ph_terms, type = "response")
          cld_df_g    <- if (length(cld_all) > 0L) cld_all[[1L]] else empty_cld_df()
          dharma_g    <- NULL
          if (requireNamespace("DHARMa", quietly = TRUE)) {
            sim_g <- tryCatch(
              DHARMa::simulateResiduals(final_model, n = 500L, plot = FALSE),
              error = function(e) NULL)
            if (!is.null(sim_g)) {
              png_g <- file.path(DIR_PLOTS, paste0("dharma_", safe_name(indicator),
                                                   "_", layout_label, ".png"))
              grDevices::png(png_g, width = 1600, height = 800, res = 150)
              try(DHARMa::plot(sim_g, title = indicator), silent = TRUE)
              grDevices::dev.off()
              dharma_g <- tryCatch(data.frame(
                outcome      = indicator,
                KS_p         = DHARMa::testUniformity(sim_g, plot = FALSE)$p.value,
                dispersion_p = DHARMa::testDispersion(sim_g, plot = FALSE)$p.value,
                stringsAsFactors = FALSE), error = function(e) NULL)
            }
          }
          r2_g <- tryCatch(MuMIn::r.squaredGLMM(final_model), error = function(e) NULL)
          cand_g <- do.call(rbind, lapply(gamma_fits, function(x)
            if (is.null(x))
              data.frame(formula = NA_character_, aicc = NA_real_, stringsAsFactors = FALSE)
            else
              data.frame(formula = x$formula, aicc = x$aicc, stringsAsFactors = FALSE)))
          # Update diag_row to reflect the Gamma family actually used
          .s5_diag_row <- if (exists("build_diag_table_row", mode = "function")) {
            tryCatch(
              build_diag_table_row(
                "fish", "Gamma GLMM (log link)", nrow(data),
                "not applicable", "not applicable",
                "Gamma GLMM (log link)", "Gaussian LMMs all failed; Gamma GLMM fallback",
                "not applicable", "not applicable",
                indicator, ""
              ),
              error = function(e) .s5_diag_row
            )
          } else .s5_diag_row
          return(list(
            indicator = indicator, status = "ok",
            note = paste0("Gamma GLMM fallback (Gaussian LMMs all failed", shift_note, ")"),
            transformation = list(chosen_transform = "gamma_glmm_shifted",
                                  decision_path    = "Gamma GLMM fallback",
                                  decision_label   = paste0("Gamma GLMM (shifted",
                                                            shift_note, ")"),
                                  decision_note    = "Gaussian LMMs all failed; Gamma GLMM used"),
            candidate_table   = cand_g,
            selected_model    = final_model,
            selected_formula  = best_g$formula,
            selected_aicc     = best_g$aicc,
            anova_table       = anova_df,
            anova_string      = anova_string,
            posthoc           = ph_list,
            cld               = cld_df_g,
            cld_all_combinations = cld_all,
            cld_tables        = cld_tabs,
            diagnostics       = list(),
            diag_row          = .s5_diag_row,
            dharma            = dharma_g,
            r2m = if (!is.null(r2_g) && "R2m" %in% colnames(r2_g)) r2_g[1, "R2m"] else NA_real_,
            r2c = if (!is.null(r2_g) && "R2c" %in% colnames(r2_g)) r2_g[1, "R2c"] else NA_real_
          ))
        }
      }
    }
    return(placeholder_result(indicator, "No converged Gaussian or Gamma mixed model"))
  }

  final_model <- best_ml$model
  anova_df <- extract_car_anova_table(final_model)
  anova_string <- build_anova_caption_multiline(anova_df)
  has_sig_non_intercept <- has_non_intercept_sig_terms(anova_df)

  if (layout_label == "all_groups") {
    posthoc_terms <- get_posthoc_terms_layout_i(anova_df)
  } else {
    posthoc_terms <- get_posthoc_terms_layout_ii(anova_df)
  }
  
  posthoc_list <- if (RUN_POSTHOC_ONLY_IF_SIG && !has_sig_non_intercept) {
    list()
  } else {
    run_posthoc_for_terms(final_model, posthoc_terms, type = "response")
  }
  
  cld_all_combinations <- if (RUN_CLD_ONLY_IF_SIG && !has_sig_non_intercept) {
    list()
  } else {
    run_cld_all_combinations(final_model, type = "response")
  }
  
  cld_tables <- if (RUN_CLD_ONLY_IF_SIG && !has_sig_non_intercept) {
    list()
  } else {
    run_cld_for_terms(final_model, posthoc_terms, type = "response")
  }
  
  cld_df <- if ("timepoint:Treatment:Group" %in% names(cld_all_combinations)) {
    cld_all_combinations[["timepoint:Treatment:Group"]]
  } else if ("Treatment:Group" %in% names(cld_all_combinations)) {
    cld_all_combinations[["Treatment:Group"]]
  } else if (length(cld_all_combinations) > 0) {
    cld_all_combinations[[1]]
  } else {
    empty_cld_df()
  }
  
  diag_paths <- save_gaussian_diag_plots(
    model = final_model,
    data = dat,
    response_name = indicator,
    prefix = paste0(safe_name(indicator), "_layout_", layout_label, "_final")
  )

  # DHARMa simulated-residual diagnostics
  dharma_summ <- NULL
  if (requireNamespace("DHARMa", quietly = TRUE)) {
    sim <- tryCatch(DHARMa::simulateResiduals(final_model, n = 500L, plot = FALSE),
                   error = function(e) NULL)
    if (!is.null(sim)) {
      png_path <- file.path(DIR_PLOTS,
                            paste0("dharma_", safe_name(indicator), "_", layout_label, ".png"))
      grDevices::png(png_path, width = 1600, height = 800, res = 150)
      try(DHARMa::plot(sim, title = indicator), silent = TRUE)
      grDevices::dev.off()
      dharma_summ <- tryCatch(
        data.frame(
          outcome      = indicator,
          KS_p         = DHARMa::testUniformity(sim, plot = FALSE)$p.value,
          dispersion_p = DHARMa::testDispersion(sim, plot = FALSE)$p.value,
          stringsAsFactors = FALSE
        ),
        error = function(e) NULL
      )
    }
  }

  # R²m / R²c (MuMIn)
  r2_vals <- tryCatch(MuMIn::r.squaredGLMM(final_model), error = function(e) NULL)
  r2m <- if (!is.null(r2_vals) && "R2m" %in% colnames(r2_vals)) r2_vals[1, "R2m"] else NA_real_
  r2c <- if (!is.null(r2_vals) && "R2c" %in% colnames(r2_vals)) r2_vals[1, "R2c"] else NA_real_

  list(
    indicator = indicator,
    status = "ok",
    note = paste0(
      "Decision path ", tf$decision_path, ": ", tf$decision_label,
      "; chosen transform = ", tf$chosen_transform,
      "; ", tf$decision_note
    ),
    transformation = tf,
    candidate_table = final_rank$table,
    selected_model = final_model,
    selected_formula = best_ml$formula,
    selected_aicc = best_ml$aicc,
    anova_table = anova_df,
    anova_string = anova_string,
    posthoc = posthoc_list,
    cld = cld_df,
    cld_all_combinations = cld_all_combinations,
    cld_tables = cld_tables,
    diagnostics = diag_paths,
    diag_row    = .s5_diag_row,
    dharma = dharma_summ,
    r2m = r2m,
    r2c = r2c
  )
}

analyze_count_indicator <- function(data, indicator, layout_label) {
  rank_pois <- rank_candidate_models(
    data, indicator,
    model_kind = "poisson",
    REML = FALSE,
    random_candidates = RANDOM_CANDIDATES,
    auto_drop_singular = TRUE
  )

  # Overdispersion diagnostic: Pearson chi-sq / residual df for the best Poisson model.
  # Ratio > 1.5 suggests overdispersion (informational only; NB is always included by AICc).
  best_pois_fit <- tryCatch(
    select_best_model(rank_pois, allow_singular = FALSE, delta_prefer_nonsingular = 2)$model,
    error = function(e) NULL
  )
  overdisp_ratio <- if (!is.null(best_pois_fit) && inherits(best_pois_fit, "merMod")) {
    tryCatch({
      pr <- sum(residuals(best_pois_fit, type = "pearson")^2, na.rm = TRUE)
      df <- df.residual(best_pois_fit)
      if (is.finite(df) && df > 0) round(pr / df, 3) else NA_real_
    }, error = function(e) NA_real_)
  } else NA_real_

  rank_nb <- rank_candidate_models(
    data, indicator,
    model_kind = "nb",
    REML = FALSE,
    random_candidates = RANDOM_CANDIDATES,
    auto_drop_singular = TRUE
  )

  all_rows <- rbind(
    transform(rank_pois$table, family = "poisson"),
    transform(rank_nb$table, family = "negative_binomial")
  )
  # Attach overdispersion ratio as a column on Poisson rows for the saved table
  all_rows$overdisp_pearson_ratio <- ifelse(
    all_rows$family == "poisson", overdisp_ratio, NA_real_
  )
  all_rows <- all_rows[order(all_rows$aicc, na.last = TRUE), , drop = FALSE]
  
  all_fits <- c(rank_pois$fits, rank_nb$fits)
  
  fit_tbl <- data.frame(
    idx = seq_along(all_fits),
    ok = vapply(all_fits, function(x) isTRUE(x$ok), logical(1)),
    singular = vapply(all_fits, function(x) isTRUE(x$singular), logical(1)),
    aicc = vapply(all_fits, function(x) x$aicc, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  fit_tbl <- fit_tbl[fit_tbl$ok & is.finite(fit_tbl$aicc), , drop = FALSE]
  
  if (nrow(fit_tbl) == 0) {
    return(placeholder_result(indicator, "No count mixed model with finite AICc"))
  }
  
  fit_tbl <- fit_tbl[order(fit_tbl$aicc), , drop = FALSE]
  fit_tbl$delta_aicc <- fit_tbl$aicc - fit_tbl$aicc[1]
  
  best_row <- fit_tbl[1, , drop = FALSE]
  
  if (!isTRUE(best_row$singular)) {
    best <- all_fits[[best_row$idx]]
  } else {
    near_tbl <- fit_tbl[fit_tbl$delta_aicc < 2, , drop = FALSE]
    near_nonsingular <- near_tbl[!near_tbl$singular, , drop = FALSE]
    
    if (nrow(near_nonsingular) > 0) {
      best <- all_fits[[near_nonsingular$idx[1]]]
    } else {
      best <- all_fits[[best_row$idx]]
    }
  }
  
  final_model <- best$model
  anova_df <- extract_car_anova_table(final_model)
  anova_string <- build_anova_caption_multiline(anova_df)
  has_sig_non_intercept <- has_non_intercept_sig_terms(anova_df)

  if (layout_label == "all_groups") {
    posthoc_terms <- get_posthoc_terms_layout_i(anova_df)
  } else {
    posthoc_terms <- get_posthoc_terms_layout_ii(anova_df)
  }
  
  posthoc_list <- if (RUN_POSTHOC_ONLY_IF_SIG && !has_sig_non_intercept) {
    list()
  } else {
    run_posthoc_for_terms(final_model, posthoc_terms, type = "response")
  }
  
  cld_all_combinations <- if (RUN_CLD_ONLY_IF_SIG && !has_sig_non_intercept) {
    list()
  } else {
    run_cld_all_combinations(final_model, type = "response")
  }
  
  cld_tables <- if (RUN_CLD_ONLY_IF_SIG && !has_sig_non_intercept) {
    list()
  } else {
    run_cld_for_terms(final_model, posthoc_terms, type = "response")
  }
  
  cld_df <- if ("timepoint:Treatment:Group" %in% names(cld_all_combinations)) {
    cld_all_combinations[["timepoint:Treatment:Group"]]
  } else if ("Treatment:Group" %in% names(cld_all_combinations)) {
    cld_all_combinations[["Treatment:Group"]]
  } else if (length(cld_all_combinations) > 0) {
    cld_all_combinations[[1]]
  } else {
    empty_cld_df()
  }
  
  diag_paths <- save_glmm_diag_plots(
    final_model,
    prefix = paste0(safe_name(indicator), "_layout_", layout_label, "_final")
  )
  
  list(
    indicator = indicator,
    status = "ok",
    note = "Decision path C: model-based alternative preferred for count outcome",
    transformation = list(
      decision_path = "C",
      decision_label = "Model-based alternative preferred",
      chosen_transform = "none"
    ),
    candidate_table = all_rows,
    selected_model = final_model,
    selected_formula = best$formula,
    selected_aicc = best$aicc,
    anova_table = anova_df,
    anova_string = anova_string,
    posthoc = posthoc_list,
    cld = cld_df,
    cld_all_combinations = cld_all_combinations,
    cld_tables = cld_tables,
    diagnostics = diag_paths
  )
}

analyze_angular_indicator <- function(data, indicator, layout_label) {
  if (!ANALYZE_ANGULAR_AS_LINEAR) {
    return(placeholder_result(
      indicator,
      "Angular indicator flagged. No inferential model run by default. Set ANALYZE_ANGULAR_AS_LINEAR <- TRUE to force an approximate Gaussian branch."
    ))
  }
  analyze_continuous_indicator(data, indicator, layout_label, tf = NULL)
}

############################
# 10) Assumption report objects
############################

make_assumption_row_continuous <- function(indicator, tf_obj) {
  data.frame(
    indicator = indicator,
    outcome_class = "continuous",
    decision_path = tf_obj$decision_path,
    decision_label = tf_obj$decision_label,
    chosen_transform = tf_obj$chosen_transform,
    raw_shapiro_W = tf_obj$raw_eval$shapiro$statistic,
    raw_shapiro_p = tf_obj$raw_eval$shapiro$p.value,
    raw_levene_stat = tf_obj$raw_eval$levene$statistic,
    raw_levene_p = tf_obj$raw_eval$levene$p.value,
    transformed_shapiro_W = tf_obj$transformed_eval$shapiro$statistic,
    transformed_shapiro_p = tf_obj$transformed_eval$shapiro$p.value,
    transformed_levene_stat = tf_obj$transformed_eval$levene$statistic,
    transformed_levene_p = tf_obj$transformed_eval$levene$p.value,
    decision_note = tf_obj$decision_note,
    stringsAsFactors = FALSE
  )
}

make_assumption_row_count <- function(indicator) {
  data.frame(
    indicator = indicator,
    outcome_class = "count",
    decision_path = "C",
    decision_label = "Model-based alternative preferred",
    chosen_transform = "none",
    raw_shapiro_W = NA_real_,
    raw_shapiro_p = NA_real_,
    raw_levene_stat = NA_real_,
    raw_levene_p = NA_real_,
    transformed_shapiro_W = NA_real_,
    transformed_shapiro_p = NA_real_,
    transformed_levene_stat = NA_real_,
    transformed_levene_p = NA_real_,
    decision_note = "Count outcome; no Shapiro/Levene transform loop applied",
    stringsAsFactors = FALSE
  )
}

make_assumption_row_angular <- function(indicator) {
  data.frame(
    indicator = indicator,
    outcome_class = "angular",
    decision_path = "C",
    decision_label = "Model-based alternative preferred / flagged",
    chosen_transform = "none",
    raw_shapiro_W = NA_real_,
    raw_shapiro_p = NA_real_,
    raw_levene_stat = NA_real_,
    raw_levene_p = NA_real_,
    transformed_shapiro_W = NA_real_,
    transformed_shapiro_p = NA_real_,
    transformed_levene_stat = NA_real_,
    transformed_levene_p = NA_real_,
    decision_note = "Angular outcome flagged; ordinary Gaussian transform pipeline not applied by default",
    stringsAsFactors = FALSE
  )
}

############################
# 11) Word report writers
############################

write_report_models <- function(model_results, layout_label, subset_label, subset_note = "") {
  doc <- .officer_read_docx()
  doc <- .officer_body_end_section_continuous(doc, a4_section)
  
  doc <- append_text(doc, paste0("GOT Activity Statistics - Model Report"), style = "heading 1")
  doc <- append_text(doc, paste0("Branch: ", layout_label, " | Subset: ", subset_label), style = "heading 2")
  doc <- append_text(doc, subset_note, style = "Normal")
  
  kept_names <- names(model_results)[vapply(
    model_results,
    function(res) {
      is.list(res) &&
        !is.null(res$status) &&
        identical(as.character(res$status)[1], "ok") &&
        !is.null(res$anova_table) &&
        is.data.frame(res$anova_table)
    },
    logical(1)
  )]
  
  if (length(kept_names) == 0) {
    doc <- append_text(doc, "No analysed indicators available for this branch/subset.", style = "Normal")
  } else {
    for (i in seq_along(kept_names)) {
      nm <- kept_names[i]
      res <- model_results[[nm]]
      
      if (i > 1) doc <- append_page_break(doc)
      
      doc <- append_text(doc, pretty_indicator_name(nm), style = "heading 2")
      doc <- append_text(doc, paste0("Model note: ", res$note), style = "Normal")

      # ---- Assumption-diagnostics table (new: from diagnostics module) ------
      if (!is.null(res$diag_row) && is.data.frame(res$diag_row) &&
          nrow(res$diag_row) > 0L) {
        doc <- append_text(doc, "Assumption diagnostics", style = "heading 3")
        doc <- append_ft(doc, res$diag_row)
      }

      sel_tbl <- data.frame(
        selected_formula = ifelse(is.null(res$selected_formula), NA_character_, res$selected_formula),
        selected_aicc = ifelse(is.null(res$selected_aicc), NA_real_, res$selected_aicc),
        stringsAsFactors = FALSE
      )
      doc <- append_text(doc, "Selected model", style = "heading 3")
      doc <- append_ft(doc, sel_tbl)
      
      doc <- append_text(doc, "Candidate model ranking", style = "heading 3")
      doc <- append_ft(doc, res$candidate_table)
      
      doc <- append_text(doc, "Type III ANOVA", style = "heading 3")
      doc <- append_ft(doc, res$anova_table)
      
      doc <- append_text(doc, "ANOVA summary", style = "heading 3")
      doc <- append_text(doc, ifelse(nzchar(res$anova_string), res$anova_string, "ANOVA caption not available."), style = "Normal")
      
      if (length(res$posthoc) > 0) {
        doc <- append_text(doc, "Post hoc tables", style = "heading 3")
        for (term in names(res$posthoc)) {
          ph <- res$posthoc[[term]]
          term_title <- if (identical(term, "timepoint:Treatment:Group")) {
            "Overall interaction post hoc (Interval × Treatment × Group)"
          } else {
            paste0("Post hoc: ", pretty_term_name(term))
          }
          doc <- append_text(doc, term_title, style = "heading 3")
          if (!is.null(ph$table) && is.data.frame(ph$table) && nrow(ph$table) > 0) {
            doc <- append_ft(doc, ph$table)
          } else {
            doc <- append_text(doc, ifelse(is.null(ph$error), "No post hoc table available.", ph$error), style = "Normal")
          }
        }
      }
      
      doc <- append_text(doc, "Compact letter display tables", style = "heading 3")
      
      printed_any_cld <- FALSE
      
      if (!is.null(res$cld_tables) && is.list(res$cld_tables) && length(res$cld_tables) > 0) {
        for (cld_term in names(res$cld_tables)) {
          cld_tbl <- res$cld_tables[[cld_term]]
          
          if (!is.data.frame(cld_tbl) || nrow(cld_tbl) == 0) next
          
          if ("cld_letter" %in% names(cld_tbl)) {
            all_empty_letters <- all(is.na(cld_tbl$cld_letter) | !nzchar(trimws(as.character(cld_tbl$cld_letter))))
            if (all_empty_letters) next
          } else {
            next
          }
          
          cld_title <- if (identical(cld_term, "timepoint:Treatment:Group")) {
            "CLD: overall interaction (Interval × Treatment × Group)"
          } else if (identical(cld_term, "Treatment:Group")) {
            "CLD: overall interaction (Treatment × Group)"
          } else {
            paste0("CLD: ", pretty_term_name(cld_term))
          }
          
          doc <- append_text(doc, cld_title, style = "heading 3")
          doc <- append_ft(doc, cld_tbl)
          # Remapped CLD table for graph-annotation transparency
          if (exists("build_remapped_report_tbl", mode = "function") &&
              is.data.frame(cld_tbl)) {
            .remap_tbl5 <- tryCatch(build_remapped_report_tbl(cld_tbl),
                                     error = function(e) NULL)
            if (!is.null(.remap_tbl5) && nrow(.remap_tbl5) > 0L) {
              doc <- append_text(doc, "Remapped CLD letters used in graphs",
                                 style = "heading 3")
              doc <- append_ft(doc, flextable::flextable(.remap_tbl5))
            }
          }
          printed_any_cld <- TRUE
        }
      }

      if (!printed_any_cld) {
        doc <- append_text(
          doc,
          paste0("No printable CLD table available. Reason: ", summarize_cld_availability(res)),
          style = "Normal"
        )
      }
      
      if (!is.null(res$diagnostics$resid_vs_fitted) && file.exists(res$diagnostics$resid_vs_fitted)) {
        doc <- append_text(doc, "Residuals vs fitted", style = "heading 3")
        doc <- append_plot(doc, res$diagnostics$resid_vs_fitted)
      }
      if (!is.null(res$diagnostics$qq) && file.exists(res$diagnostics$qq)) {
        doc <- append_text(doc, "Normal Q-Q plot", style = "heading 3")
        doc <- append_plot(doc, res$diagnostics$qq)
      }
      if (!is.null(res$diagnostics$dharma) && file.exists(res$diagnostics$dharma)) {
        doc <- append_text(doc, "DHARMa diagnostics", style = "heading 3")
        doc <- append_plot(doc, res$diagnostics$dharma)
      }
    }
  }
  
  path_out <- file.path(
    DIR_REPORTS,
    paste0("Report_models_", safe_name(layout_label), "_", safe_name(subset_label), ".docx")
  )
  .officer_print(doc, path_out)
  path_out
}

write_report_assumptions <- function(qc_tables, assumption_rows, continuous_tf_list, input_info, layout_label, subset_label) {
  doc <- .officer_read_docx()
  doc <- .officer_body_end_section_continuous(doc, a4_section)
  
  doc <- append_text(doc, "GOT Activity Statistics - Assumptions Report", style = "heading 1")
  doc <- append_text(doc, paste0("Branch: ", layout_label, " | Subset: ", subset_label), style = "heading 2")
  doc <- append_text(doc, "Experimental design and input", style = "heading 2")
  doc <- append_text(doc, paste0("Data source: ", input_info$source), style = "Normal")
  doc <- append_text(doc, paste0("Source path: ", ifelse(is.na(input_info$source_path), "NA", input_info$source_path)), style = "Normal")
  doc <- append_text(doc, "Unit of analysis: one fish observation at one timepoint within one trial.", style = "Normal")
  doc <- append_text(doc, "Repeated structure: timepoint repeated within trial_id when timepoint has more than one observed level in the subset.", style = "Normal")
  
  doc <- append_text(doc, "N by trial and timepoint", style = "heading 3")
  doc <- append_ft(doc, qc_tables$qc_counts_by_trial_timepoint)
  
  doc <- append_text(doc, "Exclusions", style = "heading 3")
  doc <- append_ft(doc, qc_tables$qc_exclusion_table)
  
  doc <- append_text(doc, "Missingness", style = "heading 3")
  doc <- append_ft(doc, qc_tables$missingness_table)
  
  doc <- append_text(doc, "Transformation decisions summary", style = "heading 3")
  doc <- append_ft(doc, assumption_rows)
  
  tf_names <- names(continuous_tf_list)
  if (length(tf_names) == 0) {
    doc <- append_text(doc, "No analysed continuous indicators available for this subset.", style = "Normal")
  } else {
    for (i in seq_along(tf_names)) {
      ind <- tf_names[i]
      tf <- continuous_tf_list[[ind]]
      if (i > 1) doc <- append_page_break(doc)
      
      doc <- append_text(doc, pretty_indicator_name(ind), style = "heading 2")
      
      doc <- append_text(doc, "Raw-scale candidate models", style = "heading 3")
      doc <- append_ft(doc, tf$raw_eval$candidate_table)
      
      raw_tbl <- data.frame(
        test = c("Shapiro-Wilk", "Levene"),
        statistic = c(tf$raw_eval$shapiro$statistic, tf$raw_eval$levene$statistic),
        p_value = c(tf$raw_eval$shapiro$p.value, tf$raw_eval$levene$p.value),
        note = c(tf$raw_eval$shapiro$note, tf$raw_eval$levene$note),
        stringsAsFactors = FALSE
      )
      doc <- append_text(doc, "Raw-scale diagnostics", style = "heading 3")
      doc <- append_ft(doc, raw_tbl)
      
      if (!is.null(tf$raw_eval$diag_paths$resid_vs_fitted) && file.exists(tf$raw_eval$diag_paths$resid_vs_fitted)) {
        doc <- append_text(doc, "Raw-scale residuals vs fitted", style = "heading 3")
        doc <- append_plot(doc, tf$raw_eval$diag_paths$resid_vs_fitted)
      }
      if (!is.null(tf$raw_eval$diag_paths$qq) && file.exists(tf$raw_eval$diag_paths$qq)) {
        doc <- append_text(doc, "Raw-scale normal Q-Q plot", style = "heading 3")
        doc <- append_plot(doc, tf$raw_eval$diag_paths$qq)
      }
      
      decision_tbl <- data.frame(
        decision_path = tf$decision_path,
        decision_label = tf$decision_label,
        chosen_transform = tf$chosen_transform,
        decision_note = tf$decision_note,
        stringsAsFactors = FALSE
      )
      doc <- append_text(doc, "Transformation decision", style = "heading 3")
      doc <- append_ft(doc, decision_tbl)
      
      doc <- append_text(doc, "Chosen-scale candidate models", style = "heading 3")
      doc <- append_ft(doc, tf$transformed_eval$candidate_table)
      
      trans_tbl <- data.frame(
        test = c("Shapiro-Wilk", "Levene"),
        statistic = c(tf$transformed_eval$shapiro$statistic, tf$transformed_eval$levene$statistic),
        p_value = c(tf$transformed_eval$shapiro$p.value, tf$transformed_eval$levene$p.value),
        note = c(tf$transformed_eval$shapiro$note, tf$transformed_eval$levene$note),
        stringsAsFactors = FALSE
      )
      doc <- append_text(doc, "Chosen-scale diagnostics", style = "heading 3")
      doc <- append_ft(doc, trans_tbl)
      
      if (!is.null(tf$transformed_eval$diag_paths$resid_vs_fitted) && file.exists(tf$transformed_eval$diag_paths$resid_vs_fitted)) {
        doc <- append_text(doc, "Chosen-scale residuals vs fitted", style = "heading 3")
        doc <- append_plot(doc, tf$transformed_eval$diag_paths$resid_vs_fitted)
      }
      if (!is.null(tf$transformed_eval$diag_paths$qq) && file.exists(tf$transformed_eval$diag_paths$qq)) {
        doc <- append_text(doc, "Chosen-scale normal Q-Q plot", style = "heading 3")
        doc <- append_plot(doc, tf$transformed_eval$diag_paths$qq)
      }
    }
  }
  
  path_out <- file.path(
    DIR_REPORTS,
    paste0("Report_assumptions_", safe_name(layout_label), "_", safe_name(subset_label), ".docx")
  )
  .officer_print(doc, path_out)
  path_out
}

validate_model_result <- function(res, indicator) {
  make_bad <- function(msg) {
    placeholder_result(indicator, paste0("Malformed result replaced: ", msg))
  }
  
  if (is.null(res)) {
    return(make_bad("result was NULL"))
  }
  
  if (!is.list(res)) {
    return(make_bad(paste0("object was not a list (class: ", paste(class(res), collapse = ", "), ")")))
  }
  
  if (is.null(res$indicator)) res$indicator <- indicator
  if (is.null(res$status)) res$status <- "ok"
  if (is.null(res$note)) res$note <- ""
  if (is.null(res$candidate_table) || !is.data.frame(res$candidate_table)) {
    res$candidate_table <- data.frame(note = "Candidate model table unavailable", stringsAsFactors = FALSE)
  }
  if (is.null(res$selected_model)) res$selected_model <- NULL
  if (is.null(res$selected_formula)) res$selected_formula <- NA_character_
  if (is.null(res$selected_aicc)) res$selected_aicc <- NA_real_
  if (is.null(res$anova_table) || !is.data.frame(res$anova_table)) {
    res$anova_table <- data.frame(note = "ANOVA table unavailable", stringsAsFactors = FALSE)
  }
  if (is.null(res$anova_string)) res$anova_string <- ""
  if (is.null(res$posthoc) || !is.list(res$posthoc)) res$posthoc <- list()
  if (is.null(res$cld) || !is.data.frame(res$cld)) res$cld <- empty_cld_df()
  if (is.null(res$cld_all_combinations) || !is.list(res$cld_all_combinations)) res$cld_all_combinations <- list()
  if (is.null(res$cld_tables) || !is.list(res$cld_tables)) res$cld_tables <- list()
  if (is.null(res$diagnostics) || !is.list(res$diagnostics)) res$diagnostics <- list()
  
  res
}

############################
# 11b) Markdown report writers
############################

# Convert a data.frame to a Markdown table string.
# Numeric columns are formatted to 3 decimal places.
# Row names are included as a "Term" column if they are non-trivial.
df_to_md_table <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0L) return("*No data.*")
  df2 <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) formatC(col, digits = 3L, format = "f") else as.character(col)
  }), stringsAsFactors = FALSE, check.names = FALSE)
  rn <- rownames(df)
  if (!is.null(rn) && any(rn != as.character(seq_len(nrow(df))))) {
    df2 <- cbind(Term = rn, df2)
  }
  header <- paste0("| ", paste(names(df2), collapse = " | "), " |")
  sep    <- paste0("| ", paste(rep("---", ncol(df2)), collapse = " | "), " |")
  rows   <- apply(df2, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# Write one Markdown report for a single layout x subset branch.
# model_results: named list (per indicator), each element has fields:
#   $status, $note, $selected_formula, $selected_aicc, $anova_table,
#   $anova_string, $posthoc (list), $cld_tables (list)
write_report_markdown_step5 <- function(model_results, layout_label, subset_label, out_dir) {
  md_dir <- file.path(out_dir, "reports_md")
  dir.create(md_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(md_dir,
    paste0("report_", layout_label, "_", subset_label, ".md"))

  kept <- names(model_results)[vapply(model_results, function(r) {
    is.list(r) && identical(r$status, "ok") &&
      !is.null(r$anova_table) && is.data.frame(r$anova_table)
  }, logical(1L))]

  lines <- c(
    "# GOT Activity Statistics -- Interval-Level Analysis", "",
    paste0("**Branch:** ", layout_label, "  "),
    paste0("**Subset:** ", subset_label, "  "),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "", "---", ""
  )

  # Machine-readable summary index
  lines <- c(lines, "## Summary Index", "", "```yaml")
  for (nm in kept) {
    r  <- model_results[[nm]]
    at <- r$anova_table
    pcol <- attr(at, "p_col")
    if (is.null(pcol)) pcol <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(at))[1]
    primary <- c("timepoint:Treatment:Group", "Treatment:Group",
                 "Group:Treatment", "Treatment:timepoint")
    rn  <- rownames(at)
    hit <- intersect(primary, rn)
    if (length(hit) == 0L) hit <- rn[nrow(at)]
    int_p <- if (!is.null(pcol) && pcol %in% names(at))
      suppressWarnings(as.numeric(at[[pcol]][match(hit[1], rn)])) else NA_real_
    sig_rn <- if (!is.null(pcol) && pcol %in% names(at))
      rn[!is.na(suppressWarnings(as.numeric(at[[pcol]]))) &
         suppressWarnings(as.numeric(at[[pcol]])) < 0.05 &
         rn != "(Intercept)"]
    else character(0L)
    lines <- c(lines,
      paste0("- indicator: ", nm),
      paste0("  label: \"", pretty_indicator_name(nm), "\""),
      paste0("  formula: \"", ifelse(is.null(r$selected_formula), "NA", r$selected_formula), "\""),
      paste0("  aicc: ", ifelse(is.null(r$selected_aicc), "NA", round(r$selected_aicc, 2L))),
      paste0("  interaction_term: ", hit[1]),
      paste0("  interaction_p: ", ifelse(is.na(int_p), "NA", round(int_p, 4L))),
      paste0("  significant_terms: \"", paste(sig_rn, collapse = "; "), "\""),
      ""
    )
  }
  lines <- c(lines, "```", "")

  # Per-indicator sections
  for (nm in kept) {
    r <- model_results[[nm]]
    lines <- c(lines, "---", "", paste0("## ", pretty_indicator_name(nm)), "")
    lines <- c(lines, paste0("**Model note:** ", r$note), "")

    if (!is.null(r$selected_formula)) {
      lines <- c(lines,
        paste0("**Selected formula:** `", r$selected_formula, "`  "),
        paste0("**AICc:** ", round(r$selected_aicc, 3L)),
        "")
    }

    # ANOVA table
    lines <- c(lines, "### Type III ANOVA", "")
    at_pr <- r$anova_table
    if (!("Term" %in% names(at_pr))) at_pr <- cbind(Term = rownames(at_pr), at_pr)
    rownames(at_pr) <- NULL
    lines <- c(lines, df_to_md_table(at_pr), "")

    # ANOVA caption
    anova_cap <- if (!is.null(r$anova_string) && nzchar(r$anova_string))
      r$anova_string else "ANOVA caption not available."
    lines <- c(lines, paste0("**ANOVA caption:** ", anova_cap), "")

    # Post-hoc
    if (!is.null(r$posthoc) && length(r$posthoc) > 0L) {
      lines <- c(lines, "### Post-hoc comparisons", "")
      for (term in names(r$posthoc)) {
        ph <- r$posthoc[[term]]
        lines <- c(lines, paste0("#### ", pretty_term_name(term)), "")
        if (!is.null(ph$table) && is.data.frame(ph$table) && nrow(ph$table) > 0L) {
          lines <- c(lines, df_to_md_table(ph$table), "")
        } else {
          lines <- c(lines, paste0("*", ifelse(is.null(ph$error), "No post-hoc table.", ph$error), "*"), "")
        }
      }
    }

    # CLD tables
    if (!is.null(r$cld_tables) && is.list(r$cld_tables) && length(r$cld_tables) > 0L) {
      lines <- c(lines, "### Compact letter display", "")
      for (term in names(r$cld_tables)) {
        cld_tbl <- r$cld_tables[[term]]
        lines <- c(lines, paste0("#### ", pretty_term_name(term)), "")
        if (!is.null(cld_tbl) && is.data.frame(cld_tbl) && nrow(cld_tbl) > 0L) {
          lines <- c(lines, df_to_md_table(cld_tbl), "")
        } else {
          lines <- c(lines, "*No CLD table.*", "")
        }
      }
    }
  }

  tryCatch(
    writeLines(lines, out_file, useBytes = FALSE),
    error = function(e) warning("Markdown report write failed: ", conditionMessage(e), call. = FALSE)
  )
  message("Markdown report written: ", out_file)
  invisible(out_file)
}

# Write Markdown report for the trial-aggregated branch.
write_report_markdown_aggregated <- function(trial_agg_results, out_dir) {
  md_dir <- file.path(out_dir, "reports_md")
  dir.create(md_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(md_dir, "report_trial_aggregated.md")

  lines <- c(
    "# GOT Activity Statistics -- Trial-Aggregated Analysis", "",
    "**Module:** Trial-aggregated (Treatment x Group; one value per trial, mean across Intervals 2-3)  ",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "", "---", ""
  )

  valid <- names(trial_agg_results)[!vapply(trial_agg_results, is.null, logical(1L))]
  for (ind in valid) {
    r <- trial_agg_results[[ind]]
    if (is.null(r)) next
    lines <- c(lines, paste0("## ", pretty_indicator_name(ind)), "")
    lines <- c(lines,
      paste0("**Formula:** `", ind, " ~ Treatment * Group`  "),
      paste0("**Model type:** ", r$model_type, "  "),
      paste0("**N trials:** ", r$n_trials),
      "")

    at <- r$anova_table
    if (is.data.frame(at) && nrow(at) > 0L) {
      at_pr <- cbind(Term = rownames(at), at)
      rownames(at_pr) <- NULL
      lines <- c(lines, "### Type III ANOVA", "", df_to_md_table(at_pr), "")
    }
    lines <- c(lines, "---", "")
  }

  tryCatch(
    writeLines(lines, out_file, useBytes = FALSE),
    error = function(e) warning("Aggregated MD report write failed: ", conditionMessage(e), call. = FALSE)
  )
  message("Aggregated Markdown report written: ", out_file)
  invisible(out_file)
}

############################
# 12) Main analysis driver
############################

analyze_layout <- function(data, layout_label = c("all_groups", "exe_choice"), precomputed_tf = NULL) {
  layout_label <- match.arg(layout_label)

  dat <- data

  if (layout_label == "exe_choice") {
    dat <- dat[dat$Group %in% c("Exercise choice", "No choice"), , drop = FALSE]
    dat <- droplevels(dat)
  }

  # OPT-2: Parallelize per-indicator LMM fits via future_lapply
  # Each indicator is processed independently with its own response column and AICc candidate set.

  # Build worker globals: functions and atomic constants from globalenv() + local data.
  # Exclude large list/data-frame objects (those are not atomic) to keep serialization lean.
  # Raise future.globals.maxSize to 4 GB so that the function closures from ggplot2/officer
  # that land in globalenv after STEP4 don't trigger a FutureLaunchError size cap.
  .all_globals <- eapply(globalenv(), identity)
  .wk_globals  <- Filter(function(x) is.function(x) || (is.atomic(x) && !is.null(x)),
                         .all_globals)
  # Locally scoped objects not in globalenv
  .wk_globals[["dat"]]            <- dat
  .wk_globals[["precomputed_tf"]] <- precomputed_tf
  .wk_globals[["layout_label"]]   <- layout_label

  # Raise globals size cap to 4 GB (default 500 MB is too small after STEP4)
  .old_maxsize <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 4 * 1024^3)
  on.exit(options(future.globals.maxSize = .old_maxsize), add = TRUE)

  .run_parallel_indicators <- function() {
    future.apply::future_lapply(
    ALL_INDICATORS,
    function(indicator) {
      # FIX-2: Per-indicator seed using hash of indicator name for DHARMa reproducibility
      # ensures reproducibility regardless of execution order
      seed_val <- 42L + sum(as.integer(charToRaw(indicator))) %% 10000L
      set.seed(seed_val)

      # BLAS thread safety for lme4/lmerTest
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1L)
      }

      message("Running ", layout_label, " :: ", indicator)

      if (!indicator %in% names(dat)) {
        return(placeholder_result(indicator, "Indicator column missing in data"))
      }

      local_dat <- dat[!is.na(dat[[indicator]]), , drop = FALSE]

      if (nrow(local_dat) < 5) {
        return(placeholder_result(indicator, "Too few non-missing rows"))
      }

      res <- tryCatch({
        if (indicator %in% CONTINUOUS_INDICATORS) {
          tf_use <- NULL

          if (!is.null(precomputed_tf) &&
              layout_label == "all_groups" &&
              indicator %in% names(precomputed_tf)) {
            tf_use <- precomputed_tf[[indicator]]
          }

          analyze_continuous_indicator(local_dat, indicator, layout_label, tf = tf_use)

        } else if (indicator %in% BETA_INDICATORS) {
          analyze_beta_indicator(local_dat, indicator, layout_label)

        } else if (indicator %in% COUNT_INDICATORS) {
          analyze_count_indicator(local_dat, indicator, layout_label)

        } else if (indicator %in% ANGULAR_INDICATORS) {
          analyze_angular_indicator(local_dat, indicator, layout_label)

        } else {
          placeholder_result(indicator, "Indicator class not mapped")
        }
      }, error = function(e) {
        placeholder_result(
          indicator,
          paste0("Analysis failed and was replaced with placeholder_result(): ", conditionMessage(e))
        )
      })

      # Strip model objects before crossing process boundary:
      # glmmTMB models contain C++ pointers that cannot be serialized over multisession sockets.
      # lme4 models are safe but large; stripping avoids RAM overhead across 6 workers.
      # All model-derived statistics (ANOVA, post-hoc, CLD, diagnostics) are already in res.
      if (is.list(res)) res$selected_model <- NULL

      validate_model_result(res, indicator)
    },
    future.packages = c(
      "lme4", "lmerTest", "emmeans", "DHARMa", "car",
      "multcomp", "multcompView", "glmmTMB", "AICcmodavg",
      "ggplot2", "officer", "flextable", "stringr", "readxl", "dplyr", "tidyr"
    ),
    future.globals = .wk_globals,
    future.seed = TRUE
  )
  }

  # Launch parallel; fall back to sequential if workers fail to start or crash.
  all_results <- tryCatch(
    .run_parallel_indicators(),
    FutureLaunchError = function(e) {
      warning("OPT-2: FutureLaunchError caught — falling back to sequential lapply. ",
              "Error was: ", conditionMessage(e))
      future::plan(future::sequential())
      lapply(ALL_INDICATORS, function(indicator) {
        message("Running (sequential fallback) :: ", indicator)
        if (!indicator %in% names(dat)) return(placeholder_result(indicator, "Indicator column missing"))
        local_dat <- dat[!is.na(dat[[indicator]]), , drop = FALSE]
        if (nrow(local_dat) < 5) return(placeholder_result(indicator, "Too few rows"))
        tryCatch({
          if (indicator %in% CONTINUOUS_INDICATORS) {
            tf_use <- if (!is.null(precomputed_tf) && layout_label == "all_groups" &&
                          indicator %in% names(precomputed_tf)) precomputed_tf[[indicator]] else NULL
            analyze_continuous_indicator(local_dat, indicator, layout_label, tf = tf_use)
          } else if (indicator %in% BETA_INDICATORS) {
            analyze_beta_indicator(local_dat, indicator, layout_label)
          } else if (indicator %in% COUNT_INDICATORS) {
            analyze_count_indicator(local_dat, indicator, layout_label)
          } else if (indicator %in% ANGULAR_INDICATORS) {
            analyze_angular_indicator(local_dat, indicator, layout_label)
          } else {
            placeholder_result(indicator, "Indicator class not mapped")
          }
        }, error = function(e2) placeholder_result(indicator, conditionMessage(e2)))
      })
    },
    error = function(e) {
      warning("OPT-2: unexpected error in parallel block — falling back to sequential. ",
              conditionMessage(e))
      future::plan(future::sequential())
      lapply(ALL_INDICATORS, function(indicator) {
        message("Running (sequential fallback) :: ", indicator)
        if (!indicator %in% names(dat)) return(placeholder_result(indicator, "Indicator column missing"))
        local_dat <- dat[!is.na(dat[[indicator]]), , drop = FALSE]
        if (nrow(local_dat) < 5) return(placeholder_result(indicator, "Too few rows"))
        tryCatch({
          if (indicator %in% CONTINUOUS_INDICATORS) {
            tf_use <- if (!is.null(precomputed_tf) && layout_label == "all_groups" &&
                          indicator %in% names(precomputed_tf)) precomputed_tf[[indicator]] else NULL
            analyze_continuous_indicator(local_dat, indicator, layout_label, tf = tf_use)
          } else if (indicator %in% BETA_INDICATORS) {
            analyze_beta_indicator(local_dat, indicator, layout_label)
          } else if (indicator %in% COUNT_INDICATORS) {
            analyze_count_indicator(local_dat, indicator, layout_label)
          } else {
            analyze_angular_indicator(local_dat, indicator, layout_label)
          }
        }, error = function(e2) placeholder_result(indicator, conditionMessage(e2)))
      })
    }
  )

  # Ensure names match ALL_INDICATORS
  names(all_results) <- ALL_INDICATORS

  stopifnot(
    length(all_results) == length(ALL_INDICATORS),
    all(names(all_results) == ALL_INDICATORS),
    !anyNA(names(all_results))
  )
  cat(sprintf("✓ Ran %d indicators (%s)\n",
              length(all_results),
              if (future::nbrOfWorkers() > 1L) paste0(future::nbrOfWorkers(), " workers") else "sequential"))

  # FIX-6: Check for placeholder results (early failure detection)
  n_failed <- sum(vapply(all_results, function(r) isTRUE(attr(r, "is_placeholder")), logical(1)))
  if (n_failed > 0) {
    warning(sprintf("OPT-2: %d/%d indicators failed; returned placeholders",
                    n_failed, length(all_results)))
  }

  # FIX-5: Return results directly; inline validation in lambda is sufficient
  # Removed redundant post-hoc lapply(validate_model_result) loop
  all_results
}
  
############################
# 13) Build summary outputs
############################

save_model_detail_tables <- function(model_results, suffix) {
  for (nm in names(model_results)) {
    res <- model_results[[nm]]
    
    if (!is.list(res) || !identical(res$status, "ok")) next
    
    if (!is.null(res$candidate_table) && is.data.frame(res$candidate_table)) {
      save_csv(
        res$candidate_table,
        file.path(DIR_TABLES, paste0("candidate_models_", suffix, "_", safe_name(nm), ".csv"))
      )
    }
    
    if (!is.null(res$anova_table) && is.data.frame(res$anova_table)) {
      save_csv(
        trim_table_numeric(res$anova_table),
        file.path(DIR_TABLES, paste0("anova_", suffix, "_", safe_name(nm), ".csv"))
      )
    }
    
    if (!is.null(res$posthoc) && length(res$posthoc) > 0) {
      for (term in names(res$posthoc)) {
        ph <- res$posthoc[[term]]
        if (!is.null(ph$table) && is.data.frame(ph$table)) {
          save_csv(
            trim_table_numeric(ph$table),
            file.path(
              DIR_TABLES,
              paste0("posthoc_", suffix, "_", safe_name(nm), "_", safe_name(term), ".csv")
            )
          )
        }
      }
    }
    
    if (!is.null(res$cld_all_combinations) && is.list(res$cld_all_combinations) && length(res$cld_all_combinations) > 0) {
      for (term in names(res$cld_all_combinations)) {
        cld_tbl <- res$cld_all_combinations[[term]]
        if (!is.null(cld_tbl) && is.data.frame(cld_tbl) && nrow(cld_tbl) > 0) {
          save_csv(
            trim_table_numeric(cld_tbl),
            file.path(
              DIR_TABLES,
              paste0("cld_all_combinations_", suffix, "_", safe_name(nm), "_", safe_name(term), ".csv")
            )
          )
        }
      }
    }
    
    if (!is.null(res$cld_tables) && is.list(res$cld_tables) && length(res$cld_tables) > 0) {
      for (term in names(res$cld_tables)) {
        cld_tbl <- res$cld_tables[[term]]
        if (!is.null(cld_tbl) && is.data.frame(cld_tbl) && nrow(cld_tbl) > 0) {
          save_csv(
            trim_table_numeric(cld_tbl),
            file.path(
              DIR_TABLES,
              paste0("cld_by_term_", suffix, "_", safe_name(nm), "_", safe_name(term), ".csv")
            )
          )
        }
      }
    }
    
    if (!is.null(res$cld) && is.data.frame(res$cld) && nrow(res$cld) > 0) {
      save_csv(
        trim_table_numeric(res$cld),
        file.path(DIR_TABLES, paste0("cld_", suffix, "_", safe_name(nm), ".csv"))
      )
    }
  }
}

make_anova_string_df <- function(model_results) {
  data.frame(
    indicator = names(model_results),
    anova_string = vapply(
      model_results,
      function(x) {
        if (is.list(x) && !is.null(x$anova_string)) x$anova_string else NA_character_
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

make_cld_list <- function(model_results) {
  out <- setNames(vector("list", length(model_results)), names(model_results))
  
  for (nm in names(model_results)) {
    res <- model_results[[nm]]
    
    if (!is.list(res)) {
      out[[nm]] <- empty_cld_df()
      next
    }
    
    cld_list_obj <- res[["cld_all_combinations"]]
    
    if (is.null(cld_list_obj) || !is.list(cld_list_obj) || length(cld_list_obj) == 0) {
      out[[nm]] <- empty_cld_df()
      next
    }
    
    cld_parts <- lapply(names(cld_list_obj), function(term_nm) {
      df <- cld_list_obj[[term_nm]]
      
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(NULL)
      }
      
      if (!"timepoint" %in% names(df)) df$timepoint <- NA_character_
      if (!"Group" %in% names(df)) df$Group <- NA_character_
      if (!"Treatment" %in% names(df)) df$Treatment <- NA_character_
      if (!"cld_letter" %in% names(df)) df$cld_letter <- NA_character_
      if (!"term" %in% names(df)) df$term <- term_nm
      
      df[, c("timepoint", "Group", "Treatment", "cld_letter", "term"), drop = FALSE]
    })
    
    cld_parts <- cld_parts[!vapply(cld_parts, is.null, logical(1))]
    
    if (length(cld_parts) == 0) {
      out[[nm]] <- empty_cld_df()
    } else {
      out[[nm]] <- do.call(rbind, cld_parts)
      rownames(out[[nm]]) <- NULL
    }
  }
  
  out
}

############################
# 14) Execute pipeline
############################

subset_time_branch <- function(data, subset_label = c("all_timepoints", "timepoints23", "timepoint1")) {
  subset_label <- match.arg(subset_label)
  
  dat <- data
  
  if (subset_label == "timepoints23") {
    dat <- dat[as.character(dat$timepoint) %in% c("2", "3"), , drop = FALSE]
  } else if (subset_label == "timepoint1") {
    dat <- dat[as.character(dat$timepoint) %in% c("1"), , drop = FALSE]
  }
  
  droplevels(dat)
}

run_one_subbranch <- function(data, layout_label, subset_label, branch_dir, input_info) {
  subbranch_dir <- file.path(branch_dir, subset_label)
  dir.create(subbranch_dir, recursive = TRUE, showWarnings = FALSE)
  set_output_dirs(subbranch_dir)
  
  sink(file.path(DIR_LOGS, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
  
  dat_sub <- subset_time_branch(data, subset_label)
  
  if (nrow(dat_sub) == 0) {
    note_tbl <- data.frame(note = paste("No rows available for", layout_label, "/", subset_label), stringsAsFactors = FALSE)
    save_csv(note_tbl, file.path(DIR_TABLES, "EMPTY_SUBSET.csv"))
    return(list(
      assumption_results = note_tbl,
      model_results = list(),
      anova_string = note_tbl,
      cld_list = list(),
      output_manifest = note_tbl
    ))
  }
  
  continuous_tf_list <- list()
  assumption_rows <- list()
  
  for (indicator in ALL_INDICATORS) {
    if (!indicator %in% names(dat_sub)) next
    local_dat <- dat_sub[!is.na(dat_sub[[indicator]]), , drop = FALSE]
    if (nrow(local_dat) < 5) next
    
    if (indicator %in% CONTINUOUS_INDICATORS) {
      tf <- choose_transformation(local_dat, indicator)
      continuous_tf_list[[indicator]] <- tf
      assumption_rows[[indicator]] <- make_assumption_row_continuous(indicator, tf)
    } else if (indicator %in% BETA_INDICATORS) {
      assumption_rows[[indicator]] <- data.frame(
        indicator = indicator,
        outcome_class = "beta",
        decision_path = "B",
        decision_label = "Beta GLMM",
        chosen_transform = "prop = raw_s / 600",
        raw_shapiro_W = NA_real_, raw_shapiro_p = NA_real_,
        raw_levene_stat = NA_real_, raw_levene_p = NA_real_,
        transformed_shapiro_W = NA_real_, transformed_shapiro_p = NA_real_,
        transformed_levene_stat = NA_real_, transformed_levene_p = NA_real_,
        decision_note = "Bounded proportion; no Shapiro/Levene transform loop applied",
        stringsAsFactors = FALSE
      )
    } else if (indicator %in% COUNT_INDICATORS) {
      assumption_rows[[indicator]] <- make_assumption_row_count(indicator)
    } else if (indicator %in% ANGULAR_INDICATORS) {
      assumption_rows[[indicator]] <- make_assumption_row_angular(indicator)
    }
  }
  
  assumption_results <- if (length(assumption_rows) > 0) do.call(rbind, assumption_rows) else {
    data.frame(note = "No assumption rows generated for this subset.", stringsAsFactors = FALSE)
  }
  save_csv(assumption_results, file.path(DIR_TABLES, "assumption_results.csv"))
  
  precomputed_tf_use <- if (layout_label == "all_groups") continuous_tf_list else NULL
  
  model_results <- analyze_layout(
    dat_sub,
    layout_label = layout_label,
    precomputed_tf = precomputed_tf_use
  )
  
  status_table <- data.frame(
    indicator = names(model_results),
    status = vapply(
      model_results,
      function(x) if (is.list(x) && "status" %in% names(x)) as.character(x$status)[1] else NA_character_,
      character(1)
    ),
    note = vapply(
      model_results,
      function(x) if (is.list(x) && "note" %in% names(x)) as.character(x$note)[1] else NA_character_,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  
  save_csv(status_table, file.path(DIR_TABLES, "analysis_status_table.csv"))
  
  anova_string <- make_anova_string_df(model_results)
  save_csv(anova_string, file.path(DIR_TABLES, paste0("anova_string_", layout_label, ".csv")))
  
  cld_list <- make_cld_list(model_results)
  for (nm in names(cld_list)) {
    save_csv(cld_list[[nm]], file.path(DIR_TABLES, paste0("cld_", layout_label, "_", safe_name(nm), ".csv")))
  }
  
  saveRDS(model_results, file.path(DIR_RDS, paste0("model_results_", layout_label, ".rds")))
  saveRDS(cld_list, file.path(DIR_RDS, paste0("cld_list_", layout_label, ".rds")))
  saveRDS(continuous_tf_list, file.path(DIR_RDS, "continuous_transformation_objects.rds"))
  
  save_model_detail_tables(model_results, layout_label)
  
  report1_path <- write_report_assumptions(
    qc_tables = list(
      qc_exclusion_table = qc$qc_exclusion_table,
      qc_counts_by_trial_timepoint = qc$qc_counts_by_trial_timepoint,
      missingness_table = qc$missingness_table
    ),
    assumption_rows = assumption_results,
    continuous_tf_list = continuous_tf_list,
    input_info = input_info,
    layout_label = layout_label,
    subset_label = subset_label
  )
  
  subset_note <- paste0(
    "Layout ", layout_label,
    "; subset = ", subset_label, "."
  )
  
  report_model_path <- write_report_models(
    model_results = model_results,
    layout_label = layout_label,
    subset_label = subset_label,
    subset_note = subset_note
  )
  
  # Results-drafting module (sourced externally via results_drafting_module.R)
  results_draft_path <- NA_character_
  if (exists("write_results_report", mode = "function")) {
    results_draft_path <- tryCatch(
      write_results_report(model_results, layout_label, subset_label),
      error = function(e) {
        warning("Results drafting failed: ", e$message)
        NA_character_
      }
    )
  }

  output_manifest <- data.frame(
    file = c(
      report1_path,
      report_model_path,
      if (!is.na(results_draft_path)) results_draft_path else NULL,
      file.path(DIR_TABLES, "assumption_results.csv"),
      file.path(DIR_RDS, paste0("model_results_", layout_label, ".rds")),
      file.path(DIR_RDS, paste0("cld_list_", layout_label, ".rds"))
    ),
    stringsAsFactors = FALSE
  )
  save_csv(output_manifest, file.path(DIR_TABLES, "output_manifest.csv"))
  
  list(
    assumption_results = assumption_results,
    model_results = model_results,
    anova_string = anova_string,
    cld_list = cld_list,
    output_manifest = output_manifest
  )
}

input_info <- locate_input_data()
fish_activity_summary <- input_info$data

required_design_cols <- c("trial_id", "timepoint", "Group", "Treatment")
check_required_columns(fish_activity_summary, c(required_design_cols, ALL_INDICATORS))

fish_activity_summary <- coerce_design_columns(fish_activity_summary)

qc <- run_qc(fish_activity_summary)
fish_activity_summary_clean <- qc$data

save_csv(fish_activity_summary_clean, file.path(OUTDIR, "fish_activity_summary_clean.csv"))
save_csv(qc$qc_exclusion_table, file.path(OUTDIR, "qc_exclusion_table.csv"))
save_csv(qc$qc_counts_by_trial_timepoint, file.path(OUTDIR, "qc_counts_by_trial_timepoint.csv"))
save_csv(qc$missingness_table, file.path(OUTDIR, "missingness_table.csv"))

all_branch_results <- list()

# OPT-2: Set parallel plan at top level for proper scoping in analyze_layout()
n_workers_step5 <- if (exists("STEP5_PARALLEL_WORKERS"))
  as.integer(STEP5_PARALLEL_WORKERS) else 0L
if (n_workers_step5 > 0L) {
  message("[OPT-2] Setting up parallel plan with ", n_workers_step5, " workers for STEP 5")
  future::plan(future::multisession(workers = n_workers_step5))
}

# Gate fish_id RE candidates based on BL-match reliability tier.
# Looks for a BL match result in the global environment (set by match BL script).
.bl_match_df <- if (exists("bl_match_result", envir = .GlobalEnv, inherits = FALSE))
  get("bl_match_result", envir = .GlobalEnv) else NULL
RANDOM_CANDIDATES <- gate_random_candidates(.bl_match_df)
message("Active RE candidates: ", paste(RANDOM_CANDIDATES, collapse = " | "))

active_layouts  <- if (isTRUE(EXCLUDE_PHARMA_GROUPS)) "exe_choice" else c("all_groups", "exe_choice")
active_subsets  <- if (isTRUE(USE_TP1)) c("all_timepoints", "timepoints23", "timepoint1") else c("timepoints23")

.got_write_run_settings(
  OUTDIR,
  EXCLUDE_PHARMA_GROUPS = EXCLUDE_PHARMA_GROUPS,
  USE_TP1               = USE_TP1,
  active_layouts        = paste(active_layouts,  collapse = ","),
  active_subsets        = paste(active_subsets,  collapse = ",")
)

for (layout_label in active_layouts) {
  branch_dir <- file.path(OUTDIR, layout_label)
  dir.create(branch_dir, recursive = TRUE, showWarnings = FALSE)

  all_branch_results[[layout_label]] <- list()

  for (subset_label in active_subsets) {
    message("Running ", layout_label, " / ", subset_label)
    all_branch_results[[layout_label]][[subset_label]] <- run_one_subbranch(
      data = fish_activity_summary_clean,
      layout_label = layout_label,
      subset_label = subset_label,
      branch_dir = branch_dir,
      input_info = input_info
    )
  }
}

saveRDS(all_branch_results, file.path(OUTDIR, "all_branch_results.rds"))

# ---- R6: also write to canonical stats_cache/ in dataset/ -------------------
# Stable, fixed-name location so STEP 4's get_anova_string_from_results() can
# find the per-indicator RDS files without having to guess the latest
# STEP5_output_<TS>/ folder.
if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
  .stats_cache_dir <- file.path(PIPELINE_DATA_DIR, "dataset", "stats_cache")
  dir.create(.stats_cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(all_branch_results,
          file.path(.stats_cache_dir, "all_branch_results.rds"))
  # Mirror the per-indicator r_objects/ tree (one file per layout) so callers
  # can do file.path(.stats_cache_dir, "<layout>", "model_results_<layout>.rds")
  for (.layout in names(all_branch_results)) {
    .layout_dir <- file.path(.stats_cache_dir, .layout)
    dir.create(.layout_dir, recursive = TRUE, showWarnings = FALSE)
    .res <- all_branch_results[[.layout]]
    if (length(.res) > 0L) {
      for (.subset in names(.res)) {
        .sub_dir <- file.path(.layout_dir, .subset)
        dir.create(.sub_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(.res[[.subset]],
                file.path(.sub_dir, "branch_results.rds"))
      }
    }
  }
  message("R6: stats cache mirrored to canonical: ", .stats_cache_dir)
  rm(.stats_cache_dir, .layout, .layout_dir, .res, .subset, .sub_dir)
}

# =============================================================================
# BH-adjusted p-values across the outcome family (informational)
# =============================================================================
# For each layout × subset, collect the p-value of the primary interaction term
# (Treatment:Group or timepoint:Treatment:Group) and apply Benjamini-Hochberg.
for (layout_nm in names(all_branch_results)) {
  for (subset_nm in names(all_branch_results[[layout_nm]])) {
    res_list <- all_branch_results[[layout_nm]][[subset_nm]]$model_results
    if (!is.list(res_list) || length(res_list) == 0) next

    p_rows <- lapply(names(res_list), function(ind) {
      r <- res_list[[ind]]
      if (!is.list(r) || !identical(r$status, "ok")) return(NULL)
      at <- r$anova_table
      if (!is.data.frame(at) || nrow(at) == 0) return(NULL)
      p_col <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(at))
      if (length(p_col) == 0) return(NULL)
      primary_terms <- c("timepoint:Treatment:Group", "Treatment:Group",
                         "Group:Treatment", "Treatment:timepoint")
      rn <- rownames(at)
      hit <- intersect(primary_terms, rn)
      if (length(hit) == 0) hit <- rn[nrow(at)]
      data.frame(indicator = ind, term = hit[1],
                 p_value = as.numeric(at[hit[1], p_col[1]]),
                 stringsAsFactors = FALSE)
    })
    p_df <- do.call(rbind, Filter(Negate(is.null), p_rows))
    if (!is.null(p_df) && nrow(p_df) > 0) {
      p_df$p_BH <- p.adjust(p_df$p_value, method = "BH")
      bh_path <- file.path(OUTDIR, layout_nm, subset_nm, "tables",
                           paste0("BH_adjusted_pvalues_", layout_nm, "_", subset_nm, ".csv"))
      dir.create(dirname(bh_path), recursive = TRUE, showWarnings = FALSE)
      save_csv(p_df, bh_path)
    }
  }
}

# -- Markdown reports (interval-level) ----------------------------------------
for (layout_nm in names(all_branch_results)) {
  for (subset_nm in names(all_branch_results[[layout_nm]])) {
    mr <- all_branch_results[[layout_nm]][[subset_nm]]$model_results
    if (is.list(mr) && length(mr) > 0L) {
      tryCatch(
        write_report_markdown_step5(mr, layout_nm, subset_nm, OUTDIR),
        error = function(e) warning("MD report failed [", layout_nm, "/", subset_nm, "]: ",
                                    conditionMessage(e), call. = FALSE)
      )
    }
  }
}

# =============================================================================
# Trial-aggregated 2-way branch (Treatment * Group, pseudoreplication-free)
# =============================================================================
# Aggregate fish × TP data to one row per trial (mean across TP2/TP3 and fish).
# Fits lm(y ~ Treatment * Group) with optional (1|trial_date) if not collinear.
# Robust to fish-identity non-preservation across timepoints.
tryCatch({
  if (exists("aggregate_trial_level", mode = "function")) {
    trial_agg_dir <- file.path(OUTDIR, "trial_aggregated_2way")
    dir.create(trial_agg_dir, recursive = TRUE, showWarnings = FALSE)

    dat_agg_full <- fish_activity_summary_clean
    if (isTRUE(EXCLUDE_PHARMA_GROUPS)) {
      dat_agg_full <- dat_agg_full[
        !dat_agg_full$Group %in% c("Morphine", "Naloxone"), , drop = FALSE]
    }

    outcome_cols_agg <- intersect(
      c(CONTINUOUS_INDICATORS, BETA_INDICATORS, COUNT_INDICATORS),
      names(dat_agg_full)
    )

    # Convert time_spent_zone_flow_s to proportion for Beta branch
    if ("time_spent_zone_flow_s" %in% outcome_cols_agg && "time_spent_zone_flow_s" %in% names(dat_agg_full)) {
      dat_agg_full$prop_in_flow <- pmax(
        pmin(dat_agg_full$time_spent_zone_flow_s / RECORDING_DURATION_SECONDS, 1 - 1e-4), 1e-4)
      outcome_cols_agg <- c(setdiff(outcome_cols_agg, "time_spent_zone_flow_s"), "prop_in_flow")
    }

    dat_trial <- aggregate_trial_level(dat_agg_full, outcome_cols_agg, tps_to_avg = c(2L, 3L))

    trial_agg_results <- lapply(outcome_cols_agg, function(ind) {
      if (!ind %in% names(dat_trial)) return(NULL)
      d <- dat_trial[is.finite(dat_trial[[ind]]), , drop = FALSE]
      if (nrow(d) < 5) return(NULL)

      # Check trial_date collinearity with trial_id
      use_date_re <- "trial_date" %in% names(d) && length(unique(d$trial_date)) > 1 &&
        tryCatch({
          cor_td <- abs(cor(as.numeric(factor(d$trial_date)),
                            as.numeric(factor(d$trial_id)), use = "complete.obs"))
          is.finite(cor_td) && cor_td < 0.95
        }, error = function(e) FALSE)

      fit <- if (use_date_re) {
        tryCatch(
          lme4::lmer(stats::as.formula(paste(ind, "~ Treatment * Group + (1|trial_date)")), data = d, REML = FALSE),
          error = function(e) NULL
        )
      } else NULL

      if (is.null(fit)) {
        fit <- tryCatch(stats::lm(stats::as.formula(paste(ind, "~ Treatment * Group")), data = d),
                        error = function(e) NULL)
      }
      if (is.null(fit)) return(NULL)

      anova_tab <- tryCatch(
        as.data.frame(car::Anova(fit, type = 3)),
        error = function(e) data.frame()
      )

      # Post-hoc assumption diagnostics for this aggregated model
      .s5_agg_diag <- if (exists("run_agg_diag", mode = "function")) {
        tryCatch({
          .s5_grp <- tryCatch(
            interaction(d$Group, d$Treatment, drop = TRUE),
            error = function(e) rep("all", nrow(d))
          )
          run_agg_diag(fit, .s5_grp, ind, n_obs = nrow(d))
        }, error = function(e) NULL)
      } else NULL

      list(indicator = ind, model = fit, anova_table = anova_tab,
           n_trials = nrow(d), model_type = class(fit)[1L],
           diagnostics = .s5_agg_diag)
    })
    names(trial_agg_results) <- outcome_cols_agg

    # Save ANOVA summary CSV
    anova_summary_agg <- do.call(rbind, lapply(names(trial_agg_results), function(ind) {
      r <- trial_agg_results[[ind]]
      if (is.null(r) || !is.data.frame(r$anova_table) || nrow(r$anova_table) == 0) return(NULL)
      at <- r$anova_table
      p_col <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(at))
      cbind(data.frame(indicator = ind, term = rownames(at),
                       n_trials = r$n_trials, model_type = r$model_type,
                       stringsAsFactors = FALSE),
            at[, union(p_col, names(at)), drop = FALSE])
    }))
    if (!is.null(anova_summary_agg) && nrow(anova_summary_agg) > 0) {
      p_col_sum <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(anova_summary_agg))
      if (length(p_col_sum) > 0) {
        anova_summary_agg$p_BH <- p.adjust(anova_summary_agg[[p_col_sum[1]]], method = "BH")
      }
      save_csv(anova_summary_agg, file.path(trial_agg_dir, "trial_aggregated_anova_summary.csv"))
    }
    save_csv(dat_trial, file.path(trial_agg_dir, "trial_aggregated_data.csv"))
    saveRDS(trial_agg_results, file.path(trial_agg_dir, "trial_aggregated_results.rds"))
    message("Trial-aggregated 2-way branch written to: ", trial_agg_dir)

    # Export to GlobalEnv so STEP6 (when run in the same session) can access
    # the results without disk I/O, and also persist to stats_cache for
    # standalone STEP6 runs.
    assign("trial_agg_results_step5", trial_agg_results, envir = .GlobalEnv)
    if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
      .sc5 <- file.path(get("PIPELINE_DATA_DIR", envir = .GlobalEnv),
                        "dataset", "stats_cache")
      dir.create(.sc5, recursive = TRUE, showWarnings = FALSE)
      tryCatch(
        saveRDS(trial_agg_results, file.path(.sc5, "trial_agg_results_step5.rds")),
        error = function(e)
          warning("stats_cache trial_agg save failed: ", e$message, call. = FALSE)
      )
      rm(.sc5)
    }

    # ── Aggregated graphs ─────────────────────────────────────────────────
    png_dir_agg <- file.path(trial_agg_dir, "interval_aggregated_png")
    dir.create(png_dir_agg, recursive = TRUE, showWarnings = FALSE)

    grp_cols_agg <- c("No choice" = "#56B4E9", "Exercise choice" = "#D55E00")
    trt_shapes_agg <- c("Control" = 17L, "Exercise" = 16L)
    grp_order_agg  <- c("No choice", "Exercise choice")
    trt_order_agg  <- c("Control", "Exercise")

    for (ind in names(trial_agg_results)) {
      r <- trial_agg_results[[ind]]
      if (is.null(r) || is.null(r$model)) next
      d_plot <- dat_trial[is.finite(dat_trial[[ind]]), , drop = FALSE]
      if (nrow(d_plot) < 3) next

      # Factor ordering
      d_plot$Group     <- factor(d_plot$Group,     levels = grp_order_agg)
      d_plot$Treatment <- factor(d_plot$Treatment, levels = trt_order_agg)

      # Build ANOVA caption
      anova_df_i    <- tryCatch(extract_car_anova_table(r$model), error = function(e) NULL)
      anova_cap_agg <- if (!is.null(anova_df_i) && nrow(anova_df_i) > 0) {
        tryCatch(build_anova_caption_multiline(anova_df_i, caption_format = "plain"),
                 error = function(e) "")
      } else ""

      n_trials_i <- r$n_trials
      caption_txt <- paste0(
        "n = ", n_trials_i, " trials (mean across Intervals 2–3)",
        if (nzchar(anova_cap_agg)) paste0("\n", anova_cap_agg) else ""
      )

      pos_jd_agg <- ggplot2::position_jitterdodge(jitter.width = 0.15, jitter.height = 0,
                                                   dodge.width = 0.4, seed = 42L)

      p_agg <- ggplot2::ggplot(d_plot,
                               ggplot2::aes(x = Group, y = .data[[ind]],
                                            colour = Group, shape = Treatment)) +
        ggplot2::geom_point(position = pos_jd_agg, size = 2, alpha = 0.8) +
        ggplot2::stat_summary(
          ggplot2::aes(group = Treatment),
          fun.data = mean_se, geom = "pointrange",
          position = ggplot2::position_dodge(width = 0.4),
          size = 0.6, colour = "black", show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(values = grp_cols_agg, breaks = grp_order_agg) +
        ggplot2::scale_shape_manual(values = trt_shapes_agg, breaks = trt_order_agg) +
        ggplot2::labs(
          y       = pretty_indicator_name(ind),
          x       = "Group",
          colour  = "Group",
          shape   = "Treatment",
          caption = caption_txt
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position  = "right",
          plot.caption     = ggplot2::element_text(size = 8, hjust = 0),
          axis.title       = ggplot2::element_text(size = 11),
          axis.text        = ggplot2::element_text(size = 10)
        )

      # --- CLD annotation (aggregated graphs) ----------------------------------
      # x = Group, dodge = Treatment (matches stat_summary position_dodge 0.4)
      if (exists("compute_agg_cld_df",    mode = "function") &&
          exists("add_cld_to_scatter",    mode = "function")) {
        .cld_agg <- tryCatch(
          compute_agg_cld_df(r$model, term = ~ Treatment * Group),
          error = function(e) NULL
        )
        if (!is.null(.cld_agg)) {
          .cld_agg <- tryCatch(
            remap_cld_letters(.cld_agg,
                              treat_col = "Treatment", grp_col = "Group",
                              tp_col    = "timepoint"),   # no tp in agg
            error = function(e) .cld_agg
          )
          p_agg <- tryCatch(
            add_cld_to_scatter(p_agg,
                               y_data       = d_plot[[ind]],
                               cld_df       = .cld_agg,
                               x_col        = "Group",
                               grp_col      = "Treatment",
                               tp_col       = NULL,
                               tp_values    = NULL,
                               dodge_width  = 0.4,
                               font_size    = 4.5,
                               fixed_colour = "black"),
            error = function(e) p_agg
          )
        }
      }

      safe_nm  <- gsub("[^A-Za-z0-9_]", "_", ind)
      png_file <- file.path(png_dir_agg, paste0(safe_nm, "_aggregated.png"))
      tryCatch(
        ggplot2::ggsave(png_file, plot = p_agg, width = 8, height = 5.5,
                        dpi = 150, bg = "white"),
        error = function(e) warning("Aggregated PNG save failed for ", ind, ": ",
                                    conditionMessage(e), call. = FALSE)
      )
    }
    message("Aggregated graphs written to: ", png_dir_agg)

    # -- Aggregated Markdown report -------------------------------------------
    tryCatch(
      write_report_markdown_aggregated(trial_agg_results, OUTDIR),
      error = function(e) warning("Aggregated MD report failed: ", conditionMessage(e), call. = FALSE)
    )
  }
}, error = function(e) {
  warning("Trial-aggregated branch failed: ", conditionMessage(e), call. = FALSE)
})

# OPT-2: Reset parallel plan to sequential so STEP6 starts with a clean slate.
if (n_workers_step5 > 0L) {
  future::plan(future::sequential())
  message("[OPT-2] Parallel plan reset to sequential after STEP5.")
}

message("Pipeline completed.")
message("Outputs saved to: ", OUTDIR)

###############################################################################
# Final objects in session
###############################################################################
# fish_activity_summary_clean
# qc
# all_branch_results
# OUTDIR
###############################################################################