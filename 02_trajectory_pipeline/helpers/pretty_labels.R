# =============================================================================
# pretty_labels.R
# -----------------------------------------------------------------------------
# Purpose : Centralised, manuscript-friendly label definitions for the GOT
#           zebrafish behavioural analysis pipeline.
#
# Usage   : source("path/to/pretty_labels.R")
#           This file is intended to be sourced by STEP4, STEP5, and STEP6.
#
# Side-effects: NONE beyond defining named vectors and functions.
#               No packages are loaded here; the calling script is assumed to
#               have already loaded any required packages.
#
# Sections:
#   1. Indicator labels
#   2. Group labels
#   3. Treatment labels
#   4. Interval labels
#   5. ANOVA term labels
#   6. Analysis type labels
#   7. ANOVA caption builder
#   8. Colour and shape constants
# =============================================================================


# =============================================================================
# 1. INDICATOR LABELS
# =============================================================================
# Canonical pretty names for every behavioural/kinematic indicator.
# STEP5 versions are used as canonical (more descriptive than STEP4 variants).
# Additional entries cover STEP6 group-dynamic indicators.

GOT_INDICATOR_LABELS <- c(
  # --- Distance ---
  distance_m                 = "Distance travelled (m)",
  distance_BL                = "Distance travelled (body lengths)",

  # --- Speed (cm/s) ---
  speed_cm_s_mean            = "Mean speed (cm/s)",
  speed_cm_s_median          = "Median speed (cm/s)",
  speed_cm_s_max             = "Maximum speed (cm/s)",

  # --- Speed (body lengths/s) ---
  speed_BL_s_mean            = "Mean speed (body lengths/s)",
  speed_BL_s_median          = "Median speed (body lengths/s)",
  speed_BL_s_max             = "Maximum speed (body lengths/s)",

  # --- Acceleration ---
  accel_cm_s2_mean           = "Mean acceleration (cm/s²)",
  accel_cm_s2_median         = "Median acceleration (cm/s²)",
  accel_cm_s2_max            = "Maximum acceleration (cm/s²)",

  # --- Zone / activity ---
  n_zone_transitions         = "Number of zone transitions",
  time_moving_TS0.5_s        = "Time moving (TS0.5; s)",
  time_moving_TS1_s          = "Time moving (TS1; s)",
  time_spent_zone_flow_s     = "Time spent in flow zone (s)",

  # --- STEP6 group-dynamic indicators ---
  mean_iid_cm                = "Mean inter-individual distance (cm)",
  mean_iid_cm_fish           = "Mean inter-individual distance (cm)",
  mean_nnd_cm                = "Mean nearest-neighbour distance (cm)",
  zone_synchrony             = "Zone synchrony (fraction of frames)",
  mean_centroid_speed_cm_s   = "Group centroid speed (cm/s)",
  mean_hull_area_cm2         = "Convex hull area (cm²)",
  mean_polarization          = "Group polarization (0–1)",
  mean_ang_momentum          = "Angular momentum (0–1)"
)

# Returns the pretty label for a single indicator name.
# Falls back to gsub("_", " ", x) if the indicator is not in the lookup table.
pretty_indicator_name <- function(x) {
  x   <- as.character(x)[1]            # scalar input expected
  lab <- unname(GOT_INDICATOR_LABELS[x])
  if (length(lab) == 1L && !is.na(lab) && nzchar(lab)) return(lab)
  gsub("_", " ", x)
}

# Accessor — returns the full named vector (useful for ggplot scale_*_discrete).
get_indicator_labels <- function() GOT_INDICATOR_LABELS


# =============================================================================
# 2. GROUP LABELS
# =============================================================================
# Groups as they appear in the data.  Both capitalisation variants of
# "Exercise choice" / "Exercise Choice" are mapped to the canonical lower-c form.

GOT_GROUP_LABELS <- c(
  "Exercise choice" = "Exercise choice",
  "Exercise Choice" = "Exercise choice",   # normalise capitalisation
  "No choice"       = "No choice",
  "No Choice"       = "No choice",         # normalise capitalisation
  "Morphine"        = "Morphine",
  "Naloxone"        = "Naloxone"
)

# Vectorised: handles character vectors of any length.
# Unknown values are passed through unchanged.
pretty_group_name <- function(x) {
  x   <- as.character(x)
  lab <- unname(GOT_GROUP_LABELS[x])
  ifelse(!is.na(lab) & nzchar(lab), lab, x)
}

get_group_labels <- function() GOT_GROUP_LABELS


# =============================================================================
# 3. TREATMENT LABELS
# =============================================================================

GOT_TREATMENT_LABELS <- c(
  "Control"  = "Control",
  "Exercise" = "Exercise"
)

# Vectorised; unknown values are passed through unchanged.
pretty_treatment_name <- function(x) {
  x   <- as.character(x)
  lab <- unname(GOT_TREATMENT_LABELS[x])
  ifelse(!is.na(lab) & nzchar(lab), lab, x)
}

get_treatment_labels <- function() GOT_TREATMENT_LABELS


# =============================================================================
# 4. INTERVAL LABELS
# =============================================================================
# The internal column is "timepoint" (integer 1 / 2 / 3).
# User-facing label should read "Interval 1", "Interval 2", "Interval 3".
# This is the KEY CHANGE: "Timepoint" -> "Interval" throughout.

GOT_INTERVAL_LABELS <- c(
  "1" = "Interval 1",
  "2" = "Interval 2",
  "3" = "Interval 3"
)

# Vectorised; unknown values fall back to paste0("Interval ", x).
pretty_interval_name <- function(x) {
  x   <- as.character(x)
  lab <- unname(GOT_INTERVAL_LABELS[x])
  ifelse(!is.na(lab) & nzchar(lab), lab, paste0("Interval ", x))
}

get_interval_labels <- function() GOT_INTERVAL_LABELS


# =============================================================================
# 5. ANOVA TERM LABELS
# =============================================================================
# Converts raw model-term strings (as returned by lmerTest / car::Anova) into
# human-readable form:
#   - Interaction operator ":"  becomes " × "
#   - The column name "timepoint" becomes "Interval"
#   - Treatment and Group capitalisation is preserved as-is

pretty_term_name <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return(NA_character_)
  x <- as.character(x)
  x <- gsub(":", " × ", x)               # : -> × (Unicode multiplication sign)
  x <- gsub("\\btimepoint\\b", "Interval", x) # timepoint -> Interval (whole-word match)
  x
}


# =============================================================================
# 6. ANALYSIS TYPE LABELS
# =============================================================================

pretty_analysis_type <- function(type) {
  switch(type,
    "interval_level" = "Interval-Level Analysis",
    "aggregated"     = "Trial-Aggregated Analysis",
    as.character(type)   # fall-through: return the raw value
  )
}


# =============================================================================
# 7. ANOVA CAPTION BUILDER
# =============================================================================
# Builds a multi-line character string suitable for use as a plot/table caption.
#
# Design decisions:
#   - The FULL interaction term (highest-order term in the model) is ALWAYS
#     included, even when its p-value exceeds alpha.  This ensures readers can
#     always see the interaction test result without having to hunt for it.
#   - All additional terms with p < alpha are appended.
#   - "(Intercept)" is always excluded from output.
#   - Format for F-tests (LMM/lmerTest):  "Term: F(numDF, denDF) = stat, p = value"
#   - Format for chi-sq tests (GLMM/glmer): "Term: χ²(df) = stat, p = value"
#
# The anova_df is expected to have the following attributes set by
# STEP5's run_anova_table():
#   attr(anova_df, "stat_col")  — name of the F / Chisq column
#   attr(anova_df, "p_col")     — name of the p-value column
#   attr(anova_df, "df_cols")   — character vector of df column name(s)
# The data frame must also have a column named "term" containing term names.

build_anova_caption_multiline <- function(anova_df,
                                          indicator_label = NA_character_,
                                          alpha = 0.05,
                                          caption_format = c("ggtext", "markdown", "plain")) {

  caption_format <- match.arg(caption_format)

  # Guard: return empty string for NULL / non-data-frame / empty input
  if (is.null(anova_df) || !is.data.frame(anova_df) || nrow(anova_df) == 0L) {
    return("")
  }

  stat_col <- attr(anova_df, "stat_col")
  p_col    <- attr(anova_df, "p_col")
  df_cols  <- attr(anova_df, "df_cols")   # may be NULL for chi-sq models

  if (is.null(stat_col) || is.null(p_col)) return("")

  # Determine whether this is an F-test or a chi-squared test
  # (stat_col names typically start with "F" or "Chisq" / "Chi.sq" etc.)
  is_f_test <- grepl("^F", stat_col, ignore.case = FALSE)

  # ---- Identify the "full interaction" term to always include ----------------
  all_terms <- as.character(anova_df$term)
  all_terms <- all_terms[!is.na(all_terms) & trimws(all_terms) != "(Intercept)"]

  # Prefer 3-way, then 2-way Treatment:Group, then the term with most colons
  full_interaction <- if ("timepoint:Treatment:Group" %in% all_terms) {
    "timepoint:Treatment:Group"
  } else if ("Treatment:Group" %in% all_terms) {
    "Treatment:Group"
  } else if (length(all_terms) > 0L) {
    # Count colons as a proxy for interaction order
    colon_counts  <- nchar(gsub("[^:]", "", all_terms))
    all_terms[which.max(colon_counts)]
  } else {
    character(0)
  }

  # ---- Collect significant terms (excluding Intercept) ----------------------
  is_intercept <- trimws(as.character(anova_df$term)) == "(Intercept)"
  p_vals       <- suppressWarnings(as.numeric(anova_df[[p_col]]))
  sig_mask     <- !is.na(p_vals) & p_vals < alpha & !is_intercept

  sig_terms <- as.character(anova_df$term)[sig_mask]

  # Union: always include full_interaction, then any significant terms
  terms_to_show <- unique(c(full_interaction, sig_terms))

  # Restore original row order from anova_df
  terms_to_show <- intersect(as.character(anova_df$term), terms_to_show)

  x <- anova_df[as.character(anova_df$term) %in% terms_to_show, , drop = FALSE]
  if (nrow(x) == 0L) return("")

  # ---- Internal formatting helpers ------------------------------------------

  # Round to `digits` decimal places, strip trailing zeros and lone decimal point.
  # Uses fixed-point format so that values like 0.00123 are shown as "0.001"
  # (max 3 decimal places) rather than "0.00123" (3 significant figures).
  fmt_stat <- function(val, digits = 3L) {
    val <- as.numeric(val)
    if (is.na(val)) return("NA")
    formatted <- formatC(round(val, digits), format = "f", digits = digits)
    # Remove trailing zeros after decimal and lonely decimal point
    formatted <- sub("\\.$", "", sub("0+$", "", trimws(formatted)))
    trimws(formatted)
  }

  # Format a p-value: very small -> "< 0.0001", otherwise 3 decimal places
  fmt_p_val <- function(p) {
    p <- as.numeric(p)
    if (is.na(p)) return("NA")
    if (p < 0.0001) return("< 0.0001")
    formatC(p, digits = 3L, format = "f")
  }

  # ---- Build one line per selected term -------------------------------------
  pieces <- apply(x, 1L, function(row_i) {
    # apply() coerces data.frame to matrix (all values become character strings).
    # as.list() preserves column names; as.numeric() below restores numeric type.
    row_i  <- as.list(row_i)
    term_i <- pretty_term_name(as.character(row_i$term))
    stat_i <- as.numeric(row_i[[stat_col]])
    p_i    <- as.numeric(row_i[[p_col]])

    # Degrees of freedom: collect from df_cols (may be 1 or 2 columns).
    # All df are rounded to nearest integer (handles fractional Satterthwaite df).
    df_txt <- ""
    if (!is.null(df_cols) && length(df_cols) > 0L) {
      df_vals <- unlist(row_i[df_cols])
      df_vals <- df_vals[!is.na(df_vals)]
      if (length(df_vals) > 0L) {
        df_fmt <- vapply(df_vals, function(v) {
          v <- as.numeric(v)
          if (is.na(v)) return("NA")
          as.character(as.integer(round(v)))   # always integer, no decimal
        }, character(1L))
        df_txt <- paste(df_fmt, collapse = ",")   # no spaces inside subscript
      }
    }
    if (!nzchar(df_txt)) df_txt <- "NA"

    stat_fmt <- fmt_stat(stat_i)
    p_fmt    <- fmt_p_val(p_i)
    p_part   <- if (startsWith(p_fmt, "<")) {
      paste0("p ", p_fmt)
    } else {
      paste0("p = ", p_fmt)
    }

    if (is_f_test) {
      stat_str <- switch(caption_format,
        ggtext   = paste0("F<sub>", df_txt, "</sub> = ", stat_fmt),
        markdown = paste0("F~", df_txt, "~ = ", stat_fmt),
        plain    = paste0("F(", gsub(",", ", ", df_txt), ") = ", stat_fmt)
      )
    } else {
      # χ² test — df is a single value for chi-sq
      stat_str <- switch(caption_format,
        ggtext   = paste0("χ²<sub>", df_txt, "</sub> = ", stat_fmt),
        markdown = paste0("χ²~", df_txt, "~ = ", stat_fmt),
        plain    = paste0("χ²(", df_txt, ") = ", stat_fmt)
      )
    }

    paste0(term_i, ": ", stat_str, ", ", p_part)
  })

  paste(pieces, collapse = "\n")
}


# =============================================================================
# 8. COLOUR AND SHAPE CONSTANTS
# =============================================================================
# Canonical colours and shapes for all interval × treatment × group graphs.
# Keys are LOWER-CASE to allow case-insensitive matching.

GOT_GROUP_COLORS <- c(
  "no choice"       = "#56B4E9",   # blue
  "exercise choice" = "#D55E00"    # vermillion
)

GOT_TREATMENT_SHAPES <- c(
  "control"  = 17L,   # filled triangle
  "exercise" = 16L    # filled circle
)

# Canonical factor orders (for scale_x_discrete / forcats::fct_relevel)
GOT_GROUP_ORDER     <- c("No choice", "Exercise choice")
GOT_TREATMENT_ORDER <- c("Control", "Exercise")
GOT_INTERVAL_ORDER  <- c("1", "2", "3")   # internal column values

# Accessor helpers
get_group_colors     <- function() GOT_GROUP_COLORS
get_treatment_shapes <- function() GOT_TREATMENT_SHAPES
get_group_order      <- function() GOT_GROUP_ORDER
get_treatment_order  <- function() GOT_TREATMENT_ORDER
get_interval_order   <- function() GOT_INTERVAL_ORDER
