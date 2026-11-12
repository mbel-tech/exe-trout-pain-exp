# =============================================================================
# cld_helpers.R
# CLD (Compact Letter Display) utilities for graph annotation and report tables
#
# Sourced by 00_run_trajectory_pipeline.R after pretty_labels.R.
# Provides:
#   remap_cld_letters()         – add cld_remapped column to a CLD data.frame
#   all_cld_same()              – TRUE when all groups share one letter
#   .maybe_expand_timepoint()   – replicate CLD rows across timepoints
#   extract_cld_for_plot()      – pick + expand + remap CLD from posthoc result
#   add_cld_to_scatter()        – annotate a ggplot scatter with CLD letters
#   build_remapped_report_tbl() – wide report table: original + remapped + order
#   compute_agg_cld_df()        – CLD from an aggregated lm (STEP5 agg graphs)
# =============================================================================

# --- Internal display-order constants ----------------------------------------
.CLD_TREAT_ORDER <- c("Control", "Exercise")
.CLD_GRP_ORDER   <- c("No choice", "Exercise choice")
.CLD_TP_ORDER    <- c("2", "3", "1")   # interval 2 first per spec

# --- Internal: normalise the CLD letter column name -------------------------
# multcomp::cld() produces a `.group` column; the pipeline may also use
# `cld_letter`.  This helper ensures `cld_letter` is always present.
.norm_cld_col <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  if (!"cld_letter" %in% names(df) && ".group" %in% names(df)) {
    names(df)[names(df) == ".group"] <- "cld_letter"
    df$cld_letter <- trimws(as.character(df$cld_letter))
  }
  df
}

# =============================================================================
# remap_cld_letters
# =============================================================================
# Add a `cld_remapped` column to cld_df by remapping letters so the first
# displayed group (closest to y-axis, leftmost Treatment × Group × timepoint)
# always receives letter 'a'.
#
# Algorithm: walk rows in display-order (Treatment → Group → timepoint).
# Build an old-character → new-character mapping from first appearances.
# Compound labels (e.g. "ab") are handled character-by-character; new chars
# are sorted alphabetically before being joined.
#
# Args:
#   cld_df      data.frame with at least a `cld_letter` column
#   treat_col   Treatment column name (default "Treatment")
#   grp_col     Group column name     (default "Group")
#   tp_col      timepoint column name (default "timepoint")
#   treat_order character vector of Treatment levels in display order
#   grp_order   character vector of Group levels in display order
#   tp_order    character vector of timepoint levels in display order
#
# Returns: cld_df with new column `cld_remapped` (original column unchanged)
remap_cld_letters <- function(cld_df,
                               treat_col   = "Treatment",
                               grp_col     = "Group",
                               tp_col      = "timepoint",
                               treat_order = .CLD_TREAT_ORDER,
                               grp_order   = .CLD_GRP_ORDER,
                               tp_order    = .CLD_TP_ORDER) {
  cld_df <- .norm_cld_col(cld_df)
  if (is.null(cld_df) || !is.data.frame(cld_df) || nrow(cld_df) == 0L ||
      !"cld_letter" %in% names(cld_df)) {
    return(cld_df)
  }

  # Build numeric sort keys (unrecognised levels get 9999 → sorted last)
  sort_treat <- if (treat_col %in% names(cld_df)) {
    match(as.character(cld_df[[treat_col]]), treat_order, nomatch = 9999L)
  } else rep(1L, nrow(cld_df))

  sort_grp <- if (grp_col %in% names(cld_df)) {
    match(as.character(cld_df[[grp_col]]), grp_order, nomatch = 9999L)
  } else rep(1L, nrow(cld_df))

  sort_tp <- if (tp_col %in% names(cld_df)) {
    match(as.character(cld_df[[tp_col]]), tp_order, nomatch = 9999L)
  } else rep(1L, nrow(cld_df))

  row_order <- order(sort_treat, sort_grp, sort_tp)

  # Infinite letter supply: a, b, ..., z, aa, ab, ...
  new_letters <- c(letters,
                   outer(letters, letters, paste0)[seq_len(26L * 26L)])
  letter_map  <- character(0L)   # named vector: old_char -> new_char
  next_idx    <- 1L

  for (i in row_order) {
    old_str <- trimws(as.character(cld_df$cld_letter[i]))
    if (is.na(old_str) || !nzchar(old_str)) next
    for (ch in strsplit(old_str, "")[[1L]]) {
      if (!ch %in% names(letter_map) && next_idx <= length(new_letters)) {
        letter_map[[ch]] <- new_letters[[next_idx]]
        next_idx <- next_idx + 1L
      }
    }
  }

  cld_df$cld_remapped <- vapply(
    trimws(as.character(cld_df$cld_letter)),
    function(old_str) {
      if (is.na(old_str) || !nzchar(old_str)) return(old_str)
      chars     <- strsplit(old_str, "")[[1L]]
      new_chars <- vapply(chars,
                          function(ch) {
                            if (ch %in% names(letter_map)) letter_map[[ch]] else ch
                          },
                          character(1L))
      paste(sort(new_chars), collapse = "")
    },
    character(1L)
  )

  cld_df
}

# =============================================================================
# all_cld_same
# =============================================================================
# Return TRUE when all groups share the same CLD letter — in that case the
# annotation conveys no information and should be suppressed on the graph.
#
# Args:
#   cld_df     data.frame with a cld_letter or cld_remapped column, or NULL
#   letter_col column to inspect; NULL → auto-detect (prefers cld_remapped)
#
# Returns: logical scalar
all_cld_same <- function(cld_df, letter_col = NULL) {
  if (is.null(cld_df) || !is.data.frame(cld_df) || nrow(cld_df) == 0L)
    return(TRUE)
  cld_df <- .norm_cld_col(cld_df)
  if (is.null(letter_col))
    letter_col <- if ("cld_remapped" %in% names(cld_df)) "cld_remapped" else "cld_letter"
  if (!letter_col %in% names(cld_df))
    return(TRUE)
  vals <- trimws(as.character(cld_df[[letter_col]]))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0L) return(TRUE)
  length(unique(vals)) == 1L
}

# =============================================================================
# .maybe_expand_timepoint (internal)
# =============================================================================
# If cld_df has no timepoint column but the plot facets by timepoint,
# replicate cld_df rows once per timepoint value so geom_text can match facets.
# Returns cld_df unchanged when the column already exists.
.maybe_expand_timepoint <- function(cld_df, tp_values, tp_col = "timepoint") {
  if (is.null(cld_df) || !is.data.frame(cld_df) || nrow(cld_df) == 0L)
    return(cld_df)
  if (tp_col %in% names(cld_df))
    return(cld_df)
  if (is.null(tp_values) || length(tp_values) == 0L)
    return(cld_df)
  rows <- lapply(as.character(tp_values), function(tp) {
    d <- cld_df
    d[[tp_col]] <- tp
    d
  })
  do.call(rbind, rows)
}

# =============================================================================
# extract_cld_for_plot
# =============================================================================
# Extract the most informative CLD table from a run_posthoc_gd() result.
# Picks among cld_tables in preference order:
#   timepoint:Treatment:Group > Treatment:Group > any other valid table
# Expands the timepoint dimension when tp_values are supplied but the picked
# table has no timepoint column (.maybe_expand_timepoint).
# Adds cld_remapped via remap_cld_letters().
#
# Args:
#   post_result  list from run_posthoc_gd(), or NULL
#   tp_values    character/integer vector of timepoints present in plot data
#   ...          forwarded to remap_cld_letters()
#
# Returns: data.frame ready for add_cld_to_scatter(), or NULL
extract_cld_for_plot <- function(post_result, tp_values = NULL, ...) {
  if (is.null(post_result)) return(NULL)

  # Gather candidate tables
  cld_tbls <- if (!is.null(post_result$cld_tables) &&
                  length(post_result$cld_tables) > 0L) {
    post_result$cld_tables
  } else if (!is.null(post_result$cld) && is.data.frame(post_result$cld)) {
    list("cld" = post_result$cld)
  } else {
    return(NULL)
  }

  # Helper: is this table non-empty and has a usable CLD letter column?
  # Accepts both 'cld_letter' (pipeline convention) and '.group' (multcomp default).
  .valid_tbl <- function(tbl) {
    if (!is.data.frame(tbl) || nrow(tbl) == 0L) return(FALSE)
    tbl <- .norm_cld_col(tbl)
    "cld_letter" %in% names(tbl) &&
      !all(is.na(tbl$cld_letter) | !nzchar(trimws(as.character(tbl$cld_letter))))
  }

  # Preferred term order for graph annotation.
  # Includes all common orderings produced by run_posthoc_gd().
  pref_terms <- c("timepoint:Treatment:Group",
                  "Treatment:Group:timepoint",
                  "Treatment:Group",
                  "Group:Treatment",
                  "timepoint:Treatment",
                  "Treatment:timepoint",
                  "Group:timepoint",
                  "timepoint:Group",
                  "Treatment",
                  "Group",
                  "timepoint")

  # Two-pass selection:
  #   Pass 1 – prefer a table whose letters are NOT all the same (informative)
  #   Pass 2 – accept any valid table (may be all-same → suppression handles it)
  pick <- NULL
  for (term in pref_terms) {
    if (term %in% names(cld_tbls)) {
      tbl_cand <- cld_tbls[[term]]
      if (.valid_tbl(tbl_cand) &&
          !all_cld_same(.norm_cld_col(tbl_cand))) {
        pick <- tbl_cand
        break
      }
    }
  }
  if (is.null(pick)) {
    for (term in pref_terms) {
      if (term %in% names(cld_tbls) && .valid_tbl(cld_tbls[[term]])) {
        pick <- cld_tbls[[term]]
        break
      }
    }
  }
  # Fallback: first valid table in list order
  if (is.null(pick)) {
    for (tbl in cld_tbls) {
      if (.valid_tbl(tbl)) { pick <- tbl; break }
    }
  }

  if (is.null(pick)) return(NULL)

  # Normalise letter column name (.group → cld_letter)
  pick <- .norm_cld_col(pick)

  # Expand timepoint rows if needed
  if (!is.null(tp_values))
    pick <- .maybe_expand_timepoint(pick, tp_values, tp_col = "timepoint")

  # Remap letters
  pick <- remap_cld_letters(pick, ...)
  pick
}

# =============================================================================
# add_cld_to_scatter
# =============================================================================
# Add CLD letter annotations to a ggplot scatter (make_scatter_gd() output).
#
# All letters share one y-height = max(raw data, error bar tops) + 10% range.
# Letters are suppressed (plot returned unchanged) when all groups have the
# same letter.
#
# The cld_df must already have been processed by extract_cld_for_plot() so that
# it contains a `cld_remapped` column and only the timepoints present in the
# plot data (pass tp_values to extract_cld_for_plot() upstream to ensure this).
#
# Args:
#   p            ggplot object (from make_scatter_gd())
#   y_data       numeric vector of finite raw y values in the plot
#   cld_df       data.frame from extract_cld_for_plot()
#   summ_df      data.frame with m (mean) and s (SEM) for error-bar-top inclusion
#   x_col        column mapped to x aesthetic   (default "Treatment")
#   grp_col      column used for position_dodge  (default "Group")
#   tp_col       facet column                    (default "timepoint")
#   tp_values    timepoints present in plot data; cld_df rows outside this set
#                are dropped so no extra empty facets are created
#   dodge_width  must match the scatter's dodge_width (default 0.75)
#   font_size    CLD letter text size in ggplot units (default 4.5)
#   letter_col   column to display; NULL → cld_remapped or cld_letter
#   fixed_colour fixed colour string for letters (NULL = colour by grp_col)
#
# Returns: modified ggplot object
add_cld_to_scatter <- function(p, y_data, cld_df, summ_df = NULL,
                                x_col        = "Treatment",
                                grp_col      = "Group",
                                tp_col       = "timepoint",
                                tp_values    = NULL,
                                dodge_width  = 0.75,
                                font_size    = 4.5,
                                letter_col   = NULL,
                                fixed_colour = NULL) {
  # Guard: nothing to annotate
  if (is.null(cld_df) || !is.data.frame(cld_df) || nrow(cld_df) == 0L)
    return(p)

  # Determine letter column
  if (is.null(letter_col))
    letter_col <- if ("cld_remapped" %in% names(cld_df)) "cld_remapped" else "cld_letter"
  if (!letter_col %in% names(cld_df)) return(p)

  # Suppression: skip when all groups share one letter
  if (all_cld_same(cld_df, letter_col)) return(p)

  # ---------- Expand missing x_col / grp_col by replicating rows ------------
  # A CLD table may cover fewer dimensions than the scatter (e.g. Group×Interval
  # CLD used on a Treatment×Group×Interval plot).  When Treatment or Group is
  # absent, replicate rows for each level so position_dodge places letters
  # correctly over every Treatment (or Group) position in the plot.

  # Helper: levels actually present in the plot, derived from summ_df if possible,
  # then from a global ORDER vector (filtered to data), then a hard default.
  .plot_levels <- function(col_nm, summ_df, fallback_obj, fallback_default, data_col = NULL) {
    # 1. From summ_df (most accurate — reflects post-filter data)
    if (!is.null(summ_df) && is.data.frame(summ_df) && col_nm %in% names(summ_df))
      return(unique(as.character(summ_df[[col_nm]])))
    # 2. From y_data companion column (data_col not usually available here; skip)
    # 3. From global ORDER vector, but preserve only levels seen in cld_df row context
    if (!is.null(fallback_obj))
      return(fallback_obj)
    fallback_default
  }

  if (!x_col %in% names(cld_df)) {
    x_fb <- if (exists("TREATMENT_ORDER", envir = .GlobalEnv, inherits = TRUE) &&
                x_col == "Treatment") {
      get("TREATMENT_ORDER", envir = .GlobalEnv)
    } else if (exists("GROUP_ORDER", envir = .GlobalEnv, inherits = TRUE) &&
               x_col == "Group") {
      get("GROUP_ORDER", envir = .GlobalEnv)
    } else NULL
    x_levels <- .plot_levels(x_col, summ_df, x_fb, c("Control", "Exercise"))
    expanded <- lapply(x_levels, function(lv) { d <- cld_df; d[[x_col]] <- lv; d })
    cld_df <- do.call(rbind, expanded)
  }

  if (!grp_col %in% names(cld_df)) {
    g_fb <- if (exists("GROUP_ORDER", envir = .GlobalEnv, inherits = TRUE) &&
                grp_col == "Group") {
      get("GROUP_ORDER", envir = .GlobalEnv)
    } else if (exists("TREATMENT_ORDER", envir = .GlobalEnv, inherits = TRUE) &&
               grp_col == "Treatment") {
      get("TREATMENT_ORDER", envir = .GlobalEnv)
    } else NULL
    g_levels <- .plot_levels(grp_col, summ_df, g_fb, c("No choice", "Exercise choice"))
    expanded <- lapply(g_levels, function(lv) { d <- cld_df; d[[grp_col]] <- lv; d })
    cld_df <- do.call(rbind, expanded)
  }

  # Required x and group columns must now exist
  if (!x_col %in% names(cld_df) || !grp_col %in% names(cld_df)) return(p)

  # ---------- Compute shared y-annotation height ----------------------------
  y_finite <- y_data[is.finite(y_data)]
  if (length(y_finite) == 0L) return(p)

  y_max_data <- max(y_finite)
  y_min_data <- min(y_finite)

  # Include error bar tops if summ_df provided
  if (!is.null(summ_df) && is.data.frame(summ_df) &&
      all(c("m", "s") %in% names(summ_df))) {
    eb_tops    <- summ_df$m + summ_df$s
    y_max_data <- max(c(y_max_data, eb_tops[is.finite(eb_tops)]))
  }

  y_range <- y_max_data - y_min_data
  y_cld   <- y_max_data + 0.10 * max(y_range, abs(y_max_data) * 0.01, 1e-6)

  # Extra expansion so the letter itself fits inside the panel
  y_expand <- y_cld + 0.05 * max(y_range, abs(y_max_data) * 0.01, 1e-6)

  # ---------- Build annotation data frame -----------------------------------
  ann_df <- cld_df
  ann_df$.cld_label <- trimws(as.character(ann_df[[letter_col]]))
  ann_df <- ann_df[!is.na(ann_df$.cld_label) & nzchar(ann_df$.cld_label), , drop = FALSE]
  if (nrow(ann_df) == 0L) return(p)

  # Filter to timepoints present in plot data (prevents empty facets)
  if (!is.null(tp_col) && tp_col %in% names(ann_df) && !is.null(tp_values)) {
    ann_df <- ann_df[as.character(ann_df[[tp_col]]) %in% as.character(tp_values), ,
                     drop = FALSE]
    if (nrow(ann_df) == 0L) return(p)
  }

  ann_df$.y_cld <- y_cld

  # Coerce x and group columns to factors (same levels as plot data)
  if (exists("TREATMENT_ORDER", envir = .GlobalEnv, inherits = TRUE)) {
    trt_ord <- get("TREATMENT_ORDER", envir = .GlobalEnv)
    if (x_col == "Treatment") ann_df[[x_col]] <- factor(ann_df[[x_col]], levels = trt_ord)
    if (grp_col == "Treatment") ann_df[[grp_col]] <- factor(ann_df[[grp_col]], levels = trt_ord)
  }
  if (exists("GROUP_ORDER", envir = .GlobalEnv, inherits = TRUE)) {
    grp_ord <- get("GROUP_ORDER", envir = .GlobalEnv)
    if (x_col == "Group") ann_df[[x_col]] <- factor(ann_df[[x_col]], levels = grp_ord)
    if (grp_col == "Group") ann_df[[grp_col]] <- factor(ann_df[[grp_col]], levels = grp_ord)
  }

  # ---------- Add geom_text layer -------------------------------------------
  pos_dodge <- ggplot2::position_dodge(width = dodge_width)

  if (!is.null(fixed_colour)) {
    # Fixed colour (e.g. "black") — no colour aesthetic
    p <- p +
      ggplot2::geom_text(
        data = ann_df,
        ggplot2::aes(x     = .data[[x_col]],
                     y     = .data$.y_cld,
                     label = .data$.cld_label,
                     group = .data[[grp_col]]),
        position    = pos_dodge,
        colour      = fixed_colour,
        size        = font_size,
        fontface    = "bold",
        show.legend = FALSE,
        inherit.aes = FALSE
      )
  } else {
    # Colour follows the group dodge variable (matches data point colours)
    p <- p +
      ggplot2::geom_text(
        data = ann_df,
        ggplot2::aes(x      = .data[[x_col]],
                     y      = .data$.y_cld,
                     label  = .data$.cld_label,
                     group  = .data[[grp_col]],
                     colour = .data[[grp_col]]),
        position    = pos_dodge,
        size        = font_size,
        fontface    = "bold",
        show.legend = FALSE,
        inherit.aes = FALSE
      )
  }

  p <- p + ggplot2::expand_limits(y = y_expand)
  p
}

# =============================================================================
# build_remapped_report_tbl
# =============================================================================
# Build a compact report table that documents the letter remapping used in
# graphs.  Columns: grouping columns + "CLD (original)" + "CLD (remapped)" +
# "Remapping order" (description in first row, blank elsewhere).
#
# Args:
#   cld_df    data.frame from extract_cld_for_plot() or remap_cld_letters()
#             Must contain `cld_letter`; `cld_remapped` is added if absent.
#   tp_col    timepoint column name (default "timepoint")
#   treat_col Treatment column name (default "Treatment")
#   grp_col   Group column name     (default "Group")
#
# Returns: data.frame or NULL (if input is invalid or has no cld_letter column)
build_remapped_report_tbl <- function(cld_df,
                                       tp_col    = "timepoint",
                                       treat_col = "Treatment",
                                       grp_col   = "Group") {
  cld_df <- .norm_cld_col(cld_df)
  if (is.null(cld_df) || !is.data.frame(cld_df) || nrow(cld_df) == 0L ||
      !"cld_letter" %in% names(cld_df))
    return(NULL)

  # Ensure remapped column exists
  if (!"cld_remapped" %in% names(cld_df))
    cld_df <- remap_cld_letters(cld_df,
                                 treat_col = treat_col,
                                 grp_col   = grp_col,
                                 tp_col    = tp_col)

  # Grouping columns present in df (preserve column order)
  grp_cols <- intersect(c(tp_col, treat_col, grp_col), names(cld_df))

  out <- cld_df[, grp_cols, drop = FALSE]
  out[["CLD (original)"]] <- trimws(as.character(cld_df$cld_letter))
  out[["CLD (remapped)"]] <- trimws(as.character(cld_df$cld_remapped))

  # Remapping-order description: first row only, rest blank
  tp_part <- if (any(.CLD_TP_ORDER %in% c("2", "3")))
    paste0(" > ", paste0("Interval ", .CLD_TP_ORDER[.CLD_TP_ORDER %in% c("2", "3")],
                         collapse = " > "))
  else ""
  remap_order_str <- paste0(
    "Display order: ",
    paste(.CLD_TREAT_ORDER, collapse = " > "),
    " > ",
    paste(.CLD_GRP_ORDER, collapse = " > "),
    tp_part
  )
  out[["Remapping order"]] <- c(remap_order_str, rep("", nrow(out) - 1L))

  out
}

# =============================================================================
# compute_agg_cld_df
# =============================================================================
# Compute a CLD data.frame directly from an aggregated lm model using emmeans
# and multcomp::cld().  Used in the STEP5 trial-aggregated graph loop where
# CLD is not pre-computed and stored in results.
#
# Args:
#   model      lm (or lmer/glmmTMB) model object, or NULL
#   term       emmeans specification formula (default ~ Treatment * Group)
#   adjust     p-value adjustment method passed to multcomp::cld()
#   treat_col  expected Treatment column name (default "Treatment")
#   grp_col    expected Group column name     (default "Group")
#
# Returns: data.frame with Treatment, Group, cld_letter, or NULL on failure
compute_agg_cld_df <- function(model,
                                term      = ~ Treatment * Group,
                                adjust    = "tukey",
                                treat_col = "Treatment",
                                grp_col   = "Group") {
  if (is.null(model)) return(NULL)
  tryCatch({
    em  <- emmeans::emmeans(model, term)
    cld_obj <- multcomp::cld(em, Letters = letters, adjust = adjust,
                              reversed = FALSE)
    df <- as.data.frame(cld_obj)
    # Rename .group → cld_letter
    if (".group" %in% names(df) && !"cld_letter" %in% names(df))
      names(df)[names(df) == ".group"] <- "cld_letter"
    if (!"cld_letter" %in% names(df)) return(NULL)
    df$cld_letter <- trimws(as.character(df$cld_letter))
    # Return only grouping columns + cld_letter
    keep <- c(treat_col, grp_col, "cld_letter")
    df[, intersect(keep, names(df)), drop = FALSE]
  }, error = function(e) {
    NULL
  })
}
