# =============================================================================
# house_style.R -- CANONICAL source of truth for GOT palette, theme,
# marks, CLD, number formatting, effect sizes, and export.
# =============================================================================
# Shared palette, theme, marks, compact letter display, number formatting,
# effect sizes and export for the 2x2 Treatment x Group design (4 cells).
# Comments cite numbered rules (HANDOFF S*, CONVENTIONS C*, RESULTS_STYLE R*)
# of the internal analysis plan; that plan is not part of this archive.
#
# DECISION (figure look): GOT keeps its own palette/theme identity rather than
# adopting CHOICE's -- only rule VIOLATIONS are fixed here (statistical ink
# black, one CLD_SIZE, one sourced theme object, bg="white" on PNG export).
# Values below are lifted verbatim from figure_3_locomotion.R where noted, so this
# file is a refactor (single source of truth), not a re-design.
#
# DO NOT RETYPE any value in this file into another script -- source it:
#   source(file.path(GOT_ROOT, "05_statistics/shared/house_style.R"))
# =============================================================================

suppressPackageStartupMessages({
  .HAS_EFFECTSIZE <- requireNamespace("effectsize", quietly = TRUE)
  .HAS_GGTEXT     <- requireNamespace("ggtext",     quietly = TRUE)
})

# =============================================================================
# GOT_ROOT -- single path anchor (HANDOFF S11: hardcoded absolute paths
# scattered across dozens of files is a named repo fragility; this constant is
# the one place a future move of the repo has to be edited).
# =============================================================================
if (!exists("GOT_ROOT")) GOT_ROOT <- "."

# =============================================================================
# A1 -- PALETTE  (CONVENTIONS.md C1, C2)
# =============================================================================
# GOT crosses Treatment x Group (2x2), so colour carries the crossed cell
# identity and shape carries Treatment alone -- NOT the C1 [DR] 6/5/3/1
# single-axis sequence, which assumes one unordered categorical variable.
# This crossed encoding is a GOT-specific [DR], recorded here rather than
# silently deviating from C1's worked example.
#
# Values lifted verbatim from figure_3_locomotion.R (this IS GOT's declared house
# palette -- "keep GOT's look" decision, 2026-08-24). Statistical ink is
# black, same as CHOICE's INK_STAT -- that part of C1 is a real invariant.

TG_COLS <- c(
  "Exercise:Exercise choice" = "#D55E00",
  "Exercise:No choice"       = "#56B4E9",
  "Control:Exercise choice"  = "#E69F00",
  "Control:No choice"        = "#0072B2"
)
TG_LBLS <- c(
  "Control:Exercise choice"  = "Control / Exercise-choice",
  "Control:No choice"        = "Control / No-choice",
  "Exercise:Exercise choice" = "Exercise / Exercise-choice",
  "Exercise:No choice"       = "Exercise / No-choice"
)
# Shape by Group: 17 = No choice (triangle), 19 = Exercise choice (circle)
GRP_SHPS <- c("Exercise choice" = 19, "No choice" = 17)
TG_SHPS <- c(
  "Control:Exercise choice"  = 19,
  "Control:No choice"        = 17,
  "Exercise:Exercise choice" = 19,
  "Exercise:No choice"       = 17
)
# Shape by Treatment: 17 = Control (triangle), 19 = Exercise (circle)
TRT_SHPS <- c("Control" = 17, "Exercise" = 19)

# Canonical cell order (dodge/legend/CLD anchor order). Control-first within
# each Group, No-choice-first between Groups -- matches figure_3_locomotion.R's
# TG_BREAKS_PLOT / AGG_INT_TG_ORDER.
TG_BREAKS_PLOT <- c(
  "Control:No choice",
  "Exercise:No choice",
  "Control:Exercise choice",
  "Exercise:Exercise choice"
)

GRP_LVLS <- c("No choice", "Exercise choice")
TRT_LVLS <- c("Control", "Exercise")

INK_STAT <- "black"   # every mean, error bar, CLD letter, bracket -- C1 [INV]

SCALE_COL <- ggplot2::scale_colour_manual(
  name   = "Treatment : Group",
  values = TG_COLS,
  breaks = TG_BREAKS_PLOT,
  labels = TG_LBLS[TG_BREAKS_PLOT]
)
SCALE_SHP_GRP <- ggplot2::scale_shape_manual(
  name   = "Group",
  values = GRP_SHPS,
  breaks = GRP_LVLS,
  guide  = "none"
)

# =============================================================================
# A2 -- THEME  (CONVENTIONS.md C3, C8/A-5)
# =============================================================================
# theme_classic(base_size = 10) is GOT's declared base (kept, per the "keep
# GOT's look" decision) -- moved here from figure_3_locomotion.R unchanged EXCEPT
# plot.caption: colour grey25 -> black and italic added, so statistical ink is
# black everywhere (C1 [INV]) and caption formatting matches C6.
#
# FONT_CAP + cairo_pdf are load-bearing: the ANOVA captions contain chi, superscript
# 2, subscript digits and x. The default grDevices pdf() device cannot encode
# these (mbcsToSbcs failures) -- confirmed as a SILENT corruption failure mode,
# not just tofu glyphs. Keep both.
FONT_CAP <- "Segoe UI Symbol"

BASE_THEME <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title.x   = ggplot2::element_blank(),
    axis.title.y   = ggplot2::element_text(size = 9.5),
    axis.text      = ggplot2::element_text(size = 8.5, colour = "black"),
    axis.line      = ggplot2::element_line(linewidth = 0.4),
    plot.caption   = ggplot2::element_text(size = 7.5, hjust = 0, lineheight = 1.3,
                                           margin = ggplot2::margin(t = 6),
                                           colour = INK_STAT, face = "italic",
                                           family = FONT_CAP),
    legend.title   = ggplot2::element_text(size = 8.5),
    legend.text    = ggplot2::element_text(size = 8),
    legend.key.size = grid::unit(0.4, "cm"),
    plot.tag       = ggplot2::element_text(size = 10, face = "bold")
  )

# One theme dialect for the whole chain (C8 A-5): theme_bw()/theme_minimal()
# local definitions in the aggression scripts are retired in favour of this
# object.
theme_aggression <- function(base_size = 12) {
  BASE_THEME + ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
    panel.grid        = ggplot2::element_blank(),
    legend.position    = "bottom",
    plot.title         = ggplot2::element_text(face = "bold", size = base_size + 1),
    plot.subtitle      = ggplot2::element_text(size = base_size - 1, colour = "grey40")
  )
}

# =============================================================================
# A3/A4 -- DATA MARKS AND GEOMETRY  (CONVENTIONS.md C4, C5/A-4)
# =============================================================================
# Geometry constants lifted verbatim from figure_3_locomotion.R -- figure sizes do
# not change under this refactor.
TG_DODGE      <- 0.38
SUMMARY_WIDTH <- 0.22
POINT_SIZE    <- 2.8

FIG45_DODGE            <- TG_DODGE * 3
FIG45_SUMMARY_WIDTH    <- SUMMARY_WIDTH * 2
FIG45_JITTER_WIDTH      <- SUMMARY_WIDTH * 2
FIG45_X_TEXT_SIZE       <- 8.5 * 3 * 0.7 * 0.5
FIG45_X_EXPAND_ADD      <- 0.65 * 2
FIG45_X_EXPAND_LOW      <- 0.85
FIG45_X_EXPAND_HIGH     <- 0.85
FIG45_TAG_SIZE          <- 15
FIG45_CAPTION_MARGIN_T  <- 24
PANEL5B_YLAB_SIZE       <- 9.5 * 2 * 0.7
INT_DODGE               <- 0.58 * 1.7
X_TICK_SIZE             <- 8.5 * 2.5 * 0.7
CAPTION_MARGIN_T        <- 24
CAPTION_SIZE_45         <- 7.5 * 2
FIG7_CAPTION_SIZE       <- 7.5 * 1.5 * 2
FIG7_TAG_SIZE           <- 10 * 2.5
FIG7_TAG_Y              <- 0.985 + 0.0175

# ADJUDICATION (GOT A-4 analogue): a single live CLD_SIZE, replacing the
# `4.2 * 1.3 * 0.7 * 1.5` literal at figure_3_locomotion.R:25. Numeric value (5.733)
# preserved so existing figures do not change size under this refactor; the
# handoff's own CHOICE value is 6 (house_style.R A3/A4) -- noted here so a
# future harmonisation pass knows both numbers exist and only one is live.
CLD_SIZE <- 4.2 * 1.3 * 0.7 * 1.5   # == 5.733

.sem <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }

# =============================================================================
# A5/A13 -- NUMBER FORMATTING  (CONVENTIONS.md C7, ADJUDICATION A-9/A-10)
# =============================================================================
# Absorbed verbatim from handoff/house_style.R (author decision, 2026-08-07;
# "< 0.001" collapse, 2026-08-09). Replaces SIX independent local copies of
# fmt3()/format_p_statement() across the GOT chain (behaviour, rosner, three
# aggression scripts) -- the exact "seven-copies" drift failure mode this file
# exists to close off (see HANDOFF.md S3/S11, house_style.R header).
#
# DECISION (p convention, 2026-08-24): the GOT-specific "< 0.0001" floor and
# the "marginal" (0.05 < p <= 0.095) verbal tier are RETIRED. p < 0.001 is the
# floor everywhere in this chain (fmt_p_house); "marginal"/"(trend)" language
# is removed at the call sites (see Phase 1 item 2 of the plan).

fmt_F <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  if (grepl(".", s, fixed = TRUE)) {
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  s
}

fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}

fmt_Fstat <- function(x) if (!is.finite(x)) "NA" else if (abs(x) < 0.001) "< 0.001" else fmt_F(x)
fmt_es3   <- function(x) if (!is.finite(x)) "NA" else if (abs(x) < 0.001) "< 0.001" else fmt3(x)

fmt_p_house <- function(p) {
  if (!is.finite(p)) return("NA")
  if (round(p, 3) == 0) return("p < 0.001")
  paste0("p = ", fmt3(p))
}

fmt_p_boot <- function(p, N) {
  floor_p <- 1 / (N + 1)
  if (!is.finite(p)) return("NA")
  if (p <= floor_p) return(paste0("p \u2264 ", fmt3(floor_p)))
  paste0("p = ", fmt3(p))
}

# =============================================================================
# A6 -- EFFECT SIZES  (HANDOFF.md S6, CONVENTIONS.md C6, RESULTS_STYLE.md R5)
# =============================================================================
# effect_size() is the single entry point every caption/results builder in the
# GOT chain must call. stat_type = "chisq" returns es = NA_real_ BY
# CONSTRUCTION so a Wald chi-square caption cannot attach a partial eta-squared
# even by accident.

eta2_p <- function(F, df1, df2) {
  if (!(is.finite(F) && is.finite(df1) && is.finite(df2))) return(NA_real_)
  (F * df1) / (F * df1 + df2)
}

omega2_p <- function(F, df1, df2) {
  if (!(is.finite(F) && is.finite(df1) && is.finite(df2))) return(NA_real_)
  max(0, (df1 * (F - 1)) / (df1 * (F - 1) + df2 + 1))
}

effect_size <- function(stat_type = c("F", "chisq", "t"), F = NA_real_, df1 = NA_real_,
                        df2 = NA_real_, t = NA_real_, df = NA_real_) {
  stat_type <- match.arg(stat_type)
  if (identical(stat_type, "chisq")) {
    return(list(es = NA_real_, es_type = NA_character_,
               note = "Wald chi-square carries no partial eta-squared by design"))
  }
  if (identical(stat_type, "F")) {
    return(list(es = eta2_p(F, df1, df2), es_type = "eta2_p",
               es_secondary = omega2_p(F, df1, df2), es_secondary_type = "omega2_p"))
  }
  # stat_type == "t": signed Hedges' g, noncentral-t CI, J correction.
  if (!.HAS_EFFECTSIZE || !is.finite(t) || !is.finite(df) || df <= 1) {
    return(list(es = NA_real_, es_type = "g", note = "effectsize unavailable or df<=1"))
  }
  J     <- 1 - 3 / (4 * df - 1)
  d_res <- tryCatch(effectsize::t_to_d(t, df, ci = 0.95), error = function(e) NULL)
  if (is.null(d_res) || !nrow(d_res)) return(list(es = NA_real_, es_type = "g"))
  list(es = d_res$d[1] * J, ci_lo = d_res$CI_low[1] * J, ci_hi = d_res$CI_high[1] * J,
      es_type = "g (total-variance standardiser -- report alongside g_resid, never pooled)")
}

# hedges_g(): convenience wrapper for a 1-df pairwise contrast from an emmeans
# / multcomp contrast row (t statistic + df). Returns the signed g and its CI.
hedges_g <- function(t, df) effect_size("t", t = t, df = df)

# pct_change(): biological effect size mandatory alongside every statistical
# one (RESULTS_STYLE.md R5). ref is the reference-group mean.
pct_change <- function(comparison_mean, ref_mean) {
  if (!is.finite(comparison_mean) || !is.finite(ref_mean) || ref_mean == 0) return(NA_real_)
  100 * (comparison_mean - ref_mean) / ref_mean
}

# =============================================================================
# A7 -- THE ANOVA STATEMENT / CAPTION STRING  (CONVENTIONS.md C6)
# =============================================================================
# Replaces format_anova_statement() (behaviour_models.R) and the ad
# hoc caption assembly in figure_3_locomotion.R / the mono-plasma scripts. Keeps
# GOT's existing Unicode-subscript rendering (to_subscript_statement) rather
# than the handoff's <sub>/<sup> HTML, because GOT's caption device path is
# cairo_pdf + Segoe UI Symbol, not ggtext::element_markdown(). Order is law:
# statistic -> p -> effect size (C6 [INV]).

to_subscript_statement <- function(x) {
  s <- as.character(x)
  for (i in 0:9)
    s <- gsub(as.character(i), intToUtf8(0x2080L + i), s, fixed = TRUE)
  s
}

# anova_statement(): term_label, stat_type ("F"|"chisq"), the statistic, dfs,
# and p. For "F" appends eta2_p (omega2_p is NOT printed in-caption -- C6
# [DR]: eta2_p in the caption, omega2_p in the supplementary anova.csv only).
# For "chisq" the effect-size clause is structurally empty (effect_size()
# returns es = NA), matching handoff/house_style.R's anova_caption().
anova_statement <- function(term_label, stat_type = c("F", "chisq"), stat, df1, df2 = NA, p) {
  stat_type <- match.arg(stat_type)
  p_txt <- fmt_p_house(p)
  if (identical(stat_type, "F")) {
    es <- effect_size("F", F = stat, df1 = df1, df2 = df2)
    stat_txt <- sprintf("F%s,%s = %s", to_subscript_statement(round(df1)),
                        to_subscript_statement(round(df2)), fmt_Fstat(stat))
    es_txt   <- if (is.finite(es$es))
      sprintf(", \u03b7\u00b2p = %s", fmt_es3(es$es)) else ""
  } else {
    es       <- effect_size("chisq")
    stat_txt <- sprintf("\u03c7\u00b2%s = %s", to_subscript_statement(round(df1)), fmt_Fstat(stat))
    es_txt   <- ""   # structurally empty -- see effect_size()
  }
  sprintf("%s %s, %s%s", term_label, stat_txt, p_txt, es_txt)
}

# =============================================================================
# A8 -- CLD PIPELINE  (CONVENTIONS.md C5)
# =============================================================================
# Four steps in MANDATORY order: multcomp::cld() -> drop redundant letters ->
# remap (control-first) -> post-condition guard. Letters are NEVER truncated.
# .GOT_CELL_ORDER anchors "control-first" on GOT's own cell order
# (TG_BREAKS_PLOT, Control before Exercise, No-choice before Exercise-choice)
# rather than CHOICE's two-level c("control","exercise choice").

.GOT_CELL_ORDER <- TG_BREAKS_PLOT

.cld_display_valid <- function(groups, same) {
  ids <- names(groups)
  lets <- lapply(groups, function(g) strsplit(trimws(g), "")[[1]])
  for (i in seq_along(ids)) for (j in seq_along(ids)) {
    if (j <= i) next
    shares <- length(intersect(lets[[i]], lets[[j]])) > 0L
    if (shares != isTRUE(same[i, j])) return(FALSE)
  }
  TRUE
}

.cld_drop_redundant_letters <- function(groups, pmat, alpha = 0.05, label = "") {
  ids <- names(groups)
  if (length(ids) < 2L) return(groups)
  same <- pmat >= alpha
  diag(same) <- TRUE
  if (!.cld_display_valid(groups, same)) {
    warning("CLD reduction skipped", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": letters inconsistent with the p-value matrix. Left unchanged.")
    return(groups)
  }
  repeat {
    dropped <- FALSE
    for (i in seq_along(ids)) {
      lets <- sort(unique(strsplit(trimws(groups[[i]]), "")[[1]]))
      if (length(lets) < 2L) next
      for (l in lets) {
        trial <- groups
        keep  <- setdiff(strsplit(trimws(trial[[i]]), "")[[1]], l)
        if (!length(keep)) next
        trial[[i]] <- paste(sort(unique(keep)), collapse = "")
        if (.cld_display_valid(trial, same)) { groups <- trial; dropped <- TRUE; break }
      }
    }
    if (!dropped) break
  }
  groups
}

.cld_remap_control_first <- function(cld_df, cell_col = "cell", group_col = ".group",
                                     cell_order = .GOT_CELL_ORDER) {
  if (is.null(cld_df) || !nrow(cld_df)) return(cld_df)
  ord <- order(match(as.character(cld_df[[cell_col]]), cell_order))
  seen <- character(0)
  for (i in ord) {
    lets <- strsplit(trimws(as.character(cld_df[[group_col]][i])), "")[[1]]
    for (l in lets) if (nzchar(l) && !(l %in% seen)) seen <- c(seen, l)
  }
  if (!length(seen) || length(seen) > length(letters)) return(cld_df)
  map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df[[group_col]] <- vapply(as.character(cld_df[[group_col]]), function(g) {
    r <- map[strsplit(trimws(g), "")[[1]]]; r <- r[!is.na(r)]
    if (!length(r)) return(trimws(g))
    paste(sort(unique(r)), collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  cld_df
}

.cld_is_uniform <- function(cld_df, group_col = ".group") {
  if (is.null(cld_df) || !nrow(cld_df) || !group_col %in% names(cld_df)) return(FALSE)
  length(unique(trimws(as.character(cld_df[[group_col]])))) == 1L
}

.cld_check_control_first <- function(cld_df, label = "", cell_col = "cell", group_col = ".group",
                                     cell_order = .GOT_CELL_ORDER) {
  if (is.null(cld_df) || !nrow(cld_df)) return(invisible(TRUE))
  i <- which(as.character(cld_df[[cell_col]]) == cell_order[1])
  if (!length(i)) return(invisible(TRUE))
  g  <- trimws(as.character(cld_df[[group_col]][i[1]]))
  ok <- nzchar(g) && substr(g, 1, 1) == "a"
  if (!ok) warning("CLD policy violated", if (nzchar(label)) paste0(" [", label, "]") else "",
                   ": ", cell_order[1], " carries '", g, "', expected a group starting with 'a'.")
  invisible(ok)
}

# cld_pipeline(): the four mandatory steps in order. `cld_df` needs a `cell`
# column (identifying the Treatment:Group cell) and a `.group` column
# (multcomp::cld() output); `pmat` a symmetric p-value matrix with dimnames
# matching `cell_ids`.
cld_pipeline <- function(cld_df, cell_ids, pmat, alpha = 0.05, label = "",
                         cell_col = "cell", cell_order = .GOT_CELL_ORDER) {
  if (.cld_is_uniform(cld_df)) return(list(cld_df = NULL, omitted = TRUE,
    reason = "every cell shares one group; no pairwise difference -- CLD omitted per house rule"))
  groups <- stats::setNames(as.list(trimws(as.character(cld_df$.group))), cell_ids)
  groups <- .cld_drop_redundant_letters(groups, pmat, alpha = alpha, label = label)
  cld_df$.group <- unlist(groups, use.names = FALSE)
  # cell_order may be a compound key (e.g. Treatment:Group:interval) supplied
  # by the caller for grids with an extra factor beyond Treatment:Group --
  # defaults to the plain 4-cell .GOT_CELL_ORDER.
  cld_df <- .cld_remap_control_first(cld_df, cell_col = cell_col, cell_order = cell_order)
  .cld_check_control_first(cld_df, label = label, cell_col = cell_col, cell_order = cell_order)
  list(cld_df = cld_df, omitted = FALSE)
}

CLD_LETTER_FRAC <- 0.15
CLD_EXPAND_FRAC <- 0.26
cld_place <- function(cell_max, data_range) cell_max + CLD_LETTER_FRAC * data_range

CLD_CAPTION_NOTE <- paste(
  "Compact letter display: cells sharing a letter do not differ at p < 0.05.",
  "Letters are relabelled so that Control/No-choice always begins with 'a'; this",
  "is a relabelling only, and which cells share a letter is unchanged. A cell",
  "carrying two letters (e.g. 'ab') belongs to two non-difference groups -- it",
  "differs from neither -- and cannot be reduced to one letter without asserting",
  "a difference the model does not support. Panels without letters have no",
  "significant pairwise differences at all: every cell falls in the same group,",
  "so letters would carry no information and are omitted rather than printed as",
  "identical marks."
)

# =============================================================================
# A10 -- EXPORT  (CONVENTIONS.md C9)
# =============================================================================
# GOT's two-tree Figures/{PDF,PNG} layout is retained (one output tree per
# figure format, not per run -- C9 [INV] "one output tree per run" is
# satisfied because figure_3_locomotion.R intentionally rebuilds these paths in
# place each run, same as CHOICE's house_save() target). bg = "white" is
# added to the PNG branch -- missing in the original ggsave() call.
house_save <- function(fig, name, out_dir_pdf, out_dir_png, w, h, dpi = 300) {
  if (is.null(fig)) { warning("house_save: NULL plot, skipped ", name); return(invisible()) }
  dir.create(out_dir_pdf, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_png, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir_pdf, paste0(name, ".pdf")), fig,
                  width = w, height = h, units = "in", device = grDevices::cairo_pdf,
                  bg = "white")
  ggplot2::ggsave(file.path(out_dir_png, paste0(name, ".png")), fig,
                  width = w, height = h, units = "in", dpi = dpi, bg = "white")
  message("  saved: ", name)
}

# house_caption_txt(): sibling caption .txt file (CONVENTIONS.md C10). Every
# saved figure gets one alongside it.
house_caption_txt <- function(out_dir, name, figure_number, title, n_per_cell,
                              unit_desc, panel_key, metric_desc, model_term_df,
                              error_bar_desc, robustness_note = NULL,
                              omission_note = NULL, cld_note = CLD_CAPTION_NOTE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lines <- c(
    sprintf("Figure %s: %s", figure_number, title),
    "",
    sprintf("N per cell: %s (unit: %s)", n_per_cell, unit_desc),
    sprintf("Panel key: %s", panel_key),
    "",
    metric_desc,
    "",
    sprintf("Model term / df convention: %s", model_term_df),
    sprintf("Error bars: %s", error_bar_desc),
    "",
    cld_note
  )
  if (!is.null(robustness_note)) lines <- c(lines, "", robustness_note)
  if (!is.null(omission_note))   lines <- c(lines, "", omission_note)
  writeLines(lines, file.path(out_dir, paste0(name, "_caption.txt")))
  invisible(TRUE)
}

# =============================================================================
# A11 -- SINGULARITY GUARD  (HANDOFF.md S5 step 4, CONVENTIONS.md C12)
# =============================================================================
# The bundle's actual rule (not "fall back to any non-singular candidate"):
# if SOME but not ALL candidates are singular, singular ones become
# ineligible; if EVERY candidate is singular, the guard keeps the full set
# rather than discarding everything. Returns the eligible subset (a
# zero-length-safe filter) plus a full aicc_selection.csv-ready table.
# =============================================================================
# A11b -- SINGULARITY REDUCTION  (parsimonious random-effect structure)
# =============================================================================
# When a fitted mixed model is SINGULAR, drop the variance component that is
# causing it and refit -- iteratively, until the fit is non-singular or no
# random effect remains. This applies REGARDLESS of whether the term is
# design-motivated.
#
# WHY (see docs/GOT_methods.md S3.2 for the manuscript-facing version):
#
#   1. A singular fit means the random-effect structure is not identified by
#      the data. The estimated covariance matrix is on the boundary of the
#      parameter space and is not positive definite; lme4's own documentation
#      (?isSingular) describes such fits as overfitted and recommends removing
#      terms. Retaining the term does not make inference conservative -- it
#      makes the likelihood surface degenerate.
#
#   2. Kenward-Roger -- the denominator-df approximation this pipeline relies
#      on -- is DERIVED under the assumption that the variance parameters lie
#      in the INTERIOR of the parameter space. At a boundary estimate the
#      adjustment is not well defined, so retaining a zero-variance term
#      corrupts precisely the machinery it was retained to protect.
#
#   3. Empirically, in this dataset, a term whose variance is estimated at
#      exactly zero still halved the reported denominator df (Group moved from
#      ~8 to ~4). A component that explains no variance should not cost degrees
#      of freedom; that df is an artefact of model specification, not of the
#      design.
#
#   4. Bates, Kliegl, Vasishth & Baayen (2015), "Parsimonious mixed models",
#      and Matuschek, Kliegl, Vasishth, Baayen & Bates (2017), J. Mem. Lang.
#      86:305-315, both show that maximal random-effect structures which the
#      data cannot support cost substantial power while buying little Type I
#      error protection, and recommend reducing to a structure the data
#      identify.
#
# Every drop is RECORDED (`dropped`, `drop_reason`) so the reduction is visible
# in aicc_selection.csv rather than silent, and so the design-motivated
# structure that was attempted first is still on the record.

# Variance components of a fitted model, as a named numeric vector of SDs,
# excluding the residual. Handles lme4/lmerTest and glmmTMB.
re_sd_components <- function(m) {
  vc <- tryCatch(
    if (inherits(m, "glmmTMB")) lme4::VarCorr(m)$cond else lme4::VarCorr(m),
    error = function(e) NULL)
  if (is.null(vc)) return(numeric(0))
  out <- vapply(vc, function(g) sqrt(max(diag(as.matrix(g)))), numeric(1))
  out <- out[names(out) != "Residual"]
  out[is.finite(out)]
}

.is_singular_any <- function(m) {
  if (inherits(m, "glmmTMB")) {
    sds <- re_sd_components(m)
    return(length(sds) > 0 && any(sds < 1e-5))
  }
  tryCatch(isTRUE(lme4::isSingular(m)), error = function(e) FALSE)
}

# reduce_singular_re(): fit_fn(fixed_str, re_str_or_NULL) -> model or NULL.
# re_terms is a character vector of INDIVIDUAL bar terms, e.g.
#   c("(1 | trial_date)", "(1 | trial_id)", "(1 | fish_id_matched_trial)")
# Returns the first non-singular fit found, with the drop trail.
reduce_singular_re <- function(fit_fn, re_terms, label = "") {
  dropped <- character(0)
  reasons <- character(0)
  current <- re_terms
  repeat {
    re_str <- if (length(current)) paste(current, collapse = " + ") else NULL
    m <- tryCatch(fit_fn(re_str), error = function(e) NULL)
    if (is.null(m)) {
      # This structure will not fit at all; drop the last term and try again.
      if (!length(current))
        return(list(model = NULL, re = NA_character_, dropped = dropped,
                    reasons = reasons, singular = NA, n_dropped = length(dropped)))
      dropped <- c(dropped, current[length(current)])
      reasons <- c(reasons, "did not converge")
      current <- current[-length(current)]
      next
    }
    if (!.is_singular_any(m) || !length(current))
      return(list(model = m, re = if (length(current)) re_str else "none",
                  dropped = dropped, reasons = reasons,
                  singular = .is_singular_any(m), n_dropped = length(dropped)))

    # Singular: identify the offending component as the one with the smallest
    # estimated SD, and map it back to the bar term naming that grouping factor.
    sds <- re_sd_components(m)
    if (!length(sds)) {
      dropped <- c(dropped, current[length(current)])
      reasons <- c(reasons, "singular, no variance components recoverable")
      current <- current[-length(current)]
      next
    }
    worst_grp <- names(sds)[which.min(sds)]
    hit <- grep(worst_grp, current, fixed = TRUE)
    if (!length(hit)) hit <- length(current)   # fall back to the last term
    dropped <- c(dropped, current[hit[1]])
    reasons <- c(reasons, sprintf("singular; SD(%s) = %.3g", worst_grp, min(sds)))
    if (nzchar(label))
      cat(sprintf("    [%s] singular -> dropping %s (SD = %.3g)\n",
                  label, current[hit[1]], min(sds)))
    current <- current[-hit[1]]
  }
}

# Split a random-effect string such as
#   "(1 | trial_date) + (1 | trial_id) + (1 | fish_id_matched_trial)"
# into its individual bar terms, so reduce_singular_re() can drop them one at
# a time. Splitting on ") + (" rather than " + " keeps any term that itself
# contains a plus sign intact.
split_re_terms <- function(re_str) {
  if (is.null(re_str) || is.na(re_str) || !nzchar(re_str)) return(character(0))
  parts <- strsplit(re_str, "\\)\\s*\\+\\s*\\(", perl = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- ifelse(grepl("^\\(", parts), parts, paste0("(", parts))
  parts <- ifelse(grepl("\\)$", parts), parts, paste0(parts, ")"))
  parts
}

# Accumulator so every reduction is visible in one artifact rather than only
# in the console log.
.RE_REDUCTION_LOG <- list()
re_reduction_log_add <- function(label, model_type, attempted, retained, dropped, reasons) {
  .RE_REDUCTION_LOG[[length(.RE_REDUCTION_LOG) + 1]] <<- data.frame(
    label = label, model_type = model_type,
    re_attempted = attempted,
    re_retained  = if (length(retained) && nzchar(retained)) retained else "none",
    n_dropped    = length(dropped),
    re_dropped   = if (length(dropped)) paste(dropped, collapse = " + ") else NA_character_,
    drop_reason  = if (length(reasons)) paste(reasons, collapse = "; ") else NA_character_,
    stringsAsFactors = FALSE)
}

# Writes re_reduction.csv into `run_dir`. Written whether or not anything was
# dropped, so an empty file means "nothing was singular", not "nobody looked".
re_reduction_log_write <- function(run_dir) {
  red_tbl <- if (length(.RE_REDUCTION_LOG)) dplyr::bind_rows(.RE_REDUCTION_LOG) else
    data.frame(label = character(0), model_type = character(0),
               re_attempted = character(0), re_retained = character(0),
               n_dropped = integer(0), re_dropped = character(0),
               drop_reason = character(0))
  out <- file.path(run_dir, "re_reduction.csv")
  write.csv(red_tbl, out, row.names = FALSE)
  cat("  wrote re_reduction.csv (", nrow(red_tbl), " reductions)\n", sep = "")
  invisible(out)
}
apply_singularity_guard <- function(cand_results, is_singular_fn) {
  singular <- vapply(cand_results, function(cr) {
    tryCatch(isTRUE(is_singular_fn(cr$model)), error = function(e) NA)
  }, logical(1))
  all_singular <- isTRUE(all(singular, na.rm = TRUE)) && !any(is.na(singular))
  eligible <- if (all_singular) rep(TRUE, length(cand_results)) else !singular
  eligible[is.na(eligible)] <- TRUE   # unknown singularity status -> don't exclude
  list(eligible = eligible, singular = singular, all_singular = all_singular)
}

# =============================================================================
# A12 -- PARAMETRIC-BOOTSTRAP LRT  (HANDOFF.md S7)
# =============================================================================
# For treatment terms where Kenward-Roger/Satterthwaite is undefined (beta
# GLMM). pbkrtest::PBmodcomp has no glmmTMB method, hence hand-rolled. N and
# seed are what Methods and the caption print -- a bootstrap without a stated
# N and seed is not a check (HANDOFF S7 [INV]).
pb_lrt <- function(fixed_str, reduced_fixed_str, re_str, data, family,
                   nsim = 1000L, seed) {
  na_out <- list(lrt_obs = NA_real_, p_pb = NA_real_, nsim_ok = 0L, N = nsim, seed = seed, method = "na")
  full_form    <- tryCatch(stats::as.formula(paste(".y ~", fixed_str, "+", re_str)), error = function(e) NULL)
  reduced_form <- tryCatch(stats::as.formula(paste(".y ~", reduced_fixed_str, "+", re_str)), error = function(e) NULL)
  if (is.null(full_form) || is.null(reduced_form)) return(modifyList(na_out, list(method = "bad_formula")))

  .fit <- function(fmla, dat) tryCatch(
    suppressWarnings(suppressMessages(glmmTMB::glmmTMB(
      fmla, data = dat, family = family,
      control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500))))),
    error = function(e) NULL)

  m_full <- .fit(full_form, data); m_red <- .fit(reduced_form, data)
  if (is.null(m_full) || is.null(m_red)) return(modifyList(na_out, list(method = "observed_fit_failed")))
  lrt_obs <- tryCatch(max(0, 2 * (as.numeric(stats::logLik(m_full)) - as.numeric(stats::logLik(m_red)))),
                      error = function(e) NA_real_)
  if (!is.finite(lrt_obs)) return(modifyList(na_out, list(method = "loglik_failed")))

  set.seed(seed)
  sim_y <- tryCatch(stats::simulate(m_red, nsim = nsim), error = function(e) NULL)
  if (is.null(sim_y)) return(modifyList(na_out, list(lrt_obs = lrt_obs, method = "simulate_failed")))

  sim_lrt <- vapply(seq_len(nsim), function(i) {
    d_sim <- data; d_sim$.y <- sim_y[[i]]
    mf <- .fit(full_form, d_sim); mr <- .fit(reduced_form, d_sim)
    if (is.null(mf) || is.null(mr)) return(NA_real_)
    v <- tryCatch(2 * (as.numeric(stats::logLik(mf)) - as.numeric(stats::logLik(mr))), error = function(e) NA_real_)
    if (!is.finite(v)) NA_real_ else max(0, v)
  }, numeric(1))

  n_ok <- sum(is.finite(sim_lrt))
  if (n_ok < nsim * 0.5)
    return(modifyList(na_out, list(lrt_obs = lrt_obs, nsim_ok = n_ok, method = "too_many_sim_failures")))
  p_pb <- (1 + sum(sim_lrt[is.finite(sim_lrt)] >= lrt_obs)) / (n_ok + 1)
  list(lrt_obs = lrt_obs, p_pb = p_pb, nsim_ok = n_ok, N = nsim, seed = seed, method = "pb_lrt")
}

# =============================================================================
# A12b -- ONE-DEGREE-OF-FREEDOM (MARGINAL) PARAMETRIC-BOOTSTRAP LRT
# =============================================================================
# pb_lrt() compares NESTED models, so dropping a main effect from a fully-crossed
# factor structure also drops its interactions: those rows are composite tests
# (4 df for a main factor) and do not answer the question the Wald / Kenward-Roger
# tests answer. The marginal question -- "is the Group MAIN effect zero, given
# the interactions are in the model?" -- is a single coefficient under
# sum-to-zero coding. It is tested by coding each two-level factor as +1 / -1
# numeric columns and removing ONLY the main-effect column from the reduced
# model. With two-level factors under contr.sum the full numeric model is the
# same model as the reported factor model (same log-likelihood), and the reduced
# model has exactly one fewer fixed parameter, so this is a genuine 1-df Type
# III-style test that reuses pb_lrt() unchanged (it takes formula strings).
#
# Add `<factor>_c` = +1 (first level) / -1 (last level) for each two-level factor.
# Refuses anything else: ±1 coding of a 3-level factor would be silently wrong.
add_sum_coded_cols <- function(data, factors = c("Treatment", "Group", "interval")) {
  for (f in factors) {
    if (!f %in% names(data)) stop("add_sum_coded_cols: column '", f, "' not found")
    fac <- droplevels(factor(data[[f]]))
    if (nlevels(fac) != 2L)
      stop("add_sum_coded_cols: '", f, "' has ", nlevels(fac),
           " levels; +/-1 coding needs exactly 2")
    data[[paste0(f, "_c")]] <- ifelse(fac == levels(fac)[1], 1, -1)
  }
  data
}

# Fixed-effect strings (over the `_c` columns) for the marginal test of `term`.
pb_1df_specs <- function(term = c("Group", "Treatment")) {
  term <- match.arg(term)
  oth <- paste0(setdiff(c("Treatment", "Group", "interval"), term), "_c")
  x <- paste0(term, "_c")
  list(full    = "Treatment_c * Group_c * interval_c",
       reduced = paste0(oth[1], " * ", oth[2], " + ", x, ":", oth[1], " + ",
                        x, ":", oth[2], " + ", x, ":", oth[1], ":", oth[2]),
       dropped = paste0(term, " main effect only (its interactions are retained)"),
       df = 1L)
}
# =============================================================================
# A14 -- PREFLIGHT CHECK  (CONVENTIONS.md App. C)
# =============================================================================
preflight_check <- function(env = globalenv()) {
  failures <- character(0)
  chk <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

  chk(identical(TG_COLS, c(
    "Exercise:Exercise choice" = "#D55E00", "Exercise:No choice" = "#56B4E9",
    "Control:Exercise choice" = "#E69F00", "Control:No choice" = "#0072B2"
  )), "1: TG_COLS does not match the declared GOT palette")
  chk(identical(INK_STAT, "black"), "1b: INK_STAT is not black")
  chk(is.null(BASE_THEME$theme$panel.grid) ||
      inherits(BASE_THEME$theme$panel.grid, "element_blank") ||
      !("panel.grid" %in% names(BASE_THEME$theme)),
      "3: BASE_THEME sets a non-blank panel.grid")
  chk(fmt_F(30) == "30", "16a: fmt_F(30) != '30' -- integer df trap")
  chk(fmt3(0.5) == "0.50", "16b: fmt3(0.5) != '0.50'")
  chk(fmt3(0.77) == "0.77", "16c: fmt3(0.77) != '0.77'")
  chk(fmt_p_house(0.0004) == "p < 0.001", "9: fmt_p_house floor incorrect")
  chk(fmt_p_boot(0.001, 200) == "p \u2264 0.005", "16d: fmt_p_boot floor incorrect")
  chk(is.na(effect_size("chisq")$es), "8/9: effect_size('chisq') must return es = NA")
  chk(length(unique(TG_COLS)) == length(TG_COLS), "1c: TG_COLS has a duplicated colour")

  if (length(failures)) {
    warning("preflight_check() FAILED:\n  - ", paste(failures, collapse = "\n  - "))
    return(invisible(FALSE))
  }
  message("preflight_check(): all checks passed.")
  invisible(TRUE)
}
