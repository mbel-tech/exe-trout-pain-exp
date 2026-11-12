################################################################################
# figure_3_locomotion.R
# ─────────────────────────────────────────────────────────────────────────────
# R-drawn figures. Figure 3 of the article (change in locomotion between the two
# observation intervals) is the file Figure_5 written by this script (working
# number). The other figures it writes (Figure_4, Figure_6 to Figure_9 and their
# variants) are working versions; the article versions of those figures are
# drawn by the Python scripts in this folder.
# Working numbering: 4 behaviour, 5 interval-individual, 6 interval-school,
#            7 aggression, 8 brain neurochemistry, 9 plasma.
# Inputs: the latest complete runs of behaviour_models.R, plasma_monoamine_models.R
# and aggression_models.R under output/.
#
# Output: PDF + PNG in output/Figures/{PDF,PNG}
################################################################################

suppressMessages({
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(patchwork)
})

# -- House style / sensitivity modules (single source of truth for palette,
#    theme, geometry constants, CLD pipeline, effect sizes, formatters).
#    Everything this file used to define locally (TG_COLS, TG_LBLS, GRP_SHPS,
#    TG_SHPS, TRT_SHPS, GRP_LVLS, TRT_LVLS, TG_BREAKS_PLOT, SCALE_COL,
#    SCALE_SHP_GRP, FONT_CAP, BASE_THEME, CLD_SIZE, and the FIG45_*/geometry
#    block) is now sourced from here instead -- see the plan's Phase 0/1.
# GOT_ROOT can be overridden from the environment (the repo has lived on both D: and E:).
GOT_ROOT <- Sys.getenv("GOT_ROOT", ".")
source(file.path(GOT_ROOT, "05_statistics/shared/house_style.R"))

# ── Output directories ─────────────────────────────────────────────────────────
# GOT_FIG_DIR redirects the output (e.g. to a scratch folder) without touching the shared Figures folder.
BASE  <- Sys.getenv("GOT_FIG_DIR", file.path(GOT_ROOT, "output/Figures"))

# ── ANOVA statements under the panels: two figure sets ────────────────────────
# STMT_MODE = "significant" (manuscript set, written to BASE/{PDF,PNG}): a panel
#   shows only the statements with p < STMT_ALPHA; a panel with none shows none.
# STMT_MODE = "all" (reference set, written to BASE/all_statements/{PDF,PNG}):
#   every term of the model is shown, significant or not.
# A plain run builds BOTH sets (the script re-sources itself once at the end).
# Set GOT_STMT_MODE=significant or =all in the environment to build only one.
# Every statement carries an effect size in both sets: partial eta-squared for
# F-tests (from the pipeline) and Cohen's w = sqrt(chi-square / N) for Wald
# chi-square tests, N = number of observations the model was fitted to.
if (!exists("STMT_MODE")) STMT_MODE <- Sys.getenv("GOT_STMT_MODE", "significant")
stopifnot(STMT_MODE %in% c("significant", "all"))
STMT_ALPHA <- 0.05
if (!exists("THIS_FILE")) THIS_FILE <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE)))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(of)) NA_character_ else normalizePath(of)
})
if (identical(STMT_MODE, "all")) BASE <- file.path(BASE, "all_statements")
cat("ANOVA statements shown:", STMT_MODE, "->", BASE, "\n")
D_PDF <- file.path(BASE, "PDF")
D_PNG <- file.path(BASE, "PNG")
for (d in c(D_PDF, D_PNG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Run discovery ─────────────────────────────────────────────────────────────
# HANDOFF.md S10 [INV]: downstream consumers AUTO-DISCOVER the newest matching
# run directory, require a completed-run sentinel, and reject backup siblings.
# A-12 records what happens when this discipline lapses: the behaviour and
# monoamine pins here were hardcoded and went stale twice, most recently
# leaving the manuscript figures one analysis vintage behind the manuscript
# text. Discovery removes the class of error rather than moving it to a newer
# path. Set GOT_PIN_BEHAV_RUN / GOT_PIN_MONO_RUN in the environment to override
# for a deliberate rebuild of a withdrawn vintage (S10 [DR]).
EASY_OUT_BASE <- file.path(GOT_ROOT, "output")

discover_run_dir <- function(base, pattern, what, override_env = NULL) {
  if (!is.null(override_env) && nzchar(Sys.getenv(override_env, ""))) {
    p <- Sys.getenv(override_env)
    cat("  [", what, "] PINNED via ", override_env, ": ", basename(p), "\n", sep = "")
    return(p)
  }
  runs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  runs <- runs[grepl(pattern, basename(runs))]
  runs <- runs[!grepl("BACKUP|BAK|PRE_|_old", basename(runs), ignore.case = TRUE)]
  runs <- runs[file.exists(file.path(runs, "_RUN_COMPLETE"))]
  if (!length(runs))
    stop("No completed run matching '", pattern, "' under ", base,
         " (a run directory must carry a _RUN_COMPLETE sentinel).", call. = FALSE)
  chosen <- runs[order(basename(runs))][length(runs)]
  cat("  [", what, "] ", basename(chosen), "\n", sep = "")
  chosen
}

cat("Discovering analysis runs:\n")

# ── Raw monoamine datasets for Figure 4 ───────────────────────────────────────
# These files are used to draw true one-point-per-sample scatterplots.
MONO_RUN_DIR <- discover_run_dir(EASY_OUT_BASE, "^rosner_[0-9]{8}_[0-9]{6}$",
                                 "monoamine/plasma", "GOT_PIN_MONO_RUN")
RAW_MONO_DIR <- file.path(MONO_RUN_DIR, "datasets")
MONO_DM_CSV  <- file.path(RAW_MONO_DIR, "mono_DM.csv")
MONO_RN_CSV  <- file.path(RAW_MONO_DIR, "mono_RN.csv")
PLASMA_CSV   <- file.path(RAW_MONO_DIR, "plasma.csv")

# Pipeline ANOVA statements (one CSV per effect; columns: key, unit, indicator,
# statement). Used to source the Figure 4 (plasma) captions.
ANOVA_STMT_DIR <- file.path(MONO_RUN_DIR, "anovas")
# Tukey CLD letters (one CSV per indicator; columns include Treatment, Group, cld).
CLD_MONO_DIR   <- file.path(MONO_RUN_DIR, "cld")

# ── Raw behaviour datasets for Figure 4 ───────────────────────────────────────
# These files are used to draw raw scatterpoints for behaviour panels.
# Figure 4A is a group-dynamics metric, so the raw observational unit is trial_id.
# Figure 4B–C are individual metrics, so rows are aggregated to one point per
# fish/trial using master_fish_id across intervals 2 and 3.
BEHAV_RUN_DIR <- discover_run_dir(EASY_OUT_BASE, "^[0-9]{8}_[0-9]{6}$",
                                  "behaviour", "GOT_PIN_BEHAV_RUN")
RAW_BEHAV_DIR <- file.path(BEHAV_RUN_DIR, "datasets")
GD_INTERVAL_CSV  <- file.path(RAW_BEHAV_DIR, "gd_interval.csv")
GD_TRIAL_CSV     <- file.path(RAW_BEHAV_DIR, "gd_trial.csv")
IND_INTERVAL_CSV <- file.path(RAW_BEHAV_DIR, "ind_interval.csv")
IND_TRIAL_CSV    <- file.path(RAW_BEHAV_DIR, "ind_trial.csv")
# Behaviour CLD letters (trial-level). All behaviour Treatment×Group interactions
# are non-significant (all cells share letter "a"), so no CLD letters will appear
# in Figure 6. read_cld() returns NULL automatically when all cells share a letter.
CLD_BEHAV_DIR      <- file.path(BEHAV_RUN_DIR, "cld")
MARGINAL_CLD_DIR   <- file.path(BEHAV_RUN_DIR, "marginal-emm", "cld")
BEHAV_ANOVA_DIR    <- file.path(BEHAV_RUN_DIR, "anovas")

# ── Aggression analysis output for Figure 7 ───────────────────────────────────
# Most recent run of aggression_models.R. Output folders are named
# YYYYMMDD_HHMMSS, so the lexicographically last subdirectory is the latest run.
# Figure 7 reads its CLD (emmeans, IRR/probability scale), post-hoc contrasts, and
# the raw crowding event table directly from there.
AGG_OUT_BASE <- file.path(EASY_OUT_BASE, "aggression_output")
.agg_runs <- list.dirs(AGG_OUT_BASE, recursive = FALSE, full.names = TRUE)
.agg_runs <- .agg_runs[grepl("/[0-9]{8}_[0-9]{6}$", .agg_runs)]
# HANDOFF.md S10 [INV]: refuse an in-progress or backup run -- require the
# completed-run sentinel, reject BACKUP|BAK|PRE_|_old paths.
.agg_runs <- .agg_runs[!grepl("BACKUP|BAK|PRE_|_old", .agg_runs, ignore.case = TRUE)]
.agg_runs <- .agg_runs[file.exists(file.path(.agg_runs, "_RUN_COMPLETE"))]
AGG_RUN_DIR     <- if (length(.agg_runs)) sort(.agg_runs)[length(.agg_runs)] else AGG_OUT_BASE
AGG_CLD_DIR     <- file.path(AGG_RUN_DIR, "cld")
AGG_POSTHOC_DIR <- file.path(AGG_RUN_DIR, "post-hoc")
AGG_DATA_DIR    <- file.path(AGG_RUN_DIR, "datasets")
AGG_ANOVA_DIR   <- file.path(AGG_RUN_DIR, "anovas")

# ── Okabe-Ito colour + shape encoding ─────────────────────────────────────────
# TG_COLS / TG_LBLS / GRP_SHPS / TG_SHPS / TRT_SHPS sourced from
# house_style.R (A1).

# ── Custom combined legend (Group colours + Treatment shapes) ─────────────────
# Native ggplot guides cannot place two colour swatches next to a single shared
# label, so the requested "Group" + "Treatment" legend block is drawn as a
# standalone plot and composed to the right of the panels with patchwork.
#   Group     : two filled circles per row, sharing one Group label.
#               No choice       -> #0072B2 (Control)  #56B4E9 (Exercise)
#               Exercise choice -> #E69F00 (Control)  #D55E00 (Exercise)
#   Treatment : the two shapes -> 17 Control (triangle), 19 Exercise (circle)
# scale    : overall size multiplier for glyphs + text (Fig 7 uses 3x).
# item_gap : vertical distance between consecutive legend items, in a FIXED
#            data frame (ylim 0.5-6.4). Default 0.5 = the requested half-spacing
#            (distance between items reduced 2x from the original 1.0).
# frame_ylim : data range mapped onto the legend panel's physical height. Glyph
#   sizes are absolute, so to make a legend look identical across panels of
#   different physical heights (e.g. Fig 6 single row vs Fig 7 two-row column),
#   widen frame_ylim proportionally to the taller panel.
build_combined_legend <- function(scale = 1, item_gap = 0.5,
                                  frame_ylim = c(0.5, 6.4)) {
  ts  <- 3.0 * scale; tts <- 3.4 * scale
  ds  <- 4.2 * scale; ss  <- 4.0 * scale
  dx        <- 0.55 * scale          # horizontal gap between the two colour dots
  lab_x     <- 1 + dx + 0.5 * scale  # x where the shared Group label starts
  shp_lab_x <- 1 + 0.5 * scale       # x where the Treatment label starts
  ctr <- mean(frame_ylim)            # items centred in the frame
  ys  <- ctr + item_gap * c(2.5, 1.5, 0.5, -0.5, -1.5, -2.5)
  #     ys: 1=Group title 2=No choice 3=Exercise choice 4=Treatment title 5=Control 6=Exercise
  dots <- data.frame(
    x   = c(1, 1 + dx, 1, 1 + dx),
    y   = c(ys[2], ys[2], ys[3], ys[3]),
    col = c("#0072B2", "#56B4E9", "#E69F00", "#D55E00"))
  dot_labels   <- data.frame(x = lab_x, y = c(ys[2], ys[3]),
                             lab = c("No Choice", "Exercise Choice"))
  shapes       <- data.frame(x = 1, y = c(ys[5], ys[6]), shp = c(17, 19))
  shape_labels <- data.frame(x = shp_lab_x, y = c(ys[5], ys[6]),
                             lab = c("Control", "Exercise"))
  titles       <- data.frame(x = 0.7, y = c(ys[1], ys[4]),
                             lab = c("Group", "Treatment"))
  ggplot() +
    geom_point(data = dots, aes(x, y), colour = dots$col, size = ds) +
    geom_text(data = dot_labels, aes(x, y, label = lab), hjust = 0, size = ts) +
    geom_point(data = shapes, aes(x, y), shape = shapes$shp,
               colour = "black", size = ss) +
    geom_text(data = shape_labels, aes(x, y, label = lab), hjust = 0, size = ts) +
    geom_text(data = titles, aes(x, y, label = lab), hjust = 0,
              fontface = "bold", size = tts) +
    coord_cartesian(xlim = c(0.4, lab_x + 5.5), ylim = frame_ylim, clip = "off") +
    theme_void()
}

# Manual panel-letter tag (top-left corner). Used instead of
# plot_annotation(tag_levels=) so the custom legend plot stays untagged.
tag_theme <- function(size = 10, y = 0.985) {
  theme(plot.tag = element_text(size = size, face = "bold"),
        plot.tag.position = c(0.015, y))
}

# x-axis ordering: No choice LEFT, Exercise choice RIGHT; within: Control inner, Exercise outer
# GRP_LVLS / TRT_LVLS / TG_BREAKS_PLOT / SCALE_COL / SCALE_SHP_GRP sourced from
# house_style.R (A1).
# Dodge order for Panel B (group × treatment × interval): Control left of Exercise,
# No choice left of Exercise choice within each treatment.
AGG_INT_TG_ORDER <- c(
  "Control:No choice",
  "Exercise:No choice",
  "Control:Exercise choice",
  "Exercise:Exercise choice"
)

# FONT_CAP / BASE_THEME sourced from house_style.R (A2) -- statistical
# ink is black there (was grey25 here; C1 [INV]).

# ── Helper: save a patchwork figure as both PDF and PNG ───────────────────────
# PDF uses cairo_pdf so Unicode glyphs in the ANOVA captions render correctly;
# PNG uses the default ragg device, which already handles them.
# ── Keep the statements under a panel inside that panel's column ──────────────
# The statements are plot captions drawn at a fixed font size, so a long line
# (e.g. a three-way interaction with its effect size) used to run into the
# neighbouring panel. fit_captions() measures the longest caption line of every
# panel against the width of the column the panel sits in and, if any is too
# long, shrinks ALL captions of the figure by one common factor (so the panels
# of a figure keep the same caption size). Nothing changes when everything fits.
.caption_width_in <- function(p) {
  cap <- p$labels$caption; el <- p$theme$plot.caption
  if (is.null(cap) || !nzchar(cap) || is.null(el$size)) return(NA_real_)
  lines <- strsplit(cap, "\n", fixed = TRUE)[[1]]
  fam <- if (is.null(el$family) || !nzchar(el$family)) "" else el$family
  ital <- !is.null(el$face) && el$face %in% c("italic", "bold.italic")
  w <- tryCatch(systemfonts::string_width(lines, family = fam, italic = ital,
                                          size = el$size, res = 72) / 72,
                error = function(e) nchar(lines) * el$size * 0.5 / 72)
  max(w)
}
.patch_cell_widths <- function(fig, w) {
  n <- length(fig$patches$plots) + 1
  lay <- fig$patches$layout
  if (!is.null(lay$design) && !is.null(lay$design$l)) {
    l <- as.integer(lay$design$l); r <- as.integer(lay$design$r); ncol <- max(r)
  } else {
    ncol <- if (!is.null(lay$ncol)) lay$ncol else ggplot2::wrap_dims(n, nrow = lay$nrow)[2]
    l <- r <- ((seq_len(n) - 1) %% ncol) + 1
  }
  wd <- lay$widths
  wd <- if (is.null(wd) || !is.numeric(wd)) rep(1, ncol) else rep_len(as.numeric(wd), ncol)
  vapply(seq_len(min(n, length(l))), function(i) w * sum(wd[l[i]:r[i]]) / sum(wd), numeric(1))
}
# Width (in) taken by tag / y title / tick labels left of the panel(s) of p.
.panel_left_in <- function(p, default = 0.8) {
  tryCatch({
    g <- ggplot2::ggplotGrob(p)
    l <- min(g$layout$l[grepl("^panel", g$layout$name)])
    if (l <= 1) return(0)
    sum(grid::convertWidth(g$widths[seq_len(l - 1)], "in", valueOnly = TRUE))
  }, error = function(e) default)
}
fit_captions <- function(fig, w, min_scale = 0.35) {
  if (!inherits(fig, "patchwork")) return(fig)
  plots <- c(fig$patches$plots, list(fig))
  cw <- .patch_cell_widths(fig, w)
  ratio <- 1
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)   # for unit conversion
  for (i in seq_along(cw)) {
    tw <- .caption_width_in(plots[[i]])
    if (!is.finite(tw) || tw <= 0) next
    p_i <- plots[[i]]; class(p_i) <- setdiff(class(p_i), "patchwork")
    ratio <- min(ratio, (cw[i] - .panel_left_in(p_i) - 0.1) / tw)
  }
  if (ratio >= 1) return(fig)
  ratio <- max(ratio * 0.98, min_scale)
  shrink <- function(p) {
    sz <- p$theme$plot.caption$size
    if (is.null(sz)) p else p + theme(plot.caption = element_text(size = sz * ratio))
  }
  for (i in seq_along(fig$patches$plots)) {
    pi <- fig$patches$plots[[i]]
    if (inherits(pi, "ggplot") && !inherits(pi, "patchwork")) fig$patches$plots[[i]] <- shrink(pi)
  }
  sz <- fig$theme$plot.caption$size
  if (!is.null(sz)) fig <- fig + theme(plot.caption = element_text(size = sz * ratio))
  message(sprintf("  captions shrunk to %.0f %% to fit the panel width", 100 * ratio))
  fig
}

save_fig <- function(fig, name, w, h) {
  fig <- fit_captions(fig, w)
  ggsave(file.path(D_PDF, paste0(name, ".pdf")), fig,
         width = w, height = h, units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(D_PNG, paste0(name, ".png")), fig,
         width = w, height = h, units = "in", dpi = 300, bg = "white")
  message("  saved: ", name)
}

# ── Helper: build an ANOVA caption from the pipeline statement CSVs ────────────
# Reads the per-effect statement files, keeps the statements for `key` whose
# p-value is below `p_max`, one per line. A colon is inserted after the effect
# name to match the Figure 4/5 caption style (e.g. "Group: F1,40 = ...").
# fixed = TRUE is used for the colon insertion because the "Treatment × Group"
# label contains a multibyte glyph that corrupts regex backreferences on Windows R.
# Every term is collected here; finish_caption() then adds the effect size and,
# in the manuscript set (STMT_MODE = "significant"), keeps only p < STMT_ALPHA.
# The complete, unselected statements (HANDOFF.md S9) are in the
# all_statements figure set. Order is main effects first, interaction beneath
# (C6 [DR]) -- `spec` already lists Group/Treatment before Treatment x Group.
anova_caption <- function(key) {
  spec <- list(
    c("anova_stmt_Group.csv",           "Group"),
    c("anova_stmt_Treatment.csv",       "Treatment"),
    c("anova_stmt_TreatmentxGroup.csv", "Treatment × Group")
  )
  out <- character(0)
  for (sp in spec) {
    df <- read.csv(file.path(ANOVA_STMT_DIR, sp[1]), stringsAsFactors = FALSE,
                   check.names = FALSE, encoding = "UTF-8")
    s <- df$statement[df$key == key]
    if (length(s) == 0) next
    s <- s[1]
    out <- c(out, sub(paste0(sp[2], " "), paste0(sp[2], ": "), s, fixed = TRUE))
  }
  finish_caption(paste(out, collapse = "\n"), .n_mono(key))
}

# ── Helpers: Figure 4–5 raw scatter + mean ± SE ──────────────────────────────
# Desired final format: raw scatterpoints (one row per sample/fish/trial unit) +
# black mean crossbar + black SE errorbar.

TG_DODGE      <- 0.38
SUMMARY_WIDTH <- 0.22
POINT_SIZE    <- 2.8

# Figure 4–5 requested geometry updates:
# - triple distance between Control/Exercise scatters sharing the same Group tick
# - double scatter spread and mean crossbar/SE errorbar width
# - x tick labels reduced again by 0.5 from the previous setting
# - x tick spacing increased by expanding the discrete x scale
FIG45_DODGE            <- TG_DODGE * 3
FIG45_SUMMARY_WIDTH    <- SUMMARY_WIDTH * 2
FIG45_JITTER_WIDTH     <- SUMMARY_WIDTH * 2
FIG45_X_TEXT_SIZE      <- 8.5 * 3 * 0.7 * 0.5
FIG45_X_EXPAND_ADD     <- 0.65 * 2
# Requested x-tick spacing change: pull the first ("No choice") group hard against
# the y-axis (small left expansion) and let the two group ticks sit farther apart
# (smaller total expansion = ticks span more of the panel width).
FIG45_X_EXPAND_LOW     <- 0.85
FIG45_X_EXPAND_HIGH    <- 0.85
FIG45_TAG_SIZE         <- 15      # panel letters increased x1.5 (10 -> 15)
FIG45_CAPTION_MARGIN_T <- 24
PANEL5B_YLAB_SIZE      <- 9.5 * 2 * 0.7

# Interval figures (Figs 6-7): increase distance between mean datapoints x1.7,
# i.e. wider dodge of the four Treatment:Group series within each interval.
INT_DODGE              <- 0.58 * 1.7

# ── Cross-figure unified settings (this round of requests) ────────────────────
# x-axis tick label size: Fig 6 default (8.5*2.5) reduced by 0.3, then applied to
# every figure's x-axis ticks.
X_TICK_SIZE      <- 8.5 * 2.5 * 0.7
# Padding between the x-axis line and the ANOVA caption (Fig 4 value, applied to all).
CAPTION_MARGIN_T <- 24
# ANOVA caption size shared by Fig 4 & Fig 5 (and never bold).
CAPTION_SIZE_45  <- 7.5 * 2
# Fig 7 ANOVA captions enlarged x2.
FIG7_CAPTION_SIZE <- 7.5 * 1.5 * 2
# Fig 7 panel letter positioning
FIG7_TAG_SIZE <- 10 * 2.5
FIG7_TAG_Y    <- 0.985 + 0.0175
# Panel border: 0.5 cm white margin + white background (Figs 6 & 7)
PANEL_BORDER  <- theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
                       plot.background = element_rect(fill = "white", colour = NA))
# Aggression (Fig 7) panels: taller TOP margin so the high-placed panel letters
# (A/B/C/D) are never clipped when the panels are assembled.
AGG_PANEL_BORDER <- theme(plot.margin = margin(0.95, 0.5, 0.5, 0.5, unit = "cm"),
                          plot.background = element_rect(fill = "white", colour = NA))
# Fig 7 y-axis-title -> axis-line gap (pt). Panel C was margin(r = 0) (ggplot-tight,
# effective base ~4 pt); "0.5x more space" -> base x1.5. Panels A, B, C share this;
# panel D keeps its own prior 3.25 pt gap, also +0.5x.
AGG_YLAB_GAP   <- 6
AGG_YLAB_GAP_D <- 3.25 * 1.5

raw_trial_plot <- function(data, value_col, ylabel, caption_text, dodge = TG_DODGE,
                           ylim = NULL, y_breaks = waiver(),
                           ylab_size = 9.5, caption_size = 7.5,
                           caption_face = "plain", cld = NULL,
                           jitter_width = SUMMARY_WIDTH, summary_width = SUMMARY_WIDTH,
                           point_size = POINT_SIZE, x_text_size = 8.5 * 3,
                           x_expand_add = 0.65, x_expand_low = NULL, x_expand_high = NULL,
                           caption_margin_t = 12,
                           show_legend = TRUE, legend_title = "Treatment : Group",
                           shape_by = c("Group", "TG")) {
  shape_by <- match.arg(shape_by)
  if (is.null(x_expand_low))  x_expand_low  <- x_expand_add
  if (is.null(x_expand_high)) x_expand_high <- x_expand_add

  d <- data %>%
    mutate(
      Group_f = factor(Group, levels = GRP_LVLS),
      Treat_f = factor(Treatment, levels = TRT_LVLS),
      TG      = paste(Treatment, Group, sep = ":"),
      # Shape now encodes Treatment (17 Control, 19 Exercise).
      ShapeKey = as.character(Treat_f)
    )

  sdat <- d %>%
    group_by(Treatment, Group, Group_f, Treat_f, TG, ShapeKey) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      SE   = sd(.data[[value_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[value_col]]))),
      n    = sum(!is.na(.data[[value_col]])),
      .groups = "drop"
    )

  colour_scale <- scale_colour_manual(
    name   = legend_title,
    values = TG_COLS,
    breaks = TG_BREAKS_PLOT,
    labels = TG_LBLS[TG_BREAKS_PLOT]
  )

  shape_scale <- scale_shape_manual(
    name   = "Treatment",
    values = TRT_SHPS,
    breaks = TRT_LVLS,
    guide  = "none"
  )

  p <- ggplot(d, aes(x = Group_f, y = .data[[value_col]],
                     colour = TG, shape = ShapeKey, group = Treat_f)) +
    geom_point(position = position_jitterdodge(jitter.width = jitter_width,
                                               dodge.width = dodge,
                                               seed = 1),
               size = point_size, alpha = 0.85) +
    geom_errorbar(data = sdat,
                  aes(x = Group_f, y = mean, ymin = mean - SE, ymax = mean + SE,
                      group = Treat_f),
                  inherit.aes = FALSE,
                  position = position_dodge(dodge),
                  width = summary_width, linewidth = 0.55, colour = "black") +
    geom_crossbar(data = sdat,
                  aes(x = Group_f, y = mean, ymin = mean, ymax = mean,
                      group = Treat_f),
                  inherit.aes = FALSE,
                  position = position_dodge(dodge),
                  width = summary_width, linewidth = 0.55, colour = "black") +
    colour_scale + shape_scale +
    scale_x_discrete(labels = c("No choice"       = "No choice",
                                "Exercise choice" = "Exercise\nchoice"),
                     expand = expansion(add = c(x_expand_low, x_expand_high))) +
    scale_y_continuous(breaks = y_breaks) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.text.x  = element_text(size = x_text_size, face = "bold", colour = "black"),
      axis.title.y = element_text(size = ylab_size, face = "bold"),
      plot.caption = element_text(size = caption_size, face = caption_face, hjust = 0,
                                  lineheight = 1.3, margin = margin(t = caption_margin_t),
                                  colour = "grey25", family = FONT_CAP)
    )

  if (!is.null(cld)) {
    cld_d <- cld %>%
      mutate(
        Group_f = factor(Group, levels = GRP_LVLS),
        Treat_f = factor(Treatment, levels = TRT_LVLS),
        TG      = paste(Treatment, Group, sep = ":")
      )
    p <- p +
      geom_text(data = cld_d,
                aes(x = Group_f, y = y, label = label, group = Treat_f),
                inherit.aes = FALSE,
                position = position_dodge(dodge),
                colour = "black", size = 4.2, fontface = "bold")
  }

  if (!show_legend) {
    p <- p + guides(colour = "none", shape = "none")
  } else {
    p <- p + labs(colour = legend_title, shape = legend_title)
  }

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

summary_for_raw_plot <- function(data, value_col) {
  data %>%
    group_by(Treatment, Group) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      SE   = sd(.data[[value_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[value_col]]))),
      n    = sum(!is.na(.data[[value_col]])),
      .groups = "drop"
    )
}

make_cld_positions <- function(summary_data, labels, prop = 0.10) {
  rng <- diff(range(c(summary_data$mean - summary_data$SE,
                      summary_data$mean + summary_data$SE), na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) rng <- max(summary_data$mean, na.rm = TRUE) * 0.10
  summary_data %>%
    mutate(label = labels,
           y = mean + SE + prop * rng)
}

make_cld_positions_from_raw_max <- function(raw_data, value_col, labels, prop = 0.10) {
  d <- raw_data %>%
    group_by(Treatment, Group) %>%
    summarise(max_y = max(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  rng <- diff(range(raw_data[[value_col]], na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) rng <- max(raw_data[[value_col]], na.rm = TRUE) * 0.10
  d %>%
    mutate(label = labels,
           y = max_y + prop * rng)
}

# ── House rule for letters: sweep the cells in plot order (left to right) ─────
# and hand out "a", "b", "c" ... in order of first appearance, so the leftmost
# cell carries "a" whenever it has a single letter (a leftmost cell with several
# letters, e.g. "ab", keeps several). Only the letter NAMES change; which cells
# share a letter is untouched. `ord` = row indices in plot order.
relabel_cld <- function(cld, ord = seq_along(cld)) {
  cld <- trimws(as.character(cld))
  map <- character(0)
  for (i in ord) {
    if (is.na(cld[i])) next
    for (ch in sort(strsplit(cld[i], "", fixed = TRUE)[[1]]))
      if (!ch %in% names(map)) map[ch] <- letters[length(map) + 1]
  }
  vapply(cld, function(x) {
    if (is.na(x)) return(NA_character_)
    paste(sort(unname(map[strsplit(x, "", fixed = TRUE)[[1]]])), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

# ── Letters only where the panel's FULL interaction is significant ────────────
# A panel gets letters only if the highest-order interaction among the factors
# it displays is significant (p < STMT_ALPHA): Treatment × Group for a
# four-cell panel, Treatment × Group × interval for an eight-point interval
# panel, Group × interval for the marginal interval panels, and
# Group × Treatment × side / × segment for the aggression panels split by side
# or interval. A significant main effect or lower-order interaction alone
# (e.g. only Interval) gives no letters, even if pairwise comparisons differ.
interaction_sig <- function(dir, file, key, key_col = "variable") {
  p <- .p_from_stmt(.stmt_lookup(dir, file, key, "statement", key_col))
  is.finite(p) && p < STMT_ALPHA
}
agg_term_sig <- function(anova_file, term) {
  df <- .read_stmt_csv(file.path(AGG_ANOVA_DIR, anova_file))
  if (is.null(df) || !all(c("term", "Pr(>Chisq)") %in% names(df))) return(FALSE)
  p <- df$`Pr(>Chisq)`[df$term == term]
  length(p) == 1 && is.finite(p) && p < STMT_ALPHA
}

# ── Source Tukey CLD letters from a pipeline `cld/<key>.csv` ───────────────────
# Returns a data frame (Treatment, Group, cld) with Group recoded to GRP_LVLS,
# or NULL if the file is missing.
read_cld <- function(cld_dir, key) {
  f <- file.path(cld_dir, paste0(key, ".csv"))
  if (!file.exists(f)) return(NULL)
  # four-cell panels: letters only if Treatment × Group is significant
  sig <- if (identical(cld_dir, CLD_MONO_DIR))
    interaction_sig(ANOVA_STMT_DIR, "anova_stmt_TreatmentxGroup.csv", key, "key")
  else
    interaction_sig(BEHAV_ANOVA_DIR, "anova_stmt_trial_TreatmentxGroup.csv", sub("^trial_", "", key))
  if (!sig) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  if (!all(c("Treatment", "Group", "cld") %in% names(df))) return(NULL)
  df$Group <- dplyr::recode(df$Group,
                            "Exercise Choice" = "Exercise choice",
                            "No Choice"       = "No choice")
  # plot order of the trial/sample panels: Group on the x axis, Treatment dodged
  df$cld <- relabel_cld(df$cld, order(match(df$Group, GRP_LVLS), match(df$Treatment, TRT_LVLS)))
  df[, c("Treatment", "Group", "cld")]
}

# Build CLD label positions for a raw scatter, sourcing letters from `cld_df`.
# Letters are placed prop (=10%) of the total raw data range above each scatter
# group's maximum point. Returns NULL when all cells share the same letter
# (i.e. no significant Tukey differences) so that no letters are drawn.
cld_from_csv <- function(raw_data, value_col, cld_df, prop = 0.10) {
  if (is.null(cld_df) || length(unique(cld_df$cld)) <= 1) return(NULL)
  d <- raw_data %>%
    mutate(Treatment = as.character(Treatment), Group = as.character(Group)) %>%
    group_by(Treatment, Group) %>%
    summarise(max_y = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    left_join(cld_df, by = c("Treatment", "Group"))
  rng <- diff(range(raw_data[[value_col]], na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) rng <- max(raw_data[[value_col]], na.rm = TRUE) * 0.10
  d %>%
    mutate(Treatment = factor(Treatment, levels = TRT_LVLS),
           Group     = factor(Group, levels = GRP_LVLS),
           label     = cld,
           y         = max_y + prop * rng)
}

# Read a Group×interval marginal CLD from the marginal-emm/cld/ subfolder.
# Returns data frame (Group, interval [integer], cld) or NULL.
read_marginal_cld <- function(key, term_str = "Group:interval") {
  safe <- gsub(":", "x", term_str)
  f <- file.path(MARGINAL_CLD_DIR, paste0("interval_", key, "_", safe, ".csv"))
  if (!file.exists(f)) return(NULL)
  # marginal Group × interval panels: letters only if Group × interval is significant
  if (identical(term_str, "Group:interval") &&
      !interaction_sig(BEHAV_ANOVA_DIR, "anova_stmt_interval_GroupxInterval.csv", key)) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  df$Group <- dplyr::recode(df$Group,
                            "Exercise Choice" = "Exercise choice",
                            "No Choice"       = "No choice")
  df$interval <- as.integer(df$interval)
  df[, c("Group", "interval", "cld")]
}

# Remap CLD so "No choice / interval 2" (leftmost in plot) anchors on "a".
# Only swaps a↔b; "ab" cells are unchanged.
remap_cld_to_leftmost <- function(cld_df) {
  if (is.null(cld_df)) return(NULL)
  cld_df$cld <- relabel_cld(cld_df$cld,
                            order(cld_df$interval, match(cld_df$Group, GRP_LVLS)))
  cld_df
}

# ANOVA caption for behaviour interval model (Interval + Group×interval terms).
# Filtered / completed by finish_caption() (see STMT_MODE).
behav_int_caption <- function(key) {
  spec <- list(
    c("anova_stmt_interval_Interval.csv",       "Interval"),
    c("anova_stmt_interval_GroupxInterval.csv", "Group × interval")
  )
  out <- character(0)
  for (sp in spec) {
    f <- file.path(BEHAV_ANOVA_DIR, sp[1])
    if (!file.exists(f)) next
    df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
    s  <- df$statement[df$variable == key]
    if (length(s) == 0 || is.na(s[1])) next
    s <- s[1]
    out <- c(out, sub(paste0(sp[2], " "), paste0(sp[2], ": "), s, fixed = TRUE))
  }
  finish_caption(paste(out, collapse = "\n"), .n_behav(key, "interval"))
}

# Group-only colours for marginal interval plots (row 2 of Figure 5 extension).
GRP_COLS_ROW2   <- c("No choice" = "#3D98D0", "Exercise choice" = "#DE8400")
GRP_SHAPES_ROW2 <- c("No choice" = 17, "Exercise choice" = 19)

# Build a simple "Group" colour-only standalone legend for row-2 panels.
# item_gap and frame_ylim intentionally match build_combined_legend defaults so
# the physical distance between consecutive dot rows is identical in both legends.
build_group_legend <- function(scale = 1, item_gap = 0.5,
                               frame_ylim = c(0.5, 6.4)) {
  ts  <- 3.0 * scale; tts <- 3.4 * scale; ds <- 4.2 * scale
  lab_x <- 1 + 0.5 * scale
  ctr   <- mean(frame_ylim)
  # Three y positions using same item_gap as build_combined_legend
  ys    <- ctr + item_gap * c(1, 0, -1)
  dots  <- data.frame(x = 1, y = c(ys[2], ys[3]),
                      col = c("#3D98D0", "#DE8400"))
  labs_d <- data.frame(x = lab_x, y = c(ys[2], ys[3]),
                       lab = c("No Choice", "Exercise Choice"))
  title  <- data.frame(x = 0.7, y = ys[1], lab = "Group")
  ggplot() +
    geom_point(data = dots, aes(x, y), colour = dots$col, size = ds) +
    geom_text(data = labs_d, aes(x, y, label = lab), hjust = 0, size = ts) +
    geom_text(data = title, aes(x, y, label = lab), hjust = 0,
              fontface = "bold", size = tts) +
    coord_cartesian(xlim = c(0.4, lab_x + 5.5), ylim = frame_ylim, clip = "off") +
    theme_void()
}

# Line plot of Group×interval MARGINAL means (raw data, averaged over Treatment).
# Colors by Group (GRP_COLS_ROW2). CLD letters from model marginal-emm CSV.
# y_mult: apply unit conversion (e.g. proportion → seconds ×600).
interval_marginal_plot <- function(raw_data, value_col, ylabel,
                                   caption_text, cld_df = NULL,
                                   y_mult = 1, ylim = NULL,
                                   y_breaks = waiver(),
                                   ylab_size = PANEL5B_YLAB_SIZE,
                                   caption_size = CAPTION_SIZE_45,
                                   caption_margin_t = CAPTION_MARGIN_T,
                                   x_text_size = X_TICK_SIZE,
                                   dodge_width = INT_DODGE * 0.5) {
  d_sum <- raw_data %>%
    filter(!is.na(.data[[value_col]]), interval %in% c(2, 3)) %>%
    group_by(Group, interval) %>%
    summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE) * y_mult,
      SE_val   = sd(.data[[value_col]], na.rm = TRUE) /
                   sqrt(sum(!is.na(.data[[value_col]]))) * y_mult,
      .groups  = "drop"
    ) %>%
    mutate(
      Group_f = factor(Group, levels = GRP_LVLS),
      Int_f   = factor(interval, levels = c(2, 3))
    )

  dg <- position_dodge(width = dodge_width)

  p <- ggplot(d_sum, aes(x = Int_f, y = mean_val,
                         colour = Group_f, group = Group_f)) +
    geom_errorbar(aes(ymin = mean_val - SE_val, ymax = mean_val + SE_val),
                  position = dg, width = 0.10, linewidth = 0.5) +
    geom_line(position = dg, linewidth = 0.65) +
    geom_point(position = dg, size = 2.4 * 2, shape = 19) +
    scale_colour_manual(name = "Group", values = GRP_COLS_ROW2, breaks = GRP_LVLS) +
    scale_x_discrete(labels = c("2" = "20–30 min", "3" = "40–50 min"),
                     expand = expansion(add = 0.55)) +
    scale_y_continuous(breaks = y_breaks) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.title.y    = element_text(size = ylab_size, face = "bold"),
      axis.text.x     = element_text(size = x_text_size, face = "bold", colour = "black"),
      plot.caption    = element_text(size = caption_size, hjust = 0, lineheight = 1.3,
                                     margin = margin(t = caption_margin_t), colour = "grey25",
                                     family = FONT_CAP),
      legend.position = "none"
    )

  if (!is.null(cld_df) && length(unique(cld_df$cld)) > 1) {
    rng <- diff(range(d_sum$mean_val, na.rm = TRUE))
    if (!is.finite(rng) || rng == 0) rng <- max(d_sum$mean_val, na.rm = TRUE) * 0.1
    cld_pos <- d_sum %>%
      left_join(cld_df, by = c("Group", "interval")) %>%
      mutate(y_cld = mean_val + SE_val + 0.10 * rng)
    p <- p +
      geom_text(data = cld_pos,
                aes(x = Int_f, y = y_cld, label = cld, group = Group_f),
                inherit.aes = FALSE, position = dg,
                colour = "black", size = 4.2, fontface = "bold")
  }

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

trial_plot <- function(data, ylabel, caption_text, dodge = TG_DODGE,
                       ylim = NULL, y_breaks = waiver(),
                       ylab_size = 9.5, caption_size = 7.5,
                       caption_face = "plain", cld = NULL) {
  d <- data %>%
    mutate(
      Group_f = factor(Group, levels = GRP_LVLS),
      Treat_f = factor(Treatment, levels = TRT_LVLS),
      TG      = paste(Treatment, Group, sep = ":")
    )

  p <- ggplot(d, aes(x = Group_f, colour = TG, shape = Group_f, group = Treat_f)) +
    geom_point(aes(y = mean),
               position = position_dodge(dodge), size = POINT_SIZE) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE),
                  position = position_dodge(dodge),
                  width = SUMMARY_WIDTH, linewidth = 0.55, colour = "black") +
    geom_crossbar(aes(y = mean, ymin = mean, ymax = mean),
                  position = position_dodge(dodge),
                  width = SUMMARY_WIDTH, linewidth = 0.55, colour = "black") +
    SCALE_COL + SCALE_SHP_GRP +
    scale_x_discrete(labels = c("No choice"       = "No choice",
                                "Exercise choice" = "Exercise\nchoice")) +
    scale_y_continuous(breaks = y_breaks) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.text.x  = element_text(size = 8.5 * 3, face = "bold", colour = "black"),
      axis.title.y = element_text(size = ylab_size, face = "bold"),
      plot.caption = element_text(size = caption_size, face = caption_face, hjust = 0,
                                  lineheight = 1.3, margin = margin(t = caption_margin_t),
                                  colour = "grey25", family = FONT_CAP)
    )

  if (!is.null(cld)) {
    cld_d <- cld %>%
      mutate(
        Group_f = factor(Group, levels = GRP_LVLS),
        Treat_f = factor(Treatment, levels = TRT_LVLS),
        TG      = paste(Treatment, Group, sep = ":")
      )
    p <- p +
      geom_text(data = cld_d,
                aes(x = Group_f, y = y, label = label, group = Treat_f),
                inherit.aes = FALSE,
                position = position_dodge(dodge),
                colour = "black", size = 4.2, fontface = "bold")
  }

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# Read the Treatment×Group×interval Tukey CLD (cld/interval_<key>.csv): one letter
# per plotted point of interval_plot(). Returns (Treatment, Group, Interval, cld)
# or NULL when the file is missing or all cells share the same letters, so that
# letters are drawn only where cells differ.
read_interval_cld <- function(key) {
  f <- file.path(CLD_BEHAV_DIR, paste0("interval_", key, ".csv"))
  if (!file.exists(f)) return(NULL)
  # eight-point panels: letters only if Treatment × Group × interval is significant
  if (!interaction_sig(BEHAV_ANOVA_DIR, "anova_stmt_interval_TreatmentxGroupxInterval.csv", key))
    return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  if (!all(c("Treatment", "Group", "interval", "cld") %in% names(df))) return(NULL)
  df$Group <- dplyr::recode(df$Group, "Exercise Choice" = "Exercise choice",
                            "No Choice" = "No choice")
  df$cld <- trimws(df$cld)
  if (length(unique(df$cld)) <= 1) return(NULL)
  # plot order of interval_plot(): interval on the x axis, the four
  # Treatment:Group series dodged in the order TG_ORDER_INT
  df$cld <- relabel_cld(df$cld, order(as.integer(df$interval),
                                      match(paste(df$Treatment, df$Group, sep = ":"), TG_ORDER_INT)))
  data.frame(Treatment = df$Treatment, Group = df$Group,
             Interval = as.integer(df$interval), cld = df$cld,
             stringsAsFactors = FALSE)
}

# Order of the four series within an interval, shared by all manuscript figures:
# Control/No Choice, Exercise/No Choice, Control/Exercise Choice, Exercise/Exercise Choice.
TG_ORDER_INT <- c("Control:No choice", "Exercise:No choice",
                  "Control:Exercise choice", "Exercise:Exercise choice")

# ── Helper: interval-level line plot  ─────────────────────────────────────────
# x = interval (2 → 3); one line per Treatment:Group combo. Points/lines are
# horizontally dodged within each interval so Group/Treatment cells do not overlap.
interval_plot <- function(data, ylabel, caption_text, ylim = NULL, y_breaks = waiver(),
                          y_mult = 1, se_mult = 1, show_legend = TRUE,
                          legend_position = "right", ylab_size = 9.5 * 2,
                          x_text_size = 8.5 * 2.5, caption_size = 7.5 * 1.5,
                          caption_margin_t = 6, legend_title = "Legend",
                          show_shape_legend = FALSE, dodge_width = 0.58,
                          cld = NULL) {
  d <- data %>%
    mutate(
      TG     = factor(paste(Treatment, Group, sep = ":"), levels = TG_ORDER_INT),
      Int_f  = factor(Interval, levels = c(2, 3)),
      Grp    = factor(Group, levels = GRP_LVLS),
      Trt    = factor(Treatment, levels = TRT_LVLS),
      mean   = mean * y_mult,
      SE     = SE * se_mult
    )

  dodge_int <- position_dodge(width = dodge_width)

  p <- ggplot(d, aes(x = Int_f, y = mean, colour = TG, shape = Trt, group = TG)) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE),
                  position = dodge_int, width = 0.10, linewidth = 0.5) +
    geom_line(position = dodge_int, linewidth = 0.65) +
    geom_point(position = dodge_int, size = 2.4 * 2) +
    scale_colour_manual(
      name   = legend_title,
      values = TG_COLS,
      breaks = TG_BREAKS_PLOT,
      labels = TG_LBLS[TG_BREAKS_PLOT]
    ) +
    scale_shape_manual(
      name   = "Treatment",
      values = TRT_SHPS,
      breaks = TRT_LVLS,
      labels = TRT_LVLS,
      guide  = if (show_shape_legend) "legend" else "none"
    ) +
    scale_x_discrete(labels = c("2" = "20–30 min", "3" = "40–50 min"),
                     expand = expansion(add = 0.55)) +
    scale_y_continuous(breaks = y_breaks) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.title.y = element_text(size = ylab_size, face = "bold"),
      axis.text.x  = element_text(size = x_text_size, face = "bold", colour = "black"),
      plot.caption = element_text(size = caption_size, hjust = 0, lineheight = 1.3,
                                  margin = margin(t = caption_margin_t), colour = "grey25",
                                  family = FONT_CAP),
      legend.position = legend_position,
      legend.title = element_text(face = "bold")
    )

  if (!show_legend) {
    p <- p + guides(colour = "none", shape = "none")
  } else if (show_shape_legend) {
    p <- p + guides(
      colour = guide_legend(title = legend_title, order = 1),
      shape  = guide_legend(title = legend_title, order = 2)
    )
  } else {
    p <- p + guides(colour = guide_legend(title = legend_title, order = 1),
                    shape = "none")
  }

  # CLD letters (cld = output of read_interval_cld(); NULL when cells do not
  # differ): one letter above each point, 6 % of the plotted range above the
  # upper SE tip, dodged like the points.
  if (!is.null(cld) && length(unique(cld$cld)) > 1) {
    rng <- diff(range(c(d$mean - d$SE, d$mean + d$SE), na.rm = TRUE))
    cld_d <- d %>%
      mutate(Interval = as.integer(as.character(Interval))) %>%
      left_join(cld, by = c("Treatment", "Group", "Interval")) %>%
      filter(!is.na(cld)) %>%
      mutate(y_cld = mean + SE + 0.06 * rng)
    p <- p + geom_label(data = cld_d,
                        aes(x = Int_f, y = y_cld, label = cld, group = TG),
                        inherit.aes = FALSE, position = dodge_int,
                        colour = "black", size = 4.2, fontface = "bold",
                        fill = scales::alpha("white", 0.75), border.colour = NA,
                        label.padding = unit(0.08, "lines"))
  }

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

################################################################################
###  DATA  (raw means ± SE — original measurement scale, post-outlier removal)
################################################################################

# ── Monoamine DM  (rosner_20260608_170823 / raw_means_SE.md) ──────────────────
d_DM_ratio <- tribble(
  ~Treatment, ~Group,             ~mean,  ~SE,
  "Control",  "Exercise choice",  0.230,  0.017,
  "Control",  "No choice",        0.219,  0.014,
  "Exercise", "Exercise choice",  0.180,  0.018,
  "Exercise", "No choice",        0.260,  0.021
)
d_DM_HT5 <- tribble(
  ~Treatment, ~Group,             ~mean,   ~SE,
  "Control",  "Exercise choice", 19.853,  2.981,
  "Control",  "No choice",       13.881,  2.128,
  "Exercise", "Exercise choice", 14.526,  2.486,
  "Exercise", "No choice",       15.523,  1.990
)

# ══ Behaviour summary cells — DERIVED from the behaviour run, never hardcoded ══
#
# These eleven tables were `tribble()` literals until 2026-08-21: numbers typed
# in by hand from a specific run (the header used to read
# "output/20260608_105522"). Because they were static,
# re-running this script could not update Figures 5 and 6 — the panels rendered
# byte-identically no matter what the pipeline produced.
#
# The transcription itself was accurate: every one of the eleven tables was
# checked against the 2026-08-19 run and reproduced to the 3 decimal places the
# literals carried. The problem was never wrong numbers, it was frozen ones —
# they could only stay correct for as long as the upstream data never changed.
# When it did change, they froze a known artifact into the figures:
# `d_int_speed`'s Exercise/Exercise-choice interval-3 cell read
# 135.846 ± 121.966 cm/s, which was fish T28_S3_F3's tracking-glitch inflation
# (LOGS.md 2026-08-20). That single frozen cell is what the "⚠ ... inflated by
# one extreme individual" caption caveat was describing. With the STEP 1
# trajectory fix the same cell is 13.319 ± 0.713 — in line with every other cell.
#
# Deriving them from BEHAV_RUN_DIR's CSVs means repointing BEHAV_RUN_DIR is now
# sufficient to update every panel, and the figures cannot disagree with the
# statistics they sit next to. Provenance below was established empirically by
# reproducing each old hardcoded table from the CSVs, not by assumption.
#
# NOTE: SE here is the plain standard error across the panel's observational
# unit (fish for ind_*, trial for gd_*) — matching how the original numbers were
# produced. It is deliberately NOT the model-based SE; caption statistics come
# from the pipeline's own ANOVA tables further down.

.behav_csv_cache <- new.env(parent = emptyenv())

.read_behav_csv <- function(fname) {
  if (!exists(fname, envir = .behav_csv_cache, inherits = FALSE)) {
    p <- file.path(RAW_BEHAV_DIR, fname)
    if (!file.exists(p))
      stop("manu_graphs: required behaviour CSV not found: ", p,
           "\n  Check BEHAV_RUN_DIR points at a complete behaviour_models.R run.",
           call. = FALSE)
    assign(fname, readr::read_csv(p, show_col_types = FALSE, progress = FALSE),
           envir = .behav_csv_cache)
  }
  get(fname, envir = .behav_csv_cache)
}

# mean +/- SE per Treatment x Group [x Interval], from one column of one CSV.
summarise_cells <- function(fname, col, by_interval = FALSE) {
  df <- .read_behav_csv(fname)
  if (!col %in% names(df))
    stop("manu_graphs: column '", col, "' not found in ", fname,
         " (has: ", paste(names(df), collapse = ", "), ")", call. = FALSE)

  grp <- c("Treatment", "Group")
  if (by_interval) {
    ivcol <- intersect(c("timepoint", "interval"), names(df))
    if (length(ivcol) == 0L)
      stop("manu_graphs: no timepoint/interval column in ", fname, call. = FALSE)
    df$Interval <- as.numeric(df[[ivcol[1]]])
    grp <- c(grp, "Interval")
  }

  out <- df %>%
    dplyr::filter(!is.na(.data[[col]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      mean = mean(.data[[col]]),
      SE   = if (dplyr::n() > 1L) stats::sd(.data[[col]]) / sqrt(dplyr::n()) else NA_real_,
      .groups = "drop"
    )

  if (by_interval) out <- dplyr::filter(out, Interval %in% c(2, 3))
  if (nrow(out) == 0L)
    stop("manu_graphs: no rows after summarising ", col, " from ", fname, call. = FALSE)
  as.data.frame(out)
}

# ── Trial-level behaviour ─────────────────────────────────────────────────────
d_tr_polar   <- summarise_cells("gd_trial.csv",     "mean_polarization")
d_tr_flow    <- summarise_cells("ind_trial.csv",    "time_spent_zone_flow_prop")
d_tr_trans   <- summarise_cells("ind_trial.csv",    "n_zone_transitions")

# ── Interval-level individual (observational unit = fish) ─────────────────────
d_int_trans  <- summarise_cells("ind_interval.csv", "n_zone_transitions",        by_interval = TRUE)
d_int_dist   <- summarise_cells("ind_interval.csv", "distance_m",                by_interval = TRUE)
d_int_mov    <- summarise_cells("ind_interval.csv", "time_moving_TS1_prop",      by_interval = TRUE)
d_int_speed  <- summarise_cells("ind_interval.csv", "speed_cm_s_mean",           by_interval = TRUE)

# ── Interval-level school / group dynamics (observational unit = trial) ───────
d_int_pol    <- summarise_cells("gd_interval.csv",  "mean_polarization",         by_interval = TRUE)
d_int_nnd    <- summarise_cells("gd_interval.csv",  "mean_nnd_cm",               by_interval = TRUE)
d_int_hull   <- summarise_cells("gd_interval.csv",  "mean_hull_area_cm2",        by_interval = TRUE)
d_int_cspeed <- summarise_cells("gd_interval.csv",  "mean_centroid_speed_cm_s",  by_interval = TRUE)

################################################################################
###  ANOVA STATEMENTS — read from the pipeline's own ANOVA tables
###  Unicode: × = ×  χ = χ  ² = ²  subscripts: ₀-₉
################################################################################
#
# These were hardcoded strings until 2026-08-21, the same freeze-drift problem as
# the summary tables above but one layer worse: a frozen *caption* misreports the
# statistics next to a panel that has itself moved on. They had already gone
# stale — Ind_distance_m's Group × interval term read F₁,₈₀ = 11.646, p = 0.001
# in the script while the run it was supposedly describing gave 8.482, p = 0.005,
# and GD_mean_nnd_cm had even changed denominator df (F₁,₇ -> F₁,₈).
#
# Everything mechanical is now read from disk. Behaviour and monoamine runs write
# pre-formatted `statement` columns, so those are used verbatim (guaranteeing the
# caption and the report cannot disagree). The aggression run writes raw
# term/Chisq/Df/Pr tables, so those are formatted here with the same rules the
# upstream formatter uses.
#
# Term ORDER is canonical (main effects, then interactions by order). Every
# builder collects every term and hands the lines to finish_caption(), which
# adds the effect size and applies STMT_MODE (significant-only for the
# manuscript set; the unselected statements of HANDOFF.md S9 are kept in the
# all_statements set).

# Strip a UTF-8 BOM from the first column name if read.csv kept it.
.debom_names <- function(df) {
  names(df) <- sub("^﻿", "", names(df))
  names(df) <- sub("^..¿", "", names(df))
  df
}

.read_stmt_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8"),
    error = function(e) NULL
  )
  if (is.null(df)) NULL else .debom_names(df)
}

# numeric p out of "... p = 0.019" / "... p < 0.0001"
.p_from_stmt <- function(s) {
  if (length(s) == 0 || is.na(s) || !nzchar(s)) return(NA_real_)
  m <- regmatches(s, regexpr("p[[:space:]]*[<=][[:space:]]*[0-9.]+", s))
  if (!length(m)) return(NA_real_)
  suppressWarnings(as.numeric(sub(".*[<=][[:space:]]*", "", m)))
}

# chi-square value out of "... χ²₁ = 18.13, p ..." (NA for F statements).
.chisq_from_stmt <- function(s) {
  i <- regexpr("\u03c7\u00b2", s, fixed = TRUE)
  if (i < 0) return(NA_real_)
  rest <- substring(s, i)
  m <- regmatches(rest, regexpr("=[[:space:]]*[0-9.]+", rest))
  if (!length(m)) return(NA_real_)          # e.g. "< 0.001": w not computed
  suppressWarnings(as.numeric(sub("=[[:space:]]*", "", m)))
}

# Final step of every caption builder. `text`: statements, one per line;
# `n`: number of observations of the model (for Cohen's w).
#  1. effect size on every statement: F statements already carry partial
#     eta-squared; chi-square statements get Cohen's w = sqrt(chi-square / n);
#  2. STMT_MODE "significant": drop the lines with p >= STMT_ALPHA (lines with
#     no p-value, e.g. prose notes, are kept). Returns "" if nothing is left.
finish_caption <- function(text, n = NA_real_) {
  if (is.null(text) || !length(text) || !nzchar(text)) return("")
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- vapply(lines, function(l) {
    chi <- .chisq_from_stmt(l)
    if (!is.finite(chi) || grepl(", w = ", l, fixed = TRUE)) return(l)
    if (!is.finite(n) || n <= 0) {
      warning("no N for Cohen's w: ", l, call. = FALSE); return(l)
    }
    paste0(l, ", w = ", fmt_es3(sqrt(chi / n)))
  }, character(1), USE.NAMES = FALSE)
  if (identical(STMT_MODE, "significant")) {
    p <- vapply(lines, .p_from_stmt, numeric(1), USE.NAMES = FALSE)
    lines <- lines[is.na(p) | p < STMT_ALPHA]
  }
  paste(lines, collapse = "\n")
}

# Number of observations behind each model (for Cohen's w).
.n_behav <- function(key, level = "interval") {
  f <- file.path(RAW_BEHAV_DIR, paste0(if (startsWith(key, "GD_")) "gd_" else "ind_", level, ".csv"))
  df <- .read_stmt_csv(f); col <- sub("^(GD|Ind)_", "", key)
  if (is.null(df) || !col %in% names(df)) return(NA_real_)
  sum(!is.na(df[[col]]))
}
.n_mono <- function(key) {
  df <- .read_stmt_csv(file.path(ANOVA_STMT_DIR, "model_decisions_assumptions.csv"))
  if (is.null(df) || !all(c("variable", "unit", "n") %in% names(df))) return(NA_real_)
  unit <- if (startsWith(key, "DM_")) "Monoamine DM" else if (startsWith(key, "RN_")) "Monoamine RN" else "Plasma"
  n <- df$n[df$variable == sub("^(DM|RN)_", "", key) & df$unit == unit]
  if (length(n)) as.numeric(n[1]) else NA_real_
}
# Aggression: rows of the table each model was fitted to.
.AGG_N_FILE <- c("ANOVA_R1a.csv" = "counts_trial_segment_side.csv",
                 "ANOVA_R1b.csv" = "counts_trial_side.csv",
                 "ANOVA_R2b.csv" = "probs_trial.csv",
                 "ANOVA_R3b.csv" = "crowding_events.csv")
.n_agg <- function(anova_file) {
  if (!anova_file %in% names(.AGG_N_FILE)) return(NA_real_)
  df <- .read_stmt_csv(file.path(AGG_DATA_DIR, .AGG_N_FILE[[anova_file]]))
  if (is.null(df)) NA_real_ else nrow(df)
}

# "Interval F subscripted terms" -> "Interval: F subscripted terms" (colon after the label).
.tidy_stmt <- function(s) sub("^(.*?) (F|χ²)", "\\1: \\2", s)

# One pre-formatted statement for `key` from a statement-table CSV.
.stmt_lookup <- function(dir, file, key, col = "statement", key_col = "variable") {
  df <- .read_stmt_csv(file.path(dir, file))
  if (is.null(df) || !all(c(key_col, col) %in% names(df))) return(NA_character_)
  s <- df[[col]][df[[key_col]] == key]
  if (length(s) == 0 || is.na(s[1]) || !nzchar(s[1])) return(NA_character_)
  s[1]
}

# Assemble a caption from an ordered list of (file, column) specs.
.build_caption <- function(dir, specs, key, key_col = "variable", extra = NULL, n = NA_real_) {
  out <- character(0)
  for (sp in specs) {
    s <- .stmt_lookup(dir, sp[[1]], key, sp[[2]], key_col)
    if (is.na(s)) next
    out <- c(out, .tidy_stmt(s))
  }
  out <- finish_caption(paste(out, collapse = "\n"), n)
  paste(c(if (nzchar(out)) out, extra), collapse = "\n")
}

# Canonical term order: main effects first, then interactions.
.BEHAV_INT_SPECS <- list(
  list("anova_stmt_interval_main_effects.csv",           "stmt_Treatment"),
  list("anova_stmt_interval_main_effects.csv",           "stmt_Group"),
  list("anova_stmt_interval_Interval.csv",               "statement"),
  list("anova_stmt_interval_TreatmentxGroup.csv",        "statement"),
  list("anova_stmt_interval_GroupxInterval.csv",         "statement"),
  list("anova_stmt_interval_TreatmentxGroupxInterval.csv","statement")
)
.BEHAV_TRI_SPECS <- list(
  list("anova_stmt_trial_main_effects.csv",    "stmt_Treatment"),
  list("anova_stmt_trial_main_effects.csv",    "stmt_Group"),
  list("anova_stmt_trial_TreatmentxGroup.csv", "statement")
)
.MONO_SPECS <- list(
  list("anova_stmt_Treatment.csv",        "statement"),
  list("anova_stmt_Group.csv",            "statement"),
  list("anova_stmt_TreatmentxGroup.csv",  "statement")
)

behav_stmt <- function(key, level = c("interval", "trial"), extra = NULL) {
  level <- match.arg(level)
  .build_caption(BEHAV_ANOVA_DIR,
                 if (level == "interval") .BEHAV_INT_SPECS else .BEHAV_TRI_SPECS,
                 key, key_col = "variable", extra = extra, n = .n_behav(key, level))
}
mono_stmt <- function(key) .build_caption(ANOVA_STMT_DIR, .MONO_SPECS, key, key_col = "key",
                                          n = .n_mono(key))

# ── Aggression: raw term/Chisq/Df/Pr tables, formatted here ───────────────────
.AGG_TERM_LABELS <- c(
  "group"                     = "Group",
  "treatment_clean"           = "Treatment",
  "side"                      = "Side",
  "segment"                   = "Interval",
  "crowding_diff_z"           = "Crowding difference",
  "group:treatment_clean"     = "Group × Treatment",
  "group:side"                = "Group × Side",
  "treatment_clean:side"      = "Treatment × Side"
)

# .to_subscript() retired -- to_subscript_statement() is sourced from
# house_style.R (A7); agg_stmt() below now routes through
# anova_statement() instead of formatting the chisq line by hand.

# Wald chi-square carries no eta2p (C6 [DR], structural via effect_size() in
# house_style.R); finish_caption() adds Cohen's w and applies STMT_MODE.
agg_stmt <- function(anova_file, terms) {
  df <- .read_stmt_csv(file.path(AGG_ANOVA_DIR, anova_file))
  if (is.null(df) || !all(c("term", "Chisq", "Df", "Pr(>Chisq)") %in% names(df)))
    return("")
  # `terms` is the hand-picked list for the manuscript set; the all_statements
  # set shows every term of the model (intercept excluded).
  if (identical(STMT_MODE, "all")) terms <- setdiff(df$term, "(Intercept)")
  out <- character(0)
  for (tm in terms) {
    row <- df[df$term == tm, , drop = FALSE]
    if (nrow(row) == 0) next
    p <- row$`Pr(>Chisq)`[1]
    lab <- if (tm %in% names(.AGG_TERM_LABELS)) .AGG_TERM_LABELS[[tm]] else
      paste(vapply(strsplit(tm, ":", fixed = TRUE)[[1]], function(x)
        if (x %in% names(.AGG_TERM_LABELS)) .AGG_TERM_LABELS[[x]] else x, character(1)),
        collapse = " × ")
    out <- c(out, anova_statement(lab, "chisq", stat = row$Chisq[1], df1 = row$Df[1], p = p))
  }
  finish_caption(paste(out, collapse = "\n"), .n_agg(anova_file))
}

# ── Monoamine (Figure 8) ──────────────────────────────────────────────────────
stmt_ratio <- mono_stmt("DM_ratio_5HIAA_5HT")
stmt_HT5   <- mono_stmt("DM_HT5")

# ── Trial-level behaviour ─────────────────────────────────────────────────────
stmt_tr_polar <- behav_stmt("GD_mean_polarization",          "trial")
stmt_tr_flow  <- behav_stmt("Ind_time_spent_zone_flow_prop", "trial")
stmt_tr_trans <- behav_stmt("Ind_n_zone_transitions",        "trial")

# ── Interval-level behaviour ──────────────────────────────────────────────────
stmt_int_trans <- behav_stmt("Ind_n_zone_transitions")
stmt_int_dist  <- behav_stmt("Ind_distance_m")
stmt_int_mov   <- behav_stmt("Ind_time_moving_TS1_prop")
# The "inflated by one extreme individual" caveat that used to be appended here
# described fish T28_S3_F3's tracking artifact, removed by the STEP 1 trajectory
# fix (LOGS.md 2026-08-20). The cell it referred to is now 13.319 +/- 0.713 cm/s,
# in line with every other cell, so the caveat is dropped rather than carried
# forward describing a problem that no longer exists.
stmt_int_speed <- behav_stmt("Ind_speed_cm_s_mean")
stmt_int_pol   <- behav_stmt("GD_mean_polarization")
stmt_int_nnd   <- behav_stmt("GD_mean_nnd_cm")
stmt_int_hull  <- behav_stmt("GD_mean_hull_area_cm2")
stmt_int_cspd  <- behav_stmt("GD_mean_centroid_speed_cm_s")

################################################################################
###  FIGURE 8  —  Monoamine DM: 5-HT (A) + 5-HIAA (B) + 5-HIAA/5-HT ratio (C)
################################################################################
cat("Building Figure 8...\n")

# Figure 4 uses raw DM monoamine rows from MONO_DM_CSV. Outliers are removed to
# match the post-Rosner summaries used in the manuscript/statistical output:
#   [5-HT]: Sample_ID 74 and 116
#   [5-HIAA]/[5-HT]: Sample_ID 31, 74 and 116
mono_DM_raw <- read.csv(MONO_DM_CSV, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Group = case_when(
      Group == "Exercise Choice" ~ "Exercise choice",
      Group == "No Choice"       ~ "No choice",
      TRUE                       ~ Group
    ),
    Treatment = factor(Treatment, levels = TRT_LVLS),
    Group     = factor(Group, levels = GRP_LVLS),
    HT5             = as.numeric(HT5),
    HIAA5           = as.numeric(HIAA5),
    ratio_5HIAA_5HT = as.numeric(ratio_5HIAA_5HT)
  )

# The dropped samples are read from the run's own outlier report, so the panels
# follow whichever monoamine run was discovered (the IDs used to be hardcoded
# for the June data and no longer matched the current data file).
mono_outliers <- function(indicator, module = "Monoamine DM") {
  f <- file.path(MONO_RUN_DIR, "reports", "outlier_report.md")
  if (!file.exists(f)) return(character(0))
  parts <- strsplit(readLines(f, encoding = "UTF-8", warn = FALSE), "|", fixed = TRUE)
  parts <- parts[lengths(parts) > 4]
  ids <- vapply(parts, function(p) {
    p <- trimws(p)
    if (p[2] == module && p[3] == indicator && grepl("^[0-9]+$", p[4])) p[4] else NA_character_
  }, character(1))
  unique(ids[!is.na(ids)])
}
fig4_HT5_raw   <- mono_DM_raw %>% filter(!Sample_ID %in% mono_outliers("5-HT"))
fig4_HIAA5_raw <- mono_DM_raw %>% filter(!Sample_ID %in% mono_outliers("5-HIAA"))
fig4_ratio_raw <- mono_DM_raw %>% filter(!Sample_ID %in% mono_outliers("5-HIAA/5-HT"))
# Fixed tick positions only where the data still fall in the range they were
# written for; otherwise ggplot picks them (fixed ticks piled up on each other
# when the concentration scale of the data file changed).
.breaks_or_auto <- function(x, breaks)
  if (max(x, na.rm = TRUE) > max(breaks) * 1.25) waiver() else breaks

# Panels: A = [5-HT], B = [5-HIAA], C = [5-HIAA]/[5-HT].
# CLD letters sourced from pipeline Tukey post-hoc (cld/<key>.csv); NULL when
# all cells share one letter. [5-HT] and [5-HIAA] are all "a" -> no letters.
# Ratio CLD: manual remap so both Control cells = "a" (Control No choice was
# "ab", Control Exercise choice was "ab"); Exercise No choice stays "b".
cld_ht5   <- cld_from_csv(fig4_HT5_raw,   "HT5",
                          read_cld(CLD_MONO_DIR, "DM_HT5"))
cld_hiaa5 <- cld_from_csv(fig4_HIAA5_raw, "HIAA5",
                          read_cld(CLD_MONO_DIR, "DM_HIAA5"))
cld_ratio_raw <- read_cld(CLD_MONO_DIR, "DM_ratio_5HIAA_5HT")
# read_cld() already relabels the letters in plot order (leftmost cell = "a"
# when it has one letter), replacing the manual a↔b swap that used to sit here.
cld_ratio <- cld_from_csv(fig4_ratio_raw, "ratio_5HIAA_5HT", cld_ratio_raw)

p4a <- raw_trial_plot(
  fig4_HT5_raw,
  value_col = "HT5",
  ylabel = "[5-HT] (ng/mg)",
  caption_text = stmt_HT5,
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  y_breaks = .breaks_or_auto(fig4_HT5_raw$HT5, seq(0, 40, by = 5)),
  ylab_size = 9.5 * 3 * 0.5,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_ht5,
  show_legend = FALSE
)

p4b <- raw_trial_plot(
  fig4_HIAA5_raw,
  value_col = "HIAA5",
  ylabel = "[5-HIAA] (ng/mg)",
  caption_text = anova_caption("DM_HIAA5"),
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  ylab_size = 9.5 * 3 * 0.5,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_hiaa5,
  show_legend = FALSE
)

p4c <- raw_trial_plot(
  fig4_ratio_raw,
  value_col = "ratio_5HIAA_5HT",
  ylabel = "[5-HIAA]/[5-HT]",
  caption_text = stmt_ratio,
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  y_breaks = seq(0, 0.40, by = 0.05),
  ylab_size = 9.5 * 3 * 0.5,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_ratio,
  show_legend = FALSE
) +
  coord_cartesian(ylim = c(NA, max(0.40,
    if (!is.null(cld_ratio)) max(cld_ratio$y, na.rm = TRUE) + 0.02 else 0.40)))

fig4 <- wrap_plots(
  p4a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p4b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  p4c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(),
  nrow = 1, widths = c(1, 1, 1, 0.5)
)
save_fig(fig4, "Figure_8", w = 17, h = 5.5)

# ── Figure_8_ver  (vertical: A / B / C stacked, legend at bottom) ────────────
fig4_ver <- wrap_plots(
  p4a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p4b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  p4c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(item_gap = 1.0),   # legend row tall enough for its six lines
  design = "A\nB\nC\nL", heights = c(1, 1, 1, 0.75)
)
save_fig(fig4_ver, "Figure_8_ver", w = 5.5, h = 16)

################################################################################
###  FIGURE 9  —  Plasma: lactate (A) + cortisol (B)
################################################################################
# NOTE: built after the monoamine block only because of code flow.
# New manuscript numbering: 4 behaviour, 5 interval-individual, 6 interval-school,
# [7 reserved for aggression], 8 brain neurochemistry, 9 plasma.
cat("Building Figure 9...\n")

# Mirrors every aesthetic of the monoamine scatter figure (raw points + black
# mean ± SE crossbar, Treatment shapes, tight x-expansion, custom legend, tags)
# but plots plasma lactate (Panel A) and plasma cortisol (Panel B). ANOVA
# captions are sourced from ANOVA_STMT_DIR, keeping statements with p < 0.09.
plasma_raw <- read.csv(PLASMA_CSV, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    Group = case_when(
      Group == "Exercise Choice" ~ "Exercise choice",
      Group == "No Choice"       ~ "No choice",
      TRUE                       ~ Group
    ),
    Treatment       = factor(Treatment, levels = TRT_LVLS),
    Group           = factor(Group, levels = GRP_LVLS),
    Lactate_mmol_L  = as.numeric(Lactate_mmol_L),
    plasma_cortisol = as.numeric(plasma_cortisol)
  )

# Sourced Tukey CLD (both are all-"a" -> NULL -> no letters drawn).
cld_lactate  <- cld_from_csv(plasma_raw, "Lactate_mmol_L",
                             read_cld(CLD_MONO_DIR, "Lactate_mmol_L"))
cld_cortisol <- cld_from_csv(plasma_raw, "plasma_cortisol",
                             read_cld(CLD_MONO_DIR, "plasma_cortisol"))

p4r_a <- raw_trial_plot(
  plasma_raw,
  value_col = "Lactate_mmol_L",
  ylabel = "Plasma lactate (mmol/L)",
  caption_text = anova_caption("Lactate_mmol_L"),
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  ylab_size = 9.5 * 3 * 0.5,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_lactate,
  show_legend = FALSE
)

p4r_b <- raw_trial_plot(
  plasma_raw,
  value_col = "plasma_cortisol",
  ylabel = "Plasma cortisol (ng/mL)",
  caption_text = anova_caption("plasma_cortisol"),
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  ylab_size = 9.5 * 3 * 0.5,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_cortisol,
  show_legend = FALSE
)

fig4_real <- wrap_plots(
  p4r_a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p4r_b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(),
  nrow = 1, widths = c(1, 1, 0.55)
)
save_fig(fig4_real, "Figure_9", w = 13, h = 5.5)

# ── Figure_9_ver  (vertical: A / B stacked, legend at bottom) ────────────────
fig4_real_ver <- wrap_plots(
  p4r_a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p4r_b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(item_gap = 1.0),   # legend row tall enough for its six lines
  design = "A\nB\nL", heights = c(1, 1, 0.55)
)
save_fig(fig4_real_ver, "Figure_9_ver", w = 5.5, h = 11)

################################################################################
###  FIGURE 4  —  Trial behaviour: Group polarisation + Flow-zone time + Transitions
################################################################################
# SKIP Figure 4 if CLD_BEHAV_DIR (behavior analysis outputs) is unavailable
if (dir.exists(CLD_BEHAV_DIR)) {
cat("Building Figure 4...\n")

# Read raw behaviour datasets.
# Panel A uses gd_trial.csv because polarisation is a group/school-level metric.
# Panels B–C use ind_interval.csv aggregated across intervals 2 and 3 to obtain
# one plotted point per fish/trial (master_fish_id).
gd_trial_raw <- read.csv(GD_TRIAL_CSV, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    trial_id = as.character(trial_id),
    Group = case_when(
      Group == "Exercise Choice" ~ "Exercise choice",
      Group == "No Choice"       ~ "No choice",
      TRUE                       ~ Group
    ),
    Treatment = factor(Treatment, levels = TRT_LVLS),
    Group     = factor(Group, levels = GRP_LVLS),
    mean_polarization = as.numeric(mean_polarization)
  )

ind_interval_raw <- read.csv(IND_INTERVAL_CSV, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    trial_id = as.character(trial_id),
    Group = case_when(
      Group == "Exercise Choice" ~ "Exercise choice",
      Group == "No Choice"       ~ "No choice",
      TRUE                       ~ Group
    ),
    Treatment = factor(Treatment, levels = TRT_LVLS),
    Group     = factor(Group, levels = GRP_LVLS),
    interval  = as.numeric(interval),
    time_spent_zone_flow_prop = as.numeric(time_spent_zone_flow_prop),
    n_zone_transitions = as.numeric(n_zone_transitions),
    master_fish_id = ifelse(is.na(master_fish_id) | master_fish_id == "",
                            as.character(fish_id), as.character(master_fish_id))
  )

fig5_polar_raw <- gd_trial_raw %>%
  transmute(
    unit_id = paste0("trial_", trial_id),
    Treatment, Group,
    polarisation = mean_polarization
  )

fig5_individual_raw <- ind_interval_raw %>%
  filter(interval %in% c(2, 3)) %>%
  group_by(Treatment, Group, trial_id, master_fish_id) %>%
  summarise(
    # Requested conversion: proportion × 2 × 600.
    # At fish level this is implemented as sum(interval proportion × 600 s)
    # across intervals 2 and 3, equivalent to mean proportion × 2 × 600.
    time_flow_s = mean(time_spent_zone_flow_prop, na.rm = TRUE) * 2 * 600,
    flow_calm_transitions = sum(n_zone_transitions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(unit_id = paste0(trial_id, "_", master_fish_id))

# Panel CLD letters sourced from the behaviour Tukey post-hoc (trial level,
# Treatment×Group). Every behaviour interaction is non-significant, so these are
# all-"a" -> NULL -> no letters drawn. CLD_BEHAV_DIR does not exist until
# behaviour_models.R is re-run (its upstream efficient_scripts/ inputs
# are currently unavailable), in which case read_cld() returns NULL.
cld_polar     <- cld_from_csv(fig5_polar_raw, "polarisation",
                              read_cld(CLD_BEHAV_DIR, "trial_GD_mean_polarization"))
cld_time_flow <- cld_from_csv(fig5_individual_raw, "time_flow_s",
                              read_cld(CLD_BEHAV_DIR, "trial_Ind_time_spent_zone_flow_prop"))
cld_trans     <- cld_from_csv(fig5_individual_raw, "flow_calm_transitions",
                              read_cld(CLD_BEHAV_DIR, "trial_Ind_n_zone_transitions"))

p5a <- raw_trial_plot(
  fig5_polar_raw,
  value_col = "polarisation",
  ylabel = "Group polarisation.",
  caption_text = stmt_tr_polar,
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  ylab_size = PANEL5B_YLAB_SIZE,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_polar,
  show_legend = FALSE
) +
  coord_cartesian(ylim = c(0.45, 0.50))

p5b <- raw_trial_plot(
  fig5_individual_raw,
  value_col = "time_flow_s",
  ylabel = "Time spent in flow zone (s)",
  caption_text = stmt_tr_flow,
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  ylab_size = PANEL5B_YLAB_SIZE,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_time_flow,
  y_breaks = seq(0, 1200, by = 200),   # max explicit y tick = 1200 (Panel 5B)
  show_legend = FALSE
) +
  coord_cartesian(
    ylim = c(0, max(1200,
      if (!is.null(cld_time_flow)) max(cld_time_flow$y, na.rm = TRUE) * 1.02 else 1200)))

p5c <- raw_trial_plot(
  fig5_individual_raw,
  value_col = "flow_calm_transitions",
  ylabel = "N° of flow-calm transitions",
  caption_text = stmt_tr_trans,
  dodge = FIG45_DODGE,
  jitter_width = FIG45_JITTER_WIDTH,
  summary_width = FIG45_SUMMARY_WIDTH,
  x_text_size = X_TICK_SIZE,
  x_expand_low = FIG45_X_EXPAND_LOW,
  x_expand_high = FIG45_X_EXPAND_HIGH,
  caption_margin_t = CAPTION_MARGIN_T,
  ylab_size = PANEL5B_YLAB_SIZE,
  caption_size = CAPTION_SIZE_45,
  caption_face = "plain",
  cld = cld_trans,
  show_legend = FALSE
)

fig5 <- wrap_plots(
  p5a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p5b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  p5c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(),
  nrow = 1, widths = c(1, 1, 1, 0.5)
)
save_fig(fig5, "Figure_4", w = 13, h = 5.5)

# ── Figure_4_ver  (vertical: A / B / C stacked, legend at bottom) ────────────
fig5_ver <- wrap_plots(
  p5a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p5b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  p5c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(item_gap = 1.0),   # legend row tall enough for its six lines
  design = "A\nB\nC\nL", heights = c(1, 1, 1, 0.75)
)
save_fig(fig5_ver, "Figure_4_ver", w = 5.5, h = 16)
} else {
  cat("Skipping Figure 4 (behavior analysis output CLD_BEHAV_DIR not available).\n")
}

################################################################################
###  FIGURE 5  —  Interval individual: transitions / distance / moving
################################################################################
if (dir.exists(CLD_BEHAV_DIR)) {
cat("Building Figure 5...\n")

# Figure 6 is displayed horizontally, with the only legend shown in panel C
# on the left. Y-axis title sizes are matched to Figure 5B.
p6a <- interval_plot(
  d_int_trans,
  ylabel = "N° of flow-calm transitions",
  caption_text = stmt_int_trans,
  show_legend = FALSE,
  ylab_size = PANEL5B_YLAB_SIZE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("Ind_n_zone_transitions")
)

p6b <- interval_plot(
  d_int_dist,
  ylabel = "Distance (m)",
  caption_text = stmt_int_dist,
  y_breaks = seq(0, 140, by = 10),
  show_legend = FALSE,
  ylab_size = PANEL5B_YLAB_SIZE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("Ind_distance_m")
)

p6c <- interval_plot(
  d_int_mov,
  ylabel = "Time spent moving\n>1 BL s\u207b\u00b9 (s)",
  caption_text = stmt_int_mov,
  y_mult = 600,
  se_mult = 600,
  show_legend = FALSE,
  ylab_size = PANEL5B_YLAB_SIZE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("Ind_time_moving_TS1_prop")
)

# Panel D / speed is temporarily dropped from Figure 6.
# p6d <- interval_plot(d_int_speed, ylabel = "Mean speed (cm/s)",
#                      caption_text = stmt_int_speed,
#                      ylim = c(0, 35))

# Panel letters matched to Figure 5 (size 15), legend to the right of panel C.
fig6 <- wrap_plots(
  p6a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p6b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  p6c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(),
  nrow = 1, widths = c(1, 1, 1, 0.45)
)
save_fig(fig6, "Figure_5", w = 16, h = 4.8)

# ── Figure_5_ext  (2-row: row 1 = existing interval scatter; row 2 = Group×interval emmeans)
# Row 2 uses marginal means averaged over Treatment, coloured by Group only.
# CLD sourced from marginal-emm/cld/, remapped so No-choice/interval-2 = "a".
# ANOVA captions: Interval + Group×interval terms, p < 0.05 only.

# Raw individual data already in ind_interval_raw (loaded in Figure 4 block).
# n_zone_transitions is per fish per interval; distance_m is per fish per interval;
# time_moving_TS1_prop is a proportion — multiply by 600 for seconds.

cld_trans_marg  <- remap_cld_to_leftmost(
  read_marginal_cld("Ind_n_zone_transitions",   "Group:interval"))
cld_dist_marg   <- remap_cld_to_leftmost(
  read_marginal_cld("Ind_distance_m",            "Group:interval"))
cld_mov_marg    <- remap_cld_to_leftmost(
  read_marginal_cld("Ind_time_moving_TS1_prop",  "Group:interval"))

p6a_emm <- interval_marginal_plot(
  ind_interval_raw, "n_zone_transitions",
  ylabel       = "N° of flow-calm transitions",
  caption_text = behav_int_caption("Ind_n_zone_transitions"),
  cld_df       = cld_trans_marg
)
p6b_emm <- interval_marginal_plot(
  ind_interval_raw, "distance_m",
  ylabel       = "Distance (m)",
  caption_text = behav_int_caption("Ind_distance_m"),
  cld_df       = cld_dist_marg
)
p6c_emm <- interval_marginal_plot(
  ind_interval_raw, "time_moving_TS1_prop",
  ylabel       = "Time spent moving\n>1 BL s\u207b\u00b9 (s)",
  caption_text = behav_int_caption("Ind_time_moving_TS1_prop"),
  cld_df       = cld_mov_marg,
  y_mult       = 600
)

fig6_ext <- wrap_plots(
  # Row 1 (A-D): interval scatter same as Figure_5; D = combined legend
  p6a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE),
  p6b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE),
  p6c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE),
  build_combined_legend(),
  # Row 2 (E-H): Group×interval marginal means; H = group-only legend
  p6a_emm + labs(tag = "D") + tag_theme(FIG45_TAG_SIZE),
  p6b_emm + labs(tag = "E") + tag_theme(FIG45_TAG_SIZE),
  p6c_emm + labs(tag = "F") + tag_theme(FIG45_TAG_SIZE),
  build_group_legend(),
  # patchwork assigns A→1st, B→2nd … H→8th element in alphabetical order:
  # row1: A=p6a | B=p6b | C=p6c | D=combined_legend
  # row2: E=p6a_emm | F=p6b_emm | G=p6c_emm | H=group_legend
  design = "ABCD\nEFGH",
  widths = c(1, 1, 1, 0.45)
)
save_fig(fig6_ext, "Figure_5_ext", w = 16, h = 9.6)
} else {
  cat("Skipping Figure 5 (behavior analysis output not available).\n")
}

################################################################################
###  FIGURE 6  —  Interval school: polarisation / NND / hull area / speed
################################################################################
if (dir.exists(CLD_BEHAV_DIR)) {
cat("Building Figure 6...\n")

p7a <- interval_plot(
  d_int_pol,
  ylabel = "Group polarisation.",
  caption_text = stmt_int_pol,
  show_legend = FALSE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  caption_size = FIG7_CAPTION_SIZE,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("GD_mean_polarization"),
  ylab_size = PANEL5B_YLAB_SIZE
)
p7b <- interval_plot(
  d_int_nnd,
  ylabel = "Mean NND (cm)",
  caption_text = stmt_int_nnd,
  show_legend = FALSE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  caption_size = FIG7_CAPTION_SIZE,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("GD_mean_nnd_cm"),
  ylab_size = PANEL5B_YLAB_SIZE
)
p7c <- interval_plot(
  d_int_hull,
  ylabel = "School area (cm²)",
  caption_text = stmt_int_hull,
  show_legend = FALSE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  caption_size = FIG7_CAPTION_SIZE,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("GD_mean_hull_area_cm2"),
  ylab_size = PANEL5B_YLAB_SIZE
)
p7d <- interval_plot(
  d_int_cspeed,
  ylabel = "School speed (cm/s)",
  caption_text = stmt_int_cspd,
  show_legend = FALSE,
  x_text_size = X_TICK_SIZE,
  caption_margin_t = CAPTION_MARGIN_T,
  caption_size = FIG7_CAPTION_SIZE,
  dodge_width = INT_DODGE,
  cld = read_interval_cld("GD_mean_centroid_speed_cm_s"),
  ylab_size = PANEL5B_YLAB_SIZE
)

# (Note: FIG7_TAG_SIZE, FIG7_TAG_Y, and PANEL_BORDER are defined at the top with other constants.)

# Panel tags at the Figure 5 size: the larger FIG7_TAG_SIZE sat on the top y tick label.
fig7 <- wrap_plots(
  p7a + labs(tag = "A") + tag_theme(FIG45_TAG_SIZE) + PANEL_BORDER,
  p7b + labs(tag = "B") + tag_theme(FIG45_TAG_SIZE) + PANEL_BORDER,
  p7c + labs(tag = "C") + tag_theme(FIG45_TAG_SIZE) + PANEL_BORDER,
  p7d + labs(tag = "D") + tag_theme(FIG45_TAG_SIZE) + PANEL_BORDER,
  build_combined_legend(frame_ylim = c(-2.4, 9.3)),
  design = "ABE\nCDE", widths = c(1, 1, 0.45)
)
save_fig(fig7, "Figure_6", w = 14, h = 9)
} else {
  cat("Skipping Figure 6 (behavior analysis output not available).\n")
}

################################################################################
###  FIGURE 7  —  Aggression: rate interaction (A) + robustness forest (B) +
###               side preference (C) + crowding mechanism (D)
###  Source: aggression_models.R latest run (AGG_RUN_DIR).
###    A  cld/R1a_cld.csv        (NB GLMM emmeans, events / 600 s, by side)
###    B  post-hoc/R1*_pairs     (No-choice Control / Exercise IRR, flow side)
###    C  cld/R2b_cld.csv        (binomial GLMM emmeans, P(flow), trial scope)
###    D  datasets/crowding_events.csv (raw side_flow vs crowding difference)
################################################################################
cat("Building Figure 7...\n")

# ── Read an aggression emmeans CLD (cld/<key>_cld.csv) ─────────────────────────
# Aggression CLD columns differ from the behaviour/mono ones:
#   group, treatment_clean, [side|], <value_col>, SE, df, asymp.LCL, asymp.UCL, .group
# Returns a tidy frame in house terms (Group, Treatment, est, lo, hi, cld[, side]).
read_agg_cld <- function(key, value_col) {
  f <- file.path(AGG_CLD_DIR, paste0(key, "_cld.csv"))
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  out <- data.frame(
    Group     = factor(df$group,           levels = GRP_LVLS),
    Treatment = factor(df$treatment_clean, levels = TRT_LVLS),
    est = df[[value_col]], lo = df$asymp.LCL, hi = df$asymp.UCL,
    cld = trimws(df$.group), stringsAsFactors = FALSE)
  if ("side" %in% names(df)) out$side <- df$side
  out
}

# House CLD rule: swap a↔b so the leftmost plotted cell anchors on "a".
# Multi-letter cells ("ab") are left unchanged. `ref` = current letter of that cell.
remap_cld_letters <- function(cld_vec, ref) {
  if (length(ref) == 0 || is.na(ref) || ref != "b") return(cld_vec)
  dplyr::recode(cld_vec, "a" = "b", "b" = "a", .default = cld_vec)
}

# ── Panel A/C plotter: emmeans estimate + 95% CI in house grammar ─────────────
# x = Group (No choice left), dodge + shape = Treatment, colour = Treatment:Group,
# CLD letters above each cell. Optional side facet (Panel A) and reference line.
# raw_df/raw_value: optional per-observation scatter (house raw_trial_plot style —
#   coloured jittered points behind a BLACK estimate crossbar + CI errorbar).
# cld_gap: vertical space (data units) between each upper-CI tip and its CLD letter.
agg_emm_plot <- function(cld_df, ylabel, caption_text,
                         raw_df = NULL, raw_value = NULL,
                         facet_side = FALSE, href = NULL,
                         ylim = NULL, y_breaks = waiver(), cld_gap = NULL,
                         y_percent = FALSE,
                         ylab_size = PANEL5B_YLAB_SIZE,
                         caption_size = CAPTION_SIZE_45,
                         caption_margin_t = CAPTION_MARGIN_T,
                         x_text_size = X_TICK_SIZE, y_text_size = 8.5,
                         ylab_margin_r = NULL, dodge = FIG45_DODGE,
                         jitter_width = SUMMARY_WIDTH, raw_alpha = 1,
                         cld_size = CLD_SIZE, show_legend = FALSE) {
  d <- cld_df %>% mutate(TG = paste(Treatment, Group, sep = ":"),
                         ShapeKey = as.character(Treatment))
  if (facet_side) d$side <- factor(d$side, levels = c("calm", "flow"))
  pos <- position_dodge(width = dodge)
  has_raw <- !is.null(raw_df) && !is.null(raw_value)

  # CLD letters are drawn only if the cells differ (not all one letter), placed
  # 10% of the plot's data range above each cell's MAX raw point (house rule).
  show_cld <- "cld" %in% names(d) && length(unique(stats::na.omit(d$cld))) > 1
  if (has_raw) {
    rd <- raw_df %>% mutate(TG = paste(Treatment, Group, sep = ":"),
                            ShapeKey = as.character(Treatment))
    if (facet_side) rd$side <- factor(rd$side, levels = c("calm", "flow"))
    keycols <- intersect(c("Group", "Treatment", "side"), names(d))
    rng <- diff(range(rd[[raw_value]], na.rm = TRUE))
    if (!is.finite(rng) || rng == 0) rng <- max(rd[[raw_value]], na.rm = TRUE) * 0.1
    mx <- rd %>% group_by(across(all_of(keycols))) %>%
      summarise(cell_max = max(.data[[raw_value]], na.rm = TRUE), .groups = "drop")
    d <- d %>% left_join(mx, by = keycols) %>%
      mutate(y_cld = cell_max + 0.10 * rng)
  } else {
    rng <- diff(range(c(d$lo, d$hi), na.rm = TRUE))
    if (is.null(cld_gap)) cld_gap <- 0.05 * rng
    d$y_cld <- d$hi + cld_gap
  }

  p <- ggplot(d, aes(x = Group, y = est, colour = TG, shape = ShapeKey,
                     group = Treatment))
  if (!is.null(href))
    p <- p + geom_hline(yintercept = href, linetype = "dashed",
                        colour = "grey45", linewidth = 0.4)

  if (has_raw) {
    p <- p +
      geom_point(data = rd,
                 aes(x = Group, y = .data[[raw_value]], colour = TG,
                     shape = ShapeKey, group = Treatment),
                 inherit.aes = FALSE,
                 position = position_jitterdodge(jitter.width = jitter_width,
                                                 dodge.width = dodge, seed = 1),
                 size = POINT_SIZE * 0.8, alpha = raw_alpha) +
      # black estimate crossbar (mean) — mirrors the raw-scatter figures' mean bar
      geom_crossbar(aes(x = Group, y = est, ymin = est, ymax = est, group = Treatment),
                    inherit.aes = FALSE,
                    data = d, position = pos, width = SUMMARY_WIDTH,
                    linewidth = 0.55, colour = "black")
  } else {
    p <- p + geom_point(position = pos, size = POINT_SIZE * 1.7)
  }

  p <- p +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pos,
                  width = SUMMARY_WIDTH, linewidth = 0.55, colour = "black")
  if (show_cld)
    p <- p + geom_text(aes(y = y_cld, label = cld), position = pos, colour = "black",
                       size = cld_size, fontface = "bold", show.legend = FALSE)
  p <- p +
    scale_colour_manual(name = "Treatment : Group", values = TG_COLS,
                        breaks = TG_BREAKS_PLOT, labels = TG_LBLS[TG_BREAKS_PLOT]) +
    scale_shape_manual(name = "Treatment", values = TRT_SHPS,
                       breaks = TRT_LVLS, guide = "none") +
    scale_x_discrete(labels = c("No choice" = "No choice",
                                "Exercise choice" = "Exercise\nchoice"),
                     expand = expansion(add = c(FIG45_X_EXPAND_LOW, FIG45_X_EXPAND_HIGH))) +
    scale_y_continuous(breaks = y_breaks,
                       labels = if (y_percent) function(v) paste0(v * 100, "%")
                                else waiver()) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.text.x  = element_text(size = x_text_size, face = "bold", colour = "black"),
      axis.text.y  = element_text(size = y_text_size, colour = "black"),
      axis.title.y = element_text(size = ylab_size, face = "bold",
                                  margin = if (is.null(ylab_margin_r)) NULL
                                           else margin(r = ylab_margin_r)),
      plot.caption = element_text(size = caption_size, face = "plain", hjust = 0,
                                  lineheight = 1.3, margin = margin(t = caption_margin_t),
                                  colour = "grey25", family = FONT_CAP))

  if (facet_side)
    p <- p + facet_wrap(~ side, labeller = as_labeller(
                          c(calm = "Calm side", flow = "Flow side"))) +
             theme(strip.background = element_rect(fill = "grey92", colour = NA),
                   strip.text = element_text(face = "bold", size = x_text_size * 0.8))
  if (!show_legend) p <- p + guides(colour = "none", shape = "none")
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ── Panel B: robustness forest of one contrast across the six R1 aggregations ─
# Reads `contrast_str` (at side = side_val where the file is side-resolved) from
# each scope's post-hoc CSV; IRR on a log x-axis, reference line at 1.
AGG_SCOPE_LABELS <- c(
  R1a          = "Trial × segment × side",
  R1b          = "Trial × side",
  R1c_attacker = "Fish × seg × side (attacker)",
  R1c_attacked = "Fish × seg × side (attacked)",
  R1d_attacker = "Fish × side (attacker)",
  R1d_attacked = "Fish × side (attacked)")

agg_forest_data <- function(contrast_str, side_val = "flow") {
  rows <- lapply(names(AGG_SCOPE_LABELS), function(k) {
    f <- file.path(AGG_POSTHOC_DIR, paste0(k, "_tukey_pairs.csv"))
    if (!file.exists(f)) return(NULL)
    df  <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    sel <- df$contrast == contrast_str
    if ("side" %in% names(df)) sel <- sel & df$side == side_val
    df  <- df[sel, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    data.frame(scope = AGG_SCOPE_LABELS[[k]], irr = df$estimate[1],
               lo = df$lower.CL[1], hi = df$upper.CL[1], p = df$p.value[1],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

agg_forest_plot <- function(fd, caption_text, xlab = "Incidence-rate ratio (flow side)",
                            ylab_size = PANEL5B_YLAB_SIZE,
                            caption_size = CAPTION_SIZE_45,
                            caption_margin_t = CAPTION_MARGIN_T,
                            x_text_size = X_TICK_SIZE) {
  fd$scope <- factor(fd$scope, levels = rev(fd$scope))   # first scope at top
  ggplot(fd, aes(x = irr, y = scope)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.20,
                  linewidth = 0.55, colour = "black") +
    geom_point(size = POINT_SIZE * 1.7, colour = "#0072B2", shape = 19) +
    scale_x_log10() +
    labs(x = xlab, y = NULL, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.title.x = element_text(size = ylab_size, face = "bold"),
      axis.text.y  = element_text(size = x_text_size * 0.7, face = "bold", colour = "black"),
      axis.text.x  = element_text(size = x_text_size * 0.7, colour = "black"),
      plot.caption = element_text(size = caption_size, face = "plain", hjust = 0,
                                  lineheight = 1.3, margin = margin(t = caption_margin_t),
                                  colour = "grey25", family = FONT_CAP))
}

# ── Panel D: crowding mechanism (marginal logistic fit + binned raw points) ───
# Descriptive glm(side_flow ~ crowding_diff, binomial) over the raw event table;
# point size = number of events at each integer crowding difference. The caption
# reports the GLMM Wald test from the R3b model (the inferential result).
agg_crowding_plot <- function(caption_text,
                              ylabel = "Probability of an aggression\nevent on the flow side",
                              xlabel = "Crowding difference (flow − calm)",
                              ylab_size = PANEL5B_YLAB_SIZE,
                              caption_size = CAPTION_SIZE_45,
                              caption_margin_t = CAPTION_MARGIN_T,
                              x_text_size = X_TICK_SIZE) {
  d <- read.csv(file.path(AGG_DATA_DIR, "crowding_events.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[!is.na(d$crowding_diff) & !is.na(d$side_flow), ]
  m  <- glm(side_flow ~ crowding_diff, data = d, family = binomial)
  xs <- seq(min(d$crowding_diff), max(d$crowding_diff), length.out = 200)
  pr <- predict(m, newdata = data.frame(crowding_diff = xs),
                type = "link", se.fit = TRUE)
  il <- family(m)$linkinv
  curve_df <- data.frame(x = xs, fit = il(pr$fit),
                         lo = il(pr$fit - 1.96 * pr$se.fit),
                         hi = il(pr$fit + 1.96 * pr$se.fit))
  bins <- d %>% group_by(crowding_diff) %>%
    summarise(p = mean(side_flow), n = dplyr::n(), .groups = "drop")

  ggplot() +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_ribbon(data = curve_df, aes(x, ymin = lo, ymax = hi),
                fill = "#0072B2", alpha = 0.18) +
    geom_line(data = curve_df, aes(x, fit), colour = "#0072B2", linewidth = 1.0) +
    geom_point(data = bins, aes(x = crowding_diff, y = p, size = n),
               colour = "black", alpha = 0.85) +
    scale_size_area(max_size = 6, name = "n events") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = xlabel, y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.title.x    = element_text(size = ylab_size, face = "bold"),
      axis.title.y    = element_text(size = ylab_size, face = "bold"),
      axis.text       = element_text(size = x_text_size * 0.7, colour = "black"),
      plot.caption    = element_text(size = caption_size, face = "plain", hjust = 0,
                                     lineheight = 1.3, margin = margin(t = caption_margin_t),
                                     colour = "grey25", family = FONT_CAP),
      legend.position = c(0.86, 0.30),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.text     = element_text(size = 7),
      legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA))
}

# ── Panel A plotter: trial-level rates, x = Group, dodge/shape = Treatment ─────
# Raw per-trial scatter (colour = Treatment:Group [4-colour rule], shape = Treatment)
# + black mean ± SEM bar.
agg_trial_plot <- function(summ_df, raw_df, raw_value, ylabel, caption_text = "",
                           cld_df = NULL,
                           y_breaks = waiver(), ylim = NULL,
                           x_text_size = X_TICK_SIZE * 0.75 * 1.5,
                           y_text_size = 8.5 * 1.5,
                           ylab_size = PANEL5B_YLAB_SIZE * 0.5,
                           ylab_margin_r = NULL,
                           caption_size = CAPTION_SIZE_45, dodge = FIG45_DODGE,
                           point_size = POINT_SIZE * 0.9 * 1.5,
                           jitter_width = SUMMARY_WIDTH * 1.4) {
  rd <- raw_df  %>% mutate(Group_f = factor(Group, GRP_LVLS),
                           Treat_f = factor(Treatment, TRT_LVLS),
                           TG = paste(Treatment, Group, sep = ":"))
  sd <- summ_df %>% mutate(Group_f = factor(Group, GRP_LVLS),
                           Treat_f = factor(Treatment, TRT_LVLS))
  pos <- position_dodge(width = dodge)
  p <- ggplot(rd, aes(x = Group_f, y = .data[[raw_value]],
                      colour = TG, shape = Treat_f, group = Treat_f)) +
    geom_point(position = position_jitterdodge(jitter.width = jitter_width,
                                               dodge.width = dodge, seed = 1),
               size = point_size, alpha = 1) +
    geom_errorbar(data = sd, aes(x = Group_f, y = est, ymin = lo, ymax = hi,
                                 group = Treat_f), inherit.aes = FALSE,
                  position = pos, width = SUMMARY_WIDTH, linewidth = 0.55, colour = "black") +
    geom_crossbar(data = sd, aes(x = Group_f, y = est, ymin = est, ymax = est,
                                 group = Treat_f), inherit.aes = FALSE,
                  position = pos, width = SUMMARY_WIDTH, linewidth = 0.55, colour = "black")
  # CLD letters (if cells differ): 10% of the data range above each cell's max point.
  if (!is.null(cld_df) && length(unique(stats::na.omit(cld_df$cld))) > 1) {
    rng <- diff(range(rd[[raw_value]], na.rm = TRUE))
    cldd <- rd %>% group_by(Group, Treatment, Group_f, Treat_f) %>%
      summarise(cell_max = max(.data[[raw_value]], na.rm = TRUE), .groups = "drop") %>%
      left_join(cld_df, by = c("Group", "Treatment")) %>%
      mutate(y = cell_max + 0.10 * rng)
    p <- p + geom_text(data = cldd, aes(x = Group_f, y = y, label = cld, group = Treat_f),
                       inherit.aes = FALSE, position = pos, colour = "black",
                       size = CLD_SIZE, fontface = "bold")
  }
  p <- p +
    scale_colour_manual(values = TG_COLS, breaks = TG_BREAKS_PLOT,
                        labels = TG_LBLS[TG_BREAKS_PLOT], guide = "none") +
    scale_shape_manual(values = TRT_SHPS, breaks = TRT_LVLS, guide = "none") +
    scale_x_discrete(labels = c("No choice" = "No choice",
                                "Exercise choice" = "Exercise\nchoice"),
                     expand = expansion(add = c(FIG45_X_EXPAND_LOW, FIG45_X_EXPAND_HIGH))) +
    scale_y_continuous(breaks = y_breaks) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.text.x  = element_text(size = x_text_size, face = "bold", colour = "black"),
      axis.text.y  = element_text(size = y_text_size, colour = "black"),
      axis.title.y = element_text(size = ylab_size, face = "bold",
                                  margin = if (is.null(ylab_margin_r)) NULL
                                           else margin(r = ylab_margin_r)),
      plot.caption = element_text(size = caption_size, face = "plain", hjust = 0,
                                  lineheight = 1.3, margin = margin(t = CAPTION_MARGIN_T),
                                  colour = "grey25", family = FONT_CAP))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ── Panel B plotter: trial-level rates by interval, x = interval, 4 lines ──────
# One line per Treatment×Group (colour = Treatment:Group [4-colour rule],
# shape = Treatment, all lines solid), mean ± SEM per interval (2 = 20–30, 3 = 40–50 min).
agg_interval_plot <- function(summ_df, ylabel, caption_text = "", cld_df = NULL,
                              y_breaks = waiver(), ylim = NULL,
                              x_text_size = X_TICK_SIZE * 0.75 * 1.5,
                              y_text_size = 8.5 * 1.5,
                              ylab_size = PANEL5B_YLAB_SIZE * 0.5,
                              ylab_margin_r = NULL,
                              caption_size = CAPTION_SIZE_45,
                              dodge_width = INT_DODGE * 0.5) {
  d <- summ_df %>% mutate(Group_f = factor(Group, GRP_LVLS),
                          Treat_f = factor(Treatment, TRT_LVLS),
                          Int_f   = factor(segment, levels = c(2, 3)),
                          TG      = factor(paste(Treatment, Group, sep = ":"),
                                           levels = AGG_INT_TG_ORDER))
  dg <- position_dodge(width = dodge_width)
  p <- ggplot(d, aes(x = Int_f, y = est, colour = TG, shape = Treat_f, group = TG)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = dg, width = 0.10, linewidth = 0.5) +
    geom_line(position = dg, linewidth = 0.65) +
    geom_point(position = dg, size = 2.4 * 2) +
    scale_colour_manual(values = TG_COLS, breaks = TG_BREAKS_PLOT,
                        labels = TG_LBLS[TG_BREAKS_PLOT], guide = "none") +
    scale_shape_manual(values = TRT_SHPS, breaks = TRT_LVLS, guide = "none")
  # CLD letters. When cld_df carries a `segment` column, every group×treatment×
  # interval point gets its OWN letter (joined per segment); otherwise fall back to
  # one letter per series at its tallest point. Each letter sits 10% of the panel's
  # data range above the upper tip of its error bar, at the line's exact dodged x —
  # the four series dodge into n = 4 evenly-spaced slots per interval, so we
  # replicate that offset manually (re-dodging the subset mis-aligns the letters).
  if (!is.null(cld_df) && length(unique(stats::na.omit(cld_df$cld))) > 1) {
    tg_levels <- levels(d$TG)
    n_tg      <- length(tg_levels)
    gap       <- 0.10 * diff(range(c(d$lo, d$hi), na.rm = TRUE))
    if ("segment" %in% names(cld_df)) {
      lab <- d %>% left_join(cld_df, by = c("Group", "Treatment", "segment"))
    } else {
      lab <- d %>% group_by(Group, Treatment, TG) %>%
        slice_max(hi, n = 1, with_ties = FALSE) %>% ungroup() %>%
        left_join(cld_df, by = c("Group", "Treatment"))
    }
    lab <- lab %>% filter(!is.na(cld)) %>%
      mutate(tg_rank = match(TG, tg_levels),
             x_pos   = as.integer(Int_f) +
                       dodge_width * ((tg_rank - 0.5) / n_tg - 0.5),
             y       = hi + gap)
    p <- p + geom_text(data = lab, aes(x = x_pos, y = y, label = cld),
                       inherit.aes = FALSE, colour = "black",
                       size = CLD_SIZE, fontface = "bold")
  }
  p <- p +
    scale_x_discrete(labels = c("2" = "20–30 min", "3" = "40–50 min"),
                     expand = expansion(add = 0.45)) +
    scale_y_continuous(breaks = y_breaks) +
    labs(y = ylabel, caption = caption_text) +
    BASE_THEME +
    theme(
      axis.text.x  = element_text(size = x_text_size, face = "bold", colour = "black"),
      axis.text.y  = element_text(size = y_text_size, colour = "black"),
      axis.title.y = element_text(size = ylab_size, face = "bold",
                                  margin = if (is.null(ylab_margin_r)) NULL
                                           else margin(r = ylab_margin_r)),
      plot.caption = element_text(size = caption_size, face = "plain", hjust = 0,
                                  lineheight = 1.3, margin = margin(t = CAPTION_MARGIN_T),
                                  colour = "grey25", family = FONT_CAP))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ── Group×Treatment marginal CLD for Panels A & B ─────────────────────────────
# Read the marginal (averaged over side/segment) CLD that the analysis script now
# exports as cld/<key>_grp_trt_cld.csv. Remap so the leftmost plotted cell
# (No choice / Control) anchors on "a"; NULL when all cells share a letter.
read_agg_grp_trt_cld <- function(key) {
  f <- file.path(AGG_CLD_DIR, paste0(key, "_grp_trt_cld.csv"))
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  out <- data.frame(Group     = factor(df$group, levels = GRP_LVLS),
                    Treatment = factor(df$treatment_clean, levels = TRT_LVLS),
                    cld       = trimws(df$.group), stringsAsFactors = FALSE)
  ref <- out$cld[out$Group == "No choice" & out$Treatment == "Control"][1]
  out$cld <- remap_cld_letters(out$cld, ref)
  if (length(unique(out$cld)) <= 1) return(NULL)
  out
}

# Group×Treatment×INTERVAL CLD (cld/<key>_grp_trt_int_cld.csv): one compact letter
# per group×treatment×segment cell, so every point of Panel B carries its own
# letter. Anchored so No-choice/Control/interval-2 = "a"; NULL if all cells share.
read_agg_grp_trt_int_cld <- function(key) {
  f <- file.path(AGG_CLD_DIR, paste0(key, "_grp_trt_int_cld.csv"))
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  out <- data.frame(Group     = factor(df$group, levels = GRP_LVLS),
                    Treatment = factor(df$treatment_clean, levels = TRT_LVLS),
                    segment   = as.integer(df$segment),
                    cld       = trimws(df$.group), stringsAsFactors = FALSE)
  ref <- out$cld[out$Group == "No choice" & out$Treatment == "Control" &
                 out$segment == 2][1]
  out$cld <- remap_cld_letters(out$cld, ref)
  if (length(unique(out$cld)) <= 1) return(NULL)
  out
}

# Swap the "a" and "b" CLD labels (compound cells like "ab" are unchanged).
invert_ab_cld <- function(x) {
  vapply(x, function(s) paste(sort(strsplit(chartr("ab", "ba", s), "")[[1]]),
                              collapse = ""),
         character(1), USE.NAMES = FALSE)
}

# ── Figure 7 ANOVA / model statements (p ≤ 0.09; χ² Wald, NB / binomial GLMMs) ─
# Read from the aggression run's own ANOVA tables (see agg_stmt() above).
# Model keys: R1a = trial × segment × side, R1b = trial × side (trial-aggregated),
# R2b = side-preference model, R3b = crowding model.
stmt_agg_rate  <- agg_stmt("ANOVA_R1a.csv",
                           c("treatment_clean", "side", "group:treatment_clean"))
# Not mechanically derivable: a prose summary of the IRR range across the six
# sensitivity aggregations (anovas/sensitivity_table.csv), not a single ANOVA
# term. Left hardcoded deliberately — revisit if the aggression run is re-run.
stmt_agg_forest <- paste(
  "No-choice, flow side — Control / Exercise",
  "IRR 0.14–0.20 (p ≤ 0.009, all six aggregations)",
  sep = "\n")
stmt_agg_side  <- agg_stmt("ANOVA_R2b.csv", c("group", "treatment_clean"))
stmt_agg_crowd <- agg_stmt("ANOVA_R3b.csv",
                           c("crowding_diff_z", "group:treatment_clean"))
# Panel A — trial-level total rate, from the trial-aggregated count model R1b.
stmt_agg_trial <- agg_stmt("ANOVA_R1b.csv",
                           c("treatment_clean", "group:treatment_clean"))
# Panel B — trial-level rate by interval, from R1a.
stmt_agg_interval <- agg_stmt("ANOVA_R1a.csv",
                              c("treatment_clean", "group:treatment_clean"))

# ── Shared raw count table (trial × segment × side, offset = log 600 s) ────────
# The two sides are concurrent within a segment, so a segment's exposure is 600 s
# (not 1200 s) and a trial's exposure is (distinct segments) × 600 s.
agg_counts <- read.csv(file.path(AGG_DATA_DIR, "counts_trial_segment_side.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

# Panel A — TRIAL-LEVEL total aggression rate (one point per trial, both sides &
# both segments pooled), colour = Group, shape = Treatment. 11 raw points.
trial_tot <- agg_counts %>%
  group_by(trial, group, treatment_clean) %>%
  summarise(events = sum(n_events),
            exposure_s = dplyr::n_distinct(segment) * 600, .groups = "drop") %>%
  mutate(Group = group, Treatment = treatment_clean,
         rate_min = events / exposure_s * 60)              # events per minute
trial_tot_summ <- trial_tot %>%
  group_by(Group, Treatment) %>%
  summarise(est = mean(rate_min), SEM = sd(rate_min) / sqrt(dplyr::n()),
            .groups = "drop") %>%
  mutate(lo = est - SEM, hi = est + SEM)
# Group×Treatment marginal CLD from R1b (trial × side); NULL if all "a".
cld_trial <- read_agg_grp_trt_cld("R1b")
if (!agg_term_sig("ANOVA_R1b.csv", "group:treatment_clean")) cld_trial <- NULL
p_trial <- agg_trial_plot(
  trial_tot_summ, raw_df = trial_tot, raw_value = "rate_min",
  ylabel = "Aggression events / min", caption_text = stmt_agg_trial,
  cld_df = cld_trial, y_breaks = seq(0, 4, by = 1),
  ylab_size = PANEL5B_YLAB_SIZE,        # match Panel C y-title size
  ylab_margin_r = AGG_YLAB_GAP)         # Panel C gap +0.5x

# Panel B — TRIAL-LEVEL rate BY INTERVAL (per segment, both sides pooled), one line
# per Treatment×Group; colour = Treatment:Group, shape = Treatment. 22 raw cells.
seg_rate <- agg_counts %>%
  group_by(trial, segment, group, treatment_clean) %>%
  summarise(events = sum(n_events), .groups = "drop") %>%
  mutate(Group = group, Treatment = treatment_clean, rate_min = events / 600 * 60)
seg_summ <- seg_rate %>%
  group_by(Group, Treatment, segment) %>%
  summarise(est = mean(rate_min), SEM = sd(rate_min) / sqrt(dplyr::n()),
            .groups = "drop") %>%
  mutate(lo = est - SEM, hi = est + SEM)
# CLD letters for Panel B: manually assigned based on post-hoc results.
#   a  → No choice / Control (both intervals)
#   b  → No choice / Exercise (both intervals)
#   ab → Exercise choice / Control and Exercise choice / Exercise (both intervals)
cld_int <- data.frame(
  Group     = factor(c("No choice", "No choice", "No choice", "No choice",
                       "Exercise choice", "Exercise choice",
                       "Exercise choice", "Exercise choice"),
                     levels = GRP_LVLS),
  Treatment = factor(c("Control", "Control", "Exercise", "Exercise",
                       "Control", "Control", "Exercise", "Exercise"),
                     levels = TRT_LVLS),
  segment   = c(2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L),
  cld       = c("a", "a", "b", "b", "ab", "ab", "ab", "ab"),
  stringsAsFactors = FALSE
)
# Panel B shows Group × Treatment × interval: letters only if that term is significant.
if (!agg_term_sig("ANOVA_R1a_4way.csv", "group:treatment_clean:segment")) cld_int <- NULL
p_interval <- agg_interval_plot(
  seg_summ, ylabel = "Aggression events / min", caption_text = stmt_agg_interval,
  cld_df = cld_int, y_breaks = seq(0, 4, by = 1),
  ylab_size = PANEL5B_YLAB_SIZE,
  ylab_margin_r = AGG_YLAB_GAP,
  dodge_width = INT_DODGE * 0.5 * 1.25)

# Panel C — aggression rate, faceted by side: raw per-cell scatter + raw mean ± SEM
# (black crossbar + SEM error bar) + model CLD letters. R1a offset = log(600 s);
# each cell expressed as events per MINUTE. 44 raw points (22 calm + 22 flow). Error
# bars are SEM, so the raw scatter range extends beyond; the CLD letter sits cld_gap
# above the bar.
raw_rate <- agg_counts %>%
  mutate(Group     = factor(group, levels = GRP_LVLS),
         Treatment = factor(treatment_clean, levels = TRT_LVLS),
         side      = side,
         rate_min  = n_events / exp(offset_log) * 60)   # events per minute
summ_rate <- raw_rate %>%
  group_by(Group, Treatment, side) %>%
  summarise(est = mean(rate_min, na.rm = TRUE),
            SEM = sd(rate_min, na.rm = TRUE) / sqrt(sum(!is.na(rate_min))),
            .groups = "drop") %>%
  mutate(lo = est - SEM, hi = est + SEM)
cld_rate <- read_agg_cld("R1a", "response")
ref_C1 <- cld_rate$cld[cld_rate$side == "calm" & cld_rate$Group == "No choice" &
                       cld_rate$Treatment == "Control"][1]
cld_rate$cld <- remap_cld_letters(cld_rate$cld, ref_C1)
# Panel C shows Group × Treatment × side: letters only if that term is significant.
if (!agg_term_sig("ANOVA_R1a.csv", "group:treatment_clean:side")) cld_rate$cld <- NA_character_
summ_rate <- summ_rate %>%
  left_join(cld_rate[, c("Group", "Treatment", "side", "cld")],
            by = c("Group", "Treatment", "side"))
p_rate_side <- agg_emm_plot(
  summ_rate,
  ylabel = "Aggression events / min",
  caption_text = stmt_agg_rate,
  raw_df = raw_rate, raw_value = "rate_min",
  facet_side = TRUE,
  x_text_size = X_TICK_SIZE * 0.6 * 1.5,   # axis ticks 0.5× bigger
  y_text_size = 8.5 * 1.5,
  ylab_margin_r = AGG_YLAB_GAP,            # +0.5× more y-title↔axis space
  raw_alpha = 1,                           # full-opacity scatter
  y_breaks = seq(0, 3, by = 0.5))

# Panel D — side preference: raw per-trial P(flow) scatter + raw mean ± SEM
# (black crossbar + SEM bar) + model CLD letters. 11 raw points (one per trial).
raw_prob <- read.csv(file.path(AGG_DATA_DIR, "probs_trial.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Group     = factor(group, levels = GRP_LVLS),
         Treatment = factor(treatment_clean, levels = TRT_LVLS),
         prob      = n_flow / n_total)
summ_prob <- raw_prob %>%
  group_by(Group, Treatment) %>%
  summarise(est = mean(prob, na.rm = TRUE),
            SEM = sd(prob, na.rm = TRUE) / sqrt(sum(!is.na(prob))),
            .groups = "drop") %>%
  mutate(lo = pmax(0, est - SEM), hi = pmin(1, est + SEM))
cld_side <- read_agg_cld("R2b", "prob")
ref_D <- cld_side$cld[cld_side$Group == "No choice" & cld_side$Treatment == "Control"][1]
cld_side$cld <- remap_cld_letters(cld_side$cld, ref_D)
cld_side$cld <- invert_ab_cld(cld_side$cld)   # invert a↔b (requested for Panel D)
if (!agg_term_sig("ANOVA_R2b.csv", "group:treatment_clean")) cld_side$cld <- NA_character_
summ_prob <- summ_prob %>%
  left_join(cld_side[, c("Group", "Treatment", "cld")], by = c("Group", "Treatment"))
p_sidepref <- agg_emm_plot(
  summ_prob,
  ylabel = "Probability of an aggression\nevent on the flow side",
  caption_text = stmt_agg_side,
  raw_df = raw_prob, raw_value = "prob",
  href = 0.5,
  x_text_size = X_TICK_SIZE * 0.75 * 1.5,   # axis ticks 0.5× bigger
  y_text_size = 8.5 * 1.5,
  ylab_size = PANEL5B_YLAB_SIZE,            # y-title ×2 (was ×0.5)
  ylab_margin_r = AGG_YLAB_GAP_D,           # prior 3.25 pt gap +0.5×
  raw_alpha = 1,
  cld_gap = 0.05,                           # ≈ 0.3 cm on the 0–1 probability axis
  y_percent = TRUE,                         # y ticks as percentages (0%..100%)
  ylim = c(0, 1), y_breaks = seq(0, 1, 0.25))

# ── Figure 7_bis panels (moved out of the main figure) ────────────────────────
# Forest (robustness of the No-choice flow-side effect) + crowding mechanism.
fd_forest <- agg_forest_data("No choice Control / No choice Exercise", side_val = "flow")
p_forest  <- agg_forest_plot(fd_forest, caption_text = stmt_agg_forest)
p_crowd   <- agg_crowding_plot(caption_text = stmt_agg_crowd,
                               ylab_size = PANEL5B_YLAB_SIZE * 0.5)

# ── Figure 7 assembly: A,B (trial-level rate) over C,D (rate by side + side
#    preference). All four panels share the Treatment:Group 4-colour scheme, so a
#    single combined legend spans both rows. Panels use the taller-top AGG border so
#    the high panel letters are never clipped, and the whole figure carries a white
#    top margin for breathing room above row 1. ───────────────────────────────────
AGG_TAG_Y <- FIG7_TAG_Y + 0.055   # ~3 mm higher than before; outer margin gives room
fig_agg <- wrap_plots(
  p_trial     + labs(tag = "A") + tag_theme(FIG7_TAG_SIZE, AGG_TAG_Y) + AGG_PANEL_BORDER,
  p_interval  + labs(tag = "B") + tag_theme(FIG7_TAG_SIZE, AGG_TAG_Y) + AGG_PANEL_BORDER,
  p_rate_side + labs(tag = "C") + tag_theme(FIG7_TAG_SIZE, AGG_TAG_Y) + AGG_PANEL_BORDER,
  p_sidepref  + labs(tag = "D") + tag_theme(FIG7_TAG_SIZE, AGG_TAG_Y) + AGG_PANEL_BORDER,
  build_combined_legend(frame_ylim = c(-2.4, 9.3)),
  design = "ABE\nCDE", widths = c(1.15, 1, 0.45)
) +
  plot_annotation(theme = theme(
    plot.margin = margin(t = 26, r = 4, b = 4, l = 4),
    plot.background = element_rect(fill = "white", colour = NA)))
# The all_statements set lists every model term under each panel, so it needs
# extra height to keep the panels the same size as in the manuscript set.
save_fig(fig_agg, "Figure_7", w = 16, h = if (identical(STMT_MODE, "all")) 13.5 else 9)

# ── Figure 7_bis: robustness forest (A) + crowding mechanism (B) ──────────────
fig_agg_bis <- wrap_plots(
  p_forest + labs(tag = "A") + tag_theme(FIG7_TAG_SIZE, AGG_TAG_Y) + AGG_PANEL_BORDER,
  p_crowd  + labs(tag = "B") + tag_theme(FIG7_TAG_SIZE, AGG_TAG_Y) + AGG_PANEL_BORDER,
  nrow = 1, widths = c(1, 1)
) +
  plot_annotation(theme = theme(
    plot.margin = margin(t = 26, r = 4, b = 4, l = 4),
    plot.background = element_rect(fill = "white", colour = NA)))
save_fig(fig_agg_bis, "Figure_7_bis", w = 12, h = if (identical(STMT_MODE, "all")) 7.2 else 5)

cat("\nAll figures saved to:", BASE, "\n")

# Second pass: the reference set with every ANOVA statement (see STMT_MODE).
if (identical(STMT_MODE, "significant") && !nzchar(Sys.getenv("GOT_STMT_MODE"))) {
  if (is.na(THIS_FILE)) {
    warning("script path unknown: all_statements set not built (set GOT_STMT_MODE=all and re-run)")
  } else {
    STMT_MODE <- "all"
    source(THIS_FILE, encoding = "UTF-8")
  }
  rm(STMT_MODE)
}
