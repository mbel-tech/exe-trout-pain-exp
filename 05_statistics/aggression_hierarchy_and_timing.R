################################################################################
# aggression_hierarchy_and_timing.R
# ─────────────────────────────────────────────────────────────────────────────
# A NEW aggression analysis, deliberately different from the count / side /
# crowding GLMMs in aggression_models.R and a fish-level aggression script (not included).
# Those scripts treat aggression as binary event COUNTS. This one instead mines
# two dimensions those analyses collapse away:
#
#   SECTION 1 — DOMINANCE HIERARCHY / SOCIAL NETWORK
#     Uses the directed attacker→attacked identity structure to reconstruct
#     per-trial (and per-segment) dominance hierarchies and asks whether group
#     (No choice vs Exercise choice) and treatment (Control vs Exercise) change
#     their steepness, linearity, and consistency.
#       • David's Scores (DS / NormDS)      de Vries, Stevens & Vervaecke 2006
#       • Steepness (+ permutation p)       de Vries et al. 2006
#       • Landau's h' (linearity)           Landau 1951 / de Vries 1995
#       • Directional Consistency (DCI)     van Hooff & Wensing 1987
#       • Triangle transitivity (Ptr, ttri) Shizuka & McDonald 2012
#       • rank consistency across intervals (Spearman ρ, segment 2 vs 3)
#       • GIVEN vs RECEIVED aggression per fish: paired test (own given vs own
#         received), net aggression (given−received) ~ group*treatment, and the
#         given↔received correlation across fish (reciprocal conflict vs despotism)
#
#   SECTION 2 — TEMPORAL DYNAMICS (point process)
#     Uses the complete event timestamps to characterise the temporal STRUCTURE
#     of aggression within each 600 s trial×segment window:
#       • inter-event intervals (IEI); CV; burstiness B (Goh & Barabási 2008)
#       • Fano factor / index of dispersion (60 s bins)
#       • latency to first attack (survival: Cox PH + Weibull AFT)
#       • within-segment escalation vs habituation (Poisson trend on 60 s bins)
#     Each dispersion metric is reported against its homogeneous-Poisson null.
#
# KEY MODELLING ASSUMPTION: for each high-intensity event the `attacker` is the
# winner / dominant individual and the `attacked` is the loser / subordinate.
# This is the standard reading of initiated aggression and the only one the data
# support (there is no separate contest-outcome field).
#
# DESIGN NOTE (drives the inference): the study is BETWEEN-TRIAL — each of the 12
# trials is a single group×treatment combo (≈3 trials per 2×2 cell). Trial-level
# hierarchy/temporal metrics are therefore small-n (n≤12, or ≤24 with segment).
# They are treated as EXPLORATORY: primary reporting is descriptive (cell means ±
# CI) backed by PERMUTATION tests; parametric models are secondary.
#
# PACKAGES: all hierarchy and burstiness indices are implemented from their
# published closed forms in base R (no EloRating/steepness/compete/igraph needed);
# only the installed `survival` is used, for latency. Input, output-folder
# convention, and I/O helpers mirror aggression_models.R. Nothing existing
# is modified — this writes to its own timestamped output root.
################################################################################


################################################################################
###  LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glmmTMB)      # secondary parametric models (Gaussian LMM, Poisson trend)
  library(car)          # Type III ANOVA
  library(survival)     # Cox PH + Weibull AFT for latency
  library(ggplot2)
  library(scales)
  library(officer)      # Word report
  library(flextable)    # pretty Word tables
})

options(contrasts = c("contr.sum", "contr.poly"))


################################################################################
###  CONFIG
################################################################################

ALPHA   <- 0.05
N_PERM  <- 5000        # permutations for trial-level metric tests
MIN_INT <- 3           # min events in a cell to compute IEI dispersion metrics
set.seed(20260705)     # reproducible permutations / bootstrap


################################################################################
###  PATHS
################################################################################

MASTER_CSV <- file.path(
  "data",
  "aggression_annotated_june2026.csv"
)

TIMESTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(
  "output/aggression_social_temporal_output",
  TIMESTAMP
)
DIR_DATASETS <- file.path(OUTPUT_DIR, "datasets")
DIR_ANOVAS   <- file.path(OUTPUT_DIR, "anovas")
DIR_POSTHOC  <- file.path(OUTPUT_DIR, "post-hoc")
DIR_REPORTS  <- file.path(OUTPUT_DIR, "reports")
DIR_FIGURES  <- file.path(OUTPUT_DIR, "figures")
for (.d in c(OUTPUT_DIR, DIR_DATASETS, DIR_ANOVAS, DIR_POSTHOC, DIR_REPORTS, DIR_FIGURES))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# -- House style / sensitivity modules (single source of truth) -------------
GOT_ROOT <- "."
source(file.path(GOT_ROOT, "05_statistics/shared/house_style.R"))
source(file.path(GOT_ROOT, "05_statistics/shared/sensitivity_checks.R"))

cat("\n=== aggression_hierarchy_and_timing.R ===\n")
cat("Input :", MASTER_CSV, "\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("  subfolders: datasets / anovas / post-hoc / reports / figures\n\n")


################################################################################
###  SMALL I/O HELPERS  (mirrors aggression_models.R)
################################################################################

write_df <- function(df, name, dir = OUTPUT_DIR) {
  path <- file.path(dir, name)
  tryCatch(data.table::fwrite(as.data.table(df), path),
           error = function(e) message("  write_df failed for ", name, ": ", e$message))
  invisible(path)
}

safe_write_lines <- function(lines, path) {
  tryCatch(writeLines(lines, path),
           error = function(e) message("  write failed ", path, ": ", e$message))
}

# p < 0.001 (fmt_p_house(), house_style.R A5) -- supersedes the old
# "<0.0001" floor.
fmt_num <- function(x, digits = 3L, na_str = "NA") {
  x <- as.numeric(x)
  result <- format(round(x, digits), trim = TRUE, scientific = FALSE)
  result[is.na(x)] <- na_str
  result
}
fmt_num_col <- function(col_name, x, digits = 3L, na_str = "NA") {
  if (grepl("p|Pr\\(", col_name, ignore.case = TRUE))
    ifelse(is.na(x), na_str, sub("^p = ", "", vapply(x, fmt_p_house, character(1))))
  else fmt_num(x, digits = digits, na_str = na_str)
}
# For the many bare scalar p-values this script prints directly in prose
# (not inside a data-frame column fmt_num_col() could dispatch on).
fmt_p_num <- function(p) ifelse(is.na(p), "NA", sub("^p = ", "", sub("^p < ", "< ", vapply(p, fmt_p_house, character(1)))))

df_to_md <- function(df, digits = 3, na = "NA") {
  df <- as.data.frame(df)
  if (is.null(df) || nrow(df) == 0L) return("_(no rows)_\n")
  cells <- as.data.frame(mapply(function(col, nm) {
    if (is.numeric(col)) fmt_num_col(nm, col, digits = digits, na_str = na)
    else { out <- as.character(col); out[is.na(out)] <- na; gsub("\\|", "\\\\|", out) }
  }, df, names(df), SIMPLIFY = FALSE), stringsAsFactors = FALSE, optional = TRUE)
  names(cells) <- names(df)
  header <- paste0("| ", paste(names(cells), collapse = " | "), " |")
  sep    <- paste0("| ", paste(rep("---", ncol(cells)), collapse = " | "), " |")
  body   <- vapply(seq_len(nrow(cells)),
                   function(i) paste0("| ", paste(unlist(cells[i, ]), collapse = " | "), " |"),
                   character(1))
  paste0(paste(c(header, sep, body), collapse = "\n"), "\n")
}

.df_to_flex <- function(df, font_size = 8L) {
  df <- as.data.frame(df)
  df[] <- mapply(function(col, nm) if (is.numeric(col)) fmt_num_col(nm, col) else as.character(col),
                df, names(df), SIMPLIFY = FALSE)
  ft <- flextable::flextable(df)
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::bg(ft, part = "header", bg = "#2E75B6")
  ft <- flextable::color(ft, part = "header", color = "white")
  ft <- flextable::fontsize(ft, part = "all", size = font_size)
  ft <- flextable::font(ft, part = "all", fontname = "Calibri")
  ft <- flextable::align(ft, part = "header", align = "center")
  ft <- flextable::set_table_properties(ft, width = 1, layout = "autofit")
  ft
}


################################################################################
###  DATA IMPORT AND VALIDATION  (mirrors aggression_models.R)
################################################################################

cat("--- Loading and validating master CSV ---\n")
if (!file.exists(MASTER_CSV)) stop("Master CSV not found: ", MASTER_CSV)

dt <- fread(MASTER_CSV, encoding = "UTF-8")
cat(sprintf("  Loaded %d rows × %d columns\n", nrow(dt), ncol(dt)))
stopifnot("Expected 442 rows" = nrow(dt) == 442, "Expected 54 columns" = ncol(dt) == 54)

dt[, trial   := as.integer(trial)]
dt[, segment := factor(segment, levels = c(2L, 3L))]
dt[, group   := factor(group, levels = c("No choice", "Exercise choice"))]
dt[, treatment_clean := factor(treatment_clean, levels = c("Control", "Exercise"))]
dt[, side    := factor(side, levels = c("calm", "flow"))]
dt[, attacker_match_accepted := as.logical(attacker_match_accepted)]
dt[, attacked_match_accepted := as.logical(attacked_match_accepted)]
dt[, date := as.Date(date)]
dt[, time_absolute_sec       := as.numeric(time_absolute_sec)]
dt[, time_within_segment_sec := as.numeric(time_within_segment_sec)]

stopifnot(
  "segment values"       = all(dt$segment %in% c(2L, 3L)),
  "side values"          = all(dt$side %in% c("calm", "flow")),
  "no missing trial"     = all(!is.na(dt$trial)),
  "valid_fish always 4"  = all(dt$valid_fish_in_window == 4L),
  "all high_intensity"   = all(dt$interaction_type == "high_intensity"),
  "timestamps present"   = all(!is.na(dt$time_within_segment_sec))
)
cat(sprintf("  Validation passed. Trials: %s\n",
            paste(sort(unique(dt$trial)), collapse = ", ")))


################################################################################
###  DATA DERIVATION
################################################################################

cat("\n--- Building derived datasets ---\n")

trial_design <- unique(dt[, .(trial, group, treatment_clean, date)])
stopifnot(nrow(trial_design) == uniqueN(dt$trial))

# events_full: matched attacker/attacked identities (as in aggression_models.R).
events_full <- copy(dt)
events_full[, fish_attacker := attacker_fish_id_matched_trial]
events_full[, fish_attacked := attacked_fish_id_matched_trial]
events_full[, match_accepted_both := attacker_match_accepted & attacked_match_accepted]

# DYADIC set for hierarchy: both fish identified AND both matches accepted, and a
# real contest (attacker ≠ attacked). This is the winner/loser event stream.
events_dyad <- events_full[
  !is.na(fish_attacker) & !is.na(fish_attacked) &
  match_accepted_both == TRUE & fish_attacker != fish_attacked]
n_excluded <- nrow(events_full) - nrow(events_dyad)
cat(sprintf("  Dyadic events (both fish matched & accepted): %d/%d  (%d excluded)\n",
            nrow(events_dyad), nrow(events_full), n_excluded))
write_df(events_dyad[, .(trial, segment, group, treatment_clean,
                         time_absolute_sec, time_within_segment_sec,
                         fish_attacker, fish_attacked, side)],
         "events_dyadic.csv", DIR_DATASETS)


################################################################################
###  SECTION 1 HELPERS — dominance-hierarchy indices (closed form, base R)
################################################################################

# Build an n×n directed win matrix M where M[i,j] = # times i attacked (beat) j.
# `ids` fixes the row/column order (individuals present in the trial).
build_win_matrix <- function(events, ids) {
  M <- matrix(0L, length(ids), length(ids), dimnames = list(ids, ids))
  if (nrow(events) == 0L) return(M)
  for (k in seq_len(nrow(events))) {
    a <- as.character(events$fish_attacker[k]); b <- as.character(events$fish_attacked[k])
    if (a %in% ids && b %in% ids) M[a, b] <- M[a, b] + 1L
  }
  M
}

# David's Scores using the chance-corrected dyadic dominance index Dij
# (de Vries, Stevens & Vervaecke 2006). Returns DS and NormDS per individual.
#   Pij = Mij / (Mij + Mji);  Dij = Pij − (Pij − 0.5)/(nij + 1);  (Dij = 0 if nij = 0)
#   w  = ΣDij ; w2 = Σ(Dij·wj) ; l = ΣDji ; l2 = Σ(Dji·lj)
#   DS = w + w2 − l − l2 ;  NormDS = (DS + n(n−1)/2) / n   ∈ [0, n−1]
davids_scores <- function(M) {
  n <- nrow(M); ids <- rownames(M)
  D <- matrix(0, n, n, dimnames = dimnames(M))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
    nij <- M[i, j] + M[j, i]
    if (nij > 0) {
      Pij <- M[i, j] / nij
      D[i, j] <- Pij - (Pij - 0.5) / (nij + 1)
    }
  }
  w  <- rowSums(D)
  l  <- colSums(D)
  w2 <- as.numeric(D %*% w)
  l2 <- as.numeric(t(D) %*% l)
  DS <- w + w2 - l - l2
  NormDS <- (DS + n * (n - 1) / 2) / n
  data.frame(fish = ids, DS = DS, NormDS = NormDS,
             attacks_given = rowSums(M), attacks_received = colSums(M),
             row.names = NULL, stringsAsFactors = FALSE)
}

# Binary dominance matrix: B[i,j] = 1 if i dominates j (Mij > Mji), 0.5 if the dyad
# is tied OR unknown (nij = 0), 0 otherwise. The 0.5-for-unknown convention makes
# the linearity index the "improved" h' of de Vries (1995).
binary_dominance <- function(M) {
  n <- nrow(M); B <- matrix(0, n, n, dimnames = dimnames(M))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
    if      (M[i, j] >  M[j, i]) B[i, j] <- 1
    else if (M[i, j] == M[j, i]) B[i, j] <- 0.5   # tie or unknown
  }
  B
}

# Landau's linearity h' from the binary dominance matrix (row sums Va).
landau_hprime <- function(M) {
  n <- nrow(M); if (n < 2) return(NA_real_)
  Va <- rowSums(binary_dominance(M))
  (12 / (n^3 - n)) * sum((Va - (n - 1) / 2)^2)
}

# Directional Consistency Index: Σ(H−L)/Σ(H+L) over unordered dyads.
# H = attacks in the more-frequent direction, L = in the less-frequent. 0..1.
dci <- function(M) {
  n <- nrow(M); H <- 0; L <- 0
  for (i in seq_len(n)) for (j in seq_len(n)) if (i < j) {
    a <- M[i, j]; b <- M[j, i]
    H <- H + max(a, b); L <- L + min(a, b)
  }
  if ((H + L) == 0) return(NA_real_)
  (H - L) / (H + L)
}

# Triangle transitivity. Over all triads whose three dyads are all "known"
# (a strict winner each), Ptr = proportion transitive; ttri = (Ptr − 0.75)/0.25
# (Shizuka & McDonald 2012). Perfect linear order → Ptr = 1, ttri = 1.
triangle_transitivity <- function(M) {
  n <- nrow(M); if (n < 3) return(c(Ptr = NA_real_, ttri = NA_real_, n_triads = 0))
  dominates <- function(i, j) M[i, j] > M[j, i]
  known     <- function(i, j) M[i, j] != M[j, i]
  tri <- utils::combn(n, 3); n_ok <- 0L; n_trans <- 0L
  for (t in seq_len(ncol(tri))) {
    ijk <- tri[, t]; i <- ijk[1]; j <- ijk[2]; k <- ijk[3]
    if (known(i, j) && known(j, k) && known(i, k)) {
      n_ok <- n_ok + 1L
      score <- c(sum(dominates(i, j), dominates(i, k)),
                 sum(dominates(j, i), dominates(j, k)),
                 sum(dominates(k, i), dominates(k, j)))
      if (setequal(sort(score), c(0, 1, 2))) n_trans <- n_trans + 1L  # transitive
    }
  }
  if (n_ok == 0L) return(c(Ptr = NA_real_, ttri = NA_real_, n_triads = 0))
  Ptr <- n_trans / n_ok
  c(Ptr = Ptr, ttri = (Ptr - 0.75) / 0.25, n_triads = n_ok)
}

# Steepness (de Vries, Stevens & Vervaecke 2006): the absolute slope of the linear
# regression of NormDS on rank (rank 1 = highest NormDS). This is the canonical,
# Dij-based (chance-corrected) steepness; its theoretical maximum approaches 1 as
# interactions per dyad → ∞ but is < 1 under finite sampling (the correction shrinks
# extreme win proportions toward 0.5). Optional permutation p reshuffles the outcome
# of each interaction within its dyad (fair-coin null) and counts steeper-or-equal.
steepness <- function(M, n_perm = 0) {
  n <- nrow(M)
  nds <- sort(davids_scores(M)$NormDS, decreasing = TRUE)
  rank_x <- seq_len(n)
  obs <- tryCatch(abs(as.numeric(coef(lm(nds ~ rank_x))[2])), error = function(e) NA_real_)
  p_perm <- NA_real_
  if (n_perm > 0 && !is.na(obs)) {
    cnt <- 0L; valid <- 0L
    for (p in seq_len(n_perm)) {
      Mp <- M
      for (i in seq_len(n)) for (j in seq_len(n)) if (i < j) {
        tot <- M[i, j] + M[j, i]
        if (tot > 0) { fwd <- rbinom(1, tot, 0.5); Mp[i, j] <- fwd; Mp[j, i] <- tot - fwd }
      }
      ndp <- sort(davids_scores(Mp)$NormDS, decreasing = TRUE)
      sp  <- tryCatch(abs(coef(lm(ndp ~ rank_x))[2]), error = function(e) NA_real_)
      if (!is.na(sp)) { valid <- valid + 1L; if (sp >= obs) cnt <- cnt + 1L }
    }
    p_perm <- if (valid > 0) (cnt + 1) / (valid + 1) else NA_real_
  }
  c(steepness = obs, p_perm = p_perm)
}

# One-shot: all hierarchy metrics for a win matrix.
hierarchy_metrics <- function(M, n_perm = 0) {
  st <- steepness(M, n_perm = n_perm)
  tt <- triangle_transitivity(M)
  data.frame(
    n_ind     = nrow(M),
    n_events  = sum(M),
    steepness = as.numeric(st["steepness"]),
    steep_p   = as.numeric(st["p_perm"]),
    landau_hprime = landau_hprime(M),
    DCI       = dci(M),
    Ptr       = as.numeric(tt["Ptr"]),
    ttri      = as.numeric(tt["ttri"]),
    n_triads  = as.numeric(tt["n_triads"]),
    row.names = NULL
  )
}


################################################################################
###  SECTION 1 SELF-TESTS  (fail fast if a closed form is wrong)
################################################################################

cat("\n--- Hierarchy metric self-tests ---\n")
# Perfect linear 4-fish order A>B>C>D: each higher beats each lower 3 times, 0 back.
.Lin <- matrix(0L, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
for (i in 1:4) for (j in 1:4) if (i < j) .Lin[i, j] <- 3L
.ds  <- davids_scores(.Lin)
.hm  <- hierarchy_metrics(.Lin, n_perm = 0)
stopifnot(
  "DS sum to 0"      = abs(sum(.ds$DS)) < 1e-8,
  "NormDS in [0,n-1]"= all(.ds$NormDS >= -1e-8 & .ds$NormDS <= 3 + 1e-8),
  "h' == 1 (linear)" = abs(.hm$landau_hprime - 1) < 1e-8,
  "DCI == 1 (linear)"= abs(.hm$DCI - 1) < 1e-8,
  "Ptr == 1 (linear)"= abs(.hm$Ptr - 1) < 1e-8,
  "ttri == 1 (linear)"= abs(.hm$ttri - 1) < 1e-8
)
cat("  PASS: linear order → h'=1, DCI=1, Ptr=1, ttri=1, ΣDS=0.\n")
cat(sprintf("  (linear-order steepness w/ 3 interactions per dyad = %.3f; →1 as nij→∞)\n",
            .hm$steepness))


################################################################################
###  SECTION 1 — DOMINANCE HIERARCHY / SOCIAL NETWORK  (analysis)
################################################################################

cat("\n\n########## SECTION 1 — DOMINANCE HIERARCHY ##########\n")

trial_ids <- sort(unique(events_dyad$trial))
# Fish present per trial = union of attacker/attacked matched IDs in the trial.
fish_in_trial <- lapply(trial_ids, function(tr) {
  e <- events_dyad[trial == tr]
  sort(unique(c(as.character(e$fish_attacker), as.character(e$fish_attacked))))
})
names(fish_in_trial) <- as.character(trial_ids)

meta_for_trial <- function(tr) trial_design[trial == tr][1]

# ── 1a. Per-trial hierarchy metrics (segments pooled) ────────────────────────
hier_trial <- rbindlist(lapply(trial_ids, function(tr) {
  ids <- fish_in_trial[[as.character(tr)]]
  M   <- build_win_matrix(events_dyad[trial == tr], ids)
  hm  <- hierarchy_metrics(M, n_perm = N_PERM)
  md  <- meta_for_trial(tr)
  # verify win-matrix totals equal per-fish given/received counts
  stopifnot(sum(M) == nrow(events_dyad[trial == tr]))
  # Keep the attacker (row) fish id as an explicit column; columns = attacked fish.
  M_out <- cbind(attacker_fish = rownames(M), as.data.frame(M))
  write_df(M_out, sprintf("win_matrix_trial%02d.csv", tr), DIR_DATASETS)
  data.table(trial = tr, group = md$group, treatment_clean = md$treatment_clean, hm)
}), fill = TRUE)
write_df(hier_trial, "hierarchy_metrics_trial.csv", DIR_ANOVAS)
cat("  Per-trial hierarchy metrics (n = ", nrow(hier_trial), " trials):\n", sep = "")
print(hier_trial[, .(trial, group, treatment_clean, n_events,
                     steepness, steep_p, landau_hprime, DCI, ttri)])

# ── 1b. Per-trial×segment hierarchy metrics ──────────────────────────────────
hier_seg <- rbindlist(lapply(trial_ids, function(tr) {
  rbindlist(lapply(c(2L, 3L), function(sg) {
    e   <- events_dyad[trial == tr & segment == sg]
    ids <- fish_in_trial[[as.character(tr)]]
    if (nrow(e) == 0L) return(NULL)
    M   <- build_win_matrix(e, ids)
    hm  <- hierarchy_metrics(M, n_perm = 0)   # permutation skipped at segment level (sparse)
    md  <- meta_for_trial(tr)
    data.table(trial = tr, segment = sg, group = md$group,
               treatment_clean = md$treatment_clean, hm)
  }))
}), fill = TRUE)
write_df(hier_seg, "hierarchy_metrics_trial_segment.csv", DIR_ANOVAS)

# ── 1c. Fish-level David's scores (segments pooled) ──────────────────────────
david_fish <- rbindlist(lapply(trial_ids, function(tr) {
  ids <- fish_in_trial[[as.character(tr)]]
  M   <- build_win_matrix(events_dyad[trial == tr], ids)
  ds  <- davids_scores(M)
  md  <- meta_for_trial(tr)
  data.table(trial = tr, group = md$group, treatment_clean = md$treatment_clean, ds)
}), fill = TRUE)
# (written after 1h, once net_aggression = given − received is added)

# ── 1d. Rank consistency across intervals (Spearman ρ of David's rank) ────────
rank_consist <- rbindlist(lapply(trial_ids, function(tr) {
  ids <- fish_in_trial[[as.character(tr)]]
  e2  <- events_dyad[trial == tr & segment == 2L]
  e3  <- events_dyad[trial == tr & segment == 3L]
  if (nrow(e2) < 2L || nrow(e3) < 2L) return(NULL)
  ds2 <- davids_scores(build_win_matrix(e2, ids))
  ds3 <- davids_scores(build_win_matrix(e3, ids))
  m   <- merge(ds2[, c("fish", "DS")], ds3[, c("fish", "DS")], by = "fish",
               suffixes = c("_s2", "_s3"))
  if (nrow(m) < 3L) return(NULL)
  rho <- suppressWarnings(stats::cor(m$DS_s2, m$DS_s3, method = "spearman"))
  md  <- meta_for_trial(tr)
  data.table(trial = tr, group = md$group, treatment_clean = md$treatment_clean,
             n_fish = nrow(m), spearman_rho = rho)
}), fill = TRUE)
if (nrow(rank_consist)) {
  write_df(rank_consist, "rank_consistency.csv", DIR_ANOVAS)
  cat(sprintf("\n  Rank consistency (segment 2 vs 3): mean ρ = %.3f over %d trials\n",
              mean(rank_consist$spearman_rho, na.rm = TRUE), nrow(rank_consist)))
}

# ── 1e. Descriptives + permutation tests of group / treatment main effects ───
# Between-trial 2×2: shuffle the group (or treatment) label across the 12 trials.
perm_main_effect <- function(y, fac, n_perm = N_PERM) {
  ok <- !is.na(y); y <- y[ok]; fac <- droplevels(as.factor(fac[ok]))
  if (nlevels(fac) != 2 || length(y) < 4) return(c(diff = NA_real_, p = NA_real_))
  lv  <- levels(fac)
  obs <- mean(y[fac == lv[2]]) - mean(y[fac == lv[1]])
  cnt <- 0L
  for (p in seq_len(n_perm)) {
    fp <- sample(fac)
    d  <- mean(y[fp == lv[2]]) - mean(y[fp == lv[1]])
    if (abs(d) >= abs(obs) - 1e-12) cnt <- cnt + 1L
  }
  c(diff = obs, p = (cnt + 1) / (n_perm + 1))
}

metric_cols <- c("steepness", "landau_hprime", "DCI", "ttri")
perm_tbl <- rbindlist(lapply(metric_cols, function(mc) {
  y <- hier_trial[[mc]]
  pe_g <- perm_main_effect(y, hier_trial$group)
  pe_t <- perm_main_effect(y, hier_trial$treatment_clean)
  data.table(
    metric              = mc,
    mean_overall        = mean(y, na.rm = TRUE),
    mean_NoChoice       = mean(y[hier_trial$group == "No choice"], na.rm = TRUE),
    mean_ExChoice       = mean(y[hier_trial$group == "Exercise choice"], na.rm = TRUE),
    group_diff          = pe_g["diff"], group_p_perm = pe_g["p"],
    mean_Control        = mean(y[hier_trial$treatment_clean == "Control"], na.rm = TRUE),
    mean_Exercise       = mean(y[hier_trial$treatment_clean == "Exercise"], na.rm = TRUE),
    treatment_diff      = pe_t["diff"], treatment_p_perm = pe_t["p"])
}), fill = TRUE)
write_df(perm_tbl, "hierarchy_permutation_tests.csv", DIR_POSTHOC)
cat("\n  Permutation tests (group & treatment main effects on trial-level metrics):\n")
print(perm_tbl[, .(metric, group_diff, group_p_perm, treatment_diff, treatment_p_perm)])

# ── 1f. Secondary parametric models (LM per metric; individual NormDS LMM) ────
lm_tbl <- rbindlist(lapply(metric_cols, function(mc) {
  d <- as.data.frame(hier_trial); d$y <- d[[mc]]
  m <- tryCatch(lm(y ~ group * treatment_clean, data = d), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  a <- tryCatch(car::Anova(m, type = 3), error = function(e) NULL)
  if (is.null(a)) return(NULL)
  ad <- as.data.frame(a); ad$term <- rownames(ad); ad$metric <- mc
  as.data.table(ad)
}), fill = TRUE)
if (nrow(lm_tbl)) write_df(lm_tbl, "hierarchy_lm_anova.csv", DIR_ANOVAS)

# Individual NormDS ~ group*treatment + (1|trial): tests whether individual
# dominance strength shifts with group/treatment (≈48 fish).
m_nds <- tryCatch(
  glmmTMB(NormDS ~ group * treatment_clean + (1 | trial), data = david_fish),
  error = function(e) { message("  NormDS LMM failed: ", e$message); NULL })
if (!is.null(m_nds)) {
  a_nds <- tryCatch(car::Anova(m_nds, type = 3), error = function(e) NULL)
  if (!is.null(a_nds)) {
    a_nds_df <- cbind(term = rownames(a_nds), as.data.frame(a_nds))
    write_df(a_nds_df, "ANOVA_individual_NormDS.csv", DIR_ANOVAS)
  }
}

# ── 1g. Sociograms (one panel per trial; node size ~ NormDS, arrows ~ attacks) ─
soc_nodes <- rbindlist(lapply(trial_ids, function(tr) {
  ids <- fish_in_trial[[as.character(tr)]]; k <- length(ids)
  ds  <- davids_scores(build_win_matrix(events_dyad[trial == tr], ids))
  ang <- seq(0, 2 * pi, length.out = k + 1)[seq_len(k)]
  data.table(trial = tr, fish = ids, x = cos(ang), y = sin(ang),
             NormDS = ds$NormDS[match(ids, ds$fish)])
}), fill = TRUE)
soc_edges <- rbindlist(lapply(trial_ids, function(tr) {
  ids <- fish_in_trial[[as.character(tr)]]
  M   <- build_win_matrix(events_dyad[trial == tr], ids)
  nd  <- soc_nodes[trial == tr]
  ee  <- list()
  for (i in ids) for (j in ids) if (i != j && M[i, j] > 0) {
    ee[[length(ee) + 1]] <- data.table(
      trial = tr, from = i, to = j, w = M[i, j],
      x = nd$x[nd$fish == i], y = nd$y[nd$fish == i],
      xend = nd$x[nd$fish == j], yend = nd$y[nd$fish == j])
  }
  if (length(ee)) rbindlist(ee) else NULL
}), fill = TRUE)
p_soc <- ggplot() +
  geom_segment(data = soc_edges,
               aes(x = x, y = y, xend = 0.82 * xend + 0.18 * x,
                   yend = 0.82 * yend + 0.18 * y, linewidth = w),
               arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
               colour = "grey45", alpha = 0.75) +
  geom_point(data = soc_nodes, aes(x, y, size = NormDS), colour = "#2E75B6") +
  geom_text(data = soc_nodes, aes(x, y, label = round(NormDS, 1)),
            colour = "white", size = 2) +
  scale_linewidth(range = c(0.3, 2.2), name = "attacks") +
  scale_size(range = c(4, 10), name = "NormDS") +
  coord_equal() + facet_wrap(~ trial, ncol = 4) +
  labs(title = "Per-trial dominance sociograms (arrow = attacker→attacked, node = NormDS)") +
  theme_void() + theme(strip.text = element_text(face = "bold"))
ggsave(file.path(DIR_FIGURES, "sociograms_by_trial.png"), p_soc,
       width = 12, height = 9, dpi = 200)
cat("  wrote figures/sociograms_by_trial.png\n")

# ── 1h. Given vs received aggression (per fish) ──────────────────────────────
# Three complementary tests on the same david_fish table (attacks_given,
# attacks_received already computed per fish per trial):
#   (i)   PAIRED test — is a fish's OWN given count different from its OWN
#         received count (Wilcoxon signed-rank, + paired t as a parametric check)?
#   (ii)  NET AGGRESSION (given − received) ~ group*treatment_clean + (1|trial) —
#         does the given:received BALANCE shift with group/treatment?
#   (iii) CORRELATION of given vs received across fish — positive ⇒ mutual/
#         reciprocal conflict (dominant attackers also get attacked back more);
#         near-zero/negative ⇒ one-sided despotism (attackers ≠ targets).
cat("\n  --- Given vs received aggression (per fish) ---\n")

# (i) Paired test.
gr_wilcox <- tryCatch(
  wilcox.test(david_fish$attacks_given, david_fish$attacks_received, paired = TRUE),
  error = function(e) NULL)
gr_ttest  <- tryCatch(
  t.test(david_fish$attacks_given, david_fish$attacks_received, paired = TRUE),
  error = function(e) NULL)
given_received_paired <- data.frame(
  n_fish            = nrow(david_fish),
  mean_given        = mean(david_fish$attacks_given),
  mean_received     = mean(david_fish$attacks_received),
  mean_diff         = mean(david_fish$attacks_given - david_fish$attacks_received),
  wilcox_V          = if (!is.null(gr_wilcox)) as.numeric(gr_wilcox$statistic) else NA_real_,
  wilcox_p          = if (!is.null(gr_wilcox)) gr_wilcox$p.value else NA_real_,
  paired_t          = if (!is.null(gr_ttest))  as.numeric(gr_ttest$statistic)  else NA_real_,
  paired_t_p        = if (!is.null(gr_ttest))  gr_ttest$p.value               else NA_real_)
write_df(given_received_paired, "given_vs_received_paired_test.csv", DIR_POSTHOC)
cat(sprintf("    Paired test (n=%d fish): mean given=%.2f, mean received=%.2f, diff=%.2f | Wilcoxon p=%s, paired-t p=%s\n",
            given_received_paired$n_fish, given_received_paired$mean_given,
            given_received_paired$mean_received, given_received_paired$mean_diff,
            fmt_p_num(given_received_paired$wilcox_p), fmt_p_num(given_received_paired$paired_t_p)))

# (ii) Net aggression (given − received) ~ group * treatment_clean + (1|trial).
david_fish[, net_aggression := attacks_given - attacks_received]
a_net_df <- NULL
m_net <- tryCatch(
  glmmTMB(net_aggression ~ group * treatment_clean + (1 | trial), data = david_fish),
  error = function(e) { message("    net-aggression LMM failed: ", e$message); NULL })
if (!is.null(m_net)) {
  a_net <- tryCatch(car::Anova(m_net, type = 3), error = function(e) NULL)
  if (!is.null(a_net)) {
    a_net_df <- cbind(term = rownames(a_net), as.data.frame(a_net))
    write_df(a_net_df, "ANOVA_net_aggression.csv", DIR_ANOVAS)
    cat("    Net aggression (given-received) ~ group*treatment ANOVA written (ANOVA_net_aggression.csv).\n")
  }
}
write_df(david_fish, "david_scores_fish.csv", DIR_ANOVAS)   # re-write with net_aggression added

# (iii) Given vs received correlation across fish (mutual conflict vs despotism).
gr_cor <- suppressWarnings(cor.test(david_fish$attacks_given, david_fish$attacks_received,
                                    method = "spearman"))
given_received_cor <- data.frame(
  n_fish = nrow(david_fish),
  spearman_rho = unname(gr_cor$estimate), p_value = gr_cor$p.value)
write_df(given_received_cor, "given_vs_received_correlation.csv", DIR_POSTHOC)
cat(sprintf("    Given↔received correlation: Spearman ρ = %.3f (p = %s) — %s\n",
            given_received_cor$spearman_rho, fmt_p_num(given_received_cor$p_value),
            if (!is.na(given_received_cor$p_value) && given_received_cor$p_value < ALPHA &&
                given_received_cor$spearman_rho > 0)
              "reciprocal/mutual conflict (dominant attackers also targeted more)"
            else "no evidence of reciprocity (given and received are not linked)"))

p_gr <- ggplot(david_fish, aes(attacks_given, attacks_received, colour = treatment_clean,
                               shape = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(data = david_fish, mapping = aes(attacks_given, attacks_received),
              method = "lm", se = TRUE, colour = "black", linewidth = 0.6,
              inherit.aes = FALSE) +
  scale_colour_manual(values = c("Control" = "#2E75B6", "Exercise" = "#D55E00")) +
  labs(x = "Attacks given (per fish, per trial)", y = "Attacks received (per fish, per trial)",
       colour = "Treatment", shape = "Group",
       title = "Given vs received aggression per fish",
       subtitle = sprintf("dashed line = equal given/received; Spearman ρ = %.2f, p = %s",
                          given_received_cor$spearman_rho, fmt_p_num(given_received_cor$p_value))) +
  BASE_THEME
ggsave(file.path(DIR_FIGURES, "given_vs_received.png"), p_gr, width = 7, height = 5.5, dpi = 200)
cat("    wrote figures/given_vs_received.png\n")


################################################################################
###  SECTION 2 — TEMPORAL DYNAMICS (point process)
################################################################################

cat("\n\n########## SECTION 2 — TEMPORAL DYNAMICS ##########\n")

SEG_LEN  <- 600      # segment window length (s)
BIN_W    <- 60       # bin width for Fano / trend (s) → 10 bins per segment
N_BINS   <- SEG_LEN / BIN_W

# Temporal structure of ALL aggression events (identity not required here).
ev_time <- events_full[, .(trial, segment, group, treatment_clean,
                           t = time_within_segment_sec)]
setorder(ev_time, trial, segment, t)

# Burstiness / dispersion metrics for one ordered event-time vector.
#   CV   = sd(IEI)/mean(IEI)            (1 under homogeneous Poisson)
#   B    = (sd-mean)/(sd+mean)          (0 under Poisson; →1 bursty; <0 regular)
#   Fano = var/mean of counts in bins   (1 under Poisson)
temporal_metrics <- function(t, seg_len = SEG_LEN, bin_w = BIN_W) {
  t <- sort(t[!is.na(t)]); n <- length(t)
  iei <- if (n >= 2) diff(t) else numeric(0)
  cv  <- if (length(iei) >= 2 && mean(iei) > 0) sd(iei) / mean(iei) else NA_real_
  B   <- if (length(iei) >= 2 && (sd(iei) + mean(iei)) > 0)
           (sd(iei) - mean(iei)) / (sd(iei) + mean(iei)) else NA_real_
  bins <- cut(t, breaks = seq(0, seg_len, by = bin_w), include.lowest = TRUE)
  cnts <- as.integer(table(factor(bins, levels = levels(cut(0, seq(0, seg_len, by = bin_w))))))
  fano <- if (mean(cnts) > 0) var(cnts) / mean(cnts) else NA_real_
  data.frame(n_events = n, mean_iei = if (length(iei)) mean(iei) else NA_real_,
             cv_iei = cv, burstiness_B = B, fano = fano)
}

# ── 2a. Per trial×segment dispersion metrics ─────────────────────────────────
disp <- rbindlist(lapply(trial_ids, function(tr) {
  rbindlist(lapply(c(2L, 3L), function(sg) {
    e <- ev_time[trial == tr & segment == sg]
    md <- meta_for_trial(tr)
    tm <- temporal_metrics(e$t)
    data.table(trial = tr, segment = sg, group = md$group,
               treatment_clean = md$treatment_clean, tm)
  }))
}), fill = TRUE)
write_df(disp, "temporal_dispersion_metrics.csv", DIR_ANOVAS)
cat(sprintf("  Dispersion metrics over %d trial×segment cells (Poisson null: CV=1, B=0, Fano=1):\n",
            nrow(disp)))
cat(sprintf("    mean CV = %.2f | mean B = %.2f | mean Fano = %.2f  (cells with ≥%d events used)\n",
            mean(disp$cv_iei, na.rm = TRUE), mean(disp$burstiness_B, na.rm = TRUE),
            mean(disp$fano, na.rm = TRUE), MIN_INT))

# One-sample tests vs the Poisson null (are events clustered overall?).
disp_ok <- disp[n_events >= MIN_INT]
null_tests <- data.table(
  metric = c("cv_iei", "burstiness_B", "fano"),
  null_value = c(1, 0, 1),
  mean_obs = c(mean(disp_ok$cv_iei, na.rm = TRUE),
               mean(disp_ok$burstiness_B, na.rm = TRUE),
               mean(disp_ok$fano, na.rm = TRUE)),
  wilcox_p = c(
    tryCatch(wilcox.test(disp_ok$cv_iei,       mu = 1)$p.value, error = function(e) NA_real_),
    tryCatch(wilcox.test(disp_ok$burstiness_B, mu = 0)$p.value, error = function(e) NA_real_),
    tryCatch(wilcox.test(disp_ok$fano,         mu = 1)$p.value, error = function(e) NA_real_)))
write_df(null_tests, "temporal_vs_poisson_null.csv", DIR_POSTHOC)
cat("  Deviation from homogeneous-Poisson null (Wilcoxon signed-rank):\n")
print(null_tests)

# Store IEIs for the record.
iei_rows <- rbindlist(lapply(trial_ids, function(tr) {
  rbindlist(lapply(c(2L, 3L), function(sg) {
    e <- ev_time[trial == tr & segment == sg]; t <- sort(e$t)
    if (length(t) < 2) return(NULL)
    md <- meta_for_trial(tr)
    data.table(trial = tr, segment = sg, group = md$group,
               treatment_clean = md$treatment_clean, iei = diff(t))
  }))
}), fill = TRUE)
if (nrow(iei_rows)) write_df(iei_rows, "iei_by_cell.csv", DIR_DATASETS)

# ── 2b. Group×treatment models on dispersion metrics (secondary, small-n) ─────
disp_anova <- rbindlist(lapply(c("cv_iei", "burstiness_B", "fano"), function(mc) {
  d <- as.data.frame(disp_ok); d$y <- d[[mc]]
  d <- d[!is.na(d$y), ]
  m <- tryCatch(glmmTMB(y ~ group * treatment_clean * segment + (1 | trial), data = d),
                error = function(e) NULL)
  if (is.null(m) || !is.finite(AIC(m)))
    m <- tryCatch(glmmTMB(y ~ group * treatment_clean + segment + (1 | trial), data = d),
                  error = function(e) NULL)
  if (is.null(m)) return(NULL)
  a <- tryCatch(car::Anova(m, type = 3), error = function(e) NULL)
  if (is.null(a)) return(NULL)
  ad <- cbind(term = rownames(a), as.data.frame(a)); ad$metric <- mc
  as.data.table(ad)
}), fill = TRUE)
if (nrow(disp_anova)) write_df(disp_anova, "temporal_dispersion_anova.csv", DIR_ANOVAS)

# ── 2c. Latency to first attack (survival) ───────────────────────────────────
# One row per trial×segment; latency = first within-segment event time; a cell with
# zero events is right-censored at SEG_LEN.
lat <- rbindlist(lapply(trial_ids, function(tr) {
  rbindlist(lapply(c(2L, 3L), function(sg) {
    e <- ev_time[trial == tr & segment == sg]; md <- meta_for_trial(tr)
    if (nrow(e) == 0L)
      data.table(trial = tr, segment = sg, group = md$group,
                 treatment_clean = md$treatment_clean, latency = SEG_LEN, event = 0L)
    else
      data.table(trial = tr, segment = sg, group = md$group,
                 treatment_clean = md$treatment_clean, latency = min(e$t), event = 1L)
  }))
}), fill = TRUE)
write_df(lat, "latency_first_attack.csv", DIR_ANOVAS)
cat(sprintf("\n  Latency to first attack: %d cells, %d observed / %d censored; median = %.1f s\n",
            nrow(lat), sum(lat$event), sum(lat$event == 0), median(lat$latency)))

cox_lat <- tryCatch(
  survival::coxph(survival::Surv(latency, event) ~ group * treatment_clean +
                  cluster(trial), data = lat),
  error = function(e) { message("  Cox failed: ", e$message); NULL })
if (!is.null(cox_lat)) {
  cx  <- summary(cox_lat)
  cxt <- cbind(term = rownames(cx$coefficients), as.data.frame(cx$coefficients))
  write_df(cxt, "latency_coxph.csv", DIR_POSTHOC)
  zph <- tryCatch(survival::cox.zph(cox_lat), error = function(e) NULL)
  if (!is.null(zph)) {
    zdf <- cbind(term = rownames(zph$table), as.data.frame(zph$table))
    write_df(zdf, "latency_coxph_zph.csv", DIR_POSTHOC)
  }
  safe_write_lines(capture.output(print(cx)),
                   file.path(DIR_REPORTS, "latency_coxph_summary.txt"))
}
# Weibull AFT as a parametric cross-check.
aft_lat <- tryCatch(
  survival::survreg(survival::Surv(latency, event) ~ group * treatment_clean,
                    data = lat, dist = "weibull"),
  error = function(e) NULL)
if (!is.null(aft_lat))
  safe_write_lines(capture.output(summary(aft_lat)),
                   file.path(DIR_REPORTS, "latency_weibull_aft_summary.txt"))

# KM curves by treatment.
km <- tryCatch(survival::survfit(survival::Surv(latency, event) ~ treatment_clean, data = lat),
               error = function(e) NULL)
if (!is.null(km)) {
  grDevices::png(file.path(DIR_FIGURES, "latency_km_by_treatment.png"),
                 width = 1000, height = 700, res = 130)
  plot(km, col = c("#2E75B6", "#D55E00"), lwd = 2, mark.time = TRUE,
       xlab = "Time within segment (s)", ylab = "P(no attack yet)",
       main = "Latency to first attack (Kaplan–Meier)")
  legend("topright", legend = levels(lat$treatment_clean),
         col = c("#2E75B6", "#D55E00"), lwd = 2, bty = "n")
  grDevices::dev.off()
  cat("  wrote figures/latency_km_by_treatment.png\n")
}

# ── 2d. Within-segment escalation vs habituation (Poisson trend on 60 s bins) ─
bin_counts <- rbindlist(lapply(trial_ids, function(tr) {
  rbindlist(lapply(c(2L, 3L), function(sg) {
    e  <- ev_time[trial == tr & segment == sg]; md <- meta_for_trial(tr)
    br <- seq(0, SEG_LEN, by = BIN_W)
    cnts <- as.integer(table(cut(e$t, breaks = br, include.lowest = TRUE,
                                 labels = FALSE, right = TRUE)))
    cnts <- c(cnts, rep(0L, N_BINS - length(cnts)))[seq_len(N_BINS)]
    data.table(trial = tr, segment = sg, group = md$group,
               treatment_clean = md$treatment_clean,
               bin = seq_len(N_BINS), count = cnts)
  }))
}), fill = TRUE)
bin_counts[is.na(count), count := 0L]
write_df(bin_counts, "within_segment_bin_counts.csv", DIR_DATASETS)

m_trend <- tryCatch(
  glmmTMB(count ~ poly(bin, 2) * group * treatment_clean + segment + (1 | trial),
          data = bin_counts, family = poisson),
  error = function(e) { message("  trend model failed: ", e$message); NULL })
if (!is.null(m_trend)) {
  a_tr <- tryCatch(car::Anova(m_trend, type = 3), error = function(e) NULL)
  if (!is.null(a_tr)) {
    a_tr_df <- cbind(term = rownames(a_tr), as.data.frame(a_tr))
    write_df(a_tr_df, "within_segment_trend_anova.csv", DIR_ANOVAS)
    cat("  Within-segment trend (Poisson, 60 s bins) — Type III ANOVA written.\n")
  }
}

# Mean event rate per bin (escalation/habituation curve).
trend_curve <- bin_counts[, .(mean_count = mean(count), sem = sd(count) / sqrt(.N)),
                          by = .(treatment_clean, bin)]
p_trend <- ggplot(trend_curve, aes(bin, mean_count, colour = treatment_clean)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = mean_count - sem, ymax = mean_count + sem), width = 0.2) +
  scale_x_continuous(breaks = 1:N_BINS,
                     labels = paste0(seq(0, SEG_LEN - BIN_W, BIN_W), "–",
                                     seq(BIN_W, SEG_LEN, BIN_W))) +
  scale_colour_manual(values = c("Control" = "#2E75B6", "Exercise" = "#D55E00")) +
  labs(x = "Within-segment 60 s bin (s)", y = "Mean aggression events / bin",
       colour = "Treatment",
       title = "Within-segment temporal trend (escalation vs habituation)") +
  BASE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_FIGURES, "within_segment_trend.png"), p_trend,
       width = 9, height = 5, dpi = 200)
cat("  wrote figures/within_segment_trend.png\n")

# Raster of event times per trial×segment.
ev_time[, cell := factor(sprintf("T%02d s%d", trial, segment))]
p_raster <- ggplot(ev_time, aes(t, cell, colour = treatment_clean)) +
  geom_point(shape = 124, size = 2.4) +
  scale_colour_manual(values = c("Control" = "#2E75B6", "Exercise" = "#D55E00")) +
  labs(x = "Time within segment (s)", y = NULL, colour = "Treatment",
       title = "Aggression event raster (each tick = one event)") +
  BASE_THEME
ggsave(file.path(DIR_FIGURES, "event_raster.png"), p_raster,
       width = 10, height = 7, dpi = 200)
cat("  wrote figures/event_raster.png\n")


################################################################################
###  SECTION 3 — CONSOLIDATED REPORT
################################################################################

cat("\n\n########## SECTION 3 — REPORT ##########\n")

# Plain-language verdicts.
verdict_line <- function(name, p, what) {
  if (is.na(p)) return(sprintf("%s: not testable (insufficient data).", name))
  sprintf("%s: %s (p = %s).", name,
          if (p < ALPHA) paste0(what, " — SUPPORTED") else paste0(what, " — not supported"),
          fmt_p_num(p))
}
verdicts <- c(
  "SECTION 1 — DOMINANCE HIERARCHY",
  verdict_line("  Steepness ~ group",      perm_tbl[metric == "steepness"]$group_p_perm,     "group changes steepness"),
  verdict_line("  Steepness ~ treatment",  perm_tbl[metric == "steepness"]$treatment_p_perm, "treatment changes steepness"),
  verdict_line("  Linearity(h') ~ group",  perm_tbl[metric == "landau_hprime"]$group_p_perm,     "group changes linearity"),
  verdict_line("  Linearity(h') ~ treatment", perm_tbl[metric == "landau_hprime"]$treatment_p_perm, "treatment changes linearity"),
  verdict_line("  DCI ~ treatment",        perm_tbl[metric == "DCI"]$treatment_p_perm, "treatment changes directional consistency"),
  sprintf("  Rank consistency seg2↔seg3: mean Spearman ρ = %s over %d trials.",
          fmt_num(mean(rank_consist$spearman_rho, na.rm = TRUE)), nrow(rank_consist)),
  verdict_line("  Given ≠ received (paired, per fish)", given_received_paired$wilcox_p,
              "a fish's given count differs from its own received count"),
  sprintf("  Mean given = %s, mean received = %s, mean diff = %s (per fish, per trial).",
          fmt_num(given_received_paired$mean_given), fmt_num(given_received_paired$mean_received),
          fmt_num(given_received_paired$mean_diff)),
  verdict_line("  Net aggression ~ group×treatment", {
    an <- a_net_df; if (is.null(an)) NA_real_ else {
      pp <- suppressWarnings(as.numeric(an[["Pr(>Chisq)"]][an$term == "group:treatment_clean"]))
      if (length(pp) == 1) pp else NA_real_ } },
    "group×treatment changes the given-minus-received balance"),
  verdict_line("  Given↔received correlation", given_received_cor$p_value,
              if (!is.na(given_received_cor$spearman_rho) && given_received_cor$spearman_rho > 0)
                "reciprocal/mutual conflict (attackers are also targeted)"
              else "one-sided despotism (attackers and targets are distinct roles)"),
  "",
  "SECTION 2 — TEMPORAL DYNAMICS",
  verdict_line("  Bursty vs Poisson (CV≠1)", null_tests[metric == "cv_iei"]$wilcox_p, "events depart from Poisson timing"),
  verdict_line("  Clustering (Fano≠1)",      null_tests[metric == "fano"]$wilcox_p,   "events cluster in time"),
  sprintf("  Mean CV = %s, mean burstiness B = %s (Poisson null: CV=1, B=0).",
          fmt_num(mean(disp_ok$cv_iei, na.rm = TRUE)), fmt_num(mean(disp_ok$burstiness_B, na.rm = TRUE))),
  sprintf("  Latency: %d/%d cells had an attack; median latency = %s s.",
          sum(lat$event), nrow(lat), fmt_num(median(lat$latency))))
safe_write_lines(verdicts, file.path(DIR_REPORTS, "verdicts.txt"))
cat(paste(verdicts, collapse = "\n"), "\n")

# Markdown + Word report.
md <- c(
  "# Aggression — dominance hierarchy & temporal dynamics",
  sprintf("_Generated %s. Input: %s. Dyadic events: %d/%d._",
          TIMESTAMP, basename(MASTER_CSV), nrow(events_dyad), nrow(events_full)),
  "", "## Section 1 — per-trial hierarchy metrics",
  df_to_md(hier_trial[, .(trial, group, treatment_clean, n_events,
                          steepness, steep_p, landau_hprime, DCI, ttri)]),
  "## Section 1 — group/treatment permutation tests",
  df_to_md(perm_tbl),
  "## Section 1 — given vs received aggression (paired test, per fish)",
  df_to_md(given_received_paired),
  "## Section 1 — given↔received correlation across fish",
  df_to_md(given_received_cor),
  "## Section 2 — temporal dispersion vs Poisson null",
  df_to_md(null_tests),
  "## Section 2 — latency (per trial×segment)",
  df_to_md(lat),
  "## Verdicts", paste0(verdicts, collapse = "\n"))
safe_write_lines(md, file.path(DIR_REPORTS, "social_temporal_report.md"))

tryCatch({
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Aggression — Dominance Hierarchy & Temporal Dynamics",
                               style = "heading 1")
  doc <- officer::body_add_par(doc, sprintf(
    "Generated %s. Input: %s. Dyadic events: %d/%d (%d excluded for unmatched identity).",
    TIMESTAMP, basename(MASTER_CSV), nrow(events_dyad), nrow(events_full), n_excluded),
    style = "Normal")
  doc <- officer::body_add_par(doc, "Section 1 — Dominance hierarchy (per trial)", style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(as.data.frame(
    hier_trial[, .(trial, group, treatment_clean, n_events, steepness, steep_p,
                   landau_hprime, DCI, ttri)]), 8L))
  doc <- officer::body_add_par(doc, "Group / treatment permutation tests", style = "heading 3")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(as.data.frame(perm_tbl), 7L))
  doc <- officer::body_add_par(doc, "Given vs received aggression (paired test, per fish)", style = "heading 3")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(given_received_paired, 8L))
  doc <- officer::body_add_par(doc, "Given↔received correlation across fish", style = "heading 3")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(given_received_cor, 8L))
  doc <- officer::body_add_par(doc, "Section 2 — Temporal dispersion vs Poisson null", style = "heading 2")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(as.data.frame(null_tests), 8L))
  doc <- officer::body_add_par(doc, "Latency to first attack", style = "heading 3")
  doc <- flextable::body_add_flextable(doc, .df_to_flex(as.data.frame(lat), 8L))
  doc <- officer::body_add_par(doc, "Verdicts", style = "heading 2")
  for (v in verdicts) doc <- officer::body_add_par(doc, v, style = "Normal")
  print(doc, target = file.path(DIR_REPORTS, "social_temporal_report.docx"))
  cat("  wrote reports/social_temporal_report.docx\n")
}, error = function(e) message("  docx report failed: ", e$message))

cat("\n=== DONE ===\n")
cat("All outputs under:", OUTPUT_DIR, "\n")

# HANDOFF.md S10 [INV]: completed-run sentinel for downstream auto-discovery.
writeLines(as.character(Sys.time()), file.path(OUTPUT_DIR, "_RUN_COMPLETE"))

