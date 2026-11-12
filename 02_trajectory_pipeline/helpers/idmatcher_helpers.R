# =============================================================================
# 00_idmatcher_helpers.R
#
# Helper functions for the STEP 1.5 idmatcher.ai integration:
#
#   find_latest_idmatcher_run()        – auto-detect the most recent Colab output
#   assign_fish_id_matched_trial()     – assign trial-scoped stable identities
#   recalculate_body_length()          – compute mean BL across matched segments
#   compute_idmatcher_reliability_tier() – per-trial tier for RE gating
#
# Sourced by step1_5_annotate_idmatcher_GOT.R (and by 00_master_pipeline_GOT.R
# so the helpers are in GlobalEnv before any step script runs).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ── Tier thresholds (tune here, not downstream) ───────────────────────────────
.IDMATCHER_TIER_CONSERVATIVE_STRONG_RATE <- 1.0
.IDMATCHER_TIER_RELIABLE_ACCEPT_RATE     <- 0.75
.IDMATCHER_TIER_BORDERLINE_ACCEPT_RATE   <- 0.50
# conservative additionally requires: accept_rate == 1 AND no rejected_conflict


# =============================================================================
# find_latest_idmatcher_run()
#
# Scans `root` for subfolders matching the YYYYMMDD_HHMMSS timestamp format
# emitted by run_idmatcherai_by_trial.py.  Returns the path to the most recent
# one (lexicographic sort == chronological sort for this format).
#
# Stops loudly when:
#   - root does not exist
#   - no timestamped subfolders are found
#   - the latest folder is missing idmatcher_identity_correspondences.csv or
#     run_report.json
# =============================================================================
find_latest_idmatcher_run <- function(
    root = file.path(get("PROJECT_ROOT", envir = .GlobalEnv, inherits = TRUE),
                     "idmatcherai_colab", "idmatcherai outputs from colab")
) {
  if (!dir.exists(root))
    stop("idmatcher Colab root not found: ", root)

  subs    <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  ts_subs <- subs[grepl("^\\d{8}_\\d{6}$", basename(subs))]

  if (length(ts_subs) == 0)
    stop("No YYYYMMDD_HHMMSS subfolders found under: ", root)

  latest <- ts_subs[order(basename(ts_subs), decreasing = TRUE)][1]

  required <- c("idmatcher_identity_correspondences.csv", "run_report.json")
  missing  <- required[!file.exists(file.path(latest, required))]
  if (length(missing))
    stop("Latest idmatcher run (", basename(latest), ") is missing: ",
         paste(missing, collapse = ", "))

  attr(latest, "has_matrices") <- dir.exists(file.path(latest, "matrices"))
  attr(latest, "has_logs")     <- dir.exists(file.path(latest, "logs"))
  latest
}


# =============================================================================
# assign_fish_id_matched_trial()
#
# Takes the annotated long table from annotate_idmatcher_scores_run() and:
#
#  1. Anchors each seg2 fish as T<trial>_MF<fish_id_master>
#     (MF = matched fish identity).
#  2. Binds accepted seg3 fish to the same fish_id_matched_trial.
#  3. Resolves rejected_conflict rows when reliable anchors exist in the trial:
#     the highest-scoring assignment (by direct_score) in the conflicting group
#     is accepted as truth and labelled "accepted_from_conflict".
#  4. Detects structural many-to-one / one-to-many duplicates after step 3 and
#     demotes losers to match_quality = "rejected_duplicate".
#  5. Emits a QC warning for any trial with conflicts or duplicates.
#
# Returns: annotated tibble with new columns:
#   fish_id_matched_trial  – e.g. "T5_MF1" (NA for unmatched fish)
#   tracker_fish_id_seg2   – paste0("T",trial_id,"_S2_F",fish_id_master)
#   tracker_fish_id_seg3   – paste0("T",trial_id,"_S3_F",fish_id)
# =============================================================================
assign_fish_id_matched_trial <- function(annotated) {

  annotated <- annotated |>
    mutate(
      trial_id = as.integer(trial),
      fish_id_master = as.integer(fish_id_master),
      fish_id        = as.integer(fish_id)
    )

  # ── Step 3: resolve rejected_conflict when reliable anchors exist ──────────
  annotated <- annotated |>
    group_by(trial_id, segment) |>
    mutate(
      n_reliable = sum(match_quality %in% c("strong", "acceptable",
                                            "by_exclusion", "by_consensus",
                                            "by_matrix", "by_bodylength"),
                       na.rm = TRUE)
    ) |>
    ungroup()

  conflict_idx <- which(annotated$match_quality == "rejected_conflict" &
                        annotated$n_reliable > 0)

  if (length(conflict_idx) > 0) {
    # Within each (trial_id, segment) conflict group, keep highest direct_score
    conflict_df <- annotated[conflict_idx, ] |>
      group_by(trial_id, segment) |>
      mutate(
        rank_in_group = rank(-direct_score, ties.method = "first"),
        match_quality = if_else(rank_in_group == 1,
                                "accepted_from_conflict", "rejected_conflict"),
        match_accepted = if_else(rank_in_group == 1, TRUE, FALSE)
      ) |>
      ungroup() |>
      select(-rank_in_group)

    annotated[conflict_idx, ] <- conflict_df

    n_resolved <- sum(conflict_df$match_quality == "accepted_from_conflict")
    warning(sprintf(
      "STEP 1.5: resolved %d rejected_conflict row(s) as 'accepted_from_conflict' — check idmatcher_diagnostics output.",
      n_resolved
    ), call. = FALSE)
  }

  # ── Step 1: anchor seg2 fish ───────────────────────────────────────────────
  seg2_anchors <- annotated |>
    distinct(trial_id, fish_id_master) |>
    mutate(
      fish_id_matched_trial = paste0("T", trial_id, "_MF", fish_id_master),
      tracker_fish_id_seg2  = paste0("T", trial_id, "_S2_F", fish_id_master)
    )

  # ── Step 2: bind accepted seg3 fish ───────────────────────────────────────
  annotated <- annotated |>
    left_join(
      seg2_anchors |> select(trial_id, fish_id_master, fish_id_matched_trial,
                             tracker_fish_id_seg2),
      by = c("trial_id", "fish_id_master")
    ) |>
    mutate(
      fish_id_matched_trial = if_else(
        !is.na(match_accepted) & match_accepted,
        fish_id_matched_trial,
        NA_character_
      ),
      tracker_fish_id_seg3 = paste0("T", trial_id, "_S", segment, "_F", fish_id)
    )

  # ── Step 4: detect structural duplicates ──────────────────────────────────
  # Many-to-one: same tracker_fish_id_seg3 claimed by > 1 fish_id_matched_trial
  dupes_seg3 <- annotated |>
    filter(!is.na(fish_id_matched_trial)) |>
    add_count(trial_id, segment, fish_id, name = "n_seg3_claims") |>
    filter(n_seg3_claims > 1)

  # One-to-many: same fish_id_matched_trial claims > 1 fish_id in seg3
  dupes_seg2 <- annotated |>
    filter(!is.na(fish_id_matched_trial)) |>
    add_count(trial_id, fish_id_matched_trial, segment, name = "n_seg2_claims") |>
    filter(n_seg2_claims > 1)

  dup_rows <- union(rownames(dupes_seg3), rownames(dupes_seg2))

  if (length(dup_rows) > 0) {
    dup_trials <- unique(c(dupes_seg3$trial_id, dupes_seg2$trial_id))
    warning(sprintf(
      "STEP 1.5: structural duplicate assignments detected in trial(s) %s — demoting losers to 'rejected_duplicate'. Check idmatcher_diagnostics output.",
      paste(dup_trials, collapse = ", ")
    ), call. = FALSE)

    # Keep highest direct_score per (trial_id, segment, fish_id) conflict
    annotated <- annotated |>
      group_by(trial_id, segment, fish_id) |>
      mutate(
        n_claims = sum(!is.na(fish_id_matched_trial)),
        dup_rank = if_else(
          n_claims > 1,
          rank(-replace(direct_score, is.na(direct_score), -Inf),
               ties.method = "first"),
          1L
        ),
        match_quality = if_else(
          n_claims > 1 & dup_rank > 1,
          "rejected_duplicate", match_quality
        ),
        match_accepted = if_else(match_quality == "rejected_duplicate", FALSE, match_accepted),
        fish_id_matched_trial = if_else(
          match_quality == "rejected_duplicate", NA_character_, fish_id_matched_trial
        )
      ) |>
      select(-n_claims, -dup_rank) |>
      ungroup()
  }

  # ── Fish absent from correspondences table (unmatched_no_correspondence) ──
  # These would be seg3 fish tracked but never appearing in the wide CSV.
  # They can't appear in `annotated` (which came from the wide CSV), so this
  # label is applied by the STEP 1.5 orchestrator after building idmatcher_id_map.

  # ── Cleanup ────────────────────────────────────────────────────────────────
  annotated |>
    select(-n_reliable) |>
    mutate(trial_id = as.integer(trial_id))
}


# =============================================================================
# recalculate_body_length()
#
# Joins per-fish BL estimates from master_bodylength_by_fish.csv onto the
# annotated (+ fish_id_matched_trial) table, computes the mean BL across
# accepted, stable segments, and returns a per-fish lookup table.
#
# Parameters:
#   annotated    tibble from assign_fish_id_matched_trial()
#   bl_csv_path  path to master_bodylength_by_fish.csv
#
# Returns: tibble with columns:
#   fish_id_matched_trial, trial_id,
#   body_length_seg2_cm, body_length_seg3_cm,
#   body_length_mean_cm, body_length_delta_cm,
#   body_length_unstable_seg2, body_length_unstable_seg3,
#   body_length_source
# =============================================================================
recalculate_body_length <- function(annotated, bl_csv_path) {

  bl <- read_csv(bl_csv_path, show_col_types = FALSE) |>
    select(trial_id, segment_id, fish_in_session,
           body_length_cm, body_length_unstable)

  # Accepted rows only (master seg2 anchor + accepted seg3)
  acc <- annotated |>
    filter(!is.na(fish_id_matched_trial)) |>
    select(trial_id, fish_id_matched_trial,
           fish_id_master, master_segment, segment, fish_id)

  # Join BL for seg2 (master)
  seg2_bl <- acc |>
    distinct(trial_id, fish_id_matched_trial, fish_id_master, master_segment) |>
    left_join(bl, by = c("trial_id",
                         "master_segment" = "segment_id",
                         "fish_id_master"  = "fish_in_session")) |>
    rename(body_length_seg2_cm     = body_length_cm,
           body_length_unstable_seg2 = body_length_unstable)

  # Join BL for seg3 (matched)
  seg3_rows <- acc |>
    filter(segment != master_segment) |>
    distinct(trial_id, fish_id_matched_trial, segment, fish_id)

  seg3_bl <- seg3_rows |>
    left_join(bl, by = c("trial_id",
                         "segment"  = "segment_id",
                         "fish_id"  = "fish_in_session")) |>
    rename(body_length_seg3_cm     = body_length_cm,
           body_length_unstable_seg3 = body_length_unstable)

  # Combine and compute mean
  bl_lookup <- seg2_bl |>
    left_join(seg3_bl |> select(trial_id, fish_id_matched_trial,
                                body_length_seg3_cm, body_length_unstable_seg3),
              by = c("trial_id", "fish_id_matched_trial")) |>
    mutate(
      # prefer stable segments; if unstable, set to NA for the first pass
      seg2_for_mean = if_else(!is.na(body_length_unstable_seg2) & body_length_unstable_seg2,
                              NA_real_, body_length_seg2_cm),
      seg3_for_mean = if_else(!is.na(body_length_unstable_seg3) & body_length_unstable_seg3,
                              NA_real_, body_length_seg3_cm),
      body_length_mean_cm = case_when(
        # at least one stable segment available — use only stable values
        !is.na(seg2_for_mean) & !is.na(seg3_for_mean) ~
          (seg2_for_mean + seg3_for_mean) / 2,
        !is.na(seg2_for_mean) ~ seg2_for_mean,
        !is.na(seg3_for_mean) ~ seg3_for_mean,
        # both segments unstable — fall back to raw (unstable) values rather than NA
        !is.na(body_length_seg2_cm) & !is.na(body_length_seg3_cm) ~
          (body_length_seg2_cm + body_length_seg3_cm) / 2,
        !is.na(body_length_seg2_cm) ~ body_length_seg2_cm,
        !is.na(body_length_seg3_cm) ~ body_length_seg3_cm,
        TRUE ~ NA_real_
      ),
      body_length_delta_cm = abs(body_length_seg2_cm - body_length_seg3_cm),
      body_length_source   = case_when(
        !is.na(seg2_for_mean) & !is.na(seg3_for_mean) ~ "seg2_seg3_mean",
        !is.na(seg2_for_mean)                          ~ "seg2_only",
        !is.na(seg3_for_mean)                          ~ "seg3_only",
        !is.na(body_length_mean_cm)                    ~ "unstable_fallback",
        TRUE                                           ~ "no_data"
      )
    ) |>
    select(fish_id_matched_trial, trial_id,
           body_length_seg2_cm, body_length_seg3_cm,
           body_length_mean_cm, body_length_delta_cm,
           body_length_unstable_seg2, body_length_unstable_seg3,
           body_length_source) |>
    select(-starts_with("seg2_for_mean"), -starts_with("seg3_for_mean"))

  # Warn about fish where both segments are unstable (BL still computed from raw values)
  all_unstable <- bl_lookup |>
    filter(body_length_unstable_seg2 == TRUE & body_length_unstable_seg3 == TRUE)
  if (nrow(all_unstable) > 0)
    warning(sprintf(
      "STEP 1.5: %d matched fish have unstable BL in BOTH seg2 and seg3 — body_length_mean_cm computed from raw (unstable) values (body_length_source = 'unstable_fallback').",
      nrow(all_unstable)
    ), call. = FALSE)

  bl_lookup
}


# =============================================================================
# compute_idmatcher_reliability_tier()
#
# Derives one reliability tier per trial from the annotated + assigned table.
#
# Tier thresholds (tunable via module-level constants):
#   conservative  accept_rate == 1 AND strong_rate == 1 AND no rejected_conflict
#   reliable      accept_rate >= 0.75 AND no rejected_conflict
#   borderline    accept_rate >= 0.50 OR (>= 0.75 AND some rejected_conflict)
#   unreliable    otherwise
#
# Returns: tibble with columns:
#   trial_id, idmatcher_tier,
#   n_total, n_accepted, n_strong, accept_rate, strong_rate, n_rejected_conflict
# =============================================================================
compute_idmatcher_reliability_tier <- function(annotated) {

  annotated |>
    group_by(trial_id) |>
    summarise(
      n_total            = n(),
      n_accepted         = sum(!is.na(match_accepted) & match_accepted),
      n_strong           = sum(match_quality %in% c("strong", "accepted_from_conflict"),
                               na.rm = TRUE),
      n_rejected_conflict = sum(match_quality == "rejected_conflict", na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      accept_rate = if_else(n_total > 0, n_accepted / n_total, NA_real_),
      strong_rate = if_else(n_total > 0, n_strong  / n_total, NA_real_),
      idmatcher_tier = dplyr::case_when(
        accept_rate == .IDMATCHER_TIER_CONSERVATIVE_STRONG_RATE &
          strong_rate == .IDMATCHER_TIER_CONSERVATIVE_STRONG_RATE &
          n_rejected_conflict == 0 ~ "conservative",
        accept_rate >= .IDMATCHER_TIER_RELIABLE_ACCEPT_RATE &
          n_rejected_conflict == 0 ~ "reliable",
        accept_rate >= .IDMATCHER_TIER_BORDERLINE_ACCEPT_RATE |
          (accept_rate >= .IDMATCHER_TIER_RELIABLE_ACCEPT_RATE &
           n_rejected_conflict > 0) ~ "borderline",
        TRUE ~ "unreliable"
      )
    ) |>
    select(trial_id, idmatcher_tier, n_total, n_accepted, n_strong,
           accept_rate, strong_rate, n_rejected_conflict)
}


# =============================================================================
# build_idmatcher_id_map()
#
# Builds the canonical per-(trial_id, segment, fish_in_session) lookup table
# used by STEP2 and the BL match script.
#
# One row per fish per segment (both seg2 anchors and seg3 matched fish).
# Unmatched seg3 fish that appear in the correspondences table but have
# fish_id_matched_trial == NA are included with match_quality = their label.
#
# Returns: tibble with columns:
#   trial_id, segment, fish_in_session, tracker_fish_id,
#   fish_id_matched_trial, body_length_mean_cm, match_quality
# =============================================================================
build_idmatcher_id_map <- function(annotated, bl_lookup) {

  # Seg3 rows from annotation
  seg3 <- annotated |>
    select(trial_id, segment, fish_in_session = fish_id,
           tracker_fish_id = tracker_fish_id_seg3,
           fish_id_matched_trial, match_quality)

  # Seg2 anchor rows (one per master fish per trial)
  seg2 <- annotated |>
    distinct(trial_id, fish_id_master, master_segment,
             tracker_fish_id_seg2, fish_id_matched_trial) |>
    filter(!is.na(fish_id_matched_trial)) |>
    transmute(
      trial_id,
      segment          = as.integer(master_segment),
      fish_in_session  = as.integer(fish_id_master),
      tracker_fish_id  = tracker_fish_id_seg2,
      fish_id_matched_trial,
      match_quality    = "(seg2 anchor)"
    )

  id_map <- bind_rows(seg2, seg3) |>
    left_join(
      bl_lookup |> select(fish_id_matched_trial, body_length_mean_cm),
      by = "fish_id_matched_trial"
    ) |>
    arrange(trial_id, segment, fish_in_session)

  id_map
}
