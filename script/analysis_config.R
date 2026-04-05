# analysis_config.R — default settings for PD + AF OAC sequential trials
# Source from script/01_build_seqcohort.R or script/02_analysis_primary.R
# Override before sourcing, e.g.: analysis_config$exposure_mode <- "current_oac_interval"

analysis_config <- list(
  # "manuscript_initiator" = Word methods (monthly trials, 1-mo initiation window,
  #   prevalent OAC exclusion, expo from cohort date.oac).
  # "current_oac_interval" = aligns with script/run.R (daily index, Rx interval overlap).
  exposure_mode = "manuscript_initiator",

  # "month" = 191 trials (after drop_first_index); "day" = daily index dates.
  index_by = "month",

  study_start = as.Date("2007-01-01"),
  study_end = as.Date("2022-12-31"),

  # If TRUE, first calendar index (2007-01-01) is dropped (matches run.R index_date[-1]).
  drop_first_index = TRUE,

  oac_end_impute_days = 90L,
  verbose_trial_msg = FALSE,

  # Cap follow-up in obs_end (pmin); manuscript primary outcome text uses 3 years.
  # script/run.R uses 5 — override here to match manuscript vs legacy script.
  obs_end_cap_years = 3,

  # Paths relative to project root (set working directory to project root before run).
  path_cohort_rds = "data/pd_oac-cohort-20240813.RDS",
  path_dx_rds = "data/PD_OAC-dx.RDS",
  path_drug_rds = "data/PD_OAC-drug.RDS",
  path_demo_rds = "data/PD_OAC-demo.RDS",
  path_protocol_xlsx = "documents/PD OAC-Protocol.xlsx",
  protocol_covar_sheet = "Dx.Cov",

  # Output seqcohort path (01_build_seqcohort.R); derived if NULL.
  path_seqcohort_rds = NULL,

  # MatchIt (02_analysis_primary.R)
  matchit_seed = 456L,
  matchit_ratio = 10L,
  matchit_caliper = 0.2,
  matchit_replace = TRUE,

  path_tableone_before_csv = "out/tableone_before_matching.csv",
  path_tableone_after_csv = "out/tableone_after_matching.csv",
  path_person_trial_rds = "out/person_trial_primary.rds"
)

load_analysis_config <- function(root = ".") {
  cfg <- analysis_config
  if (is.null(cfg$path_seqcohort_rds)) {
    suf <- sprintf("%s_%s", cfg$exposure_mode, cfg$index_by)
    cfg$path_seqcohort_rds <- file.path(
      root, "data", paste0("seqcohort_", suf, ".RDS")
    )
  }
  cfg
}
