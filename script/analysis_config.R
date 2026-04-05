# analysis_config.R — default settings for PD + AF OAC sequential trials
# Source from script/01, 02, or 03 (project root as working directory).
# Override before sourcing, e.g.: analysis_config$exposure_mode <- "current_oac_interval"

analysis_config <- list(
  # "manuscript_initiator" = Word / script/run_incident.R style: monthly trials,
  #   prevalent OAC excluded (Rx ending <= index - 1 month), new-user expo =
  #   (index - oac_initiation_lookback_months, index] per run_incident line 274.
  # "current_oac_interval" = aligns with script/run.R (daily index, Rx interval overlap).
  exposure_mode = "manuscript_initiator",

  # Prevalent OAC: exclude if date.oac <= index - this many months (run_incident: 1).
  oac_prevalent_exclusion_months = 1L,
  # New-user window for expo==1: date.oac in (index - m months, index] (run_incident: 2).
  oac_initiation_lookback_months = 2L,
  # If TRUE, trial rows with dx.stroke_embo==1 are dropped in 01. run_incident.R
  # does NOT do this; it relies on pt.allstroke==0 in 02. Keep FALSE for Word parity.
  cohort_exclude_dx_stroke_embo = FALSE,

  # "month" = one trial per calendar month (~191 after drop_first_index).
  # "day" = one trial per calendar day (large; use for daily person-time in 03).
  index_by = "day",

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
  path_person_trial_rds = "out/person_trial_primary.rds",

  # 03_regression_pooled_logistic.R
  # If both FALSE, 03 stops. Default: speedglm coefs + cluster bootstrap RD/RR (run_incident).
  regression_run_speedglm = TRUE,
  regression_run_bootstrap = TRUE,
  bootstrap_r = 500L,
  bootstrap_seed = 456L,
  path_est_pooled_expo_csv = "out/est_pooled_logistic_expo.csv",
  path_est_pooled_all_csv = "out/est_pooled_logistic_all_terms.csv",
  path_bootstrap_risk_csv = "out/bootstrap_pooled_risk_run_incident.csv",
  # Long-format TableOne used for SMD in 03 (person-month rows).
  path_tableone_long_smd_csv = "out/tableone_long_format_for_smd.csv",

  # Monthly person-month outcome `out` (see script/run_incident.R out_indicator == 1)
  # "pmin_composite" — first event time t = pmin(ische, haem, death); out = 1 if
  #   t >= obs_date and month(t) == month(obs_date). Matches run_incident primary.
  # "any_event_same_month" — out = 1 if any of ischemic, haem, or death falls in
  #   the same calendar month as obs_date (differs from pmin when multiple events).
  outcome_definition = "pmin_composite",

  # "default" = script/02 as originally written (stroke in PS, caliper, replace=T).
  # "run_incident" = replicate script/run_incident.R: obs_end uses ische/haem/death;
  #   MatchIt: fu in PS, no dx.stroke_embo in PS, replace=F, exact ~ trial_id, no caliper;
  #   Table 1 before match omits dx.stroke_embo; 03 uses SMD from long-format data.
  pipeline_style = "run_incident"
)

load_analysis_config <- function(root = ".") {
  cfg <- analysis_config
  if (is.null(cfg$path_seqcohort_rds)) {
    suf <- sprintf("%s_%s", cfg$exposure_mode, cfg$index_by)
    cfg$path_seqcohort_rds <- file.path(
      root, "data", paste0("seqcohort_", suf, ".RDS")
    )
  }
  # Defaults if an older analysis_config list omits new keys
  if (is.null(cfg$oac_prevalent_exclusion_months)) {
    cfg$oac_prevalent_exclusion_months <- 1L
  }
  if (is.null(cfg$oac_initiation_lookback_months)) {
    cfg$oac_initiation_lookback_months <- 2L
  }
  if (is.null(cfg$cohort_exclude_dx_stroke_embo)) {
    cfg$cohort_exclude_dx_stroke_embo <- FALSE
  }
  if (is.null(cfg$path_tableone_long_smd_csv)) {
    cfg$path_tableone_long_smd_csv <- "out/tableone_long_format_for_smd.csv"
  }
  if (is.null(cfg$path_bootstrap_risk_csv)) {
    cfg$path_bootstrap_risk_csv <- "out/bootstrap_pooled_risk_run_incident.csv"
  }
  if (is.null(cfg$regression_run_speedglm)) {
    cfg$regression_run_speedglm <- TRUE
  }
  if (is.null(cfg$regression_run_bootstrap)) {
    cfg$regression_run_bootstrap <- TRUE
  }
  if (is.null(cfg$bootstrap_r)) {
    cfg$bootstrap_r <- 500L
  }
  if (is.null(cfg$bootstrap_seed)) {
    cfg$bootstrap_seed <- 456L
  }
  if (is.null(cfg$index_by)) {
    cfg$index_by <- "day"
  }
  if (!cfg$index_by %in% c("day", "month")) {
    stop('analysis_config: index_by must be "day" or "month".')
  }
  cfg
}
