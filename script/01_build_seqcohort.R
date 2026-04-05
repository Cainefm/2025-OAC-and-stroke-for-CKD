# 01_build_seqcohort.R — sequential target-trial cohort (PD + AF + OAC)
# Run from project root:  Rscript script/01_build_seqcohort.R
# Requires: data/*.RDS, documents/PD OAC-Protocol.xlsx
# Settings: script/analysis_config.R (override before sourcing if needed)

options(scipen = 6L, digits = 4L)
if (.Platform$OS.type == "windows") {
  try(memory.limit(30000000L), silent = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(lubridate)
  library(pbapply)
  library(readxl)
  library(tibble)
})

source("script/00_functions_pd_oac.R")
root <- project_root()

source(file.path(root, "script/analysis_config.R"))
cfg <- load_analysis_config(root = root)

cohort <- setDT(readRDS(file.path(root, cfg$path_cohort_rds)))
ymd_cols <- grep("date", names(cohort), value = TRUE)
cohort[, (ymd_cols) := lapply(.SD, ymd), .SDcols = ymd_cols]
cohort[, date.PD.end := pmin(date.transplant, date.PD.end, date.hd, na.rm = TRUE)]

covar <- setDT(
  readxl::read_xlsx(
    file.path(root, cfg$path_protocol_xlsx),
    sheet = cfg$protocol_covar_sheet
  )
)
dx <- setDT(readRDS(file.path(root, cfg$path_dx_rds)))
dx[, Reference_Date := ymd(Reference_Date)]
dx <- dx[grepl(covar[, paste(Regex, collapse = "|")], All_Diagnosis_Code_ICD9)]
dx <- dx[, .(Reference_Key, Reference_Date, All_Diagnosis_Code_ICD9)][,
         .SD[Reference_Date == min(Reference_Date)],
         .(Reference_Key, All_Diagnosis_Code_ICD9)]

rx <- readRDS(file.path(root, cfg$path_drug_rds))
setDT(rx)

drug_oac <- if (cfg$exposure_mode == "current_oac_interval") {
  prepare_drug_oac(rx, cfg$oac_end_impute_days)
} else {
  NULL
}

collapse_ap <- (cfg$index_by == "month")
drug_antiplatelet <- prepare_antiplatelet(rx, collapse_to_month = collapse_ap)

index_date <- build_index_dates(cfg)
message(
  "Building seqcohort: mode=", cfg$exposure_mode,
  ", index_by=", cfg$index_by,
  ", n_index=", length(index_date),
  ", oac_init_months=", cfg$oac_initiation_lookback_months,
  ", exclude_dx_stroke_embo=", cfg$cohort_exclude_dx_stroke_embo,
  ", output=", cfg$path_seqcohort_rds
)

run_seq <- make_run_seq(cfg, dx, covar, drug_oac)
seqcohort <- pblapply(index_date, run_seq, cohort)

out_seq <- file.path(root, cfg$path_seqcohort_rds)
saveRDS(seqcohort, out_seq)
message("Saved: ", out_seq)

ap_path <- sub("\\.RDS$", "_antiplatelet.RDS", cfg$path_seqcohort_rds)
saveRDS(drug_antiplatelet, file.path(root, ap_path))
message("Saved: ", file.path(root, ap_path))
