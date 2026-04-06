# 02_analysis_primary.R — person_trial, Table 1, MatchIt (matches core of script/run.R)
# Run from project root after 01_build_seqcohort.R:
#   Rscript script/02_analysis_primary.R
# Uses script/analysis_config.R for paths and follow-up cap.
#
# pipeline_style "run_incident" mirrors script/run_incident.R (out_indicator == 1):
#   pt.allstroke == 0, obs_end includes date.ische/haem/death, Table 1 omits
#   dx.stroke_embo before match, MatchIt: fu in PS, no stroke_embo, replace=F,
#   exact ~ trial_id, no caliper.

options(scipen = 6L, digits = 4L)
if (.Platform$OS.type == "windows") {
  try(memory.limit(30000000L), silent = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(lubridate)
  library(tableone)
  library(labelled)
  library(MatchIt)
  library(readxl)
})

source("script/00_functions_pd_oac.R")
root <- project_root()

source(file.path(root, "script/analysis_config.R"))
cfg <- load_analysis_config(root = root)
message("pipeline_style: ", cfg$pipeline_style)

seqcohort <- read_seqcohort_with_trial_ids(
  file.path(root, cfg$path_seqcohort_rds)
)

person_trial <- rbindlist(seqcohort, fill = T)
person_trial[, id := paste(trial_id, Reference_Key, sep = "-")]

demo <- as.data.table(readRDS(file.path(root, cfg$path_demo_rds)))
person_trial <- merge(
  person_trial,
  unique(demo[, .(Reference_Key, Sex, dob = Date_of_Birth_yyyymmdd)]),
  by = "Reference_Key"
)
setkey(person_trial, NULL)
person_trial[, id := paste(trial_id, Reference_Key, sep = "_")]

# Past stroke exclusion (script/run_incident.R): prefer cohort column pt.allstroke.
if (cfg$pipeline_style == "run_incident") {
  if ("pt.allstroke" %in% names(person_trial)) {
    n_ex <- person_trial[pt.allstroke == 1L, .N]
    message("run_incident: excluding pt.allstroke == 1 (n = ", n_ex, ")")
    person_trial <- person_trial[pt.allstroke == 0L]
  } else if ("dx.stroke_embo" %in% names(person_trial)) {
    n_ex <- person_trial[dx.stroke_embo == 1L, .N]
    message(
      "run_incident: cohort has no pt.allstroke; excluding dx.stroke_embo == 1 ",
      "(n = ", n_ex, ") as stroke-history proxy — for exact run_incident parity, ",
      "include pt.allstroke on the cohort RDS."
    )
    person_trial <- person_trial[dx.stroke_embo == 0L]
  } else {
    warning(
      "run_incident style: neither pt.allstroke nor dx.stroke_embo found; ",
      "stroke history not excluded (differs from script/run_incident.R)."
    )
  }
}

cap_years <- cfg$obs_end_cap_years
study_end <- as.Date(cfg$study_end)
if (cfg$pipeline_style == "run_incident") {
  person_trial[, obs_end := pmin(
    study_end,
    date.ische,
    date.haem,
    date.death,
    date.hd,
    date.transplant,
    indx_date %m+% years(cap_years),
    date.PD.end,
    na.rm = TRUE
  )]
} else {
  person_trial[, obs_end := pmin(
    study_end,
    date.stroke,
    date.death,
    date.hd,
    date.transplant,
    indx_date %m+% years(cap_years),
    date.PD.end,
    na.rm = TRUE
  )]
}

person_trial[, age_indx := suppressWarnings(as.numeric(
  (as.Date(indx_date) - as.Date(lubridate::ymd(dob))) / 365.25
))]
person_trial[, time_af_index := as.numeric(indx_date - date.af)]
compute_chadsvas_score(person_trial)

person_trial[, fu := as.numeric(obs_end - indx_date)]

drug_antiplatelet <- load_or_build_antiplatelet_monthly(cfg, root)
person_trial <- merge_antiplatelet_baseline_at_index(
  person_trial,
  drug_antiplatelet
)

vars_default <- c(
  "fu", "Sex", "age_indx", "time_af_index", "chadsvas_score",
  "dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
  "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
  "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
  "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
  "dx.dm_com1", "dx.htn", "dx.constipation", "antiplatelet.baseline"
)
vars_run_incident <- c(
  "fu", "Sex", "age_indx", "time_af_index", "chadsvas_score",
  "dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
  "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
  "dx.liver_modsev", "dx.ulcers", "dx.asthma",
  "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
  "dx.dm_com1", "dx.htn", "dx.constipation", "antiplatelet.baseline"
)
vars <- if (cfg$pipeline_style == "run_incident") vars_run_incident else vars_default

covar <- setDT(readxl::read_xlsx(
  file.path(root, cfg$path_protocol_xlsx),
  sheet = cfg$protocol_covar_sheet
))
apply_protocol_var_labels(person_trial, covar)
set_core_var_labels(person_trial)
flatten_scalar_list_cols(person_trial)

notnormal <- intersect(c("age_indx", "time_af_index", "fu"), vars)
tabone <- CreateTableOne(
  vars = vars,
  factorVars = tableone_factor_vars(vars),
  strata = "expo",
  data = prep_for_tableone(person_trial, "expo", vars),
  test = FALSE
)
tableone_b4_matching <- as.data.table(
  print(tabone, smd = TRUE, varLabels = TRUE, nonnormal = notnormal),
  keep.rownames = TRUE
)

dir.create(file.path(root, "out"), showWarnings = FALSE)
write.csv(
  tableone_b4_matching,
  file.path(root, cfg$path_tableone_before_csv),
  row.names = FALSE
)

set.seed(cfg$matchit_seed)
if (cfg$pipeline_style == "run_incident") {
  match_formula <- as.formula(paste(
    "expo ~ Sex + age_indx + dx.chf + dx.mi + dx.pvd + dx.cbd + dx.copd +",
    "dx.dementia + dx.paralysis + dx.dm + dx.crf + dx.liver_mild +",
    "dx.liver_modsev + dx.ulcers + dx.asthma + dx.ra + dx.aids + dx.cancer +",
    "dx.cancer_mets + dx.dm_com0 + dx.dm_com1 + dx.htn + dx.constipation +",
    "time_af_index + chadsvas_score + fu"
  ))
  m.out <- matchit(
    match_formula,
    data = person_trial,
    replace = FALSE,
    method = "nearest",
    ratio = cfg$matchit_ratio,
    exact = ~ trial_id
  )
} else {
  match_formula <- as.formula(paste(
    "expo ~ Sex + age_indx + dx.chf + dx.mi + dx.pvd + dx.cbd + dx.copd +",
    "dx.dementia + dx.paralysis + dx.dm + dx.crf + dx.liver_mild +",
    "dx.liver_modsev + dx.ulcers + dx.stroke_embo + dx.asthma + dx.ra +",
    "dx.aids + dx.cancer + dx.cancer_mets + dx.dm_com0 + dx.dm_com1 +",
    "dx.htn + dx.constipation + time_af_index + chadsvas_score"
  ))
  m.out <- matchit(
    match_formula,
    data = person_trial,
    replace = cfg$matchit_replace,
    method = "nearest",
    caliper = cfg$matchit_caliper,
    ratio = cfg$matchit_ratio
  )
}
print(summary(m.out))

person_trial <- as.data.table(match.data(m.out))
flatten_scalar_list_cols(person_trial)
tabone2 <- CreateTableOne(
  vars = vars,
  factorVars = tableone_factor_vars(vars),
  strata = "expo",
  data = prep_for_tableone(person_trial, "expo", vars),
  test = FALSE
)
tableone_after_matching <- as.data.table(
  print(tabone2, smd = TRUE, varLabels = TRUE, nonnormal = notnormal),
  keep.rownames = TRUE
)

write.csv(
  tableone_after_matching,
  file.path(root, cfg$path_tableone_after_csv),
  row.names = FALSE
)

saveRDS(person_trial, file.path(root, cfg$path_person_trial_rds))
message("Saved: ", cfg$path_person_trial_rds)
message("N rows (matched): ", nrow(person_trial))
message("expo table: ")
print(person_trial[, .N, expo])
