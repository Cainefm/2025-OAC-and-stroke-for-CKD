# 02_analysis_primary.R — person_trial, Table 1, MatchIt (matches core of script/run.R)
# Run from project root after 01_build_seqcohort.R:
#   Rscript script/02_analysis_primary.R
# Uses script/analysis_config.R for paths and follow-up cap.

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

root <- if (dir.exists("data")) "." else stop("Set working directory to project root.")

source(file.path(root, "script/analysis_config.R"))
cfg <- load_analysis_config(root = root)
source(file.path(root, "script/00_functions_pd_oac.R"))

seqcohort <- readRDS(file.path(root, cfg$path_seqcohort_rds))
seqcohort <- lapply(seq_along(seqcohort), function(i) {
  df <- seqcohort[[i]]
  df$trial_id <- i
  df
})

person_trial <- rbindlist(seqcohort)
person_trial[, id := paste(trial_id, Reference_Key, sep = "-")]

demo <- as.data.table(readRDS(file.path(root, cfg$path_demo_rds)))
person_trial <- merge(
  person_trial,
  unique(demo[, .(Reference_Key, Sex, dob = Date_of_Birth_yyyymmdd)]),
  by = "Reference_Key"
)
setkey(person_trial, NULL)
person_trial[, id := paste(trial_id, Reference_Key, sep = "_")]

cap_years <- cfg$obs_end_cap_years
person_trial[, obs_end := pmin(
  ymd("2022-12-31"),
  date.stroke,
  date.death,
  date.hd,
  date.transplant,
  indx_date %m+% years(cap_years),
  date.PD.end,
  na.rm = TRUE
)]

person_trial[, age_indx := as.numeric((indx_date - ymd(dob)) / 365.25)]
person_trial[, time_af_index := as.numeric(indx_date - date.af)]
person_trial[, chadsvas_score := dplyr::case_when(
  age_indx < 65 ~ 0,
  age_indx >= 65 & age_indx < 75 ~ 1,
  age_indx > 75 ~ 2
) + ifelse(Sex == "F", 1L, 0L) + dx.chf + dx.htn + dx.cbd +
  dx.stroke_embo * 2L + dx.pvd + dx.dm]

person_trial[, fu := as.numeric(obs_end - indx_date)]

ap_path <- sub("\\.RDS$", "_antiplatelet.RDS", cfg$path_seqcohort_rds)
if (!file.exists(file.path(root, ap_path))) {
  message("Rebuilding antiplatelet from drug RDS (run 01 to cache).")
  rx <- readRDS(file.path(root, cfg$path_drug_rds))
  setDT(rx)
  collapse_ap <- (cfg$index_by == "month")
  drug_antiplatelet <- prepare_antiplatelet(rx, collapse_to_month = collapse_ap)
} else {
  drug_antiplatelet <- readRDS(file.path(root, ap_path))
}

person_trial <- merge(
  person_trial,
  drug_antiplatelet,
  by.x = c("Reference_Key", "indx_date"),
  by.y = c("Reference_Key", "Date"),
  all.x = TRUE
)
setnames(person_trial, "antiplatelet", "antiplatelet.baseline")
person_trial[, antiplatelet.baseline := ifelse(
  is.na(antiplatelet.baseline), 0L, antiplatelet.baseline)]

vars <- c(
  "fu", "Sex", "age_indx", "time_af_index", "chadsvas_score",
  "dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
  "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
  "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
  "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
  "dx.dm_com1", "dx.htn", "dx.constipation", "antiplatelet.baseline"
)

covar <- setDT(readxl::read_xlsx(
  file.path(root, cfg$path_protocol_xlsx),
  sheet = cfg$protocol_covar_sheet
))
for (i in seq_len(nrow(covar))) {
  vn <- covar[i, Name]
  vl <- covar[i, Description]
  if (vn %in% names(person_trial)) {
    var_label(person_trial[[vn]]) <- vl
  }
}
var_label(person_trial[["age_indx"]]) <- "Age at index date (Years)"
var_label(person_trial[["time_af_index"]]) <- "Time from AF to Index date (Days)"
var_label(person_trial[["chadsvas_score"]]) <- "CHA₂DS₂-VASc Score"
var_label(person_trial[["antiplatelet.baseline"]]) <- "Antiplatelet user history"
var_label(person_trial[["fu"]]) <- "Follow-up period"
var_label(person_trial)

notnormal <- c("age_indx", "time_af_index", "fu")
tabone <- CreateTableOne(
  vars = vars,
  factorVars = setdiff(vars, c(
    "fu", "trial_id", "chadsvas_score",
    "age_indx", "time_af_index"
  )),
  strata = "expo",
  data = person_trial,
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
print(summary(m.out))

person_trial <- match.data(m.out)
tabone2 <- CreateTableOne(
  vars = vars,
  factorVars = setdiff(vars, c(
    "fu", "trial_id", "chadsvas_score",
    "age_indx", "time_af_index"
  )),
  strata = "expo",
  data = person_trial,
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
