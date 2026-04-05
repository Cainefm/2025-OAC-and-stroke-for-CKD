# 03_regression_pooled_logistic.R — monthly person-time + pooled logistic + bootstrap
# Run after 02_analysis_primary.R.
#
# pipeline_style "run_incident": long-format SMD, speedglm coefs (optional), and/or
# bootstrap RD/RR like script/run_incident.R (~700–772).
#
# Run: Rscript script/03_regression_pooled_logistic.R
#
# Config (analysis_config.R):
#   regression_run_speedglm — SMD-adjusted speedglm coef tables (optional).
#   regression_run_bootstrap — cluster bootstrap RD/RR (run_incident-style; optional).
# At least one must be TRUE.

options(scipen = 6L, digits = 4L)
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(lubridate)
  library(tableone)
  library(broom)
  library(boot)
})

source("script/00_functions_pd_oac.R")
root <- project_root()

source(file.path(root, "script/analysis_config.R"))
cfg <- load_analysis_config(root = root)
run_sg <- isTRUE(cfg$regression_run_speedglm)
run_boot <- isTRUE(cfg$regression_run_bootstrap)
if (!run_sg && !run_boot) {
  stop(
    "Set regression_run_speedglm = TRUE and/or regression_run_bootstrap = TRUE ",
    "in script/analysis_config.R (both FALSE would run no models)."
  )
}
if (run_sg || run_boot) {
  suppressPackageStartupMessages(library(speedglm))
}
if (run_boot) {
  suppressPackageStartupMessages(library(splitstackshape))
}
message("pipeline_style: ", cfg$pipeline_style)
message(
  "regression_run_speedglm: ", run_sg,
  " | regression_run_bootstrap: ", run_boot
)

path_pt <- file.path(root, cfg$path_person_trial_rds)
if (!file.exists(path_pt)) {
  stop("Run script/02_analysis_primary.R first; missing: ", path_pt)
}

person_trial <- as.data.table(readRDS(path_pt))
message(
  "Matched person-trial rows (Word Table 1 / step 02): n = ",
  nrow(person_trial), " | ",
  paste(person_trial[, .N, expo][, paste0("expo=", expo, ": ", N)], collapse = "; ")
)

person_trial <- ensure_antiplatelet_baseline_wide(person_trial, cfg, root)

# --- Wide-format TableOne for SMD (default pipeline only) -------------------
if (cfg$pipeline_style != "run_incident" && run_sg) {
  vars_tab <- c(
    "fu", "Sex", "age_indx", "time_af_index", "chadsvas_score",
    "dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
    "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
    "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
    "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
    "dx.dm_com1", "dx.htn", "dx.constipation", "antiplatelet.baseline"
  )
  vars_tab <- intersect(vars_tab, names(person_trial))
  if (!"antiplatelet.baseline" %in% vars_tab) {
    message("Warning: antiplatelet.baseline not in person_trial; Table 1 omits it.")
  }

  tabone <- CreateTableOne(
    vars = vars_tab,
    factorVars = intersect(
      tableone_factor_vars(vars_tab),
      names(person_trial)
    ),
    strata = "expo",
    data = person_trial,
    test = FALSE
  )
  notnormal_tab <- intersect(c("age_indx", "time_af_index", "fu"), vars_tab)
  smd_lines <- tryCatch(
    tableone_smd_to_dt(tabone, nonnormal = notnormal_tab, var_labels = FALSE),
    error = function(e) NULL
  )
  if (is.null(smd_lines) || !"SMD" %in% names(smd_lines)) {
    message("Could not parse TableOne SMD; using expo + time + timesqr only.")
    unbalanced_variable <- character(0)
  } else {
    raw_rn <- smd_lines[suppressWarnings(as.numeric(SMD)) > 0.1, rn]
    unbalanced_variable <- clean_unbalanced_from_tableone_rn(
      raw_rn,
      names(person_trial)
    )
    unbalanced_variable <- setdiff(unbalanced_variable, c("trial_id", ""))
    if (length(unbalanced_variable) == 0L) {
      message("No covariates with SMD > 0.1; using expo + time + timesqr only.")
      unbalanced_variable <- character(0)
    }
  }

  person_trial[, antiplatelet_baseline := as.integer(antiplatelet.baseline)]
  unbalanced_variable <- gsub(
    "antiplatelet.baseline",
    "antiplatelet_baseline",
    unbalanced_variable,
    fixed = TRUE
  )
} else if (cfg$pipeline_style != "run_incident" && !run_sg) {
  unbalanced_variable <- character(0)
}

# --- Monthly person-time (same structure as script/run_incident.R) ------------
long_format <- person_trial[, .(id, indx_date, obs_end)][,
  .(obs_date = seq.Date(indx_date, obs_end, by = "1 month")),
  by = id]
long_format[, time := seq_len(.N) - 1L, by = id]
long_format[, timesqr := time^2]
person_trial <- merge(person_trial, long_format, by = "id", all.x = TRUE)

if (!"dob" %in% names(person_trial)) {
  stop("Need dob on person_trial (from PD_OAC-demo merge in 02).")
}
person_trial[, age_indx := as.numeric((obs_date - ymd(dob) + 1L) / 365.25)]
compute_chadsvas_score(person_trial)
person_trial[, time_af_index := interval(date.af, obs_date) %/% months(1L)]
person_trial[, fu := as.numeric(obs_date - indx_date)]

message("outcome_definition = ", cfg$outcome_definition)
if (!all(c("date.ische", "date.haem", "date.death") %in% names(person_trial))) {
  stop("person_trial needs date.ische, date.haem, date.death from seqcohort.")
}

if (cfg$outcome_definition == "pmin_composite") {
  if (cfg$pipeline_style == "run_incident") {
    person_trial[, out := ifelse(
      pmin(date.ische, date.haem, date.death, na.rm = TRUE) >= obs_date &
        substring(
          pmin(date.ische, date.haem, date.death, na.rm = TRUE), 1L, 7L
        ) == substring(obs_date, 1L, 7L),
      1, 0
    )]
    person_trial[, out := ifelse(is.na(out), 0, out)]
    person_trial[, out := as.integer(out)]
  } else {
    person_trial[, t_evt := pmin(date.ische, date.haem, date.death, na.rm = TRUE)]
    person_trial[, out := fifelse(
      is.na(t_evt) | !is.finite(t_evt),
      0L,
      as.integer(
        t_evt >= obs_date &
          substring(t_evt, 1L, 7L) == substring(obs_date, 1L, 7L)
      )
    )]
    person_trial[, t_evt := NULL]
  }
} else if (cfg$outcome_definition == "any_event_same_month") {
  person_trial[, out := 0L]
  person_trial[, out := fifelse(
    !is.na(date.ische) &
      format(date.ische, "%Y-%m") == format(obs_date, "%Y-%m") &
      date.ische >= obs_date, 1L, out)]
  person_trial[, out := fifelse(
    !is.na(date.haem) &
      format(date.haem, "%Y-%m") == format(obs_date, "%Y-%m") &
      date.haem >= obs_date, 1L, out)]
  person_trial[, out := fifelse(
    !is.na(date.death) &
      format(date.death, "%Y-%m") == format(obs_date, "%Y-%m") &
      date.death >= obs_date, 1L, out)]
} else {
  stop("Unknown outcome_definition: ", cfg$outcome_definition)
}

rx <- readRDS(file.path(root, cfg$path_drug_rds))
setDT(rx)
drug_ap <- prepare_antiplatelet(rx, collapse_to_month = FALSE)
person_trial <- merge(
  person_trial,
  drug_ap,
  by.x = c("Reference_Key", "obs_date"),
  by.y = c("Reference_Key", "Date"),
  all.x = TRUE
)
person_trial[, antiplatelet := ifelse(is.na(antiplatelet), 0L, 1L)]

# --- SMD / unbalanced (run_incident: after tv antiplatelet) -------------------
if (cfg$pipeline_style == "run_incident" && (run_sg || run_boot)) {
  message(
    "SMD for pooled models uses long-format data. Person-month rows: ",
    nrow(person_trial), ". See ", cfg$path_tableone_long_smd_csv
  )
  ut <- person_trial[, .(n_unique_ids = uniqueN(id)), by = expo]
  message(
    "Unique participants (id) after expansion: ",
    paste(ut[, paste0("expo=", expo, ": ", n_unique_ids)], collapse = "; ")
  )
  vars_long <- c(
    "fu", "Sex", "age_indx", "time_af_index", "chadsvas_score",
    "dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
    "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
    "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
    "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
    "dx.dm_com1", "dx.htn", "dx.constipation", "antiplatelet"
  )
  vars_long <- intersect(vars_long, names(person_trial))
  notnormal_long <- intersect(c("age_indx", "time_af_index", "fu"), vars_long)
  tabone_long <- CreateTableOne(
    vars = vars_long,
    factorVars = tableone_factor_vars(vars_long),
    strata = "expo",
    data = person_trial,
    test = FALSE
  )
  smd_df <- tableone_smd_to_dt(
    tabone_long,
    nonnormal = notnormal_long,
    var_labels = TRUE
  )
  if (run_sg) {
    write.csv(
      smd_df,
      file.path(root, cfg$path_tableone_long_smd_csv),
      row.names = FALSE
    )
    message("Saved long-format TableOne (SMD): ", cfg$path_tableone_long_smd_csv)
  }
  raw_rn <- smd_df[suppressWarnings(as.numeric(SMD)) > 0.1, rn]
  unbalanced_variable <- clean_unbalanced_from_tableone_rn(
    raw_rn,
    names(person_trial)
  )
  unbalanced_variable <- setdiff(unbalanced_variable, c("trial_id", ""))
  if (length(unbalanced_variable) == 0L) {
    message("No covariates with SMD > 0.1; using expo + time + timesqr only.")
  }
} else if (cfg$pipeline_style != "run_incident") {
  unbalanced_variable <- intersect(unbalanced_variable, names(person_trial))
}

# --- Formulas (run_incident: speedglm and/or bootstrap) ----------------------
if (cfg$pipeline_style == "run_incident" && (run_sg || run_boot)) {
  ri_adj <- function(excl = character(0)) {
    u <- setdiff(
      gsub(".baseline", "", unbalanced_variable, fixed = TRUE),
      c("fu", excl)
    )
    u <- intersect(u, names(person_trial))
    if (length(u) == 0L) {
      return("")
    }
    paste0("+", paste(u, collapse = "+"))
  }
  formula_rr_main <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr)", ri_adj()
  ))
  formula_rr_sex <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr)", ri_adj("Sex")
  ))
  formula_rr_age <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr)", ri_adj("age_indx")
  ))
  formula_rr_inci <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr)", ri_adj("dx.stroke_embo")
  ))
  formula_rr_sex_inci <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr)",
    ri_adj(c("Sex", "dx.stroke_embo"))
  ))

  cov_no_fu <- paste0(
    setdiff(gsub(".baseline", "", unbalanced_variable, fixed = TRUE), "fu"),
    collapse = "+"
  )
  cov_no_sex_fu <- paste0(
    setdiff(
      gsub(".baseline", "", unbalanced_variable, fixed = TRUE),
      c("Sex", "fu")
    ),
    collapse = "+"
  )
  cov_no_age_fu <- paste0(
    setdiff(
      gsub(".baseline", "", unbalanced_variable, fixed = TRUE),
      c("age_indx", "fu")
    ),
    collapse = "+"
  )
  formula_rr_main_cum_plot <- as.formula(paste0(
    "out~expo+expo*Sex+I(time)+I(timesqr)+I(time*expo)+I(timesqr*expo)",
    if (nzchar(cov_no_fu)) paste0("+", cov_no_fu) else ""
  ))
  formula_rr_sex_cum_plot <- as.formula(paste0(
    "out~expo+I(time)+I(timesqr)+I(time*expo)+I(timesqr*expo)",
    if (nzchar(cov_no_sex_fu)) paste0("+", cov_no_sex_fu) else ""
  ))
  formula_rr_age_cum_plot <- as.formula(paste0(
    "out~expo+I(time)+I(timesqr)+I(time*expo)+I(timesqr*expo)",
    if (nzchar(cov_no_age_fu)) paste0("+", cov_no_age_fu) else ""
  ))
} else if (run_sg) {
  rhs_extra <- if (length(unbalanced_variable) > 0L) {
    paste0("+", paste(unbalanced_variable, collapse = "+"))
  } else {
    ""
  }
  formula_rr_main <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr) ", rhs_extra
  ))
  formula_rr_sex <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr) ",
    if (length(setdiff(unbalanced_variable, "Sex")) > 0L) {
      paste0("+", paste(setdiff(unbalanced_variable, "Sex"), collapse = "+"))
    } else {
      ""
    }
  ))
  formula_rr_inci <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr) ",
    if (length(setdiff(unbalanced_variable, "dx.stroke_embo")) > 0L) {
      paste0(
        "+",
        paste(setdiff(unbalanced_variable, "dx.stroke_embo"), collapse = "+")
      )
    } else {
      ""
    }
  ))
  formula_rr_sex_inci <- as.formula(paste0(
    "out ~ expo + I(time) + I(timesqr) ",
    if (length(setdiff(unbalanced_variable, c("Sex", "dx.stroke_embo"))) > 0L) {
      paste0(
        "+",
        paste(
          setdiff(unbalanced_variable, c("Sex", "dx.stroke_embo")),
          collapse = "+"
        )
      )
    } else {
      ""
    }
  ))
}

fit_one <- function(f, dat, label) {
  fit <- tryCatch(
    speedglm(f, data = dat, family = binomial(link = "logit")),
    error = function(e) {
      warning("Model failed (", label, "): ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) {
    return(NULL)
  }
  out <- as.data.table(broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE))
  out[, model := label]
  out
}

if (run_sg) {
  model.primary <- fit_one(formula_rr_main, person_trial, "pri")
  model.F <- fit_one(formula_rr_sex, person_trial[as.character(Sex) == "F"], "F")
  model.M <- fit_one(formula_rr_sex, person_trial[as.character(Sex) == "M"], "M")
  model.inci <- fit_one(formula_rr_inci, person_trial[dx.stroke_embo == 0L], "inci")
  model.F.inci <- fit_one(
    formula_rr_sex_inci,
    person_trial[as.character(Sex) == "F" & dx.stroke_embo == 0L],
    "F_inc"
  )
  model.M.inci <- fit_one(
    formula_rr_sex_inci,
    person_trial[as.character(Sex) == "M" & dx.stroke_embo == 0L],
    "M_inc"
  )
  res_list <- Filter(Negate(is.null), list(
    model.primary, model.F, model.M, model.inci, model.F.inci, model.M.inci
  ))
  all_pool <- rbindlist(res_list, fill = TRUE)
  est_expo <- all_pool[term == "expo"]
} else {
  all_pool <- NULL
  est_expo <- NULL
}

# --- Bootstrap RD/RR (script/run_incident.R ~700–772) -------------------------
year <- cfg$obs_end_cap_years
K <- year * 12L
R_boot <- cfg$bootstrap_r
out_boot <- NULL

if (cfg$pipeline_style == "run_incident" && run_boot) {
  source(file.path(root, "script/bootstrap_pooled_run_incident.R"))
  set.seed(cfg$bootstrap_seed)

  run_one_boot <- function(elig_ids, label_rr, label_rd, model_type) {
    formula_rr <<- label_rr
    formula_rd <<- label_rd
    current_bootstrap <<- 0L
    pb <<- txtProgressBar(min = 0, max = R_boot, style = 3L)
    res <- boot::boot(data = elig_ids, statistic = std.boot, R = R_boot)
    close(pb)
    orgnize_ci(res, model_type, elig_ids)
  }

  elig_main <- data.table(id = unique(person_trial$id))
  main_risk_results_o <- run_one_boot(
    elig_main, formula_rr_main, formula_rr_main_cum_plot, "main"
  )

  elig_m <- data.table(id = unique(person_trial[as.character(Sex) == "M", id]))
  male_risk_results_o <- run_one_boot(
    elig_m, formula_rr_sex, formula_rr_sex_cum_plot, "male"
  )

  elig_f <- data.table(id = unique(person_trial[as.character(Sex) == "F", id]))
  female_risk_results_o <- run_one_boot(
    elig_f, formula_rr_sex, formula_rr_sex_cum_plot, "female"
  )

  elig_y <- data.table(id = unique(person_trial[age_indx < 65, id]))
  younger_risk_results_o <- run_one_boot(
    elig_y, formula_rr_age, formula_rr_age_cum_plot, "younger"
  )

  elig_e <- data.table(id = unique(person_trial[age_indx >= 65, id]))
  elder_risk_results_o <- run_one_boot(
    elig_e, formula_rr_age, formula_rr_age_cum_plot, "elder"
  )

  out_boot <- as.data.frame(rbind(
    main_risk_results_o,
    male_risk_results_o,
    female_risk_results_o,
    younger_risk_results_o,
    elder_risk_results_o
  ))
  write.csv(
    out_boot,
    file.path(root, cfg$path_bootstrap_risk_csv),
    row.names = FALSE
  )
  message("Saved bootstrap RD/RR (run_incident-style): ", cfg$path_bootstrap_risk_csv)
} else if (run_boot) {
  message(
    "Skipping bootstrap (only implemented for pipeline_style = \"run_incident\")."
  )
}

dir.create(file.path(root, "out"), showWarnings = FALSE)
if (run_sg) {
  write.csv(est_expo, file.path(root, cfg$path_est_pooled_expo_csv), row.names = FALSE)
  message("Saved (speedglm expo): ", cfg$path_est_pooled_expo_csv)
  write.csv(all_pool, file.path(root, cfg$path_est_pooled_all_csv), row.names = FALSE)
  message("Saved (speedglm all terms): ", cfg$path_est_pooled_all_csv)
}

message("Rows (person-months): ", nrow(person_trial))
if (run_sg) {
  message("Exposure OR — speedglm (SMD-adjusted formulas):\n")
  print(est_expo)
}
if (!is.null(out_boot)) {
  message("Bootstrap RD/RR (g-computation style per run_incident.R):\n")
  print(out_boot)
}
