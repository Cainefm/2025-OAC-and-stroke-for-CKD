# 00_functions_pd_oac.R — helpers for sequential PD + AF + OAC emulation
# Sourced by 01_build_seqcohort.R, 02_analysis_primary.R, 03_regression_pooled_logistic.R

#' Working directory must be project root (folder `data/` exists).
project_root <- function() {
  if (!dir.exists("data")) {
    stop("Set working directory to project root (folder data/ missing).")
  }
  "."
}

#' CHA2DS2-VASc score (age bands match prior code: gap at exactly 75 years).
compute_chadsvas_score <- function(dt) {
  dt[, age_indx := suppressWarnings(as.numeric(as.vector(age_indx)))]
  dt[, chadsvas_score := data.table::fcase(
    is.na(age_indx), NA_integer_,
    age_indx < 65, 0L,
    age_indx >= 65 & age_indx < 75, 1L,
    age_indx > 75, 2L,
    default = NA_integer_
  ) + data.table::fifelse(as.character(Sex) == "F", 1L, 0L) +
    as.integer(dx.chf) + as.integer(dx.htn) + as.integer(dx.cbd) +
    as.integer(dx.stroke_embo) * 2L + as.integer(dx.pvd) + as.integer(dx.dm)]
  invisible(dt)
}

#' Factor columns for TableOne (continuous vars excluded).
tableone_factor_vars <- function(vars) {
  setdiff(vars, c("fu", "trial_id", "chadsvas_score", "age_indx", "time_af_index"))
}

#' In-place: turn mistaken list-of-scalars columns into atomic (Excel / labels).
flatten_scalar_list_cols <- function(tab) {
  n <- nrow(tab)
  for (nm in names(tab)) {
    v <- tab[[nm]]
    if (!is.list(v) || inherits(v, "data.frame") || inherits(v, "POSIXlt")) {
      next
    }
    if (length(v) != n) {
      next
    }
    if (!all(lengths(v) == 1L)) {
      stop("Column ", nm, ": list cells must be length 1 (see protocol xlsx).")
    }
    val <- unlist(v, use.names = FALSE)
    if (inherits(tab, "data.table")) {
      tab[, (nm) := val]
    } else {
      tab[[nm]] <- val
    }
  }
  invisible(tab)
}

#' data.frame for CreateTableOne: strip haven labels and list wrappers.
prep_for_tableone <- function(dt, strata, vars) {
  cols <- intersect(unique(c(strata, vars)), names(dt))
  sub <- if (inherits(dt, "data.table")) {
    dt[, cols, with = FALSE]
  } else {
    dt[, cols, drop = FALSE]
  }
  out <- as.data.frame(sub)
  if (requireNamespace("haven", quietly = TRUE)) {
    out <- haven::zap_labels(out)
  } else {
    for (j in seq_len(ncol(out))) {
      if (inherits(out[[j]], "labelled")) {
        out[[j]] <- as.vector(out[[j]])
      }
    }
  }
  flatten_scalar_list_cols(out)
  out
}

#' Apply protocol xlsx Name / Description as variable labels.
apply_protocol_var_labels <- function(person_trial, covar) {
  for (i in seq_len(nrow(covar))) {
    nm <- trimws(as.character(covar[["Name"]][i]))
    desc <- paste(
      as.character(unlist(covar[["Description"]][i], use.names = FALSE)),
      collapse = " "
    )
    if (nzchar(nm) && nm %in% names(person_trial)) {
      var_label(person_trial[[nm]]) <- desc
    }
  }
  invisible(person_trial)
}

#' Override labels for core Table-1 columns.
set_core_var_labels <- function(person_trial) {
  core <- c(
    age_indx = "Age at index date (Years)",
    time_af_index = "Time from AF to Index date (Days)",
    chadsvas_score = "CHA₂DS₂-VASc Score",
    antiplatelet.baseline = "Antiplatelet user history",
    fu = "Follow-up period"
  )
  for (nm in names(core)) {
    if (nm %in% names(person_trial)) {
      var_label(person_trial[[nm]]) <- core[[nm]]
    }
  }
  invisible(person_trial)
}

#' Read seqcohort list and attach `trial_id` per trial.
read_seqcohort_with_trial_ids <- function(path) {
  seqcohort <- readRDS(path)
  lapply(seq_along(seqcohort), function(i) {
    df <- seqcohort[[i]]
    df$trial_id <- i
    df
  })
}

#' Cached monthly antiplatelet at index (same file as 01 saves next to seqcohort).
load_or_build_antiplatelet_monthly <- function(cfg, root) {
  ap_path <- sub("\\.RDS$", "_antiplatelet.RDS", cfg$path_seqcohort_rds)
  ap_full <- file.path(root, ap_path)
  if (file.exists(ap_full)) {
    return(readRDS(ap_full))
  }
  message("Rebuilding antiplatelet from drug RDS (run 01 to cache).")
  rx <- readRDS(file.path(root, cfg$path_drug_rds))
  data.table::setDT(rx)
  prepare_antiplatelet(rx, collapse_to_month = (cfg$index_by == "month"))
}

#' Merge antiplatelet at index date; result column `antiplatelet.baseline` 0/1.
merge_antiplatelet_baseline_at_index <- function(person_trial, drug_ap) {
  out <- merge(
    person_trial,
    drug_ap,
    by.x = c("Reference_Key", "indx_date"),
    by.y = c("Reference_Key", "Date"),
    all.x = TRUE
  )
  data.table::setnames(out, "antiplatelet", "antiplatelet.baseline")
  out[, antiplatelet.baseline := ifelse(
    is.na(antiplatelet.baseline), 0L, as.integer(antiplatelet.baseline)
  )]
  out
}

#' Ensure `antiplatelet.baseline` on wide `person_trial` (for 02 / 03).
ensure_antiplatelet_baseline_wide <- function(person_trial, cfg, root) {
  if ("antiplatelet" %in% names(person_trial) &&
        !"antiplatelet.baseline" %in% names(person_trial)) {
    person_trial[, antiplatelet.baseline := as.integer(antiplatelet)]
    person_trial[, antiplatelet := NULL]
  }
  if (!"antiplatelet.baseline" %in% names(person_trial)) {
    message("Adding antiplatelet.baseline from drug file (merge at indx_date).")
    rx_fix <- readRDS(file.path(root, cfg$path_drug_rds))
    data.table::setDT(rx_fix)
    drug_ap_idx <- prepare_antiplatelet(
      rx_fix,
      collapse_to_month = (cfg$index_by == "month")
    )
    person_trial <- merge(
      person_trial,
      drug_ap_idx,
      by.x = c("Reference_Key", "indx_date"),
      by.y = c("Reference_Key", "Date"),
      all.x = TRUE
    )
    person_trial[, antiplatelet.baseline := ifelse(is.na(antiplatelet), 0L, 1L)]
    if ("antiplatelet" %in% names(person_trial)) {
      person_trial[, antiplatelet := NULL]
    }
  }
  person_trial[, antiplatelet.baseline := ifelse(
    is.na(antiplatelet.baseline), 0L, as.integer(antiplatelet.baseline)
  )]
  person_trial
}

#' Map one TableOne row name (`rn`) to a column name in `data_names`.
#' Strips factor levels and continuous summaries so formulas do not parse
#' `var (mean (SD))` as a call to `var`.
sanitize_tableone_rn <- function(rn, data_names) {
  if (length(rn) != 1L) {
    return(NA_character_)
  }
  s <- as.character(rn)[1L]
  if (is.na(s) || !nzchar(trimws(s))) {
    return(NA_character_)
  }
  s <- gsub(" = [01] \\(%\\)", "", s)
  s <- gsub(" = [MN] \\(%\\)", "", s)
  s <- gsub(" \\(mean \\(SD\\)\\)", "", s)
  s <- gsub(" \\(median.*\\)", "", s, perl = TRUE)
  s <- sub("\\s*\\(.*$", "", s)
  s <- trimws(s)
  if (!nzchar(s)) {
    return(NA_character_)
  }
  if (s %in% data_names) {
    return(s)
  }
  s2 <- gsub(".baseline", "", s, fixed = TRUE)
  if (s2 %in% data_names) {
    return(s2)
  }
  NA_character_
}

#' Vector of `rn` values -> unique valid `data_names` (SMD > threshold picks).
clean_unbalanced_from_tableone_rn <- function(rn_vec, data_names) {
  if (length(rn_vec) == 0L) {
    return(character(0))
  }
  rn_vec <- as.character(rn_vec)
  out <- vapply(
    rn_vec,
    sanitize_tableone_rn,
    character(1L),
    data_names = data_names,
    USE.NAMES = FALSE
  )
  out <- out[!is.na(out) & nzchar(out)]
  unique(out)
}

#' TableOne SMD matrix as `data.table` (console table suppressed).
tableone_smd_to_dt <- function(tab_obj, nonnormal = NULL, var_labels = TRUE) {
  args <- list(
    x = tab_obj,
    smd = TRUE,
    varLabels = var_labels,
    printToggle = FALSE
  )
  if (length(nonnormal) > 0L) {
    args$nonnormal <- nonnormal
  }
  smd_mat <- do.call(print, args)
  as.data.table(
    as.data.frame(smd_mat, stringsAsFactors = FALSE),
    keep.rownames = TRUE
  )
}

#' `seq.Date` step for long person-time in 03 (`"1 day"` or `"1 month"`).
long_person_time_seq_by <- function(cfg) {
  if (!cfg$index_by %in% c("day", "month")) {
    stop('analysis_config index_by must be "day" or "month".')
  }
  if (cfg$index_by == "day") {
    "1 day"
  } else {
    "1 month"
  }
}

#' Horizon length for g-computation in bootstrap (matches `obs_end_cap_years`).
#' Monthly index: 12 steps per year; daily index: ~365.25 days per year.
gcomputation_n_steps <- function(cfg) {
  if (cfg$index_by == "day") {
    as.integer(floor(365.25 * cfg$obs_end_cap_years))
  } else {
    as.integer(cfg$obs_end_cap_years * 12L)
  }
}

#' Build calendar index dates from config
build_index_dates <- function(cfg) {
  by_unit <- if (cfg$index_by == "month") "month" else "day"
  d <- seq.Date(cfg$study_start, cfg$study_end, by = paste0("1 ", by_unit))
  if (isTRUE(cfg$drop_first_index)) {
    d[-1L]
  } else {
    d
  }
}

#' OAC prescription intervals for current-use exposure
prepare_drug_oac <- function(rx, oac_end_impute_days) {
  drug_name_oac <- paste0(
    "WARFARIN|COUMADIN|ACENOCOUMAROL|PHENINDIONE|PHENPROCOUMON|",
    "RIVAROXABAN|XARELTO|",
    "APIXABAN|ELIQUIS|",
    "DABIGATRAN|PRADAXA|",
    "EDOXABAN|LIXIANA|SAVAYSA"
  )
  out <- rx[grepl(drug_name_oac, Drug_Name, ignore.case = TRUE),
            .(Reference_Key,
              Prescription_Start_Date = ymd(Prescription_Start_Date),
              Prescription_End_Date = ymd(Prescription_End_Date))]
  out <- out[!is.na(Prescription_Start_Date)]
  out[is.na(Prescription_End_Date),
      Prescription_End_Date := Prescription_Start_Date + oac_end_impute_days]
  out[Prescription_Start_Date <= Prescription_End_Date]
}

#' Antiplatelet person-days; optional collapse to month start for monthly index
prepare_antiplatelet <- function(rx, collapse_to_month) {
  drug_name_antiplatelet <- paste0(
    "ASPIRIN|CARDIPRIN|CARTIA|PROPIRIN|ASPI-COR|ASPILETS|AGGRENOX|",
    "CLOPIE?DOG?(ER|RE)L|PLAVIX|PRASUGREL|EFFIENT|TICAGREC?LOR|BRILINTA|",
    "TICLOPIDINE|TICLID|CANGRELOR|KENGREX?AL|D?IPYRIDAMOLE|PERSANTIN|",
    "PROCARDIN|AGGRENOX|CILOSTAZOL|PLETAAL|VORAPAXAR|ZONTIVITY|",
    "ABCIXIMAB|REOPRO|EPTIFIBATIDE|INTEGRILIN"
  )
  out <- rx[grepl(drug_name_antiplatelet, Drug_Name, ignore.case = TRUE),
            .(Reference_Key,
              Prescription_Start_Date = ymd(Prescription_Start_Date),
              Prescription_End_Date = ymd(Prescription_End_Date))]
  out <- out[!is.na(Prescription_End_Date)]
  out <- out[Prescription_Start_Date <= Prescription_End_Date]
  # Same-day vs multi-day: grouped seq() can mix Date storage (int/double) in j
  same_day <- out[Prescription_Start_Date == Prescription_End_Date,
                  .(Reference_Key,
                    Date = as.Date(Prescription_Start_Date),
                    antiplatelet = 1L)]
  multi_day <- out[Prescription_Start_Date < Prescription_End_Date,
                   .(Date = as.Date(seq(
                     Prescription_Start_Date[1L],
                     Prescription_End_Date[1L],
                     by = "day"
                   ))),
                   by = .(Reference_Key, Prescription_Start_Date, Prescription_End_Date)]
  multi_day[, antiplatelet := 1L]
  multi_day <- multi_day[, .(Reference_Key, Date, antiplatelet)]
  out <- unique(rbindlist(list(same_day, multi_day), use.names = TRUE))
  if (isTRUE(collapse_to_month)) {
    out[, Date := ymd(paste0(substring(Date, 1L, 7L), "-01"))]
    out <- unique(out)
  }
  out
}

assign_oac_expo_current_use <- function(dt, drug_oac, ind) {
  if (nrow(dt) == 0L) {
    return(dt)
  }
  keys <- unique(dt$Reference_Key)
  sub <- drug_oac[Reference_Key %in% keys &
                  Prescription_Start_Date <= ind &
                  Prescription_End_Date >= ind]
  on_key <- unique(sub$Reference_Key)
  dt[, expo := as.integer(Reference_Key %in% on_key)]
  dt
}

getpx <- function(refkey, regex, data) {
  as.integer(nrow(data[Reference_Key == refkey &
    grepl(regex, All_Diagnosis_Code_ICD9), ]) > 0L)
}

search_outcomes <- function(cohort, each_trial_dx, icd, outcomename) {
  tmp <- merge(
    cohort,
    each_trial_dx[grepl(icd, All_Diagnosis_Code_ICD9),
                  .(Reference_Key, Reference_Date, out = 1L)][,
                  .(Reference_Date = min(Reference_Date), out = 1L),
                  by = Reference_Key],
    by = "Reference_Key",
    all.x = TRUE
  )
  tmp[, out := ifelse(is.na(out), 0L, out)]
  setnames(
    tmp,
    c("Reference_Date", "out"),
    c(paste0("date.", outcomename), paste0("out.", outcomename))
  )
  tmp
}

#' Factory returning run_seq(index_date, cohort) aligned with analysis_config
make_run_seq <- function(cfg, dx, covar, drug_oac) {
  function(index_date, cohort) {
    if (isTRUE(cfg$verbose_trial_msg)) {
      message("index_date: ", index_date)
    }
    each_trial <- data.table::copy(cohort)
    each_trial[, indx_date := index_date]
    each_trial_PD <- each_trial[date.PD <= indx_date & date.PD.end >= indx_date]
    each_trial_PD_AF <- each_trial_PD[date.af <= indx_date]
    each_trial_PD_AF <- each_trial_PD_AF[date.af >= date.PD]

    if (cfg$exposure_mode == "manuscript_initiator") {
      m_prev <- cfg$oac_prevalent_exclusion_months
      id_prev_oac <- each_trial_PD_AF[
        date.oac <= indx_date %m-% months(m_prev), Reference_Key]
      base_af <- each_trial_PD_AF[!Reference_Key %in% id_prev_oac]
    } else if (cfg$exposure_mode == "current_oac_interval") {
      base_af <- each_trial_PD_AF
    } else {
      stop("Unknown exposure_mode: ", cfg$exposure_mode)
    }

    id_prev_out <- base_af[
      date.death <= indx_date | date.transplant <= indx_date |
        date.hd <= indx_date,
      Reference_Key]
    base2 <- base_af[!Reference_Key %in% id_prev_out]
    each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx <- data.table::copy(base2)

    if (cfg$exposure_mode == "manuscript_initiator") {
      m_new <- cfg$oac_initiation_lookback_months
      each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[, expo := as.integer(
        (date.oac > indx_date %m-% months(m_new) & date.oac <= indx_date) &
          !is.na(date.oac))]
    } else {
      if (is.null(drug_oac) || nrow(drug_oac) == 0L) {
        stop("current_oac_interval requires non-empty drug_oac from prepare_drug_oac()")
      }
      each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx <-
        assign_oac_expo_current_use(
          each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx,
          drug_oac,
          index_date
        )
    }

    each_trial_dx <- dx[
      Reference_Key %in%
        each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key &
        Reference_Date < index_date]

    px <- mapply(
      function(x) {
        sapply(
          each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key,
          function(y) getpx(y, x, each_trial_dx)
        )
      },
      as.list(tibble::deframe(covar[, .(Name, Regex)])),
      USE.NAMES = TRUE
    )

    each_trial_px <- merge(
      each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx,
      data.table::as.data.table(px, keep.rownames = TRUE)[,
        setnames(.SD, "rn", "Reference_Key")],
      by = "Reference_Key"
    )

    if (isTRUE(cfg$cohort_exclude_dx_stroke_embo)) {
      each_trial_px <- each_trial_px[dx.stroke_embo == 0L]
    }

    each_trial_dx <- dx[
      Reference_Key %in%
        each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key &
        Reference_Date >= index_date]

    each_trial_px <- search_outcomes(
      each_trial_px, each_trial_dx, "^43[0-6]|^444", "stroke")
    icd_is <- "^43[3-6]|^444"
    each_trial_px <- search_outcomes(
      each_trial_px, each_trial_dx, icd_is, "ische")
    icd_hs <- "^43[0-2]"
    each_trial_px <- search_outcomes(
      each_trial_px, each_trial_dx, icd_hs, "haem")
    icd_gi <- paste0(
      "^456.0^|456.20|^530.21|^530.7|^530.82|^531.[0246]|^532.[0246]|",
      "^533.[0246]|^534.[0246]|^537.8[34]|^562.0[23]|^562.1[23]|",
      "^569.3|^569.85|^578")
    each_trial_px <- search_outcomes(
      each_trial_px, each_trial_dx, icd_gi, "GIbleed")
    icd_ob <- paste0(
      "^287.[89]|^596.7|^784.8|^599.7|^627.1|^459.0|^719.1|^786.3|",
      "^363.6|^376.32|^377.42|^729.92|^423.0|^801.[2378]|^803.[2378]|",
      "^804.[2378]|^800.[2378]|^85[23]")
    each_trial_px <- search_outcomes(
      each_trial_px, each_trial_dx, icd_ob, "Otherbleed")

    each_trial_px
  }
}
