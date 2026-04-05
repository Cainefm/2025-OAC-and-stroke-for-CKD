# 00_functions_pd_oac.R — helpers for sequential PD + AF + OAC emulation
# Sourced by 01_build_seqcohort.R and 02_analysis_primary.R

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
  out <- out[, .(
    Date = seq(Prescription_Start_Date, Prescription_End_Date, by = "day")
  ), by = .(Reference_Key, Prescription_Start_Date, Prescription_End_Date)]
  out <- unique(out[, .(Reference_Key, Date, antiplatelet = 1L)])
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
      id_prev_oac <- each_trial_PD_AF[
        date.oac <= indx_date %m-% months(1L), Reference_Key]
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
      each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[, expo := as.integer(
        (date.oac > indx_date %m-% months(1L) & date.oac <= indx_date) &
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

    each_trial_px <- each_trial_px[dx.stroke_embo == 0L]

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
