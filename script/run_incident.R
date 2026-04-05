## Script name:  run.R
## Purpose of script: crosscheck for Franco about OAC and PD patients
## Author: FAN Min
## Date Created: 2024-11-14
## Copyright (c) FAN Min, 2024
## Email: minfan@connect.hku.hk
##
## Notes:
##   OAC (oral anticoagulation substantially reduce teh risk of ischemic stroke and embolism in patients
##      with atrial fibrillation (AF)
##   OAC vs placebo,
##   outcome stroke and bleeding
##   among PD patients
## Ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC7538212/
##      https://academic-oup-com.eproxy.lib.hku.hk/ckj/article/17/10/sfae270/7811318


# Package loading and environment setting ----------------------------------
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
options(scipen = 6, digits = 4)
memory.limit(30000000)
library(renv)
renv::init()
library(dplyr)
library(data.table)
library(lubridate)
library(pbapply)
library(survival)
library(broom)
library(tableone)
library(survey)
library(ggplot2)
library(MatchIt)
library(labelled)
library(speedglm)
library(splitstackshape)
library(boot)
library(data.table)

# functions

save_plot <-function(filename, plot, width_in, height_in, dpi =300) {
  if(tools::file_ext(filename) =="pdf") {
    ggsave(filename, plot = plot, width = width_in, height = height_in, units ="in")
  }else{
    ggsave(filename, plot = plot, width = width_in * dpi, height = height_in * dpi, units ="px", dpi = dpi)
  }
}


std.boot <- function(data, indices) {
  #* ---- Progress Tracking Optimization ----
  current_bootstrap <<- current_bootstrap + 1
  setTxtProgressBar(pb, current_bootstrap)
  
  #* ---- Bootstrap Sampling Optimization ----
  boot_ids <- data.table(id = data$id[indices])[, bid := .I]  # Generate indexed DT
  
  #* ---- Data Merging Optimization ----
  d <- merge(boot_ids, person_trial,
             by = "id",
             all.x = TRUE,
             allow.cartesian = TRUE)[, bid_new := .GRP, by = .(bid, trial_id)]  # Replace interaction
  
  d <- left_join(boot_ids, person_trial, by="id",relationship="many-to-many") 
  d$bid_new <- interaction(d$bid, d$trial_id) 
  #* ---- Risk Calculation Optimization ----
  risk_diff <- calculate_risk(d, formula_rd, "diff")
  risk_ratio <- calculate_risk(d, formula_rr, "ratio")
  if(all(is.na(risk_diff)) | all(is.na(risk_ratio))) return(c(NA,NA,NA,NA,NA))
  return(c(
    risk_diff[time_0==K, risk0],    # Control group risk
    risk_diff[time_0==K, risk1],    # Treatment group risk
    risk_diff[time_0==K, rd],       # Risk difference
    risk_ratio[time_0==K, rr],       # Risk ratio
    risk_diff[time_0==K, rr]       # Risk ratio from rd model
  ))
}

#* ---- Risk Calculation Helper Function ----
calculate_risk <- function(dt, formula, type) {
  #! Model Fitting Optimization
  model <- tryCatch({
    speedglm(formula, family = binomial(link = "logit"), data = dt)
  }, error = function(e) {
    message("Model failed at bootstrap ", current_bootstrap, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(model)) return(list(risk0=NA, risk1=NA, rd=NA, rr=NA))
  print(model)
  # Create dataset with all time points for each individual under each treatment level
  
  #* ---- Counterfactual Prediction ----
  base_data <- expandRows(dt[time==0,], count=year*12, count.is.col=F)[, 
                                                                       time := rep(0:(year*12-1), dt[time==0,.N])][, 
                                                                                                                   timesqr := time^2]
  # Control scenario prediction
  control <- copy(base_data)[, expo := 0][, 
                                          p_event := predict(model, .SD, type = "response")]
  
  # Treatment scenario prediction
  treatment <- copy(base_data)[, expo := 1][, 
                                            p_event := predict(model, .SD, type = "response")]
  
  #* ---- Survival Analysis ----
  process_survival <- function(dt) {
    dt[, surv := cumprod(1 - p_event), by = bid_new][, 
                                                     risk := 1 - surv]
  }
  
  control_risk <- process_survival(control)
  treatment_risk <- process_survival(treatment)
  
  #* ---- Risk Aggregation ----
  aggregate_risk <- function(dt) {
    dt[, .(risk = mean(risk)), by = .(time, expo)]
  }
  # View(control_risk[id=="100_1299894"])
  # View(treatment_risk[id=="155_100617"])
  
  graph_pred <- merge(
    aggregate_risk(control_risk)[, .(time, risk0 = risk)],
    aggregate_risk(treatment_risk)[, .(time, risk1 = risk)],
    by = "time"
  )[, `:=`(time_0 = time + 1,
           rd = risk1 - risk0,
           rr = risk1 / risk0)]
  
  #! Initialize time=0 data
  zero_row <- data.table(time=0, risk0=0, risk1=0, rd=0, rr=0, time_0=0)
  rbindlist(list(zero_row, graph_pred), fill=TRUE)
}

#* ---- Results Processing ----
extract_ci <- function(results, index) {
  boot.ci(results, conf = 0.95, type = "perc", index = index)$percent[4:5]
}

#* ---- orgnize ci  ----
orgnize_ci <- function(x,model_type,ids){
  est <- dcast(data.table(
    arm = c("risk_nonuser", "risk_user", "risk_diff", "risk_ratio","risk_ratio_rdmodel"),
    est = round(x$t0, 4),
    lower = sapply(1:5, function(i) round(extract_ci(x, i)[1],4)),
    upper = sapply(1:5, function(i) round(extract_ci(x, i)[2],4))
  )[,estimation:=paste0(round(est,2)," (",round(lower,2),", ",round(upper,2),")")][,.(arm,estimation)],.~arm)[,.(risk_nonuser,risk_user,risk_diff,risk_ratio,risk_ratio_rdmodel)]
  
  # return(person_trial[id %in% id,.N,expo][,num.of.p:="person-months"][])
  # return(person_trial[id %in% id,][,unique(.SD),.SDcols=c("Reference_Key","expo")][,.N,expo][,num.of.p:="individuals"][])
  # return(rbind(
    # person_trial[id %in% id,.N,expo][,num.of.p:="person-months"],
    # person_trial[id %in% id,][,unique(.SD),.SDcols=c("Reference_Key","expo")][,.N,expo][,num.of.p:="individuals"]))
  num <- transpose(
    as.data.table(rbind(
      as.data.frame(person_trial[id %in% ids$id,.N,expo][,num.of.p:="person-months"]),
      as.data.frame(person_trial[id %in% ids$id,][,unique(.SD),.SDcols=c("Reference_Key","expo")][,.N,expo][,num.of.p:="individuals"])))[,.(num.of.p=paste0(num.of.p," ",ifelse(expo==0,"non-users","users")),N)],make.names = "num.of.p" )
  return(cbind(model_type,num,est))
}


# -------------------------------------------------------------------------

# future feature define exposure / exclusion
# any AF after during PD
#
# 2009-11-01 / first index date /
#   2009-12-01 / past hx of AF-> eligible /
#
#   4:
#   148677 2022-05-04 2022-06-17 2022-06-20 WARF02
#
# 2022-06-01 / Not past hx of AF / not eligible
# 2022-07-01 / eligible past hx of PD/ past hx of AF / oac after AF time window? exposure group
# 2022-08-01 / eligible
# 2010-01-01 / past hx of AF/

# to do
# af should happened within the PD period: fix
# oac should happened after the AF (noAF): fix
# censoring: death, transplant, HD, PD end (focusing on PD patients only)

# Generate study start to end date by sequence (180 months):
index_date <- seq.Date(ymd("2007-1-1"),ymd("2022-12-31"),by = "1 month")

# Load full cohort ------------------------------------------------------
cohort <- setDT(readRDS("data/pd_oac-cohort-20240813.RDS")) # 1632 patients/ records

ymd_cols <- grep("date",colnames(cohort),value=T) # change to date format
cohort[,(ymd_cols):=lapply(.SD, ymd),.SDcols=ymd_cols]

cohort[,date.PD.end:=pmin(date.transplant,date.PD.end,date.hd,na.rm = T)]

# Load covariate and dx data ----------------------------------------------
covar <- setDT(readxl::read_xlsx("documents/PD OAC-Protocol.xlsx",sheet = "Dx.Cov"))
icd_ischemicstroke <- "^43[3-6]|^444"
icd_HaemorrhagicStroke <- "^43[0-2]"
icd_GI_Bleed <-  "^456.0^|456.20|^530.21|^530.7|^530.82|^531.[0246]|^532.[0246]|^533.[0246]|^534.[0246]|^537.8[34]|^562.0[23]|^562.1[23]|^569.3|^569.85|^578"
icd_Other_Bleed <- "^287.[89]|^596.7|^784.8|^599.7|^627.1|^459.0|^719.1|^786.3|^363.6|^376.32|^377.42|^729.92|^423.0|^801.[2378]|^803.[2378]|^804.[2378]|^800.[2378]|^85[23]"
icd_MI <- "^410"
icd_fracture <- "^8[012][0-9]"

dx <- setDT(readRDS("data/PD_OAC-dx.RDS"))
dx[,Reference_Date:=ymd(Reference_Date)]
all_needed_icd <- paste0(covar[,paste(Regex,collapse = "|")],
                         icd_ischemicstroke,
                         icd_HaemorrhagicStroke,
                         icd_GI_Bleed,
                         icd_Other_Bleed,
                         icd_MI,
                         icd_fracture,
                         collapse = "|")
dx <- dx[grepl(all_needed_icd,All_Diagnosis_Code_ICD9)]
dx <- dx[,.(Reference_Key,Reference_Date,All_Diagnosis_Code_ICD9)][,.SD[Reference_Date==min(Reference_Date)],.(Reference_Key,All_Diagnosis_Code_ICD9)]

rx <- readRDS("data/PD_OAC-drug.RDS")
setDT(rx)

drug_name_antiplatelet <- "ASPIRIN|CARDIPRIN|CARTIA|PROPIRIN|ASPI-COR|ASPILETS|AGGRENOX|CLOPIE?DOG?(ER|RE)L|PLAVIX|PRASUGREL|EFFIENT|TICAGREC?LOR|BRILINTA|TICLOPIDINE|TICLID|CANGRELOR|KENGREX?AL|D?IPYRIDAMOLE|PERSANTIN|PROCARDIN|AGGRENOX|CILOSTAZOL|PLETAAL|VORAPAXAR|ZONTIVITY|ABCIXIMAB|REOPRO|EPTIFIBATIDE|INTEGRILIN"
drug_antiplatelet <- rx[grepl(drug_name_antiplatelet,Drug_Name,ignore.case = T),.(Reference_Key,Prescription_Start_Date=ymd(Prescription_Start_Date),Prescription_End_Date=ymd(Prescription_End_Date))]
drug_antiplatelet <- drug_antiplatelet[!is.na(Prescription_End_Date)]
drug_antiplatelet_sameday <- drug_antiplatelet[Prescription_Start_Date==Prescription_End_Date][,.(Reference_Key,Date=Prescription_Start_Date,antiplatelet=1L)]
drug_antiplatelet <- drug_antiplatelet[!Prescription_Start_Date==Prescription_End_Date, .(
  Date = seq(Prescription_Start_Date, Prescription_End_Date, by = "day")
), by = .(Reference_Key, Prescription_Start_Date, Prescription_End_Date)][,.(Reference_Key,Date,antiplatelet=1)]
drug_antiplatelet <- rbind(drug_antiplatelet, drug_antiplatelet_sameday)
drug_antiplatelet <- unique(drug_antiplatelet)
drug_antiplatelet[,Date:=ymd(paste0(substring(Date,1,7),'-01'))]
drug_antiplatelet <- unique(drug_antiplatelet)
# sequential trail ---------------------------------------------------------

# 1.eligibility -----------------------------------------------------------
num_delete <- list()
run_seq <- function(index_date,cohort){
  message("-------------------------\n\n")
  message("index_date: ",index_date)
  each_trial <- copy(cohort)
  each_trial[,indx_date:=index_date]
  num_delete[["N. patients in the trial"]] <<- append(num_delete[["N. patients in the trial"]],nrow(each_trial))
  
  num_delete[["N. patients stop dialysis before index"]] <<- append(num_delete[["N. patients stop dialysis before index"]],each_trial[!(date.PD <= indx_date & date.PD.end >=indx_date),.N])
  each_trial_PD <- each_trial[date.PD <= indx_date & date.PD.end >=indx_date,] # 1.patients had a history of Peritoneal dialysis (PD) before the index date and still under dialysis
  message("Still with PD: ",nrow(each_trial_PD))
  
  num_delete[["N. patients donot have a af before index"]] <<- append(num_delete[["N. patients donot have a af before index"]],each_trial_PD[date.af > indx_date,.N])
  each_trial_PD_AF <- each_trial_PD[date.af <= indx_date] # 3.patients had a history of AF before the index date
  
  num_delete[["N. patients have AF before PD"]] <<- append(num_delete[["N. patients have AF before PD"]], each_trial_PD_AF[date.af < date.PD,.N])
  each_trial_PD_AF <- each_trial_PD_AF[date.af >= date.PD] # 3.AF after the index date
  message("PD patients with AF: ",nrow(each_trial_PD_AF))
  
  id_with_prevalant_oac <- each_trial_PD_AF[date.oac <= indx_date %m-% months(1),Reference_Key]
  num_delete[["N. patients have a history of OAC before index"]] <<- append(num_delete[["N. patients have a history of OAC before index"]],length(id_with_prevalant_oac))
  each_trial_PD_AF_rmpastOAC <- each_trial_PD_AF[!Reference_Key %in% id_with_prevalant_oac] # 4. remove prevalent user of OAC
  
  # each_trial_PD_AF_rmpastOAC <- each_trial_PD_AF_rmpastOAC[(date.oac >= date.af) | is.na(date.oac)]
  message("PD patients with no OAC one-month before index date: ",nrow(each_trial_PD_AF_rmpastOAC))
  rm(id_with_prevalant_oac)
  id_with_prevalant_outcome <- each_trial_PD_AF_rmpastOAC[date.death <= indx_date |
                                                            date.transplant <= indx_date |
                                                            date.hd <= indx_date,Reference_Key]
  num_delete[["N. patients have a history of outcome before index"]] <<- append(num_delete[["N. patients have a history of outcome before index"]],length(id_with_prevalant_outcome))
  each_trial_PD_AF_rmpastOAC_NoOutcome <-  each_trial_PD_AF_rmpastOAC[!Reference_Key %in% id_with_prevalant_outcome]
  rm(id_with_prevalant_outcome)
  message("No HD/transplant/disconti PD/death before index date: ",nrow(each_trial_PD_AF_rmpastOAC_NoOutcome))

  # id_with_past_hx <-   merge(each_trial_PD_AF_rmpastOAC_NoOutcome[,.(Reference_Key,indx_date)],
  #                            dx[grepl(covar[13,Regex],All_Diagnosis_Code_ICD9)],by="Reference_Key")[Reference_Date < indx_date,Reference_Key]
  # message("No past hx of hemo and stroke before index date: ",length(unique(id_with_past_hx)))
  # each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx <-  each_trial_PD_AF_rmpastOAC_NoOutcome[!Reference_Key %in% id_with_past_hx]
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx <- each_trial_PD_AF_rmpastOAC_NoOutcome
  message("\nNumber of patients:",each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[,uniqueN(Reference_Key)])
  
  # 2.identifying exposure --------------------------------------------------
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[,expo:=ifelse((date.oac > indx_date %m-% months(2) & date.oac <=indx_date) & !is.na(date.oac),1,0)]
  message("Expo: ",sum(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$expo),"\nNon-exp: ",each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[,.N-sum(expo)])
  # each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[,expo:=ifelse((date.oac > date.af & date.oac <=indx_date) & !is.na(date.oac),1,0)]
  # 3.covarates -------------------------------------------------------------
  each_trial_dx <- dx[Reference_Key %in% each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key & Reference_Date < index_date]
  px <- mapply(function(x) sapply(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key,function(y) getpx(y,x,each_trial_dx)),as.list(tibble::deframe(covar[,.(Name,Regex)])),USE.NAMES = T)
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- merge(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx,as.data.table(px,keep.rownames = T)[,setnames(.SD,"rn","Reference_Key")],by="Reference_Key")
  message("\nIncident: ",each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px[dx.stroke_embo==0,.N])
  # num_delete[["N. patients have a history of stroke before index"]] <<- append(num_delete[["N. patients have a history of stroke before index"]],each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px[dx.stroke_embo==1,.N])
  # each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px[dx.stroke_embo==0,]
  
  # 4. outcomes  ----------------------------------------------------------------
  each_trial_dx <- dx[Reference_Key %in% each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key & Reference_Date >=index_date ]
    # 4.1 all stroke ----------------------------------------------------------
  
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,"^43[0-6]|^444","stroke")[]
  
  # 4.2 all bleeding ----------------------------------------------------------------
  #"^43[0-6]|^444"
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_ischemicstroke,"ische")[]
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_HaemorrhagicStroke,"haem")[]
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_GI_Bleed,"GIbleed")[]
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_Other_Bleed,"Otherbleed")[]
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_MI,"mi")[]
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_fracture,"fracture")[]
  # each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px[date.stroke>=index_date | is.na(date.stroke)]
  return(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px)
}
getpx <- function(refkey,regex,data){
  temp <- ifelse(nrow(data[Reference_Key==refkey & grepl(regex,All_Diagnosis_Code_ICD9),])>0,1,0)
  return(temp)
}
search_outcomes <- function(cohort,each_trial_dx,icd,outcomename){
  
  temp <- merge(cohort,
                each_trial_dx[grepl(icd,All_Diagnosis_Code_ICD9),
                              .(Reference_Key,
                                Reference_Date,
                                out=1)][,.(Reference_Date=min(Reference_Date),out=1),Reference_Key],
                by="Reference_Key",all.x=T)
  temp[,out:=ifelse(is.na(out),0,out)]
  setnames(temp,c("Reference_Date","out"),c(paste0("date.",outcomename),paste0("out.",outcomename)))
  return(temp)
}

seqcohort <- pblapply(index_date[-1],run_seq,cohort)

# outcome definition ------------------------------------------------------
seqcohort <- lapply(seq_along(seqcohort), function(i) {
  df <- seqcohort[[i]]
  df$trial_id <- i  # Add the ID column with the index
  return(df)
})
saveRDS(seqcohort,"data/cleaned_20250508.RDS")




# 1：composite outcome (ischemic, hemorrhagic, mortality) 2. composite outcome (GI, other, Hemo)
# 3. ischemic; 4:hemorrhagic; 5: mortality; 6: mi 7:fracture
out_indicator <- 1
year <- 3
K <- year *12 
person_trial <- rbindlist(seqcohort)
person_trial[,id:=paste(trial_id,Reference_Key,sep = "-")]
demo <- as.data.table(readRDS("data/PD_OAC-demo.RDS"))
person_trial <- merge(person_trial,
                      unique(demo[,.(Reference_Key,Sex,dob=Date_of_Birth_yyyymmdd)]))


setkey(person_trial,NULL)
person_trial[,id:=paste(trial_id,Reference_Key,sep = "_")]
# primary outcome - stroke ------------------------------------------------
if(out_indicator %in% c(1)){
  message("Number of people with past hx of stroke ",person_trial[pt.allstroke==1,.N])
  person_trial <- person_trial[pt.allstroke==0]
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.ische,
                              date.haem,
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)]  
}else if(out_indicator==3){
  message("Number of people with past hx of ischemic ",person_trial[dx.stroke_embo==1,.N])
  person_trial <- person_trial[dx.stroke_embo==0]
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.ische,
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)] 
}else if(out_indicator==4){
  message("Number of people with past hx of hemorrhagic ",person_trial[pt.heam==1,.N])
  person_trial <- person_trial[pt.heam==0]
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.haem,
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)] 
}else if(out_indicator==5){
  message("Number of people with past hx of mortality 0")
  # person_trial <- person_trial[pt.heam==0]
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)] 
}else if(out_indicator == 6){
  message("Number of people with past hx of mi ",person_trial[pt.mi==1,.N])
  person_trial <- person_trial[pt.mi==0]
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.mi,
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)]
}else if(out_indicator==2){
  message("Number of people with past hx of bleeding",person_trial[pt.bleeding==1,.N])
  person_trial <- person_trial[pt.bleeding==0]
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.GIbleed,
                              date.Otherbleed,
                              date.haem,
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)]  
}else if(out_indicator==7){
  person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                              date.fracture,
                              date.death,
                              date.hd,
                              date.transplant,
                              indx_date %m+% years(year),
                              date.PD.end,na.rm = T)]
}

person_trial[,summary(as.numeric(obs_end-indx_date)/365.25)]
person_trial[,.N,expo]

# covaraites for ps 
person_trial[,age_indx:=as.numeric((indx_date-ymd(dob))/365.25)]
person_trial[,time_af_index:=as.numeric(indx_date-date.af)]
person_trial[,chadsvas_score:=dplyr::case_when(
  age_indx < 65 ~ 0,
  age_indx >= 65 & age_indx <75 ~ 1,
  age_indx > 75 ~ 2)+ifelse(Sex=="F",1,0)+dx.chf+dx.htn+dx.cbd+dx.stroke_embo*2+dx.pvd+dx.dm]

person_trial[,fu:=as.numeric(obs_end-indx_date)]
person_trial <- merge(person_trial,drug_antiplatelet,by.x=c("Reference_Key","indx_date"),by.y=c("Reference_Key","Date"),all.x=T)
setnames(person_trial,"antiplatelet","antiplatelet.baseline")
person_trial[,antiplatelet.baseline:=ifelse(is.na(antiplatelet.baseline),0,antiplatelet.baseline)]

if(out_indicator ==1){
  # table one ---------------------------------------------------------------
  vars <- c("fu","Sex", "age_indx","time_af_index", "chadsvas_score","dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
            "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
            "dx.liver_modsev", "dx.ulcers", "dx.asthma",
            "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
            "dx.dm_com1", "dx.htn", "dx.constipation","antiplatelet.baseline"
  )
  
  {
    
    for (i in 1:nrow(covar)) {
      var_name <- covar[i, Name]
      var_label <- covar[i, Description]
      
      if (var_name %in% names(person_trial)) {
        var_label(person_trial[[var_name]]) <- var_label
      }
    }
    
    var_label(person_trial[["age_indx"]]) <- "Age at index date (Years)"
    var_label(person_trial[["time_af_index"]]) <- "Time from AF to Index date (Days)"
    var_label(person_trial[["chadsvas_score"]]) <- "CHA₂DS₂-VASc Score"
    var_label(person_trial[["antiplatelet.baseline"]]) <- "Antiplatelet user history"
    var_label(person_trial[["fu"]]) <- "Follow-up period"
    labelled::var_label(person_trial) 
  }
  tabone <- CreateTableOne(vars = vars, 
                           factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                           strata = "expo", 
                           data = person_trial, test = F)
  
  notnormal <- c("age_indx","time_af_index","fu")
  tableone_b4_matching <- as.data.table(print(tabone, smd = T,varLabels=T,nonnormal = notnormal),keep.rownames = T)
  tableone_b4_matching[,rn:=gsub(" = 1","",rn)]
  tableone_b4_matching[,rn:=gsub(" = 0","",rn)]
  setnames(tableone_b4_matching,c("rn","before matching-non-user","before matching-user","before matching-SMD"))
  
  # matching step -----------------------------------------------------------
  set.seed(456)
  m.out <- matchit(expo~Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+
                     dx.ulcers+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation+time_af_index+
                     chadsvas_score+fu,
                   data = person_trial,replace = F,method = "nearest",ratio = 10,exact = ~trial_id
  )
  summary(m.out)
  plot(summary(m.out))
  plot(m.out, type = "jitter", interactive = FALSE)
  person_trial <- match.data(m.out)
  
  {
    
    for (i in 1:nrow(covar)) {
      var_name <- covar[i, Name]
      var_label <- covar[i, Description]
      
      if (var_name %in% names(person_trial)) {
        var_label(person_trial[[var_name]]) <- var_label
      }
    }
    
    var_label(person_trial[["age_indx"]]) <- "Age at index date (Years)"
    var_label(person_trial[["time_af_index"]]) <- "Time from AF to Index date (Days)"
    var_label(person_trial[["chadsvas_score"]]) <- "CHA₂DS₂-VASc Score"
    var_label(person_trial[["antiplatelet.baseline"]]) <- "Antiplatelet user history"
    var_label(person_trial[["fu"]]) <- "Follow-up period"
    labelled::var_label(person_trial) 
  }
  tabone <- CreateTableOne(vars = vars, 
                           factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                           strata = "expo", data = person_trial, test = F)
  tableone_after_matching <- as.data.table(print(tabone, smd = T,varLabels=T,nonnormal = notnormal),keep.rownames = T)
  tableone_after_matching[,rn:=gsub(" = 1","",rn)]
  tableone_after_matching[,rn:=gsub(" = 0","",rn)]
  setnames(tableone_after_matching,c("rn","after matching-non-user","after matching-user","after matching-SMD"))
  
  writexl::write_xlsx(merge(tableone_b4_matching,tableone_after_matching,by="rn"),"out/tableone.xlsx")
}


long_format <- person_trial[,.(id,indx_date,obs_end)
][,.(obs_date=seq.Date(indx_date,obs_end,by="1 months")),id]
long_format[,time:=0:(.N-1),id]
long_format[,timesqr:=time^2]
person_trial <- merge(person_trial,long_format,by="id")

person_trial[,age_indx := as.numeric((obs_date-ymd(dob)+1)/365.25)]
# https://www.mdcalc.com/calc/801/cha2ds2-vasc-score-atrial-fibrillation-stroke-risk#next-steps
person_trial[,chadsvas_score:=dplyr::case_when(
  age_indx < 65 ~ 0,
  age_indx >= 65 & age_indx <75 ~ 1,
  age_indx > 75 ~ 2)+ifelse(Sex=="F",1,0)+dx.chf+dx.htn+dx.cbd+dx.stroke_embo*2+dx.pvd+dx.dm]


person_trial[,time_af_index:=interval(date.af,obs_date) %/% months(1)]
person_trial[,fu:=as.numeric(obs_date-indx_date)]

if(out_indicator==3){
  # Ischemic stroke
  name <- paste0("ischemic","_",year,"_Yrs")
  person_trial[,out:=ifelse(date.ische >= obs_date & substring(date.ische,1,7)== substring(obs_date,1,7) ,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
}else if(out_indicator==4){
  # hemorrhagic stroke
  name <- paste0("hemorrhagic","_",year,"_Yrs")
  person_trial[,out:=ifelse(date.haem >= obs_date & substring(date.haem,1,7)== substring(obs_date,1,7) ,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
  
}else if(out_indicator==5){
  # all-cause mortality
  name <- paste0("mortality","_",year,"_Yrs")
  person_trial[,out:=ifelse(date.death >= obs_date & substring(date.death,1,7)== substring(obs_date,1,7) ,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
}else if(out_indicator==6){
  # Mi
  name <- paste0("mi","_",year,"_Yrs")
  person_trial[,out:=ifelse(date.mi >= obs_date & substring(date.mi,1,7)== substring(obs_date,1,7) ,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
}else if(out_indicator==2){
  name <- paste0("all_bleed","_",year,"_Yrs")
  person_trial[,out:=ifelse(pmin(date.haem,date.GIbleed,date.Otherbleed,na.rm = T) >= obs_date & substring(pmin(date.haem,date.GIbleed,date.Otherbleed,na.rm = T),1,7)== substring(obs_date,1,7) ,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
}else if(out_indicator==1){
  # Ischemic or hemorrhagic or death
  name <- paste0("ischhemodeath","_",year,"_Yrs")
  person_trial[,out:=ifelse(pmin(date.ische,date.haem,date.death,na.rm = T) >= obs_date & 
                              substring(pmin(date.ische,date.haem,date.death,na.rm = T),1,7)== substring(obs_date,1,7) ,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
}else if(out_indicator==7){
  # Ischemic or hemorrhagic or death
  name <- paste0("fracture","_",year,"_Yrs")
  person_trial[,out:=ifelse(date.fracture >= obs_date & substring(date.fracture,1,7)== substring(obs_date,1,7) ,1,0)]
}


# add antiplatete
person_trial <- merge(person_trial,drug_antiplatelet,by.x=c("Reference_Key","obs_date"),by.y=c("Reference_Key","Date"),all.x=T)
person_trial[,antiplatelet:=ifelse(is.na(antiplatelet),0,1)]

person_trial[,.N,expo]
person_trial[,.N,out]

# remove those records already have events
# person_trial <- person_trial[order(id,time),][,tempi:=cumsum(out),.(id)][tempi <= 1]


if(out_indicator==1){
  # plot for number of eligible people --------------------------------------
  test <- 
    person_trial[,uniqueN(Reference_Key),trial_id][order(trial_id)] |>
    ggplot(aes(trial_id,V1))+
    geom_bar(stat = "identity", fill = "skyblue")+
    labs(x = "Trial ID",
         y = "Number of eligibal patients")+
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
  
  save_plot("out/number of eligible.pdf", test, width_in = 10, height_in =3)
  test <- person_trial[,uniqueN(Reference_Key),.(trial_id,expo=factor(expo,levels=c(1,0),labels=c("Users","Non-users")))][order(trial_id,expo)] |>
    ggplot(aes(trial_id,V1, fill = expo))+
    geom_bar(stat = "identity",position="stack")+
    labs(x = "Trial ID",
         y = "Number of eligibal patients") +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
  save_plot("out/number of eligible by users or not.pdf", test, width_in = 10, height_in =3)
  
  # table one ----------------------------------------------------------------
  vars <- c("fu","Sex", "age_indx","time_af_index", "chadsvas_score","dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
            "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
            "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
            "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
            "dx.dm_com1", "dx.htn", "dx.constipation","antiplatelet"
  )
  {
    
    for (i in 1:nrow(covar)) {
      var_name <- covar[i, Name]
      var_label <- covar[i, Description]
      
      if (var_name %in% names(person_trial)) {
        var_label(person_trial[[var_name]]) <- var_label
      }
    }
    
    var_label(person_trial[["age_indx"]]) <- "Age at index date (Years)"
    var_label(person_trial[["time_af_index"]]) <- "Time from AF to Index date (Days)"
    var_label(person_trial[["chadsvas_score"]]) <- "CHA₂DS₂-VASc Score"
    var_label(person_trial[["antiplatelet.baseline"]]) <- "Antiplatelet user history"
    var_label(person_trial[["fu"]]) <- "Follow-up period"
    labelled::var_label(person_trial) 
  }
  tabone <- CreateTableOne(vars = vars, 
                           factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                           strata = "expo", data = person_trial, test = F)
  
  tabone <- as.data.table(print(tabone, smd = T,varLabels=T,nonnormal = notnormal),keep.rownames = T)
  labelled::var_label(person_trial) 
  notnormal <- c("age_indx","time_af_index","fu")
  # tabone[,rn:=gsub(" = 1","",rn)]
  # tabone[,rn:=gsub(" = 0","",rn)]
  print(tabone, smd = T,varLabels=T,nonnormal = notnormal)
  
  tabone <- CreateTableOne(vars = vars, 
                           factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                           strata = "expo", data = person_trial, test = F)
  
  unbalanced_variable <- as.data.table(as.data.frame(print(tabone,smd=T)),keep.rownames = T)[as.numeric(SMD) > 0.1,gsub(" = [1M] \\(%\\)| \\(mean \\(SD\\)\\)","",rn)]
  unbalanced_variable <- setdiff(unbalanced_variable,c("trial_id"))
  formula_rr_main <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(setdiff(gsub(".baseline","",unbalanced_variable),c("fu")),collapse = "+")))
  formula_rr_sex <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(setdiff(gsub(".baseline","",unbalanced_variable),c("Sex","fu")),collapse = "+")))
  formula_rr_age <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(setdiff(gsub(".baseline","",unbalanced_variable),c("age_indx","fu")),collapse = "+")))
  formula_rr_main_cum_plot <- formula(paste("out~expo+expo*Sex+I(time)+I(timesqr)+I(time*expo)+I(timesqr*expo)+",paste0(setdiff(gsub(".baseline","",unbalanced_variable),c("fu")),collapse = "+")))
  formula_rr_sex_cum_plot <- formula(paste("out~expo+I(time)+I(timesqr)+I(time*expo)+I(timesqr*expo)+",paste0(setdiff(gsub(".baseline","",unbalanced_variable),c("Sex","fu")),collapse = "+")))
  formula_rr_age_cum_plot <- formula(paste("out~expo+I(time)+I(timesqr)+I(time*expo)+I(timesqr*expo)+",paste0(setdiff(gsub(".baseline","",unbalanced_variable),c("age_indx","fu")),collapse = "+")))
}



# models for risk ratio (pooled logistic as HR)  ------------------------------------------------------------------


rm(main_risk_results,male_risk_results,female_risk_results,younger_risk_results,elder_risk_results)
rm(main_risk_results_o,male_risk_results_o,female_risk_results_o,younger_risk_results_o,elder_risk_results_o)

R <- 500
# main analysis -----------------------------------------------------------
set.seed(456)
current_bootstrap <- 0
elig_ids <- data.table(id = unique(person_trial$id))
formula_rr <- formula_rr_main
formula_rd <- formula_rr_main_cum_plot
pb <- txtProgressBar(min=0, max=R, style=3)
main_risk_results <- boot(
  data = elig_ids,
  statistic = std.boot,
  R = R
)
main_risk_results_o <- orgnize_ci(main_risk_results,"main",elig_ids)


# male analysis -----------------------------------------------------------
set.seed(456)
current_bootstrap <- 0
elig_ids <- data.table(id = unique(person_trial[Sex=="M",id]))
formula_rr <- formula_rr_sex
formula_rd <- formula_rr_sex_cum_plot
pb <- txtProgressBar(min=0, max=R, style=3)
male_risk_results <- boot(
  data = elig_ids,
  statistic = std.boot,
  R = R
)
male_risk_results_o <- orgnize_ci(male_risk_results,"male",elig_ids)


# female analysis -----------------------------------------------------------
set.seed(456)
current_bootstrap <- 0
elig_ids <- data.table(id = unique(person_trial[Sex=="F",id]))
formula_rr <- formula_rr_sex
formula_rd <- formula_rr_sex_cum_plot
pb <- txtProgressBar(min=0, max=R, style=3)
female_risk_results <- boot(
  data = elig_ids,
  statistic = std.boot,
  R = R
)
female_risk_results_o <- orgnize_ci(female_risk_results,"female",elig_ids)


# younger analysis -----------------------------------------------------------
set.seed(456)
current_bootstrap <- 0
elig_ids <- data.table(id = unique(person_trial[age_indx<65,id]))
formula_rr <- formula_rr_age
formula_rd <- formula_rr_age_cum_plot
pb <- txtProgressBar(min=0, max=R, style=3)
younger_risk_results <- boot(
  data = elig_ids,
  statistic = std.boot,
  R = R
)
younger_risk_results_o <- orgnize_ci(younger_risk_results,"younger",elig_ids)

# elder analysis -----------------------------------------------------------
set.seed(456)
current_bootstrap <- 0
elig_ids <- data.table(id = unique(person_trial[age_indx>=65,id]))
formula_rr <- formula_rr_age
formula_rd <- formula_rr_age_cum_plot
pb <- txtProgressBar(min=0, max=R, style=3)
elder_risk_results <- boot(
  data = elig_ids,
  statistic = std.boot,
  R = R
)
elder_risk_results_o <- orgnize_ci(elder_risk_results,"elder",elig_ids)

out_est <- rbind(main_risk_results_o,
                 male_risk_results_o,
                 female_risk_results_o,
                 younger_risk_results_o,
                 elder_risk_results_o) %>% mutate(out=name,.before="model_type")

if(out_indicator==1){
  pool_results <- list()
  pool_results[[name]]$main <- main_risk_results
  pool_results[[name]]$main_o <- main_risk_results_o
  pool_results[[name]]$male <- male_risk_results
  pool_results[[name]]$male_o <- male_risk_results_o
  pool_results[[name]]$female <- female_risk_results
  pool_results[[name]]$female_o <- female_risk_results_o
}else{
  pool_results[[name]]$main <- main_risk_results
  pool_results[[name]]$main_o <- main_risk_results_o
  pool_results[[name]]$male <- male_risk_results
  pool_results[[name]]$male_o <- male_risk_results_o
  pool_results[[name]]$female <- female_risk_results
  pool_results[[name]]$female_o <- female_risk_results_o
}


# openxlsx::write.xlsx(out_est,
#                     "stroke_and_oac_results.xlsx",
#                     sheetName = paste0(name,"-",format(Sys.time(), '%d-%m-%Y')),
#                     append = T)

file_path <- "out/stroke_and_oac_results_20250508.xlsx"
library('openxlsx')
# Load existing workbook or create a new one if it doesn't exist
if (file.exists(file_path)) {
  wb <- loadWorkbook(file_path)
} else {
  wb <- createWorkbook()
}
new_sheet_name <- paste0(name,"-",format(Sys.time(), '%d-%m-%Y'))
addWorksheet(wb, sheetName = new_sheet_name)
writeData(
  wb, 
  sheet = new_sheet_name, 
  x = out_est,
  startCol = 1,  # Optional: specify starting column/row
  startRow = 1
)
saveWorkbook(wb, file_path, overwrite = TRUE)



# https://rpubs.com/mbounthavong/sample_size_power_analysis_R
library("pwr")
### alpha = sig.level option and is equal to 0.05
### power = 0.80
### p1 = 0.60 
### p2 = 0.50
plot(pwr.2p.test(h = 0.5, sig.level = 0.05, power = .80,alternative = "greater"))

plot(pwr.2p.test(h = 0.5, sig.level = 0.05, power = .80,alternative = "greater"))


power1 <-pwr.2p.test(h = ES.h(p1 = 2/(2+82), p2 = 280/(8487+280)), sig.level = 0.05, power = .80)
power1
plot(power1)

person_trial[fu==0][,.N,keyby=.(expo,out)]
# 1:     0     0  8487
# 2:     0     1   280
# 3:     1     0    82
# 4:     1     1     2

pwr.2p2n.test(h = ES.h(p1 = 3/(3+81), p2 = 32/(32+8735)),n1 = (3+81),n2=(32+8735),
              sig.level = 0.05) #0.07367

person_trial[,.N,keyby=.(expo,out)]
# 1:     0     0 138740
# 2:     0     1   4196
# 3:     1     0   1263
# 4:     1     1     34
pwr.2p2n.test(h = ES.h(p1 = 34/(34+1263), p2 = 4196/(4196+138740)),n1 = (34+1263),n2=(4196+138740),
              sig.level = 0.05) #0.1053


pwr.2p2n.test(h = ES.h(p1 = 2/(2+82), p2 = 280/(8487+280)),n1 = (2+82),n2=(8487+280),
              sig.level = 0.05, power=0.8)

plot(pwr.2p2n.test(h = ES.h(p1 = 2/(2+82), p2 = 280/(8487+280)),n1 = (2+82),n2=(8487+280),
              sig.level = 0.05))


plot(pwr.2p2n.test(h = ES.h(p1 = 2/(2+82), p2 = 280/(8487+280)),n1 = (2+82),n2=(8487+280),
                   sig.level = 0.05, power = 0.8))

p1 <- seq(0.5, 1.0, 0.05)
power1 <-pwr.2p.test(h = ES.h(p1 = p1, p2 = 0.50),
                     n = 388,
                     sig.level = 0.05)
powerchange <- data.frame(p1, power = power1$power * 100)
plot(powerchange$p1, 
     powerchange$power, 
     type = "b", 
     xlab = "Proportion of Responders in Treatment A", 
     ylab = "Power (%)")

ischemic_stroke_5y <- rbind(
  orgnize_ci(main_risk_results,"main"),
  orgnize_ci(male_risk_results,"male"),
  orgnize_ci(female_risk_results,"female")) %>% mutate(out="ischemic",.before="model")

person_trial[,out:=ifelse(date.haem >= obs_date & substring(date.haem,1,7)== substring(obs_date,1,7) ,1,0)]
person_trial[,out:=ifelse(is.na(out),0,out)]

haemorraghic_stroke_5y <- rbind(
  orgnize_ci(main_risk_results,"main"),
  orgnize_ci(male_risk_results,"male"),
  orgnize_ci(female_risk_results,"female")) %>% mutate(out="haemorraghic",.before="model")

writexl::write_xlsx(rbind(ischemic_stroke_5y,haemorraghic_stroke_5y),
                    "ischemic_haemorraghic_incident_stroke_results_5y.xlsx")



rbind(
  orgnize_ci(younger_risk_results,"younger"),
  orgnize_ci(elder_risk_results,"elder")) %>% mutate(out=name,.before="model") %>% 
  mutate(estimate=paste0( sprintf("%.2f",round(est,2))," (",
                          sprintf("%.2f",round(lower,2)),", ",
                          sprintf("%.2f",round(upper,2)),")")) %>% 
  dplyr::select(out,model,arm,estimate) %>% 
  dcast(out+model~arm) %>% View()



# old without ps matching
# # main
# 0.39565  0.34130 -0.05435  0.86263
# ( 0.2854,  0.5133 )
# ( 0.1708,  0.6105 )  
# (-0.1909,  0.2368 ) 
# ( 0.5021,  1.6834 )  
# 
# #male
# 0.39634  0.34195 -0.05439  0.86278
# ( 0.2779,  0.5171 )  
# ( 0.1881,  0.5475 )
# (-0.2059,  0.1537 ) 
# ( 0.4836,  1.4133 )  
# 
# #female
# 0.3692  0.3009 -0.0683  0.8150
# ( 0.2538,  0.5244 ) 
# ( 0.1333,  0.7407 )  
# (-0.2024,  0.3226 )  
# ( 0.4382,  1.6660 )  
# 
# #incidence
# 0.35754 0.39515 0.03762 1.10522
# ( 0.2065,  0.4366 )  
# ( 0.1282,  0.9717 ) 
# (-0.1894,  0.6552 )  
# ( 0.423,  3.766 ) 
# 
# #recurrent
# 0.6596  0.2725 -0.3872  0.4130
# ( 0.2065,  0.4366 ) 
# ( 0.1282,  0.9717 )  
# (-0.1894,  0.6552 )  
# ( 0.423,  3.766 )  




library(ggplot2)
library(ggview)
# rbind(rbindlist(lapply(names(all_sub_cox),function(x) all_sub_cox[[x]]$estimation[1,] %>% mutate(model=x))) %>% mutate(regression="cox") %>% select(-robust.se),
#       rbindlist(lapply(names(all_sub_logistic),function(x) all_sub_logistic[[x]]$estimation %>% filter(term=="expo") %>% mutate(model=x))) %>% mutate(regression="logistic") ) %>%
#   filter(model %in% c("main","M","F",1,0)) %>%
#   mutate(model=factor(model,levels=c("main","M","F","0","1","M.0","M.1","F.0","F.1"),
#                       labels=c("Main analysis","Male","Female","Incident cases",
#                                "Prevelant cases",
#                                "Male incident cases","Male prevelant cases",
#                                "Female incident cases","Female prevalent cases"))) %>% 
#   ggplot(.,aes(x=estimate,y=model)) +
#   geom_point(size=0.6)+
#   geom_errorbar(aes(xmin=conf.low,xmax=conf.high),width=0.3)+
#   scale_y_discrete(limits=rev)+
#   facet_grid(.~regression)+
#   geom_vline(xintercept=1, linetype="dashed", 
#              color = "red")
# 
# ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
#   geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1, position=pd) +
#   geom_line(position=pd) +
#   geom_point(position=pd)


# abstract with pooled logistic  ------------------------------------------

abstract <- function(outcome,sub="All"){
  person_trial <- rbindlist(seqcohort)
  demo <- as.data.table(readRDS("data/PD_OAC-demo.RDS"))
  person_trial <- merge(person_trial,
                        unique(demo[,.(Reference_Key,Sex,dob=Date_of_Birth_yyyymmdd)]))
  setkey(person_trial,NULL)
  person_trial[,id:=paste(trial_id,Reference_Key,sep = "_")]
  if(outcome=="date.stroke"){
    person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                                date.stroke,
                                date.death,
                                date.hd,
                                date.transplant,
                                date.PD.end,na.rm = T)]
  }else if(outcome=="date.ische"){
    person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                                date.ische,
                                date.death,
                                date.hd,
                                date.transplant,
                                date.haem,
                                date.PD.end,na.rm = T)]
  }else{
    person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                                date.haem,
                                date.death,
                                date.hd,
                                date.transplant,
                                date.ische,
                                date.PD.end,na.rm = T)]
  }
  long_format <- person_trial[,.(id,indx_date,obs_end)
  ][,.(obs_date=seq.Date(indx_date,obs_end,by="1 months")),id]
  long_format[,time:=1:.N,id]
  person_trial <- merge(person_trial,long_format,by="id")
  
  person_trial[,age_indx := as.numeric((obs_date-ymd(dob)+1)/365.25)]
  # https://www.mdcalc.com/calc/801/cha2ds2-vasc-score-atrial-fibrillation-stroke-risk#next-steps
  person_trial[,chadsvas_score:=dplyr::case_when(
    age_indx < 65 ~ 0,
    age_indx >= 65 & age_indx <75 ~ 1,
    age_indx > 75 ~ 2)+ifelse(Sex=="F",1,0)+dx.chf+dx.htn+dx.cbd+dx.stroke_embo*2+dx.pvd+dx.dm]
  
  
  person_trial[,time_af_index:=as.numeric(obs_date-date.af)+1]
  person_trial[,fu:=as.numeric(obs_end-obs_date+1)]
  person_trial[,out:=ifelse(substr(get(outcome),1,7) == substr(obs_date,1,7),1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
  
  pool.logistic.model <- glm(out ~ expo + Sex + age_indx + dx.mi + dx.stroke_embo + I(time) + I(time^2) , data = person_trial,family =binomial)
  pool.logistic.model.M <- glm(out ~ expo + age_indx + dx.mi + dx.stroke_embo + I(time) + I(time^2), data = person_trial[Sex=="M"],family =binomial)
  pool.logistic.model.F <- glm(out ~ expo + age_indx + dx.mi + dx.stroke_embo + I(time) + I(time^2), data = person_trial[Sex=="F"],family =binomial)
  pool.logistic.model.I <- glm(out ~ expo + Sex + age_indx + dx.mi + I(time) + I(time^2) , data = person_trial[dx.stroke_embo==0],family =binomial)
  pool.logistic.model.P <- glm(out ~ expo + Sex + age_indx + dx.mi + I(time) + I(time^2), data = person_trial[dx.stroke_embo==1],family =binomial)
  pool.logistic.model.IM <- glm(out ~ expo + age_indx + dx.mi + I(time) + I(time^2), data = person_trial[dx.stroke_embo==1 & Sex=="M"],family =binomial)
  pool.logistic.model.IF <- glm(out ~ expo + age_indx + dx.mi + I(time) + I(time^2), data = person_trial[dx.stroke_embo==1 & Sex=="F"],family =binomial)
  # return(person_trial[id=="136_3569443",.(id,indx_date,obs_end,obs_date,month,out,expo,date.stroke)])
  return(list(pri=pool.logistic.model,
              male=pool.logistic.model.M,
              femal=pool.logistic.model.F,
              inci=pool.logistic.model.I,
              prev=pool.logistic.model.P,
              inci_m=pool.logistic.model.IM,
              inci_f=pool.logistic.model.IF))
}
results_abs <- list(abstract("date.stroke"),
                    abstract("date.ische"),
                    abstract("date.haem"))

results_abs_tidy <- data.table()
for(i in results_abs){
  for(j in i)
    results_abs_tidy <- rbind(results_abs_tidy,
                              tidy(j,exp=T,conf.int = T) %>% filter(term=="expo"))
}
results_abs_tidy$out <- rep(c("stroke","ische","haem"),each=7)
results_abs_tidy$model <- rep(c("pri","M","F","I","R","IM","IF"),times=3)
results_abs_tidy
results_abs_tidy[model=="pri"]
results_abs_tidy[model %in% c("I","IM","IF")]
results_abs_tidy[out=="stroke"]





person_trial <- rbindlist(seqcohort)
demo <- as.data.table(readRDS("data/PD_OAC-demo.RDS"))
person_trial <- merge(person_trial,
                      unique(demo[,.(Reference_Key,Sex,dob=Date_of_Birth_yyyymmdd)]))
setkey(person_trial,NULL)
person_trial[,id:=paste(trial_id,Reference_Key,sep = "_")]
person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                            date.ische,
                            date.death,
                            date.hd,
                            date.transplant,
                            date.haem,
                            date.PD.end,na.rm = T)]

long_format <- person_trial[,.(id,indx_date,obs_end)
][,.(obs_date=seq.Date(indx_date,obs_end,by="1 months")),id]
long_format[,month:=1:.N,id]
person_trial <- merge(person_trial,long_format,by="id")

person_trial[,age_indx := as.numeric((obs_date-ymd(dob)+1)/365.25)]
# https://www.mdcalc.com/calc/801/cha2ds2-vasc-score-atrial-fibrillation-stroke-risk#next-steps
person_trial[,chadsvas_score:=dplyr::case_when(
  age_indx < 65 ~ 0,
  age_indx >= 65 & age_indx <75 ~ 1,
  age_indx > 75 ~ 2)+ifelse(Sex=="F",1,0)+dx.chf+dx.htn+dx.cbd+dx.stroke_embo*2+dx.pvd+dx.dm]


person_trial[,time_af_index:=as.numeric(obs_date-date.af)+1]
person_trial[,fu:=as.numeric(obs_end-obs_date+1)]
# person_trial[,out:=ifelse(get("date.ische") >= obs_date & get("date.ische") - obs_date <= fu,1,0)]
person_trial[,out:=ifelse(substr(get("date.ische"),1,7) == substr(obs_date,1,7),1,0)]
person_trial[,out:=ifelse(is.na(out),0,out)]

pool.logistic.model.M.I <- glm(out ~ expo + age_indx + dx.mi + I(month) + I(month^2) , data = person_trial[Sex=="M" &dx.stroke_embo==0],family =binomial)
pool.logistic.model.F.I <- glm(out ~ expo + age_indx + dx.mi + I(month) + I(month^2) , data = person_trial[Sex=="F" & dx.stroke_embo==0],family =binomial)
tidy(pool.logistic.model.M.I,exponentiate = T,conf.int = T)
tidy(pool.logistic.model.F.I,exponentiate = T,conf.int = T)
# reccurent case-> confounder by indication

# for abstract only incident cases ----------------------------------------

abstract <- function(outcome,sub="All"){
  person_trial <- rbindlist(seqcohort)
  demo <- as.data.table(readRDS("data/PD_OAC-demo.RDS"))
  person_trial <- merge(person_trial,
                        unique(demo[,.(Reference_Key,Sex,dob=Date_of_Birth_yyyymmdd)]))
  person_trial[,age_indx := as.numeric((indx_date-ymd(dob)+1)/365.25)]
  
  person_trial[,.N,expo]
  if(outcome=="date.stroke"){
    person_trial[,fu:=pmin(ymd("2022-12-31"),
                           get(outcome),
                           date.death,
                           date.hd,
                           date.transplant,
                           date.PD.end,na.rm = T)-indx_date+1]
  }else if(outcome=="date.ische"){
    person_trial[,fu:=pmin(ymd("2022-12-31"),
                           get(outcome),
                           date.death,
                           date.hd,
                           date.transplant,
                           date.haem,
                           date.PD.end,na.rm = T)-indx_date+1]
  }else{
    person_trial[,fu:=pmin(ymd("2022-12-31"),
                           get(outcome),
                           date.death,
                           date.hd,
                           date.transplant,
                           date.ische,
                           date.PD.end,na.rm = T)-indx_date+1]
  }
  
  # person_trial[,fu:=pmin(ymd("2022-12-31"),
  #                        get(outcome),
  #                        date.death,
  #                        date.hd,
  #                        date.transplant,
  #                        date_adjusted,
  #                        date.haem,
  #                        date.ische,
  #                        date.PD.end,na.rm = T)-indx_date+1]
  person_trial[,time_af_index:=as.numeric(indx_date-date.af)+1]
  person_trial[,out:=ifelse(get(outcome)>=indx_date & get(outcome) - indx_date <= fu,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
  person_trial[,.N,out]
  
  person_trial[,fu:=as.numeric(fu)]
  
  if(sub=="All"){
    cox.model <- coxph(Surv(fu, out) ~ expo + Sex+age_indx+dx.mi+dx.stroke_embo+time_af_index, 
                       data = person_trial[dx.stroke_embo==0],robust = T)
  }else if(sub=="M"){
    cox.model <- coxph(Surv(fu, out) ~ expo +age_indx+dx.mi+dx.stroke_embo+time_af_index, 
                       data = person_trial[dx.stroke_embo==0 & Sex=="M"],robust = T)
  }else if(sub=="F"){
    cox.model <- coxph(Surv(fu, out) ~ expo +age_indx+dx.mi+dx.stroke_embo+time_af_index, 
                       data = person_trial[dx.stroke_embo==0 & Sex=="F"],robust = T)
  }
  tidy(cox.model,exp = TRUE,conf.int=T) %>% 
    mutate_if(is.numeric, round, 5)
}
abstract("date.stroke")
abstract("date.ische")
abstract("date.haem")

abstract("date.stroke",sub = "M")
abstract("date.stroke",sub = "F")


outcome <- "date.ische"
{
  person_trial <- rbindlist(seqcohort)
  demo <- as.data.table(readRDS("data/PD_OAC-demo.RDS"))
  person_trial <- merge(person_trial,
                        unique(demo[,.(Reference_Key,Sex,dob=Date_of_Birth_yyyymmdd)]))
  person_trial[,age_indx := as.numeric((indx_date-ymd(dob)+1)/365.25)]
  
  person_trial[,.N,expo]
  person_trial[,fu:=pmin(ymd("2022-12-31"),
                         get(outcome),
                         date.death,
                         date.hd,
                         date.transplant,
                         date.haem,
                         date.PD.end,na.rm = T)-indx_date+1]
  person_trial[,time_af_index:=as.numeric(indx_date-date.af)+1]
  person_trial[,out:=ifelse(get(outcome)>=indx_date & get(outcome) - indx_date <= fu,1,0)]
  person_trial[,out:=ifelse(is.na(out),0,out)]
  person_trial[,.N,out]
  
  person_trial[,fu:=as.numeric(fu)]
  
  cox.model <- coxph(Surv(fu, out) ~ expo + Sex+age_indx+dx.mi+dx.stroke_embo+time_af_index, 
                     data = person_trial[dx.stroke_embo==0],robust = T)
  #cox.model <- coxph(Surv(fu, out) ~ expo + Sex+age_indx+dx.mi+dx.stroke_embo, data = person_trial,robust = T)
  #cox.zph(cox.model)
  tidy(cox.model,exp = TRUE,conf.int=T) %>% 
    mutate_if(is.numeric, round, 5)
}





# Propensity score ---------------------------------------------------------
ps.model <- glm(expo~Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+dx.ulcers+dx.stroke_embo+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation+time_af_index+chadsvas_score,family = binomial, data = person_trial)


# Unstabilized weight
person_trial$ps_score <- predict(ps.model, type = "response")
person_trial$usweight <- with(person_trial, ifelse(expo==1,
                                                   1/ps_score, 1/(1-ps_score)))

round(summary(person_trial$usweight), 2)

# Stabilized weight
person_trial$sweight <- with(person_trial, ifelse(expo==1,
                                                  mean(expo==1)/ps_score,
                                                  (1-mean(expo==1))/(1-ps_score)))
round(summary(person_trial$sweight), 2)


# Truncating unstabilized weight
person_trial <- person_trial %>%
  mutate(usweight_t = pmin(pmax(usweight, quantile(usweight, 0.01)),
                           quantile(usweight, 0.99)))
summary(person_trial$usweight_t)

# Truncating stabilized weight
person_trial <- person_trial %>%
  mutate(sweight_t = pmin(pmax(sweight, quantile(sweight, 0.01)),
                          quantile(sweight, 0.99)))
summary(person_trial$sweight_t)


# Covariates
vars <- c("expo", "dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
          "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
          "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
          "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
          "dx.dm_com1", "dx.htn", "dx.constipation", "trial_id",
          "fu", "Sex", "dob", "age_indx","chadsvas_score")

# Design with truncated unstabilized weight
design.unstab <- svydesign(ids = ~Reference_Key, weights = ~usweight_t, data = person_trial)
design.stab <- svydesign(ids = ~Reference_Key, weights = ~sweight_t, data = person_trial)

tab.stab <- svyCreateTableOne(vars = vars, strata = "expo", data = design.stab, test = F)

as.data.table(print(tab.stab, smd = T),keep.rownames = T)[as.numeric(SMD)<0.1,rn]

model<-svycoxph(Surv(fu, out)~expo+dx.chf+dx.pvd+dx.paralysis+dx.liver_mild+dx.ulcers+dx.asthma+dx.dm_com1+dx.constipation,design=design.stab)
tidy(model,exp=T,conf.int=T)




tidy(coxph(Surv(fu, out) ~ expo + Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+dx.ulcers+dx.stroke_embo+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation, data = person_trial,subset=dx.cbd==1 |dx.stroke_embo ==1 ,robust = T),exp=TRUE,conf.int=T)

tidy(coxph(Surv(fu, out) ~ expo + Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+dx.ulcers+dx.stroke_embo+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation, data = person_trial,subset= dx.cbd==0 & dx.stroke_embo==0 ,robust = T),exp=TRUE,conf.int=T)


# library(survey)
# svycox.design <- svydesign(ids=~1, weights = 1, data=person_trial)
# svycox.model <- svycoxph(formula = Surv(fu, out) ~ expo+Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+dx.ulcers+dx.stroke_embo+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation, design = svycox.design)



# plot --------------------------------------------------------------------

# Create dataset with all time points for each individual under each treatment level
K=140
treat0 <- splitstackshape::expandRows(person_trial[which(person_trial$month==1),], count=K, count.is.col=F) 
treat0$month <- rep(seq(0, K-1), nrow(person_trial[which(person_trial$month==1),]))
# Under flu 
treat0$expo <- as.numeric(as.character(treat0$expo))
treat0$expo <- 0
# Under covid 
treat1 <- treat0
treat1$expo <- 1


# Extract predicted values from pooled logistic regression model for each person-time row
# Predicted values correspond to discrete-time hazards
treat0$p.event0 <- predict(pool.logistic.model, treat0, type="response")
treat1$p.event1 <- predict(pool.logistic.model, treat1, type="response")

# Obtain predicted survival probabilities from discrete-time hazards
treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
treat0.surv$risk0 <- 1 - treat0.surv$surv0
treat1.surv$risk1 <- 1 - treat1.surv$surv1

# Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
risk0 <- aggregate(treat0.surv[c("expo", "month", "risk0")], by=list(treat0.surv$month), FUN=mean)[c("expo", "month", "risk0")] 
treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
risk1 <- aggregate(treat1.surv[c("expo", "month", "risk1")], by=list(treat1.surv$month), FUN=mean)[c("expo", "month", "risk1")] 

# Prepare data
graph.pred <- merge(risk0, risk1, by=c("month"))
# Edit data frame to reflect that risks are estimated at the END of each interval
graph.pred$time_0 <- graph.pred$month + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred)

### Use pooled logistic regression estimates to compute causal estimates ###

# 30-days risk in flu
risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
risk0.plr

# 30-days risk in COVID
risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
risk1.plr

# 30-days risk difference
rd.plr <- risk1.plr - risk0.plr
rd.plr

# 30-days risk ratio
rr.plr <- risk1.plr / risk0.plr
rr.plr

ggplot(graph, 
       aes(x=time_0, y=risk)) + # set x and y axes
  geom_line(aes(y = risk1, # create line for COVID group
                color = "OAC"),
            size = 1.5) + 
  geom_line(aes(y = risk0, # create line for no FLU group
                color = "NonUsers"),
            size = 1.5) +
  xlab("days") + # label x axis
  # scale_x_continuous(limits = c(0, 30), # format x axis
  # breaks=c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24,26,28,30)) + 
  ylab("Cumulative Incidence (%)") + # label y axis
  # scale_y_continuous(limits=c(0, 0.1), # format y axis
  # breaks=c(0, 0.02, 0.04, 0.06, 0.08, 0.1),
  # labels=c("0%","2%",  "4%", "6%", "8%","10%")) + 
  theme_minimal()+ # set plot theme elements
  # theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
  #       axis.line = element_line(colour = "black"),
  #       legend.title = element_blank(),
  #       panel.grid.major.x = element_blank(),
  #       panel.grid.minor.x = element_blank(),
  #       panel.grid.minor.y = element_blank(),
  #       panel.grid.major.y = element_blank())+
  # font("xlab",size=14)+
  # font("ylab",size=14)+
  # font("legend.text",size=10)+
  # scale_color_manual(values=c("#E7B800","#2E9FDF"), # set colors
#                    breaks=c('Flu', 'COVID')) +
theme(plot.title = element_text(face = "bold",
                                size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      # panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
)
