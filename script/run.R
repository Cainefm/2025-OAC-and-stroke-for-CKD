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
options(scipen = 6, digits = 4)
memory.limit(30000000)

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

# Sequential trial calendar: daily index + direction-2 exposure --------------
# Direction 2: expo = 1 if any OAC prescription covers indx_date (current use
# vs not). No prevalent-user exclusion (everyone on vs off OAC at index).
# Impute missing Prescription_End_Date as start + 90 days (edit if needed).
STUDY_START <- ymd("2007-01-01")
STUDY_END <- ymd("2022-12-31")
OAC_END_IMPUTE_DAYS <- 90L
VERBOSE_TRIAL_MSG <- FALSE # set TRUE to print every trial day (very slow I/O)

index_date <- seq.Date(STUDY_START, STUDY_END, by = "day")
# index_date <- seq.Date(STUDY_START, STUDY_END, by = "week") # lighter load
# index_date <- index_date[1:120] # debug: subset of days

# Load full cohort ------------------------------------------------------
cohort <- setDT(readRDS("data/pd_oac-cohort-20240813.RDS")) # 1632 patients/ records

ymd_cols <- grep("date",colnames(cohort),value=T) # change to date format
cohort[,(ymd_cols):=lapply(.SD, ymd),.SDcols=ymd_cols]

cohort[,date.PD.end:=pmin(date.transplant,date.PD.end,date.hd,na.rm = T)]

# Load covariate and dx data ----------------------------------------------
covar <- setDT(readxl::read_xlsx("documents/PD OAC-Protocol.xlsx",sheet = "Dx.Cov"))
dx <- setDT(readRDS("data/PD_OAC-dx.RDS"))
dx[,Reference_Date:=ymd(Reference_Date)]
dx <- dx[grepl(covar[,paste(Regex,collapse = "|")],All_Diagnosis_Code_ICD9)]
dx <- dx[,.(Reference_Key,Reference_Date,All_Diagnosis_Code_ICD9)][,.SD[Reference_Date==min(Reference_Date)],.(Reference_Key,All_Diagnosis_Code_ICD9)]

rx <- readRDS("data/PD_OAC-drug.RDS")
setDT(rx)

# OAC prescriptions (intervals) for direction-2 exposure at index ------------
drug_name_oac <- paste0(
  "WARFARIN|COUMADIN|ACENOCOUMAROL|PHENINDIONE|PHENPROCOUMON|",
  "RIVAROXABAN|XARELTO|",
  "APIXABAN|ELIQUIS|",
  "DABIGATRAN|PRADAXA|",
  "EDOXABAN|LIXIANA|SAVAYSA"
)
drug_oac <- rx[grepl(drug_name_oac, Drug_Name, ignore.case = TRUE),
               .(Reference_Key,
                 Prescription_Start_Date = ymd(Prescription_Start_Date),
                 Prescription_End_Date = ymd(Prescription_End_Date))]
drug_oac <- drug_oac[!is.na(Prescription_Start_Date)]
drug_oac[is.na(Prescription_End_Date),
         Prescription_End_Date := Prescription_Start_Date + OAC_END_IMPUTE_DAYS]
drug_oac <- drug_oac[Prescription_Start_Date <= Prescription_End_Date]

drug_name_antiplatelet <- "ASPIRIN|CARDIPRIN|CARTIA|PROPIRIN|ASPI-COR|ASPILETS|AGGRENOX|CLOPIE?DOG?(ER|RE)L|PLAVIX|PRASUGREL|EFFIENT|TICAGREC?LOR|BRILINTA|TICLOPIDINE|TICLID|CANGRELOR|KENGREX?AL|D?IPYRIDAMOLE|PERSANTIN|PROCARDIN|AGGRENOX|CILOSTAZOL|PLETAAL|VORAPAXAR|ZONTIVITY|ABCIXIMAB|REOPRO|EPTIFIBATIDE|INTEGRILIN"
drug_antiplatelet <- rx[grepl(drug_name_antiplatelet,Drug_Name,ignore.case = T),.(Reference_Key,Prescription_Start_Date=ymd(Prescription_Start_Date),Prescription_End_Date=ymd(Prescription_End_Date))]
drug_antiplatelet <- drug_antiplatelet[!is.na(Prescription_End_Date)]
drug_antiplatelet <- drug_antiplatelet[, .(
  Date = seq(Prescription_Start_Date, Prescription_End_Date, by = "day")
), by = .(Reference_Key, Prescription_Start_Date, Prescription_End_Date)][,.(Reference_Key,Date,antiplatelet=1)]
drug_antiplatelet <- unique(drug_antiplatelet)
# Daily dates (needed for daily indx_date merges at baseline and follow-up)
# sequential trail ---------------------------------------------------------

# Direction-2: OAC active on indx_date (any interval overlaps index day)
assign_oac_expo <- function(dt, drug_oac, ind) {
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

# 1.eligibility -----------------------------------------------------------
run_seq <- function(index_date, cohort) {
  if (VERBOSE_TRIAL_MSG) {
    message("-------------------------\n\nindex_date: ", index_date)
  } else if (format(index_date, "%d") == "01") {
    message("index_date (1st of month): ", index_date)
  }
  each_trial <- copy(cohort)
  each_trial[, indx_date := index_date]
  each_trial_PD <- each_trial[date.PD <= indx_date & date.PD.end >= indx_date, ]
  each_trial_PD_AF <- each_trial_PD[date.af <= indx_date]
  each_trial_PD_AF <- each_trial_PD_AF[date.af >= date.PD]
  # No prevalent-OAC exclusion (direction 2: compare current use vs not)
  id_with_prevalant_outcome <- each_trial_PD_AF[date.death <= indx_date |
    date.transplant <= indx_date |
    date.hd <= indx_date, Reference_Key]
  each_trial_PD_AF_NoOutcome <- each_trial_PD_AF[
    !Reference_Key %in% id_with_prevalant_outcome]
  rm(id_with_prevalant_outcome)
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx <- copy(each_trial_PD_AF_NoOutcome)
  if (VERBOSE_TRIAL_MSG) {
    message("Eligible N: ", nrow(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx))
  }
  # 2. exposure: current OAC prescription covering indx_date ----------------
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx <- assign_oac_expo(
    each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx,
    drug_oac,
    index_date
  )
  if (VERBOSE_TRIAL_MSG) {
    message("Expo: ", sum(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$expo),
            "\nNon-exp: ",
            each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx[, .N - sum(expo)])
  }
  # 3.covarates -------------------------------------------------------------
  each_trial_dx <- dx[Reference_Key %in% each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key & Reference_Date < index_date]
  
  px <- mapply(function(x) sapply(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key,function(y) getpx(y,x,each_trial_dx)),as.list(tibble::deframe(covar[,.(Name,Regex)])),USE.NAMES = T)
  
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- merge(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx,as.data.table(px,keep.rownames = T)[,setnames(.SD,"rn","Reference_Key")],by="Reference_Key")
  if (VERBOSE_TRIAL_MSG) {
    message("\nIncident: ", each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px[dx.stroke_embo == 0, .N])
  }
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px[dx.stroke_embo==0,]
  
  # 4. outcomes  ----------------------------------------------------------------
  each_trial_dx <- dx[Reference_Key %in% each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx$Reference_Key & Reference_Date >=index_date ]
  
  # 4.1 all stroke ----------------------------------------------------------
  
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,"^43[0-6]|^444","stroke")[]
  
  # 4.2 all bleeding ----------------------------------------------------------------
  #"^43[0-6]|^444"
  icd_ischemicstroke <- "^43[3-6]|^444"
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_ischemicstroke,"ische")[]
  
  icd_HaemorrhagicStroke <- "^43[0-2]"
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_HaemorrhagicStroke,"haem")[]
  
  icd_GI_Bleed <-  "^456.0^|456.20|^530.21|^530.7|^530.82|^531.[0246]|^532.[0246]|^533.[0246]|^534.[0246]|^537.8[34]|^562.0[23]|^562.1[23]|^569.3|^569.85|^578"
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_GI_Bleed,"GIbleed")[]
  
  icd_Other_Bleed <- "^287.[89]|^596.7|^784.8|^599.7|^627.1|^459.0|^719.1|^786.3|^363.6|^376.32|^377.42|^729.92|^423.0|^801.[2378]|^803.[2378]|^804.[2378]|^800.[2378]|^85[23]"
  each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px <- search_outcomes(each_trial_PD_AF_rmpastOAC_NoOutcome_NoPastHx_px,each_trial_dx,icd_Other_Bleed,"Otherbleed")[]
  
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

seqcohort <- pblapply(index_date[-1], run_seq, cohort)
saveRDS(seqcohort, "data/seqcohort_daily_current_oac_at_index.RDS")

# a <- pblapply(indx_date[136],run_seq,cohort)

# outcome definition ------------------------------------------------------
seqcohort <- lapply(seq_along(seqcohort), function(i) {
  df <- seqcohort[[i]]
  df$trial_id <- i  # Add the ID column with the index
  return(df)
})


person_trial <- rbindlist(seqcohort)
person_trial[,id:=paste(trial_id,Reference_Key,sep = "-")]
demo <- as.data.table(readRDS("data/PD_OAC-demo.RDS"))
person_trial <- merge(person_trial,
                      unique(demo[,.(Reference_Key,Sex,dob=Date_of_Birth_yyyymmdd)]))

setkey(person_trial,NULL)
person_trial[,id:=paste(trial_id,Reference_Key,sep = "_")]
# primary outcome - stroke ------------------------------------------------
person_trial[,obs_end:=pmin(ymd("2022-12-31"),
                            date.stroke,
                            date.death,
                            date.hd,
                            date.transplant,
                            indx_date %m+% years(5),
                            date.PD.end,na.rm = T)]

person_trial[,summary(as.numeric(obs_end-indx_date)/365.25)]


person_trial[,.N,expo]
# person_trial[,Reference_Key]

# ps 
person_trial[,age_indx:=as.numeric((indx_date-ymd(dob))/365.25)]
person_trial[,time_af_index:=as.numeric(indx_date-date.af)]
person_trial[,chadsvas_score:=dplyr::case_when(
  age_indx < 65 ~ 0,
  age_indx >= 65 & age_indx <75 ~ 1,
  age_indx > 75 ~ 2)+ifelse(Sex=="F",1,0)+dx.chf+dx.htn+dx.cbd+dx.stroke_embo*2+dx.pvd+dx.dm]
# ps.model <- glm(expo~Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+dx.ulcers+dx.stroke_embo+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation+time_af_index+chadsvas_score,family = binomial, data = person_trial)
# # Unstabilized weight
# person_trial$ps_score <- predict(ps.model, type = "response")
# person_trial$usweight <- with(person_trial, ifelse(expo==1,
#                                                    1/ps_score, 1/(1-ps_score)))
# round(summary(person_trial$usweight), 2)
# # Stabilized weight
# person_trial$sweight <- with(person_trial, ifelse(expo==1,
#                                                   mean(expo==1)/ps_score,
#                                                   (1-mean(expo==1))/(1-ps_score)))
# round(summary(person_trial$sweight), 2)
# # Truncating unstabilized weight
# person_trial <- person_trial %>%
#   mutate(usweight_t = pmin(pmax(usweight, quantile(usweight, 0.01)),
#                            quantile(usweight, 0.99)))
# summary(person_trial$usweight_t)
# # Truncating stabilized weight
# person_trial <- person_trial %>%
#   mutate(sweight_t = pmin(pmax(sweight, quantile(sweight, 0.01)),
#                           quantile(sweight, 0.99)))
# summary(person_trial$sweight_t)
# 
# labs <- paste("DOAC:", c("Users", "None users"))
# person_trial %>%
#   mutate(user = ifelse(expo == 1, labs[1], labs[2])) %>%
#   ggplot(aes(x = ps_score)) +
#   geom_histogram(color = "white") +
#   facet_wrap(~user) +
#   xlab("Probability of going to Catholic school") +
#   theme_bw()

# table one ---------------------------------------------------------------
person_trial[,fu:=as.numeric(obs_end-indx_date)]
person_trial <- merge(person_trial,drug_antiplatelet,by.x=c("Reference_Key","indx_date"),by.y=c("Reference_Key","Date"),all.x=T)
setnames(person_trial,"antiplatelet","antiplatelet.baseline")
person_trial[,antiplatelet.baseline:=ifelse(is.na(antiplatelet.baseline),0,antiplatelet.baseline)]
vars <- c("fu","Sex", "age_indx","time_af_index", "chadsvas_score","dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
          "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
          "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
          "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
          "dx.dm_com1", "dx.htn", "dx.constipation","antiplatelet.baseline"
)

library(labelled)

# 遍历covar表的每一行，为person_trial的对应列添加标签
for (i in 1:nrow(covar)) {
  var_name <- covar[i, Name]
  var_label <- covar[i, Description]
  
  # 如果变量存在，则设置标签
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

tabone <- CreateTableOne(vars = vars, 
                         factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                         strata = "expo", 
                         data = person_trial, test = F)

notnormal <- c("age_indx","time_af_index","fu")
tableone_b4_matching <- as.data.table(print(tabone, smd = T,varLabels=T,nonnormal = notnormal),keep.rownames = T)


# matching step -----------------------------------------------------------
set.seed(456)
m.out <- matchit(expo~Sex+age_indx+dx.chf+dx.mi+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+dx.dm+dx.crf+dx.liver_mild+dx.liver_modsev+dx.ulcers+dx.stroke_embo+dx.asthma+dx.ra+dx.aids+dx.cancer+dx.cancer_mets+dx.dm_com0+dx.dm_com1+dx.htn+dx.constipation+time_af_index+chadsvas_score,
                    data = person_trial,replace = T,method = "nearest",caliper = 0.2,ratio = 10)
summary(m.out)
plot(summary(m.out))
person_trial <- match.data(m.out)
tabone <- CreateTableOne(vars = vars, 
                         factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                         strata = "expo", data = person_trial, test = F)
tableone_after_matching <- as.data.table(print(tabone, smd = T,varLabels=T,nonnormal = notnormal),keep.rownames = T)
# writexl::write_xlsx(list(before=tableone_b4_matching,after=tableone_after_matching),"tableone.xlsx")

# for(i in 1:person_trial[expo==1,.N]){
#   cat(i,"\n")
#   tmp_cas <- person_trial[expo==1,][i]
#   tmp_ctl <- person_trial[,tmpdiff:=abs(ps_score - tmp_cas$ps_score)][order(tmpdiff)][expo==0][1:10,id]
#   tmp_ctl <- person_trial[id %in% tmp_ctl]
#   tmp_ctl[,tmpdiff:=NULL]
#   tmp <- rbind(tmp_cas,tmp_ctl)
#   tmp[,match:=i]
#   person_trial$tmpdiff <- NULL
#   if(i==1){
#     out <- tmp
#   }else{
#     out <- rbind(out,tmp)
#   }
# }
# 
# person_trial <- person_trial[id %in% out$id]

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

person_trial[,out:=ifelse(date.ische >= obs_date & substring(date.ische,1,7)== substring(obs_date,1,7) ,1,0)]
person_trial[,out:=ifelse(is.na(out),0,out)]

person_trial[,out:=ifelse(date.haem >= obs_date & substring(date.haem,1,7)== substring(obs_date,1,7) ,1,0)]
person_trial[,out:=ifelse(is.na(out),0,out)]

# 
# person_trial[,out:=ifelse(pmin(date.death,date.ische) >= obs_date & substring(date.death,1,7)== substring(obs_date,1,7) ,1,0)]
# person_trial[,out:=ifelse(is.na(out),0,out)]

# add antiplatete


person_trial <- merge(person_trial,drug_antiplatelet,by.x=c("Reference_Key","obs_date"),by.y=c("Reference_Key","Date"),all.x=T)
person_trial[,antiplatelet:=ifelse(is.na(antiplatelet),0,1)]

person_trial[,.N,expo]
person_trial[,.N,out]

# remove those records already have events
# person_trial <- person_trial[order(id,time),][,tempi:=cumsum(out),.(id)][tempi <= 1]



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

save_plot("number of eligible.pdf", test, width_in = 10, height_in =3)
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
save_plot("number of eligible by users or not.pdf", test, width_in = 10, height_in =3)


# description of sequential trail  ----------------------------------------

# Create a baseline dataset without duplicate IDs
no_dups <- person_trial[!duplicated(person_trial$Reference_Key),]

# Number of unique individuals 
length(no_dups$Reference_Key)

# Number of unique events
stroke <- person_trial[which(person_trial$out==1),]
no_dups_stroke <- stroke[!duplicated(stroke$Reference_Key),]
table(no_dups_stroke$out)

# Number of non-unique individuals
length(person_trial$Reference_Key)

# Number of non-unique events
table(person_trial$out[which(person_trial$fu==0)])

# Compute average number of trials each individual contributed to
avg_trials <- person_trial %>% dplyr::group_by(Reference_Key) %>%
  dplyr::summarise(trials_count = length(Reference_Key))

# Print average number of trials
summary(avg_trials$trials_count)


# table one ----------------------------------------------------------------
vars <- c("fu","Sex", "age_indx","time_af_index", "chadsvas_score","dx.chf", "dx.mi", "dx.pvd", "dx.cbd", "dx.copd",
          "dx.dementia", "dx.paralysis", "dx.dm", "dx.crf", "dx.liver_mild",
          "dx.liver_modsev", "dx.ulcers", "dx.stroke_embo", "dx.asthma",
          "dx.ra", "dx.aids", "dx.cancer", "dx.cancer_mets", "dx.dm_com0",
          "dx.dm_com1", "dx.htn", "dx.constipation","antiplatelet"
)
tabone <- CreateTableOne(vars = vars, 
                         factorVars = setdiff(vars,c("fu","trial_id","chadsvas_score","age_indx","time_af_index")),
                         strata = "expo", data = person_trial, test = F)
labelled::var_label(person_trial) 
notnormal <- c("age_indx","time_af_index","fu")
print(tabone, smd = T,varLabels=T,nonnormal = notnormal)

unbalanced_variable <- as.data.table(as.data.frame(print(tabone,smd=T)),keep.rownames = T)[as.numeric(SMD) > 0.1,gsub(" = [1M] \\(%\\)| \\(mean \\(SD\\)\\)","",rn)]
unbalanced_variable <- setdiff(unbalanced_variable,c("trial_id"))
formula_rr_main <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(unbalanced_variable,collapse = "+")))
formula_rr_sex <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(setdiff(unbalanced_variable,"Sex"),collapse = "+")))
formula_rr_inci <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(setdiff(unbalanced_variable,"dx.stroke_embo"),collapse = "+")))
formula_rr_sex_inci <- formula(paste("out~expo+I(time)+I(timesqr)+",paste0(setdiff(unbalanced_variable,c("Sex","dx.stroke_embo")),collapse = "+")))

library(speedglm)
library(splitstackshape)
library(boot)


# models for risk ratio (pooled logistic as HR)  ------------------------------------------------------------------
# only those without balancing 

# Fit pooled logistic regression model with covariates -> time timesqr -> average HR druign the whole up
model.primary <- tidy(speedglm(formula_rr_main , data = person_trial,family =binomial(link = 'logit')),exponentiate = T,conf.int = T) %>% mutate(model="pri")
model.F <- tidy(speedglm(formula_rr_sex , data = person_trial[Sex=="F"],family =binomial(link = 'logit')),exponentiate = T,conf.int = T)%>% mutate(model="F")
model.M <- tidy(speedglm(formula_rr_sex, data = person_trial[Sex=="M"],family =binomial(link = 'logit')),exponentiate = T,conf.int = T)%>% mutate(model="M")
model.inci <- tidy(speedglm(formula_rr_inci , data = person_trial[dx.stroke_embo==0],family =binomial(link = 'logit')),exponentiate = T,conf.int = T)%>% mutate(model="inci")
# model.prev <- tidy(speedglm(formula_rr_inci , data = person_trial[dx.stroke_embo==1],family =binomial(link = 'logit')),exponentiate = T,conf.int = T)%>% mutate(model="prev")
model.F.inci <- tidy(speedglm(formula_rr_sex_inci , data = person_trial[Sex=="F" & dx.stroke_embo==0],family =binomial(link = 'logit')),exponentiate = T,conf.int = T)%>% mutate(model="F_inc")
model.M.inci <- tidy(speedglm(formula_rr_sex_inci , data = person_trial[Sex=="M" & dx.stroke_embo==0],family =binomial(link = 'logit')),exponentiate = T,conf.int = T)%>% mutate(model="M_inc")

all_pool_logi <- do.call(rbind, list(model.primary,model.F, model.M, model.inci,model.F.inci, model.M.inci))
rbindlist(lapply(unique(all_pool_logi$model),
                 function(x) all_pool_logi %>% filter(term=="expo" & model==x))) -> est_result_noboostrap
# models for risk difference (pooled logistic as HR)  ------------------------------------------------------------------

{
current_bootstrap <- 0
elig_ids <- data.frame(id = unique(person_trial$id))
std.boot <- function(data,indices){
  current_bootstrap <<- current_bootstrap + 1
  if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
  # Select individuals into each bootstrapped sample
  ids <- data$id 
  boot.ids <- data.frame(id = ids[indices]) 
  # boot.ids <- data.frame(id = ids[sample(length(ids))])
  boot.ids$bid <- 1:nrow(boot.ids) 
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
  d$bid_new <- interaction(d$bid, d$trial_id) 
  
  { # for risk difference
    # Fit pooled logistic model to estimate discrete hazards
    fit.pool1 <- tryCatch({
      fit.pool1 <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
                             chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                             dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                           family=binomial(link="logit"),data=d)
    }, error = function(e) {
      message("Risk difference model failed: ", e$message)
      return(NULL)
    })
    if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量

    # Create dataset with all time points for each individual under each treatment level
    treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
    treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
    treat0$timesqr <- treat0$time^2
    
    # Create "treat_b" variable under no baseline vaccination
    treat0$expo <- 0
    
    # Create "treat_b" variable under baseline CROWN vaccination
    treat1 <- treat0
    treat1$expo <- 1
    
    # Extract predicted values from pooled logistic regression model for each person-time row
    # Predicted values correspond to discrete-time hazards
    treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
    treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
    # The above creates a person-time dataset where we have predicted discrete-time hazards
    # For each person-time row in the dataset
    
    # Obtain predicted survival probabilities from discrete-time hazards
    treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
    treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
    
    # Estimate risks from survival probabilities
    # Risk = 1 - S(t)
    treat0.surv$risk0 <- 1 - treat0.surv$surv0
    treat1.surv$risk1 <- 1 - treat1.surv$surv1
    
    # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
    risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
    risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
    
    # Prepare data
    graph.pred <- merge(risk0, risk1, by=c("time"))
    # Edit data frame to reflect that risks are estimated at the END of each interval
    graph.pred$time_0 <- graph.pred$time + 1
    zero <- data.frame(cbind(0,0,0,1,0,0))
    zero <- setNames(zero,names(graph.pred))
    graph <- rbind(zero, graph.pred)
    
    graph$rd <- graph$risk1-graph$risk0
    graph$rr <- graph$risk1/graph$risk0
    rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
  }
  { # for risk ratio
    # Fit pooled logistic model to estimate discrete hazards
    fit.pool1 <- tryCatch({
      fit.pool1 <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
                             chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                             dx.ra + dx.dm_com0 + dx.htn + antiplatelet,
                           family=binomial(link="logit"),data=d)
    }, error = function(e) {
      message("Risk difference model failed: ", e$message)
      return(NULL)
    })
    if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
    
    # Create dataset with all time points for each individual under each treatment level
    treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
    treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
    treat0$timesqr <- treat0$time^2
    
    # Create "treat_b" variable under no baseline vaccination
    treat0$expo <- 0
    
    # Create "treat_b" variable under baseline CROWN vaccination
    treat1 <- treat0
    treat1$expo <- 1
    
    # Extract predicted values from pooled logistic regression model for each person-time row
    # Predicted values correspond to discrete-time hazards
    treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
    treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
    # The above creates a person-time dataset where we have predicted discrete-time hazards
    # For each person-time row in the dataset
    
    # Obtain predicted survival probabilities from discrete-time hazards
    treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
    treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
    
    # Estimate risks from survival probabilities
    # Risk = 1 - S(t)
    treat0.surv$risk0 <- 1 - treat0.surv$surv0
    treat1.surv$risk1 <- 1 - treat1.surv$surv1
    
    # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
    risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
    risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
    
    # Prepare data
    graph.pred <- merge(risk0, risk1, by=c("time"))
    # Edit data frame to reflect that risks are estimated at the END of each interval
    graph.pred$time_0 <- graph.pred$time + 1
    zero <- data.frame(cbind(0,0,0,1,0,0))
    zero <- setNames(zero,names(graph.pred))
    graph1 <- rbind(zero, graph.pred)
    
    graph1$rd <- graph1$risk1-graph1$risk0
    graph1$rr <- graph1$risk1/graph1$risk0
  }
  return(c(graph$risk0[which(graph$time==59)],
           graph$risk1[which(graph$time==59)],
           graph$rd[which(graph$time==59)],
           graph1$rr[which(graph$time==59)]))
}


# Run bootstrap samples
set.seed(456)
risk.results <- boot(data = elig_ids,
                     statistic = std.boot,
                     R=500)

est_result <- 
as.data.frame(matrix(c(risk.results$t0[1],
                       # 95% CI for risk in no drug arm
                       boot.ci(risk.results,
                               conf = 0.95, 
                               type = "perc", 
                               index = 1)$percent[c(4,5)],
                       # 95% CI for risk in drug arm
                       risk.results$t0[2],
                       boot.ci(risk.results,
                               conf = 0.95, 
                               type = "perc",
                               index = 2)$percent[c(4,5)],
                       # 95% CI for risk difference
                       risk.results$t0[3],
                       boot.ci(risk.results,
                               conf = 0.95, 
                               type = "perc", 
                               index = 3)$percent[c(4,5)],
                       # 95% CI for risk ratio
                       risk.results$t0[4],
                       boot.ci(risk.results,
                               conf = 0.95,
                               type = "perc", 
                               index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
setnames(est_result,c("est",'lower','upper'))
rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
est_result$model <- "full"
main_est_result <- est_result
}

saveRDS(risk.results,"data/main_boostrap_500.rds")
# est   lower  upper model
# risk_nonuser    0.1893  0.1036 0.4161  full
# risk_user       0.3284  0.1380 0.7022  full
# risk difference 0.1391 -0.1521 0.4490  full
# risk ratio      2.4875  1.3861 3.9819  full


# hemorraghic stroke
# est    lower  upper model
# risk_nonuser     0.15508  0.06895 0.4119  full
# risk_user        0.13462  0.04987 0.2467  full
# risk difference -0.02046 -0.30054 0.1026  full
# risk ratio       1.89957  0.79947 3.2877  full


# est    lower  upper model
# risk_nonuser     0.32350  0.14231 0.8382  full
# risk_user        0.22905  0.07747 0.6959  full
# risk difference -0.09445 -0.54477 0.2104  full
# risk ratio       1.33567  0.55664 2.4893  full



#subgroup by sex

{
  current_bootstrap <- 0
  elig_ids <- data.frame(id = person_trial[Sex=="M",unique(id)])
  std.boot <- function(data,indices){
    current_bootstrap <<- current_bootstrap + 1
    if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
    # Select individuals into each bootstrapped sample
    ids <- data$id 
    boot.ids <- data.frame(id = ids[indices]) 
    # boot.ids <- data.frame(Reference_Key = ids[sample(length(ids))])
    boot.ids$bid <- 1:nrow(boot.ids) 
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
    d$bid_new <- interaction(d$bid, d$trial_id) 
    
    { # for risk difference
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Risk difference model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph <- rbind(zero, graph.pred)
      
      graph$rd <- graph$risk1-graph$risk0
      graph$rr <- graph$risk1/graph$risk0
      rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
    }
    { # for risk ratio
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  ,
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Risk difference model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph1 <- rbind(zero, graph.pred)
      
      graph1$rd <- graph1$risk1-graph1$risk0
      graph1$rr <- graph1$risk1/graph1$risk0
    }
    return(c(graph$risk0[which(graph$time==59)],
             graph$risk1[which(graph$time==59)],
             graph$rd[which(graph$time==59)],
             graph1$rr[which(graph$time==59)]))
  }
  
  
  # Run bootstrap samples
  set.seed(456)
  risk.results <- boot(data = elig_ids,
                       statistic = std.boot,
                       R=500)
  est_result <- 
    as.data.frame(matrix(c(risk.results$t0[1],
                           # 95% CI for risk in no drug arm
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 1)$percent[c(4,5)],
                           # 95% CI for risk in drug arm
                           risk.results$t0[2],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc",
                                   index = 2)$percent[c(4,5)],
                           # 95% CI for risk difference
                           risk.results$t0[3],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 3)$percent[c(4,5)],
                           # 95% CI for risk ratio
                           risk.results$t0[4],
                           boot.ci(risk.results,
                                   conf = 0.95,
                                   type = "perc", 
                                   index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
  setnames(est_result,c("est",'lower','upper'))
  rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
  est_result$model <- "male"
  male_est_result <- est_result
}

saveRDS(risk.results,"data/male_boostrap_500.rds")
# est    lower  upper model
# risk_nonuser    0.08905  0.03652 0.3857  male
# risk_user       0.16801  0.06232 0.9874  male
# risk difference 0.07896 -0.16969 0.8459  male
# risk ratio      2.58281  0.97321 4.7756  male

# est          lower  upper model
# risk_nonuser     0.14740  0.07274576411 0.3311  male
# risk_user        0.11734  0.00000000921 1.0000  male
# risk difference -0.03006 -0.23548442492 0.8561  male
# risk ratio       1.09352  0.00000006082 2.8954  male


# Female
{
  current_bootstrap <- 0
  elig_ids <- data.frame(id = person_trial[Sex=="F",unique(id)])
  std.boot <- function(data,indices){
    current_bootstrap <<- current_bootstrap + 1
    if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
    # Select individuals into each bootstrapped sample
    ids <- data$id 
    boot.ids <- data.frame(id = ids[indices]) 
    # boot.ids <- data.frame(Reference_Key = ids[sample(length(ids))])
    boot.ids$bid <- 1:nrow(boot.ids) 
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
    d$bid_new <- interaction(d$bid, d$trial_id) 
    
    { # for risk difference
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Risk difference model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph <- rbind(zero, graph.pred)
      
      graph$rd <- graph$risk1-graph$risk0
      graph$rr <- graph$risk1/graph$risk0
      rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
    }
    { # for risk ratio
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        fit.pool1 <-speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                               chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                               dx.ra + dx.dm_com0 + dx.htn + antiplatelet,
                             family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Risk difference model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
  
     # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph1 <- rbind(zero, graph.pred)
      
      graph1$rd <- graph1$risk1-graph1$risk0
      graph1$rr <- graph1$risk1/graph1$risk0
    }
    return(c(graph$risk0[which(graph$time==59)],
             graph$risk1[which(graph$time==59)],
             graph$rd[which(graph$time==59)],
             graph1$rr[which(graph$time==59)]))
  }
  
  
  # Run bootstrap samples
  set.seed(456)
  risk.results <- boot(data = elig_ids,
                       statistic = std.boot,
                       R=500)
  est_result <- 
    as.data.frame(matrix(c(risk.results$t0[1],
                           # 95% CI for risk in no drug arm
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 1)$percent[c(4,5)],
                           # 95% CI for risk in drug arm
                           risk.results$t0[2],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc",
                                   index = 2)$percent[c(4,5)],
                           # 95% CI for risk difference
                           risk.results$t0[3],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 3)$percent[c(4,5)],
                           # 95% CI for risk ratio
                           risk.results$t0[4],
                           boot.ci(risk.results,
                                   conf = 0.95,
                                   type = "perc", 
                                   index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
  setnames(est_result,c("est",'lower','upper'))
  rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
  est_result$model <- "female"
  female_est_result <- est_result
  
}

saveRDS(risk.results,"data/female_boostrap_500.rds")
# est   lower  upper  model
# risk_nonuser    0.2302  0.1091 0.4928 female
# risk_user       0.6748  0.1048 0.9637 female
# risk difference 0.4445 -0.2266 0.7643 female
# risk ratio      2.4960  1.1415 4.3685 female

# est    lower  upper  model
# risk_nonuser    0.12256  0.02331 0.7245 female
# risk_user       0.14699  0.04941 0.3501 female
# risk difference 0.02443 -0.54662 0.2223 female
# risk ratio      3.53141  1.02318 9.8051 female

# incident cases
{
  current_bootstrap <- 0
  elig_ids <- data.frame(id = person_trial[dx.stroke_embo==0,unique(id)])
  std.boot <- function(data,indices){
    current_bootstrap <<- current_bootstrap + 1
    if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
    # Select individuals into each bootstrapped sample
    ids <- data$id 
    boot.ids <- data.frame(id = ids[indices]) 
    # boot.ids <- data.frame(Reference_Key = ids[sample(length(ids))])
    boot.ids$bid <- 1:nrow(boot.ids) 
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
    d$bid_new <- interaction(d$bid, d$trial_id) 
    
    { # for risk difference
      # Fit pooled logistic model to estimate discrete hazards
      
      fit.pool1 <- tryCatch({
        fit.pool1 <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
                               chadsvas_score + dx.cbd + dx.copd + dx.dementia  + 
                               dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                             family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Risk difference model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph <- rbind(zero, graph.pred)
      
      graph$rd <- graph$risk1-graph$risk0
      graph$rr <- graph$risk1/graph$risk0
      rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
    }
    { # for risk ratio
      # Fit pooled logistic model to estimate discrete hazards
      
      fit.pool1 <- tryCatch({
        fit.pool1 <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
                               chadsvas_score + dx.cbd + dx.copd + dx.dementia  + 
                               dx.ra + dx.dm_com0 + dx.htn + antiplatelet ,
                             family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Risk difference model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph1 <- rbind(zero, graph.pred)
      
      graph1$rd <- graph1$risk1-graph1$risk0
      graph1$rr <- graph1$risk1/graph1$risk0
    }
    return(c(graph$risk0[which(graph$time==59)],
             graph$risk1[which(graph$time==59)],
             graph$rd[which(graph$time==59)],
             graph1$rr[which(graph$time==59)]))
  }
  
  
  # Run bootstrap samples
  set.seed(456)
  risk.results <- boot(data = elig_ids,
                       statistic = std.boot,
                       R=500)
  est_result <- 
    as.data.frame(matrix(c(risk.results$t0[1],
                           # 95% CI for risk in no drug arm
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 1)$percent[c(4,5)],
                           # 95% CI for risk in drug arm
                           risk.results$t0[2],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc",
                                   index = 2)$percent[c(4,5)],
                           # 95% CI for risk difference
                           risk.results$t0[3],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 3)$percent[c(4,5)],
                           # 95% CI for risk ratio
                           risk.results$t0[4],
                           boot.ci(risk.results,
                                   conf = 0.95,
                                   type = "perc", 
                                   index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
  setnames(est_result,c("est",'lower','upper'))
  rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
  est_result$model <- "inci"
  inci_est_result <- est_result
  
  saveRDS(risk.results,"data/inci_boostrap_500.rds")
}
# est    lower  upper model
# risk_nonuser    0.12036  0.06108 0.3292  inci
# risk_user       0.18281  0.06998 0.4159  inci
# risk difference 0.06245 -0.14327 0.2626  inci
# risk ratio      1.93670  0.75203 3.6133  inci

# est    lower  upper model
# risk_nonuser     0.14247  0.05317 0.5619  inci
# risk_user        0.11859  0.04746 0.2671  inci
# risk difference -0.02388 -0.39273 0.1122  inci
# risk ratio       1.84130  0.77384 3.8093  inci




# recurrent cases
{
  current_bootstrap <- 0
  elig_ids <- data.frame(id = person_trial[dx.stroke_embo==1,unique(id)])
  std.boot <- function(data,indices){
    current_bootstrap <<- current_bootstrap + 1
    if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
    # Select individuals into each bootstrapped sample
    ids <- data$id 
    boot.ids <- data.frame(id = ids[indices]) 
    # boot.ids <- data.frame(id = ids[sample(length(ids))])
    boot.ids$bid <- 1:nrow(boot.ids) 
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
    d$bid_new <- interaction(d$bid, d$trial_id) 
    
    { # for risk difference
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + Sex+ age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia +  
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph <- rbind(zero, graph.pred)
      
      graph$rd <- graph$risk1-graph$risk0
      graph$rr <- graph$risk1/graph$risk0
      rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
    }
    { # for risk ratio
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + Sex+ age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia +  
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet,
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph1 <- rbind(zero, graph.pred)
      
      graph1$rd <- graph1$risk1-graph1$risk0
      graph1$rr <- graph1$risk1/graph1$risk0
    }
    return(c(graph$risk0[which(graph$time==59)],
             graph$risk1[which(graph$time==59)],
             graph$rd[which(graph$time==59)],
             graph1$rr[which(graph$time==59)]))
  }
  
  # Run bootstrap samples
  set.seed(456)
  risk.results <- boot(data = elig_ids,
                       statistic = std.boot,
                       R=500)
  est_result <- 
    as.data.frame(matrix(c(risk.results$t0[1],
                           # 95% CI for risk in no drug arm
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 1)$percent[c(4,5)],
                           # 95% CI for risk in drug arm
                           risk.results$t0[2],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc",
                                   index = 2)$percent[c(4,5)],
                           # 95% CI for risk difference
                           risk.results$t0[3],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 3)$percent[c(4,5)],
                           # 95% CI for risk ratio
                           risk.results$t0[4],
                           boot.ci(risk.results,
                                   conf = 0.95,
                                   type = "perc", 
                                   index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
  setnames(est_result,c("est",'lower','upper'))
  rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
  est_result$model <- "rec"
  rec_est_result <- est_result
  saveRDS(risk.results,"data/recurrent_boostrap_500.rds")
}
# est   lower  upper model
# risk_nonuser     NA  0.1399 0.8589   rec
# risk_user        NA  0.2064 1.0000   rec
# risk difference  NA -0.6141 0.8109   rec
# risk ratio       NA  1.1272 6.0167   rec


# incident cases -male
{
  current_bootstrap <- 0
  elig_ids <- data.frame(id = person_trial[dx.stroke_embo==0 & Sex =="M",unique(id)])
  std.boot <- function(data,indices){
    current_bootstrap <<- current_bootstrap + 1
    if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
    # Select individuals into each bootstrapped sample
    ids <- data$id 
    boot.ids <- data.frame(id = ids[indices]) 
    # boot.ids <- data.frame(Reference_Key = ids[sample(length(ids))])
    boot.ids$bid <- 1:nrow(boot.ids) 
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
    d$bid_new <- interaction(d$bid, d$trial_id) 
    
    { # for risk difference
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia +  
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph <- rbind(zero, graph.pred)
      
      graph$rd <- graph$risk1-graph$risk0
      graph$rr <- graph$risk1/graph$risk0
      rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
    }
    { # for risk ratio
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia +  
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph1 <- rbind(zero, graph.pred)
      
      graph1$rd <- graph1$risk1-graph1$risk0
      graph1$rr <- graph1$risk1/graph1$risk0
    }
    return(c(graph$risk0[which(graph$time==59)],
             graph$risk1[which(graph$time==59)],
             graph$rd[which(graph$time==59)],
             graph1$rr[which(graph$time==59)]))
  }
  
  
  # Run bootstrap samples
  set.seed(456)
  risk.results <- boot(data = elig_ids,
                       statistic = std.boot,
                       R=500)
  est_result <- 
    as.data.frame(matrix(c(risk.results$t0[1],
                           # 95% CI for risk in no drug arm
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 1)$percent[c(4,5)],
                           # 95% CI for risk in drug arm
                           risk.results$t0[2],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc",
                                   index = 2)$percent[c(4,5)],
                           # 95% CI for risk difference
                           risk.results$t0[3],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 3)$percent[c(4,5)],
                           # 95% CI for risk ratio
                           risk.results$t0[4],
                           boot.ci(risk.results,
                                   conf = 0.95,
                                   type = "perc", 
                                   index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
  setnames(est_result,c("est",'lower','upper'))
  rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
  est_result$model <- "inc_male"
  inc_male_est_result <- est_result
  
  saveRDS(risk.results,"data/inci_male_boostrap_500.rds")
}

# est          lower   upper    model
# risk_nonuser    0.03078  0.01216068040 0.08693 inc_male
# risk_user       0.04970  0.00000001745 0.12697 inc_male
# risk difference 0.01892 -0.05891922943 0.08272 inc_male
# risk ratio      1.61470  0.00000039026 3.83130 inc_male
# 
# est           lower  upper    model
# risk_nonuser     0.10661  0.041221797120 0.3585 inc_male
# risk_user        0.09372  0.000000002888 0.3249 inc_male
# risk difference -0.01288 -0.240475924188 0.2057 inc_male
# risk ratio       0.87916  0.000000040910 3.6886 inc_male



# incident cases -female
{
  current_bootstrap <- 0
  elig_ids <- data.frame(id = person_trial[dx.stroke_embo==0 & Sex =="F",unique(id)])
  std.boot <- function(data,indices){
    current_bootstrap <<- current_bootstrap + 1
    if(current_bootstrap %% 100 ==0){cat("Processing Bootstrap sample:", current_bootstrap, "\n")}
    # Select individuals into each bootstrapped sample
    ids <- data$id 
    boot.ids <- data.frame(id = ids[indices]) 
    # boot.ids <- data.frame(Reference_Key = ids[sample(length(ids))])
    boot.ids$bid <- 1:nrow(boot.ids) 
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- left_join(boot.ids, person_trial, by="id",relationship="many-to-many") 
    d$bid_new <- interaction(d$bid, d$trial_id) 
    
    { # for risk difference
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia +  
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph <- rbind(zero, graph.pred)
      
      graph$rd <- graph$risk1-graph$risk0
      graph$rr <- graph$risk1/graph$risk0
      rm(fit.pool1,treat0,treat1,risk0,risk1,treat0.surv,treat1.surv,graph.pred,zero)
    }
    { # for risk ratio
      # Fit pooled logistic model to estimate discrete hazards
      fit.pool1 <- tryCatch({
        speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                   chadsvas_score + dx.cbd + dx.copd + dx.dementia +  
                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet  + I(time*expo) + I(timesqr*expo),
                 family=binomial(link="logit"),data=d)
      }, error = function(e) {
        message("Model failed: ", e$message)
        return(NULL)
      })
      if (is.null(fit.pool1)) return(rep(NA, 4))  # 返回 NA 向量
      # Create dataset with all time points for each individual under each treatment level
      treat0 <- expandRows(d[which(d$time==0),], count=5*12, count.is.col=F) 
      treat0$time <- rep(seq(0, 59), nrow(d[which(d$time==0),]))
      treat0$timesqr <- treat0$time^2
      
      # Create "treat_b" variable under no baseline vaccination
      treat0$expo <- 0
      
      # Create "treat_b" variable under baseline CROWN vaccination
      treat1 <- treat0
      treat1$expo <- 1
      
      # Extract predicted values from pooled logistic regression model for each person-time row
      # Predicted values correspond to discrete-time hazards
      treat0$p.event0 <- predict(fit.pool1, treat0, type="response")
      treat1$p.event1 <- predict(fit.pool1, treat1, type="response")
      # The above creates a person-time dataset where we have predicted discrete-time hazards
      # For each person-time row in the dataset
      
      # Obtain predicted survival probabilities from discrete-time hazards
      treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup() 
      treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup() 
      
      # Estimate risks from survival probabilities
      # Risk = 1 - S(t)
      treat0.surv$risk0 <- 1 - treat0.surv$surv0
      treat1.surv$risk1 <- 1 - treat1.surv$surv1
      
      # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
      risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
      risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
      
      # Prepare data
      graph.pred <- merge(risk0, risk1, by=c("time"))
      # Edit data frame to reflect that risks are estimated at the END of each interval
      graph.pred$time_0 <- graph.pred$time + 1
      zero <- data.frame(cbind(0,0,0,1,0,0))
      zero <- setNames(zero,names(graph.pred))
      graph1 <- rbind(zero, graph.pred)
      
      graph1$rd <- graph1$risk1-graph1$risk0
      graph1$rr <- graph1$risk1/graph1$risk0
    }
    return(c(graph$risk0[which(graph$time==59)],
             graph$risk1[which(graph$time==59)],
             graph$rd[which(graph$time==59)],
             graph1$rr[which(graph$time==59)]))
  }
  
  
  # Run bootstrap samples
  set.seed(456)
  risk.results <- boot(data = elig_ids,
                       statistic = std.boot,
                       R=500)
  
  est_result <- 
    as.data.frame(matrix(c(risk.results$t0[1],
                           # 95% CI for risk in no drug arm
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 1)$percent[c(4,5)],
                           # 95% CI for risk in drug arm
                           risk.results$t0[2],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc",
                                   index = 2)$percent[c(4,5)],
                           # 95% CI for risk difference
                           risk.results$t0[3],
                           boot.ci(risk.results,
                                   conf = 0.95, 
                                   type = "perc", 
                                   index = 3)$percent[c(4,5)],
                           # 95% CI for risk ratio
                           risk.results$t0[4],
                           boot.ci(risk.results,
                                   conf = 0.95,
                                   type = "perc", 
                                   index = 4)$percent[c(4,5)]),ncol=3,byrow = T))
  setnames(est_result,c("est",'lower','upper'))
  rownames(est_result) <- c("risk_nonuser","risk_user","risk difference","risk ratio")
  est_result$model <- "inc_female"
  inc_female_est_result <- est_result
  
  saveRDS(risk.results,"data/inci_female_boostrap_500.rds")
}
# est    lower  upper      model
# risk_nonuser    0.1868  0.07846 0.5085 inc_female
# risk_user       0.5965  0.06492 1.0000 inc_female
# risk difference 0.4098 -0.26695 0.8003 inc_female
# risk ratio      3.1942  0.19238 8.2380 inc_female

# est     lower   upper      model
# risk_nonuser     0.153252  0.006678  0.8076 inc_female
# risk_user        0.148166  0.035725  0.6317 inc_female
# risk difference -0.005086 -0.485176  0.3056 inc_female
# risk ratio       0.966811  0.164153 22.1077 inc_female


writexl::write_xlsx(do.call(rbind,list(main_est_result,male_est_result,female_est_result,inci_est_result,inc_male_est_result,inc_female_est_result)),
      "hemo_stroke_results.xlsx")


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
