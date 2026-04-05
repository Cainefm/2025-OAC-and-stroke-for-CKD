# 2025-OAC-and-stroke-for-CKD

## Oral anticoagulants in PD patients with atrial fibrillation (target trial emulation)

## Relation to `script/run.R`

`script/run.R` was **not** modified. As it stands, it implements **daily** index dates and **current OAC** exposure (prescription intervals on the index day). That **does not** match the submitted Word manuscript, which describes **191 monthly** sequential trials and **OAC initiation within one month** before index (initiators vs non-initiators), with numbers such as 8,230 eligible patient-trials and 81 initiators.

The new pipeline below defaults to **`manuscript_initiator` + `month`** so the **methods align with the Word file**. Exact numerical reproduction also depends on the same input `.RDS` files and on analysis choices (for example follow-up cap).

## New pipeline (recommended)

| File | Role |
|------|------|
| `script/analysis_config.R` | `analysis_config` list + `load_analysis_config()` |
| `script/00_functions_pd_oac.R` | Index dates, drug/OAC helpers, `make_run_seq()` |
| `script/01_build_seqcohort.R` | Build and save sequential trial list |
| `script/02_analysis_primary.R` | `person_trial`, Table 1, 1:10 MatchIt, CSV + RDS outputs |

From the project root (where `data/` and `documents/` live):

```bash
Rscript script/01_build_seqcohort.R
Rscript script/02_analysis_primary.R
```

Outputs (by default):

- `data/seqcohort_manuscript_initiator_month.RDS` — trial list  
- `data/seqcohort_manuscript_initiator_month_antiplatelet.RDS` — antiplatelet table for merging  
- `out/tableone_before_matching.csv`, `out/tableone_after_matching.csv`  
- `out/person_trial_primary.rds` — matched cohort for downstream models  

## Configuration

Edit entries in `script/analysis_config.R` before running, or override in R:

```r
source("script/analysis_config.R")
analysis_config$exposure_mode <- "current_oac_interval"  # Rx-based current use
analysis_config$index_by <- "day"
analysis_config$path_seqcohort_rds <- NULL  # auto name from mode + index_by
cfg <- load_analysis_config()
```

- **`manuscript_initiator`**: prevalent OAC exclusion (`date.oac` ≤ index − 1 month); `expo` from first OAC in (index − 1 month, index]; uses cohort `date.oac`.  
- **`current_oac_interval`**: no prevalent exclusion; `expo` from OAC prescription intervals (`prepare_drug_oac()`). Matches the **logic** of the current `script/run.R` when index is daily.

`obs_end_cap_years` defaults to **3** in `analysis_config` to align with the manuscript’s three-year risk description. `script/run.R` uses **five** years in `obs_end`; set `analysis_config$obs_end_cap_years <- 5` to mirror that script’s follow-up cap.

## Dependencies

R packages used in 01/02 include: `data.table`, `dplyr`, `lubridate`, `pbapply`, `readxl`, `tibble`, `tableone`, `labelled`, `MatchIt`. Install as needed.

## Git

This folder may be initialized as a git repository separately; the peritonitis project uses `renv` for locking packages—consider `renv` here if you need strict reproducibility.
