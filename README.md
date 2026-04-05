# 2025-OAC-and-stroke-for-CKD

## Oral anticoagulants in PD patients with atrial fibrillation (target trial emulation)

## Relation to `script/run.R`

`script/run.R` was **not** modified. As it stands, it implements **daily** index dates and **current OAC** exposure (prescription intervals on the index day). That **does not** match the submitted Word manuscript, which describes **191 monthly** sequential trials and **OAC initiation within one month** before index (initiators vs non-initiators), with numbers such as 8,230 eligible patient-trials and 81 initiators.

The pipeline defaults to **`manuscript_initiator`** exposure and **`index_by = "day"`** (daily sequential trials and daily long person-time in 03). Set **`index_by = "month"`** to match the Word manuscript’s ~191 monthly trials and **`script/run_incident.R`** month-based person-time. Exact counts depend on input `.RDS` files and choices such as follow-up cap.

## Sequential pipeline (01–03)

Run each step from the **project root** (where `data/` and `documents/` live). Scripts call `project_root()` after sourcing `script/00_functions_pd_oac.R`; the working directory must be the root so `data/` is found.

| File | Role |
|------|------|
| `script/analysis_config.R` | `analysis_config` list + `load_analysis_config()` |
| `script/00_functions_pd_oac.R` | Shared helpers: index dates, OAC/antiplatelet prep, `make_run_seq()`, CHA₂DS₂-VASc (`compute_chadsvas_score()`), baseline antiplatelet merge, TableOne factor/SMD helpers, `project_root()` |
| `script/01_build_seqcohort.R` | Build and save the sequential trial list + cached antiplatelet at index (month- or day-level per `index_by`) |
| `script/02_analysis_primary.R` | Wide `person_trial`, Table 1, MatchIt (1:10 by default), CSV + RDS outputs |
| `script/03_regression_pooled_logistic.R` | Long person-time (`index_by`: day or month), time-varying antiplatelet, SMD-adjusted `speedglm` (optional), cluster bootstrap RD/RR (optional) |
| `script/bootstrap_pooled_run_incident.R` | Sourced by 03: `std.boot`, `calculate_risk`, `orgnize_ci` (same logic as `script/run_incident.R` ~51–158) |

```bash
Rscript script/01_build_seqcohort.R
Rscript script/02_analysis_primary.R
Rscript script/03_regression_pooled_logistic.R
```

Default outputs (paths can be changed in `script/analysis_config.R`):

- `data/seqcohort_<exposure_mode>_<index_by>.RDS` — trial list (name derived when `path_seqcohort_rds` is `NULL`)
- `data/seqcohort_<...>_antiplatelet.RDS` — antiplatelet table for merging at index
- `out/tableone_before_matching.csv`, `out/tableone_after_matching.csv`
- `out/person_trial_primary.rds` — matched cohort for 03
- `out/est_pooled_logistic_expo.csv`, `out/est_pooled_logistic_all_terms.csv` — **`speedglm`** coefficients (SMD-adjusted; when `regression_run_speedglm = TRUE`)
- `out/bootstrap_pooled_risk_run_incident.csv` — **bootstrap** risk estimates (main, male, female, younger, elder) aligned with **`script/run_incident.R`** (~700–772), when `pipeline_style = "run_incident"` and `regression_run_bootstrap = TRUE`
- `out/tableone_long_format_for_smd.csv` — long-format TableOne with SMD (written when `regression_run_speedglm = TRUE` and `pipeline_style = "run_incident"`)

## Configuration

Edit entries in `script/analysis_config.R` before running, or override in R:

```r
source("script/analysis_config.R")
analysis_config$exposure_mode <- "current_oac_interval"  # Rx-based current use
analysis_config$index_by <- "day"
analysis_config$path_seqcohort_rds <- NULL  # auto name from mode + index_by
cfg <- load_analysis_config()
```

- **`manuscript_initiator`**: prevalent OAC exclusion (`date.oac` ≤ index − 1 month); `expo` from first OAC in (index − `oac_initiation_lookback_months`, index]; uses cohort `date.oac`.
- **`current_oac_interval`**: no prevalent exclusion; `expo` from OAC prescription intervals (`prepare_drug_oac()`). Matches the **logic** of the current `script/run.R` when index is daily.

**`index_by`** (`"day"` default, or `"month"`):

- **`day`**: `01` builds one trial per calendar day between `study_start` and `study_end` (after `drop_first_index`); `03` expands each matched person-trial to **one row per day** from index to `obs_end`. **Much larger** row counts and runtime than monthly.
- **`month`**: one trial per month (~191 after dropping the first index); `03` uses **one row per calendar month** (aligned with **`script/run_incident.R`** person-time).

For **`index_by = "day"`**, step **03** (especially **`regression_run_bootstrap`**) can be **very slow**; reduce **`bootstrap_r`** while testing.

**`pipeline_style`** (in `analysis_config`):

- **`default`**: MatchIt includes `dx.stroke_embo` in the PS, uses caliper and `replace` as configured.
- **`run_incident`**: Mirrors `script/run_incident.R` (e.g. `obs_end` uses ischemic/haemorrhagic/death where applicable; MatchIt with `fu` in PS, no `dx.stroke_embo` in PS, `replace = FALSE`, `exact ~ trial_id`, no caliper; Table 1 before match omits `dx.stroke_embo`). Script 03 can use long-format SMD when **`regression_run_speedglm`** is **`TRUE`**.

`obs_end_cap_years` defaults to **3** in `analysis_config` to align with the manuscript’s three-year risk description. `study_start` / **`study_end`** bound the calendar study window and feed **`obs_end`** (via `pmin` with **`study_end`**) in 02. `script/run.R` uses **five** years in `obs_end`; set `analysis_config$obs_end_cap_years <- 5` to mirror that script’s follow-up cap.

## Regression (`script/03_regression_pooled_logistic.R`)

- **`script/02_analysis_primary.R`** stops after MatchIt and Table 1 (matched wide cohort). It does **not** run outcome regression.
- **`script/03_regression_pooled_logistic.R`** expands to **long person-time** (day or month steps per `index_by`) from `indx_date` to `obs_end`, merges **time-varying antiplatelet**, defines **`out`**, then:
  - **`speedglm`** (optional): pooled logistic coefficients with `expo` + `time` + `timesqr` and SMD-based adjustment — when **`regression_run_speedglm = TRUE`**.
  - **Cluster bootstrap RD/RR** (optional): same pattern as **`script/run_incident.R`** (~700–772) — resample **`id`**, fit **`speedglm`** with RR and RD formulas, g-computation over follow-up — when **`pipeline_style = "run_incident"`** and **`regression_run_bootstrap = TRUE`**. Replicates **`bootstrap_r`** (default **500**) and **`bootstrap_seed`** (default **456**).

At least one of **`regression_run_speedglm`** or **`regression_run_bootstrap`** must be **`TRUE`**. Defaults are both **`TRUE`**; set **`regression_run_bootstrap = FALSE`** to skip the long bootstrap (e.g. only need **`speedglm`** tables).

**Outcome `out`** (`outcome_definition` in `script/analysis_config.R`):

- **`pmin_composite`** (default): aligned with **`script/run_incident.R`** when `out_indicator == 1`. Let `t = pmin(date.ische, date.haem, date.death)`. Then `out = 1` if `t` is not missing, `t >= obs_date`, and `t` falls in the **same calendar month** as `obs_date`. This is **not** the same as counting each event type separately in that month.
- **`any_event_same_month`**: `out = 1` if **any** of ischemic stroke, haemorrhagic stroke, or death occurs in the same month as `obs_date` (older 03 behaviour).

**`script/run.R`** (unchanged) is a separate legacy pipeline. The **bootstrap + g-computation** path for primary RD/RR is now reproduced in **03** (via **`bootstrap_pooled_run_incident.R`**) when **`pipeline_style = "run_incident"`**.

## Optional checks

- **`script/compare_run_incident_vs_pipeline.R`**: compare artifacts between **`run_incident.R`** and this pipeline (when both have been run with matching settings).

## Dependencies

R packages: **`data.table`**, **`dplyr`**, **`lubridate`**, **`pbapply`**, **`readxl`**, **`tibble`** (01); **`tableone`**, **`labelled`**, **`MatchIt`** (02); **`tableone`**, **`broom`**, **`boot`**, **`speedglm`**, **`splitstackshape`** (03). The project uses **`renv`**; run **`renv::restore()`** (or install missing packages) before running scripts.

## Git

The repository can use **`renv`** for locked dependencies. Initialize or extend git as needed for your workflow.
