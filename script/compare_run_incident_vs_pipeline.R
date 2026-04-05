# compare_run_incident_vs_pipeline.R — pairwise checks: script/run_incident.R vs 01/02/03
#
# Run from project root:
#   Rscript script/compare_run_incident_vs_pipeline.R
#
# Requires the same data/ files used by run_incident. Place run_incident outputs next
# to pipeline outputs (or set RI_* paths below) to see diffs.
#
# Typical results when both pipelines use the same seqcohort:
#   - cleaned_20250508.RDS vs seqcohort_*.RDS: row counts and expo counts MATCH.
#   - Extra columns on RI seqcohort: MI/fracture outcomes if search_outcomes was run.
#   - tableone.xlsx vs CSV: cell text may differ if row labels differ (varLabels) or
#     one file is from an older 02 run; compare the 'n' row and after-matching counts.
#   - Pooled OR: RI workbook = bootstrap; PL csv = speedglm (not comparable numbers).

options(scipen = 6L, digits = 6L)
suppressPackageStartupMessages({
  library(data.table)
})

root <- if (dir.exists("data")) "." else stop("Working directory = project root.")

source(file.path(root, "script/analysis_config.R"))
cfg <- load_analysis_config(root = root)

# --- Paths: run_incident (RI) vs pipeline (PL) --------------------------------
ri_paths <- list(
  seqcohort_rds = file.path(root, "data/cleaned_20250508.RDS"),
  tableone_xlsx = file.path(root, "out/tableone.xlsx"),
  stroke_results_xlsx = file.path(root, "out/stroke_and_oac_results_20250508.xlsx"),
  eligible_pdf = file.path(root, "out/number of eligible.pdf"),
  eligible_by_expo_pdf = file.path(root, "out/number of eligible by users or not.pdf")
)

pl_paths <- list(
  seqcohort_rds = file.path(root, cfg$path_seqcohort_rds),
  tableone_before_csv = file.path(root, cfg$path_tableone_before_csv),
  tableone_after_csv = file.path(root, cfg$path_tableone_after_csv),
  tableone_long_smd_csv = file.path(root, cfg$path_tableone_long_smd_csv),
  person_trial_rds = file.path(root, cfg$path_person_trial_rds),
  est_pooled_expo_csv = file.path(root, cfg$path_est_pooled_expo_csv),
  est_pooled_all_csv = file.path(root, "out/est_pooled_logistic_all_terms.csv")
)

cat("\n=== Pairing (run_incident -> modular pipeline) ===\n\n")
cat(
  "1. Seqcohort RDS:\n   RI: data/cleaned_20250508.RDS\n",
  "  PL: ", cfg$path_seqcohort_rds, "\n\n",
  sep = ""
)
cat(
  "2. Table 1 (pre/post match):\n   RI: out/tableone.xlsx (merged before/after)\n",
  "  PL: ", cfg$path_tableone_before_csv, " + ", cfg$path_tableone_after_csv, "\n\n",
  sep = ""
)
cat(
  "3. Matched cohort for regression:\n   RI: in-memory only (no RDS); built from RI seqcohort\n",
  "  PL: ", cfg$path_person_trial_rds, "\n\n",
  sep = ""
)
cat(
  "4. Pooled logistic / bootstrap:\n   RI: out/stroke_and_oac_results_20250508.xlsx (boot)\n",
  "  PL: ", cfg$path_est_pooled_expo_csv, " (speedglm point estimates)\n\n",
  sep = ""
)
cat(
  "5. Long-format SMD table (adjustment):\n   RI: inside run_incident.R console / not saved as CSV\n",
  "  PL: ", cfg$path_tableone_long_smd_csv, "\n\n",
  sep = ""
)

# --- Helpers ------------------------------------------------------------------
exists_msg <- function(path) {
  if (file.exists(path)) {
    paste0("OK (", round(file.size(path) / 1024, 1), " KiB)")
  } else {
    "MISSING"
  }
}

summarise_seqcohort <- function(path) {
  sc <- readRDS(path)
  pt <- rbindlist(sc, fill = TRUE)
  if (!"trial_id" %in% names(pt)) {
    pt[, trial_id := rep(seq_along(sc), times = vapply(sc, nrow, 1L))]
  }
  list(
    path = path,
    n_trials = length(sc),
    n_rows = nrow(pt),
    n_expo1 = sum(pt$expo == 1L, na.rm = TRUE),
    n_expo0 = sum(pt$expo == 0L, na.rm = TRUE),
    n_ref = uniqueN(pt$Reference_Key),
    cols = sort(names(pt))
  )
}

compare_atomic <- function(label, a, b) {
  if (identical(a, b)) {
    cat("  ", label, ": MATCH\n", sep = "")
  } else {
    cat("  ", label, ": DIFF\n", sep = "")
    cat("    A: ", paste(format(a, trim = TRUE), collapse = " "), "\n", sep = "")
    cat("    B: ", paste(format(b, trim = TRUE), collapse = " "), "\n", sep = "")
  }
}

# --- 1. Seqcohort RDS --------------------------------------------------------
cat("=== 1. Seqcohort RDS ===\n")
cat("  RI: ", exists_msg(ri_paths$seqcohort_rds), "\n", sep = "")
cat("  PL: ", exists_msg(pl_paths$seqcohort_rds), "\n", sep = "")

if (file.exists(ri_paths$seqcohort_rds) && file.exists(pl_paths$seqcohort_rds)) {
  s_ri <- summarise_seqcohort(ri_paths$seqcohort_rds)
  s_pl <- summarise_seqcohort(pl_paths$seqcohort_rds)
  cat("\n  Summary run_incident (cleaned_20250508.RDS):\n")
  cat(
    "    rows=", s_ri$n_rows, " expo0=", s_ri$n_expo0, " expo1=", s_ri$n_expo1,
    " unique_ref=", s_ri$n_ref, " trials=", s_ri$n_trials, "\n",
    sep = ""
  )
  cat("\n  Summary pipeline (seqcohort from 01):\n")
  cat(
    "    rows=", s_pl$n_rows, " expo0=", s_pl$n_expo0, " expo1=", s_pl$n_expo1,
    " unique_ref=", s_pl$n_ref, " trials=", s_pl$n_trials, "\n",
    sep = ""
  )
  compare_atomic("n_rows stacked", s_ri$n_rows, s_pl$n_rows)
  compare_atomic("n_expo0", s_ri$n_expo0, s_pl$n_expo0)
  compare_atomic("n_expo1", s_ri$n_expo1, s_pl$n_expo1)
  only_ri <- setdiff(s_ri$cols, s_pl$cols)
  only_pl <- setdiff(s_pl$cols, s_ri$cols)
  if (length(only_ri) || length(only_pl)) {
    cat("  Column names: DIFF\n")
    if (length(only_ri)) {
      cat("    only in RI: ", paste(only_ri, collapse = ", "), "\n", sep = "")
    }
    if (length(only_pl)) {
      cat("    only in PL: ", paste(only_pl, collapse = ", "), "\n", sep = "")
    }
  } else {
    cat("  Column names: same set\n")
  }
} else {
  cat("  Skip detailed seqcohort compare (need both RDS files).\n")
}

# --- 2. Table 1 before / after (RI xlsx vs PL csv) -----------------------------
# run_incident strips " = 1" / " = 0" anywhere in rn (see run_incident.R ~469–470)
normalize_rn <- function(x) {
  x <- gsub(" = 1", "", x, fixed = TRUE)
  x <- gsub(" = 0", "", x, fixed = TRUE)
  trimws(x)
}

cell_diff <- function(a, b) {
  sa <- gsub("\\s+", " ", trimws(as.character(a)))
  sb <- gsub("\\s+", " ", trimws(as.character(b)))
  na_a <- is.na(a) | sa == ""
  na_b <- is.na(b) | sb == ""
  (na_a != na_b) | (!na_a & !na_b & sa != sb)
}

cat("\n=== 2. Table 1 before matching ===\n")
cat("  RI xlsx: ", exists_msg(ri_paths$tableone_xlsx), "\n", sep = "")
cat("  PL csv:  ", exists_msg(pl_paths$tableone_before_csv), "\n", sep = "")

if (file.exists(ri_paths$tableone_xlsx) && file.exists(pl_paths$tableone_before_csv)) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    cat("  Install readxl to compare tableone.xlsx\n")
  } else {
    sheets <- readxl::excel_sheets(ri_paths$tableone_xlsx)
    ri_merged <- as.data.table(
      readxl::read_xlsx(ri_paths$tableone_xlsx, sheet = sheets[[1L]])
    )
    pl_before <- fread(pl_paths$tableone_before_csv, header = TRUE, encoding = "UTF-8")
    if (!"rn" %in% names(pl_before) && all(grepl("^V[0-9]+$", names(pl_before)))) {
      setnames(pl_before, c("rn", "0", "1", "SMD")[seq_len(ncol(pl_before))])
    }
    cat(
      "  RI merged xlsx: ", nrow(ri_merged), " rows x ", ncol(ri_merged),
      " cols (sheet ", sheets[[1L]], ")\n",
      sep = ""
    )
    cat("  PL before csv: ", nrow(pl_before), " rows x ", ncol(pl_before), " cols\n", sep = "")

    # RI: merge(tableone_b4, tableone_after, by=rn) ->
    #   rn, before-non-user, before-user, before-SMD, after-...
    b0 <- "before matching-non-user"
    b1 <- "before matching-user"
    bs <- "before matching-SMD"
    pl_nm <- names(pl_before)
    has_pl012 <- length(pl_nm) >= 4L
    if (all(c(b0, b1, bs) %in% names(ri_merged)) && has_pl012) {
      ri_b <- ri_merged[, .(rn, ri0 = get(b0), ri1 = get(b1), riS = get(bs))]
      cn_b <- names(pl_before)
      pl_b <- pl_before[, .(
        rn = get(cn_b[1L]),
        pl0 = get(cn_b[2L]),
        pl1 = get(cn_b[3L]),
        plS = get(cn_b[4L])
      )]
      ri_b[, rn2 := normalize_rn(as.character(rn))]
      pl_b[, rn2 := normalize_rn(as.character(rn))]
      # Prefer merge on identical normalized labels; then report unmatched keys
      cmp <- merge(
        ri_b[, .(rn2, ri0, ri1, riS)],
        pl_b[, .(rn2, pl0, pl1, plS)],
        by = "rn2",
        all = TRUE
      )
      cat(
        "  Unmatched keys (RI only): ", sum(!is.na(cmp$ri0) & is.na(cmp$pl0)),
        "; (PL only): ", sum(is.na(cmp$ri0) & !is.na(cmp$pl0)), "\n",
        sep = ""
      )
      diff0 <- cmp[cell_diff(ri0, pl0)]
      diff1 <- cmp[cell_diff(ri1, pl1)]
      cat("  Before match: merged rows by normalized rn = ", nrow(cmp), "\n", sep = "")
      ncmp <- cmp[rn2 == "n"]
      if (nrow(ncmp) == 1L) {
        cat(
          "  'n' row (before match) RI expo0 / PL expo0: ",
          trimws(as.character(ncmp$ri0[1L])), " vs ",
          trimws(as.character(ncmp$pl0[1L])), "\n",
          sep = ""
        )
        cat(
          "  'n' row (before match) RI expo1 / PL expo1: ",
          trimws(as.character(ncmp$ri1[1L])), " vs ",
          trimws(as.character(ncmp$pl1[1L])), "\n",
          sep = ""
        )
      }
      cat("  Cells differing (expo=0 stratum): ", nrow(diff0), " rows\n", sep = "")
      if (nrow(diff0) > 0L && nrow(diff0) <= 30L) {
        print(diff0[, .(rn2, ri0, pl0)])
      }
      cat("  Cells differing (expo=1 stratum): ", nrow(diff1), " rows\n", sep = "")
      if (nrow(diff1) > 0L && nrow(diff1) <= 30L) {
        print(diff1[, .(rn2, ri1, pl1)])
      }
    } else {
      cat("  Could not map RI before-match columns; RI names: ", paste(names(ri_merged), collapse = ", "), "\n", sep = "")
    }

    pl_after <- pl_paths$tableone_after_csv
    if (file.exists(pl_after) && all(c("after matching-non-user", "after matching-user", "after matching-SMD") %in% names(ri_merged))) {
      pa <- fread(pl_after, header = TRUE, encoding = "UTF-8")
      if (!"rn" %in% names(pa) && all(grepl("^V[0-9]+$", names(pa)))) {
        setnames(pa, c("rn", "0", "1", "SMD")[seq_len(ncol(pa))])
      }
      a0 <- "after matching-non-user"
      a1 <- "after matching-user"
      asmd <- "after matching-SMD"
      ri_a <- ri_merged[, .(rn, ri0 = get(a0), ri1 = get(a1), riS = get(asmd))]
      cn_a <- names(pa)
      pl_a <- pa[, .(
        rn = get(cn_a[1L]),
        pl0 = get(cn_a[2L]),
        pl1 = get(cn_a[3L]),
        plS = get(cn_a[4L])
      )]
      ri_a[, rn2 := normalize_rn(as.character(rn))]
      pl_a[, rn2 := normalize_rn(as.character(rn))]
      cmp2 <- merge(
        ri_a[, .(rn2, ri0, ri1, riS)],
        pl_a[, .(rn2, pl0, pl1, plS)],
        by = "rn2",
        all = TRUE
      )
      diff_a0 <- cmp2[cell_diff(ri0, pl0)]
      diff_a1 <- cmp2[cell_diff(ri1, pl1)]
      diff_as <- cmp2[cell_diff(riS, plS)]
      cat("\n=== 2b. Table 1 after matching ===\n")
      ncmp2 <- cmp2[rn2 == "n"]
      if (nrow(ncmp2) == 1L) {
        cat(
          "  'n' row (after match) RI expo0 / PL expo0: ",
          trimws(as.character(ncmp2$ri0[1L])), " vs ",
          trimws(as.character(ncmp2$pl0[1L])), "\n",
          sep = ""
        )
        cat(
          "  'n' row (after match) RI expo1 / PL expo1: ",
          trimws(as.character(ncmp2$ri1[1L])), " vs ",
          trimws(as.character(ncmp2$pl1[1L])), "\n",
          sep = ""
        )
      }
      cat("  Cells differing (expo=0): ", nrow(diff_a0), "\n", sep = "")
      cat("  Cells differing (expo=1): ", nrow(diff_a1), "\n", sep = "")
      cat("  Cells differing (SMD): ", nrow(diff_as), "\n", sep = "")
      if (nrow(diff_a0) > 0L && nrow(diff_a0) <= 15L) {
        print(diff_a0[, .(rn2, ri0, pl0)])
      }
    }
  }
} else {
  cat("  Skip (need RI out/tableone.xlsx and PL tableone_before csv).\n")
}

# --- 3. person_trial (PL only) -----------------------------------------------
cat("\n=== 3. Matched person_trial (pipeline 02 only) ===\n")
cat("  ", exists_msg(pl_paths$person_trial_rds), " ", pl_paths$person_trial_rds, "\n", sep = "")
if (file.exists(pl_paths$person_trial_rds)) {
  pt <- readRDS(pl_paths$person_trial_rds)
  setDT(pt)
  cat(
    "  rows=", nrow(pt), " | expo: ",
    paste(pt[, .N, expo][, paste0(expo, "=", N)], collapse = ", "),
    "\n",
    sep = ""
  )
}

# --- 4. Pooled estimates CSV (PL) --------------------------------------------
cat("\n=== 4. Pooled expo estimates (pipeline 03) ===\n")
cat("  ", exists_msg(pl_paths$est_pooled_expo_csv), "\n", sep = "")
if (file.exists(pl_paths$est_pooled_expo_csv)) {
  est <- fread(pl_paths$est_pooled_expo_csv)
  print(est)
}
cat(
  "\n  Note: run_incident primary bootstrap CIs live in stroke_and_oac_results",
  "_20250508.xlsx — not the same object as speedglm ORs above.\n"
)

# --- 5. Other RI files (presence only) -----------------------------------------
cat("\n=== 5. Other run_incident outputs (presence) ===\n")
for (nm in names(ri_paths)) {
  if (nm %in% c("seqcohort_rds", "tableone_xlsx")) {
    next
  }
  p <- ri_paths[[nm]]
  cat("  ", nm, ": ", exists_msg(p), "\n", sep = "")
}

cat("\nDone.\n")
