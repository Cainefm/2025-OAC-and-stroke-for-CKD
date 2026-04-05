# bootstrap_pooled_run_incident.R — pooled RD/RR bootstrap (script/run_incident.R ~51–158)
# Globals: person_trial, formula_rr, formula_rd, n_gcomp, K, current_bootstrap, pb
# (n_gcomp = horizon steps: 12 * years if monthly index, ~365.25 * years if daily.)

std.boot <- function(data, indices) {
  current_bootstrap <<- current_bootstrap + 1L
  setTxtProgressBar(pb, current_bootstrap)

  boot_ids <- data.table(id = data$id[indices])[, bid := .I]

  d <- merge(
    boot_ids,
    person_trial,
    by = "id",
    all.x = TRUE,
    allow.cartesian = TRUE
  )[, bid_new := .GRP, by = .(bid, trial_id)]

  d <- dplyr::left_join(
    boot_ids,
    person_trial,
    by = "id",
    relationship = "many-to-many"
  )
  d$bid_new <- interaction(d$bid, d$trial_id)

  risk_diff <- calculate_risk(d, formula_rd, "diff")
  risk_ratio <- calculate_risk(d, formula_rr, "ratio")
  if (is.null(risk_diff) || is.null(risk_ratio)) {
    return(c(NA, NA, NA, NA, NA))
  }
  c(
    risk_diff[time_0 == K, risk0],
    risk_diff[time_0 == K, risk1],
    risk_diff[time_0 == K, rd],
    risk_ratio[time_0 == K, rr],
    risk_diff[time_0 == K, rr]
  )
}

calculate_risk <- function(dt, formula, type) {
  model <- tryCatch(
    speedglm::speedglm(formula, family = binomial(link = "logit"), data = dt),
    error = function(e) {
      message(
        "Model failed at bootstrap ", current_bootstrap, ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  if (is.null(model)) {
    return(NULL)
  }
  if (isTRUE(getOption("bootstrap_pooled_verbose", FALSE))) {
    print(model)
  }

  base_data <- splitstackshape::expandRows(
    dt[time == 0, ],
    count = n_gcomp,
    count.is.col = FALSE
  )[, time := rep(0:(n_gcomp - 1L), dt[time == 0, .N])][, timesqr := time^2]

  control <- data.table::copy(base_data)[, expo := 0][,
    p_event := predict(model, .SD, type = "response")
  ]

  treatment <- data.table::copy(base_data)[, expo := 1][,
    p_event := predict(model, .SD, type = "response")
  ]

  process_survival <- function(dt) {
    dt[, surv := cumprod(1 - p_event), by = bid_new][, risk := 1 - surv]
  }

  control_risk <- process_survival(control)
  treatment_risk <- process_survival(treatment)

  aggregate_risk <- function(dt) {
    dt[, .(risk = mean(risk)), by = .(time, expo)]
  }

  graph_pred <- merge(
    aggregate_risk(control_risk)[, .(time, risk0 = risk)],
    aggregate_risk(treatment_risk)[, .(time, risk1 = risk)],
    by = "time"
  )[, `:=`(time_0 = time + 1L, rd = risk1 - risk0, rr = risk1 / risk0)]

  zero_row <- data.table(time = 0, risk0 = 0, risk1 = 0, rd = 0, rr = 0, time_0 = 0)
  rbindlist(list(zero_row, graph_pred), fill = TRUE)
}

extract_ci <- function(results, index) {
  tcol <- results$t[, index]
  if (!any(is.finite(tcol))) {
    return(c(NA_real_, NA_real_))
  }
  tryCatch(
    boot::boot.ci(results, conf = 0.95, type = "perc", index = index)$percent[4:5],
    error = function(e) c(NA_real_, NA_real_)
  )
}

orgnize_ci <- function(x, model_type, ids) {
  est <- dcast(
    data.table(
      arm = c(
        "risk_nonuser",
        "risk_user",
        "risk_diff",
        "risk_ratio",
        "risk_ratio_rdmodel"
      ),
      est = round(x$t0, 4),
      lower = sapply(1:5, function(i) round(extract_ci(x, i)[1], 4)),
      upper = sapply(1:5, function(i) round(extract_ci(x, i)[2], 4))
    )[, estimation := paste0(
      round(est, 2), " (", round(lower, 2), ", ", round(upper, 2), ")"
    )][, .(arm, estimation)],
    . ~ arm,
    value.var = "estimation"
  )[, .(risk_nonuser, risk_user, risk_diff, risk_ratio, risk_ratio_rdmodel)]

  num <- transpose(
    as.data.table(rbind(
      as.data.frame(
        person_trial[id %in% ids$id, .N, expo][,
          num.of.p := getOption("bootstrap_person_time_label", "person-periods")
        ]
      ),
      as.data.frame(
        person_trial[id %in% ids$id, ][
          ,
          unique(.SD),
          .SDcols = c("Reference_Key", "expo")
        ][, .N, expo][, num.of.p := "individuals"]
      )
    ))[, .(
      num.of.p = paste0(
        num.of.p, " ",
        ifelse(expo == 0, "non-users", "users")
      ),
      N
    )],
    make.names = "num.of.p"
  )
  cbind(model_type, num, est)
}
