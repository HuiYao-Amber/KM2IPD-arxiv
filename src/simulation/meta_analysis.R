# Trial collections, missing-HR masks, and random-effects meta-analysis.
meta_require_pkgs <- function(pkgs, stop_if_missing = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L && stop_if_missing) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      "\nPlease install with install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(length(missing) == 0L)
}

meta_require_pkgs(c("data.table", "survival", "metafor"), stop_if_missing = TRUE)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(metafor)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

source_single_trial_core <- function(core_path = "src/simulation/single_trial.R") {
  if (!file.exists(core_path)) {
    stop(
      "Cannot find single-trial core at: ", core_path, "\n",
      "Run the simulation entry point from the repository root.",
      call. = FALSE
    )
  }
  source(core_path)
  invisible(TRUE)
}

meta_parse_methods <- function(x) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(c("QP", "IPDfromKM", "KMtoIPD", "iKM"))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

meta_parse_numeric <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(default)
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

meta_parse_character <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

meta_safe <- function(expr, default = NULL) {
  tryCatch(expr, error = function(e) default)
}

meta_fit_rma <- function(logHR_vec, se_vec, re_method = "REML") {
  ok <- is.finite(logHR_vec) & is.finite(se_vec) & se_vec > 0
  yi <- logHR_vec[ok]
  sei <- se_vec[ok]
  if (length(yi) < 2L) {
    return(list(ok = FALSE, logHR = NA_real_, se = NA_real_, tau2 = NA_real_, I2 = NA_real_, k = length(yi), error = "fewer_than_two_studies"))
  }
  tryCatch({
    fit <- metafor::rma(yi = yi, sei = sei, method = re_method)
    list(
      ok = TRUE,
      logHR = as.numeric(fit$b),
      se = as.numeric(fit$se),
      tau2 = as.numeric(fit$tau2),
      I2 = as.numeric(fit$I2),
      k = as.integer(fit$k),
      error = NA_character_
    )
  }, error = function(e) {
    list(ok = FALSE, logHR = NA_real_, se = NA_real_, tau2 = NA_real_, I2 = NA_real_, k = length(yi), error = conditionMessage(e))
  })
}

meta_pool_true_ipd <- function(ipd_list) {
  if (length(ipd_list) < 1L) return(list(ok = FALSE, logHR = NA_real_, se = NA_real_, error = "empty_ipd_list"))
  ipd_all <- data.table::rbindlist(lapply(seq_along(ipd_list), function(i) {
    dt <- data.table::as.data.table(ipd_list[[i]])
    dt[, trial_id := factor(i)]
    dt
  }), fill = TRUE)
  if (!all(c("time", "status", "arm", "trial_id") %in% names(ipd_all))) {
    return(list(ok = FALSE, logHR = NA_real_, se = NA_real_, error = "bad_ipd_columns"))
  }
  ipd_all[, arm := factor(as.character(arm), levels = c("Control", "Treatment"))]
  tryCatch({
    fit <- survival::coxph(survival::Surv(time, status) ~ arm + strata(trial_id), data = ipd_all, ties = "efron")
    cf <- stats::coef(fit)
    vc <- stats::vcov(fit)
    beta <- unname(cf["armTreatment"])
    se <- sqrt(unname(vc["armTreatment", "armTreatment"]))
    list(ok = TRUE, logHR = beta, se = se, error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, logHR = NA_real_, se = NA_real_, error = conditionMessage(e))
  })
}

make_meta_grid <- function(
    preset = c("pilot", "main", "stress", "custom"),
    K_vec = NULL,
    n_per_arm_vec = NULL,
    p_miss_vec = NULL,
    tau_logHR_vec = NULL,
    tau_label_vec = NULL,
    censoring_vec = NULL,
    target_censoring_vec = NULL,
    digitize_sd_S_vec = NULL,
    weib_shape_vec = NULL,
    mean_time = 10
) {
  preset <- match.arg(preset)

  if (preset == "pilot") {
    K_vec <- K_vec %||% c(10)
    n_per_arm_vec <- n_per_arm_vec %||% c(100)
    p_miss_vec <- p_miss_vec %||% c(0, 0.4, 0.6, 0.8)
    tau_logHR_vec <- tau_logHR_vec %||% c(0.05, 0.25)
    tau_label_vec <- tau_label_vec %||% c("low", "high")
    censoring_vec <- censoring_vec %||% c("random", "front")
    target_censoring_vec <- target_censoring_vec %||% c(0.30)
    digitize_sd_S_vec <- digitize_sd_S_vec %||% c(0.01)
    weib_shape_vec <- weib_shape_vec %||% c(1.0)
  } else if (preset == "main") {
    K_vec <- K_vec %||% c(30)
    n_per_arm_vec <- n_per_arm_vec %||% c(250)
    p_miss_vec <- p_miss_vec %||% c(0, 0.2, 0.4, 0.6, 0.8)
    tau_logHR_vec <- tau_logHR_vec %||% c(0.05, 0.10, 0.25)
    tau_label_vec <- tau_label_vec %||% c("low", "middle", "high")
    censoring_vec <- censoring_vec %||% c("random", "front", "back")
    target_censoring_vec <- target_censoring_vec %||% c(0.30)
    digitize_sd_S_vec <- digitize_sd_S_vec %||% c(0.01)
    weib_shape_vec <- weib_shape_vec %||% c(1.0)
  } else if (preset == "stress") {
    K_vec <- K_vec %||% c(10, 15)
    n_per_arm_vec <- n_per_arm_vec %||% c(100)
    p_miss_vec <- p_miss_vec %||% c(0.6, 0.8)
    tau_logHR_vec <- tau_logHR_vec %||% c(0.25)
    tau_label_vec <- tau_label_vec %||% c("high")
    censoring_vec <- censoring_vec %||% c("random", "front", "back", "informative")
    target_censoring_vec <- target_censoring_vec %||% c(0.60)
    digitize_sd_S_vec <- digitize_sd_S_vec %||% c(0.04)
    weib_shape_vec <- weib_shape_vec %||% c(1.0)
  } else {
    K_vec <- K_vec %||% c(10)
    n_per_arm_vec <- n_per_arm_vec %||% c(250)
    p_miss_vec <- p_miss_vec %||% c(0, 0.4, 0.8)
    tau_logHR_vec <- tau_logHR_vec %||% c(0.05, 0.25)
    tau_label_vec <- tau_label_vec %||% paste0("tau", seq_along(tau_logHR_vec))
    censoring_vec <- censoring_vec %||% c("random")
    target_censoring_vec <- target_censoring_vec %||% c(0.30)
    digitize_sd_S_vec <- digitize_sd_S_vec %||% c(0.01)
    weib_shape_vec <- weib_shape_vec %||% c(1.0)
  }

  if (length(tau_label_vec) != length(tau_logHR_vec)) {
    tau_label_vec <- paste0("tau", seq_along(tau_logHR_vec))
  }
  tau_dt <- data.table(tau_label = tau_label_vec, tau_logHR = tau_logHR_vec)

  base_grid <- CJ(
    K = as.integer(K_vec),
    n_per_arm = as.integer(n_per_arm_vec),
    censoring = censoring_vec,
    target_censoring = as.numeric(target_censoring_vec),
    digitize_sd_S = as.numeric(digitize_sd_S_vec),
    weib_shape = as.numeric(weib_shape_vec),
    tau_idx = seq_len(nrow(tau_dt))
  )
  base_grid[, `:=`(
    tau_label = tau_dt$tau_label[tau_idx],
    tau_logHR = tau_dt$tau_logHR[tau_idx],
    mean_time = mean_time
  )]
  base_grid[, tau_idx := NULL]

  grid <- base_grid[, .(p_miss = as.numeric(p_miss_vec)), by = .(
    K, n_per_arm, censoring, target_censoring, digitize_sd_S,
    weib_shape, tau_label, tau_logHR, mean_time
  )]

  grid[, scenario_id := sprintf(
    "K%d_n%d_miss%02d_tau%s_cens%s_tc%02d_dig%03d_shape%s",
    K, n_per_arm, as.integer(round(100 * p_miss)), tau_label, censoring,
    as.integer(round(100 * target_censoring)), as.integer(round(1000 * digitize_sd_S)),
    gsub("\\.", "p", as.character(weib_shape))
  )]

  setcolorder(grid, c("scenario_id", "K", "n_per_arm", "p_miss", "tau_label", "tau_logHR",
                      "censoring", "target_censoring", "digitize_sd_S", "weib_shape", "mean_time"))
  grid[]
}

make_trial_specific_scenario_lib <- function(hr, shape = 1.0, mean_time = 10, id = "meta_trial_dynamic") {
  scn <- make_weibull_ph_scenario(id = id, shape = shape, hr = hr, mean_time = mean_time)
  x <- list(scn)
  names(x) <- id
  x
}

simulate_meta_base_world <- function(
    K,
    mu_HR = 0.75,
    tau_logHR,
    n_per_arm,
    censoring,
    target_censoring,
    digitize_sd_S,
    weib_shape = 1.0,
    mean_time = 10,
    methods = c("QP", "IPDfromKM", "KMtoIPD", "iKM"),
    digitize_mode = "event_times",
    n_grid = 200,
    re_method = "REML",
    rep_id = 1L,
    seed = NULL,
    method_timeout_sec = 180
) {
  if (!is.null(seed)) set.seed(seed)
  methods <- unique(methods)
  mu_logHR <- log(mu_HR)

  logHR_gen <- stats::rnorm(K, mean = mu_logHR, sd = tau_logHR)
  HR_gen <- exp(logHR_gen)

  trial_rows <- vector("list", K)
  trial_method_rows <- list()
  ipd_true_list <- vector("list", K)
  idx_tm <- 1L

  for (k in seq_len(K)) {
    trial_seed <- if (!is.null(seed)) as.integer(seed + 1000L * k) else NULL
    scenario_id <- sprintf("meta_trial_%03d", k)
    scn_lib <- make_trial_specific_scenario_lib(
      hr = HR_gen[k], shape = weib_shape, mean_time = mean_time, id = scenario_id
    )

    sim <- simulate_single_trial(
      scenario_id = scenario_id,
      scenario_lib = scn_lib,
      n_control = n_per_arm,
      n_treatment = n_per_arm,
      censoring = censoring,
      target_censoring = target_censoring,
      admin_time = NULL,
      fixed_x_end = NULL,
      seed = trial_seed
    )
    ipd_true_list[[k]] <- sim$ipd

    est_true <- .estimate_loghr_cox(sim$ipd)
    tau_eval <- max(sim$risk_table$time, na.rm = TRUE)

    digi <- digitize_km_numeric(
      ipd_true = sim$ipd,
      risk_table = sim$risk_table,
      mode = digitize_mode,
      n_grid = n_grid,
      sd_S = digitize_sd_S,
      x_end_plot = sim$axis$x_end_plot,
      seed = if (!is.null(trial_seed)) trial_seed + 10L else NULL
    )

    trial_rows[[k]] <- data.table(
      rep = rep_id,
      trial = k,
      HR_gen = HR_gen[k],
      logHR_gen = logHR_gen[k],
      logHR_true_hat = est_true$logHR,
      se_true_hat = est_true$se,
      true_ok = isTRUE(est_true$ok),
      achieved_censor_prop = sim$censor_prop,
      true_events = sim$ipd[, sum(status == 1L)],
      true_censors = sim$ipd[, sum(status == 0L)],
      tau_eval = tau_eval
    )

    for (m in methods) {
      t0 <- proc.time()[[3]]
      rec <- tryCatch({
        .with_optional_timeout({
          reconstruct_ipd_twoarm_final(
            digi = digi,
            risk_table = sim$risk_table,
            method = m,
            curve_map = c(Control = 1, Treatment = 2),
            total_events = sim$total_events,
            censor_marks_by_arm = sim$censor_marks_by_arm,
            kmtoipd_use_marks = TRUE,
            drop_tol = 1e-8
          )
        }, timeout_sec = method_timeout_sec)
      }, error = function(e) {
        list(method = m, ipd = NULL, error = conditionMessage(e))
      })
      runtime_sec <- proc.time()[[3]] - t0

      if (is.null(rec$ipd)) {
        tm <- data.table(
          rep = rep_id, trial = k, recon_method = m,
          recon_ok = FALSE,
          recon_error = rec$error %||% "unknown reconstruction error",
          logHR_recon_hat = NA_real_, se_recon_hat = NA_real_,
          runtime_sec = as.numeric(runtime_sec),
          km_rmse_mean = NA_real_, abs_cox_logHR_diff = NA_real_, abs_rmst_diff = NA_real_,
          nrisk_max_abs = NA_real_, abs_event_diff_total = NA_real_, event_diff_total = NA_real_
        )
      } else {
        est_rec <- .estimate_loghr_cox(rec$ipd)
        met <- tryCatch({
          evaluate_reconstruction(
            true_ipd = sim$ipd,
            rec_ipd = rec$ipd,
            digi_dt = digi$data,
            risk_table = sim$risk_table,
            tau_eval = tau_eval
          )
        }, error = function(e) {
          data.table(eval_ok = FALSE, km_rmse_mean = NA_real_, abs_cox_logHR_diff = NA_real_,
                     abs_rmst_diff = NA_real_, nrisk_max_abs = NA_real_,
                     event_diff_total = NA_real_)
        })
        tm <- data.table(
          rep = rep_id, trial = k, recon_method = m,
          recon_ok = isTRUE(est_rec$ok),
          recon_error = ifelse(isTRUE(est_rec$ok), NA_character_, est_rec$error %||% "cox_failed"),
          logHR_recon_hat = est_rec$logHR,
          se_recon_hat = est_rec$se,
          runtime_sec = as.numeric(runtime_sec),
          km_rmse_mean = met$km_rmse_mean[1],
          abs_cox_logHR_diff = abs(est_rec$logHR - est_true$logHR),
          abs_rmst_diff = met$abs_rmst_diff[1],
          nrisk_max_abs = met$nrisk_max_abs[1],
          abs_event_diff_total = abs(met$event_diff_total[1]),
          event_diff_total = met$event_diff_total[1]
        )
      }
      trial_method_rows[[idx_tm]] <- tm
      idx_tm <- idx_tm + 1L
    }
  }

  trials <- rbindlist(trial_rows, fill = TRUE)
  trial_methods <- rbindlist(trial_method_rows, fill = TRUE)
  pooled <- meta_pool_true_ipd(ipd_true_list)

  list(
    trials = trials,
    trial_methods = trial_methods,
    pooled_ipd = pooled,
    ipd_true_list = ipd_true_list
  )
}

fit_meta_strategies_for_missingness <- function(
    trials,
    trial_methods,
    pooled_ipd,
    p_miss,
    u_missing,
    methods,
    mu_HR = 0.75,
    re_method = "REML"
) {
  trials <- copy(as.data.table(trials))
  trial_methods <- copy(as.data.table(trial_methods))
  trials[, is_missing := u_missing < p_miss]
  mu_logHR <- log(mu_HR)

  out <- list()
  idx <- 1L

  ad <- trials[is_missing == FALSE & true_ok == TRUE & is.finite(logHR_true_hat) & is.finite(se_true_hat) & se_true_hat > 0]
  fit_ad <- meta_fit_rma(ad$logHR_true_hat, ad$se_true_hat, re_method = re_method)
  out[[idx]] <- data.table(
    meta_strategy = "AD_only",
    recon_method = NA_character_,
    ok = fit_ad$ok,
    logHR_hat = fit_ad$logHR,
    se_hat = fit_ad$se,
    tau2_hat = fit_ad$tau2,
    I2_hat = fit_ad$I2,
    k_used = fit_ad$k,
    n_missing = sum(trials$is_missing),
    fit_error = fit_ad$error
  ); idx <- idx + 1L

  full <- trials[true_ok == TRUE & is.finite(logHR_true_hat) & is.finite(se_true_hat) & se_true_hat > 0]
  fit_full <- meta_fit_rma(full$logHR_true_hat, full$se_true_hat, re_method = re_method)
  out[[idx]] <- data.table(
    meta_strategy = "Full_AD_oracle",
    recon_method = NA_character_,
    ok = fit_full$ok,
    logHR_hat = fit_full$logHR,
    se_hat = fit_full$se,
    tau2_hat = fit_full$tau2,
    I2_hat = fit_full$I2,
    k_used = fit_full$k,
    n_missing = sum(trials$is_missing),
    fit_error = fit_full$error
  ); idx <- idx + 1L

  out[[idx]] <- data.table(
    meta_strategy = "Pooled_IPD_oracle",
    recon_method = NA_character_,
    ok = pooled_ipd$ok,
    logHR_hat = pooled_ipd$logHR,
    se_hat = pooled_ipd$se,
    tau2_hat = NA_real_,
    I2_hat = NA_real_,
    k_used = nrow(trials),
    n_missing = sum(trials$is_missing),
    fit_error = pooled_ipd$error
  ); idx <- idx + 1L

  for (m in methods) {
    tm <- trial_methods[recon_method == m]
    tm <- merge(tm, trials[, .(rep, trial, is_missing, logHR_true_hat, se_true_hat, true_ok)], by = c("rep", "trial"), all.x = TRUE)

    hy_true <- tm[is_missing == FALSE & true_ok == TRUE & is.finite(logHR_true_hat) & is.finite(se_true_hat) & se_true_hat > 0,
                  .(logHR = logHR_true_hat, se = se_true_hat)]
    hy_rec <- tm[is_missing == TRUE & recon_ok == TRUE & is.finite(logHR_recon_hat) & is.finite(se_recon_hat) & se_recon_hat > 0,
                 .(logHR = logHR_recon_hat, se = se_recon_hat)]
    hy <- rbindlist(list(hy_true, hy_rec), use.names = TRUE, fill = TRUE)
    fit_hy <- meta_fit_rma(hy$logHR, hy$se, re_method = re_method)
    out[[idx]] <- data.table(
      meta_strategy = paste0("Hybrid_", m),
      recon_method = m,
      ok = fit_hy$ok,
      logHR_hat = fit_hy$logHR,
      se_hat = fit_hy$se,
      tau2_hat = fit_hy$tau2,
      I2_hat = fit_hy$I2,
      k_used = fit_hy$k,
      n_missing = sum(trials$is_missing),
      fit_error = fit_hy$error
    ); idx <- idx + 1L

    rm <- tm[recon_ok == TRUE & is.finite(logHR_recon_hat) & is.finite(se_recon_hat) & se_recon_hat > 0]
    fit_rm <- meta_fit_rma(rm$logHR_recon_hat, rm$se_recon_hat, re_method = re_method)
    out[[idx]] <- data.table(
      meta_strategy = paste0("Recon_only_", m),
      recon_method = m,
      ok = fit_rm$ok,
      logHR_hat = fit_rm$logHR,
      se_hat = fit_rm$se,
      tau2_hat = fit_rm$tau2,
      I2_hat = fit_rm$I2,
      k_used = fit_rm$k,
      n_missing = sum(trials$is_missing),
      fit_error = fit_rm$error
    ); idx <- idx + 1L
  }

  res <- rbindlist(out, fill = TRUE)
  res[, `:=`(
    p_miss = p_miss,
    err = logHR_hat - mu_logHR,
    HR_hat = exp(logHR_hat)
  )]
  res[]
}

run_meta_simulation <- function(
    grid,
    n_rep = 20,
    methods = c("QP", "IPDfromKM", "KMtoIPD", "iKM"),
    mu_HR = 0.75,
    seed_base = 20260619L,
    digitize_mode = "event_times",
    n_grid = 200,
    re_method = "REML",
    method_timeout_sec = 180,
    verbose = TRUE
) {
  grid <- as.data.table(grid)
  methods <- unique(methods)

  base_worlds <- unique(grid[, .(
    K, n_per_arm, tau_label, tau_logHR, censoring,
    target_censoring, digitize_sd_S, weib_shape, mean_time
  )])
  base_worlds[, base_id := .I]

  methods_out <- list()
  trials_out <- list()
  trial_methods_out <- list()
  idx <- 1L
  idx_trials <- 1L
  idx_tm <- 1L

  for (b in seq_len(nrow(base_worlds))) {
    bw <- base_worlds[b]
    p_values <- sort(unique(grid[
      K == bw$K & n_per_arm == bw$n_per_arm & tau_label == bw$tau_label &
        censoring == bw$censoring & target_censoring == bw$target_censoring &
        digitize_sd_S == bw$digitize_sd_S & weib_shape == bw$weib_shape,
      p_miss
    ]))

    for (r in seq_len(n_rep)) {
      seed_world <- as.integer(seed_base + b * 100000L + r)
      if (verbose && (r == 1L || r %% 10L == 0L)) {
        message(sprintf("Base %d/%d | rep %d/%d | K=%d n=%d tau=%s censor=%s targetC=%.2f digit=%.3f",
                        b, nrow(base_worlds), r, n_rep, bw$K, bw$n_per_arm, bw$tau_label,
                        bw$censoring, bw$target_censoring, bw$digitize_sd_S))
      }

      one <- simulate_meta_base_world(
        K = bw$K,
        mu_HR = mu_HR,
        tau_logHR = bw$tau_logHR,
        n_per_arm = bw$n_per_arm,
        censoring = bw$censoring,
        target_censoring = bw$target_censoring,
        digitize_sd_S = bw$digitize_sd_S,
        weib_shape = bw$weib_shape,
        mean_time = bw$mean_time,
        methods = methods,
        digitize_mode = digitize_mode,
        n_grid = n_grid,
        re_method = re_method,
        rep_id = r,
        seed = seed_world,
        method_timeout_sec = method_timeout_sec
      )

      trials <- one$trials
      trials[, `:=`(
        base_id = bw$base_id, seed = seed_world,
        K = bw$K, n_per_arm = bw$n_per_arm,
        tau_label = bw$tau_label, tau_logHR = bw$tau_logHR,
        censoring = bw$censoring, target_censoring = bw$target_censoring,
        digitize_sd_S = bw$digitize_sd_S, weib_shape = bw$weib_shape,
        mean_time = bw$mean_time
      )]
      set.seed(seed_world + 7777L)
      trials[, u_missing := runif(.N)]

      tm <- one$trial_methods
      tm[, `:=`(
        base_id = bw$base_id, seed = seed_world,
        K = bw$K, n_per_arm = bw$n_per_arm,
        tau_label = bw$tau_label, tau_logHR = bw$tau_logHR,
        censoring = bw$censoring, target_censoring = bw$target_censoring,
        digitize_sd_S = bw$digitize_sd_S, weib_shape = bw$weib_shape,
        mean_time = bw$mean_time
      )]

      trials_out[[idx_trials]] <- trials; idx_trials <- idx_trials + 1L
      trial_methods_out[[idx_tm]] <- tm; idx_tm <- idx_tm + 1L

      for (pm in p_values) {
        res <- fit_meta_strategies_for_missingness(
          trials = trials,
          trial_methods = tm,
          pooled_ipd = one$pooled_ipd,
          p_miss = pm,
          u_missing = trials$u_missing,
          methods = methods,
          mu_HR = mu_HR,
          re_method = re_method
        )
        res[, `:=`(
          base_id = bw$base_id, rep = r, seed = seed_world,
          K = bw$K, n_per_arm = bw$n_per_arm,
          tau_label = bw$tau_label, tau_logHR = bw$tau_logHR,
          censoring = bw$censoring, target_censoring = bw$target_censoring,
          digitize_sd_S = bw$digitize_sd_S, weib_shape = bw$weib_shape,
          mean_time = bw$mean_time,
          scenario_id = sprintf("K%d_n%d_miss%02d_tau%s_cens%s_tc%02d_dig%03d_shape%s",
                                bw$K, bw$n_per_arm, as.integer(round(100 * pm)), bw$tau_label,
                                bw$censoring, as.integer(round(100 * bw$target_censoring)),
                                as.integer(round(1000 * bw$digitize_sd_S)),
                                gsub("\\.", "p", as.character(bw$weib_shape)))
        )]
        methods_out[[idx]] <- res
        idx <- idx + 1L
      }
    }
  }

  list(
    meta_methods_long = rbindlist(methods_out, fill = TRUE),
    meta_trials_long = rbindlist(trials_out, fill = TRUE),
    meta_trial_methods_long = rbindlist(trial_methods_out, fill = TRUE),
    grid = grid,
    base_worlds = base_worlds
  )
}

summarise_meta_performance <- function(meta_methods_long, mu_HR = 0.75, alpha = 0.05) {
  dt <- copy(as.data.table(meta_methods_long))
  target <- log(mu_HR)
  z <- stats::qnorm(1 - alpha / 2)
  dt[, err := logHR_hat - target]
  dt[, ci_lower := logHR_hat - z * se_hat]
  dt[, ci_upper := logHR_hat + z * se_hat]
  dt[, cover := is.finite(ci_lower) & is.finite(ci_upper) & ci_lower <= target & ci_upper >= target]
  dt[, abs_err := abs(err)]

  perf <- dt[, .(
    n_rep = .N,
    success_rate = mean(ok & is.finite(logHR_hat) & is.finite(se_hat), na.rm = TRUE),
    mean_k_used = mean(k_used, na.rm = TRUE),
    median_k_used = stats::median(k_used, na.rm = TRUE),
    bias_logHR = mean(err, na.rm = TRUE),
    abs_bias_logHR = mean(abs_err, na.rm = TRUE),
    median_abs_error_logHR = stats::median(abs_err, na.rm = TRUE),
    p90_abs_error_logHR = as.numeric(stats::quantile(abs_err, 0.90, na.rm = TRUE)),
    rmse_logHR = sqrt(mean(err^2, na.rm = TRUE)),
    mean_se = mean(se_hat, na.rm = TRUE),
    coverage = mean(cover, na.rm = TRUE),
    mean_tau2_hat = mean(tau2_hat, na.rm = TRUE),
    mean_I2_hat = mean(I2_hat, na.rm = TRUE)
  ), by = .(
    scenario_id, K, n_per_arm, p_miss, tau_label, tau_logHR,
    censoring, target_censoring, digitize_sd_S, weib_shape,
    meta_strategy, recon_method
  )]

  data.table::setorder(perf, censoring, tau_logHR, p_miss, meta_strategy)
  perf[]
}

summarise_trial_reconstruction <- function(meta_trial_methods_long) {
  dt <- copy(as.data.table(meta_trial_methods_long))
  dt[, abs_trial_logHR_err := abs(abs_cox_logHR_diff)]
  summ <- dt[, .(
    n_trial_methods = .N,
    recon_success_rate = mean(recon_ok, na.rm = TRUE),
    median_trial_abs_logHR_error = stats::median(abs_cox_logHR_diff, na.rm = TRUE),
    p90_trial_abs_logHR_error = as.numeric(stats::quantile(abs_cox_logHR_diff, 0.90, na.rm = TRUE)),
    median_km_rmse = stats::median(km_rmse_mean, na.rm = TRUE),
    p90_km_rmse = as.numeric(stats::quantile(km_rmse_mean, 0.90, na.rm = TRUE)),
    median_abs_rmst_diff = stats::median(abs_rmst_diff, na.rm = TRUE),
    p90_abs_rmst_diff = as.numeric(stats::quantile(abs_rmst_diff, 0.90, na.rm = TRUE)),
    median_event_diff = stats::median(abs_event_diff_total, na.rm = TRUE),
    median_nrisk_max_abs = stats::median(nrisk_max_abs, na.rm = TRUE),
    median_runtime_sec = stats::median(runtime_sec, na.rm = TRUE)
  ), by = .(
    K, n_per_arm, tau_label, tau_logHR, censoring,
    target_censoring, digitize_sd_S, weib_shape, recon_method
  )]
  summ[]
}

compute_hybrid_gain <- function(perf) {
  dt <- copy(as.data.table(perf))
  ad <- dt[meta_strategy == "AD_only", .(
    scenario_id, rmse_AD_only = rmse_logHR,
    bias_AD_only = bias_logHR,
    coverage_AD_only = coverage,
    mean_k_AD_only = mean_k_used
  )]
  hy <- dt[grepl("^Hybrid_", meta_strategy)]
  out <- merge(hy, ad, by = "scenario_id", all.x = TRUE)
  out[, `:=`(
    rmse_gain_vs_AD = rmse_AD_only - rmse_logHR,
    abs_bias_gain_vs_AD = abs(bias_AD_only) - abs(bias_logHR),
    coverage_diff_vs_AD = coverage - coverage_AD_only,
    k_gain_vs_AD = mean_k_used - mean_k_AD_only
  )]
  out[]
}
