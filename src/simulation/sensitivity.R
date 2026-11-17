# Illustrative designs for running the sensitivity workflows without external inputs.
example_sensitivity_design <- function(part) {
  phase_one <- build_phase1_design()
  baseline <- phase_one[scenario_id == "W_shape1_HR067" & censoring == "random"][1]
  if (part == "input-ablation") {
    baseline[, `:=`(design_id = 1L, regime = "Example")]
    return(baseline)
  }
  noise_levels <- unique(phase_one$digitize_sd_S)
  design <- baseline[rep(1L, length(noise_levels))]
  design[, `:=`(design_id = seq_along(noise_levels), digitize_sd_S = noise_levels,
                sensitivity_block = "digitization_noise", sensitivity_level = as.character(noise_levels),
                sensitivity_value_num = noise_levels, design_label = "Example digitization noise")]
  design[]
}

run_single_factor <- function(design, reps, out_dir) {
  rows <- list()
  options(km2ipd.reconstructkm_timeout = benchmark_config$sensitivity_ikm_timeout_sec)
  for (i in seq_len(nrow(design))) for (r in seq_len(reps)) {
    message(sprintf("Single-factor design %d/%d, replicate %d/%d", i, nrow(design), r, reps))
    for (m in benchmark_config$methods) {
      timeout <- switch(m, iKM = 15, QP = benchmark_config$sensitivity_qp_timeout_sec,
                        benchmark_config$timeout_sec)
      rows[[length(rows) + 1L]] <- run_one_single_trial_rep(
        design[i], r, methods = m,
        seed_base = benchmark_config$seeds$single_factor,
        method_timeout_sec = timeout
      )
    }
    data.table::fwrite(data.table::rbindlist(rows, fill = TRUE),
                      file.path(out_dir, "results.csv"))
  }
  invisible(data.table::rbindlist(rows, fill = TRUE))
}

input_variants <- function(qp = NULL) {
  if (is.null(qp)) {
    settings <- benchmark_config$method_options
    qp <- data.table::data.table(
      variant = "QP_baseline", configuration = settings$km2ipd.cifresolve_configuration,
      information = "Baseline risk table and censor marks",
      strict_tick = settings$km2ipd.cifresolve_strict_tick,
      strict_dec = settings$km2ipd.cifresolve_strict_dec,
      constr_tol = settings$km2ipd.cifresolve_constr_tol
    )
  }
  qp <- data.table::copy(qp)
  qp[, method := "QP"]
  other <- data.table::data.table(
    variant = c("Guyot_style_iKM", "IPDfromKM_no_total_events", "IPDfromKM_total_events",
                "KMtoIPD_exact_marks", "KMtoIPD_no_marks"),
    method = c("iKM", "IPDfromKM", "IPDfromKM", "KMtoIPD", "KMtoIPD"),
    information = c("Risk table only", "No total-event input", "True total events",
                    "Exact censor marks", "No censor marks")
  )
  data.table::rbindlist(list(other, qp), fill = TRUE)
}

run_auxiliary_inputs <- function(design, reps, out_dir, variants = input_variants()) {
  data.table::fwrite(variants, file.path(out_dir, "variants.csv"))
  rows <- list()
  options(km2ipd.reconstructkm_timeout = benchmark_config$sensitivity_ikm_timeout_sec)
  for (i in seq_len(nrow(design))) for (r in seq_len(reps)) {
    d <- design[i]
    seed <- as.integer(benchmark_config$seeds$input_ablation + 10000L * d$design_id + r)
    sim <- simulate_single_trial(
      scenario_id = d$scenario_id, n_control = d$n_control, n_treatment = d$n_treatment,
      censoring = d$censoring, target_censoring = d$target_censoring, seed = seed
    )
    digi <- digitize_km_numeric(sim$ipd, sim$risk_table, mode = "event_times",
                              sd_S = d$digitize_sd_S, x_end_plot = sim$axis$x_end_plot,
                              seed = seed + 100L)
    for (j in seq_len(nrow(variants))) {
      v <- variants[j]
      message(sprintf("Auxiliary-input design %d/%d, replicate %d/%d, variant %s",
                      i, nrow(design), r, reps, v$variant))
      options(benchmark_config$method_options)
      options(km2ipd.reconstructkm_timeout = benchmark_config$sensitivity_ikm_timeout_sec)
      timeout <- switch(v$method, iKM = 15, QP = benchmark_config$sensitivity_qp_timeout_sec,
                        benchmark_config$timeout_sec)
      start <- proc.time()[[3]]
      rec <- tryCatch(.with_optional_timeout({
        if (v$method == "QP") {
          options(km2ipd.cifresolve_strict_tick = v$strict_tick,
                  km2ipd.cifresolve_strict_dec = v$strict_dec,
                  km2ipd.cifresolve_constr_tol = v$constr_tol)
          reconstruct_twoarm_cifresolve(
            digi, sim$risk_table, total_events = sim$total_events,
            censor_marks_by_arm = sim$censor_marks_by_arm, configuration = v$configuration
          )
        } else {
          marks <- if (v$variant == "KMtoIPD_no_marks") {
            list(Control = numeric(), Treatment = numeric())
          } else sim$censor_marks_by_arm
          reconstruct_ipd_twoarm_final(
            digi, sim$risk_table, method = v$method,
            total_events = if (v$variant == "IPDfromKM_total_events") sim$total_events else NULL,
            censor_marks_by_arm = marks, kmtoipd_use_marks = TRUE
          )
        }
      }, timeout_sec = timeout), error = function(e) list(ipd = NULL, error = conditionMessage(e)))
      base <- cbind(data.table::copy(d), data.table::copy(v))
      base[, `:=`(rep = r, seed = seed, runtime_sec = proc.time()[[3]] - start,
                  ok = !is.null(rec$ipd), eval_ok = FALSE)]
      if (is.null(rec$ipd)) {
        base[, error := rec$error %||% "Unknown reconstruction error"]
      } else {
        metrics <- tryCatch(evaluate_reconstruction(
          sim$ipd, rec$ipd, digi$data, sim$risk_table, max(sim$risk_table$time)
        ), error = function(e) data.table::data.table(eval_ok = FALSE, error_eval = conditionMessage(e)))
        base[, eval_ok := NULL]
        base <- cbind(base, metrics)
        base[, relative_event_deviation_pct := if ("event_diff_total" %in% names(base)) {
          100 * abs(event_diff_total) / max(sum(sim$ipd$status == 1L), 1L)
        } else NA_real_]
      }
      base[, valid_success := ok %in% TRUE & eval_ok %in% TRUE]
      rows[[length(rows) + 1L]] <- base
    }
    data.table::fwrite(data.table::rbindlist(rows, fill = TRUE), file.path(out_dir, "results.csv"))
  }
  options(benchmark_config$method_options)
  invisible(data.table::rbindlist(rows, fill = TRUE))
}
