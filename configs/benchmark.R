# Frozen benchmark settings. Run entry points from the repository root.
benchmark_config <- list(
  methods = c("iKM", "IPDfromKM", "KMtoIPD", "QP"),
  seeds = list(single = 20260611L, meta = 20260619L,
               single_factor = 20260711L, input_ablation = 202607210L),
  repetitions = list(single = 1000L, meta = 300L,
                     single_factor = 300L, input_ablation = 200L),
  timeout_sec = 180,
  sensitivity_ikm_timeout_sec = 10,
  sensitivity_qp_timeout_sec = 30,
  meta_mean_hr = 0.75,
  method_options = list(
    km2ipd.reconstructkm_input_mode = "direct_augmented",
    km2ipd.reconstructkm_timeout = 30,
    km2ipd.cifresolve_configuration = "risk_ticks",
    km2ipd.cifresolve_strict_tick = TRUE,
    km2ipd.cifresolve_strict_dec = FALSE,
    km2ipd.cifresolve_constr_tol = 1e-8,
    km2ipd.cifresolve_cen_penalty = 1e-3
  )
)
