#!/usr/bin/env Rscript
# Run a small example by default; select --full=true for the complete design.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(value)) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", tail(value, 1))
}
part <- option("part", "single")
full <- identical(tolower(option("full", "false")), "true")
reps <- as.integer(option("reps", "1"))
stopifnot(reps > 0L, part %in% c("single", "meta"))
source("src/load_methods.R")
methods <- benchmark_config$methods
out <- option("output", file.path("outputs", paste0("simulation_", part, if (full) "_full" else "")))
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty. Choose a new directory with --output=PATH.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (part == "single") {
  design <- build_single_trial_design("all")
  # The replicate seed uses the design index within each survival-model phase.
  design[, design_id := seq_len(.N), by = phase]
  if (!full) design <- design[1]
  result <- run_single_trial_benchmark(design,
    n_rep = reps, methods = methods,
    out_dir = out, seed_base = benchmark_config$seeds$single,
    method_timeout_sec = benchmark_config$timeout_sec
  )
} else {
  source("src/simulation/meta_analysis.R")
  preset <- option("preset", "main")
  stopifnot(preset %in% c("main", "stress"))
  design <- if (full) {
    make_meta_grid(preset)
  } else {
    make_meta_grid(
      preset = "custom", K_vec = 3, n_per_arm_vec = 40, p_miss_vec = c(0, .8),
      tau_logHR_vec = .05, tau_label_vec = "low", censoring_vec = "random",
      target_censoring_vec = .30, digitize_sd_S_vec = .01, weib_shape_vec = 1
    )
  }
  result <- run_meta_simulation(
    grid = design, n_rep = reps, methods = methods,
    mu_HR = benchmark_config$meta_mean_hr,
    seed_base = benchmark_config$seeds$meta, verbose = TRUE
  )
  for (name in names(result)) {
    if (is.data.frame(result[[name]])) {
      data.table::fwrite(result[[name]], file.path(out, paste0(name, ".csv")))
    }
  }
}
cat("Simulation outputs:", out, "\n")
