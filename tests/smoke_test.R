#!/usr/bin/env Rscript
# Exercise all four implementations and both simulation parts.
source("src/load_methods.R")
source("src/simulation/meta_analysis.R")
design <- build_single_trial_design("all")
design[, design_id := seq_len(.N), by = phase]
stopifnot(nrow(design) == 61L, identical(design[phase == "phase2"]$design_id, 1:26))
stopifnot(nrow(make_meta_grid("main")) == 45L, nrow(make_meta_grid("stress")) == 16L)
single <- run_one_single_trial_rep(design[1], 1L, methods = benchmark_config$methods,
                                 seed_base = benchmark_config$seeds$single)
stopifnot(nrow(single) == 4L, all(single$ok), all(single$eval_ok),
          identical(sort(single$method), sort(benchmark_config$methods)))
grid <- make_meta_grid(preset = "custom", K_vec = 3, n_per_arm_vec = 40,
                      p_miss_vec = c(0, 0.8), tau_logHR_vec = 0.05,
                      tau_label_vec = "low", censoring_vec = "random",
                      target_censoring_vec = 0.30, digitize_sd_S_vec = 0.01,
                      weib_shape_vec = 1)
meta <- run_meta_simulation(grid, n_rep = 1L, methods = benchmark_config$methods,
                            seed_base = benchmark_config$seeds$meta, verbose = FALSE)
stopifnot(nrow(meta$meta_methods_long) == 22L,
          nrow(meta$meta_trial_methods_long) == 12L,
          !anyDuplicated(meta$meta_methods_long[, .(base_id, rep, p_miss, meta_strategy)]))
dir.create("outputs/tests", recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(single, "outputs/tests/single_trial.csv")
data.table::fwrite(meta$meta_methods_long, "outputs/tests/meta_methods.csv")
writeLines(capture.output(sessionInfo()), "outputs/tests/sessionInfo.txt")
cat("PASS: all four methods, phase-specific design IDs, and both simulation parts.\n")
