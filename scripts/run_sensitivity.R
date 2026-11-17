#!/usr/bin/env Rscript
source("src/load_methods.R")
source("src/simulation/sensitivity.R")
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  values <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(values)) default else sub(paste0("^--", name, "="), "", tail(values, 1))
}
read_config <- function(flag, default_path, required, columns) {
  path <- option(flag, if (file.exists(default_path)) default_path else "")
  if (!nzchar(path)) {
    if (required) stop("Supply --", flag, "=PATH for a full run.")
    return(NULL)
  }
  if (!file.exists(path)) stop("File not found: ", path)
  input <- data.table::fread(file = path)
  missing <- setdiff(columns, names(input))
  if (length(missing)) stop(path, " is missing columns: ", paste(missing, collapse = ", "))
  if (!nrow(input)) stop("No input rows in ", path)
  input
}

part <- option("part", "single-factor")
stopifnot(part %in% c("single-factor", "input-ablation"))
full <- identical(tolower(option("full", "false")), "true")
key <- if (part == "single-factor") "single_factor" else "input_ablation"
reps <- as.integer(option("reps", if (full) benchmark_config$repetitions[[key]] else 1L))
stopifnot(length(reps) == 1L, !is.na(reps), reps > 0L)
out <- option("output", file.path("outputs", paste0(key, if (full) "_full" else "_example")))
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty. Choose a new directory with --output=PATH.")
}
columns <- c("design_id", "scenario_id", "n_control", "n_treatment", "target_censoring",
             "censoring", "digitize_sd_S")
if (part == "single-factor") {
  columns <- c(columns, "phase", "design_label", "admin_time", "fixed_x_end", "digitize_mode")
}
design <- read_config("design", file.path("configs", paste0(key, "_design.csv")), full, columns)
if (is.null(design)) {
  design <- example_sensitivity_design(part)
  message("Using an illustrative sensitivity design.")
} else if (!full) {
  design <- design[1]
}
if (part == "input-ablation") {
  qp <- read_config("qp-variants", "configs/qp_input_variants.csv", full,
                    c("variant", "configuration", "information", "strict_tick", "strict_dec", "constr_tol"))
  variants <- input_variants(qp)
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(design, file.path(out, "design.csv"))
if (part == "single-factor") {
  run_single_factor(design, reps, out)
} else {
  run_auxiliary_inputs(design, reps, out, variants)
}
cat("Sensitivity outputs:", out, "\n")
