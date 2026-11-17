# Shared initialization for entry points executed from the repository root.
if (!file.exists("configs/benchmark.R")) {
  stop("Run this command from the repository root.")
}
if (dir.exists(".Rlib")) .libPaths(c(normalizePath(".Rlib"), .libPaths()))
source("configs/benchmark.R")
options(stringsAsFactors = FALSE)
