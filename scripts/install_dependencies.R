#!/usr/bin/env Rscript
# Install the direct dependencies of the four simulation workflows.
source("scripts/bootstrap.R")
dir.create(".Rlib", showWarnings = FALSE)
lib <- normalizePath(".Rlib")
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c(
  data.table = "1.17.8", dplyr = "1.2.1", survival = "3.8-3",
  R.utils = "2.13.0", metafor = "5.0-1",
  reconstructKM = "0.4.0", IPDfromKM = "0.1.10"
)
github <- c(
  KMtoIPD = "rogula/KMtoIPD@80e9ba067492aabc6c0ac37409965186c13224de",
  CIFresolve = "andrewtitman/CIFresolve@2295c0d32b3bb5595fb38a32030f773daa47091d"
)
versions <- c(cran, KMtoIPD = "1.0.0", CIFresolve = "0.1.1")

matches <- function(pkg) {
  description <- tryCatch(packageDescription(pkg), warning = function(w) NULL)
  if (is.null(description) || package_version(description$Version) != package_version(versions[[pkg]])) {
    return(FALSE)
  }
  if (pkg %in% names(github) && grepl("@", github[[pkg]], fixed = TRUE) && !is.null(description$RemoteSha)) {
    return(identical(description$RemoteSha, sub(".*@", "", github[[pkg]])))
  }
  TRUE
}

needed <- names(versions)[!vapply(names(versions), matches, logical(1))]
if (length(needed) && !requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = lib)
}
for (pkg in needed) {
  if (pkg %in% names(github)) {
    remotes::install_github(github[[pkg]], lib = lib, upgrade = "never", dependencies = NA)
  } else {
    remotes::install_version(pkg, version = cran[[pkg]], lib = lib, upgrade = "never", dependencies = NA)
  }
}

loadable <- function(pkg) matches(pkg) && requireNamespace(pkg, quietly = TRUE)
failed <- names(versions)[!vapply(names(versions), loadable, logical(1))]
if (length(failed)) stop("Could not load the required package versions: ", paste(failed, collapse = ", "))
cat("All simulation dependencies are available.\n")
