# Load the four reconstruction workflows and their benchmark settings.
source('scripts/bootstrap.R')
required <- c('data.table','survival','IPDfromKM','KMtoIPD','reconstructKM','CIFresolve','R.utils')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
source('src/simulation/single_trial.R')
source('src/methods/ikm.R')
source('src/methods/qp.R')
options(benchmark_config$method_options)
