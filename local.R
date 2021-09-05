#' Load the parameterRecovery function and run it.
source("parameterRecovery.R")
parameterRecovery(
  savePath = "data/thesis-parameter-estimation-full.rda",
  nShuffles = 9,
  includeSplits = T,
  fixedParameters = c(NA, NA),
  customFilter = \(x) x,
  nCores = parallel::detectCores() - 4,
  installPackages = F,
  verbosity = 2
)
