#' Load the parameterRecovery function and run it.
source("parameterRecovery.R")
parameterRecovery(
  savePath = "data/thesis-parameter-estimation",
  nShuffles = 9,
  includeSplits = T,
  fixedParameters = c(NA, NA),
  customFilter = \(x) x,
  nCores = parallel::detectCores(),
  installPackages = F,
  verbosity = 2
)
