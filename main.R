#' Load the parameterRecovery function and run it.
source("/home/wolf5224/arc_thesis-parameter-estimation/parameterRecovery.R")
parameterRecovery(
  savePath = "/data/xpsy-acc/wolf5224/thesis-parameter-estimation.rda",
  nShuffles = 9,
  fixedParameters = c(NA, NA),
  customFilter = \(x) x,
  nCores = parallel::detectCores(),
  installPackages = TRUE,
  verbosity = 2
)