library(testthat)
context("parameterRecovery")

source("parameterRecovery.R")

testthat("parameterRecovery runs", {
  fName <- "data/test.rda"
  unlink(fName)
  expect_warning(parameterRecovery(
    savePath = fName,
    nShuffles = 1,
    fixedParameters = c(NA, NA),
    customFilter = \(x) head(x, 2),
    nCores = parallel::detectCores() - 4,
    installPackages = F,
    verbosity = 2
  ))
  expect_equal(file.exists(fName), TRUE)
})