library(testthat)
context("parameterRecovery")

source("parameterRecovery.R")

test_that("parameterRecovery runs", {
  fName <- "data/test.rda"
  unlink(fName)
  expect_warning(parameterRecovery(
    savePath = fName,
    nShuffles = 1,
    includeSplits = T,
    fixedParameters = c(NA, NA),
    customFilter = \(x) x[c(1, 400), ],
    nCores = parallel::detectCores() - 4,
    installPackages = F,
    verbosity = 2
  ))
  expect_equal(file.exists(fName), TRUE)
})
