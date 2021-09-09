library(testthat)
context("parameterRecovery")

source("parameterRecovery.R")

test_that("parameterRecovery runs", {
  fName <- "data/test"
  unlink(fName, recursive = T)
  parameterRecovery(
    savePath = fName,
    nShuffles = 1,
    includeSplits = T,
    fixedParameters = c(NA, NA),
    customFilter = \(x) x[c(1, 400), ],
    nCores = parallel::detectCores() - 4,
    installPackages = F,
    verbosity = 2,
    remote_URL = ""
  )
  expect_equal(dir.exists(fName), TRUE)
})
