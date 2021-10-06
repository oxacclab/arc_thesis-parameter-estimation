# ARC Thesis Parameter Estimation

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/360914235.svg)](https://zenodo.org/badge/latestdoi/360914235)
<!-- badges: end -->

These scripts are used to run the parameter recovery components of the [Exploring Social Metacognition](https://github.com/mjaquiery/oxforddown) DPhil thesis.
The main work is done within the `parameterRecovery()` function, created in `parameterRecovery.R`.
This (mammoth) function is invoked by other scripts.
For local computer running, I used `local.R`, and for running on the Advanced Research Computing (ARC) cluster, I used `main.R`, scheduled using `slurm.sh`.
I used `analysis.R` to check results while developing the code, and `test.R` to test it.

The `parameterRecovery()` function sends the names of the files being processed to a remote server (if requested) to enable coordination across multiple instances. 
When all the data have been gathered up and placed in the data folder (and appropriately named), they can be boiled down to a summary using `collate.R`.
