library(dplyr)

tmpEnv <- new.env()
recovered_parameters <- NULL
shuffles <- NULL
status <- NULL
splits <- NULL

for (f in list.files(
  "data/thesis-parameter-estimation", 
  full.names = T,
  pattern = "\\.Rdata$"
)) {
  load(f, envir = tmpEnv)
  recovered_parameters <- bind_rows(recovered_parameters, tmpEnv$recovered_parameters)
  shuffles <- bind_rows(shuffles, tmpEnv$shuffles)
  status <- bind_rows(status, tmpEnv$status)
  splits <- bind_rows(splits, tmpEnv$splits)
}

save(
  recovered_parameters, shuffles, status, splits, 
  file = "data/thesis-parameter-estimation.rda"
)
