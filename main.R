libPath <- '/home/wolf5224/R_libs4'
# Libraries ---------------------------------------------------------------
install.packages(c(
  'tibble', 'dplyr', 'tidyr', 'lubridate', 'purrr', 'parallel', 'remotes'
), repos = 'http://cran.irsn.fr', lib = libPath)
remotes::install_github('oxacclab/adviseR', lib = libPath)
remotes::install_github('oxacclab/esmData', lib = libPath)

library(tibble)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(parallel)

# Functions ---------------------------------------------------------------

#' Perform a gradient descent search on a participant's data
#' @param x tbl of data
#' @param start_coords c(weightedSelection, trustUpdate) starting coordinates.
#'  Selected randomly if not specified.
#' @param limits $max and $min are both tuples of c(weightedSelection, 
#'  trustUpdate) values. Governs the starting limits rather than the space
#'  explored.
#' @param max_steps maximum number of iterations
#' @param step_size size of steps to begin with (normalised across dimensions)
#' @param min_step_size minimum step size after which to stop
#' @param stall_count number of consecutive repetitions of values to trigger
#'  early termination
gradientDescent <- function(
  x, 
  start_coords = NULL, 
  starting_limits = tibble::tibble(min = c(-5, -.5), max = c(30, .5)), 
  max_steps = 200,
  step_size = 2,
  min_step_size = .005,
  stall_count = max_steps
) {
  out <- NULL
  # Initialize state
  i <- 1
  E <- c(
    `Advisor choice mean squared error` = Inf, 
    `Advice-taking mean squared error` = Inf
  )
  coords <- if (is.null(start_coords)) {
    # randomly selected start position
    runif(2, min = starting_limits$min, max = starting_limits$max)  
  } else { 
    start_coords
  }
  
  # Gradient descent
  while (i <= max_steps) {
    out <- bind_rows(
      out, 
      tibble(
        step = i, 
        weightedSelection = coords[1],
        trustUpdate = coords[2], 
        weightedSelection_error = E[1],
        trustUpdate_error = E[2]
      )
    )
    # Finished if the last stall_count steps are identical or flip between 2
    if (i > stall_count) {
      if (length(unique(tail(out$trustUpdate, stall_count))) < 3 &
          length(unique(tail(out$weightedSelection, stall_count))) < 3)
        return(out)
    }
    
    # Find new coordinates
    okay <- F
    while (!okay) {
      coords_new <- coords + step_size
      E_new <- adviseR::simulateFromData(
        x, tibble(w = coords_new[1], LR = coords_new[2])
      )
      E_diff <- E - E_new
      # Step in each direction according to the amount of error reduction
      if (is.infinite(sum(E_diff))) {
        coords_new <- c(
          coords[1] + step_size / 2 * sign(E_diff[1]),
          coords[2] + step_size / 2 * sign(E_diff[2])
        )
      } else {
        # If we're in a minimum we don't need to move
        if (E_diff[1] == 0 & E_diff[2] == 0) {
          return(out)
        }
        coords_new <- c(
          coords[1] + sign(E_diff[1]) * step_size * abs(E_diff[1]) / sum(abs(E_diff)),
          coords[2] + sign(E_diff[2]) * step_size * abs(E_diff[2]) / sum(abs(E_diff))
        )
      } 
      # Calculate new error
      E_new <- adviseR::simulateFromData(
        x, tibble(w = coords_new[1], LR = coords_new[2])
      )
      
      # Check we got better with the movement
      if (E[1] < E_new[1] & E[2] < E_new[2]) {
        if (step_size < min_step_size) {
          return(out)
        }
        step_size <- step_size / 2
      } else {
        okay <- T
      }
    }
    coords <- coords_new
    E <- E_new
    
    i <- i + 1
  }
  
  if (stall_count < max_steps)
    warning(glue('GradientDescent failed to converge after { max_steps } steps'))
  out
}

#' Parellelised implementation of gradientDescent over nRuns runs
#' @param x dataframe to perform descent on 
#' @param nRuns number of initial start positions to try
#' @param nCores number of cores to use
#' @dotparams gradientDescent
doGradientDescent <- function(
  x, 
  nRuns = 100,
  nCores = parallel::detectCores() / 2, 
  ...
) {
  if (!has_name(x, 'uid'))
    x <- mutate(x, uid = 'sim')
  if (!has_name(x, 'data'))
    x <- nest(x, data = -uid)
  gd <- x %>%
    crossing(tibble(run = 1:nRuns)) %>%
    mutate(core = rep(1:nCores, ceiling(nRuns / nCores))[1:nRuns]) %>%
    group_by(core) %>%
    group_split() 
  
  cl <- makeCluster(nCores)
  clusterExport(cl, c('gradientDescent', 'gradientDescentSummary'))
  on.exit(stopCluster(cl), add = T)
  f <- function(d, ...) {
    library(tidyr); library(dplyr); library(purrr)
    d %>%
      mutate(gd = map(data, function(z) gradientDescent(z, ...) %>% gradientDescentSummary())) %>%
      unnest(cols = gd)
  }
  
  gd <- parLapply(cl = cl, X = gd, fun = f, ...)
  bind_rows(gd)
}

gradientDescentSummary <- function(gd) {
  tibble::tibble(
    ws_start = gd$weightedSelection[1],
    tu_start = gd$trustUpdate[1],
    ws_end = gd$weightedSelection[nrow(gd)],
    tu_end = gd$trustUpdate[nrow(gd)],
    ws_error = gd$weightedSelection_error[nrow(gd)],
    tu_error = gd$trustUpdate_error[nrow(gd)],
    step_count = nrow(gd)
  )
}

#' Return err1 and err2 scaled for appropriate combination
#' @param err1 numeric vector
#' @param err2 numeric vector of length length(err1)
#' @return numeric vector of length length(err1) with combined error values
scaledError <- function(err1, err2) {
  x <- scale(err1); y <- scale(err2)
  x <- x - min(x); y <- y - min(y)
  log(x + y)
}

# Variables ---------------------------------------------------------------

nShuffles <- 9
nCores <- as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
cacheFile <- '/data/xpsy-acc/wolf5224/thesis-parameter-estimation.rda'
# cacheFile <- 'data/thesis-parameter-estimation.rda'
tryCatch(
  load(cacheFile),
  error = function(e) {
    recovered_parameters <- NULL
    shuffles <- NULL
    status <- NULL
  }
)

# Script ------------------------------------------------------------------

# Load up some dots task data. Like all of it.
esmData::select_experiment('dotstask')
trials <- trials %>% 
  mutate(uid = factor(paste0(studyId, studyVersion, ' p', pid)))

d <- trials %>% 
  nest(data = -uid) %>%
  mutate(okay = map_lgl(data, ~ any(.$hasChoice, na.rm = T))) %>%
  filter(okay)

ids_left <- d$uid[!(d$uid %in% status$uid)]

if (length(ids_left) > 0) {
  # Cycle through remaining ids and calculate the model fit and shuffle position
  
  while (length(ids_left) > 0) {
    t1 <- Sys.time()
    id <- ids_left[1]
    
    tryCatch({
      # Format how the C++ code wants the data:
      # * @param trials a data frame of trials with 5 columns (names may vary):
      # * initialConfidence - raw initial confidence rating, all values +ve
      # * advisorIndex - index of the advisor chosen (0 or 1; NA if no choice made)
      # * choice0 - index of the first advisor offered in the choice (NA if no choice offered)
      # * choice1 - index of the second advisor offered in the choice (NA if no choice offered)
      # * advisorAgrees - whether the chosen advisor agrees (NA if no advice provided)
      # * finalConfidence - raw final confidence rating, -ve values indicate change-of-mind
      x <- d %>%
        filter(uid == id) %>%
        mutate(
          data = map(data, ~select(
            .,
            initialConfidence, 
            advisorId, 
            choice0,
            choice1,
            advisorAgrees, 
            finalConfidence
          ) %>%
            mutate(
              finalConfidence = if_else(
                sign(initialConfidence) == sign(finalConfidence),
                abs(finalConfidence), 
                -1 * abs(finalConfidence)
              ),
              initialConfidence = abs(initialConfidence),
              advisorIndex = as.numeric(factor(advisorId)), # avoid 0-indexing
              across(
                .cols = matches('choice[01]'), 
                ~ as.numeric(factor(., levels = levels(factor(advisorId))))
              )
            )
          )
        )
      
      f <- function(x, ...) {
        doGradientDescent(x, nCores = nCores, ...) %>% 
          filter(
            scaledError(ws_error, tu_error) == min(scaledError(ws_error, tu_error))
          )
      }
      
      gd <- f(x)
      
      recovered_parameters <- bind_rows(recovered_parameters, gd)
      
      shuffle_list <- crossing(
        gd %>%
          select(uid, data, weightedSelection = ws_end, trustUpdateRate = tu_end), 
        shuffle_run = 1:nShuffles
      ) %>%
        mutate(data = map(data, function(x) {
          x$advisorAgrees <- sample(x$advisorAgrees)
          x
        }))
      # Run gradient descent on shuffles
      ### WARNING: this takes much longer than the original run because there are 
      ### multiple shuffles!
      dList <- list()
      for (i in 1:nrow(shuffle_list))
        dList[[i]] <- shuffle_list[i, ]
      
      shuffle_list <- lapply(dList, f)
      shuffles <- bind_rows(
        shuffles,
        bind_rows(shuffle_list)
      )
      
      status <- bind_rows(
        status,
        tibble(
          uid = id, 
          time_start = t1, 
          time_end = Sys.time(), 
          error = NA_character_, 
          cores = nCores
        )
      )
    },
    error = function(e) {
      status <- bind_rows(
        status,
        tibble(
          uid = id, 
          time_start = t1, 
          time_end = Sys.time(), 
          error = e$message, 
          cores = nCores
        )
      )
    })
    
    # Write to the cache 
    save(
      recovered_parameters, 
      shuffles,
      status,
      file = cacheFile
    )
    
    ids_left <- sample(ids_left[ids_left != id])
  }
}
