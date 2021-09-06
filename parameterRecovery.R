#' Fit participant data using the adviseR:: model.
#'
#' @param savePath Character string of filename to save data. This file will be
#'   loaded at the beginning so that progress can be picked up from where it
#'   left off.
#' @param nShuffles Number of times to shuffle and refit for each participant
#' @param includeSplits Whether to include a dataframe of split estimates for
#'   participants who saw multiple pairs of advisors
#' @param fixedParameters Length 2 vector of NA/number for whether
#'   weighted_selection and trust_update_rate should be fixed
#' @param customFilter function to apply to the data (can filter by UID,
#'   slice_sample, head, etc)
#' @param nCores Number of cores to use for parallel processing
#' @param installPackages Whether or not to install packages (useful for
#'   cluster/virtual environments)
#' @param remote_url URL to query with a participant UID to allow coordinatino
#'   of work across multiple computers
#' @param verbosity Number indicating how much output to print (higher = more
#'   verbose)
#'
#' @details The function returns nothing, but will save a file to
#'   \code{savePath} which contains three dataframes: recovered_parameters,
#'   shuffles, and status. These dataframes contain the recovered parameters
#'   (_end columns) for each participant, the results of the shuffling, and
#'   information on the timings and/or errors related to processing each
#'   participant's data.
#'
#'   Because of how the restarting works, participants are not guaranteed to
#'   have only a single row in the dataframes.
parameterRecovery <- function(
  savePath,
  nShuffles = 9,
  includeSplits = TRUE,
  fixedParameters = c(NA, NA),
  customFilter = function(x) x,
  nCores = parallel::detectCores(),
  installPackages = FALSE,
  remote_URL = "https://acclab.psy.ox.ac.uk/~mj221/ARC/",
  verbosity = 2
) {
  if (installPackages) {
    libPath <- '/home/wolf5224/R_libs4'
    # Libraries ---------------------------------------------------------------
    pkgs <- c(
      'tibble', 'dplyr', 'tidyr', 'lubridate', 'purrr', 'parallel', 'remotes'
    )
    pkgs <- pkgs[!(pkgs %in% installed.packages()[,1])]
    install.packages(pkgs, repos = 'http://cran.irsn.fr', lib = libPath)
    remotes::install_github('oxacclab/adviseR', lib = libPath)
    remotes::install_github('oxacclab/esmData', lib = libPath)
  }
  
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(stringr)
  library(parallel)
  library(esmData)
  library(adviseR)
  
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
    step_size = c(5, 1),
    min_step_size = c(.25, .005),
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
            coords[1] + step_size[1] / 2 * sign(E_diff[1]),
            coords[2] + step_size[2] / 2 * sign(E_diff[2])
          )
        } else {
          # If we're in a minimum we don't need to move
          if (all(E_diff == 0)) {
            return(out)
          }
          coords_new <- c(
            coords[1] + sign(E_diff[1]) * step_size[1] * abs(E_diff[1]) / sum(abs(E_diff)),
            coords[2] + sign(E_diff[2]) * step_size[2] * abs(E_diff[2]) / sum(abs(E_diff))
          )
        } 
        # Calculate new error
        E_new <- adviseR::simulateFromData(
          x, tibble(w = coords_new[1], LR = coords_new[2])
        )
        
        # Check we got better with the movement
        if (all(E_new < E)) {
          okay <- T
        } else {
          if (all(step_size <= min_step_size)) {
            return(out)
          }
          step_size[E_new >= E] <- pmax(
            step_size[E_new >= E] / 2,
            min_step_size[E_new >= E]
          )
        }
        
        if (all(E <= E_new)) {
          if (all(step_size <= min_step_size)) {
            return(out)
          }
          step_size <- pmax(step_size / 2, min_step_size)
        } else {
          okay <- T
        }
      }
      coords <- coords_new
      E <- E_new
      
      i <- i + 1
    }
    
    if (stall_count < max_steps)
      warning(paste0('GradientDescent failed to converge after ', max_steps, ' steps'))
    out
  }
  
  #' Perform gradient descent with some coordinates fixed
  gdFixed <- function(
    x, 
    start_coords = NULL, 
    starting_limits = tibble::tibble(min = c(-5, -.5), max = c(30, .5)), 
    max_steps = 200,
    step_size = c(5, 1),
    min_step_size = c(.25, .005),
    stall_count = max_steps, 
    fixed = c(NA, NA)
  ) {
    f <- function(x) ifelse(is.na(fixed), x, fixed)
    if (!is.null(start_coords))
      start_coords = f(start_coords)
    starting_limits = tibble::tibble(
      min = f(starting_limits$min), 
      max = f(starting_limits$max)
    )
    step_size = ifelse(is.na(fixed), step_size, 0)
    min_step_size = ifelse(is.na(fixed), min_step_size, 0)
    gradientDescent(
      x = x, 
      start_coords = start_coords,
      starting_limits = starting_limits,
      max_steps = max_steps,
      step_size = step_size,
      min_step_size = min_step_size,
      stall_count = stall_count
    )
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
    
    f <- function(d, ...) {
      library(tidyr); library(dplyr); library(purrr)
      d %>%
        mutate(gd = map(data, function(z) gdFixed(z, ...) %>% gradientDescentSummary())) %>%
        unnest(cols = gd)
    }
    
    if (nCores > 1) {
      cl <- makeCluster(nCores)
      clusterExport(
        cl, 
        c('gradientDescent', 'gdFixed', 'gradientDescentSummary'),
        parent.env(environment())
      )
      on.exit(stopCluster(cl), add = T)
      
      gd <- parLapply(cl = cl, X = gd, fun = f, ...)
    } else {
      gd <- lapply(gd, f, ...)
    }
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
  
  # Pick up where we left off
  tryCatch(load(savePath), error = function(e) {})
  
  # Script ------------------------------------------------------------------
  
  # Load up some dots task data. Like all of it.
  select_experiment('dotstask')
  trials <- trials %>% 
    mutate(uid = factor(paste0(studyId, studyVersion, ' p', pid)))
  
  d <- trials %>% 
    nest(data = -uid) %>%
    mutate(okay = map_lgl(data, ~ any(.$hasChoice, na.rm = T))) %>%
    filter(okay) %>%
    customFilter()
  
  if (verbosity) {
    print(paste0('Using ', nCores, ' cores.'))
    print(paste0('Found ', nrow(d), ' participants\'s data.'))
  }
  
  if ('recovered_parameters' %in% ls()) {
    if (verbosity) print(paste0(nrow(status), ' ids complete.'))
    
    ids_left <- d$uid[!(d$uid %in% status$uid)]
    
    if (verbosity) print(paste0(length(ids_left), ' ids left'))
  } else {
    ids_left <- d$uid
    recovered_parameters <- NULL
    shuffles <- NULL
    status <- NULL
    splits <- NULL
  }
  
  if (length(ids_left) > 0) {
    # Cycle through remaining ids and calculate the model fit and shuffle position
    
    while (length(ids_left) > 0) {
      t1 <- Sys.time()
      id <- ids_left[1]
      ids_left <- sample(ids_left[ids_left != id])
      
      # Lookup id
      if (nchar(remote_URL) > 0) {
        req <- curlGetHeaders(
          paste0(remote_URL, "data/", URLencode(as.character(id))), 
          verify = F
        )
        if (attr(req, "status") == 200) {  # already exists
          if (verbosity > 1) print(paste0('Skipping id ', id, ' (done remotely)'))
          next()
        } else { # record that we're working on it
          req <- curlGetHeaders(
            paste0(remote_URL, "?f=", URLencode(as.character(id))), 
            verify = F
          )
          if (attr(req, "status") == 200) {
            if (verbosity > 1) print(paste0('Reserved id ', id))
          } else {
            grr <- paste0(
              "Failed to reserve id ", id, " at ", 
              remote_URL, "?f=", URLencode(as.character(id))
            )
            warning(grr)
            if (verbosity > 0) print(grr)
            if (verbosity > 1) print(req)
          }
        }
      } else
        if (verbosity > 1) print(paste0('Processing id ', id))
      
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
                ),
                okay = case_when(
                  is.na(choice0) ~ T,
                  is.na(choice1) ~ T,
                  advisorIndex %in% c(choice0, choice1) ~ T,
                  T ~ F
                )
              )
            ),
            choices_okay = map_lgl(data, ~ all(.$okay))
          )
        
        if (!x$choices_okay) {
          stop("AdvisorIndices and choices are mismatched.")
        }
        
        x <- x %>% 
          mutate(data = map(data, ~ select(., -okay))) %>% 
          select(-choices_okay)
        
        f <- function(x, ...) {
          gd <- doGradientDescent(x, nCores = nCores, ...) 
          if (!is.na(fixedParameters[1])) {
            filter(gd, tu_error == min(tu_error))
          } else if (!is.na(fixedParameters[2])) {
            filter(gd, ws_error == min(ws_error))
          } else 
            gd %>% 
            filter(
              scaledError(ws_error, tu_error) == min(scaledError(ws_error, tu_error))
            )
        }
        
        gd <- f(x, fixed = fixedParameters)
        
        recovered_parameters <- bind_rows(recovered_parameters, gd)
        
        if (nShuffles > 0) {
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
        }
        
        # Splits
        if (includeSplits) {
          choices <- x$data[[1]] %>%
            group_by(choice0, choice1) %>%
            summarise(.groups = "drop") %>%
            mutate(
              choice = paste0(
                pmin(choice0, choice1), "_", pmax(choice0, choice1)
              )
            ) %>% 
            drop_na() %>%
            select(choice) %>%
            unique() %>%
            rowid_to_column("choiceNum") %>%
            mutate(
              c1 = str_extract(choice, "^[0-9]"),
              c2 = str_extract(choice, "[0-9]$")
            ) %>%
            pivot_longer(c(c1, c2)) %>%
            select(choiceNum, advisorIndex = value) %>%
            mutate(across(everything(), as.numeric))
          
          y <- x %>%
            mutate(
              data = map(
                data, 
                ~ left_join(. , choices, by = "advisorIndex")
              )
            ) %>%
            unnest(cols = data) %>%
            mutate(uid = paste0(uid, "_", choiceNum)) %>% 
            filter(!is.na(choiceNum)) %>%
            nest(data = -c(uid, okay))
          
          if (nrow(y) > 1) {
            for (i in unique(y$uid)) {
              splits <- bind_rows(
                splits, 
                f(
                  y %>% 
                    filter(uid == i) %>% 
                    mutate(data = map(data, ~ select(., -choiceNum))), 
                  fixed = fixedParameters
                )
              )
            }
          }
        }
        
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
        if (verbosity) print(e$message)
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
      
      if (verbosity > 1) print(status[nrow(status),])
      
      # Write to the cache 
      save(
        recovered_parameters, 
        shuffles,
        status,
        splits,
        file = savePath
      )
      
      if (verbosity > 1) print('saved.')
    }
  }
}