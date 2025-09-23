#' Helper functions for binary variable regression sensitivity analysis
#'
#' These functions perform core computations such as flipping observations,
#' calculating extreme bounds, and running permutation tests.
#' They are intended for internal use and are not exported.
#'
#' @keywords internal, helper, utility
#' @name utils
NULL

# Remove explicit require call - dependencies handled via DESCRIPTION

##### Helper Functions ####

#' Validate inputs for misclassification sensitivity analysis
#'
#' This function performs comprehensive validation of inputs to ensure
#' they meet the requirements for sensitivity analysis functions.
#'
#' @param dat Data frame containing the data
#' @param outcome Character. Name of outcome variable
#' @param treatment Character. Name of treatment variable
#' @param binary Character. Name of binary variable subject to misclassification
#' @param m Model object
#' @param nsims Integer. Number of simulations (optional)
#' @param tol Numeric. Tolerance value (optional)
#' @param R_vect Numeric vector. Reclassification rates (optional)
#' @keywords internal
validate_misclass_inputs <- function(dat, outcome, treatment, binary, m,
                                    nsims = NULL, tol = NULL, R_vect = NULL) {
  # Check data frame
  if(!inherits(dat, "data.frame")) {
    stop("dat must be a data.frame or tibble")
  }

  if(nrow(dat) == 0) {
    stop("dat cannot be empty")
  }

  # Check variable names
  if(!is.character(outcome) || length(outcome) != 1) {
    stop("outcome must be a single character string")
  }
  if(!outcome %in% names(dat)) {
    stop("outcome variable '", outcome, "' not found in data")
  }

  if(!is.character(treatment) || length(treatment) != 1) {
    stop("treatment must be a single character string")
  }
  if(!treatment %in% names(dat)) {
    stop("treatment variable '", treatment, "' not found in data")
  }

  if(!is.character(binary) || length(binary) != 1) {
    stop("binary must be a single character string")
  }
  if(!binary %in% names(dat)) {
    stop("binary variable '", binary, "' not found in data")
  }

  # Check binary variable is actually binary
  binary_vals <- unique(dat[[binary]])
  binary_vals <- binary_vals[!is.na(binary_vals)]
  if(!all(binary_vals %in% c(0, 1))) {
    stop("binary variable '", binary, "' must contain only 0s and 1s")
  }

  # Check model object
  if(!inherits(m, c("lm", "fixest"))) {
    stop("m must be a model object of class 'lm' or 'fixest'")
  }

  # Check optional parameters
  if(!is.null(nsims)) {
    if(!is.numeric(nsims) || length(nsims) != 1 || nsims < 1 || nsims != round(nsims)) {
      stop("nsims must be a positive integer")
    }
  }

  if(!is.null(tol)) {
    if(!is.numeric(tol) || length(tol) != 1 || tol <= 0 || tol >= 1) {
      stop("tol must be a numeric value between 0 and 1")
    }
  }

  if(!is.null(R_vect)) {
    if(!is.numeric(R_vect) || any(R_vect < 0) || any(R_vect > 1)) {
      stop("R_vect values must be numeric between 0 and 1")
    }
  }

  invisible(TRUE)
}

#' Validate inputs for contour plot functions
#'
#' This function validates inputs specific to contour plotting functions.
#'
#' @param dat Data frame containing the data
#' @param m Model object
#' @param treatment Character. Name of treatment variable
#' @param outcome Character. Name of outcome variable
#' @param binary Character. Name of binary variable
#' @param ground_truth Character. Name of ground truth variable (optional)
#' @param rho_M_Y Numeric. Correlation parameter (optional)
#' @param rho_M_D Numeric. Correlation parameter (optional)
#' @param nsims Integer. Number of simulations (optional)
#' @keywords internal
validate_contour_inputs <- function(dat, m, treatment, outcome, binary,
                                   ground_truth = NULL, rho_M_Y = NULL,
                                   rho_M_D = NULL, nsims = NULL) {
  # Use base validation first
  validate_misclass_inputs(dat, outcome, treatment, binary, m, nsims)

  # Check ground truth if provided
  if(!is.null(ground_truth)) {
    if(!is.character(ground_truth) || length(ground_truth) != 1) {
      stop("ground_truth must be a single character string")
    }
    if(!ground_truth %in% names(dat)) {
      stop("ground_truth variable '", ground_truth, "' not found in data")
    }
    # Check it's binary
    gt_vals <- unique(dat[[ground_truth]])
    gt_vals <- gt_vals[!is.na(gt_vals)]
    if(!all(gt_vals %in% c(0, 1))) {
      stop("ground_truth variable '", ground_truth, "' must contain only 0s and 1s")
    }
  }

  # Check correlation parameters
  if(!is.null(rho_M_Y)) {
    if(!is.numeric(rho_M_Y) || length(rho_M_Y) != 1 || abs(rho_M_Y) > 1) {
      stop("rho_M_Y must be a numeric value between -1 and 1")
    }
  }

  if(!is.null(rho_M_D)) {
    if(!is.numeric(rho_M_D) || length(rho_M_D) != 1 || abs(rho_M_D) > 1) {
      stop("rho_M_D must be a numeric value between -1 and 1")
    }
  }

  invisible(TRUE)
}

#' Extract control variables from a model object
#'
#' This function standardizes the extraction of control variables from different
#' model types (lm, fixest) to ensure consistency across the package.
#'
#' @param m Model object (lm, fixest, etc.)
#' @param treatment Character. Name of the treatment variable to exclude
#' @return Character vector of control variable names
#' @keywords internal
extract_controls <- function(m, treatment) {
  if(inherits(m, 'lm')) {
    rhs <- variable.names(m)
  } else if(inherits(m, 'fixest')) {
    # For fixest models, extract variables from formula
    rhs <- all.vars(formula(m))[-1]  # Remove outcome variable (first element)
  } else {
    rhs <- variable.names(m)
    warning('Untested model class: ', class(m)[1], '. Using variable.names() method.')
  }

  # Remove intercept and treatment variable, using exact matching for treatment
  return(rhs[!grepl(paste0('^(\\()?Intercept(\\))?$|^', treatment, '$'), rhs)])
}

##### Progress Bar Functions ####

#' Progress bar version of sapply
#'
#' This function creates a self-updating progress bar akin to that of `pbsapply` to
#' indicate percentage and time passed of long iterations. Results are simplified
#' to array format.
#'
#' @param X Vector or list to apply function over
#' @param FUN Function to apply to each element of X
#' @param ... Additional arguments passed to FUN
#' @return Array of results from applying FUN to each element of X
#' @keywords internal
prog_bar_sapply <- function(X, FUN, ...) {
  n <- length(X)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  start_time <- Sys.time()
  result <- vector("list", n)

  for (i in seq_along(X)) {
    result[[i]] <- FUN(X[[i]], ...)
    # Update progress bar less frequently for large datasets
    if (i %% max(1, n %/% 100) == 0 || i == n) {
      setTxtProgressBar(pb, i)
    }
  }

  close(pb)
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  cat(sprintf(" elapsed = %.2f seconds\n", elapsed))

  simplify2array(result, higher = TRUE)
}

#' Progress bar version of lapply
#'
#' This function creates a self-updating progress bar akin to that of `pblapply` to
#' indicate percentage and time passed of long iterations. Results are returned
#' as a list.
#'
#' @param X Vector or list to apply function over
#' @param FUN Function to apply to each element of X
#' @param ... Additional arguments passed to FUN
#' @return List of results from applying FUN to each element of X
#' @keywords internal
prog_bar_lapply <- function(X, FUN, ...) {
  n <- length(X)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  start_time <- Sys.time()
  result <- vector("list", n)

  for (i in seq_along(X)) {
    result[[i]] <- FUN(X[[i]], ...)
    # Update progress bar less frequently for large datasets
    if (i %% max(1, n %/% 100) == 0 || i == n) {
      setTxtProgressBar(pb, i)
    }
  }

  close(pb)
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  cat(sprintf(" elapsed = %.2f seconds\n", elapsed))

  result
}


##### Sensitivity Analysis #####

#### Extreme Bounds Flip ####

#' Find extreme coefficient after k flips
#'
#' This function identifies which observations will return the most extreme positive
#' or negative relationship between a continuous variable (cont) and a binary variable (bin)
#' by strategically flipping k observations to maximize or minimize the treatment coefficient.
#'
#' @param cont Numeric vector. Continuous variable (residuals from model)
#' @param bin Numeric vector. Binary variable (treatment from model)
#' @param k Integer. Number of observations to flip
#' @param treatment Character. Name of treatment variable in model
#' @param dat Data frame. The dataset
#' @param setting Character. Whether variable is treatment ('d') or dependent ('y')
#' @param direction Character. Direction to optimize ('max' or 'min')
#' @param m Model object. The fitted model
#' @return List containing the extreme coefficient and flip information
#' @keywords internal
extreme_coef_after_k_flips <- function(cont, # continuous var (residuals from model)
                                       bin, # binary var (treatment from model)
                                       k, # k number of (obs) flips
                                       treatment, # character (name treatment in model)
                                       dat, # data (dataframe)
                                       setting, # whether var is treatment or dependent
                                       direction = 'max',
                                       m # model
                                       )
  {

  stopifnot(length(cont) == length(bin)) # variables need to have same n

  # Define what "better" means, depending on whether we want the maximum or minimum extreme
  better <- if (direction == "max") `>` else `<`

  # Split continuous data
  c0 <- cont[bin == 0] # Continuous values for binary = 0
  c1 <- cont[bin == 1] # Continuous values for binary = 1

  # Record the number of treated and control observations
  n0 <- length(c0) # number of control obs
  n1 <- length(c1) # number of treated obs

  # Sort the vectors according to direction
  if(direction == 'max') { # max
    c0_sorted <- sort(c0, decreasing = TRUE)    # big positives in control
    c1_sorted <- sort(c1, decreasing = FALSE)   # big negatives in treated
  } else { # min
    c0_sorted <- sort(c0, decreasing = FALSE)   # big negatives in control
    c1_sorted <- sort(c1, decreasing = TRUE)    # big positives in treated
  }

  # Best scenarios
  ## The initial most extreme (i.e,. "best") beta is just what we observe in the data
  best_beta  <- coef(m)[treatment]
  betas <- NULL # create object for loop later
  ## Keeps track of how many obs have been flipped in each direction (treat/control):
  best_split <- c(0L, 0L)
  ### c (first element) = k flips from control to treat
  ### t (second element) = k flips from treat to control

  ## Loop = all possible ways to flip up to ⁠k observations between the two group
  ## to find the flips that produce the most extreme regression coefficient

  # NB: We are flipping based on cont, which is the *residuals* of the model for each obs

  for (c in 0:min(k, n0)) { ## cannot flip more than n0 control observations or more than k total
    t <- k - c # total flips k = t + c
    if (t > n1) next # Skip if trying to flip more treated observations than exist
    new_split <- c(c, t) # Store the current flip counts
    ## order direction flip
    if(direction == 'max') {
      which_to_flip = c(
        # Flip the control observations with the largest continuous values
        which(bin == 0)[order(c0, decreasing = TRUE)][seq_len(new_split[1])],
        # Flip the treated observations with the smallest continuous values
        which(bin == 1)[order(c1, decreasing = FALSE)][seq_len(new_split[2])]
      )
    } else {
      which_to_flip = c(
        # Flip the control observations with the smallest continuous values
        which(bin == 0)[order(c0, decreasing = FALSE)][seq_len(new_split[1])],
        # Flip the treated observations with the largest continuous values
        which(bin == 1)[order(c1, decreasing = TRUE)][seq_len(new_split[2])]
      )
    }
    ## Create a new binary variable and
    ## flip the selected observations by changing 0 → 1 or 1 → 0.
    dat$b_new <- bin
    dat$b_new[which_to_flip] <- 1 - bin[which_to_flip]

    ## Recalculate the regression coefficient
    ### If the binary variable is on the right-hand size (i.e., the treatment),
    ### update the regression accordingly
    if(setting == 'd') {
      tryCatch({
        model_updated <- suppressMessages(
          update(m, as.formula(paste0('. ~ . + b_new -', treatment)),
                                                 data = dat))
        beta_new <- coef(model_updated)['b_new']
        # Check if coefficient is NA or missing (indicates multicollinearity)
        if(is.na(beta_new) || is.null(beta_new)) next
      }, error = function(e) {
        # Skip this iteration if model update fails
        next
      })
    } else {
      ### Otherwise, the binary variable is the outcome,
      ### requiring a different approach to updating the regression
      tryCatch({
        model_updated <- suppressMessages(
          update(m,as.formula(paste0('b_new ~ .')),
                 data = dat))
        beta_new <- coef(model_updated)[treatment]
        # Check if coefficient is NA or missing (indicates multicollinearity)
        if(is.na(beta_new) || is.null(beta_new)) next
      }, error = function(e) {
        # Skip this iteration if model update fails
        next
      })
    }
    ## Store the new coefficient and the indices flipped and append to existing betas
    betas <- tibble(beta = beta_new,
                    which_to_flip = list(which_to_flip)) %>%
      as_tibble() %>%
      bind_rows(betas)

    ## Updated totals
    n1_new <- n1 + c - t
    n0_new <- n0 - c + t

    ## Skip illegal splits that zero out a group
    if (n1_new == 0 || n0_new == 0) next

    ## If the new beta is *more extreme* than the current "best" beta, update
    if (better(beta_new,best_beta)) {
      best_beta  <- beta_new
      best_split <- new_split
      best_flips <- which_to_flip
    }
    # Loop over different number of flips to get most extreme ("best") coef
  }

  # If the observed coefficient is the most extreme,
  # then give the second most extreme (edge case)

  ## if no flips were made that improved the coef, keep 0,0
  if(all(best_split == c(0,0))) {
    ## if no better coefficient was found, pick the second most extreme from the betas
    if(direction == 'min') {
      ##  pick the flip set with the minimum beta
      best_flips <- (betas %>%
                       slice(which.min(beta)) %>%
                       pull(which_to_flip))[[1]]
    } else {
      ##  pick the flip set with the maximum beta
      best_flips <- (betas %>%
                       slice(which.max(beta)) %>%
                       pull(which_to_flip))[[1]]
    }
    ## Flip the observations
    dat$b_new <- bin
    dat$b_new[best_flips] <- 1 - bin[best_flips]

    ## Recalculate the coefficient for this second-best flip set
    if(setting == 'd') {
      tryCatch({
        model_updated <- suppressMessages(update(m,as.formula(paste0('. ~ . + b_new -',treatment)),data = dat))
        best_beta <- coef(model_updated)['b_new']
        # If coefficient is still NA, use the original coefficient
        if(is.na(best_beta) || is.null(best_beta)) {
          best_beta <- coef(m)[treatment]
        }
      }, error = function(e) {
        # If model update fails, use original coefficient
        best_beta <- coef(m)[treatment]
      })
    } else {
      tryCatch({
        model_updated <- suppressMessages(update(m,as.formula(paste0('b_new ~ .')),data = dat))
        best_beta <- coef(model_updated)[treatment]
        # If coefficient is still NA, use the original coefficient
        if(is.na(best_beta) || is.null(best_beta)) {
          best_beta <- coef(m)[treatment]
        }
      }, error = function(e) {
        # If model update fails, use original coefficient
        best_beta <- coef(m)[treatment]
      })
    }
  }

  # Return value
  list(
    max_beta = best_beta,
    flips_from_control = best_split[1],
    flips_from_treated = best_split[2],
    ## identify the actual row indices to flip
    ## (this is the only part of this we really need)
    which_to_flip = best_flips
  )
}

#### Extreme Bound Prob Flip ####
# This function:
#   1) applies the extreme bounds function (extreme_coef_after_k_flips)
#   2) estimates the predicted probability of being selected for the extreme value
flip_beta_greedy <- function(dat,
                             treatment,
                             outcome,
                             m,
                             misclassRate,
                             binary,
                             direction = c("max", "min")) {

  # How many observations do we want to flip for a given mis-/re-classification rate?
  ## rounds percentage of misclassified rows to nearest integer
  n_flips <- round(misclassRate * nrow(dat))

  ## if the binary variable equals to the treatment character entry, label as independent
  if(binary == treatment) {
    setting <- 'd'
    ## if not, label as dependent variable
  } else {
    setting <- 'y'
  }

  # If the binary variable is on the right-hand-side (i.e., treatment),
  # get the extreme bounds accordingly
  if(binary == treatment) {
    res <- extreme_coef_after_k_flips(cont = dat$resid,
                                      bin = dat[[treatment]], ## treatment var
                                      k = n_flips,
                                      m,
                                      direction = direction,
                                      treatment = treatment,
                                      dat = dat,
                                      setting = setting)
  } else {
    res <- extreme_coef_after_k_flips(cont = dat$resid,
                                      bin = dat[[outcome]], ## dependent var
                                      k = n_flips,
                                      m,
                                      direction = direction,
                                      treatment = treatment,
                                      dat = dat,
                                      setting = setting)
  }

  # If the binary variable is on the right-hand-side (i.e., treatment),
  # estimate the new coefficient and calculate the predicted probabilities accordingly
  if(setting == 'd') {

    # Generate a new treatment variable whose observations will be reclassified to achieve the extreme bound
    dat$D_flipped <- dat[[treatment]]
    dat$D_flipped[res$which_to_flip] <- 1 - dat$D_flipped[res$which_to_flip]

    # Calculate the extreme bound with error handling
    tryCatch({
      model <- suppressMessages(update(m, as.formula(paste0('.~. + D_flipped -',treatment)),
                                       data = dat))
      extreme_coef <- coef(model)['D_flipped']

      # Check if coefficient is NA or missing (indicates multicollinearity)
      if(is.na(extreme_coef) || is.null(extreme_coef)) {
        # Return original coefficient if flipped variable causes multicollinearity
        extreme_coef <- coef(m)[treatment]
        # Use original variable for predictions
        dat$D_flipped <- dat[[treatment]]
      }
    }, error = function(e) {
      # If model update fails due to multicollinearity, use original coefficient
      extreme_coef <- coef(m)[treatment]
      # Reset to original variable
      dat$D_flipped <- dat[[treatment]]
    })

    # Model the probability of being selected to be flipped for the extreme bound
    preds <- ranger::ranger(as.formula(paste0('M ~ ',outcome,' + ',treatment)),
                            data = dat %>%
                              mutate(M = as.numeric(get(treatment) != D_flipped)))
    probs_rang <- preds$predictions

    # Calculate the correlation between the outcome and the reclassified variables (this is old)
    differential <- dat %>%
      mutate(M = get(treatment) != D_flipped) %>%
      group_by(get(treatment)) %>%
      summarise(diff = cor(get(outcome),M)) %>%
      pull(diff)

  } else {
    # Generate a new outcome variable whose observations will be reclassified to achieve the extreme bound
    dat$Y_flipped <- dat[[outcome]]
    dat$Y_flipped[res$which_to_flip] <- 1 - dat$Y_flipped[res$which_to_flip]

    # Calculate the extreme bound with error handling
    tryCatch({
      model <- suppressMessages(update(m,as.formula(paste0('Y_flipped~.')),
                                       data = dat))
      extreme_coef <- coef(model)[treatment]

      # Check if coefficient is NA or missing (indicates multicollinearity)
      if(is.na(extreme_coef) || is.null(extreme_coef)) {
        # Return original coefficient if flipped variable causes multicollinearity
        extreme_coef <- coef(m)[treatment]
        # Use original variable for predictions
        dat$Y_flipped <- dat[[outcome]]
      }
    }, error = function(e) {
      # If model update fails due to multicollinearity, use original coefficient
      extreme_coef <- coef(m)[treatment]
      # Reset to original variable
      dat$Y_flipped <- dat[[outcome]]
    })

    # Model the probability of being selected to be flipped for the extreme bound
    preds <- ranger::ranger(as.formula(paste0('M ~ ',outcome,' + ',treatment)),
                            data = dat %>%
                              mutate(M = as.numeric(get(outcome) != Y_flipped)))

    # Fudge to ensure we don't have perfect separation
    probs_rang <- preds$predictions

    # Calculate the correlation between the outcome and the reclassified variables (this is old)
    differential <- dat %>%
      mutate(M = get(outcome) != Y_flipped) %>%
      group_by(get(treatment)) %>%
      summarise(diff = cor(get(treatment),M)) %>%
      pull(diff)
  }



  return(list(extreme_coef = extreme_coef, # Extreme bound
              differential = differential, # Differential correlation (not needed)
              probs_rang = probs_rang, # Probability
              rawRes = res)) # Raw results (not needed)
}

#### Random Flip w/ Probs ####
# This function:
#   randomly samples observations to be flipped using the predicted probabilities
#   from the extreme bound
permutation_test_sens <- function(dat, # The data frame
                                  treatment,  # Which variable is the treatment
                                  outcome, # Which variable is the outcome
                                  M_extreme,
                                  nsims,
                                  prop_sense = .05,
                                  setting, # Whether binary variable is treatment or outcome
                                  m) # The model object
  {

  perm_sens_res <- NULL
  multicollinearity_count <- 0  # Track multicollinearity issues
  total_simulations <- 0  # Track total simulations

  for(prop in prop_sense) {
    alts <- setdiff(1:nrow(dat),M_extreme)
    size <- ceiling(prop*length(M_extreme))

    coefs <- NULL
    for(i in 1:nsims) {
      total_simulations <- total_simulations + 1
      perm_flips <- sample(alts,size = size,replace = F)
      extr_flips <- sample(M_extreme,size = size,replace = F)

      M_perm <- M_extreme
      M_perm[which(M_perm %in% extr_flips)] <- perm_flips

      # Initialize perm_beta
      perm_beta <- NA

      if(setting == 'd') {
        dat$D_flipped <- dat[[treatment]]
        dat$D_flipped[M_perm] <- 1 - dat$D_flipped[M_perm]

        perm_beta <- tryCatch({
          # Suppress warnings about multicollinearity
          model_updated <- suppressMessages(update(m, as.formula(paste0('. ~ . + D_flipped -', treatment)),
                                                   data = dat))
          beta_result <- coef(model_updated)['D_flipped']
          # Check if coefficient is NA or missing (indicates multicollinearity)
          if(is.na(beta_result) || is.null(beta_result)) {
            multicollinearity_count <- multicollinearity_count + 1
            NA  # Skip this simulation
          } else {
            beta_result
          }
        }, error = function(e) {
          multicollinearity_count <- multicollinearity_count + 1
          NA  # Skip this simulation if model fails
        })
      } else {
        dat$Y_flipped <- dat[[outcome]]
        dat$Y_flipped[M_perm] <- 1 - dat$Y_flipped[M_perm]

        perm_beta <- tryCatch({
          # Suppress warnings about multicollinearity
          model_updated <- suppressMessages(update(m,as.formula(paste0('Y_flipped ~ .')),data = dat))
          beta_result <- coef(model_updated)[treatment]
          # Check if coefficient is NA or missing (indicates multicollinearity)
          if(is.na(beta_result) || is.null(beta_result)) {
            multicollinearity_count <- multicollinearity_count + 1
            NA  # Skip this simulation
          } else {
            beta_result
          }
        }, error = function(e) {
          multicollinearity_count <- multicollinearity_count + 1
          NA  # Skip this simulation if model fails
        })
      }
      coefs <- c(coefs,perm_beta)
    }
    perm_sens_res <- data.frame(coefs = coefs,prop = prop) %>%
      bind_rows(perm_sens_res) %>%
      as_tibble()
  }

  # Report multicollinearity summary
  if(multicollinearity_count > 0) {
    cat(sprintf('\nMulticollinearity in %d of %d simulations\n',
                multicollinearity_count, total_simulations))
  }

  return(perm_sens_res)
}

### Random Flip w/ Pred Probs ###
# This function:
#   randomly samples observations to be flipped using the predicted probabilities
#   from the extreme bound
permutation_test_prob <- function(dat, # The data frame
                                  probVect, # The probabilities
                                  ## (from the ranger function inside flip_beta_greedy())
                                  R_vect, # The vector of reclassification rates
                                  nsims, # The number of simulations
                                  treatment, # Which variable is the treatment
                                  m, # The model object
                                  outcome, # Which variable is the outcome
                                  tol, # How close to the target reclassification rate are we willing to accept?
                                  # (Due to the random nature of the new M vector,
                                  # it might not be exactly the right proportion)
                                  setting) { # Is the binary variable the treatment or the outcome?

  ests <- NULL
  Mask <- list()
  for(zz in 1:length(probVect)) {
    MNew <- rbinom(n = nrow(dat), size = 1, prob = as.numeric(probVect[[zz]])) # Get a new indicator for which observations will be flipped
    while(abs(R_vect[zz] - mean(MNew)) > tol) { # Re-run this until we are within tol of the target reclassification rate
      MNew <- rbinom(n = nrow(dat), size = 1, prob = as.numeric(probVect[[zz]]))
    }
    if(setting == 'd') {
      # Generate a new treatment variable whose observations will be reclassified to approach the extreme bound
      dat$D_perm <- dat[[treatment]]
      dat$D_perm[which(MNew == 1)] <- 1 - dat[[treatment]][which(MNew == 1)]

      # Calculate the differential bound and save the result
      tryCatch({
        model <- suppressMessages(update(m,as.formula(paste0('.~. + D_perm -',treatment)),data = dat))
        est <- coef(model)['D_perm']
        # Check if coefficient is NA or missing (indicates multicollinearity)
        if(is.na(est) || is.null(est)) {
          est <- coef(m)[treatment]  # Use original coefficient
        }
      }, error = function(e) {
        est <- coef(m)[treatment]  # Use original coefficient if model fails
      })
    } else {
      # Generate a new outcome variable whose observations will be reclassified to approach the extreme bound
      dat$Y_perm <- dat[[outcome]]
      dat$Y_perm[which(MNew == 1)] <- 1 - dat[[outcome]][which(MNew == 1)]

      # Calculate the differential bound and save the result
      tryCatch({
        model <- suppressMessages(update(m,as.formula(paste0('Y_perm~.')),data = dat))
        est <- coef(model)[treatment]
        # Check if coefficient is NA or missing (indicates multicollinearity)
        if(is.na(est) || is.null(est)) {
          est <- coef(m)[treatment]  # Use original coefficient
        }
      }, error = function(e) {
        est <- coef(m)[treatment]  # Use original coefficient if model fails
      })
    }
    ests <- c(ests,est)
    Mask[[zz]] <- MNew
  }
  return(list(ests = ests,
              Mask = Mask))
}

### Naive bounds (random w/o probs) ###
# This function: randomly samples observations to be flipped
basic_flips <- function(binary, # Which is the binary variable (character)
                        misclassRate, # What proportion do we want to reclassify?
                        m, # This is the original model
                        setting, # Is the binary variable the outcome or the treatment?
                        treatment, # Which is the treatment variable (character)
                        outcome, # Which is the outcome variable (character)
                        dat) {
  n <- nrow(dat)
  # Loop over the desired reclassification rates
  ests <- NULL
  for(mcr in misclassRate) {
    nSamp <- round(mcr*n)
    # Should respect the underlying distribution though
    ## How many observations will we switch from 1 to 0?
    nSamp1s <- round(length(which(dat[[binary]] == 1)) * mcr)
    ## # How many observations will we switch from 0 to 1?
    nSamp0s <- round(length(which(dat[[binary]] == 0)) * mcr)

    # Naively select these observations at random
    inds1 <- sample(which(dat[[binary]] == 1),
                    size = nSamp1s,
                    replace = F)
    inds0 <- sample(which(dat[[binary]] == 0),
                    size = nSamp0s,
                    replace = F)
    inds <- c(inds1,inds0)

    # Generate a new binary variable whose observations will be reclassified at random
    dat$b_simp <- dat[[binary]]
    dat$b_simp[inds] <- 1 - dat$b_simp[inds]
    if(setting == 'd') {
      tryCatch({
        model <- suppressMessages(update(m,as.formula(paste0('.~. + b_simp -',treatment)),
                                         data = dat))
        est <- coef(model)['b_simp']
        # Check if coefficient is NA or missing (indicates multicollinearity)
        if(is.na(est) || is.null(est)) {
          est <- coef(m)[treatment]  # Use original coefficient
        }
        ests <- c(ests, est)
      }, error = function(e) {
        # Use original coefficient if model fails
        ests <- c(ests, coef(m)[treatment])
      })
    } else {
      tryCatch({
        model <- suppressMessages(update(m,as.formula(paste0('b_simp~.')),
                                         data = dat))
        est <- coef(model)[treatment]
        # Check if coefficient is NA or missing (indicates multicollinearity)
        if(is.na(est) || is.null(est)) {
          est <- coef(m)[treatment]  # Use original coefficient
        }
        ests <- c(ests, est)
      }, error = function(e) {
        # Use original coefficient if model fails
        ests <- c(ests, coef(m)[treatment])
      })
    }
  }
  return(ests)
}

##### Contour Plots #####

### Generate Mask Variable ###
# This function generates a binary "mask" variable M that achieves the following:
#   1) It has some proportion of 1s given by the user (p)
#   2) It is correlated with a continuous variable Y, given by the user (r_target)
# It achieves this by first choosing the observations that maximize the correlation
#   between the mask vector M and the continuous variable Y. It then randomly flips
#   some of these chosen observations until the correlation declines to within a tolerance
#   parameter (provided by the user) of the desired correlation.
corrupt_M_to_target <- function(Y,
                                p = 0.1,
                                r_target = 0.2,
                                tol = 0.01,
                                max_iter = 1000) {
  n <- length(Y)
  stopifnot(p > 0 && p < 1)

  # Step 1: generate max correlation
  if (r_target >= 0) {
    cutoff <- quantile(Y, probs = 1 - p)
    M <- as.numeric(Y >= cutoff)
    if(cor(M,Y) < r_target) {
      r_target <- cor(M,Y)
      warning('The desired correlation is greater than the maximum possible given the desired proportion.')
    }
  } else {
    cutoff <- quantile(Y, probs = p)
    M <- as.numeric(Y < cutoff)
    if(cor(M,Y) > r_target) {
      r_target <- cor(M,Y)
      warning('The desired correlation is less than the minimum possible given the desired proportion.')
    }
  }

  cor_history <- numeric()
  iter <- 0

  # Greedy swapping loop
  while ((r_target >= 0 && cor(M, Y) > r_target + tol) ||
         (r_target <  0 && cor(M, Y) < r_target - tol)) {

    iter <- iter + 1
    ones <- which(M == 1)
    zeros <- which(M == 0)

    # Randomly flip observations until the correlation decays to the desired value
    flip_1 <- sample(ones, 1)
    flip_0 <- sample(zeros, 1)

    M_new <- M
    M_new[flip_1] <- 0
    M_new[flip_0] <- 1

    # Accept swap only if it brings us closer to target
    old_diff <- abs(cor(M, Y) - r_target)
    new_diff <- abs(cor(M_new, Y) - r_target)

    if (new_diff < old_diff) {
      M <- M_new
    }

    if (iter %% 1000 == 0) cor_history <- c(cor_history, cor(M, Y))
    if (iter >= max_iter) break
  }

  return(list(
    M = M,
    final_cor = cor(M, Y),
    iterations = iter,
    final_prevalence = mean(M),
    trace = cor_history
  ))
}

### Data-generating-process w/ multiple controls ###
me_Complex <- function(N = 100,            # number of observations
                       propD = .5,         # proportion D = 1
                       a = 0,              # alpha
                       bD = 1,             # beta_D
                       bX1 = 0,            # beta_X1
                       bX2 = 0,            # beta_X2
                       prop_ME = .1,       # proportion measured with error
                       rho_ME_Y = NULL,    # correlation between measurement error and Y
                       rho_ME_D = NULL,    # correlation between measurement error and D
                       rho_X1_D = 0,       # correlation between X1 and D
                       rho_X2_D = 0,       # correlation between X2 and D
                       rho_X2_X1 = 0,      # correlation between X2 and X1
                       seed = 123,         # seed for consistency
                       noise               # disturbance term
) {
  set.seed(seed)

  # Set D and D_star to be identical to start, subject to the number of observations
  #   desired (N) and the proportion of observations for which D_star = 1 (propD).
  D    <- D_star <- rbinom(n = N,size = 1,prob = propD)

  # Generating X1 and X2 and D with a variance-covariance matrix
  # Note: mvtnorm functions are imported via package imports

  # Define the overall mean vector and covariance matrix
  # (Assuming a multivariate normal distribution with D, X1, X2)
  overall_mean <- c(propD, 10, -2) # Mean for D, X1, X2
  overall_sigma <- matrix(c(
    1, rho_X1_D, rho_X2_D,
    rho_X1_D, 1, rho_X2_X1,
    rho_X2_D, rho_X2_X1, 1
  ), nrow = 3) # Covariance matrix for D, X1, X2

  # Separate indices for fixed and target variables
  target_ind <- c(2, 3) # Indices for X1 and X2
  given_ind <- 1          # Index for D

  # Define the conditional value of D (e.g., D = 1)
  D_given <- 1

  # Calculate the conditional mean and variance-covariance matrix
  # Note: condMVNorm functions are imported via package imports
  conditional_params_D1 <- condMVNorm::rcmvnorm(n = sum(D_star),
                                    mean = overall_mean,
                                    sigma = overall_sigma,
                                    dependent.ind = target_ind,
                                    given.ind = given_ind,
                                    X.given = D_given
  )

  D_given <- 0
  conditional_params_D0 <- condMVNorm::rcmvnorm(n = sum(D_star == 0),
                                    mean = overall_mean,
                                    sigma = overall_sigma,
                                    dependent.ind = target_ind,
                                    given.ind = given_ind,
                                    X.given = D_given
  )

  # Build the right hand side of the dataset
  dat_sim <- data.frame(conditional_params_D0) %>%
    mutate(D_star = 0) %>%
    bind_rows(data.frame(conditional_params_D1) %>%
                mutate(D_star = 1)) %>%
    as_tibble()


  # Generate a continuous Y as a function of the user-defined parameters
  dat_sim$Y <- a + bD * dat_sim$D_star + bX1 * dat_sim$X1 + bX2*dat_sim$X2 +
    rnorm(n = N,mean = 0,sd = noise)


  # Generate the misclassified observations subject to some amount of differential
  #   measurement error. Currently this cannot be a function of more than one variable.
  if(!is.null(rho_ME_Y)) {
    var <- 'Y'
    r_target <- rho_ME_Y
  } else if(!is.null(rho_ME_D)) {
    var <- 'D_star'
    r_target <- rho_ME_D
  } else if(!is.null(rho_ME_X1)) {
    var <- 'X1'
    r_target <- rho_ME_X1
  } else if(!is.null(rho_ME_X2)) {
    var <- 'X2'
    r_target <- rho_ME_X2
  } else {
    var <- 'Y'
    r_target <- rho_ME_Y
  }

  # Generate a mask vector M that achieve the desired proportion to be misclassified
  #   and the desired differential measurement error.
  res <- corrupt_M_to_target(Y = dat_sim[[var]],
                             p = prop_ME,
                             r_target = r_target)
  inds <- which(res$M == 1)

  # Misclassify D
  dat_sim$D <- dat_sim$D_star
  dat_sim$D[inds] <- 1-dat_sim$D_star[inds]

  dat_sim <- dat_sim %>%
    mutate(propD = propD,
           prop_ME = prop_ME,
           rho_X1_D = rho_X1_D,
           rho_X2_D = rho_X2_D,
           rho_X2_X1 = rho_X2_X1,
           a = a,
           bD = bD,
           bX1 = bX1,
           bX2 = bX2)

  return(dat_sim)
}

##### Compute lambda ####
compute_lambda <- function(propD,
                           prop_ME,
                           rho_M_Y,
                           rho_M_D,
                           rho_X1_D,
                           rho_X2_D,
                           rho_X2_X1,
                           seed = seed,
                           noiseObs = 0.1, # default
                           method = 'Bias') {

  if(method == 'Bias') {
    simDat <- me_Complex(N = 1000,
                         propD = propD,
                         prop_ME = prop_ME,
                         rho_ME_Y = rho_M_Y,
                         rho_ME_D = rho_M_D,
                         rho_X1_D = rho_X1_D,
                         rho_X2_D = rho_X2_D,
                         rho_X2_X1 = rho_X2_X1,
                         seed = seed,
                         noise = 0.01)
    model <- lm(Y ~ D + X1 + X2,simDat)
    return(coef(model)[2])
  } else if(method == 'Coverage') {
    simDat <- me_Complex(N = 1000,
                         propD = propD,
                         prop_ME = prop_ME,
                         rho_ME_Y = rho_M_Y,
                         rho_ME_D = rho_M_D,
                         rho_X1_D = rho_X1_D,
                         rho_X2_D = rho_X2_D,
                         rho_X2_X1 = rho_X2_X1,
                         seed = seed,
                         noise = noiseObs)
    model <- lm(Y ~ D + X1 + X2,simDat)
    confint_ci <- confint(model)[2, ]
    return(1 >= confint_ci[1] && 1 <= confint_ci[2])
  }
}

compute_lambda_vect <- Vectorize(compute_lambda)

##### Gradient fill colors ####
# For `contour_plots`, this function calculates the number of filled areas in the contour
# plot and assigns color codes to indicate a clear break to/around 0
generate_palette <- function(sim_lambda) {
  ## Check range of sim_lambda
  range_sim_lambda <- range(sim_lambda, na.rm = TRUE)
  ## Breakdown pos/zero/neg
  pos_breaks <- if (any(range_sim_lambda >= 0.1)) {
    seq(0.1, floor(range_sim_lambda[2]*10)/10, by = 0.1)
  } else {
    numeric(0)
  }
  neg_breaks <- if (range_sim_lambda[1] <= -0.1) {
    seq(ceiling(range_sim_lambda[1]*10)/10, -0.1, by = 0.1)
  } else {
    numeric(0)
  }
  zero_break <- if (range_sim_lambda[1] <= 0 && range_sim_lambda[2] >= 0)
    0 else numeric(0) ## Zero break if zero is in range
  ## Depending on range, get number of colors
  neg_colors <- colorRampPalette(c("#EA0A23", "#1A18FF"))(length(neg_breaks)+1)
  pos_colors <- colorRampPalette(c("#FFB662", "#FEF281", "#3DA100"))(length(pos_breaks)+1)
  fill_colors <- c(neg_colors, pos_colors)
  return(fill_colors)
}
