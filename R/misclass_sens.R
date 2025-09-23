#' Sensitivity analysis for measurement error in binary regression coefficients
#'
#' This function performs sensitivity analysis by simulating misclassification
#' in a binary variable and estimating bounds on regression coefficients.
#'
#' @param dat Data frame containing the data.
#' @param outcome Character. Name of the outcome variable.
#' @param treatment Character. Name of the treatment variable.
#' @param binary Character. Name of the binary variable subject to misclassification.
#' @param nsims Integer. Number of simulations to run.
#' @param tol Numeric. Tolerance for reclassification rate.
#' @param max_iter Integer. Maximum iterations for joint restrictions.
#' @param m Model object. The fitted regression model.
#' @param R_vect Numeric vector. Proportions of observations to reclassify.
#'
#' @return A list containing bounds and simulation results.
#'
#' @importFrom ranger ranger
#' @import dplyr
#'
#' @examples
#' # Example usage with your data and model
#' results <- misclass_sens(dat = your_data, outcome = "Y", treatment = "D",
#'                          binary = "D", nsims = 100, tol = 0.01, m = your_model)
#' @export

# Let's pretend this is the real data supplied by the user.
#   We ask them to define how much measurement error they have,
#   and the magnitude of the correlation between ME and one of the variables.
misclass_sens <- function(dat,               # Data frame: either data.frame or tibble
                          outcome,           # Name of outcome: character
                          treatment,         # Name of treatment: character
                          binary,            # Name of variable for reclassification: character
                          nsims = 30,        # Number of simulations: integer
                          tol = 0.01,        # Tolerance for corrupt_M_to_target(): double
                          max_iter = 1000,   # Maximum number of iterations to solve joint restriction
                          m,                 # Need the model for flexibility
                          R_vect = seq(.01,.5,by = .01)
                          # How many values of the proportion reclassified should we check?
) {

  # Validate all inputs
  validate_misclass_inputs(dat = dat, outcome = outcome, treatment = treatment,
                          binary = binary, m = m, nsims = nsims, tol = tol,
                          R_vect = R_vect)

  if(inherits(m, 'lm')) {
    dat <- model.frame(formula(m), data = dat, na.action = na.omit)
  } else if(inherits(m, 'fixest')) {
    # For fixest models, keep only the observations that were used in the model
    # obsRemoved contains the indices of observations that were actually used
    if(!is.null(m$obs_selection$obsRemoved) && length(m$obs_selection$obsRemoved) > 0) {
      dat <- dat %>% slice(m$obs_selection$obsRemoved)
    }
    # Alternative: if obsRemoved actually contains removed obs, use:
    # dat <- dat %>% slice(-m$obs_selection$obsRemoved)
  } else {
    dat <- model.frame(formula(m), data = dat, na.action = na.omit)
    warning('Not yet tested for ', class(m)[1])
  }

  # Get residuals
  ## (this is what we use instead of the raw continuous value if there are controls)
  if(treatment == binary) {
    # Calculate residualized Y if D = 0
    resid0 <- resid(update(m, as.formula(paste0(outcome,'~. -',treatment)),
                           data = dat %>%
                             filter(get(treatment) == 0)))
    # Calculate residualized Y if D = 1
    resid1 <- resid(update(m, as.formula(paste0(outcome,'~. -',treatment)),
                           data = dat %>%
                             filter(get(treatment) == 1)))
    # Add the residuals to the data frame
    dat <- dat %>%
      filter(get(treatment) == 0) %>%
      mutate(resid = resid0) %>%
      bind_rows(dat %>%
                  filter(get(treatment) == 1) %>%
                  mutate(resid = resid1))
  } else {
    # Calculate residualized D if Y = 0
    resid0 <- resid(update(m, as.formula(paste0(treatment,'~. -',treatment)),
                           data = dat %>%
                             filter(get(outcome) == 0)))
    # Calculate residualized D if Y = 1
    resid1 <- resid(update(m, as.formula(paste0(treatment,'~. -',treatment)),
                           data = dat %>%
                             filter(get(outcome) == 1)))

    # Add the residuals to the data frame
    dat <- dat %>%
      filter(get(outcome) == 0) %>%
      mutate(resid = resid0) %>%
      bind_rows(dat %>%
                  filter(get(outcome) == 1) %>%
                  mutate(resid = resid1))
  }

  # Calculate the maximum extreme bound
  cat('Calculating extreme maximum bound\n')
  extMax <- prog_bar_sapply(R_vect, function(x) # Iterate over different proportions (R_vect)
    flip_beta_greedy(dat = dat, # The data frame
                     treatment = treatment, # Which variable is the treatment (character)
                     outcome = outcome, # Which variable is the outcome (character)
                     binary = binary, # Which variable is the binary variable (character)
                     misclassRate = x, # What proportion do we want to reclassify
                     m = m, # The original model
                     direction = 'max')) # Which direction is the bound

  # Calculate the minimum extreme bound
  cat('Calculating extreme minimum bound\n')
  extMin <- prog_bar_sapply(R_vect, function(x)
    flip_beta_greedy(dat = dat,
                     treatment = treatment,
                     outcome = outcome,
                     binary = binary,
                     misclassRate = x,
                     m = m,
                     direction = 'min'))

  # Inference: First by observed differential reclassification rates
  if(binary == treatment) {
    setting <- 'd'
    diffVar <- outcome
  } else {
    setting <- 'y'
    diffVar <- treatment
  }
  cat('Calculating differential 95% bound for maximum\n')

  # Filter out NULL results from extreme bounds that failed due to multicollinearity
  extMax_valid_indices <- which(!sapply(extMax[4,], is.null) &
                                  !is.na(as.numeric(extMax[1,])))
  extMin_valid_indices <- which(!sapply(extMin[4,], is.null) &
                                  !is.na(as.numeric(extMin[1,])))

  permMax <- prog_bar_lapply(extMax_valid_indices, function(x)
    permutation_test_sens(dat = dat,
                          treatment = treatment,
                          outcome = outcome,
                          M_extreme = extMax[4,][[x]]$which_to_flip,
                          nsims = nsims,
                          prop_sense = .05, # default value
                          setting = setting,
                          m = m))

  permMaxProb <- prog_bar_sapply(1:nsims, function(x)
    permutation_test_prob(dat = dat, # The data frame
                          probVect = extMax[3,],
                          # This is the vector of probabilities
                          R_vect = R_vect,
                          # This is the vector of reclassification rates
                          treatment = treatment,
                          # Which variable is the treatment (character)
                          outcome = outcome,
                          # Which variable is the outcome (character)
                          tol = tol,
                          m = m,
                          # How much deviation from the desired
                          # reclassification rate are we willing to accept?
                          nsims = nsims,
                          setting = setting)) # Whether binary variable is treatment or  outcome

  # permMax[2,][[1]]
  # apply(do.call(rbind,permMax[1,]),2,mean,na.rm=T)

  cat('Calculating differential 5% bound for minimum\n')
  # permMin <- prog_bar_sapply(1:nsims,function(x) permutation_test(dat = dat,
  #                                                          probVect = extMin[4,],
  #                                                          R_vect = R_vect,
  #                                                          treatment = treatment,
  #                                                          outcome = outcome,
  #                                                          tol = tol,
  #                                                          setting = setting))
  permMin <- prog_bar_lapply(extMin_valid_indices, function(x)
    permutation_test_sens(dat = dat,
                          treatment = treatment,
                          outcome = outcome,
                          M_extreme = extMin[4,][[x]]$which_to_flip,
                          nsims = nsims,
                          prop_sense = .05, # default value
                          setting = setting,
                          m = m))

  permMinProb <- prog_bar_sapply(1:nsims, function(x)
    permutation_test_prob(dat = dat,
                          # The data frame
                          probVect = extMin[3,],
                          # This is the vector of probabilities
                          R_vect = R_vect,
                          # This is the vector of reclassification rates
                          treatment = treatment,
                          # Which variable is the treatment (character)
                          outcome = outcome,
                          # Which variable is the outcome (character)
                          tol = tol,
                          m = m,
                          # How much deviation from desired reclassification rate we can accept?
                          nsims = nsims,
                          setting = setting))

  # Save the extreme bounds
  boundsRes <- data.frame(beta_max = as.numeric(extMax[1,]),
                          beta_min = as.numeric(extMin[1,])) %>%
    as_tibble()

  # Add in the true coefficient and the vector of reclassification rates
  boundsRes <- boundsRes %>%
    mutate(obsCoef = coef(m)[treatment],
           misclassRate = R_vect)

  # Identify valid indices (non-NA coefficients)
  valid_indices <- which(complete.cases(boundsRes) &
                           !is.na(boundsRes$beta_max) &
                           !is.na(boundsRes$beta_min))

  # Filter bounds results
  boundsRes <- boundsRes %>%
    filter(complete.cases(.), !is.na(beta_max), !is.na(beta_min))

  cat('Calculating naive bounds\n')
  combSim <- prog_bar_sapply(1:nsims, function(x)
    basic_flips(binary = binary,
                misclassRate = boundsRes$misclassRate,
                setting = setting,
                treatment = treatment,
                outcome = outcome,
                m = m,
                dat = dat))

  perm_sens_max_sum <- map(permMax, ~.x %>%
                             group_by(prop) %>%
                             summarise(coef = mean(coefs, na.rm=T))) %>%
    bind_rows(.id = 'source') %>%
    spread(prop,coef,sep = '_') %>%
    mutate(source = as.numeric(source)) %>%
    arrange(source) %>%
    select(-source)

  perm_sens_min_sum <- map(permMin, ~.x %>%
                             group_by(prop) %>%
                             summarise(coef = mean(coefs, na.rm=T))) %>%
    bind_rows(.id = 'source') %>%
    spread(prop,coef,sep = '_') %>%
    mutate(source = as.numeric(source)) %>%
    arrange(source) %>%
    select(-source)

  # Pull everything together
  boundsRes <- boundsRes %>%
    # There are two sets of differential bounds: those associated with the extreme maximum
    bind_cols(perm_sens_max_sum %>%
                rename_all(function(x) paste0(x,'_max'))) %>%
    bind_cols(perm_sens_min_sum %>%
                rename_all(function(x) paste0(x,'_min'))) %>%
    bind_cols(tryCatch({
      perm_prob_max_data <- do.call(rbind, permMaxProb[1,])
      if(is.null(perm_prob_max_data) || nrow(perm_prob_max_data) == 0) {
        data.frame(perm_prob_max = rep(NA, nrow(boundsRes)))
      } else {
        data.frame(perm_prob_max = apply(perm_prob_max_data[valid_indices, , drop=FALSE], 2,
                                         function(x) mean(x,na.rm=T)))
      }
    }, error = function(e) {
      data.frame(perm_prob_max = rep(NA, nrow(boundsRes)))
    }) %>% as_tibble()) %>%
    bind_cols(tryCatch({
      perm_prob_min_data <- do.call(rbind, permMinProb[1,])
      if(is.null(perm_prob_min_data) || nrow(perm_prob_min_data) == 0) {
        data.frame(perm_prob_min = rep(NA, nrow(boundsRes)))
      } else {
        data.frame(perm_prob_min = apply(perm_prob_min_data[valid_indices, , drop=FALSE], 2,
                                         function(x) mean(x,na.rm=T)))
      }
    }, error = function(e) {
      data.frame(perm_prob_min = rep(NA, nrow(boundsRes)))
    }) %>% as_tibble()) %>%
    as_tibble() %>%
    # The naive bounds is just a single distribution from which we grab the upper and lower 95% interval values
    bind_cols(tryCatch({
      combSim_filtered <- combSim[valid_indices, , drop=FALSE]
      if(is.null(combSim_filtered) || nrow(combSim_filtered) == 0 || ncol(combSim_filtered) == 0) {
        data.frame(sim_max = rep(NA, nrow(boundsRes)), sim_min = rep(NA, nrow(boundsRes)))
      } else {
        data.frame(t(apply(combSim_filtered, 1, function(x) quantile(x, probs = c(.025,.975), na.rm=T)))) %>%
          rename(sim_max = X2.5., sim_min = X97.5.)
      }
    }, error = function(e) {
      data.frame(sim_max = rep(NA, nrow(boundsRes)), sim_min = rep(NA, nrow(boundsRes)))
    }) %>% as_tibble()) %>%
    as_tibble()

  # Spot checking
  # extMax[4,1]
  # which(permMax[2,][[2]][[1]] == 1)
  # permMax[1,][[2]][1]
  return(list(bounds = boundsRes,
              rawExtremeMax = extMax,
              rawExtremeMin = extMin,
              rawPermSensMax = permMax,
              rawPermSensMin = permMin,
              rawPermProbMax = permMaxProb,
              rawPermProbMin = permMinProb,
              rawNaive = combSim))
}
