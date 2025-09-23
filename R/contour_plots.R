#' Generate Contour Plot Data for Misclassification Analysis
#'
#' This function generates the data needed to create contour plots that visualize the
#' impact of misclassification rates and proportions of treated units on the bias of
#' treatment effect estimates.
#'
#' @param dat A dataframe containing the data.
#' @param m A fitted model object (e.g., from `lm`) used to extract variables.
#' @param treatment A string specifying the name of the treatment variable in `dat`.
#' @param outcome A string specifying the name of the outcome variable in `dat`.
#' @param binary A string specifying the name of the observed binary treatment variable
#' subject to misclassification (could be `treatment` or `outcome`).
#' @param ground_truth Optional string specifying the name of the ground truth treatment
#' variable in `dat`.
#' If provided, it is used to calculate observed misclassification rates and correlations.
#' @param confMat Optional 2x2 confusion matrix for the binary treatment variable.
#' Used to estimate misclassification rate if `ground_truth` is not provided.
#' @param R_obs Optional numeric value specifying the user's best guess of the
#' misclassification rate if neither `ground_truth` nor `confMat` is provided.
#' @param rho_M_Y Optional numeric correlation between misclassification indicator and
#' outcome variable. Required if `ground_truth` is not provided.
#' @param rho_M_D Optional numeric correlation between misclassification indicator and
#' treatment variable. Required if `ground_truth` is not provided.
#' @param nsims Number of simulations to run for estimating lambda values. Default is 100.
#' @param R_vect Numeric vector of misclassification rates to evaluate.
#' Default is `seq(0.025, 0.5, by = 0.025)`.
#' @param pi_vect Numeric vector of proportions of treated units to evaluate.
#' Default is `seq(0.05, 0.95, by = 0.05)`.
#'
#' @return A list with components:
#' \describe{
#'   \item{toplot}{A data frame containing the grid of misclassification rates, proportions,
#'   and simulated lambda values.}
#'   \item{R_obs}{Observed or user-supplied misclassification rate.}
#'   \item{pi_obs}{Observed or user-supplied proportion of treated units.}
#' }
#'
#' @details
#' The function calculates observed misclassification rates and correlations if ground
#' truth data is available. Otherwise, it relies on user-supplied values or confusion
#' matrices. It then simulates data under various misclassification and treatment
#' proportion scenarios to estimate the bias or coverage of treatment effect estimates.
#'
#' @examples
#' \dontrun{
#' contour_data <- contour_plot_data(dat = mydata,
#'                                  m = lm(Y ~ D + X1 + X2, data = mydata),
#'                                  treatment = "D",
#'                                  outcome = "Y",
#'                                  binary = "D_obs",
#'                                  ground_truth = "D_true",
#'                                  nsims = 50)
#' }
#'
#' @export
contour_plot_data <- function(dat,
                              m,
                              treatment,
                              outcome,
                              binary,
                              ground_truth = NULL,
                              confMat = NULL,
                              R_obs = NULL,
                              rho_M_Y = NULL,
                              rho_M_D = NULL,
                              nsims = 100,
                              R_vect = seq(0.025,.5,by = .025),
                              pi_vect = seq(.05,.95,by = .05)) {

  # Validate inputs
  validate_contour_inputs(dat = dat, m = m, treatment = treatment,
                         outcome = outcome, binary = binary,
                         ground_truth = ground_truth, rho_M_Y = rho_M_Y,
                         rho_M_D = rho_M_D, nsims = nsims)

  # Extract parameters of interest from the data / model
  # Skew:
  # If the ground truth variable is provided, use this
  if(!is.null(ground_truth)) {
    pi_obs <- mean(dat[[ground_truth]], na.rm=T)
  } else { # otherwise, use the observed binary variable
    pi_obs <- mean(dat[[binary]],na.rm=T)
  }

  # Misclassification Rate:
  # If the ground truth variable is provided, use this
  if(!is.null(ground_truth)) {
    R_obs <- mean(dat[[binary]] != dat[[ground_truth]],na.rm=T)
  } else if(!is.null(confMat)) { # otherwise, if the confusion matrix is provided, use this
    R_obs <- (confMat[1,2] + confMat[2,1]) / sum(confMat)
  } else if(is.null(R_obs)) { # otherwise the user needs to provide their best guess
    R_obs <- 0
    warning("Either a ground truth variable needs to be provided (ground_truth), or a confusion matrix (confMat), or the user's best guess about the proportion of data that are incorrectly labeled (R_obs). If not, it is assumed that there is no misclassification.")
  }

  # Differential measurement error:
  # If the ground truth variable is provided, use this
  if(!is.null(ground_truth)) {
    dat$M <- dat[[binary]] != dat[[ground_truth]]
    rho_M_Y <- cor(dat$M, dat[[outcome]],
                   use = 'pairwise.complete.obs')
    rho_M_D <- cor(dat$M, dat[[treatment]],
                   use = 'pairwise.complete.obs')
  } else if(is.null(rho_M_Y) & is.null(rho_M_D)) {
    # otherwise the user needs to specify either rho_M_Y or rho_M_D
    rho_M_Y = 0
    rho_M_D = 0
    warning("Either a ground truth variable needs to be provided (ground truth), or a user-provided correlation between the misclassified observations and Y (rho_M_Y) or D (rho_M_D) needs to be provided.")
  }

  # Correlated controls
  # Extract control variables from model using standardized helper
  ctrls <- extract_controls(m, treatment)

  # Handle different numbers of control variables
  if(length(ctrls) != 0) {
    sigma <- var(dat[,c(treatment, ctrls)])
    rho_X1_D <- as.matrix(sigma)[treatment, ctrls[1]]

    # Handle cases with different numbers of control variables
    if(length(ctrls) >= 2) {
      # Use the two controls most correlated with treatment if more than 2
      if(length(ctrls) > 2) {
        # Calculate correlations with treatment for all controls
        corr_with_treatment <- abs(as.matrix(sigma)[treatment, ctrls])
        # Sort by correlation magnitude and take top 2
        top_ctrls <- names(sort(corr_with_treatment, decreasing = TRUE))[1:2]
        rho_X1_D <- as.matrix(sigma)[treatment, top_ctrls[1]]
        rho_X2_D <- as.matrix(sigma)[treatment, top_ctrls[2]]
        rho_X2_X1 <- as.matrix(sigma)[top_ctrls[1], top_ctrls[2]]
      } else {
        # Exactly 2 controls
        rho_X2_D <- as.matrix(sigma)[treatment, ctrls[2]]
        rho_X2_X1 <- as.matrix(sigma)[ctrls[1], ctrls[2]]
      }
    } else {
      # Only 1 control: set second control correlations to 0
      rho_X2_D <- 0
      rho_X2_X1 <- 0
    }
  } else {
    # No controls: set all correlations to 0
    rho_X1_D <- 0
    rho_X2_D <- 0
    rho_X2_X1 <- 0
  }

  # This is the vector of possible combinations of the misclassification rate (R_vect)
  # and the skew (pi_vect)
  grid <- expand.grid(R = R_vect,
                      pi = pi_vect)

  # One doesn't need too many simulations in the current set-up, since the noise
  # in the simulated data for the contour plots is almost zero.
  comb <- prog_bar_sapply(1:nsims, function(x)
    compute_lambda_vect(propD = grid$pi,
                        prop_ME = grid$R,
                        rho_M_Y = rho_M_Y,
                        rho_M_D = rho_M_D,
                        rho_X1_D = rho_X1_D,
                        rho_X2_D = rho_X2_D,
                        rho_X2_X1 = rho_X2_X1,
                        seed = x))

  # Average over the simulations to get even smoother contours.
  grid$sim_lambda <- apply(comb,1,mean)

  return(list(toplot = grid,
              R_obs = R_obs,
              pi_obs = pi_obs))
}

#' Contour Plot for Misclassification Impact Visualization
#'
#' This function creates a contour plot visualizing the impact of misclassification rates
#' and treatment proportions on the bias or coverage of treatment effect estimates,
#' based on simulated data.
#'
#' @param dat A data frame containing the dataset.
#' @param m A fitted model object (e.g., from `lm`) used to extract control variables.
#' @param treatment A string specifying the name of the treatment variable in `dat`.
#' @param outcome A string specifying the name of the outcome variable in `dat`.
#' @param binary A string specifying the name of the observed binary treatment variable
#' subject to misclassification.
#' @param ground_truth Optional string specifying the name of the ground truth
#' treatment variable in `dat`.
#' @param confMat Optional 2x2 confusion matrix for the binary treatment variable.
#' @param rho_M_Y Optional numeric correlation between misclassification indicator
#' and outcome variable.
#' @param rho_M_D Optional numeric correlation between misclassification indicator
#' and treatment variable.
#' @param nsims Number of simulations to run for estimating lambda values.
#' Default is 100.
#' @param R_vect Numeric vector of misclassification rates to evaluate.
#' Default is `seq(0.025, 0.5, by = 0.025)`.
#' @param pi_vect Numeric vector of proportions of treated units to evaluate.
#' Default is `seq(0.05, 0.95, by = 0.05)`.
#' @param R_obs Optional numeric value specifying the observed or user-supplied misclassification rate.
#' @param pi_obs Optional numeric value specifying the observed or user-supplied proportion of treated units.
#'
#' @return A list with components:
#' \describe{
#'   \item{p}{A ggplot2 object representing the contour plot.}
#'   \item{contour_dat}{The list output from `contour_plot_data` containing the data used for plotting.}
#' }
#'
#' @details
#' This function internally calls `contour_plot_data` to generate the data grid and then plots contour lines
#' representing the simulated lambda values. The plot highlights the observed or user-supplied misclassification
#' rate and treatment proportion with a diamond marker.
#'
#' @examples
#' \dontrun{
#' plot_res <- contour_plot(dat = mydata,
#'                          m = lm(Y ~ D + X1 + X2, data = mydata),
#'                          treatment = "D",
#'                          outcome = "Y",
#'                          binary = "D_obs",
#'                          ground_truth = "D_true",
#'                          nsims = 50)
#' plot_res$p  # to display the plot
#' }
#'
#' @import ggplot2
#' @importFrom tidyr drop_na
#' @importFrom directlabels geom_dl
#' @importFrom scales percent
#'
#' @export
contour_plot <- function(dat,
                         m,
                         treatment,
                         outcome,
                         binary,
                         ground_truth = NULL,
                         confMat = NULL,
                         rho_M_Y = NULL,
                         rho_M_D = NULL,
                         nsims = 100,
                         R_vect = seq(0.025,.5,by = .025),
                         pi_vect = seq(.05,.95,by = .05),
                         R_obs = NULL,
                         pi_obs = NULL) {

  cat("... Generating data")

  # First, generate the contour data
  contour_dat <- contour_plot_data(dat = dat,
                                   m = m,
                                   treatment = treatment,
                                   outcome = outcome,
                                   binary = binary,
                                   ground_truth = ground_truth,
                                   confMat = confMat,
                                   R_obs = R_obs,
                                   rho_M_Y = rho_M_Y,
                                   rho_M_D = rho_M_D,
                                   nsims = nsims,
                                   R_vect = R_vect,
                                   pi_vect = pi_vect)

  ### Get number of gradient areas
  #### Check range of sim_lambda
  sim_lambda = contour_dat$toplot$sim_lambda
  range_sim_lambda <- range(sim_lambda, na.rm = TRUE)
  #### Breakdown pos/zero/neg
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
  range_sim_lambda_rd <- round(range_sim_lambda, 1)
  grad_labels <- seq(range_sim_lambda_rd[1], range_sim_lambda_rd[2], by = 0.1)

  cat("... Plotting")

  # Draw the plot
  p <- contour_dat$toplot %>%
    drop_na(sim_lambda) %>%
    ggplot(aes(x = R,
               y = pi,
               z = sim_lambda)) +

    ## Filled contours
    geom_contour_filled() +

    ## Positive contour lines
    {if(length(pos_breaks) > 0) geom_contour(breaks = pos_breaks,
                                             linewidth = 0.1,
                                             linetype = 'dashed',
                                             color = 'grey30')} +
    ## Zero contour line
    {if(length(zero_break) > 0) geom_contour(breaks = zero_break,
                                             linewidth = 1,
                                             linetype = 'solid',
                                             color = 'black')} +
    ## Negative contour lines
    {if(length(neg_breaks) > 0) geom_contour(breaks = neg_breaks,
                                             linewidth = 0.1,
                                             linetype = 'dashed',
                                             color = 'grey30')} +
    ## Observed point
    geom_point(x = contour_dat$R_obs,
               y = contour_dat$pi_obs,
               shape = 'diamond', size = 3.5) +
    geom_point(x = contour_dat$R_obs,
               y = contour_dat$pi_obs,
               shape = 23, size = 1.5, fill = 'white') +
    ## Labels
    directlabels::geom_dl(aes(label = after_stat(level)),
                          method = list("smart.grid", cex = 0.8),
                          stat = "contour",
                          breaks = grad_labels[grad_labels > 0],
                          color = 'black',
                          size = 1) +
    directlabels::geom_dl(aes(label = after_stat(level)),
                          method = list("smart.grid", cex = 0.8),
                          stat = "contour",
                          breaks = grad_labels[grad_labels <= 0],
                          color = 'grey30',
                          size = 1) +
    ## Legend and axes
    theme_minimal() +
    theme(legend.position = 'none') +
    labs(x = 'Misclassification Rate',
         y = 'Proportion of 1s') +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    scale_fill_manual(values = generate_palette(sim_lambda))

  return(list(p = p,
              contour_dat = contour_dat))
}
#' Black-and-White Contour Plot for Misclassification Impact Visualization
#'
#' This function creates a black-and-white contour plot visualizing the impact of
#' misclassification rates and treatment proportions on the bias or coverage of treatment
#' effect estimates, based on simulated data.
#'
#' @param dat A data frame containing the dataset.
#' @param m A fitted model object (e.g., from `lm`) used to extract control variables.
#' @param treatment A string specifying the name of the treatment variable in `dat`.
#' @param outcome A string specifying the name of the outcome variable in `dat`.
#' @param binary A string specifying the name of the observed binary treatment variable
#' subject to misclassification.
#' @param ground_truth Optional string specifying the name of the ground truth
#' treatment variable in `dat`.
#' @param confMat Optional 2x2 confusion matrix for the binary treatment variable.
#' @param rho_M_Y Optional numeric correlation between misclassification indicator
#' and outcome variable.
#' @param rho_M_D Optional numeric correlation between misclassification indicator
#' and treatment variable.
#' @param nsims Number of simulations to run for estimating lambda values.
#' Default is 100.
#' @param R_vect Numeric vector of misclassification rates to evaluate.
#' Default is `seq(0.025, 0.5, by = 0.025)`.
#' @param pi_vect Numeric vector of proportions of treated units to evaluate.
#' Default is `seq(0.05, 0.95, by = 0.05)`.
#' @param R_obs Optional numeric value specifying the observed or user-supplied
#' misclassification rate.
#' @param pi_obs Optional numeric value specifying the observed or user-supplied
#' proportion of treated units.
#'
#' @return A list with components:
#' \describe{
#'   \item{p}{A ggplot2 object representing the contour plot.}
#'   \item{contour_dat}{The list output from `contour_plot_data_bw` containing the data
#'   used for plotting.}
#' }
#'
#' @details
#' This function internally calls `contour_plot_data` to generate the data grid and then
#'  plots contour lines representing the simulated lambda values. The plot highlights the
#' observed or user-supplied misclassification rate and treatment proportion with a diamond marker.
#'
#' @examples
#' \dontrun{
#' plot_res <- contour_plot_bw(dat = mydata,
#'                          m = lm(Y ~ D + X1 + X2, data = mydata),
#'                          treatment = "D",
#'                          outcome = "Y",
#'                          binary = "D_obs",
#'                          ground_truth = "D_true",
#'                          nsims = 50)
#' plot_res$p  # to display the plot
#' }
#'
#' @import ggplot2
#' @importFrom tidyr drop_na
#' @importFrom directlabels geom_dl
#' @importFrom scales percent
#'
#' @export
contour_plot_bw <- function(dat,
                            m,
                            treatment,
                            outcome,
                            binary,
                            ground_truth = NULL,
                            confMat = NULL,
                            rho_M_Y = NULL,
                            rho_M_D = NULL,
                            nsims = 100,
                            R_vect = seq(0.025, .5, by = .025),
                            pi_vect = seq(.05, .95, by = .05),
                            R_obs = NULL,
                            pi_obs = NULL) {

  cat("... Generating data")
  contour_dat <- contour_plot_data(dat = dat,
                                   m = m,
                                   treatment = treatment,
                                   outcome = outcome,
                                   binary = binary,
                                   ground_truth = ground_truth,
                                   confMat = confMat,
                                   R_obs = R_obs,
                                   rho_M_Y = rho_M_Y,
                                   rho_M_D = rho_M_D,
                                   nsims = nsims,
                                   R_vect = R_vect,
                                   pi_vect = pi_vect)

  sim_lambda <- contour_dat$toplot$sim_lambda
  range_sim_lambda <- range(sim_lambda, na.rm = TRUE)

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
  zero_break <- if (range_sim_lambda[1] <= 0 && range_sim_lambda[2] >= 0) 0 else numeric(0)
  range_sim_lambda_rd <- round(range_sim_lambda, 1)
  grad_labels <- seq(range_sim_lambda_rd[1], range_sim_lambda_rd[2], by = 0.1)

  cat("... Plotting")
  p <- contour_dat$toplot %>%
    drop_na(sim_lambda) %>%
    ggplot(aes(x = R, y = pi, z = sim_lambda)) +

    ## Contour lines only, no fill
    {if (length(pos_breaks) > 0) geom_contour(breaks = pos_breaks,
                                              color = "black",
                                              linetype = "dashed",
                                              linewidth = 0.5)} +
    {if (length(zero_break) > 0) geom_contour(breaks = zero_break,
                                              color = "black",
                                              linetype = "solid",
                                              linewidth = 1)} +
    {if (length(neg_breaks) > 0) geom_contour(breaks = neg_breaks,
                                              color = "black",
                                              linetype = "dotted",
                                              linewidth = 0.5)} +

    ## Observed point
    geom_point(x = contour_dat$R_obs,
               y = contour_dat$pi_obs,
               shape = 23, size = 3.5, fill = "black") +

    ## Labels
    directlabels::geom_dl(aes(label = after_stat(level)),
                          method = list("smart.grid", cex = 0.8),
                          stat = "contour",
                          breaks = grad_labels[grad_labels > 0],
                          color = "black",
                          size = 3) +
    directlabels::geom_dl(aes(label = after_stat(level)),
                          method = list("smart.grid", cex = 0.8),
                          stat = "contour",
                          breaks = grad_labels[grad_labels <= 0],
                          color = "black",
                          size = 3) +

    ## Theme and axes
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Misclassification Rate",
         y = "Proportion of 1s")

  return(list(p = p,
              contour_dat = contour_dat))
}
