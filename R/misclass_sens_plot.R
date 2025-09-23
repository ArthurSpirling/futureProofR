#' Plot for misclassification sensitivity analysis results
#'
#' This function creates a plot showing the bounds on regression coefficients
#' under different levels of misclassification in a binary variable.
#'
#' @param misclass_results List. The output from misclass_sens() function.
#' @param title Character. Main title for the plot. Default: "Sensitivity Analysis"
#' @param subtitle Character. Optional subtitle. If NULL, automatically generates
#'   based on model formula and number of variables.
#' @param point_size Numeric. Size of the observed coefficient point. Default: 3
#' @param m Model object. The fitted regression model.
#' @param treatment Character. Name of the treatment variable.
#'
#' @return A ggplot object showing the sensitivity analysis bounds
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import scales
#'
#' @examples
#' # Run sensitivity analysis first
#' results <- misclass_sens(dat = your_data, outcome = "Y", treatment = "D",
#'                          binary = "D", m = your_model)
#' # Create plot
#' misclass_sens_plot(results)
#'
#' @export
misclass_sens_plot <- function(misclass_results,
                               m,
                               treatment,
                               title = "Sensitivity Analysis",
                               subtitle = NULL,
                               point_size = 3) {

  # Get number of control variables using standardized helper
  ctrls <- extract_controls(m, treatment)
  n_ctrls <- length(ctrls)

  # Create the plot data
  plot_data <- misclass_results$bounds %>%
    select(misclassRate, matches('max|min')) %>%
    pivot_longer(cols = -misclassRate, names_to = "metric", values_to = "value") %>%
    mutate(coef = str_extract(metric, 'max|min'),
           bound = gsub('_max|_min', '', metric)) %>%
    filter(bound != 'perm_prob') %>%
    mutate(bound = gsub('_0\\.\\d', '', bound))

  # Get observed coefficient for the point
  obs_coef <- misclass_results$bounds$obsCoef[1]

  # Create the plot
  p <- plot_data %>%
    ggplot(aes(x = misclassRate, y = value,
               linetype = bound, group = metric)) +
    geom_line() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    theme_bw() +
    scale_linetype_manual(values = c('solid', 'dotted', 'dashed'),
                          labels = c('Extreme', 'Permutation', 'Naive')) +
    annotate(geom = 'point', x = 0, y = obs_coef,
             size = point_size, pch = 19) +
    scale_x_continuous(labels = scales::percent) +
    # coord_cartesian(ylim = ylim) +
    labs(x = '% Reclassified',
         y = 'Estimated Coefficient',
         linetype = 'Bound',
         title = title,
         subtitle = paste0("Binary Outcome ~ Treatment + ", n_ctrls, " controls")) +
    theme(legend.position = "bottom")

  return(p)
}

