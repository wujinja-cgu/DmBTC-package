#' Plot the forest plot
#'
#' @param bugs.out The output data from DmBTC model.
#' @return A list containing the results of the model fitting.
#' @export

#plot_forest <- function(data) {
#
#  require(meta4diag)
#  require(meta)
  # Extract individual study estimates (e.g., sensitivity) from the bugs.out object
#  bugs.out = data
#  sensitivity_estimates <- bugs.out$sims.list$pool.se
#  lower_ci <- bugs.out$summary["pool.se", "2.5%"]
#  upper_ci <- bugs.out$summary["pool.se", "97.5%"]

  # Assuming data contains study-level sensitivity estimates
#  study_data <- data.frame(
#    study = 1:length(sensitivity_estimates),
#    sensitivity = sensitivity_estimates,
#    lower_ci = lower_ci,
#    upper_ci = upper_ci
#  )

#  forestplot <- meta::forest(
#    x = study_data$sensitivity,
#    ci.lb = study_data$lower_ci,
#    ci.ub = study_data$upper_ci,
#    studlab = paste("Study", study_data$study),
#    xlab = "Sensitivity",
#    sm = "Sensitivity"
#  )

#  return(forestplot)
#}

plot_forest <- function(bugs.out) {

  # Load required packages
  require(meta)

  # Example: Extract individual study estimates (assuming 'sens' stores study-level data)
  sensitivity_estimates <- bugs.out[[10]][1]
  lower_ci <- bugs.out[[10]][3]
  upper_ci <- bugs.out[[10]][7]

  # Create a dataframe with study data
  study_data <- data.frame(
    study = 1:length(sensitivity_estimates),
    sensitivity = sensitivity_estimates,
    lower_ci = lower_ci,
    upper_ci = upper_ci
  )

  # Plot the forest plot using meta::forest
  forestplot <- meta::forest(
    x = study_data$sensitivity,
    ci.lb = study_data$lower_ci,
    ci.ub = study_data$upper_ci,
    studlab = paste("Study", study_data$study),
    xlab = "Sensitivity",
    sm = "Sensitivity"
  )

  return(forestplot)
}


# Install necessary packages
install.packages("forestplot")
library(forestplot)

# Example data: Replace this with your actual sensitivity and specificity data
studies <- c("Study 1", "Study 2", "Study 3", "Pooled")
sensitivities <- c(0.85, 0.70, 0.90, 0.80)  # Point estimates for sensitivity
specificities <- c(0.95, 0.85, 0.88, 0.90)  # Point estimates for specificity

# Lower and upper bounds for 95% credible intervals
ci_sens_lower <- c(0.75, 0.60, 0.80, 0.78)
ci_sens_upper <- c(0.95, 0.80, 1.00, 0.85)
ci_spec_lower <- c(0.90, 0.70, 0.75, 0.85)
ci_spec_upper <- c(1.00, 0.90, 0.95, 0.92)

# Create table for forest plot
tabletext <- cbind(
  c("Study", studies),
  c("Sensitivity", paste0(sensitivities, " (", ci_sens_lower, "-", ci_sens_upper, ")")),
  c("Specificity", paste0(specificities, " (", ci_spec_lower, "-", ci_spec_upper, ")"))
)

# Sensitivity plot
forestplot(labeltext = tabletext,
           mean = cbind(NA, sensitivities),  # Mean point estimates for sensitivity
           lower = cbind(NA, ci_sens_lower),  # Lower bound for CI
           upper = cbind(NA, ci_sens_upper),  # Upper bound for CI
           is.summary = c(TRUE, rep(FALSE, length(studies)-1)),  # Summary row for pooled
           xlab = "Sensitivity",
           col = forestplot::fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
           zero = 1,  # Line at 1
           boxsize = 0.1)

# Specificity plot
forestplot(labeltext = tabletext,
           mean = cbind(NA, specificities),  # Mean point estimates for specificity
           lower = cbind(NA, ci_spec_lower),  # Lower bound for CI
           upper = cbind(NA, ci_spec_upper),  # Upper bound for CI
           is.summary = c(TRUE, rep(FALSE, length(studies)-1)),  # Summary row for pooled
           xlab = "Specificity",
           col = forestplot::fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
           zero = 1,  # Line at 1
           boxsize = 0.1)
