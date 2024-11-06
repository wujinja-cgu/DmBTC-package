# Example of defining the DmBTC function in R/DmBTC.R
#' Run the DmBTC Model
#'
#' @param data The data to be used in the model.
#' @param parameters The parameters to be estimated.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of burn-in iterations.
#' @param n.chains Number of Markov chains.
#' @param n.thin Thinning interval.
#' @param digits Number of decimal places to display.
#' @param debug Whether to run in debug mode.
#' @param inits Initial values for the parameters.
#' @return A list containing the results of the model fitting.
#' @export

DmBTCmeta <- function(data,
                      parameters,
                      n.iter,
                      n.burnin,
                      n.chains,
                      n.thin,
                      digits,
                      debug,
                      inits = NULL) {

  ## load packages
  require(rstan)
  require(meta4diag)
  require(R2OpenBUGS)
  require(devtools)
  require(bamdit)
  require(meta)
  library(coda)

  rm()

  # Define the path to the WinBUGS model
  model.file <- system.file("extdata", "DmBTC_model.txt", package = "DmBTC")

  if (model.file == "") {
    stop("Model file does not exist in the specified path.")
  }

  # Run the WinBUGS model
  bugs.out <- bugs(data = data,
                   parameters = parameters,
                   model.file = model.file,
                   n.iter = n.iter,
                   n.burnin = n.burnin,
                   n.chains = n.chains,
                   n.thin = n.thin,
                   digits = digits,
                   DIC = TRUE,
                   codaPkg = TRUE,
                   debug = debug,
                   inits = NULL)
  #library(coda)
  #mcmc.output <- as.mcmc(bugs.out)
  #traceplot <- traceplot(mcmc.output)
  # Read the MCMC samples from WinBUGS output
  #mcmc.output <- read.bugs(bugs.out)

  # Plot the trace plot for all parameters
  #traceplot <- traceplot(mcmc.output)

  # Return both the model output and the plotting functions
  return(bugs.out)
}
