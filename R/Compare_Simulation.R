#' Run the DmBTC, Bamdit, and SCS Models
#'
#' @param data The data to be used in the model.
#' @return A list containing the results of the model fitting.
#' @export

Compare_Simulation <- function(data) {

  if (!is.data.frame(data)) {
    stop("Input data must be a data frame")
  }

  # Ensure required columns are present
  required_columns <- c("TP", "FP", "TN", "FN")

  if (!all(required_columns %in% colnames(data))) {
    stop("Data frame must contain the following columns: TP, FP, TN, FN")
  }

  mydata = data

  mydata_DmBTC=cbind(mydata$ID,mydata$TP,mydata$FN,mydata$TN,mydata$FP)

  colnames(mydata_DmBTC)=c('ID','TP','FN','TN','FP')

  mydata_DmBTC=as.data.frame(mydata_DmBTC)

  R <- structure(.Data = c(1, 0, 0, 1), .Dim = c(2, 2))

  bugs.data <- list(
    'R' = R,
    'tp' = mydata_DmBTC$TP,
    'fp' = mydata_DmBTC$FP,
    'n1' = mydata_DmBTC$TP + mydata_DmBTC$FN,
    'n2' = mydata_DmBTC$TN + mydata_DmBTC$FP,
    'n' = nrow(mydata_DmBTC)
  )

  bugs.parameters <- c("pool.se", "pool.sp", "new.se", "new.sp", "mu")

  fit_DmBTC <- tryCatch({
    DmBTCmeta(
      data = bugs.data,
      parameters = bugs.parameters,
      n.iter = 100000,
      n.burnin = 50000,
      n.chains = 4,
      n.thin = 1,
      digits = 3,
      debug = FALSE,
      inits = NULL
    )
  }, error = function(e) {
    message("Error in DmBTC: ", e$message)
    NULL
  })

  fit_SCS <- tryCatch({
    SCSmeta(
      tp = mydata_DmBTC$TP,
      fp = mydata_DmBTC$FP,
      tn = mydata_DmBTC$TN,
      fn = mydata_DmBTC$FN
    )
  }, error = function(e) {
    message("Error in SCSmeta: ", e$message)
    NULL
  })

  tp <- mydata_DmBTC$TP
  fp <- mydata_DmBTC$FP
  n1 <- mydata_DmBTC$TP + mydata_DmBTC$FN
  n2 <- mydata_DmBTC$FP + mydata_DmBTC$TN
  mydata2 <- cbind(tp, n1, fp, n2)

  fit_bamdit <- tryCatch({
    metadiag(
      mydata2,
      re = "normal",
      link = "logit",
      mean.mu.D = 0,
      mean.mu.S = 0,
      sd.mu.D = 10,
      sd.mu.S = 10,
      sigma.D.upper = 10,
      sigma.S.upper = 10,
      mean.Fisher.rho = 0,
      sd.Fisher.rho = 1/sqrt(2),
      df = 4,
      split.w = FALSE,
      n.1.new = 50,
      n.2.new = 50,
      nr.chains = 4,
      nr.iterations = 100000,
      nr.adapt = 1000,
      nr.burnin = 50000,
      nr.thin = 1,
      be.quiet = FALSE,
      r2jags = TRUE
    )
  }, error = function(e) {
    message("Error in metadiag: ", e$message)
    NULL
  })

  fit_results <- array(0, dim = c(4, 2))

  fit_results[1,1]=mydata[1,1]
  fit_results[1,2]=mydata[1,2]

  if (!is.null(fit_bamdit)) {
    fit_results[2, 1] <- fit_bamdit[[2]][10]$summary[5 + nrow(mydata_DmBTC) + 2]
    fit_results[2, 2] <- fit_bamdit[[2]][10]$summary[5 + nrow(mydata_DmBTC) + 4 + nrow(mydata_DmBTC) + 2]
  }

  if (!is.null(fit_SCS)) {
    fit_results[3, 1] <- fit_SCS[1]
    fit_results[3, 2] <- fit_SCS[2]
  }

  if (!is.null(fit_DmBTC)) {
    fit_results[4, 1] <- fit_DmBTC[[10]][1]
    fit_results[4, 2] <- fit_DmBTC[[10]][2]
  }

  colnames(fit_results) <- c('pooled sensitivity', 'pooled specificity')
  rownames(fit_results) <- c('TRUE','Bamdit', 'SCS', 'DmBTC')

  # Clean up objects to free memory
  #rm(fit_DmBTC, fit_SCS, fit_bamdit, bugs.data, bugs.parameters, tp, fp, n1, n2, mydata2)
  #gc()  # Trigger garbage collection

  return(fit_results)
}
