#' Toy dataset generated through [`simul.pdColG`] function
#'
#' Data simulated by the function [`simul.pdColG`] with a call:
#' `toy_data <- pdRCON.simulate(20, type=c("v", "i","a), force.symm=NULL, dens=0.3)`
#'
#' @format ## `toy_data`
#' A list containing three elements:
#' * pdColG, a \eqn{p \times p} matrix describing the coloured graphical model,
#' * K, the concentration matrix associated to the model,
#' * sample.data, a data frame with 60 rows and 20 columns.
"toy_data"
