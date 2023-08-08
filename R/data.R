#' Toy dataset generated through [`pdRCON.simulate`] function
#'
#' Data simulated by the function [`pdRCON.simulate`] with parameters:
#' * `p=20`
#' * `type=c("v", "i","a")`
#' * `force.symm=NULL`
#' * `dens=0.5`
#' * `Sigma` a \eqn{p \times p} matrix with unitary diagonal and 0.5 on off-diagonal elements.
#'
#' @format ## `toy_data`
#' A list with the following components:
#'
#' * `pdColG` a matrix representing a coloured graph for paired data; see [`pdglasso-package`] for details.
#' * `K`, a concentration matrix adapted to `pdColG`.
#' * `sample.data`, a data frame with 250 rows and 20 columns from a multivariate normal distribution
#' with zero mean vector and concentration matrix `K`.
#'
"toy_data"
