
#' ADMM graphical lasso algorithm for coloured GGMs for paired data.
#'
#' By providing a covariance matrix `S` and values for `lambda1` and `lambda2`,
#' this function estimates a concentration matrix `X` within the pdRCON
#' submodel class, identified by the arguments `type` and `force.symm`, based on the pdglasso method (Ranciati & Roverato, 2024) using an
#' (adaptive) ADMM algorithm. The output is the matrix and a list of internal
#' parameters used by the function, together with the specific call with the
#' relevant pdRCON submodel class.
#'
#' @param S a sample covariance (or correlation) matrix with the block structure
#'   described in [`pdglasso-package`].
#' @param lambda1 a non-negative scalar (or vector) penalty that encourages
#'   sparsity in the concentration matrix.
#' @param lambda2 a non-negative scalar (or vector) penalty that encourages
#'   equality constraints in the concentration matrix.
#' @param type,force.symm two subvectors of `c("vertex", "inside.block.edge",
#'   "across.block.edge")` which identify the pdRCON submodel class of interest; see
#'   [`pdglasso-package`] for details.
#' @param X.init (optional) a \eqn{p \times p} initial guess for the
#'   concentration matrix and/or starting solution for the ADMM algorithm.
#' @param rho1 a scalar; tuning parameter of the ADMM algorithm to be used for
#'   the outer loop. It must be strictly positive.
#' @param rho2 a scalar; tuning parameter of the ADMM algorithm to be used for
#'   the inner loop. It must be strictly positive.
#' @param varying.rho1 a logical; if `TRUE` the parameter rho1 is updated
#'   iteratively to speed-up convergence.
#' @param varying.rho2 a logical; if `TRUE` the parameter rho2 is updated
#'   iteratively to speed-up convergence.
#' @param max_iter an integer; maximum number of iterations to be run in case
#'   the algorithm does not converge.
#' @param eps.abs a scalar; the absolute precision required for the computation
#'   of primal and dual residuals of the ADMM algorithm.
#' @param eps.rel a scalar; the relative precision required for the computation
#'   of primal and dual residuals of the ADMM algorithm.
#' @param rcpp a logical; if `TRUE`, computations are performed using the Rcpp (C++) implementation;
#' if `FALSE`, a pure R implementation is used.
#' @param verbose a logical; if `TRUE` the progress (and internal convergence of
#'   inner loop) is shown in the console while the algorithm is running.
#' @param print.type a logical; if `TRUE` the pdRCON submodel class considered, as
#'   specified by the arguments `type` and `force.symm`, is returned as printed
#'   output in the console. This option is only used when `rcpp = FALSE`; when
#'   `rcpp = TRUE`, no output is printed for efficiency.
#'
#' @return A list with the following components:
#'
#' * `X` the estimated concentration matrix under the pdRCON submodel class
#'   considered and the values of `lambda1` and `lambda2`.
#'
#' * `acronyms` a vector of strings identifying the pdRCON submodel class
#'  considered as identified  by the arguments `type` and `force.symm`.
#'
#' * `internal.par` a list of internal parameters passed to the function at the
#' call, as well as convergence information.
#'
#' @export
#'
#' @references Ranciati, S., Roverato, A., (2024). On the application of Gaussian graphical models to paired data problems. *Statistics and Computing*, *34*(6), 1-19. \url{https://doi.org/10.1007/s11222-024-10513-6}
#' @references Ranciati, S., Roverato, A., (2023). On the application of Gaussian graphical models to paired data problems. *arXiv pre-print*. \url{https://arxiv.org/abs/2307.14160}
#' @examples
#' S <- cov(toy_data$sample.data)
#' admm.pdglasso(S)

admm.pdglasso <- function(S,
                          lambda1      = 1,
                          lambda2      = 1e-4,
                          type         = c("vertex", "inside.block.edge", "across.block.edge"),
                          force.symm   = NULL,
                          X.init       = NULL,
                          rho1         = 1,
                          rho2         = 1,
                          varying.rho1 = TRUE,
                          varying.rho2 = TRUE,
                          max_iter     = 5000,
                          eps.abs      = 1e-6,
                          eps.rel      = 1e-6,
                          rcpp         = TRUE,
                          print.type   = TRUE,
                          verbose      = FALSE) {
  
  stopifnot(is.matrix(S), nrow(S) == ncol(S))
  stopifnot(is.numeric(lambda1), all(lambda1 >= 0))
  stopifnot(is.numeric(lambda2), all(lambda2 >= 0))
  stopifnot(rho1 > 0, rho2 > 0)
  
  p <- nrow(S)
  if (p %% 2 != 0) stop("p must be even")
  q <- p/2
  
  time.start <- Sys.time()
  
  out.make.a <- make.acronyms(type, force.symm, print.type = print.type)
  acr.type  <- out.make.a$acronym.of.type
  acr.force <- out.make.a$acronym.of.force
  
  if (!is.null(acr.force))
    lambda2 <- lambda2.force.symm(p, lambda2, acr.type, acr.force)
  
  n.row.F <- get.n.row.F(q, acr.type)
  
  if (is.null(X.init)) X <- diag(1, p) else X <- X.init
  
  if (rcpp){
    
    res <- admm_pdglasso_internal(
      S           = S,
      X           = X,
      lambda1     = lambda1,
      lambda2     = lambda2,
      rho1        = rho1,
      rho2        = rho2,
      varying_rho1= varying.rho1,
      varying_rho2= varying.rho2,
      max_iter    = max_iter,
      eps_abs     = eps.abs,
      eps_rel     = eps.rel,
      acr_type    = acr.type,
      n_row_F     = n.row.F
    )
    
  }else{
    
    res <- admm_pdglasso_internal_r(
      S           = S,
      X           = X,
      lambda1     = lambda1,
      lambda2     = lambda2,
      rho1        = rho1,
      rho2        = rho2,
      varying.rho1= varying.rho1,
      varying.rho2= varying.rho2,
      max_iter    = max_iter,
      eps.abs     = eps.abs,
      eps.rel     = eps.rel,
      acr.type    = acr.type,
      n.row.F     = n.row.F,
      verbose     = verbose
    )
    
  }
  
  time.diff <- Sys.time() - time.start
  
  internal.par <- list(
    execution.time = time.diff,
    res.primal = res$res_primal,
    res.dual = res$res_dual,
    lambda1 = lambda1,
    lambda2 = unique(lambda2),
    n.iter = res$n_iter,
    n.iter.rho1_update_last = res$n_iter_rho1_update_last,
    last.rho1 = res$last_rho1,
    eps.primal = res$eps_primal,
    eps.dual = res$eps_dual,
    eps.abs = eps.abs,
    eps.rel = eps.rel,
    max_iter = max_iter,
    converged = res$converged
  )
  
  acronyms <- list(acronym.of.type = acr.type,
                  acronym.of.force = acr.force)
  
  X = res$X
  dimnames(X) <- dimnames(S)
  
  return(list(
    X = X,
    acronyms = acronyms,
    internal.par = internal.par
  ))
}
