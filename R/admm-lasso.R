
#' ADMM graphical lasso algorithm for pdRCON models
#'
#' Estimate of a concentration matrix within the pdRCON 
#' submodel class identified by the arguments `type` and `force.symm`, 
#' based on the pdglasso method with penalty terms `lambda1` and `lambda2`
#' (Ranciati and Roverato, 2024). Optimization is implemented by an
#' ADMM algorithm (see Boyd *et. al.* 2011). The output is an object of 
#' class `ADMMoutput` and the user may then call [`pdColG.get`] to obtain the Coloured Graph
#' for Paired Data (pdColG) representing the selected model.
#'
#' @param S a sample covariance matrix with the block structure
#'   described in [`pdglasso-package`].
#' @param lambda1 a non-negative scalar penalty that encourages
#'   sparsity in the concentration matrix.
#' @param lambda2 a non-negative scalar penalty that encourages
#'   equality constraints in the concentration matrix.
#' @param type,force.symm two subvectors of `c("vertex", "inside.block.edge",
#'   "across.block.edge")` which identify the pdRCON submodel class of interest; see
#'   [`pdglasso-package`] for details.
#' @param X.init a \eqn{p \times p}
#'   concentration matrix to be used as starting point of the ADMM. If `NULL` a default value is used. 
#' @param rho1 a postive scalar; tuning parameter of the ADMM algorithm to be used for
#'   the outer loop.
#' @param rho2 a positive scalar; tuning parameter of the ADMM algorithm to be used for
#'   the inner loop. 
#' @param varying.rho1 a logical; if `TRUE` the parameter `rho1` is updated
#'   iteratively to speed-up convergence.
#' @param varying.rho2 a logical; if `TRUE` the parameter `rho2` is updated
#'   iteratively to speed-up convergence.
#' @param max_iter an integer; maximum number of iterations for convergence. 
#' @param eps.abs a  positive scalar; the absolute tollerance for the computation
#'   of primal and dual feasibility tollerances of the ADMM.
#' @param eps.rel a positive scalar; the relative tollerance for the computation
#'   of primal and dual feasibility tollerances of the ADMM.
#' @param rcpp a logical; if `TRUE`, computations are performed using the Rcpp (C++) implementation;
#' if `FALSE`, a pure R implementation is used.
#' @param verbose a logical; if `TRUE` the progress (and internal convergence of
#'   inner loop) is shown in the console while the algorithm is running.
#' @param print.type a logical; if `TRUE` the pdRCON submodel class considered, as
#'   specified by the arguments `type` and `force.symm`, is returned as printed
#'   output in the console. This option is only used when `rcpp = FALSE`; when
#'   `rcpp = TRUE`, no output is printed for efficiency.
#'
#' @return A object of class `ADMMoutput` that is a list with the following components:
#'
#' * `X` the estimated concentration matrix.
#'
#' * `acronyms` a vector of strings identifying the pdRCON submodel class
#'  considered as identified  by the arguments `type` and `force.symm`.
#'
#' * `internal.par` a list of internal parameters passed to the function at the
#' call, as well as convergence information.
#'
#' @export
#'
#' @references Ranciati, S. and Roverato, A., (2024). On the application of Gaussian graphical models to paired data problems. *Statistics and Computing*, *34*(6), 1-19.
#' @references Boyd, S., Parikh, N., Chu, E., Borja, P. and Eckstein, J. (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers. *Foundations and Trends in Machine learning*, *3*(1), 1-122.
#' @examples
#' 
#' # computation of the sample covariance
#' S <- cov(toy_data$sample.data)
#' 
#' # model with all types of symmetries allowed and no full symmetry required
#' admm.out <- admm.pdglasso(S, lambda1=4, lambda2=0.7)
#' G <- pdColG.get(admm.out)
#' summary(G)
#' 
#' # model with no across-block symmetries allowed and full vertex-symmetry required
#' admm.out <-admm.pdglasso(S, , lambda1=4, lambda2=0.7, type=c("v", "i"), force.symm = c("v"))
#' G <- pdColG.get(admm.out)
#' summary(G)
#' 

admm.pdglasso <- function(S,
                          lambda1,
                          lambda2,
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
  
  return(as.ADMMoutput(
    list(
      X = X,
      acronyms = acronyms,
      internal.par = internal.par
        )
    ))
}
