#' Fit and select a coloured GGM for paired data according to eBIC
#' criterion.
#'
#' Performs a sequence of calls to [`admm.pdglasso`] providing two grids of
#' values for `lambda1` and `lambda2`. First, a grid search conditional on
#' `lambda2`=0 is run to select the best `lambda1` value among the candidates
#' (according to eBIC); conditional on the best `lambda1`, a similar search is
#' performed for `lambda2`. The output is the select model, given by the
#' estimated concenration matrix.
#' The user may call [`pdColG.get`] to obtain the corresponding Coloured Graph
#' for Paired Data (pdColG) from the selected model.
#'
#' @inheritParams admm.pdglasso
#' @param n the sample size of the data used to compute the sample covariance
#'   matrix S.
#' @param lams a 2x4 matrix; first row refers to `lambda1` and second row to
#'   `lambda2`; for each row, values are (i) minimum value for the grid; (ii)
#'   maximum value for the grid; (iii) number of points in the grid; (iv) if a
#'   logarithmic spacing (TRUE) is desired for the grid or not (FALSE); if `NULL` defaulta values are
#'   used, i.e. from max.lams to max.lams/20, and 20 grid points, log spacing for both.
#' @param gamma.eBIC the parameter for the eBIC computation. gamma=0 is equivalent to BIC.
#' @param progress a logical value; if `TRUE` provides a visual update in the console about the grid search over `lambda1` and `lambda2`
#' @param mle a logical; if `TRUE`, estimates the eBIC via the MLE (more accurate but slower computation);
#' if `FALSE`, uses the pdglasso (ADMM) estimator (faster but potentially less precise).
#'
#' @return A list with the following components:
#'
#' * `model` selected model; object resulting from the function [`admm.pdglasso`] with `best.lambdas`.
#' * `lambda.grid` the grid of values used for `lambda1` and `lambda2.`
#' * `best.lambdas` the selected values of `lambda1` and `lambda2` according to eBIC criterion.
#' * `l1.path` a matrix containing the grid values for `lambda1` as well as quantities used in eBIC computation.
#' * `l2.path` a matrix containing the grid values for `lambda2` as well as quantities used in eBIC computation.
#' * `time.exec` total execution time for the called function.
#'
#' A warning is produced if at least one run of the algorithm for the grid
#' searches has resulted in non-convergence (status can be checked by inspecting
#' `l1.path` and `l2.path`).
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' sel.mod <- pdRCON.fit(S, n=60)
#' sel.mod$l1.path
#' sel.mod$l2.path
#' pdRCON.check(sel.mod$model)
#' pdColG.get(sel.mod$model)

pdRCON.fit <- function(S,
                       n,
                       lams        = NULL,
                       gamma.eBIC  = 0.5,
                       type        = c("vertex", "inside.block.edge", "across.block.edge"),
                       force.symm  = NULL,
                       X.init      = NULL,
                       rho1        = 1,
                       rho2        = 1,
                       varying.rho1= TRUE,
                       varying.rho2= TRUE,
                       max_iter    = 5000,
                       eps.abs     = 1e-08,
                       eps.rel     = 1e-08,
                       verbose     = FALSE,
                       progress    = TRUE,
                       print.type  = TRUE,
                       mle.estimate = TRUE){
  start.time <- Sys.time()
  ## Max values for lambda_1 and lambda_2 according to theorems; only needed for the grid search
  
  if(is.null(lams)){
    lams <- matrix(0,2,4)
    lams[,2] <- lams.max(S)
    lams[,1] <- lams[,2]/20
    lams[,3] <- c(20,20)
    lams[,4] <- c(TRUE,TRUE)
  }
  rownames(lams) <- c("l1","l2")
  colnames(lams) <- c("min","max","n.pts","log.spacing")
  
  ## Prepare temp objects
  eBIC.l1 <-  matrix(0,lams[1,3],4)
  eBIC.l2 <-  matrix(0,lams[2,3],4)
  
  ### First grid search for lambda_1, with lambda_2=0
  if(lams[1,4]==1) l1.vec <- exp(seq(log(lams[1,1]),log(lams[1,2]), length.out=lams[1,3])) else l1.vec <- seq(lams[1,1], lams[1,2], length.out=lams[1,3])
  l1.vec <- sort(l1.vec, decreasing=TRUE)
  for(i in 1:lams[1,3]){
    if(progress==TRUE) cat("Searching over lambda1 grid (",i,"/",lams[1,3],").\n", sep="")
    mod.out <- admm.pdglasso(S,
                             lambda1=l1.vec[i],
                             lambda2=0,
                             type=type,
                             force.symm=force.symm,
                             X.init=X.init,
                             rho1=rho1,
                             rho2=rho2,
                             varying.rho1=varying.rho1,
                             varying.rho2=varying.rho2,
                             max_iter=max_iter,
                             eps.abs=eps.abs,
                             eps.rel=eps.rel,
                             verbose=FALSE,
                             print.type=FALSE)
    eBIC.l1[i,1:3] <- compute.eBIC(S=S,
                                   mod=mod.out,
                                   n=n,
                                   gamma.eBIC=gamma.eBIC,
                                   max_iter=max_iter,
                                   mle = mle.estimate)
    eBIC.l1[i,4] <- mod.out$internal.par$converged+0
    if(eBIC.l1[i,4]==0) cat("Convergence not achieved for this value of lambda1! \n")
  }
  best.l1 <- l1.vec[which.min(eBIC.l1[,1])]
  if(length(best.l1)==0) stop("Grid search of lambda1 failed!")
  
  if(progress==TRUE) cat("--- \n", sep="")
  
  ### Second grid search for lambda_2, with lambda_1=best.l1
  if(lams[2,4]==1) l2.vec <- exp(seq(log(lams[2,1]),log(lams[2,2]), length.out=lams[2,3])) else l2.vec <-  seq(lams[2,1], lams[2,2], length.out=lams[2,3])
  l2.vec <- sort(l2.vec, decreasing=TRUE)
  for(i in 1:lams[2,3]){
    if(progress==TRUE) cat("Searching over lambda2 grid (",i,"/",lams[2,3],").\n", sep="")
    mod.out <- admm.pdglasso(S,
                             lambda1=best.l1,
                             lambda2=l2.vec[i],
                             type=type,
                             force.symm=force.symm,
                             X.init=X.init,
                             rho1=rho1,
                             rho2=rho2,
                             varying.rho1=varying.rho1,
                             varying.rho2=varying.rho2,
                             max_iter=max_iter,
                             eps.abs=eps.abs,
                             eps.rel=eps.rel,
                             verbose=FALSE,
                             print.type=FALSE)
    eBIC.l2[i,1:3] <- compute.eBIC(S=S,
                                   mod=mod.out,
                                   n=n,
                                   gamma.eBIC=gamma.eBIC,
                                   max_iter=max_iter,
                                   mle = mle.estimate)
    eBIC.l2[i,4] <- mod.out$internal.par$converged+0
    if(eBIC.l2[i,4]==0) cat("Convergence not achieved for this value of lambda2! \n")
  }
  ### adding eBIC value/results and l2 value to path
  ### already estimated from the first grid.search where lam2=0
  eBIC.l2 <- rbind(eBIC.l2, eBIC.l1[which.min(eBIC.l1[,1]),])
  l2.vec <- c(l2.vec,0)
  best.l2 <- l2.vec[which.min(eBIC.l2[,1])]
  if(length(best.l2)==0) stop("Grid search of lambda2 failed!")
  
  ## Fit final model
  mod.out <- admm.pdglasso(S,
                           lambda1=best.l1,
                           lambda2=best.l2,
                           type=type,
                           force.symm=force.symm,
                           X.init=X.init,
                           rho1=rho1,
                           rho2=rho2,
                           varying.rho1=varying.rho1,
                           varying.rho2=varying.rho2,
                           max_iter=max_iter,
                           eps.abs=eps.abs,
                           eps.rel=eps.rel,
                           verbose=verbose,
                           print.type=print.type)
  
  l1.path=cbind(l1.vec,eBIC.l1)
  l2.path=cbind(l2.vec,eBIC.l2)
  colnames(l1.path) <- c("lambda1.grid", "eBIC     ","  log-Likelihood  ","num. of params", "converged (1=TRUE)")
  colnames(l2.path) <- c("lambda2.grid", "eBIC     ","  log-Likelihood  ","num. of params", "converged (1=TRUE)")
  
  time.exec <- Sys.time()-start.time
  return(list(model=mod.out,
              best.lambdas=c(best.l1,best.l2),
              lambda.grid=lams,
              l1.path=l1.path,
              l2.path=l2.path,
              time.exec=time.exec))
  
}




#' Check orders of magnitude of entries of the estimated concentration matrix X under the pdRCON submodel class considered.
#'
#' This function produces different plots of values of X in [`log10`] scale,
#' depending on the submodel class identified by the acronym stored in the input
#' object `mod.out`. The user might want to call this function to identify
#' what values to pass as arguments `th1` and `th2` to a call to the function
#' [`pdColGet`].
#'
#' @param mod.out An object of list type, that is the output of a call to the [`admm.pdglasso`] function. This can also be the output of a call to [`pdRCON.fit`].
#'
#' @return Depending on the acronym stored inside the input object, one or more plots depicting:
#'
#' * off-diagonal elements,
#' * vertices,
#' * inside-block,
#' * across-block.
#' @importFrom graphics par title abline legend
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' mod.out <- pdRCON.fit(S,n=60)$model
#' pdRCON.check(mod.out)
pdRCON.check <- function(mod.out){
  acronyms <- mod.out$acronyms$acronym.of.type
  th.rel <- log10(mod.out$internal.par$eps.rel)
  th.primal <- log10(mod.out$internal.par$eps.primal)
  th.dual <- log10(mod.out$internal.par$eps.dual)
  p <- ncol(mod.out$X)
  q <- p/2
  
  # Prepare canvas
  def.pars <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2))
  par(mar = c(2, 2, 1.5, 1))
  #
  # off-diagonal
  off.d <- mod.out$X[upper.tri(mod.out$X, diag=FALSE)]
  off.d <- log10(abs(off.d))
  off.d <- off.d[is.finite(off.d)]
  plot(off.d, ylim=c(min(off.d), 0))
  title(main="Off-diagonal elements")
  abline(h=c(th.rel, th.primal, th.dual), lty=2, col=c("red", "green", "blue"), lwd=2)
  #
  # vertices
  if(grepl("V",acronyms,fixed=T)){
    vertx <- diag(mod.out$X[1:q, 1:q]-mod.out$X[(q+1):p,(q+1):p])
    vertx <- log10(abs(vertx))
    vertx <- vertx[is.finite(vertx)]
    plot(vertx, ylim=c(min(off.d), 0))
    title(main="Diff. of 'Vertices' elements")
    abline(h=c(th.rel, th.primal, th.dual), lty=2, col=c("red", "green", "blue"), lwd=2)
  }
  #
  # inside-block edges
  if(grepl("I",acronyms,fixed=T)){
    inside.edges <- LL.block(mod.out$X)-RR.block(mod.out$X)
    inside.edges <- inside.edges[upper.tri(inside.edges, diag=FALSE)]
    inside.edges <- log10(abs(inside.edges))
    inside.edges <- inside.edges[is.finite(inside.edges)]
    plot(inside.edges, ylim=c(min(off.d), 0))
    title(main="Diff. of 'Inside' elements")
    abline(h=c(th.rel, th.primal, th.dual), lty=2, col=c("red", "green", "blue"), lwd=2)
  }
  #
  # across-block edges
  if(grepl("A",acronyms,fixed=T)){
    across.edges <- across.block(mod.out$X)-t(across.block(mod.out$X))
    across.edges <- across.edges[upper.tri(across.edges, diag=FALSE)]
    across.edges <- log10(abs(across.edges))
    across.edges <- across.edges[is.finite(across.edges)]
    plot(across.edges, ylim=c(min(off.d), 0))
    title(main="Diff. of 'Across' elements")
    abline(h=c(th.rel, th.primal, th.dual), lty=2, col=c("red", "green", "blue"), lwd=2)
  }
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("Rel. Eps.", "Primal Res.", "Dual Res."),
         col=c("red", "green", "blue"), lty=2, lwd=1.5, cex=0.8, horiz = TRUE)
  
  # Reset canvas
  on.exit(par(def.pars))
  #
}




#' Maximum likelihood estimation.
#'
#' Computes the maximum likelihood estimate of the concentration matrix of a pdRCON model.
#'
#' @param S a sample covariance matrix with the block structure described in [`pdglasso-package`].
#' @param pdColG pdColG a matrix representing a coloured graph for paired data; see [`pdglasso-package`] for details.
#' @inheritParams admm.pdglasso
#'
#' @return Either a matrix, that is the maximum likelihood estimate of
#'   \eqn{K=\Sigma^{-1}} under the pdRCON model represented by `pdColG`, or
#'   `NULL` if the maximum likelihood estimate does not exist.
#'
#' @details If the sample covariance matrix is not full-rank, then it is
#'   possible that the maximum likelihood estimate does not exist. The maximum
#'   likelihood estimate is computed by running the function [`admm.pdglasso`]
#'   with suitable penalties and, if it does not exist, then the ADMM algorithm
#'   fails to converge, a warning is produced and a `NULL` is returned.
#'
#' @export
#'
#' @examples
#' S <- var(toy_data$sample.data)
#' K.hat <- pdRCON.mle(S, toy_data$pdColG)

pdRCON.mle <- function(S, pdColG,
                       eps.rel=1e-6,
                       eps.abs=1e-6,
                       max_iter=5000){
  # make vector lambda1
  lambda1 <- (mat2vec(pdColG)==0)
  lambda1[lambda1] <- Inf
  
  # make vector lambda2
  p <- nrow(pdColG)
  q <- p/2
  m2v <- c(diag(pdColG[1:q,1:q]), half_vec(pdColG[1:q,1:q]), half_vec(pdColG[1:q,(q+1):p]))
  lambda2 <- (m2v==2)
  lambda2[lambda2] <- Inf
  
  # run SGL algorithm
  out.admm <- admm.pdglasso(S, lambda1 = lambda1,
                            lambda2 = lambda2,
                            type=c("V","I","A"),
                            eps.rel = eps.rel,
                            eps.abs = eps.abs,
                            max_iter= max_iter,
                            print.type=FALSE)
  K.hat <- NULL
  if(out.admm$internal.par$converged) K.hat = out.admm$X
  return(K.hat)
}




#' Check if a concentration matrix satisfies the likelihood equations.
#'
#' Checks if a matrix is the maximum likelihood estimate of the concentration
#' matrix of a pdRCON model represented by the colored graph for paired data
#' `pdColG`, computed from the covariance matrix `S`.
#'
#' @param K.mle candidate mle of \eqn{K} to be checked.
#' @param pdColG matrix representing a colored graph for paired data;
#' see [`pdglasso-package`] for details.
#' @param S a sample covariance matrix with the block structure describe
#' in [`pdglasso-package`].
#' @param toll threshold to check whether a value is equal to zero.
#' @param print.checks a logical (default `TRUE`).
#'
#' @return a logical equal to `TRUE` if `K.mle` satisfies the likelihood equations
#' and `FALSE` otherwise. If `print.checks = TRUE` the values used for the
#' checking are printed.
#'
#' @details
#' In order to check if `K.mle` is the maximum likelihood estimate, computed
#' from `S`, under the model represented by `pdColG` the following quantities
#' are considered:
#'
#' * the values of `K.mle` that are expected to be zero;
#'
#' * the differences between the entries of `K.mle` which are expected to be equal;
#'
#' * the differences between the relevant sums of the entries of `solve(K.mle)`
#'   and `S` which are expected to be equal, so as to satisfy the
#'   likelihood equations.
#'
#' * all the absolute quantities of the previous three points are divided
#'   by the \eqn{\ell_2} norm of the non-zero entries of `K.mle` so as to
#'   obtain relative quantities.
#'
#'   Then if all the above quantities are smaller than `toll` the check is
#'   considered successful and a `TRUE` is returned.
#' @keywords internal
#'
## @examples
## S <- var(toy_data$sample.data)
## K.hat <- pdRCON.mle(S, toy_data$pdColG)
## is.pdRCON.mle(K.hat, toy_data$pdColG, S)
##
is.pdRCON.mle <- function(K.mle, pdColG, S, toll=1e-8, print.checks=TRUE){
  S.mle <- solve(K.mle)
  G <- pdColG
  p <- dim(G)[1]
  q <- p/2
  l <- 1:q
  r <- (q+1):p
  #
  nz.ms <- sqrt(mean(K.mle[G!=0]^2))
  # missing edges
  Bz <- (G==0)
  #
  # uncoloured vertices and edges
  Bu <- (G==1)
  #
  # coloured vertices and inside-block coloured edges
  Bci <- (G[l,l]==2)
  #
  # across-block coloured edges
  Bca <- (G[l,r]==2)
  #
  # missing edges: check zero concentrations
  if(any(Bz)){
    am <- max(abs(K.mle[G==0]))
  }else{
    am <- 0
  }
  v  <- c()
  sK <- c()
  #
  # coloured vertices and inside-block coloured edges:
  # check equality constraints in K.mle
  # check likelihood equations in S.mle
  if(any(Bci)){
    sK[1] <- max(abs((K.mle[l,l]-K.mle[r,r])[Bci]))
    v[1]  <- max(abs((S.mle[l,l]+S.mle[r,r]-S[l,l]-S[r,r])[Bci]))
  }else{
    sK[1] <- 0
    v[1]  <- 0
  }
  ##
  # across-block coloured edges:
  # check equality constraints in K.mle
  # check likelihood equations in S.mle
  if(any(Bca)){
    sK[2] <- max(abs((K.mle[l,r]-K.mle[r,l])[Bca]))
    v[2]  <- max(abs((S.mle[l,r]+S.mle[r,l]-S[l,r]-S[r,l])[Bca]))
  }else{
    sK[2] <- 0
    v[2]  <- 0
  }
  #
  # uncoloured vertices and edges: check likelihood equations
  if(any(Bu)){
    v[3] <- max(abs((S.mle-S)[G==1]))
  }else{
    v[3] <- 0
  }
  if(print.checks){
    cat("\nABSOLUTE QUANTITIES:\n")
    cat("Max zero concentration           : ", am, "\n")
    cat("Max diff. of equal concentrations: ", max(sK), "\n")
    cat("Max error in likelihood equations: ", max(v), "\n")
    #
    cat("\nRELATIVE QUANTITIES:\n")
    cat("Max zero concentration           : ", am/nz.ms, "\n")
    cat("Max diff. of equal concentrations: ", max(sK)/nz.ms, "\n")
    cat("Max error in likelihood equations: ", max(v)/nz.ms, "\n\n")
  }
  max.all <- max(c(am, max(sK), max(v), am/nz.ms, max(sK)/nz.ms, max(v)/nz.ms))
  
  if(max.all < toll){
    check.result <- TRUE
  }else{
    check.result <- FALSE
  }
  return(check.result)
}




#' Compute the extended Bayesian Information Criterion (eBIC).
#'
#' This function computes the value of the eBIC for a given model and gamma value, for the purpose
#' of model selection (see Eq.1 of Feygel & Drton, 2010).
#'
#' @param S a sample covariance (or correlation) matrix with the block structure described in [`pdglasso-package`].
#' @param mod a list, the output object of a call to [`admm.pdglasso`].
#' @param n the sample size of the data used to compute the sample covariance matrix S.
#' @param gamma.eBIC a parameter governing the magnitude of the penalization term inside the criterion; it ranges from 0 to 1, where 0 makes the eBIC equivalent to BIC, and 0.5 being the suggested default value.
#' @param max_iter an integer; maximum number of iterations to be run in case
#'   the algorithm does not converge; passed to [`pdRCON.mle`].
#' @param mle a logical; if `TRUE`, uses the MLE (more accurate but slower computation);
#'   if `FALSE`, uses the pdglasso (ADMM) estimator (faster but potentially less precise).
#'   
#' @return A vector containing three elements:
#' * the value of the eBIC,
#' * the log-likelihood,
#' * and the number of parameters.
#' @export
#'
#' @references Foygel, R., Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. *Advances in neural information processing systems*, 23. \url{https://proceedings.neurips.cc/paper/2010/file/072b030ba126b2f4b2374f342be9ed44-Paper.pdf}
#' @examples
#' S <- cov(toy_data$sample.data)
#' mod <- admm.pdglasso(S, lambda1=1, lambda2=0.5)
#' compute.eBIC(S, mod, n=60, gamma.eBIC=0.5)

compute.eBIC <- function(S,mod,n,
                         gamma.eBIC=0.5,
                         max_iter=5000,
                         mle = TRUE){
  G <- pdColG.get(mod)
  
  if(mle){
    K <- pdRCON.mle(S,
                    pdColG = G$pdColG,
                    eps.rel = mod$internal.par$eps.rel,
                    eps.abs = mod$internal.par$eps.abs,
                    max_iter = mod$internal.par$max_iter)
  }else{
    K <- mod$X
  }
  
  
  if(is.null(K)){
    out.vec <- rep(NA,3)
  }else{
    S <- S*(n-1)/n
    p <- dim(S)[1]
    n.par <- G$n.par
    #log.lik <- log(det(K))-sum(S*K)
    log.lik <- as.numeric(determinant(K, logarithm=TRUE)$modulus) - sum(S*K)
    
    eBIC <- -n*log.lik+log(n)*n.par+4*n.par*gamma.eBIC*log(p)
    #### corrected with n/2 instead of 1/2 in front of the loglik, so -2*(n/2)*loglik=-n*loglik
    out.vec <- c(eBIC,log.lik,n.par)
  }
  names(out.vec) <- c("eBIC     ","  log-Likelihood  ","num. of params")
  return(out.vec)
}




######### Secondary Functions

#' Create the model acronym.
#'
#' Creates the model acronym from the arguments `types` and `force.symm`.
#'
#' @param type a character vector.
#' @param force.symm either a character vector or `NULL`.
#' @param print.type logical (default `TRUE`) indicating whether the model details should be printed.
#'
#' @return A list with two character strings named `acronym.of.type` and `acronym.of.force`, the latter is  `NULL` if `force.symm` is   `NULL`.
#' @noRd

make.acronyms <- function(type, force.symm, print.type=TRUE){
  # internal function which actually makes the acronym
  make.a <- function(opt.str){
    opt.str <- tolower(opt.str)
    choice  <- match.arg(opt.str, c("vertex", "inside.block.edge", "across.block.edge") , several.ok = TRUE)
    acronym <- ""
    if (any(choice=="vertex"))              acronym <- paste(acronym, "V", sep="")
    if (any(choice=="inside.block.edge"))   acronym <- paste(acronym, "I", sep="")
    if (any(choice=="across.block.edge"))   acronym <- paste(acronym, "A", sep="")
    choice.print <- sort(toupper(unique(choice)), decreasing = TRUE)
    choice.print <- paste(choice.print, collapse = ", ")
    choice.print <- paste("[", choice.print, "]", sep="")
    return(list(acronym=acronym, choice.print=choice.print, choice=choice))
  }
  #
  acr.type <- make.a(type)
  #
  if(!is.null(force.symm)){
    acr.force     <- make.a(force.symm)
    is.contained  <- acr.force$choice %in% acr.type$choice
    if(!all(is.contained)){
      cat("\n")
      warning("'force.symm' must be a subvector of 'type'.\nSome entries of 'force.symm' not contained in 'type' have been ignored. ", call.=FALSE, immediate. = TRUE)
      if(any(is.contained)){
        force.symm <- acr.force$choice[is.contained]
        acr.force <- make.a(force.symm)
      }else{
        force.symm = NULL
      }
    }
  }
  if(is.null(force.symm)) acr.force <- list("acronym"=NULL, "choice.print"="NONE")
  if(print.type){
    cat("\nCall:\nColoured GGM for paired data with:\nallowed types of coloured symmetry = ",  acr.type$choice.print, "\n", sep="")
    cat("forced coloured symmetry = ",  acr.force$choice.print, "\n\n", sep="")
  }
  return(list(acronym.of.type=acr.type$acronym, acronym.of.force=acr.force$acronym))
}

#' Computes the number of constraints, i.e. number of rows of the matrix F.
#'
#' @param q p/2
#' @param acr.type type of acronym.
#'
#' @return .
#' @noRd

get.n.row.F <- function(q, acr.type){
  dim.hs <- q*(q-1)/2
  switch(acr.type,
         V   = q,
         I   = dim.hs,
         A   = dim.hs,
         VI  = q+dim.hs,
         VA  = q+dim.hs,
         IA  = 2*dim.hs,
         VIA = q+2*dim.hs)
}




#' Imposes lambda values such that symmetry is forced.
#'
#' @param p an even number.
#' @param lambda2 a scalar.
#' @param acr.type type of acronym.
#' @param acr.force type of acronym for forced symmetries.
#'
#' @return .
#' @noRd

lambda2.force.symm <- function(p, lambda2, acr.type, acr.force){
  acronym <- paste("t", acr.type, "_f", acr.force, sep="")
  q       <- p/2
  dim.hs  <- q*(q-1)/2
  lambda.max <- Inf
  vmax.q  <- rep(lambda.max, q)
  vmax.hs <- rep(lambda.max, dim.hs)
  v2.q    <- rep(lambda2, q)
  v2.hs   <- rep(lambda2, dim.hs)
  switch(acronym,
         tV_fV = lambda.max,
         tI_fI = lambda.max,
         tA_fA = lambda.max,
         #
         tVI_fV  = c(vmax.q, v2.hs),
         tVI_fI  = c(v2.q, vmax.hs),
         tVI_fVI = lambda.max,
         #
         tVA_fV  = c(vmax.q, v2.hs),
         tVA_fA  = c(v2.q, vmax.hs),
         tVA_fVA = lambda.max,
         #
         tIA_fI  = c(vmax.hs, v2.hs),
         tIA_fA  = c(v2.hs, vmax.hs),
         tIA_fIA = lambda.max,
         #
         tVIA_fV   = c(vmax.q, v2.hs, v2.hs),
         tVIA_fI   = c(v2.q, vmax.hs, v2.hs),
         tVIA_fA   = c(v2.q, v2.hs, vmax.hs),
         tVIA_fVI  = c(vmax.q, vmax.hs, v2.hs),
         tVIA_fVA  = c(vmax.q, v2.hs, vmax.hs),
         tVIA_fIA  = c(v2.q, vmax.hs, vmax.hs),
         tVIA_fVIA = lambda.max,
  )
}










