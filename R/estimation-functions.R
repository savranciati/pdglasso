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
#' @param n the sample size of the data used to compute the sample covariance matrix S.
#' @param lams a 2x3 matrix; first row refers to `lambda1` and second row to
#'   `lambda2`; for each row, values are (i) minimum value for the grid; (ii)
#'   maximum value for the grid; (iii) number of points in the grid; if `NULL`
#'   defaulta values are used, i.e. from max.lams to max.lams/10, and 10 grid
#'   points.
#' @param gamma.eBIC the parameter for the eBIC computation. gamma=0 is equivalent to BIC.
#' @param progress a logical value; if `TRUE` provides a visual update in the console about the grid search over `lambda1` and `lambda2`
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
#' sel.mod <- pdRCON.fit(S,n=60)
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
                       print.type  = TRUE){
  start.time <- Sys.time()
  ## Max values for lambda_1 and lambda_2 according to theorems; only needed for the grid search

  if(is.null(lams)){
    lams <- matrix(0,2,3)
    lams[,2] <- lams.max(S)
    lams[,1] <- lams[,2]/10
    lams[,3] <- c(10,10)
  }
  rownames(lams) <- c("l1","l2")
  colnames(lams) <- c("min","max","n.pts")

  ## Prepare temp objects
  eBIC.l1 <-  matrix(0,lams[1,3],4)
  eBIC.l2 <-  matrix(0,lams[2,3],4)

  ### First grid search for lambda_1, with lambda_2=0
  # l1.vec <- exp(seq(min.ls[1],max.ls[1], length.out=n.l1))
  l1.vec <- seq(lams[1,1], lams[1,2], length.out=lams[1,3])
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
                                   max_iter=max_iter)
    eBIC.l1[i,4] <- mod.out$internal.par$converged+0
    if(eBIC.l1[i,4]==0) cat("Convergence not achieved for this value of lambda1! \n")
  }
  best.l1 <- l1.vec[which.min(eBIC.l1[,1])]
  if(length(best.l1)==0) stop("Grid search of lambda1 failed!")

  if(progress==TRUE) cat("--- \n", sep="")

  ### Second grid search for lambda_2, with lambda_1=best.l1
  l2.vec <- seq(lams[2,1],lams[2,2], length.out=lams[2,3])
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
                                   max_iter=max_iter)
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
#'
#' @export
#'
#' @examples
#'
#' S <- cov(toy_data$sample.data)
#' mod.out <- pdRCON.fit(S,n=60)
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


#' ADMM graphical lasso algorithm for coloured GGMs for paired data.
#'
#' By providing a covariance matrix `S` and values for `lambda1` and `lambda2`,
#' this function estimates a concentration matrix X within the pdRCON
#' submodel class, identified by the arguments `type` and `force.symm`, based on the pdglasso method (Ranciati & Roverato, 2023) using an
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
#' @param verbose a logical; if `TRUE` the progress (and internal convergence of
#'   inner loop) is shown in the console while the algorithm is running.
#' @param print.type a logical; if `TRUE` the pdRCON submodel class considered, as
#'   specified by the arguments `type` and `force.symm` - is returned as printed
#'   output in the console.
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
#' @references Ranciati, S., Roverato, A., (2023). On the application of Gaussian graphical models to paired data problems. *arXiv pre-print*. \url{https://arxiv.org/abs/2307.14160}
#' @examples
#'
#' S <- cov(toy_data$sample.data)
#' admm.pdglasso(S)
admm.pdglasso <- function(S,
                     lambda1     = 1,
                     lambda2     = 0.0001,
                     type        = c("vertex", "inside.block.edge", "across.block.edge"),
                     force.symm  = NULL,
                     X.init      = NULL,
                     rho1        = 1,
                     rho2        = 1,
                     varying.rho1= TRUE,
                     varying.rho2= TRUE,
                     max_iter    = 5000,
                     eps.abs     = 1e-6,
                     eps.rel     = 1e-6,
                     verbose     = FALSE,
                     print.type  = TRUE) {
  #
  time.start    <- Sys.time()
  #
  # initializations
  #
  converged <- TRUE
  out.make.a <- make.acronyms(type, force.symm, print.type=print.type)
  acr.type <- out.make.a$acronym.of.type
  acr.force <- out.make.a$acronym.of.force
  p <- dim(S)[1]
  q <- p/2
  if(!is.null(acr.force)) lambda2 <- lambda2.force.symm(p, lambda2, acr.type, acr.force)
  n.row.F <- get.n.row.F(q, acr.type)

  if(is.null(X.init)) X <-  diag(1,p) else X <- X.init
  U <- Z <- matrix(0,p,p)
  #
  # required for the computation of residuals,
  # in this case is is not necessary to use mat2vec
  z.v <- Z[upper.tri(Z, diag = TRUE)]
  #
  # set parameters for varying.rho1
  # (suggested tau.inc <- tau.dec <- 2, mu <- 10)
  #
  mu      <- 10
  tau.inc <-  2
  tau.dec <-  2
  n.iter.rho1_update_last <- 1
  #
  # Outer ADMM
  #
  for (k in 1:max_iter){
    #
    # Update X
    #
    #
    decomp  <- eigen(rho1*(Z-U)-S, symmetric = TRUE)
    Q       <- decomp$vectors
    X_tilde <- diag((decomp$values+(decomp$values^2+4*rho1)^0.5)/(rho1*2))
    X       <- Q%*%X_tilde%*%t(Q)
    #
    # Update Z - inner ADMM
    #
    Z <- admm.inner(X=X, U=U, rho1=rho1, lambda1=lambda1, lambda2=lambda2, rho2=rho2,
                       verbose_int=verbose, varying.rho2=varying.rho2, n.row.F=n.row.F,
                       acr.type=acr.type, eps.abs=eps.abs, eps.rel=eps.rel)
    #
    # Update U
    #
    U <- U+X-Z
    #
    # compute norm of residuals
    #
    last.z.v <- z.v
    z.v <- Z[upper.tri(Z, diag = TRUE)]
    x.v <- X[upper.tri(X, diag = TRUE)]
    r.kk <- sqrt(sum((x.v-z.v)^2))
    s.kk <- rho1*sqrt(sum((z.v-last.z.v)^2))
    #
    #
    # compute epsilon values
    #
    u.v <- U[upper.tri(U, diag = TRUE)]
    eps.pri  <- sqrt(p*(p+1)/2)*eps.abs+eps.rel*max(sqrt(sum(x.v^2)), sqrt(sum(z.v^2)))
    eps.dual <- sqrt(p*(p+1)/2)*eps.abs+eps.rel*sqrt(sum((rho1*u.v)^2))
    #
    # Check convergence
    #
    if(r.kk<eps.pri & s.kk<eps.dual){
      #cat("\nConverged:\niteration number: ",k,"\nr.kk=",r.kk, "\ns.kk=", s.kk, "\nrho1= ", rho1, "\n")
      break
    }
    #
    # varying.rho1 procedure
    #
    scale.factor <- rho1# in Boyd et al., scale.factor=1
    #
    if(varying.rho1){
      if(r.kk>(mu*s.kk/scale.factor)){
        rho1 <- rho1*tau.inc
        U <- U/tau.inc
        n.iter.rho1_update_last <- k
      }
      if((s.kk/scale.factor)>(mu*r.kk)){
        rho1 <- rho1/tau.dec
        U <- U*tau.dec
        n.iter.rho1_update_last <- k
      }
    }

  }
  if(k==max_iter){
    converged <- FALSE
    #warning(paste("Convergence not achieved; iterations performed: ",max_iter,".", sep=""))
  }
  time.diff <- Sys.time()-time.start
  internal.par <- list(execution.time=time.diff, res.primal=r.kk, res.dual=s.kk,
                       lambda1 = lambda1, lambda2=unique(lambda2), n.iter=k,
                       n.iter.rho1_update_last=n.iter.rho1_update_last, last.rho1=rho1,
                       eps.primal=eps.pri, eps.dual=eps.dual, eps.abs=eps.abs, eps.rel=eps.rel,
                       max_iter = max_iter, converged = converged)
  acronyms=list(acronym.of.type=acr.type, acronym.of.force=acr.force)
  dimnames(X) <- dimnames(S)
  return(list(X=X, acronyms=acronyms, internal.par=internal.par))
}
#

#' Inner ADMM loop called by the main function [`admm.pdglasso`].
#'
#' This inner ADMM loop is called by the outer loop (and main function) [`admm.pdglasso`].
#' It inherits most of the parameters, together with the two quantities X and U which are
#' related to the steps of the ADMM.
#'
#'
#' @inheritParams admm.pdglasso
#' @param X the concentration matrix at the current iteration (in matrix form).
#' @param U the scaled dual variable of the optimization problem (in matrix form).
#' @param max_iter_int maximum number of iterations for this inner loop.
#' @param verbose_int a `TRUE/FALSE` to switch on/off a console message during the iterations.
#' @param n.row.F the number of rows of the constraint matrix F.
#' @param acr.type the acronym of the called model, summarizing which penalties to use.
#' @return Results of inner ADMM loop in matrix form.
#'
#' @noRd
admm.inner <- function(X,
                        U,
                        rho1,
                        lambda1,
                        lambda2,
                        rho2        = 1,
                        varying.rho2= TRUE,
                        max_iter_int= 5000,
                        eps.abs   = 1e-12,
                        eps.rel   = 1e-8,
                        verbose_int = FALSE,
                        n.row.F,
                        acr.type){
  #
  # initializations
  #
  p <- dim(X)[1]
  q <- p/2
  b <- mat2vec(X+U)
  d <- length(b)
  v <- t <- rep(0, n.row.F)
  lambda <- lambda2/rho1
  x <- rep(1,d)
  #
  # set parameters for varying rho2
  # (suggested in Boyd; tau.inc <- tau.dec <- 2, mu <- 10)
  #
  mu      <- 10
  tau.inc <-  2
  tau.dec <-  2
  #
  #
  for (kk in 1:max_iter_int){
    #
    # Update x
    #
    x <- rho2/(1+2*rho2)*tF.by.vec(v-t-F.by.vec(b, p, acr.type), p, acr.type)+b
    #
    #
    prod.Fx.tmp <- F.by.vec(x, p, acr.type)
    prod.Fx <- prod.Fx.tmp+t
    #
    # Update v
    #
    last.v <- v
    #
    v <- pmax(prod.Fx-lambda/rho2,0)-pmax(-prod.Fx-lambda/rho2,0)
    #
    # Update t
    #
    t <- prod.Fx-v
    #
    # compute norm of residuals
    #
    r.kk.f <- sqrt(sum((prod.Fx.tmp - v)^2))
    s.kk.f <- sqrt(sum((rho2*tF.by.vec(v-last.v, p, acr.type))^2))
    #
    # compute epsilon values
    #
    eps.pri  <- sqrt(n.row.F)*eps.abs+eps.rel*max(sqrt(sum(prod.Fx.tmp^2)), sqrt(sum(v^2)))
    eps.dual <- sqrt(p*(p+1)/2)*eps.abs+eps.rel* sqrt(sum((tF.by.vec(rho2*t, p, acr.type))^2))
    #
    # Check convergence
    #
    if(r.kk.f<eps.pri & s.kk.f<eps.dual){
      if(verbose_int==TRUE) cat(kk,".", sep="")
      break
    }
    if(varying.rho2){
      if(r.kk.f>(mu*s.kk.f)){
        rho2 <- rho2*tau.inc
        t <- t/tau.inc
      }
      if(s.kk.f>(mu*r.kk.f)){
        rho2 <- rho2/tau.dec
        t <- t*tau.dec
      }
    }
  }
  z.temp <- sign(x)*pmax(abs(x)-lambda1/rho1,0)
  #
  if(kk==max_iter_int){
    warning(paste("Internal loop convergence not achieved; iterations performed: ",max_iter_int,".", sep=""))
  }
  #
  return(vec2mat(z.temp))
}



#' Maximum likelihood estimation
#'
#'
#' Computes the maximum likelihood estimate of the concentration matrix of a pdRCON model.
#'
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
#'
#' S <- var(toy_data$sample.data)
#' K.hat <- pdRCON.mle(S, toy_data$pdColG)
#'
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
  m2v <- c(diag(pdColG[1:q,1:q]), half.vec(pdColG[1:q,1:q]), half.vec(pdColG[1:q,(q+1):p]))
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



#' Check if a concentration matrix satisfies the likelihood equations
#'
#' Checks if a matrix is the maximum likelihood estimate of the concentration
#' matrix of a pdRCON model represented by the colored graph for paired data
#' `pdColG`, computed from the covariance matrix `S`.
#'
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
#'
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
#' compute.eBIC(S,mod,n=60,gamma.eBIC=0.5)
compute.eBIC <- function(S,mod,n,
                         gamma.eBIC=0.5,
                         max_iter=5000){
  G <- pdColG.get(mod)
  K <- pdRCON.mle(S,
                  pdColG = G$pdColG,
                  eps.rel = mod$internal.par$eps.rel,
                  eps.abs = mod$internal.par$eps.abs,
                  max_iter = mod$internal.par$max_iter)
  if(is.null(K)){
  out.vec <- rep(NA,3)
  }else{
    S <- S*(n-1)/n
    p <- dim(S)[1]
    n.par <- G$n.par
    log.lik <- log(det(K))-sum(S*K)
    eBIC <- -n*log.lik+log(n)*n.par+4*n.par*gamma.eBIC*log(p)
    #### corrected with n/2 instead of 1/2 in front of the loglik, so -2*(n/2)*loglik=-n*loglik
    out.vec <- c(eBIC,log.lik,n.par)
  }
  names(out.vec) <- c("eBIC     ","  log-Likelihood  ","num. of params")
  return(out.vec)
}


######### Secondary Functions

#' Create the model acronym
#'
#' Creates the model acronym from the arguments `types` and `force.symm`.
#'
#'
#'
#' @param type a character vector.
#' @param force.symm either a character vector or `NULL`.
#' @param print.type logical (default `TRUE`) indicating whether the model details should be printed.
#'
#' @return A list with two character strings named `acronym.of.type` and `acronym.of.force`, the latter is  `NULL` if `force.symm` is   `NULL`.
#' @noRd
#'
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

#' Computes the number of constraints, i.e. number of rows of the matrix F
#'
#' @param q p/2
#' @param acr.type type of acronym.
#'
#' @return .
#' @noRd
#'
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



#' Computes F%*%v
#'
#' @param v a vector
#' @param p number of rows/columns of S
#' @param acr.type type of acronym
#'
#' @return .
#' @noRd
#'
F.by.vec <- function(v, p, acr.type){
  q <- p/2
  dim.hs <- q*(q-1)/2
  i1 <- 0
  i2 <- i1+q
  i3 <- i2+q
  i4 <- i3+dim.hs
  i5 <- i4+dim.hs
  i6 <- i5+dim.hs
  i7 <- i6+dim.hs
  v1  <- v[(i1+1):i2]
  v1p <- v[(i2+1):i3]
  v2  <- v[(i3+1):i4]
  v2p <- v[(i4+1):i5]
  v3  <- v[(i5+1):i6]
  v3p <- v[(i6+1):i7]
  switch(acr.type,
         V   = v1-v1p,
         I   = v2-v2p,
         A   = v3-v3p,
         VI  = c(v1, v2)-c(v1p, v2p),
         VA  = c(v1, v3)-c(v1p, v3p),
         IA  = c(v2, v3)-c(v2p, v3p),
         VIA = c(v1, v2, v3)-c(v1p, v2p, v3p))
}


#' Computes t(F)%*%v
#'
#' @param v a vector
#' @param p number of rows/columns of S
#' @param acr.type type of acronym
#'
#' @return .
#' @noRd
#'
tF.by.vec <- function(v, p, acr.type){
  q <- p/2
  dim.hs <- q*(q-1)/2
  dim.h  <- p*(p+1)/2
  i1 <- 0
  switch(acr.type,
         V   = { i2 <- i1+q
                 v1 <- v[(i1+1):i2]
                 return(c(v1, -v1, rep(0, dim.h-2*q)))
               },
         I   = { i2 <- i1+dim.hs
                 v2 <- v[(i1+1):i2]
                 return(c(rep(0, 2*q), v2, -v2,  rep(0, dim.h-dim.hs*2-2*q)))
               },
         A   = { i2 <- i1+dim.hs
                 v3 <- v[(i1+1):i2]
                 return(c(rep(0, 2*q+2*dim.hs), v3, -v3, rep(0, q)))
               },
         VI  = { i2 <- i1+q
                 i3 <- i2+dim.hs
                 v1 <- v[(i1+1):i2]
                 v2 <- v[(i2+1):i3]
                 return(c(v1, -v1, v2, -v2,  rep(0, dim.h-dim.hs*2-2*q)))
               },
         VA  = { i2 <- i1+q
                 i3 <- i2+dim.hs
                 v1 <- v[(i1+1):i2]
                 v3 <- v[(i2+1):i3]
                 return(c(v1, -v1, rep(0, 2*dim.hs),v3, -v3, rep(0, q)))
               },
         IA  = { i2 <- i1+dim.hs
                 i3 <- i2+dim.hs
                 v2  <- v[(i1+1):i2]
                 v3  <- v[(i2+1):i3]
                 return(c(rep(0, 2*q), v2, -v2, v3, -v3, rep(0, q)))
         },
         VIA = { i2 <- i1+q
                 i3 <- i2+dim.hs
                 i4 <- i3+dim.hs
                 v1  <- v[(i1+1):i2]
                 v2  <- v[(i2+1):i3]
                 v3 <- v[(i3+1):i4]
                 return(c(v1, -v1, v2, -v2, v3, -v3, rep(0, q)))
               })
}

#' Imposes lambda values such that symmetry is forced
#'
#' @param p .
#' @param lambda2 .
#' @param acr.type type of acronym
#' @param acr.force type of acronym for forced symmetries
#'
#' @return .
#' @noRd
#'
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
