require(compiler)
setCompilerOptions(optimize=3)

#' Fit and select a coloured graphical models for paired data according to eBIC
#' criterion.
#'
#' Performs a sequence of calls to [admm.pdglasso()] providing two grids of
#' values for lambda_1 and lambda_2. First, a grid search conditional on
#' lambda_2=0 is run to select the best lambda_1 value among the candidates
#' (according to eBIC); conditional on the best lambda_1, a similar search is
#' performed for lambda_2. The output is the select model, given by the
#' estimated concenration matrix and corresponding graph.
#'
#' @inheritParams admm.pdglasso
#' @param n the sample size of the data used to compute the sample covariance matrix S.
#' @param n.l1 the number of values in the grid of candidates for lambda_1.
#' @param n.l2 the number of values in the grid of candidates for lambda_2.
#' @param gamma the parameter for the eBIC computation. gamma=0 is equivalent to BIC.
#'
#' @return a list:
#' * model, the final model;
#' * pdColG, the associated coloured graph;
#' * l1.path, a matrix containing the grid values for lambda_1 as well as quantities used in eBIC computation;
#' * l2.path, a matrix containing the grid values for lambda_2 as well as quantities used in eBIC computation.
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' fit.pdColG(S,n=60)
fit.pdColG <- function(S,
                       n,
                       n.l1        = 15,
                       n.l2        = 15,
                       gamma.eBIC  = 0.5,
                       type        = c("vertex", "inside.block.edge", "across.block.edge"),
                       force.symm  = NULL,
                       X.init      = NULL,
                       rho1        = 1,
                       rho2        = 1,
                       varying.rho1= TRUE,
                       varying.rho2= TRUE,
                       max_iter    = 1000,
                       eps.abs     = 1e-12,
                       eps.rel     = 1e-12,
                       verbose     = FALSE,
                       print.type  = TRUE,
                       ...){

  ## Max values for lambda_1 and lambda_2 according to theorems; only needed for the grid search
  max.ls <- max.lams(S)

  ## Prepare temp objects
  eBIC.l1 <- eBIC.l2 <-  matrix(0,n.l1,3)

  ### First grid search for lambda_1, with lambda_2=0
  l1.vec <- exp(seq(log(min(abs(S))),log(max.ls[1]), length.out=n.l1-1))
  l1.vec <- c(0, l1.vec)
  for(i in 1:n.l1){
    mod.out <- admm.pdglasso(S,
                             lambda1=l1.vec[i],
                             lambda2=0,
                             print.type=FALSE,
                             ...)
    eBIC.l1[i,] <- compute.eBIC(S, n, mod.out, gamma.eBIC=gamma.eBIC)
  }
  best.l1 <- l1.vec[which.min(eBIC.l1[,1])]

  ### Second grid search for lambda_2, with lambda_1=best.l1
  l2.vec <- exp(seq(log(min(abs(S))),log(max.ls[2]), length.out=n.l1-1))
  l2.vec <- c(0, l2.vec)
  for(i in 1:n.l2){
    mod.out <- admm.pdglasso(S,
                             lambda1=best.l1,
                             lambda2=l2.vec[i],
                             print.type=FALSE,
                             ...)
    eBIC.l2[i,] <- compute.eBIC(S, n, mod.out, gamma.eBIC=gamma.eBIC)
  }
  best.l2 <- l2.vec[which.min(eBIC.l2[,1])]

  ## Fit final model
  mod.out <- admm.pdglasso(S,
                      lambda1=best.l1,
                      lambda2=best.l2)

  G <- get.pdColG(mod.out)

  l1.path=cbind(l1.vec,eBIC.l1)
  l2.path=cbind(l2.vec,eBIC.l2)
  names(l1.path) <- c("lambda1.grid", "BIC     ","  log-Likelihood  ","DF (estimated.)")
  names(l2.path) <- c("lambda2.grid", "BIC     ","  log-Likelihood  ","DF (estimated.)")

  return(list(model=mod.out,
              pdColG=G,
              l1.path=l1.path,
              l2.path=l2.path))

}


#' Estimate a concentration matrix under the pdColG model using (adaptive) ADMM
#' graphical lasso algorithm.
#'
#' By providing a covariance matrix S and values for lambda_1 and lambda_2, this
#' function estimates a concentration matrix X under the coloured graphical
#' model for paired data, using the (adaptive) ADMM algorithm. The output is the
#' matrix and a list of internal parameters used by the function, together with
#' the specific call in terms of symmetries and penalties required by the user.
#'
#' @param S A \eqn{p \times p} covariance (or correlation) matrix.
#' @param lambda1 A non-negative scalar (or vector) penalty that encourages
#'   sparsity in the concentration matrix. If a vector is provided, it should
#'   match the appropriate length, i.e.
#' @param lambda2 A non-negative scalar (or vector) penalty that encourages
#'   equality constraints in the concentration matrix. If a vector is provided,
#'   it should match the appropriate length, i.e.
#' @param type A string or vector of strings for the type of equality
#'   constraints to be imposed; zero, one or more available options can be
#'   selected among: * "vertex", symmetries are imposed on the diagonal entries
#'   of the concentration matrix. * "inside.block.edge", symmetries are imposed
#'   between elements of the LL and RR block the concentration matrix. *
#'   "across.block.edge", symmetries are imposed between elements of the LR and
#'   RL block the concentration matrix. Shortened forms are accepted too, i.e.
#'   "V" or "vert" for "vertex".
#' @param force.symm  A string or vector of strings to impose forced symmetry on
#'   the corresponding block of the concentration matrix. Same options as
#'   "type".
#' @param X.init (optional) A \eqn{p \times p} initial guess for the
#'   concentration matrix and/or starting solution for the ADMM algorithm.
#' @param rho1 A scalar; tuning parameter of the ADMM algorithm to be used for
#'   the outer loop. It must be strictly positive.
#' @param rho2 A scalar; tuning parameter of the ADMM algorithm to be used for
#'   the inner loop. It must be strictly positive.
#' @param varying.rho1 A boolean value; if `TRUE` the parameter rho1 is updated
#'   iteratively to speed-up convergence.
#' @param varying.rho2 A boolean value; if `TRUE` the parameter rho2 is updated
#'   iteratively to speed-up convergence.
#' @param max_iter An integer; maximum number of iterations to be run in case
#'   the algorithm does not converge.
#' @param eps.abs A scalar; the absolute precision required for the computation
#'   of primal and dual residuals of the ADMM algorithm.
#' @param eps.rel A scalar; the relative precision required for the computation
#'   of primal and dual residuals of the ADMM algorithm.
#' @param verbose A boolean value; if `TRUE` the progress (and internal
#'   convergence of inner loop) is shown in the console while the algorithm is
#'   running.
#' @param print.type A boolean value; if `TRUE` the acronym used for the model -
#'   which penalties - is returned as printed output in the console.
#'
#' @return A list, whose element are:
#' * `X`, the estimated concentration matrix
#'   under the pdglasso model; the model is identified by the values of lambda1
#'   and lambda 2, together with the type of penalization imposed.
#' * `acronyms`, a vector of strings for the type of penalties and forced symmetries imposed
#'   when calling the function.
#' * `internal.par`, a list of internal parameters
#'   passed to the function at the call, as well as convergence information.
#' @export
#'
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
                     max_iter    = 1000,
                     eps.abs     = 1e-12,
                     eps.rel     = 1e-12,
                     verbose     = FALSE,
                     print.type  = TRUE) {
  #
  time.start    <- Sys.time()
  #
  # initializations
  #
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
  for (k in 1:max_iter) {
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
    Z <- admm.inner_C(X=X, U=U, rho1=rho1, lambda1=lambda1, lambda2=lambda2, rho2=rho2,
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
      cat("\nConverged:\niteration number: ",k,"\nr.kk=",r.kk, "\ns.kk=", s.kk, "\nrho1= ", rho1, "\n")
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
  time.diff <- Sys.time()-time.start
  internal.par <- list(execution.time=time.diff, res.primal=r.kk, res.dual=s.kk,
                       lambda1 = lambda1, lambda2=unique(lambda2), n.iter=k,
                       n.iter.rho1_update_last=n.iter.rho1_update_last, last.rho1=rho1,
                       eps.primal=eps.pri, eps.dual=eps.dual, eps.abs=eps.abs, eps.rel=eps.rel)
  acronyms=list(acronym.of.type=acr.type, acronym.of.force=acr.force)
  return(list(X=X, acronyms=acronyms, internal.par=internal.par))
}
admm.pdglasso_C<-cmpfun(admm.pdglasso)

#

#' Inner ADMM loop called by the main function [admm.pdglasso()].
#'
#' This inner ADMM loop is called by the outer loop (and main function) [admm.pdglasso()].
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
                        max_iter_int= 1000,
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
  if(kk==max_iter_int) cat("\n Internal: not converged \n")
  #
  return(vec2mat(z.temp))
}
admm.inner_C<-cmpfun(admm.inner)



#' Maximum likelihood estimate
#'
#' Computes the m.l.e. of the concentration matrix of a colured graphical model
#' for paired data.
#'
#' @param S a sample variance and covariance matrix.
#' @param pdColG a coloured graph for paired data.
#'
#' @return the m.l.e. of the concentration matrix \eqn{\Sigma^{-1}}.
#' @export
#'
#' @examples
#' #
#'
pdColG.mle <- function(S, pdColG){

  # make vector lambda1
  lambda1 <- (mat2vec(pdColG)==0)
  lambda1[lambda1] <- Inf

  # make vector lambda 2
  p <- nrow(pdColG)
  q <- p/2
  m2v <- c(diag(pdColG[1:q,1:q]), half.vec(pdColG[1:q,1:q]), half.vec(pdColG[1:q,(q+1):p]))
  lambda2 <- (m2v==2)
  lambda2[lambda2] <- Inf

  # run SGL algorithm
  K.hat <- admm.pdglasso(S, lambda1 = lambda1, lambda2 = lambda2, print.type=FALSE)$X

  return(K.hat)
}


#' Title
#'
#' @param S
#' @param mod
#' @param n
#' @param gamma
#' @param max_iter
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
compute.eBIC <- function(S,mod,n,
                         gamma.eBIC=0.5){
  G <- get.pdColG(mod)
  K <- pdColG.mle(S,G)
  S <- S*(n-1)/n
  p <- dim(S)[1]
  dof <- G$dof
  log.lik <- log(det(K))-sum(S*K)
  #### corrected with n/2 instead of 1/2 in front of the loglik, so -2*(n/2)*loglik=-n*loglik
  out.list <- c(-n*log.lik+log(n)*dof+4*dof*gamma.eBIC*log(p),log.lik,dof)
  names(out.list) <- c("BIC     ","  log-Likelihood  ","DF (estimated.)")
  return(out.list)
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
    cat("\nCall:\nColoured graph for paired data with:\nallowed type of coloured symmetry = ",  acr.type$choice.print, "\n", sep="")
    cat("forced coloured symmetry = ",  acr.force$choice.print, "\n\n", sep="")
  }
  return(list(acronym.of.type=acr.type$acronym, acronym.of.force=acr.force$acronym))
}



# Computes the number of constraints, i.e. number of rows of the matrix F
#
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

# # Computes F%*%v
#
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


# Computes t(F)%*%v
#

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

# Imposes lambda values such that symmetry is forced

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
