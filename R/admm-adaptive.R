require(compiler)
setCompilerOptions(optimize=3)
### WRAPPER FUNCTION for fitting
### creates grids of lambdas, fits model
### selects according to eBIC






#' Estimates a concentration matrix under the pdglasso model using adaptive ADMM
#' algorithm.
#'
#' Description here.
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
#'
#' @return AAA
#' @export
#'
#' @examples
#' !!! Create fake dataset
#' S <- cov(toy.data)
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
                     verbose     = FALSE) {
  #
  time.start    <- Sys.time()
  #
  # initializations
  #
  out.make.a <- make.acronims(type, force.symm)
  acr.type <- out.make.a$acronim.of.type
  acr.force <- out.make.a$acronim.of.force
  p <- dim(S)[1]
  q <- p/2
  if(!is.null(acr.force)) lambda2 <- lambda2.force.symm(p, lambda2, acr.type, acr.force, S)
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
  acronims=list(acronim.of.type=acr.type, acronim.of.force=acr.force)
  return(list(X=X, acronims=acronims, internal.par=internal.par))
}
admm.pdglasso_C<-cmpfun(admm.pdglasso)

# Inner ADMM loop for the main function
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


### Secondary Functions

# Produces acronym for the model
#
make.acronims <- function(type, force.symm){
  # internal function which actually makes the acronim
  make.a <- function(opt.str){
    opt.str <- tolower(opt.str)
    choice  <- match.arg(opt.str, c("vertex", "inside.block.edge", "across.block.edge") , several.ok = TRUE)
    acronim <- ""
    if (any(choice=="vertex"))              acronim <- paste(acronim, "V", sep="")
    if (any(choice=="inside.block.edge"))   acronim <- paste(acronim, "I", sep="")
    if (any(choice=="across.block.edge"))   acronim <- paste(acronim, "A", sep="")
    choice.print <- sort(toupper(unique(choice)), decreasing = TRUE)
    choice.print <- paste(choice.print, collapse = ", ")
    choice.print <- paste("(", choice.print, ")", sep="")
    return(list(acronim=acronim, choice.print=choice.print, choice=choice))
  }
  #
  acr.type <- make.a(type)
  cat("\nCall:\nRCON model for paired data with constraints\ntype of symmetry = ",  acr.type$choice.print, "\n", sep="")
  #
  if(!is.null(force.symm)){
    acr.force     <- make.a(force.symm)
    is.contained  <- acr.force$choice %in% acr.type$choice
    if(!all(is.contained)){
      warning("'force.symm' must be a subvector of 'type'.\nSome entries of 'force.symm' not contained in 'type' have been ignored. ", call.=FALSE, immediate. = TRUE)
      if(any(is.contained)){
        force.symm <- acr.force$choice[is.contained]
        acr.force <- make.a(force.symm)
      }else{
        force.symm = NULL
      }
    }
  }
  if(is.null(force.symm)) acr.force <- list("acronim"=NULL, "choice.print"="NONE")
  cat("forced  symmetry = ",  acr.force$choice.print, "\n\n", sep="")
  return(list(acronim.of.type=acr.type$acronim, acronim.of.force=acr.force$acronim))
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

lambda2.force.symm <- function(p, lambda2, acr.type, acr.force, S){
  acronim <- paste("t", acr.type, "_f", acr.force, sep="")
  q       <- p/2
  dim.hs  <- q*(q-1)/2
  lambda.max <- ceiling(max(abs(S))*10)
  vmax.q  <- rep(lambda.max, q)
  vmax.hs <- rep(lambda.max, dim.hs)
  v2.q    <- rep(lambda2, q)
  v2.hs   <- rep(lambda2, dim.hs)
  switch(acronim,
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



### Half-vectorization operator,  without the diagonal
# input is a symmetric matrix M
half.vec <- function(M){
  return(M[upper.tri(M, diag=FALSE)])
}

### Inverse of Half-vectorization operator, without the diagonal
# input is an m vector; returns a triangular (upper) matrix
# for the strict==T option, the diagonal is filled with zero
inv_half.vec<- function(m){
  mat.dim <- (1+(1+8*length(m))^0.5)/2
  M <- matrix(0, mat.dim, mat.dim)
  M[upper.tri(M, diag=FALSE)] <- m
  return(M)
}

# Transforms a symmetric matrix into a vector applying the
# v() operator defined in Section ??? of [2]. This function
# is called by admm.inner()
#
mat2vec <- function(M){
  p  <- dim(M)[1]
  q  <- p/2
  return(c(diag(M[1:q,1:q]),
           diag(M[(q+1):p,(q+1):p]),
           half.vec(M[1:q,1:q]),
           half.vec(M[(q+1):p,(q+1):p]),
           half.vec(M[1:q,(q+1):p]),
           half.vec(M[(q+1):p,1:q]),
           diag(M[1:q,(q+1):p]))
         )
}

# Inverse operation with respect to mat2vec().
# This function is called by admm.inner()
#
vec2mat <- function(m){
  p <- (-1+sqrt(1+8*length(m)))/2
  q <- p/2
  dim.hs <- q*(q-1)/2
  i1 <- 0
  i2 <- i1+q
  i3 <- i2+q
  i4 <- i3+dim.hs
  i5 <- i4+dim.hs
  i6 <- i5+dim.hs
  i7 <- i6+dim.hs
  i8 <- i7+q
  M <- matrix(0,p,p)
  M[1:q,1:q]               <- inv_half.vec(m[(i3+1):i4])
  diag(M[1:q,1:q])         <- m[(i1+1):i2]
  M[(q+1):p,(q+1):p]       <- inv_half.vec(m[(i4+1):i5])
  diag(M[(q+1):p,(q+1):p]) <- m[(i2+1):i3]
  M[1:q,(q+1):p]           <- inv_half.vec(m[(i5+1):i6])
  M[1:q,(q+1):p]           <- M[1:q,(q+1):p] + t(inv_half.vec(m[(i6+1):i7]))
  diag(M[1:q,(q+1):p])     <- m[(i7+1):i8]
  M[lower.tri(M, diag=FALSE)] <- t(M)[lower.tri(M, diag=FALSE)]
  return(M)
}


