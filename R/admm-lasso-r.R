#' Half-vectorization operator, without the diagonal.
#'
#' @param M a symmetric matrix.
#'
#' @return A vector.
#' @noRd
#'
half_vec_r <- function(M){
  return(M[upper.tri(M, diag=FALSE)])
}




#' Inverse of Half-vectorization operator, without the diagonal.
#'
#' @param m a vector.
#'
#' @return A triangular (upper) matrix.
#' @noRd
#'
inv_half_vec_r<- function(m){
  mat.dim <- (1+(1+8*length(m))^0.5)/2
  M <- matrix(0, mat.dim, mat.dim)
  M[upper.tri(M, diag=FALSE)] <- m
  return(M)
}




#' Transforms a symmetric matrix into a vector applying the
# v() operator.
#'
#' @param M a matrix.
#'
#' @return A vector.
#' @noRd
#'
mat2vec_r <- function(M){
  p  <- dim(M)[1]
  q  <- p/2
  return(c(diag(M[1:q,1:q]),
           diag(M[(q+1):p,(q+1):p]),
           half_vec_r(M[1:q,1:q]),
           half_vec_r(M[(q+1):p,(q+1):p]),
           half_vec_r(M[1:q,(q+1):p]),
           half_vec_r(M[(q+1):p,1:q]),
           diag(M[1:q,(q+1):p]))
  )
}




#' Inverse operation with respect to mat2vec_r().
#'
#' @param m a vector.
#'
#' @return A matrix.
#' @noRd
#'
vec2mat_r <- function(m){
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
  M[1:q,1:q]               <- inv_half_vec_r(m[(i3+1):i4])
  diag(M[1:q,1:q])         <- m[(i1+1):i2]
  M[(q+1):p,(q+1):p]       <- inv_half_vec_r(m[(i4+1):i5])
  diag(M[(q+1):p,(q+1):p]) <- m[(i2+1):i3]
  M[1:q,(q+1):p]           <- inv_half_vec_r(m[(i5+1):i6])
  M[1:q,(q+1):p]           <- M[1:q,(q+1):p] + t(inv_half_vec_r(m[(i6+1):i7]))
  diag(M[1:q,(q+1):p])     <- m[(i7+1):i8]
  M[lower.tri(M, diag=FALSE)] <- t(M)[lower.tri(M, diag=FALSE)]
  return(M)
}




#' Computes F%*%v.
#'
#' @param v a vector.
#' @param p number of rows/columns of S.
#' @param acr.type type of acronym.
#'
#' @return .
#' @noRd
#'
F_by_vec_r <- function(v, p, acr.type){
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




#' Computes t(F)%*%v.
#'
#' @param v a vector.
#' @param p number of rows/columns of S.
#' @param acr.type type of acronym.
#'
#' @return .
#' @noRd
#'
tF_by_vec_r <- function(v, p, acr.type){
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




#' R version, inner loop: ADMM graphical lasso algorithm for coloured GGMs for paired data.
#'
#' @param X a matrix.
#' @param U a matrix.
#' @param rho1 a positive scalar.
#' @param lambda1 a non-negative scalar (or vector).
#' @param lambda2 a non-negative scalar (or vector).
#' @param rho2 a positive scalar.
#' @param varying.rho2 a logical.
#' @param max_iter_int an integer.
#' @param eps.abs a scalar.
#' @param eps.rel a scalar.
#' @param verbose_int a logical.
#' @param n.row.F a scalar.
#' @param acr.type a string.
#'
#' @return A list.
#' @noRd

admm_inner_r <- function(X,
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
  b <- mat2vec_r(X+U)
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
    x <- rho2/(1+2*rho2)*tF_by_vec_r(v-t-F_by_vec_r(b, p, acr.type), p, acr.type)+b
    #
    #
    prod.Fx.tmp <- F_by_vec_r(x, p, acr.type)
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
    s.kk.f <- sqrt(sum((rho2*tF_by_vec_r(v-last.v, p, acr.type))^2))
    #
    # compute epsilon values
    #
    eps.pri  <- sqrt(n.row.F)*eps.abs+eps.rel*max(sqrt(sum(prod.Fx.tmp^2)), sqrt(sum(v^2)))
    eps.dual <- sqrt(p*(p+1)/2)*eps.abs+eps.rel* sqrt(sum((tF_by_vec_r(rho2*t, p, acr.type))^2))
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
  return(vec2mat_r(z.temp))
}




#' R version: ADMM graphical lasso algorithm for coloured GGMs for paired data.
#' @param S a matrix.
#' @param lambda1 a non-negative scalar (or vector).
#' @param lambda2 a non-negative scalar (or vector).
#' @param type,force.symm two subvectors of `c("vertex", "inside.block.edge",
#'   "across.block.edge")`.
#' @param X.init a matrix.
#' @param rho1 a positive scalar.
#' @param rho2 a positive scalar.
#' @param varying.rho1 a logical.
#' @param varying.rho2 a logical.
#' @param max_iter an integer.
#' @param eps.abs a scalar.
#' @param eps.rel a scalar.
#' @param verbose a logical.
#' @param print.type a logical.
#'
#' @return A list.
#' @noRd

admm_pdglasso_internal_r <- function(S,
                                     X,
                                     lambda1,
                                     lambda2,
                                     rho1,
                                     rho2,
                                     varying.rho1,
                                     varying.rho2,
                                     max_iter,
                                     eps.abs,
                                     eps.rel,
                                     acr.type,
                                     n.row.F,
                                     verbose) {
  
  converged <- TRUE
  
  p <- nrow(S)
  q <- p/2
  
  U <- Z <- matrix(0,p,p)
  #
  # required for the computation of residuals,
  # in this case is is not necessary to use mat2vec_r
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
    Z <- admm_inner_r(X=X, U=U, rho1=rho1, lambda1=lambda1, lambda2=lambda2, rho2=rho2,
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
  
  internal.par <- list(X = X,
                       res_primal = r.kk,
                       res_dual = s.kk,
                       n.iter = k,
                       n_iter_rho1_update_last = n.iter.rho1_update_last,
                       last_rho1 = rho1,
                       eps_primal = eps.pri,
                       eps_dual = eps.dual,
                       converged = converged)
  
  return(internal.par)
}


