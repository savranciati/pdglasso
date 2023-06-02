#' Half-vectorization operator,  without the diagonal
#'
#' @param a symmetric matrix M
#'
#' @return x
#' @noRd
#'
half.vec <- function(M){
  return(M[upper.tri(M, diag=FALSE)])
}

#' Inverse of Half-vectorization operator, without the diagonal
#'
#' @param m a vector
#'
#' @return a triangular (upper) matrix
#' @noRd
#'
inv_half.vec<- function(m){
  mat.dim <- (1+(1+8*length(m))^0.5)/2
  M <- matrix(0, mat.dim, mat.dim)
  M[upper.tri(M, diag=FALSE)] <- m
  return(M)
}

#
#' Transforms a symmetric matrix into a vector applying the
# v() operator
#'
#' @param M a matrix
#'
#' @return a vector
#' @noRd
#'
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

#' Inverse operation with respect to mat2vec().
#'
#' @param m a vector
#'
#' @return a matrix
#' @noRd
#'
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

#' Extracts the LL block from a matrix
#'
#' @param X a matrix
#' @param new.val optional
#'
#' @return a matrix
#' @noRd
#'
LL.block <- function(X, new.val=NULL){
  p   <- dim(X)[1]
  q   <- p/2
  if(is.null(new.val)){
    return(X[1:q,1:q])
  }else{
    X[1:q,1:q] <- new.val
    return(X)
  }
}

#' Extracts the RR block from a matrix
#'
#' @param X a matrix
#' @param new.val optional
#'
#' @return a matrix
#' @noRd
#'
RR.block <- function(X, new.val=NULL){
  p   <- dim(X)[1]
  q   <- p/2
  if(is.null(new.val)){
    return(X[(q+1):p,(q+1):p])
  }else{
    X[(q+1):p,(q+1):p] <- new.val
    return(X)
  }
}

#' Extracts the LR block from a matrix
#'
#' @param X a matrix
#' @param new.val optional
#'
#' @return a matrix
#' @noRd
#'
across.block <- function(X, new.val=NULL){
  p   <- dim(X)[1]
  q   <- p/2
  if(is.null(new.val)){
    return(X[1:q,(q+1):p])
  }else{
    X[1:q,(q+1):p] <- new.val
    return(X)
  }
}

#' Computes maximum theoretical values for lambda_1 and lambda_2
#'
#' @param S a covariance matrix.
#'
#' @return a vector of two elements.
#' @noRd
#'
max.lams <- function(S){
  max.l1 <- max(abs(S[upper.tri(S, diag=FALSE)]))
  diff.inside <- abs(LL.block(S)-RR.block(S))/2
  diff.across <- abs(across.block(S)-t(across.block(S)))/2
  max.l2 <- max(max(diff.inside),max(diff.across))
  return(c(max.l1,max.l2))
}
