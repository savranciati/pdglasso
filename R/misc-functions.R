
#' Extracts the LL block from a matrix.
#'
#' @param X a matrix.
#' @param new.val optional.
#'
#' @return A matrix.
#' @noRd

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




#' Extracts the RR block from a matrix.
#'
#' @param X a matrix.
#' @param new.val optional.
#'
#' @return A matrix.
#' @noRd

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




#' Extracts the LR block from a matrix.
#'
#' @param X a matrix.
#' @param new.val optional.
#'
#' @return A matrix.
#' @noRd

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




#' Compute maximum theoretical values for `lambda1` and `lambda2`.
#'
#' Computes the maximum values for `lambda1` and `lambda2` such that:
#' * if max of `lambda1` is used, the estimated concentration matrix will be diagonal;
#' * if max of `lambda2` is used, the estimated concentration matrix will be fully symmetric.
#'
#' @param S a covariance matrix.
#'
#' @return A vector of two elements.
#'
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' lams.max(S)

lams.max <- function(S){
  max.l1 <- max(abs(S[upper.tri(S, diag=FALSE)]))
  diff.inside <- abs(LL.block(S)-RR.block(S))/2
  diff.across <- abs(across.block(S)-t(across.block(S)))/2
  max.l2 <- max(max(diff.inside),max(diff.across))
  return(c(max.l1,max.l2))
}