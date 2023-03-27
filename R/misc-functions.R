#' Summary statistics for coloured graphs for paired data
#'
#' This function
#'
#' @param pdColG a coloured graph for paired data.
#' @param print.summary logical  (default `TRUE`) indicating whether a summary should be printed.
#'
#' @return a list
#' @export
#'
#' @examples
#' #
#' pdColG.summarize(toy_data$pdColG)
pdColG.summarize <- function(pdColG, print.summary=TRUE){
  X <- G.split(pdColG)
  G <- X$G
  p <- nrow(G)
  q <- p/2
  if(is.null(X$G.sym))    X$G.sym <- matrix(0, nrow=q, ncol=q)
  if(is.null(X$G.across)) X$G.across <- matrix(0, nrow=q, ncol=q)
  vertex.sym <- diag(X$G.sym)
  diag(X$G.sym) <- 0

  # summary statistics for vertices
  n.col.vertices <- 2*sum(vertex.sym)

  # summary statistics for inside block edges
  n.UNcol.symm.inside.edges <- 2*sum(X$G[1:q, 1:q]*X$G[(q+1):p, (q+1):p])
  n.col.inside.edges   <- 2*sum(X$G.sym)
  n.inside.edges       <- sum(X$G[1:q, 1:q])+sum(X$G[(q+1):p, (q+1):p])+n.col.inside.edges

  # summary statistics for across block edges
  Gtmp <- X$G[1:q, (q+1):p]
  diag(Gtmp) <- 0
  n.UNcol.symm.across.edges <- sum(Gtmp*t(Gtmp))
  n.col.across.edges  <- 2*sum(X$G.across)
  n.across.edges       <- sum(X$G[1:q, (q+1):p])+n.col.across.edges

  # overall
  n.edges <- sum(X$G)+n.col.inside.edges+n.col.across.edges

  #
  if(print.summary){
    cat("\nOVERALL\n")
    cat("Number of vertices: ", p, " \n")
    cat("Number of edges: ", n.edges, " \n")
    cat("Graph density  : ", round(n.edges/(p*(p-1)/2), 4) , " \n \n")
    #
    cat("VERTICES\n")
    cat("number of coloured vertices: ", n.col.vertices, " \n \n", sep="")
    #
    cat("INSIDE BLOCK EDGES\n")
    cat("number of edges: ", n.inside.edges, "\n", sep="")
    cat("number of uncoloured symmetric  edges: ", n.UNcol.symm.inside.edges, "\n", sep="")
    cat("number of coloured  (symmetric) edges: ", n.col.inside.edges, " \n \n", sep="")
    #
    cat("ACROSS BLOCK EDGES\n")
    cat("number of edges: ", n.across.edges, "\n", sep="")
    cat("number of uncoloured symmetric  edges: ", n.UNcol.symm.across.edges, "\n", sep="")
    cat("number of coloured  (symmetric) edges: ", n.col.across.edges, " \n \n", sep="")
  }
  overall <- list(n.vertices=p, n.edges=n.edges)
  vertices <- list(n.col.vertices=n.col.vertices)
  inside   <- list(n.edges=n.inside.edges, n.UNcol.symm.edges=n.UNcol.symm.inside.edges, n.col.edges=n.col.inside.edges)
  across   <- list(n.edges=n.across.edges, n.UNcol.symm.edges=n.UNcol.symm.across.edges, n.col.edges=n.col.across.edges)
  return(list(overall=overall, vertices=vertices, inside=inside, across=across))
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






## Extracts the LL block from a matrix
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
## Extracts the RR block from a matrix
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
## Extracts an LR block from a matrix
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

## Computes maximum theoretical values for lambda_1 and lambda_2
max.lams <- function(S){
  max.l1 <- max(abs(S))
  diff.inside <- abs(LL.block(S)-RR.block(S))
  diff.diag <- diag(diff.inside)/2
  diag(diff.inside) <- 0
  diff.across <- abs(across.block(S)-t(across.block(S)))
  diag(diff.across) <- 0
  max.l2 <- max(max(diff.inside),max(diff.diag),max(diff.across))
  return(c(max.l1,max.l2))
}
