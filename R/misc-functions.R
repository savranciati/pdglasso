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
