#' Build a graph from the output of a call to [`admm.pdglasso`].
#'
#' Description here.
#'
#' @param admm.out An object of list type, that is the output of a call to the [`admm-pdglasso`] function.
#' @param th1 (optional) A scalar, the threshold to identify edges in the graph; it must be non-negative.
#' @param th2 (optional) A scalar, the threshold to identify coloured edges in the graph; it must be non-negative.
#' @param verbose (optional) if `TRUE` provides summary statistics of the graph.
#'
#' @return a list, containing:
#' * g, the graph in matrix form.
#' * dof, the degrees of freedom corresponding to the graph build under the pdglasso model provided.
#' @export
#'
#' @examples
get.graph <- function(admm.out,
                      th1=NULL,
                      th2=NULL,
                      verbose=FALSE){
  # Prepare output object
  out <- list()
  # Store passed acronims used
  acronims <- admm.out$acronims$acronim.of.type
  # Store passed estimated concentration matrix
  X <- admm.out$X
  p <- dim(X)[1]
  q <- p/2

  if(is.null(th1)) th1 <- admm.out$internal.par$eps.rel
  if(is.null(th2)) th2 <- admm.out$internal.par$eps.rel

  # prepare matrix of edges (non-zero elements)
  # empty except diagonal
  mat_graph <- diag(1,p)
  # check which elements of precision matrix are above the threshold
  # set from 0 to 1 if above threshold (there is an edge)
  # diagonal not affected
  mat_graph[which(abs(X)>th1, arr.ind=TRUE)] <- 1
  # check if the resulting matrix is still symmetrical
  if(isSymmetric(mat_graph)==FALSE) warning(paste("Non-symmetric mat_graph"))
  # mat_graph now contains only information on edges

  # prepare matrix for storing symmetries
  # empty be default = no symmetric concentrations
  mat_sym <- matrix(0,p,p)

  ### between symmetries (LL, RR)
  if(grepl("V",acronims,fixed=T)){
    # computes the vertices differences (in absolute values)
    # | diag(LL)-diag(RR) |
    diff_vertex  <- abs(diag(X)[1:q]-diag(X)[(q+1):p])
    temp_vertex <- (!(diff_vertex>th2))+0
    temp_vertex <- c(temp_vertex,temp_vertex)
    diag(mat_sym) <- temp_vertex
  }

  ### between symmetries (LL, RR)  TOGLIERE LA DIAGONALE DAL BLOCCO I
  if(grepl("I",acronims,fixed=T)){
    # checks if there is at least one edge at each (LL_ij, RR_ij)
    block_inside <- pmax(LL.block(mat_graph), RR.block(mat_graph))
    # computes the inside differences (in absolute values)
    # | LL-RR |
    diff_inside  <- abs(LL.block(X)-RR.block(X))

    # creates temporary (q x q) matrix where all concentrations are the same
    temp_inside <- matrix(1,q,q)
    # any difference | LL_ij - RR_ij | > th2 is too large and
    # thus LL_ij != RR_ij
    temp_inside[which(diff_inside>th2,arr.ind=TRUE)] <- 0
    temp_inside <- temp_inside*block_inside
    diag(temp_inside) <- 0 # vertices equalities are checked inside V "if" above
    mat_sym <- LL.block(mat_sym,temp_inside)
    mat_sym <- RR.block(mat_sym,temp_inside)
  }

  ### across symmetries (LR, RL)
  if(grepl("A",acronims,fixed=T)){

    # checks if there is at least one edge at each (LR_ij, RL_ij)
    block_across <- pmax(across.block(mat_graph), t(across.block(mat_graph)))

    # computes the between differences (in absolute values)
    # | LR-RL |
    diff_across  <- abs(across.block(X)-t(across.block(X)))
    # creates temporary (q x q) matrix where all concentrations are the same
    temp_across     <- matrix(1,q,q)
    # any difference | LR_ij - RL_ij | > th2 is too large and
    # thus LL_ij != RR_ij
    temp_across[which(diff_across>th2, arr.ind=TRUE)] <- 0
    # diagonals of LR and RL are never penalized so
    # they cannot be set as equal concentrations
    diag(temp_across) <- 0
    temp_across  <- temp_across*block_across
    # enforcing symmetry on temp2
    temp_across  <- pmax(temp_across, t(temp_across))

    # update mat_sym with information on equal across concentrations
    mat_sym <- across.block(mat_sym,temp_across)
    mat_sym <- pmax(mat_sym, t(mat_sym))
  }


  ### Force coherence between edges and symmetric concentration values
  mat_graph   <- pmax(mat_graph,mat_sym)

  # computing degrees of freedom (+p before diving because math_graph has 1s on diagonal)
  tot.dof  <- (sum(mat_graph)+p)/2
  # number of pairs of symmetric off-diagonal concentrations
  nsym_offdiag  <- sum(LL.block(mat_sym)[upper.tri(LL.block(mat_sym),diag=F)]) +
                   sum(across.block(mat_sym)[upper.tri(across.block(mat_sym),diag=F)])
  # number of pairs of symmetric diagonal concentrations
  nsym_diag <- sum(diag(mat_sym))/2
  dof <- tot.dof - nsym_offdiag - nsym_diag
  #
  out <- list()
  out$g   <- mat_graph+mat_sym
  out$dof <- dof
  if(verbose){
    graph.stats(G.split(out$g))
    cat("\n\n")
  }
  return(out)
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
