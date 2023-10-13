#' Build a graph from the output of a call to [`admm.pdglasso`].
#'
#' This function returns a Coloured Graph for Paired Data (pdColG) from the
#' output of a call to [`admm.pdglasso`]. We refer to [`pdglasso-package`] both for the description of the matrix
#' encoding a coloured graph for paired and the available submodel classes.
#'
#' @param admm.out An object of list type, that is the output of a call to the [`admm.pdglasso`] function.
#' @param th1 (optional) A scalar, the threshold to identify edges in the graph; it must be non-negative.
#' @param th2 (optional) A scalar, the threshold to identify coloured edges in the graph; it must be non-negative.
#' @param print.summary (optional) if `TRUE` provides summary statistics of the graph.
#'
#' @return a list with the following components:
#'
#' * `pdColG` a matrix representing a coloured graph for paired data; see [`pdglasso-package`] for details.
#'
#' * `n.par` the number of parameters of the pdRCON model represented by `pdColg`.
#' @export
#'
#' @examples
#'
#' S <- cov(toy_data$sample.data)
#' mod.out <- admm.pdglasso(S)
#' pdColG.get(mod.out)
pdColG.get <- function(admm.out,
                      th1=NULL,
                      th2=NULL,
                      print.summary=FALSE){
  # Prepare output object
  out <- list()
  # Store passed acronyms used
  acronyms <- admm.out$acronyms$acronym.of.type
  # Store passed estimated concentration matrix
  X <- admm.out$X
  p <- dim(X)[1]
  q <- p/2

  if(is.null(th1)) th1 <- admm.out$internal.par$eps.rel*100
  if(is.null(th2)) th2 <- admm.out$internal.par$eps.rel*100

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

  ### between symmetries (LL, RR)  TOGLIERE LA DIAGONALE DAL BLOCCO I
  if(grepl("I",acronyms,fixed=T)){
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
  if(grepl("A",acronyms,fixed=T)){

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


  ### between symmetries (LL, RR)
  if(grepl("V",acronyms,fixed=T)){
    # computes the vertices differences (in absolute values)
    # | diag(LL)-diag(RR) |
    diff_vertex  <- abs(diag(X)[1:q]-diag(X)[(q+1):p])
    # create temporary (1:q) half-diag where all vertex values are equal
    temp_vertex <- rep(1,q)
    # any | LL_ii - RR_i'i' | > th2 is too large and
    # thus LL_ii != RR_i'i'
    temp_vertex[diff_vertex>th2] <- 0
    temp_vertex <- c(temp_vertex,temp_vertex)
    diag(mat_sym) <- temp_vertex
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
  n.par <- tot.dof - nsym_offdiag - nsym_diag
  #
  out <- list()
  out$pdColG   <- mat_graph+mat_sym
  dimnames(out$pdColG) <- dimnames(admm.out$X)
  out$n.par <- n.par
  if(print.summary){
    pdColG.summarize(out$pdColG)
    cat("\n\n")
  }
  return(out)
}


#' Structural properties of a coloured graph for paired data
#'
#' This function returns some summary statistics relative to the structural properties of a
#' coloured graph for paired data \eqn{\mathcal{G}}. We refer to  [`pdglasso-package`] both for the description of the
#' matrix encoding a coloured graph for paired and details on how the structural quantities are defined and computed.
#'
#' @param pdColG a matrix representing a coloured graph for paired data; see [`pdglasso-package`] for details.
#' @param print.summary a logical  (default `TRUE`) indicating whether a summary should be printed.
#'
#' @return
#'
#' An invisible list with the following components:
#'
#' * `overall` a list with the number of vertices and edges of  \eqn{\mathcal{G}}.
#'
#' * `vertex`  a list with the number of coloured vertices of  \eqn{\mathcal{G}}.
#'
#' * `inside`  a list with the number of inside-block edges, the number of uncolored symmetric and coloured inside block edges of \eqn{\mathcal{G}}.
#'
#' * `across`  a list with the number of across-block edges, the number of uncolored symmetric and coloured across block edges of \eqn{\mathcal{G}}.
#'
#' If `print.summary=TRUE` some summary statistics are also printed on the console
#'
#' @export
#'
#' @examples
#'
#' pdColG.summarize(toy_data$pdColG)
#'
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

  # summary statistics for inside-block edges
  n.UNcol.symm.inside.edges <- 2*sum(X$G[1:q, 1:q]*X$G[(q+1):p, (q+1):p])
  n.col.inside.edges   <- 2*sum(X$G.sym)
  n.inside.edges       <- sum(X$G[1:q, 1:q])+sum(X$G[(q+1):p, (q+1):p])+n.col.inside.edges

  # summary statistics for across-block edges
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
    cat("number of vertices: ", p, " \n")
    cat("number of edges: ", n.edges, " \n")
    cat("graph density  : ", round(n.edges/(p*(p-1)/2), 4) , " \n \n")
    #
    cat("VERTICES\n")
    cat("number of coloured vertices: ", n.col.vertices, " \n \n", sep="")
    #
    cat("INSIDE-BLOCK EDGES\n")
    cat("number of edges: ", n.inside.edges, "\n", sep="")
    cat("number of uncoloured symmetric edges: ", n.UNcol.symm.inside.edges, "\n", sep="")
    cat("number of coloured  symmetric edges: ", n.col.inside.edges, " \n \n", sep="")
    #
    cat("ACROSS-BLOCK EDGES\n")
    cat("number of edges: ", n.across.edges, "\n", sep="")
    cat("number of uncoloured symmetric edges: ", n.UNcol.symm.across.edges, "\n", sep="")
    cat("number of coloured  symmetric edges: ", n.col.across.edges, " \n \n", sep="")
  }
  overall <- list(n.vertices=p, n.edges=n.edges)
  vertex <- list(n.col.vertices=n.col.vertices)
  inside   <- list(n.edges=n.inside.edges, n.UNcol.symm.edges=n.UNcol.symm.inside.edges, n.col.edges=n.col.inside.edges)
  across   <- list(n.edges=n.across.edges, n.UNcol.symm.edges=n.UNcol.symm.across.edges, n.col.edges=n.col.across.edges)
  return(invisible(list(overall=overall, vertex=vertex, inside=inside, across=across)))
}


#' Visualize a coloured graph for paired data
#'
#' This functions produces a heatmap-style graphical representation of a graph object produced by
#' a call to [`pdColG.get`], or a matrix compatible with the coloured graph specification
#' used in this package.
#'
#' @param g a symmetric \eqn{p \times p} matrix coding a coloured graphical
#'   model; only values allowed are the integers \eqn{\{0, 1, 2\}}; usually, an
#'   object produced by a call to the function [`pdColG.get`]
#'
#' @param which.sym a string or a vector of strings; specifies which kind of symmetries are plotted (structural, parametric or both).
#' @param block a string; allows the user to specify which portion of the graph (a block) is to be plotted: possible mutually exclusive options are "left", "right", "across", "everything";  default is "everything", which means the entire graph is plotted.
#' @param asym.edges a logical value; if TRUE, Asymmetric edges are plotted otherwise their symbol is suppressed
#' @param export.plot a logical value; if TRUE, a .pdf file is produced in the working director instead of plotting the graph as a new panel/window.
#' @param fancy a logical value; if TRUE, symbols in the plot are used from the
#'   Latin1 encoding for characters; set to FALSE if characters are not properly
#'   displayed, so symbols are reverted to latin letters.
#'
#' @return either a plot within the running R session or a .pdf file saved in the working directory if the option export.plot is set to TRUE.
#' @export
#'
#' @examples
#'
#' pdColG.plot(toy_data$pdColG)
pdColG.plot <- function(G,
                        block = "everything",
                        which.sym = c("structural","parametric"),
                        asym.edges = TRUE,
                        export.plot = FALSE,
                        fancy = TRUE){
  # retain graph only (discard n.par) if G from pdColG.get is passed
  if(is.list(G)) G <- G$pdColG
  p <- dim(G)[1]
  q <- p/2
  strblock <- tolower(block)
  strblock <-  match.arg(strblock, c("left", "right","across","everything") , several.ok = FALSE)

  # define plot symbols
  G.symbols <- list()
  if(fancy){
    G.symbols$asym_edge <- "~"
    G.symbols$struct_sym <- "o"
    G.symbols$param_sym <- "·"
  }else{
    G.symbols$asym_edge <- "-"
    G.symbols$struct_sym <- "o"
    G.symbols$param_sym <- "x"
  }
  # check if asymmetric edges are to be plotted
  if(!asym.edges) G.symbols$asym_edge <- NA

  strsym <- tolower(which.sym)
  choice <-  match.arg(strsym, c("structural", "parametric") , several.ok = TRUE)
  strsym <- ""
  if (any(choice=="structural"))    strsym <- paste(strsym, "S", sep="")
  if (any(choice=="parametric"))   strsym <- paste(strsym, "P", sep="")
  switch(strsym,
         #
         S = G.symbols$param_sym <- NA,
         P = G.symbols$struct_sym <- NA,
         SP = G.symbols <- G.symbols
  )

  # convert content of G to symbols
  cond.inside <- which( (G[1:q,1:q]*G[(q+1):p,(q+1):p])==1 , arr.ind = T)
  cond.across <- which((G[1:q, (q+1):p]*t(G[1:q, (q+1):p]))==1, arr.ind = T)
  cond.across[,2] <- cond.across[,2]+q

  G[cond.inside] <- G.symbols$struct_sym ### structural symmetries LL RR
  #given the output of array.ind=TRUE only produces positions for LL, need to overwrite the RR block too
  G[cond.inside+q] <- G.symbols$struct_sym ### structural symmetries LL RR
  #given the output of array.ind=TRUE only produces positions for LL, need to overwrite adjust for LR block
  G[cond.across] <- G.symbols$struct_sym ### structural symmetries LR
  #diag(G[1:q, (q+1):p]) <- G.symbols$asym_edge

  G[G==0] <- NA
  G[G==1] <- G.symbols$asym_edge ### asymmetric edges
  G[G==2] <- G.symbols$param_sym ### parametric symmetries
  G[lower.tri(G)] <- NA


  # selects portion or total graph
  switch(strblock,
         left = G <- G[1:q, 1:q],
         right = G <- G[(q+1):p, (q+1):p],
         across = G <- G[1:q, (q+1):p],
         everything = G <- G
  )
  p <- dim(G)[1]
  q <- p/2

  # adjust graphical parameters
  cex.matrix <- matrix(1,p,p)
  adj.matrix <- matrix(0.5,p,p)
  if(fancy){
    cex.matrix[G==G.symbols$param_sym] <- 4
    adj.matrix[G==G.symbols$param_sym] <- 0.55
  }
  # store original graphical parameters
  op <- par(no.readonly = TRUE)

  # check if plot is to be saved as a pdf
  if(export.plot) pdf(paste("graph_",strblock,".pdf",sep=""), width = 7, height = 7, onefile = TRUE)

  par(oma=c(0,0,3,1.5))
  par(mai=c(0.1,0.1,0.2,0.4))

  image(
    x = 1:p,        # X-axis values (columns)
    y = 1:p,        # Y-axis values (rows)
    z = matrix(0,p,p),                # Data matrix
    col = "white",
    xlab = "",
    ylab = "",
    main = "",
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    asp=1
  )
  text(
    x = rep(1:p, each = p),
    y = rep(p:1, times = p),
    labels = G,
    cex = cex.matrix,                         # Adjust label size if needed
    adj = adj.matrix,                           # Position text below each cell
    col = "black"                      # Text color
  )

  mtext( text = colnames(G), at = 1:p, side = 3, line=1, las = 2, outer=FALSE)
  mtext( text = rownames(G), at = p:1, side = 4, line=1, las = 2, outer=FALSE)

  # guides
  if(strblock=="everything"){
    segments(y0=p/2 + 0.5, y1=p/2 + 0.5, x0=1-0.25, x1=p+0.25, col=1, lwd=2, lty=2)
    segments(x0=p/2 + 0.5, x1=p/2 + 0.5, y0=1-0.25, y1=p+0.25, col=1, lwd=2, lty=2)
  }
  if(strblock=="across"){
    segments(y0=p + 0.5, y1=1 + 0.0 , x0=1-0.25, x1=p+0.25, col=1, lwd=1, lty=2)
    # 0.5 shift
    segments(y0=p + 0.0, y1=1 - 0.5 , x0=1-0.25, x1=p+0.25, col=1, lwd=1, lty=2)
  }

  par(op)
  if(export.plot) dev.off()
  on.exit(dev.flush())
}

#' Conversion from the single matrix representation of the model to the multiple matrix representation.
#'
#' This is the inverse of the function [G.merge()], i.e. g is equal to `G.merge(G.split(g))`.
#'
#' @param g is a pXp symmetric matrix with entries 0, 1, and 2
#'
#' @return a list with three upper triangular matrices: G, G.sym and G.across with entries 0 and 1, any of G.sym and G.across may be NULL
#'
#' @noRd
#'
#' @examples
#'
#' # see example in function [G.merge()].
G.split <- function(g) {
  p <- dim(g)[1]
  q <- p / 2
  G <- (g == 1) * 1
  G[lower.tri(G, diag = TRUE)] <- 0

  # matrix "sym"
  if (any(g[1:q, 1:q] == 2)) {
    G.sym <- (g[1:q, 1:q] == 2) * 1
    G.sym[lower.tri(G.sym, diag = FALSE)] <- 0
  }else{
    G.sym <- NULL
  }

  # matrix "across"
  if (any(g[1:q, (q + 1):p] == 2)) {
    G.across <- (g[1:q, (q + 1):p] == 2) * 1
    G.across[lower.tri(G.across, diag = TRUE)] <- 0
  } else{
    G.across <- NULL
  }
  return(list(
    G = G,
    G.sym = G.sym,
    G.across = G.across
  ))
}



#' Conversion from the multiple matrix representation of the model to the single matrix representation.
#'
#' This is the inverse of the function [G.split], i.e. X is equal to `G.split(G.merge(X))`.
#'
#' @param X list with three upper triangular matrices with entries 0 and 1: G, G.sym and G.across, any of G.sym and G.across may be NULL
#'
#' @return a pXp symmetric matrix with entries 0, 1, and 2
#' @noRd
#'
#' @examples
#'
#' # random generation of a list(G=G, G.sym=G.sym, G.across=G.across)
#'
#' q <- 5 # this can be any integer
#' p <- q*2
#'
#' g.p <- matrix(sample(c(0,1), size=p^2, replace=TRUE), nrow=p, ncol=p)
#' g.p[lower.tri(g.p, diag=TRUE)] <- 0
#'
#' g.q1 <- matrix(sample(c(0,1), size=q^2, replace=TRUE), nrow=q, ncol=q)
#' g.q1[lower.tri(g.q1, diag = FALSE)] <- 0
#'
#' g.q2 <- matrix(sample(c(0,1), size=q^2, replace=TRUE), nrow=q, ncol=q)
#' g.q2[lower.tri(g.q2, diag = TRUE)] <- 0
#' g.q2sym <- g.q2+t(g.q2)
#'
#' g.p[1:q, 1:q] <- g.p[1:q, 1:q] *(1-g.q1)
#' g.p[(q+1):p, (q+1):p] <- g.p[(q+1):p, (q+1):p] *(1-g.q1)
# '
#' g.p[1:q, (q+1):p] <- g.p[1:q, (q+1):p] * (1-g.q2sym)
#'
#' # list obtained
#'
#' X <- list(G=g.p, G.sym=g.q1, G.across=g.q2)
#'
#' g  <- G.merge(X)
#' gs <- G.split(g)
#'
#' identical(X, gs)
#'
#' X <- list(G=g.p, G.sym=NULL, G.across=NULL)
#'
#' g  <- G.merge(X)
#' gs <- G.split(g)
#'
#' identical(X, gs)
#'
G.merge <- function(X) {
  p <- nrow(X$G)
  q <- p / 2
  G <- X$G + t(X$G) + diag(1, p)
  # G.sym
  if (!is.null(X$G.sym)) {
    diag(G) <- 0
    G.sym <- X$G.sym + t(X$G.sym)
    G[1:q, 1:q] <- G[1:q, 1:q] + 2 * G.sym
    G[(q + 1):p, (q + 1):p] <- G[(q + 1):p, (q + 1):p] + 2 * G.sym
    diag(G[1:q, 1:q]) <- diag(G[(q + 1):p, (q + 1):p]) <- diag(X$G.sym) + 1
  }
  # G.across
  if (!is.null(X$G.across)) {
    G.across <- X$G.across + t(X$G.across)
    G[1:q, (q + 1):p] <- G[1:q, (q + 1):p] + 2 * G.across
    G[(q + 1):p, 1:q] <- G[(q + 1):p, 1:q] + 2 * t(G.across)
  }
  return(G)
}

