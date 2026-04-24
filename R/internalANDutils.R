
#' Visual representation of a coloured graph for paired data
#'
#' Heatmap-style graphical representation of an object of class `pdColG`,
#' such as that raturned by a call to  [`pdColG.get`].
#'
#' @param x a matrix of class `pdColG`; see [`pdglasso-package`] for details.
#' @param y a string; allows the user to specify which portion of the graph (a block) is to be plotted: possible mutually exclusive options are "left", "right", "across", "everything";  default is "everything", which means the entire graph is plotted.
#' @param which.sym a string or a vector of strings; specifies which kind of symmetries are plotted (structural, parametric or both).
#' @param asym.edges a logical; if `TRUE`, Asymmetric edges are plotted otherwise their symbol is suppressed
#' @param export.plot a logical; if `TRUE`, a .pdf file is produced in the working director instead of plotting the graph as a new panel/window.
#' @param fancy a logical value; if `TRUE`, symbols in the plot are used from the
#'   Latin1 encoding for characters; set to FALSE if characters are not properly
#'   displayed, so symbols are reverted to latin letters.
#'
#' @return Either a plot within the running R session or a .pdf file saved in the working directory if the option export.plot is set to `TRUE`.
#' @importFrom grDevices pdf dev.off dev.flush
#' @importFrom graphics par image text mtext segments
#' @export
#'
#' @examples
#' plot(toy_data$pdColG)

plot.pdColG <- function(x, 
                        y = "everything",
                        which.sym = c("structural","parametric"),
                        asym.edges = TRUE,
                        export.plot = FALSE,
                        fancy = TRUE, ...){
  # retain graph only (discard n.par) if G from pdColG.get is passed
  G <- x
  block <- y
  p <- dim(G)[1]
  q <- p/2
  strblock <- tolower(block)
  strblock <-  match.arg(strblock, c("left", "right","across","everything") , several.ok = FALSE)
  
  # define plot symbols
  G.symbols <- list()
  if(fancy){
    G.symbols$asym_edge <- "~"
    G.symbols$struct_sym <- "o"
    G.symbols$param_sym <- "."
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
  if (any(choice=="structural"))   strsym <- paste(strsym, "S", sep="")
  if (any(choice=="parametric"))   strsym <- paste(strsym, "P", sep="")
  switch(strsym,
         #
         S  = G.symbols$param_sym <- NA,
         P  = G.symbols$struct_sym <- NA,
         SP = G.symbols <- G.symbols
  )
  
  # Store LR diag for temporary manipulation
  LR.diag <- diag(G[1:q, (q+1):p])
  diag(G[1:q, (q+1):p]) <- rep(0,q)
  # convert content of G to symbols
  cond.inside <- which( (G[1:q,1:q]*G[(q+1):p,(q+1):p])==1 , arr.ind = T)
  cond.across <- which((G[1:q, (q+1):p]*t(G[1:q, (q+1):p]))==1, arr.ind = T)
  cond.across[,2] <- cond.across[,2]+q
  
  G[cond.inside] <- G.symbols$struct_sym ### structural symmetries LL RR
  #given the output of array.ind=TRUE only produces positions for LL, need to overwrite the RR block too
  G[cond.inside+q] <- G.symbols$struct_sym ### structural symmetries LL RR
  #given the output of array.ind=TRUE only produces positions for LL, need to overwrite adjust for LR block
  G[cond.across] <- G.symbols$struct_sym ### structural symmetries LR
  diag(G[1:q, (q+1):p]) <- LR.diag
  
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






#' Structural properties of a coloured graph for paired data
#'
#' Summary statistics relative to the structural properties 
#' of an object of class `pdColG`, such as that raturned by a call to  [`pdColG.get`].
#'
#' @param object a matrix of class `pdColG`; see [`pdglasso-package`] for details.
#' @param print.summary a logical, if `TRUE` the summary is printed on the console.
#'
#' @return An invisible list with the following components:
#'
#' * `overall` a list with the number of vertices and edges.
#'
#' * `vertex`  a list with the number of coloured vertices.
#'
#' * `inside`  a list with the number of inside-block edges, the number of uncolored symmetric and coloured inside block edges.
#'
#' * `across`  a list with the number of across-block edges, the number of uncolored symmetric and coloured across block edges.
#'
#' @export
#'
#' @examples
#' summary(toy_data$pdColG)
#' 
summary.pdColG <- function(object, print.summary=TRUE, ...){
  pdColG <- object
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


#' Diagnostic plot for the output of the ADMM
#'
#' A plot, with different panels, obtained from the output of [`admm.pdglasso`], which is helpful to check the convergence of the ADMM.  
#' This function can be used to visualize the default threshold used by [`pdColG.get`] as well to identify specific `th1` and `th2` threshold values to be used in
#' place of the default one. 
#' 
#'
#' @param x an object of class `ADMMoutput` such as the output of a call to the [`admm.pdglasso`] function; see also [`pdRCON.select`].
#' @param y a logical, if `TRUE` the default threshold used in [`pdColG.get`] is represented (in log10-scale) in the plot.
#' @param th1,th2 two positive scalars as in [`pdColG.get`], if not `NULL` they are represented (in log10-scale) in the plot. 
#' @param logarithm10 a logical, if `TRUE` the values of `th1` and `th2` are expected to be provided in log10-scale. This facilitate  
#' interaction with the plot whose y-axis is in log10 scale.   
#'
#'
#' @return A plot is produced with different panels depending on the model `type`. More specifically: 
#'
#' * a panel with the values of the off-diagonal concentrations is always provided,
#' 
#' * a panel with the differences of homologous diagonal concentration values is provided if vertex symmetries are allowed,
#' 
#' * a panel with the differences of homologous inside-block concentration values is provided if inside-block symmetries are allowed,
#' 
#' * a panel with the differences of homologous across-block concentration values is provided if across-block symmetries are allowed
#' 
#' According to the values of the arguments `y`, `th1` and `th2`, horizontal dashed lines with the corresponding thresholds are included in the panels. 
#' 
#' The function returns an invisible vector with the values of `th1`, `th2`
#' and of the default threshold. The quantities `th1` and `th2`
#' can directly given as input for the same arguments of [`pdColG.get`]. Consistently with the
#' behavior of [`pdColG.get`], a `NULL` value of `th1` or`th2` is replaced by the default threshold. 
#' 
#' @importFrom graphics par title abline legend
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' mod.out <- pdRCON.select(S,n=60)$model
#' plot(mod.out)

plot.ADMMoutput <- function(x, y=TRUE, th1=NULL, th2=NULL, logarith10=TRUE, ...){
  mod.out <- x
  acronyms <- mod.out$acronyms$acronym.of.type
  th.default <- log10(mod.out$internal.par$eps.rel*10)
  
  if(!is.null(th1) & !logarith10) th1 <- log10(th1)
  if(!is.null(th2) & !logarith10) th2 <- log10(th2)
  
  p <- ncol(mod.out$X)
  q <- p/2
  
  # Prepare canvas
  def.pars <- par(no.readonly = TRUE)
  #par(mfrow = c(3, 2))
  #par(mar = c(2, 2, 1.5, 1))
  if (nchar(acronyms)==3) layout(matrix(c(1, 2, 3, 4, 5, 5), nrow=3, byrow=TRUE), heights = c(3,3,1),  respect=FALSE)
  if (nchar(acronyms)==2) layout(matrix(c(1, 1, 2, 3, 4, 4), nrow=3, byrow=TRUE), heights = c(3,3,1),  respect=FALSE)
  if (nchar(acronyms)==1) layout(matrix(c(1, 2, 3, 3), nrow=2, byrow=TRUE), heights = c(3,1),  respect=FALSE)
  #
  # off-diagonal
  off.d <- mod.out$X[upper.tri(mod.out$X, diag=FALSE)]
  off.d <- log10(abs(off.d))
  off.d <- off.d[is.finite(off.d)]
  plot(off.d, ylim=c(min(off.d), 0), xlab="", ylab="")
  title(main="off-diagonal concentration values")
  if (y) abline(h=th.default, lty=2, col="red", lwd=2)
  if (!is.null(th1)) abline(h=th1, lty=2, col="green", lwd=2)
  #
  # vertices
  if(grepl("V",acronyms,fixed=T)){
    vertx <- diag(mod.out$X[1:q, 1:q]-mod.out$X[(q+1):p,(q+1):p])
    vertx <- log10(abs(vertx))
    vertx <- vertx[is.finite(vertx)]
    plot(vertx, ylim=c(min(off.d), 0), xlab="", ylab="")
    title(main="diff. of paired 'vertex' concentrations")
    if (y) abline(h=th.default, lty=2, col="red", lwd=2)
    if (!is.null(th2)) abline(h=th2, lty=2, col="blue", lwd=2)
  }
  #
  # inside-block edges
  if(grepl("I",acronyms,fixed=T)){
    inside.edges <- LL.block(mod.out$X)-RR.block(mod.out$X)
    inside.edges <- inside.edges[upper.tri(inside.edges, diag=FALSE)]
    inside.edges <- log10(abs(inside.edges))
    inside.edges <- inside.edges[is.finite(inside.edges)]
    plot(inside.edges, ylim=c(min(off.d), 0), xlab="", ylab="")
    title(main="diff. of paired 'inside-block' concentrations")
    if (y) abline(h=th.default, lty=2, col="red", lwd=2)
    if (!is.null(th2)) abline(h=th2, lty=2, col="blue", lwd=2)
  }
  #
  # across-block edges
  if(grepl("A",acronyms,fixed=T)){
    across.edges <- across.block(mod.out$X)-t(across.block(mod.out$X))
    across.edges <- across.edges[upper.tri(across.edges, diag=FALSE)]
    across.edges <- log10(abs(across.edges))
    across.edges <- across.edges[is.finite(across.edges)]
    plot(across.edges, ylim=c(min(off.d), 0), xlab="", ylab="")
    title(main="diff. of paired 'across-block' concentrations")
    if (y) abline(h=th.default, lty=2, col="red", lwd=2)
    if (!is.null(th2)) abline(h=th2, lty=2, col="blue", lwd=2)
  }
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  #plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("center",inset = 0,
         legend = c("log10( default th. )   ", "log10( th1 )", "log10( th2 )"),
         col=c("red", "green", "blue"), lty=2, lwd=1.5, cex=0.8, horiz = TRUE)
  
  # Reset canvas
  on.exit(layout(1), add=TRUE)
  on.exit(par(def.pars), add=TRUE)
  #
  if(is.null(th1)) th1 <- th.default
  if(is.null(th2)) th2 <- th.default
  return(invisible(c(th1=10^th1, th2=10^th2, th.default=10^th.default)))
}




#' Maximum theoretical values of `lambda1` and `lambda2`
#' 
#' Computes the maximum theoretical values `lambda1.max` and `lambda2.max` of `lambda1` and `lambda2`, respectively. 
#' For every `lambda1` greater than `lambda1.max` the pdglasso estimator is a diagonal matrix whereas 
#' for every `lambda2` greater than `lambda2.max` the pdglasso estimator is fully symmetric. 
#'
#' @param S a sample covariance matrix.
#'
#' @return A vector of two elements, `lambda1.max` and `lambda2.max`.
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
  return(c(lambda1.max=max.l1, lambda2.max=max.l2))
}








#' Assign the class `pdColGraph` to a pdColG matrix 
#'
#' @param x a pdColG matrix; see [`pdglasso-package`] for details. 
#' 
#'
#' @return the input but with the addition of the class `pdColG`
#' @noRd

as.pdColG <- function(x) {
  stopifnot(inherits(x, "matrix"))
  structure(x, class = c("pdColG", class(x)))
}


#' Assign the class `ADMMoutput` to the list produced by a call to either [`admm.pdglasso`] or [`pdRCON.select`]
#'
#' @param x and object as in the description. 
#' 
#'
#' @return the input but with the addition of the class `ADMMoutput`
#' @noRd

as.ADMMoutput <- function(x) {
  stopifnot(inherits(x, "list"))
  structure(x, class = c("ADMMoutput", class(x)))
}


#' Extracts the LL block from a matrix
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



#' Extracts the RR block from a matrix
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




#' Extracts the LR block from a matrix
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



#' Conversion from the single matrix representation of the model to the multiple matrix representation
#'
#' This is the inverse of the function `G.merge`, i.e. g is equal to `G.merge(G.split(g))`.
#'
#' @param g a p x p symmetric matrix with entries 0, 1, and 2.
#'
#' @return A list with three upper triangular matrices: G, G.sym and G.across with entries 0 and 1, any of G.sym and G.across may be NULL.
#' 
#' @noRd

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




#' Conversion from the multiple matrix representation of the model to the single matrix representation
#'
#' This is the inverse of the function `G.split`, i.e. X is equal to `G.split(G.merge(X))`.
#'
#' @param X list with three upper triangular matrices with entries 0 and 1: G, G.sym and G.across, any of G.sym and G.across may be NULL.
#'
#' @return a pXp symmetric matrix with entries 0, 1, and 2
#' @noRd
#'
#' @examples
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

#' Create the model acronym
#'
#' Creates the model acronym from the arguments `types` and `force.symm`.
#'
#' @param type a character vector.
#' @param force.symm either a character vector or `NULL`.
#' @param print.type logical (default `TRUE`) indicating whether the model details should be printed.
#'
#' @return A list with two character strings named `acronym.of.type` and `acronym.of.force`, the latter is  `NULL` if `force.symm` is   `NULL`.
#' @noRd

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




#' Imposes lambda values such that symmetry is forced
#'
#' @param p an even number.
#' @param lambda2 a scalar.
#' @param acr.type type of acronym.
#' @param acr.force type of acronym for forced symmetries.
#'
#' @return .
#' @noRd

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



