
#' Visual representation of a coloured graph for paired data
#'
#' Heatmap-style graphical representation of an object of class `pdColG`,
#' such as that returned by a call to  [`pdColG.get`].
#'
#' @param x a matrix of class `pdColG`; see [`pdglasso-package`] for details.
#' @param y not used. Included for compatibility with the generic `plot()` method.
#' @param uncol.sym a logical; if `TRUE`, structural symmetries are highlighted in the heatmap.
#' @param col.sym a logical; if `TRUE`, coloured structural symmetries are highlighted in the heatmap.
#' @param ... unused. Included for compatibility with the generic `plot()` method.
#' 
#' @return A plot visualizing the coloured graph enconded in the input matrix `x'. Each square of the heatmap can either be empty (white) for the diagonal and missing edges or represent a type of edge:
#' 
#' * a `darkgrey` square represents an edge in the graph;
#' 
#' * a `deepskyblue2` square represents an uncoloured structural symmetry;
#' 
#' * a `purple2` square represents either a vertex symmetry or a coloured structural symmetry.
#'
#' @export
#' @importFrom stats heatmap
#' 
#' @examples
#' S <- var(toy_data$sample.data)
#' mod.out <- admm.pdglasso(S, lambda1=4, lambda2=0.6)
#' G <- pdColG.get(mod.out)
#' plot(G)
#' plot(G, uncol.sym=TRUE)
#' plot(G, uncol.sym=TRUE, col.sym=TRUE)

plot.pdColG <- function(x, 
                        y=NULL,
                        uncol.sym = FALSE,
                        col.sym = FALSE,
                        ...){
  p <- ncol(x)
  q <- p/2
  if(uncol.sym | col.sym){
    # store LR diagonal for manipulation purposes
    LR.diag <- diag(x[1:q, (q+1):p])
    diag(x[1:q, (q+1):p]) <- rep(0,q)
    # find structural symmetries inside blocks (LL, RR)
    cond.inside <- which((x[1:q,1:q]*x[(q+1):p,(q+1):p])==1 , arr.ind = T)
    # find structural symmetries across blocks (LR,RL)
    cond.across <- which((x[1:q, (q+1):p]*t(x[1:q, (q+1):p]))==1, arr.ind = T)
    cond.across[,2] <- cond.across[,2]+q
    cond.across <- rbind(cond.across,cond.across[,c(2,1)])
    # readjust LR diagonal
    diag(x[1:q, (q+1):p]) <- LR.diag
    
    if(uncol.sym & col.sym){
      # change 2 to 3 to represent coloured structural symmetries
      x[x==2] <- 3
      # change 1 to 2 to represent uncoloured structural symmetries
      ## inside
      x[cond.inside] <- 2
      x[cond.inside+q] <- 2
      ## and across
      x[cond.across] <- 2
      diag(x)[diag(x)<3] <- 0
      
      # plot heatmap
      heatmap(x,
              Rowv = NA,
              Colv = NA,
              revC = TRUE,
              symm = TRUE,
              add.expr = {
                legend(
                  "bottomleft",
                  legend = c("Uncoloured symmetries",
                             "Coloured symmetries"),
                  fill = c("deepskyblue2","purple2"),
                  border = "black",
                  bty = "o",
                  cex = 1,
                  text.width = NULL
                )
                abline(h = q + 0.5, lwd = 2, col="black")
                abline(v = q + 0.5, lwd = 2, col="black")
              },
              scale = "none",
              col=c("white","darkgrey","deepskyblue2","purple2"),
      )
    }else{
      if(uncol.sym){
        # treat 2 values as normal edges
        x[x==2] <- 1
        # and only highlight uncoloured structural symmetries
        ## inside
        x[cond.inside] <- 2
        x[cond.inside+q] <- 2
        ## and across
        x[cond.across] <- 2
        diag(x) <- 0
        
        # plot heatmap
        heatmap(x,
                Rowv = NA,
                Colv = NA,
                revC = TRUE,
                symm = TRUE,
                add.expr = {
                  legend(
                    "bottomleft",
                    legend = c("Uncoloured symmetries"),
                    fill = c("deepskyblue2"),
                    border = "black",
                    bty = "o",
                    cex = 1,
                    text.width = NULL
                  )
                  abline(h = q + 0.5, lwd = 2, col="black")
                  abline(v = q + 0.5, lwd = 2, col="black")
                },
                scale = "none",
                col=c("white","darkgrey","deepskyblue2"),
        )
      }else{
        diag(x)[diag(x)<2] <- 0
        
        # plot heatmap 
        heatmap(x,
                Rowv = NA,
                Colv = NA,
                revC = TRUE,
                symm = TRUE,
                add.expr = {
                  legend(
                    "bottomleft",
                    legend = c("Coloured symmetries"),
                    fill = c("purple2"),
                    border = "black",
                    bty = "o",
                    cex = 1,
                    text.width = NULL
                  )
                  abline(h = q + 0.5, lwd = 2, col="black")
                  abline(v = q + 0.5, lwd = 2, col="black")
                },
                scale = "none",
                col=c("white","darkgrey","purple2"),
        )
      }
    }
  }else{
    # Color 2 values as 1, treating uncoloured and coloured in the
    # same way (both arguments set to FALSE)
    diag(x) <- 0
    
    heatmap(x,
            Rowv = NA,
            Colv = NA,
            revC = TRUE,
            symm = TRUE,
            add.expr = {
              abline(h = q + 0.5, lwd = 2, col="black")
              abline(v = q + 0.5, lwd = 2, col="black")
            },
            scale = "none",
            col=c("white","darkgrey","darkgrey")
    )
  }
}



#' Structural properties of a coloured graph for paired data
#'
#' Summary statistics relative to the structural properties 
#' of an object of class `pdColG`, such as that raturned by a call to  [`pdColG.get`].
#'
#' @param object a matrix of class `pdColG`; see [`pdglasso-package`] for details.
#' @param print.summary a logical; if `TRUE` the summary is printed on the console.
#' @param ... not used. Included for compatibility with the generic `summary()` method.
#'
#' @return An invisible list with the following components:
#'
#' * `overall` a list with the number of vertices, edges and free parameters of the associated pdRCON model.
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
#' model.summary <- summary(toy_data$pdColG, print.summary=FALSE)
#' model.summary
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
  
  # number of parameters
  n.par <- n.edges + p -(n.col.vertices+n.col.inside.edges+n.col.across.edges)/2
  
  #
  if(print.summary){
    cat("\nOVERALL\n")
    cat("number of vertices   :", p, " \n")
    cat("number of edges      :", n.edges, " \n")
    cat("graph density        :", round(n.edges/(p*(p-1)/2), 4) , " \n")
    cat("number of parameters :", n.par , " \n \n")
    #
    cat("VERTICES\n")
    cat("number of coloured vertices: ", n.col.vertices, " \n \n", sep="")
    #
    cat("INSIDE-BLOCK EDGES\n")
    cat("number of edges:", n.inside.edges, "\n", sep="")
    cat("number of uncoloured symmetric edges: ", n.UNcol.symm.inside.edges, "\n", sep="")
    cat("number of  coloured  symmetric edges: ", n.col.inside.edges, " \n \n", sep="")
    #
    cat("ACROSS-BLOCK EDGES\n")
    cat("number of edges: ", n.across.edges, "\n", sep="")
    cat("number of uncoloured symmetric edges: ", n.UNcol.symm.across.edges, "\n", sep="")
    cat("number of  coloured  symmetric edges: ", n.col.across.edges, " \n \n", sep="")
  }
  overall <- list(n.vertices=p, n.edges=n.edges, n.par=n.par)
  vertex <- list(n.col.vertices=n.col.vertices)
  inside   <- list(n.edges=n.inside.edges, n.UNcol.symm.edges=n.UNcol.symm.inside.edges, n.col.edges=n.col.inside.edges)
  across   <- list(n.edges=n.across.edges, n.UNcol.symm.edges=n.UNcol.symm.across.edges, n.col.edges=n.col.across.edges)
  return(invisible(list(overall=overall, vertex=vertex, inside=inside, across=across)))
}


#' Diagnostic plot for the output of the ADMM
#'
#' A plot, with different panels, obtained from the output of [`admm.pdglasso`], which is helpful to check the convergence of the ADMM.  
#' This function can also be used to visualize the default threshold used by [`pdColG.get`] as well to identify specific `th1` and `th2` threshold values to be used in
#' place of the default one. 
#' 
#'
#' @param x an object of class `ADMMoutput` such as the output of a call to the [`admm.pdglasso`] function; see also [`pdRCON.select`].
#' @param y  not used. Included for compatibility with the generic `plot()` method.
#' @param add.default.th a logical; if `TRUE` the default threshold used in [`pdColG.get`] is represented (in log10-scale) in the plot.
#' @param th1,th2 two positive scalars as in [`pdColG.get`]; if not `NULL` they are represented (in log10-scale) in the plot. 
#' @param logarithm10 a logical; if `TRUE` the values of `th1` and `th2` are expected to be provided in log10-scale. This facilitate  
#' interaction with the plot whose y-axis is in log10 scale. 
#' @param ... not used. Included for compatibility with the generic `plot()` method.
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
#' According to the values of the arguments `add.default.th`, `th1` and `th2`, horizontal dashed lines with the corresponding thresholds are included in the panels. 
#' 
#' The function returns an invisible vector with the values of `th1`, `th2`
#' and of the default threshold. The quantities `th1` and `th2`
#' can directly be given as input for the same arguments of [`pdColG.get`]. Consistently with the
#' behavior of [`pdColG.get`], a `NULL` value of `th1` or `th2` is replaced by the default threshold. 
#' 
#' @importFrom graphics par title abline legend
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' mod.out <- admm.pdglasso(S, lambda1=4, lambda2=0.3)
#' 
#' # full plot with four panels
#' plot(mod.out)
#' 
#' # model without across-block symmetries
#' mod.out <- admm.pdglasso(S, lambda1=4, lambda2=0.3, type=c("v", "i"))
#' plot(mod.out)

plot.ADMMoutput <- function(x, y=NULL, add.default.th=TRUE, th1=NULL, th2=NULL, logarithm10=TRUE, ...){
  y <- add.default.th
  mod.out <- x
  acronyms <- mod.out$acronyms$acronym.of.type
  th.default <- log10(mod.out$internal.par$eps.rel*10)
  
  if(!is.null(th1) & !logarithm10) th1 <- log10(th1)
  if(!is.null(th2) & !logarithm10) th2 <- log10(th2)
  
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




#' Maximum theoretical values of lambda1 and lambda2
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



