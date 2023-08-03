require(MASS)


#' Random simulation of pdRCON models
#'
#' Randomly generates a coloured graph for paired data \eqn{\mathcal{G}}, a
#' concentration matrix \eqn{K} adapted to \eqn{\mathcal{G}} and a random sample
#' from a multivariate normal distribution with zero mean vector and covariance
#' matrix \eqn{\Sigma=K^{-1}}.
#'
#'
#' @param p an even integer, that is the number of vertices of the generated
#'   coloured graph for paired data `pdColG`.
#' @param concent.mat a logical (default `TRUE`) indicating whether a
#'   concentration matrix `K` adapted to `pdColG` should be generated.
#' @param sample a logical (default `TRUE`) indicating whether a sample from a
#'   normal distribution with zero mean vector and concentration matrix `K`
#'   should be generated.
#' @param Sigma a \eqn{p\times p} positive definite matrix. This the matrix
#'   argument of the [`rWishart`] function, which is used as starting point in
#'   the random generation of the concentration matrix, as described in the
#'   details section below. The default `NULL` is equivalent to the identity
#'   matrix `Sigma=diag(p)`.
#' @param sample.size size of the randomly generated sample. The default `NULL`
#'   is equivalent to `sample.size=3*p`.
#' @param type,force.symm two subvectors of `c("vertex", "inside.block.edge", "across.block.edge")` which
#' identify the pdRCON submodel class of interest; see [`pdglasso-package`] for details.
#' @param dens,dens.vertex,dens.inside,dens.across four values between zero and
#'   one used to specify the sparsity degree of the generated graph, as
#'   described in the details section below. The default `dens.vertex=NULL` is
#'   equivalent to `dens.vertex=dens`, and similarly for `dens.inside` and
#'   `dens.across`.
#'
#' @details
#'
#'
#'
#'   **Details on the sparsity degree of the generated graph**
#'
#'   The argument `dens.vertex` specifies the proportion of coloured vertices
#'   among the \eqn{p} vertices. This is used if the string "vertex" is a
#'   component of `type` but not of `force.symm`. The string "vertex" not being
#'   a component of `type` is equivalent to `dens.vertex=0` whereas the string
#'   "vertex" being a component of both `type` and `force.symm` is equivalent to
#'   `dens.vertex=1`.
#'
#'   The argument `dens.inside` specifies the proportion of coloured symmetric
#'   inside-block edges among the \eqn{q(q-1)/2}, with \eqn{q=p/2}, inside-block
#'   edges. This is used if the string "inside.block.edge" is a component of
#'   `type`, otherwise it is equivalent to `dens.inside=0`. The overall density
#'   of inside-block edges is obtained by the sum of the densities of coloured
#'   and uncoloured inside-block edges. This is  a value between `dens.inside`
#'   and `dens.inside+dens`. Furthermore, it is exactly equal to `dens` if
#'   "inside.block.edge" is not a component of `type` and to `dens.inside` if
#'   "inside.block.edge" is a component of both `type` and `force.symm`.
#'
#'   The argument `dens.across` specifies the proportion of coloured symmetric across-block edges among the potentially coloured \eqn{q(q-1)/2} across-
#'   block edges. This is used if the string "across.block.edge" is a component
#'   of `type`, otherwise it is equivalent to `dens.across=0`. The overall
#'   density of across-block edges is obtained by the sum of the densities of
#'   coloured and uncoloured across-block edges. This is  a value between
#'   `dens.across` and `dens.across+dens`. Furthermore, it is exactly equal  to
#'   `dens` if "across.block.edge" is not a component of `type`.
#'
#'   The argument `dens` specifies the density of uncoloured edges. Note that
#'   the algorithm generates uncoloured edges first, which may be overwritten by
#'   coloured edges. For this reason the actual density of uncoloured edges is
#'   typically smaller than `dens`.
#'
#'   **Details on the generating process of the concentration matrix**
#'
#'   The concentration matrix is obtained by first generating a random Wishart
#'   matrix with matrix parameter `Sigma` and \eqn{p} degrees of freedom, which
#'   represents an initial unconstrained covariance matrix. This is inverted and
#'   adapted to a suitable coloured graph for paired data with sparsity degree
#'   according to the  `dens.xxx` arguments.
#'
#'
#' @return A list with the following components:
#'
#'  * `pdColG` a randomly generated a matrix representing a coloured graph for
#'  paired data on \eqn{p} vertices; see [`pdglasso-package`] for details.
#'
#' * `K` a randomly generated concentration matrix adapted to `pdColG`.
#'
#' * `sample.data` a randomly generated sample form a multivariate normal distribution
#'   with mean vector zero and concentration matrix `K`.
#'
#'
#' Note that the variable in \eqn{L} are are named `L1,...,Lq` and variables in
#' \eqn{R} are are named `R1,...,Rq` where  `Li` is homologous to  `Ri` for
#' every i=1,...,q.
#'
#'
#'
#' @examples
#'
#' # generates a pdRCON model on 10 variables in the form of a pdColG matrix
#'
#' set.seed(123)
#' pdRCON.model <- pdRCON.simulate(10, concent=FALSE, sample=FALSE, dens=0.25)$pdColG
#'
#' # generates a pdRCON model on 20 variables, a concentration matrix
#' # for this model and a sample of size 50
#' # all vertices are coloured and no coloured across-block edge is allowed
#'
#' set.seed(123)
#' GenMod <- pdRCON.simulate(20, type=c("v", "i"), force.symm=c("v"), sample.size=50, dens=0.20)
#'
#' @export
#'
pdRCON.simulate <- function(p,
                            concent.mat = TRUE,
                            sample = TRUE,
                            Sigma = NULL,
                            sample.size = NULL,
                            type = c("vertex", "inside.block.edge", "across.block.edge"),
                            force.symm = NULL,
                            dens = 0.1,
                            dens.vertex = NULL,
                            dens.inside = NULL,
                            dens.across = NULL){
  if((p %% 2) != 0) stop("The number of variables p must be even")
  if(concent.mat) {
    if(is.null(Sigma)) Sigma <- diag(1, p)
    if(p != nrow(Sigma)) stop("The dimension of Sigma is not pXp")
    S <- rWishart(1, df=p, Sigma=Sigma)[,,1]
    S.inv <- solve(S)
  }else{
    S.inv <- NULL
  }
  PDCG <- make.pdColG(p = p, K = S.inv, type = type, force.symm = force.symm, dens = dens,
                      dens.vertex = dens.vertex, dens.inside = dens.vertex, dens.across = dens.vertex)
  if(concent.mat){
    K <- pdRCON.mle(S, PDCG)
    dimnames(K) <- dimnames(PDCG)
    if(sample){
      if (is.null(sample.size)) sample.size <- 3*p
      sample.data <- MASS::mvrnorm(n = sample.size, mu=rep(0, p), Sigma=solve(K))
      colnames(sample.data) <- rownames(K)
      sample.data <- as.data.frame(sample.data)
    }else{
      sample.data <- NULL
    }
  }else{
    if(sample) warning("A random sample cannot be generated because concent.mat=FALSE")
    K <- NULL
    sample.data <- NULL
  }
  return(list(pdColG = PDCG, K = K, sample.data = sample.data))
}


#' Random generation of a pdColG matrix
#'
#' @param K the inverse of a variance matrix
#' @param p,type,force.symm,dens,dens.vertex,dens.inside,dens.across the same as
#'   in [`pdRCON.simulate`]
#'
#' @details
#'
#' This function is called by [`pdRCON.simulate`]
#'
#' @return this function returns a `pdColG` matrix, see [`pdglasso-package`] for
#'   details. If a matrix `K` is passed then this function produced a pdColG
#'   graph that tries to mimic the "zero" structure of `K` in the sense that
#'   missing edges correspond to the smallest entries of `K`. Similarly the
#'   symmetric structure is somehow read from `K`. See [`symm.structure.gen`]
#'   for details.
#'
#'   If `K=NULL` then the pdColG is generated randomly.
#'
#' @noRd
#'
make.pdColG <- function(p=NULL,
                        K=NULL,
                        type = c("vertex", "inside.block.edge", "across.block.edge"),
                        force.symm=NULL,
                        dens=0.1,
                        dens.vertex=NULL,
                        dens.inside=NULL,
                        dens.across=NULL){

  # initialization and checks
  if(is.null(dens.vertex)) dens.vertex <- dens
  if(is.null(dens.inside)) dens.inside <- dens
  if(is.null(dens.across)) dens.across <- dens
  dens.v <- c(dens, dens.vertex, dens.inside, dens.across)
  if (any(dens.v>1) | any(dens.v<0)) stop("all densities must be values in the interval [0; 1]")
  #
  out.make.a <- make.acronyms(type, force.symm)
  acr.type <- strsplit(out.make.a$acronym.of.type, split = "")[[1]]
  acr.force <- out.make.a$acronym.of.force
  #
  if(is.null(K)){
    q <- p/2

    # make overall graph
    G <- symm.structure.gen(p=p, dens=dens)

    # inside.block.edge symmetries
    if(any(acr.type=="I")){
      G.sym <- symm.structure.gen(p=q, dens=dens.inside)
      # remove coloured symmetric edges from overall graph
      G[1:q,1:q] <- G[1:q,1:q]*(1-G.sym)
      G[(q+1):p,(q+1):p] <- G[(q+1):p,(q+1):p]*(1-G.sym)
    }else{
      G.sym <- NULL
    }

    # across.block.edge symmetries
    if(any(acr.type=="A")){
      G.across <- symm.structure.gen(p=q, dens=dens.across)
      # remove coloured symmetric edges from overall graph
      G[1:q,(q+1):p] <- G[1:q,(q+1):p]*(1-G.across)
      G[1:q,(q+1):p] <- G[1:q,(q+1):p]*t(1-G.across)
    }else{
      G.across <- NULL
    }

    # vertex symmetries
    if(any(acr.type=="V")){
      m <- floor(dens.vertex*q)
      idx <- sample(1:q, size=m, replace = FALSE)
      v.symm <- rep(0, q)
      v.symm[idx] <- 1
      if(is.null(G.sym)){
        G.sym <- diag(v.symm)
      }else{
        diag(G.sym) <- v.symm
      }
    }
  }else{
    p <- nrow(K)
    q <- p/2

    # overall graph
    Kvec <- abs(K[lower.tri(K, diag=FALSE)])
    threshold <- quantile(Kvec, 1-dens)
    G <- (abs(K)>=threshold)
    G[lower.tri(G, diag=TRUE)] <- 0

    # inside.block.edge symmetries
    if(any(acr.type=="I")){
      G.sym <- symm.structure.gen(K[1:q, 1:q], K[(q+1):p, (q+1):p], dens=dens.inside)
      # remove coloured symmetric edges from overall graph
      G[1:q,1:q] <- G[1:q,1:q]*(1-G.sym)
      G[(q+1):p,(q+1):p] <- G[(q+1):p,(q+1):p]*(1-G.sym)
    }else{
      G.sym <- NULL
    }

    # across.block.edge symmetries
    if(any(acr.type=="A")){
      G.across <- symm.structure.gen(K[1:q, (q+1):p], t(K[1:q, (q+1):p]), dens=dens.across)
      # remove coloured symmetric edges from overall graph
      G[1:q,(q+1):p] <- G[1:q,(q+1):p]*(1-G.across)
      G[1:q,(q+1):p] <- G[1:q,(q+1):p]*t(1-G.across)
    }else{
      G.across <- NULL
    }

    # vertex symmetries
    if(any(acr.type=="V")){
      ab.vec <- abs(diag(K[1:q,1:q])-diag(K[(q+1):p, (q+1):p]))
      n.symmetries <- floor(dens.vertex*q)
      threshold <- sort(ab.vec)[max(q-n.symmetries, 1)]
      v.symm <- (ab.vec>=threshold)*1
      if(is.null(G.sym)){
        G.sym <- diag(v.symm)
      }else{
        diag(G.sym) <- v.symm
      }
    }
  }
  # force full symmetry if required
  if(!is.null(acr.force)){
    acr.force <- strsplit(acr.force, split = "")[[1]]
    if(any(acr.force=="I")){
      G[1:q,1:q] <- 0
      G[(q+1):p, (q+1):p] <- 0
    }
    if(any(acr.force=="A")){
      G[1:q,(q+1):p] <- 0
    }
    if(any(acr.force=="V")){
      diag(G.sym) <- 1
    }
  }
  PDCG <- G.merge(list(G = G, G.sym = G.sym, G.across = G.across))
  LR.lab <- c(paste("L", 1:(p/2), sep=""), paste("R", 1:(p/2), sep=""))
  dimnames(PDCG) <- list(LR.lab, LR.lab)
  return(PDCG)
}



# This function is called by [`make.pdColG`]
# Generates a symmetric structure
#
# A,B  = two (concentration) matrices with same dimension
#
# This function returns a binary (0/1) upper triangular matrix
# with ones correspond to nonzero symmetries obtained as follows:
#
# The entries (i.e. concentrations) in A are paired with those of B and the average is
# computed. The largest abs(averaged.value) pairs are eligible for symmetry
#
# dens = between 0 and 1 is the density of the resulting symmetry. It is  computed
#        as the total number of edges with symmetry (i.e. the number of pairs multiplied
#        by 2) divided by the total numer of edges that is q(q-1)/2 + q(q-1)/2
#


#' Generate a symmetric structure
#'
#' @param A,B two (concentration) matrices with same dimension
#' @param p number of variables
#' @param dens density of the symmetric structure
#'
#' @details
#'
#' This function is called by [`make.pdColG`]
#'
#'
#' @return
#'
#' This function returns a binary (0/1) upper triangular matrix with ones
#' correspond to nonzero symmetries obtained as follows:
#'
#' If `A` and  `B` are given then `p=nrow(A)`. The entries (i.e. concentrations)
#' in A are paired with those of B and the average is computed. The largest
#' abs(averaged.value) pairs are eligible for symmetry
#'
#' dens = between 0 and 1 is the density of the resulting symmetry. In the
#' resulting graph computed as the total number of edges with symmetry (i.e. the
#' number of pairs multiplied by 2) divided by the total number of edges that is
#' q(q-1)/2 + q(q-1)/2
#'
#' if `A` and  `B` are `NULL` then the 0/1 structure is generated randomly and
#' `p` gives the number of vertices.
#'
#' @noRd
#'
symm.structure.gen <- function(A=NULL, B=NULL, p=NULL, dens){
  if(is.null(p)){
    p <- nrow(A)
    M <- abs((A+B)/2)
    M[lower.tri(M, diag=TRUE)] <- 0
    ab.vec <- M[upper.tri(M, diag=FALSE)]
    n.symmetries <- floor(dens*(p*(p-1)/2))
    threshold <- sort(ab.vec)[max((p*(p-1)/2)-n.symmetries, 1)]
    G <- (M>=threshold)*1
  }else{
    G <- matrix(0, nrow=p, ncol=p)
    k <- (p*(p-1)/2)
    m <- floor(dens*k)
    idx <- sample(1:k, size=m, replace = FALSE)
    v <- rep(0, k)
    v[idx] <- 1
    G[upper.tri(G, diag=FALSE)] <- v
  }
  return(G)
}
