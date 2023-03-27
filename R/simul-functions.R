require(MASS)


#' Title
#'
#' @param p number of variables.
#' @param concent.mat a logical (default `TRUE`) indicating whether a
#'   concentration matrix should be generated.
#' @param sample a logical (default `TRUE`) indicating whether a sample form a
#'   normal distribution with zero mean vector and the generated concentration
#'   matrix should be generated.
#' @param Sigma x
#' @param sample.size sample size with default value equal to \eqn{3p}.
#' @param type x
#' @param force.symm x
#' @param dens x
#' @param dens.vertex x
#' @param dens.inside x
#' @param dens.across x
#'
#' @return this function
#' @export
#'
#' @examples
#' #
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
    K <- pdColG.mle(S, PDCG)
    LR.lab <- c(paste("L", 1:(p/2), sep=""), paste("R", 1:(p/2), sep=""))
    dimnames(K) <- list(LR.lab, LR.lab)
    if(sample){
      if (is.null(sample.size)) sample.size <- 3*p
      sample.data <- MASS::mvrnorm(n = sample.size, mu=rep(0, p), Sigma=solve(K))
      colnames(sample.data) <- LR.lab
      sample.data <- as.data.frame(sample.data)
    }else{
      sample.data <- NULL
    }
  }else{
    if(sample)warning("concent.mat=FALSE no data are generated")
    K <- NULL
    sample.data <- NULL
  }
  return(list(pdColG = PDCG, K = K, sample.data = sample.data))
}


#' Title
#'
#' @param p x
#' @param K x
#' @param type x
#' @param force.symm x
#' @param dens x
#' @param dens.vertex x
#' @param dens.inside x
#' @param dens.across x
#'
#' @return this function
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
  return(PDCG)
}



# This function is called by generate.pdColG()
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


#' Title
#'
#' @param A x
#' @param B x
#' @param p x
#' @param dens x
#'
#' @return this function
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
