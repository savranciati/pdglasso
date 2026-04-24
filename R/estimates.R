#' Selection and estimate of a pdRCON model according to eBIC
#'
#' Performs a sequence of calls to [`admm.pdglasso`] on a sequence of
#' values for `lambda1` and `lambda2`. More specifically, firstly a path of `lambda1` values, with 
#' `lambda2`=0, is considered to select a "best" `lambda1` value according to eBIC. 
#' Next, a path of `lambda2` values is considered, with `lambda1` fixed to its previously selected value, 
#' and the final model is selected according to eBIC.
#' The user may then call [`pdColG.get`] to obtain the Coloured Graph
#' for Paired Data (pdColG) representing the selected model.
#'
#' @inheritParams admm.pdglasso
#' @param n the sample size.
#' @param lams a 2x4 numeric matrix to specify the path of lambda values, if `NULL` default values are
#'   used; see the "Details" section below for additional information.
#' @param gamma.eBIC parameter of eBIC with `gamma.eBIC=0` corresponding to the classical BIC; see [`compute.eBIC`] for details..
#' @param progress a logical; if `TRUE` provides a visual update in the console about the grid search over `lambda1` and `lambda2`
#' @param mle a logical; if `TRUE`, compute eBIC via the MLE,
#' if `FALSE` the pdglasso estimator is used; see [`compute.eBIC`] for details.
#'
#' @details **Details on the specification of the argument `lams`**
#' 
#' The argument `lams` is a 2x4 matrix used to specified the grid of penalty terms.
#' The first row of `lams` refers to `lambda1` and second row to
#' `lambda2`. For each row of the `lams` matrix entries are: the minimum and maximum value of the path (columns 1 and 2); 
#'  the number of points of the path (column 3); `0` for linear spacing and  
#' `1` logarithmic spacing (column 4). If `lams=NULL` the default values are computed as follows: maximum lambda values are computed from [`lams.max`]
#' and paths of 20 points are used form max.lams/20 to max.lams, with logarithm spacing. 
#' 
#'
#' @return A list with the following components:
#'
#' * `model` selected model; object of class `ADMMoutput` output of [`admm.pdglasso`].
#' 
#' * `lambda.grid` the grid of values used for `lambda1` and `lambda2.`
#' 
#' * `best.lambdas` the selected values of `lambda1` and `lambda2` according to eBIC.
#' 
#' * `l1.path` a matrix containing the path values for `lambda1` as well as relevant quantities of the eBIC.
#' 
#' * `l2.path` a matrix containing the path values for `lambda2` as well as relevant quantities of the eBIC.
#' 
#' * `time.exec` total execution time for the called function.
#'
#' A warning is produced if at least one run of the algorithm for the grid
#' searches has resulted in non-convergence (status can be checked by inspecting
#' `l1.path` and `l2.path`).
#' @export
#'
#' @examples
#' S <- cov(toy_data$sample.data)
#' sel.mod <- pdRCON.select(S, n=60)
#' sel.mod$l1.path
#' sel.mod$l2.path
#' plot(sel.mod$model)
#' pdColG.get(sel.mod$model)

pdRCON.select <- function(S,
                       n,
                       lams        = NULL,
                       gamma.eBIC  = 0.5,
                       type        = c("vertex", "inside.block.edge", "across.block.edge"),
                       force.symm  = NULL,
                       X.init      = NULL,
                       rho1        = 1,
                       rho2        = 1,
                       varying.rho1= TRUE,
                       varying.rho2= TRUE,
                       max_iter    = 5000,
                       eps.abs     = 1e-08,
                       eps.rel     = 1e-08,
                       verbose     = FALSE,
                       progress    = TRUE,
                       print.type  = TRUE,
                       mle.estimate = TRUE){
  start.time <- Sys.time()
  ## Max values for lambda_1 and lambda_2 according to theorems; only needed for the grid search
  
  if(is.null(lams)){
    lams <- matrix(0,2,4)
    lams[,2] <- lams.max(S)
    lams[,1] <- lams[,2]/20
    lams[,3] <- c(20,20)
    lams[,4] <- c(1,1)
  }
  rownames(lams) <- c("l1","l2")
  colnames(lams) <- c("min","max","n.pts","log.spacing")
  
  ## Prepare temp objects
  eBIC.l1 <-  matrix(0,lams[1,3],4)
  eBIC.l2 <-  matrix(0,lams[2,3],4)
  
  ### First grid search for lambda_1, with lambda_2=0
  if(lams[1,4]==1) l1.vec <- exp(seq(log(lams[1,1]),log(lams[1,2]), length.out=lams[1,3])) else l1.vec <- seq(lams[1,1], lams[1,2], length.out=lams[1,3])
  l1.vec <- sort(l1.vec, decreasing=TRUE)
  for(i in 1:lams[1,3]){
    if(progress==TRUE) cat("Searching over lambda1 grid (",i,"/",lams[1,3],").\n", sep="")
    mod.out <- admm.pdglasso(S,
                             lambda1=l1.vec[i],
                             lambda2=0,
                             type=type,
                             force.symm=force.symm,
                             X.init=X.init,
                             rho1=rho1,
                             rho2=rho2,
                             varying.rho1=varying.rho1,
                             varying.rho2=varying.rho2,
                             max_iter=max_iter,
                             eps.abs=eps.abs,
                             eps.rel=eps.rel,
                             verbose=FALSE,
                             print.type=FALSE)
    eBIC.l1[i,1:3] <- compute.eBIC(S=S,
                                   admm.out=mod.out,
                                   n=n,
                                   gamma.eBIC=gamma.eBIC,
                                   mle = mle.estimate)
    eBIC.l1[i,4] <- mod.out$internal.par$converged+0
    if(eBIC.l1[i,4]==0) cat("Convergence not achieved for this value of lambda1! \n")
  }
  best.l1 <- l1.vec[which.min(eBIC.l1[,1])]
  if(length(best.l1)==0) stop("Grid search of lambda1 failed!")
  
  if(progress==TRUE) cat("--- \n", sep="")
  
  ### Second grid search for lambda_2, with lambda_1=best.l1
  if(lams[2,4]==1) l2.vec <- exp(seq(log(lams[2,1]),log(lams[2,2]), length.out=lams[2,3])) else l2.vec <-  seq(lams[2,1], lams[2,2], length.out=lams[2,3])
  l2.vec <- sort(l2.vec, decreasing=TRUE)
  for(i in 1:lams[2,3]){
    if(progress==TRUE) cat("Searching over lambda2 grid (",i,"/",lams[2,3],").\n", sep="")
    mod.out <- admm.pdglasso(S,
                             lambda1=best.l1,
                             lambda2=l2.vec[i],
                             type=type,
                             force.symm=force.symm,
                             X.init=X.init,
                             rho1=rho1,
                             rho2=rho2,
                             varying.rho1=varying.rho1,
                             varying.rho2=varying.rho2,
                             max_iter=max_iter,
                             eps.abs=eps.abs,
                             eps.rel=eps.rel,
                             verbose=FALSE,
                             print.type=FALSE)
    eBIC.l2[i,1:3] <- compute.eBIC(S=S,
                                   admm.out=mod.out,
                                   n=n,
                                   gamma.eBIC=gamma.eBIC,
                                   mle = mle.estimate)
    eBIC.l2[i,4] <- mod.out$internal.par$converged+0
    if(eBIC.l2[i,4]==0) cat("Convergence not achieved for this value of lambda2! \n")
  }
  ### adding eBIC value/results and l2 value to path
  ### already estimated from the first grid.search where lam2=0
  eBIC.l2 <- rbind(eBIC.l2, eBIC.l1[which.min(eBIC.l1[,1]),])
  l2.vec <- c(l2.vec,0)
  best.l2 <- l2.vec[which.min(eBIC.l2[,1])]
  if(length(best.l2)==0) stop("Grid search of lambda2 failed!")
  
  ## Fit final model
  mod.out <- admm.pdglasso(S,
                           lambda1=best.l1,
                           lambda2=best.l2,
                           type=type,
                           force.symm=force.symm,
                           X.init=X.init,
                           rho1=rho1,
                           rho2=rho2,
                           varying.rho1=varying.rho1,
                           varying.rho2=varying.rho2,
                           max_iter=max_iter,
                           eps.abs=eps.abs,
                           eps.rel=eps.rel,
                           verbose=verbose,
                           print.type=print.type)
  
  l1.path=cbind(l1.vec,eBIC.l1)
  l2.path=cbind(l2.vec,eBIC.l2)
  colnames(l1.path) <- c("lambda1.grid", "eBIC     ","  log-Likelihood  ","num. of params", "converged (1=TRUE)")
  colnames(l2.path) <- c("lambda2.grid", "eBIC     ","  log-Likelihood  ","num. of params", "converged (1=TRUE)")
  
  time.exec <- Sys.time()-start.time
  return(list(model=mod.out,
              best.lambdas=c(best.l1,best.l2),
              lambda.grid=lams,
              l1.path=l1.path,
              l2.path=l2.path,
              time.exec=time.exec))
  
}




#' Maximum likelihood estimate of a pdRCON model
#'
#' Computes the maximum likelihood estimate of the concentration matrix of a pdRCON model.
#'
#' @param S a sample covariance matrix with the block structure described in [`pdglasso-package`].
#' @param pdColG a pdColG matrix representing a coloured graph for paired data; see [`pdglasso-package`] for details.
#' @inheritParams admm.pdglasso
#'
#' @return Either a matrix, that is the maximum likelihood estimate of
#'   \eqn{K=\Sigma^{-1}} under the pdRCON model represented by `pdColG`, or
#'   `NULL` if the maximum likelihood estimate does not exist.
#'
#' @details If the sample covariance matrix is not full-rank, then it is
#'   possible that the maximum likelihood estimate does not exist. The maximum
#'   likelihood estimate is computed by running the function [`admm.pdglasso`]
#'   with suitable penalties and, if it does not exist, then the ADMM algorithm
#'   fails to converge, a warning is produced and a `NULL` is returned.
#'
#' @export
#'
#' @examples
#' S <- var(toy_data$sample.data)
#' K.hat <- pdRCON.mle(S, toy_data$pdColG)

pdRCON.mle <- function(S, pdColG,
                       eps.rel=1e-6,
                       eps.abs=1e-6,
                       max_iter=5000){
  # make vector lambda1
  lambda1 <- (mat2vec(pdColG)==0)
  lambda1[lambda1] <- Inf
  
  # make vector lambda2
  p <- nrow(pdColG)
  q <- p/2
  m2v <- c(diag(pdColG[1:q,1:q]), half_vec(pdColG[1:q,1:q]), half_vec(pdColG[1:q,(q+1):p]))
  lambda2 <- (m2v==2)
  lambda2[lambda2] <- Inf
  
  # run SGL algorithm
  out.admm <- admm.pdglasso(S, lambda1 = lambda1,
                            lambda2 = lambda2,
                            type=c("V","I","A"),
                            eps.rel = eps.rel,
                            eps.abs = eps.abs,
                            max_iter= max_iter,
                            print.type=FALSE)
  K.hat <- NULL
  if(out.admm$internal.par$converged) K.hat = out.admm$X
  return(K.hat)
}




#' Check if a concentration matrix satisfies the likelihood equations
#'
#' Checks if a matrix is the maximum likelihood estimate of the concentration
#' matrix of a pdRCON model represented by the colored graph for paired data
#' `pdColG`, computed from the covariance matrix `S`.
#'
#' @param K.mle candidate mle of \eqn{K} to be checked.
#' @param pdColG matrix representing a colored graph for paired data;
#' see [`pdglasso-package`] for details.
#' @param S a sample covariance matrix with the block structure described
#' in [`pdglasso-package`].
#' @param tollerance threshold to check whether a value is equal to zero.
#' @param print.checks a logical (default `TRUE`).
#'
#' @return a logical equal to `TRUE` if `K.mle` satisfies the likelihood equations
#' and `FALSE` otherwise. If `print.checks = TRUE` the values used for the
#' checking are printed.
#'
#' @details
#' In order to check if `K.mle` is the maximum likelihood estimate, computed
#' from `S`, under the model represented by `pdColG` the following quantities
#' are considered:
#'
#' * the values of `K.mle` that are expected to be zero;
#'
#' * the differences between the entries of `K.mle` which are expected to be equal;
#'
#' * the differences between the relevant sums of the entries of `solve(K.mle)`
#'   and `S` which are expected to be equal, so as to satisfy the
#'   likelihood equations.
#'
#' * all the absolute quantities of the previous three points are divided
#'   by the \eqn{\ell_2} norm of the non-zero entries of `K.mle` so as to
#'   obtain relative quantities.
#'
#'   Then if all the above quantities are smaller than `toll` the check is
#'   considered successful and a `TRUE` is returned.
#' 
#' @noRd
## @examples
## S <- var(toy_data$sample.data)
## K.hat <- pdRCON.mle(S, toy_data$pdColG)
## is.pdRCON.mle(K.hat, toy_data$pdColG, S)
##
is.pdRCON.mle <- function(K.mle, pdColG, S, toll=1e-8, print.checks=TRUE){
  S.mle <- solve(K.mle)
  G <- pdColG
  p <- dim(G)[1]
  q <- p/2
  l <- 1:q
  r <- (q+1):p
  #
  nz.ms <- sqrt(mean(K.mle[G!=0]^2))
  # missing edges
  Bz <- (G==0)
  #
  # uncoloured vertices and edges
  Bu <- (G==1)
  #
  # coloured vertices and inside-block coloured edges
  Bci <- (G[l,l]==2)
  #
  # across-block coloured edges
  Bca <- (G[l,r]==2)
  #
  # missing edges: check zero concentrations
  if(any(Bz)){
    am <- max(abs(K.mle[G==0]))
  }else{
    am <- 0
  }
  v  <- c()
  sK <- c()
  #
  # coloured vertices and inside-block coloured edges:
  # check equality constraints in K.mle
  # check likelihood equations in S.mle
  if(any(Bci)){
    sK[1] <- max(abs((K.mle[l,l]-K.mle[r,r])[Bci]))
    v[1]  <- max(abs((S.mle[l,l]+S.mle[r,r]-S[l,l]-S[r,r])[Bci]))
  }else{
    sK[1] <- 0
    v[1]  <- 0
  }
  ##
  # across-block coloured edges:
  # check equality constraints in K.mle
  # check likelihood equations in S.mle
  if(any(Bca)){
    sK[2] <- max(abs((K.mle[l,r]-K.mle[r,l])[Bca]))
    v[2]  <- max(abs((S.mle[l,r]+S.mle[r,l]-S[l,r]-S[r,l])[Bca]))
  }else{
    sK[2] <- 0
    v[2]  <- 0
  }
  #
  # uncoloured vertices and edges: check likelihood equations
  if(any(Bu)){
    v[3] <- max(abs((S.mle-S)[G==1]))
  }else{
    v[3] <- 0
  }
  if(print.checks){
    cat("\nABSOLUTE QUANTITIES:\n")
    cat("Max zero concentration           : ", am, "\n")
    cat("Max diff. of equal concentrations: ", max(sK), "\n")
    cat("Max error in likelihood equations: ", max(v), "\n")
    #
    cat("\nRELATIVE QUANTITIES:\n")
    cat("Max zero concentration           : ", am/nz.ms, "\n")
    cat("Max diff. of equal concentrations: ", max(sK)/nz.ms, "\n")
    cat("Max error in likelihood equations: ", max(v)/nz.ms, "\n\n")
  }
  max.all <- max(c(am, max(sK), max(v), am/nz.ms, max(sK)/nz.ms, max(v)/nz.ms))
  
  if(max.all < toll){
    check.result <- TRUE
  }else{
    check.result <- FALSE
  }
  return(check.result)
}




#' Extended Bayesian Information Criterion (eBIC) for pdRCON models
#'
#' This function computes the value of the eBIC for a pdRCON model  
#' from the output of [`admm.pdglasso`]; see Eq.1 of Feygel and Drton (2010).
#'
#' @param S a sample covariance matrix with the block structure described in [`pdglasso-package`].
#' @param admm.out an object of class `ADMMoutput`, such as the output of a call to [`admm.pdglasso`]; see also [`pdRCON.select`].
#' @param n the sample size.
#' @param gamma.eBIC a parameter governing the magnitude of the penalization term inside the criterion; 
#' it ranges from 0 to 1, where 0 makes the eBIC equivalent to BIC, with 0.5 being the value 
#' suggested by Feygel and Drton (2010).
#' @param mle a logical; if `TRUE`, the MLE are used otherwise 
#'    the pdglasso estimator is used; see the "Details" section below. 
#'   
#' @details **Details on the specification of the argument `mle`**
#' 
#' The extended Bayesian Information Criterion implemented in this function requires that the parameters are estimated by maximum likelihood. 
#' However, the computation of MLEs may considerably slow down the procedure and, more seriously, in the high-dimensional setting   
#' MLEs may even not exist. For this reason, we provide the option that the eBIC is, naively, computed by using the pdglasso estimate.
#' 
#' @return A vector containing three elements:
#' 
#' * the value of the eBIC,
#' 
#' * the log-likelihood,
#' 
#' * the number of parameters.
#' 
#' @export
#'
#' @references Foygel, R. and Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. *Advances in neural information processing systems*, 23. 
#' @examples
#' S <- cov(toy_data$sample.data)
#' admm.out <- admm.pdglasso(S, lambda1=4, lambda2=0.7)
#' compute.eBIC(S, admm.out, n=60, gamma.eBIC=0.5)
#' 
#' compute.eBIC(S, admm.out, n=60, gamma.eBIC=0.5, mle=FALSE)

compute.eBIC <- function(S, admm.out,n,
                         gamma.eBIC=0.5,
                         mle = TRUE){
  G <- pdColG.get(admm.out, model.dimension=TRUE)
  
  if(mle){
    K <- pdRCON.mle(S,
                    pdColG = G$pdColG,
                    eps.rel = admm.out$internal.par$eps.rel,
                    eps.abs = admm.out$internal.par$eps.abs,
                    max_iter = admm.out$internal.par$max_iter)
  }else{
    K <- admm.out$X
  }
  
  
  if(is.null(K)){
    out.vec <- rep(NA,3)
  }else{
    S <- S*(n-1)/n
    p <- dim(S)[1]
    n.par <- G$n.par
    #log.lik <- log(det(K))-sum(S*K)
    log.lik <- as.numeric(determinant(K, logarithm=TRUE)$modulus) - sum(S*K)
    
    eBIC <- -n*log.lik+log(n)*n.par+4*n.par*gamma.eBIC*log(p)
    #### corrected with n/2 instead of 1/2 in front of the loglik, so -2*(n/2)*loglik=-n*loglik
    out.vec <- c(eBIC,log.lik,n.par)
  }
  names(out.vec) <- c("eBIC     ","  log-Likelihood  ","num. of params")
  return(out.vec)
}



#' pdColG matrix from the output of a call to [`admm.pdglasso`]
#'
#' This function returns the pdColG matrix representing the pdRCON model 
#' resulting from a call to [`admm.pdglasso`]; see also [`pdRCON.select`]. 
#'
#' @param admm.out an object of class `ADMMoutput`.
#' @param th1,th2 two positive scalars, the thresholds to identify edges  
#' (`th1`) and coloured edges (`th2`) in the graph, if `NULL` 
#' a default value is used; as describe in the "Details" section below. .
#' @param model.dimension a logical, if `TRUE` the number of parameters of the model is returned.
#'
#' @details **Details of the specification of the argument `th1` and `th2`**
#' 
#' Due to finite-precision arithmetic in computing, the fitted concentration matrix obtained from the ADMM has no exact zeros and no exact identical 
#' pairs of concentration values. Hence, `th1` and `th2` provide tollerances for values and differences, respectively, to be regarded as close enough to 
#' zero. If  `NULL` tollerances are automatically derived from tollerance values used to check convergence of the ADMM. 
#' The function [`plot.ADMMoutput`] allows to visualize the role played by the default threshold as well as to specify
#' facilitate the specification of suitable values for the personalized thresholds `th1` and `th2`.
#'
#' @return Either an object of class `pdColG`, if `model.dimension=FALSE`, or a list with the following components if `model.dimension=TRUE`:
#'
#' * `pdColG` an object of class `pdColG`, that is a matrix representing a coloured graph for paired data.
#'
#' * `n.par` the number of parameters of the pdRCON model.
#' 
#' @export
#'
#' @examples
#'
#' S <- cov(toy_data$sample.data)
#' admm.out <- admm.pdglasso(S, lambda1=4, lambda2=0.7)
#' pdColG.get(mod.out)

pdColG.get <- function(admm.out,
                       th1 = NULL,
                       th2 = NULL,
                       model.dimension = FALSE
                       ){
  # Prepare output object
  # Store passed acronyms used
  acronyms <- admm.out$acronyms$acronym.of.type
  # Store passed estimated concentration matrix
  X <- admm.out$X
  p <- dim(X)[1]
  q <- p/2
  
  if(is.null(th1)) th1 <- admm.out$internal.par$eps.rel*10
  if(is.null(th2)) th2 <- admm.out$internal.par$eps.rel*10
  
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
  if(grepl("I",acronyms,fixed=TRUE)){
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
  G <- mat_graph+mat_sym
  dimnames(G) <- dimnames(admm.out$X)
  G <- as.pdColG(G)
  #
  if(model.dimension){
    # computing degrees of freedom (+p before diving because math_graph has 1s on diagonal)
    tot.dof  <- (sum(mat_graph)+p)/2
    # number of pairs of symmetric off-diagonal concentrations
    nsym_offdiag  <- sum(LL.block(mat_sym)[upper.tri(LL.block(mat_sym),diag=F)]) +
    sum(across.block(mat_sym)[upper.tri(across.block(mat_sym),diag=F)])
    # number of pairs of symmetric diagonal concentrations
    nsym_diag <- sum(diag(mat_sym))/2
    n.par <- tot.dof - nsym_offdiag - nsym_diag
    out <- list()
    out$pdColG <- G
    out$n.par <- n.par
  }else{
    out <- G
  }
  return(out)
}








