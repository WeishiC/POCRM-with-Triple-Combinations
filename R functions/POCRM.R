######################################################
######################################################
##                                                  ##
##  Implementation and simulation of 2-stage POCRM  ##
##  Weishi Chen                                     ##
##  Last update: 17 Oct 2025                        ##
##                                                  ##
######################################################
######################################################
#
#  The minimum version of the two-stage POCRM implementation function
#  ordering selection, no overdose control.
#
pocrm.imp <- function(y.data, n.data, p.skel, ttr, cohortsize, ncohort, mprior.tox=NULL){
  nord.tox <- nrow(p.skel)                 # number of orderings
  n.dose <- length(y.data)                 # number of combinations
  
  # log-likelihood
  loglik <- function(a, p, y, n) {
    lik <- 0
    for (j in 1:length(p)) {
      pj <- p[j]^exp(a)
      lik <- lik + dbinom(y[j], size=n[j], prob=pj, log=TRUE)
    }
    return(lik)
  }
  
  est.tox <- c()
  marginal.tox <- rep(0, nord.tox)                # marginal likelihood of the data under each ordering
  for (m in 1:nord.tox) {
    try <- optim(par=0, fn=loglik, p = p.skel[m, ], y = y.data, n = n.data, method="Brent", lower=-5, upper=10, 
                 control=list(maxit=1000, fnscale=-1))
    est.tox[m] <- try$par
    if(try$convergence!=0) print("not converge")
    marginal.tox[m] <- loglik(a=est.tox[m], p=p.skel[m, ], y=y.data, n=n.data)
  }
  p.order <- numeric(nord.tox)
  for (m in 1:nord.tox) {
    p.order[m] <- mprior.tox[m] / (mprior.tox[m] + sum(exp(marginal.tox[-m]-marginal.tox[m])*mprior.tox[-m]))
  }
  # ordering selection
  ord.sel <- which.is.max(p.order)
  ptox.hat <- p.skel[ord.sel,]^exp(est.tox[m])
  # combination selection
  comb.best <- which.is.max(-abs(ptox.hat-ttr))
  
  return(list(comb.best=comb.best, ptox.hat=ptox.hat))
}

#
#  Simulate 1 trial
#
pocrm <- function(p0, p.skel, ttr, cohortsize, ncohort, x0, mprior.tox=NULL){
  if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel))
  nord.tox <- nrow(p.skel)                 # M
  n.dose <- length(p0)
  if(is.null(mprior.tox)) mprior.tox <- rep(1/nord.tox, nord.tox)  # prior probabilitiy of orderings
  y <- n <- rep(0, n.dose) 
  #  start the trail at a dose level pre-specified by start.comb
  comb.curr <- x0[1]
  ptox.hat <- numeric(n.dose)  # \hat{R}(d_j)
  comb.select <- rep(0, n.dose)
  y.tox <- numeric(ncohort)
  i <- 1      # index for number of cohorts
  x <- c()
  x[i] <- comb.curr
  
  # stage1: follow the path defined by x0
  stage1 <- c(x0, rep(n.dose, ncohort-length(x0)))
  while (i <= ncohort) {
    y.curr <- rbinom(1, cohortsize, p0[comb.curr])
    y.tox[i] <- y.curr
    y[comb.curr] <- y[comb.curr] + y.curr
    n[comb.curr] <- n[comb.curr] + cohortsize
    if (sum(y)==sum(n)) {
      # if all DLT=1, de-escalate if not already at d1
      comb.curr <- ifelse(comb.curr==1, comb.curr, comb.curr - 1)
    }
    else if (sum(y)==0) {
      # if all DLT=0, escalate follow the pre-defined path if not already at d_k
      comb.curr <- ifelse(comb.curr==n.dose, comb.curr, stage1[i + 1])
    }
    else {
      break
    }
    i <- i + 1
  }
  
  #  stage 2
  while (i <= ncohort) {
    y.curr <- rbinom(1, cohortsize, p0[comb.curr])
    y.tox[i] <- y.curr
    y[comb.curr] <- y[comb.curr] + y.curr
    n[comb.curr] <- n[comb.curr] + cohortsize
    
    # implement
    imp <- pocrm.imp(y.data=y, n.data=n, p.skel=p.skel, ttr=ttr, cohortsize=cohortsize,
                     ncohort=ncohort, mprior.tox=mprior.tox)
    comb.best <- imp$comb.best
    ptox.hat <- imp$ptox.hat
    
    comb.curr <- comb.best
    x[i] <- comb.curr
    i <- i + 1
  }
  RMSE <- sqrt(mean((p0-ptox.hat)^2))
  if(comb.curr!=0) comb.select[comb.curr] <- comb.select[comb.curr] + 1
  return(list(comb.select=comb.select, tox.data=y, n.tox=sum(y), pt.allocation=n, 
              ptox.hat=ptox.hat, x=x, RMSE=RMSE, p.order = p.order))
}

#
#  Simulate multiple trials
#
pocrm.sim <- function(n.sim, p0, p.skel, ttr, cohortsize, ncohort, x0, mprior.tox) {
  sim <- replicate(n.sim, pocrm(p0=p0, p.skel=p.skel, ttr=ttr, cohortsize=cohortsize, 
                                ncohort=ncohort, x0=x0, mprior.tox=mprior.tox))
  comb.select_mat <- tox.data_mat <- ptox.hat_mat <- matrix(nrow=n.sim, ncol=length(p0))
  n.tox_vec <- RMSE_vec <- n.back_vec <- n.tox.back_vec <- rep()
  for (r in 1:n.sim) {
    comb.select_mat[r,] <- sim[,r]$comb.select
    tox.data_mat[r,] <- sim[,r]$tox.data
    ptox.hat_mat[r,] <- sim[,r]$ptox.hat
    n.tox_vec[r] <- sim[,r]$n.tox
    RMSE_vec[r] <- sim[,r]$RMSE
  }
  return(list(comb.select = colMeans(comb.select_mat),
              comb.select.full = comb.select_mat,
              RMSE = mean(RMSE_vec),
              RMSE.full = RMSE_vec,
              n.tox = mean(n.tox_vec),
              p.hat = colMeans(ptox.hat_mat)))
}
