############################################################
############################################################
##                                                        ##
##  The "Refining" step in the Adding-Refining algorithm  ##
##  Weishi Chen                                           ##
##  Last update: 17 Oct 2025                              ##
##                                                        ##
############################################################
############################################################
Refining <- function(max.sim, all.scen, all.orderings, S, order.fix=NULL, TTL=0.25) {
  if(is.vector(all.orderings)) all.orderings <- matrix(all.orderings, nrow=1)
  if(is.vector(all.scen)) all.scen <- matrix(all.scen, nrow=1)
  # relabelling the combinations
  relabel <- function(ordering, label) {
    new.ordering <- numeric()
    for (i in 1:length(ordering)) {
      new.ordering[which(ordering==i)] <- label[i]
    }
    new.ordering
  }
  # ordering groups
  group.models <- function(u, MTC) {
    # u: k-vector, one complete ordering
    nu <- which(u==MTC)
    if(nu>1) {
      w <- sum(u[1:(nu-1)]>MTC) # number of more toxic dose ordered before MTC
    } else {
      w <- 0
    }
    if(nu<k) {
      n <- sum(u[(nu+1):k]<MTC) # number of less toxic dose ordered after MTC
    } else {
      n <- 0
    }
    c(w,n)
  }
  k <- ncol(all.orderings)
  M <- nrow(all.orderings)
  G <- nrow(all.scen)
  out <- numeric(S)
  # select orderings
  Continue <- TRUE
  sim <- 1
  while (Continue) {
    # if pre-specified orderings, include those
    if(is.null(order.fix)) {
      sel <- sample(1:M, S, replace = FALSE)
    } else {
      S1 <- length(order.fix)
      sel <- c(order.fix, sample((1:M)[-order.fix], S-S1, replace=FALSE))
    }
    # selected orderings
    ordering.sel <- all.orderings[sel,]
    if(is.vector(ordering.sel)) ordering.sel <- matrix(ordering.sel, nrow=1)
    # check at least one in the correct group under each scenario
    y <- rep(TRUE, G)
    consis_c <- c()
    for (c in 1:G) {
      new.label <- order(order(all.scen[c,]))
      ordering.new <- apply(ordering.sel, 1, relabel, label=new.label) %>% t()
      MTC_j <- apply(ordering.new, 1, group.models, MTC=which(sort(all.scen[c,])==TTL))
      rownames(MTC_j) <- c("w", "n")
      MTC_j <- MTC_j %>% t() %>% as_tibble() %>% mutate(model=sel, .before=w)
      consis_c[c] <- MTC_j %>% filter(w==0, n==0) %>% nrow()
      if(consis_c[c]==0) y[c] <- FALSE
    }
    if(all(y)==TRUE) {
      out <- sel
      Continue <- FALSE
    } 
    if(sim>=max.sim) {
      Continue <- FALSE
    } 
    sim <- sim + 1
    if((sim %% 50)==0) print(sim)
  }
  if(sum(out)==0) out <- NULL
  return(list(Orderings=all.orderings[out,], n.consis=sum(consis_c)))
}



#
#  This function can be used to (1) refining the orderings selected by the "Adding" step
#                               (2) obtain prior weights as discussed in Section 6 in the paper
#

# Performing the Refining step
# *OrderScen* and *OrderAdding* are outputs from the Adding step
# set *max.sim* to a large number
OrderRefine <- Refining(max.sim=10000, OrderScen, OrderAdding, S=9, order.fix=NULL, TTL=0.25)$Orderings

# Calculating prior weights
# *scen* are the scenarios of interest
# *orderings* are the specified orderings
# set *max.sim* and *S* to 1
n.consis <- c()
for (m in 1:nrow(orderings)) {
  n.consis[m] <- Refining(max.sim=1, scen, orderings[m,], S=1, order.fix=NULL, TTL=0.25)$n.consis
}
mprior.tox <- n.consis/sum(n.consis)
