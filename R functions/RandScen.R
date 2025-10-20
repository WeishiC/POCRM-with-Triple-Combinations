############################################################################
############################################################################
##                                                                        ##
##  Simulate random scenarios according to Clertant and O’Quigley (2017)  ##
##  Weishi Chen                                                           ##
##  Last update: 20 Oct 2025                                              ##
##                                                                        ##
############################################################################
############################################################################
library(pocrm)
set.seed(1999)
load("148Orderings.RData")

L <- 12
TTL <- 0.25

RandScen.sim <- function(L, TTL, orderings) {
  u <- numeric(L)
  # position of MTC
  MTC <- sample(1:L, size=1)
  u[MTC] <- TTL
  # generate upper bound
  M <- rbeta(1, max(L-MTC, 0.5), 1)
  Bs <- TTL + (1-TTL)*M
  # simulate ordered values
  if(MTC>1) u[1:(MTC-1)] <- sort(runif(MTC-1, 0, TTL))
  if(MTC<L) u[(MTC+1):L] <- sort(runif(L-MTC, TTL, Bs))
  
  # randomly select an ordering
  Om <- orderings[sample(1:nrow(orderings), 1),]
  u <- getwm(matrix(Om, nrow=1), matrix(u, nrow=1))
  u <- as.vector(u)
  return(u)
}


n.scen <- 1e+4
RandScen <- t(replicate(n.scen, RandScen.sim(L=L, TTL=TTL, orderings=all.orderings)))
