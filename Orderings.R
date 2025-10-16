library(pocrm); library(dfcrm)  # for CRM-related functions
library(nnet)
library(combinat)               # for listing orderings
library(dplyr); library(tidyr)
library(ggplot2)

#  Motivating example 1 (used in the paper)
#  12 combinations

#
# Listing all 148 orderings
#
dim1 <- 3
dim2 <- 4
dim3 <- 2
nsims<-10000

orders.store <- t(sapply(1:nsims, function(sim){
  y <- array(dim=c(dim1, dim2, dim3))
  y[1,1,1] <- 0
  y[1,1,2:dim3] <- sort(runif(dim3-1, min=y[1,1,1]))
  y[1,2:dim2,1] <- sort(runif(dim2-1, min=y[1,1,1]))
  y[2:dim1,1,1] <- sort(runif(dim1-1, min=y[1,1,1]))
  for (i in 2:dim1) {
    for (j in 2:dim2) {
      y[i,j,1] <- runif(1, min=max(y[i, j-1, 1], y[i-1, j, 1]))
    }
  }
  for (index3 in 2:dim3) {
    for (j in 2:dim2) {
      y[1,j,index3] <- runif(1, min=max(y[1,j-1,index3], y[1,j,index3-1]))
    }
    for (i in 2:dim1) {
      y[i,1,index3] <- runif(1, min=max(y[i-1,1,index3], y[i,1,index3-1]))
    }
    for (i in 2:dim1) {
      for (j in 2:dim2) {
        y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3], y[i, j, index3-1]))
      }
    }
  }
  y <- as.vector(y)
  y <- y[c(1, 4, 5, 8, 9, 11, 12, 17, 20, 21, 23, 24)]
  order(as.vector(t(y)))
}))
all.orderings <- unique(orders.store)
M <- nrow(all.orderings)
L <- ncol(all.orderings)
colnames(all.orderings) <- paste0("d", 1:L)
all.orderings <- all.orderings %>% as_tibble() %>%
  arrange(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12)
all.orderings <- all.orderings %>% as.matrix()
