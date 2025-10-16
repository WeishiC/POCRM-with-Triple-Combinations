#
#  List out all possible orderings
#   Motivating example 1 (used in the paper) with 12 combinations
#
source("SimScen.R")
dim1 <- 3
dim2 <- 4
dim3 <- 2
nsims <- 10000
orders.store <- t(sapply(1:nsims, function(i) {
  y <- sim.scen(sim=1, dim=c(dim1,dim2,dim3))
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
