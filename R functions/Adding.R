#
#  The Adding step of the Adding-Refining algorithm
#
Adding <- function(nsims, MTC_list, labels, dim, TTL, L=12) {
  # MTC_list: matrix, each row gives the location of the MTC
  # labels: vector of length=nrow(MTC_list), labels of the MTCs
  all.scen <- list()
  order_selected <- list()
  order_selected[[1]] <- 1:L
  u <- 1
  for(c in 1:nrow(MTC_list)){
    print(c)
    # simulate nsims scenarios with MTC at specified place  
    all.scen[[c]] <- list()
    scen.store <- t(sapply(1:nsims, function(i) {
      y <- sim.scen(sim=1, MTC=MTC_list[c,], dim=dim, TTL=TTL)
      y <- y[c(1, 4, 5, 8, 9, 11, 12, 17, 20, 21, 23, 24)]
      y
    }))
    colnames(scen.store) <- paste0("d", 1:L)
    # labels of the MTCs
    label.store <- apply(scen.store, 1, function(u) which(sort(u)==TTL))
    label.unique <- unique(label.store)
    scen_curr <- scen.store %>% as_tibble() %>% mutate(label=label.store)
    scen.index <- 1
    for (i in label.unique) {
      # simulated scenarios with the specified label
      scen_temp <- scen_curr %>% mutate(n.scen=1:nsims) %>% filter(label==i) %>% as.matrix()
      # needed scenario, i.e. one for each below MTC set
      scen.needed <- t(apply(matrix(scen_temp[,1:L], ncol=L), 1, function(u) order(u)))[,1:i] %>%
        unique() 
      if(is.vector(scen.needed)) scen.needed <- matrix(scen.needed, nrow=1)
      # the order among the below MTC set doesn't matter
      if(i>1) {
        scen.needed <- cbind(t(apply(matrix(scen.needed[,1:(i-1)], ncol=(i-1)), 1, sort)) %>% unique(), labels[c])
        colnames(scen.needed) <- NULL
      }
      # add one order-scenario for each below MTC set
      scen_include <- numeric()
      for (j in 1:nrow(scen.needed)) {
        scen.orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed[j,])))
        scen.include_temp <- scen_temp %>% as_tibble() %>% mutate(include=scen.orders.needed) %>%
          filter(include==TRUE) %>% select(n.scen) %>% unlist() %>% as.vector()
        scen_include[j] <- scen.include_temp[1]
      }
      all.scen[[c]][[scen.index]] <- scen_curr[scen_include,1:L]
      scen.index <- scen.index + 1
      # under each below MTC set
      for (j in 1:nrow(scen.needed)) {
        # test if the selected orderings has already included one from correct group
        consis <- mapply(function(x, y){
          isTRUE(all.equal(x[1:i], y[1:i]))
        }, order_selected, list(scen.needed))
        # if not, add one from the correct group
        if(sum(consis)==0) {
          u <- u+1
          orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed[j,])))
          include_temp <- scen_temp %>% as_tibble() %>% mutate(include=orders.needed) %>%
            filter(include==TRUE) %>% select(-label, -n.scen, -include) %>% as.matrix()
          order_selected[[u]] <- order(include_temp[1,1:L])
        }
      }
    }
  }
  # collect results
  # orderings added
  order.selected.mat <- order_selected %>% unlist() %>% matrix(byrow=TRUE, ncol=L)
  # order scenarios
  all.scen.temp <- all.scen %>% unlist(recursive = FALSE)
  all.scen.mat <- matrix(NA, nrow=1, ncol=L)
  colnames(all.scen.mat) <- paste0("d", 1:L)
  for (i in 1:length(all.scen.temp)) {
    all.scen.mat <- rbind(all.scen.mat, all.scen.temp[[i]])
  }
  all.scen.mat <- all.scen.mat[-1,]
  all.scen.mat <- all.scen.mat %>% as.matrix()
  colnames(all.scen.mat) <- NULL
  # number of order-scenario by MTC location
  n.scen <- sapply(all.scen, function(list){
    u <- 0
    for (i in 1:length(list)) {
      u <- u + nrow(list[[i]])
    }
    u
  })
  # total number of order-scenario
  n <- sum(n.scen)
  return(list(OrderScen=all.scen.mat, order.adding=order.selected.mat, n.scen=n.scen, n=n))
}




#
# For the 12 combination example
#
#MTC_list <- matrix(c(1, 1, 1,
#                     1, 2, 1, 
#                     2, 2, 1, 
#                     2, 3, 1, 
#                     3, 3, 1,
#                     2, 4, 1,
#                     3, 4, 1,
#                     2, 2, 2, 
#                     2, 3, 2,
#                     3, 3, 2, 
#                     2, 4, 2, 
#                     3, 4, 2), byrow=TRUE, ncol=3)
#labels <- 1:12
#dim <- c(3, 4, 2)
#TTL <- 0.25
#L <- 12
#nsims <- 8000

# the order-scenarios stored in *OrderScen*
# the orderings included by the adding step stored in *OrderAdding*
#adding.out <- Adding(nsims=8000, MTC_list, labels=1:12, dim, TTL=0.25, L=12)
#adding.out$n.scen
#adding.out$n
#OrderScen <- adding.out$OrderScen
#OrderAdding <- adding.out$order.adding

                                 
