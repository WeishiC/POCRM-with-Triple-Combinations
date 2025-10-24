library(pocrm); library(dfcrm)  # for CRM-related functions
library(nnet)
library(combinat)               # for listing orderings
library(dplyr); library(tidyr)
library(ggplot2)

#
#  Simulate scenarios for three dimensional grids of any dimension
#
Adding <- function(nsims, MTC_list, labels, dim, TTL, L, subset=NULL, symmetric=FALSE) {
  # MTC_list: matrix, each row gives the location of the MTC
  # labels: vector of length=nrow(MTC_list), labels of the MTCs
  all.scen <- list()
  order_selected <- list()
  order_selected[[1]] <- 1:L
  if(is.vector(MTC_list)) MTC_list <- matrix(MTC_list, nrow=1)
  u <- 1
  for(c in 1:nrow(MTC_list)){
    # simulate nsims scenarios with MTC at specified place  
    all.scen[[c]] <- list()
    scen.store <- t(sapply(1:nsims, function(i) {
      y <- sim.scen(sim=1, MTC=MTC_list[c,], dim=dim, TTL=TTL)
      if(!is.null(subset)) y <- y[subset]
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
          filter(include==TRUE) %>% dplyr::select(n.scen) %>% unlist() %>% as.vector()
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
            filter(include==TRUE) %>% dplyr::select(-label, -n.scen, -include) %>% as.matrix()
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
  
  out <- list(OrderScen=all.scen.mat, order.adding=order.selected.mat, n.scen=n.scen, n=n)
  
  if(symmetric) {
    OrderRev <- function(scen_mat, TTL) {
      if(is.vector(scen_mat)) scen_mat <- matrix(scen_mat, nrow=1)
      L <- ncol(scen_mat)
      B <- B_rev <- list()
      for (i in 1:nrow(scen_mat)) {
        scen <- scen_mat[i,]
        B[[i]] <- sort(which(scen<TTL))
        B_rev[[i]] <- sort(L-which(scen>TTL)+1) %>% as.integer()
      }
      B <- unique(B); B_rev <- unique(B_rev)
      return(list(B=B, B_rev=B_rev))
    }
    
    SymSecn <- function(scen, B) {
      u <- numeric(length(scen))
      MTC <- L-which(scen==TTL)+1
      u[MTC] <- TTL
      u[B] <- sort(runif(length(B), 0, TTL))
      u[-c(B, MTC)] <- sort(runif(L-length(B)-1, TTL, 1))
      u
    }
    
    ord.scen <- list()
    ord.scen[[1]] <- matrix(all.scen.mat[1,],nrow=1); ord.scen[[L]] <- matrix(all.scen.mat[nrow(all.scen.mat),],nrow=1)
    
    extend_ord <- rep(NULL, L)
    for (l in 2:(L/2)) {
      begin <- all.scen.mat[(sum(n.scen[1:(l-1)])+1):sum(n.scen[1:l]),] %>% OrderRev(TTL=TTL)
      end <- all.scen.mat[(sum(n.scen[1:(L-l)])+1):sum(n.scen[1:(L-l+1)]),] %>% OrderRev(TTL=TTL)
      begin_scen <- all.scen.mat[(sum(n.scen[1:(l-1)])+1):sum(n.scen[1:l]),]
      end_scen <- all.scen.mat[(sum(n.scen[1:(L-l)])+1):sum(n.scen[1:(L-l+1)]),]
      if(length(end$B_rev)>length(begin$B)) {
        include <- sapply(1:length(end$B_rev), function(i) Position(function(x) identical(x, end$B_rev[[i]]), begin$B, nomatch = 0) > 0)
        extend_scen <- rep(NULL,L)
        for (i in which(include==FALSE)) {
          extend_scen <- rbind(extend_scen, SymSecn(end_scen[i,], end$B_rev[[i]]))
          extend_ord <- rbind(extend_ord, order(SymSecn(end_scen[i,], end$B_rev[[i]])))
        }
        ord.scen[[l]] <- rbind(begin_scen, extend_scen)
        ord.scen[[L-l+1]] <- end_scen
      }
      if(length(end$B_rev)==length(begin$B)) {
        ord.scen[[l]] <- begin_scen
        ord.scen[[L-l+1]] <- end_scen
      } 
      if(length(end$B_rev)<length(begin$B)) {
        include <- sapply(1:length(begin$B), function(i) Position(function(x) identical(x, begin$B[[i]]), end$B_rev, nomatch = 0) > 0)
        extend_scen <- rep(NULL,L)
        for (i in which(include==FALSE)) {
          extend_scen <- rbind(extend_scen, SymSecn(begin_scen[i,], begin$B_rev[[i]]))
          extend_ord <- rbind(extend_ord, order(SymSecn(begin_scen[i,], begin$B_rev[[i]])))
        }
        ord.scen[[l]] <- begin_scen
        ord.scen[[L-l+1]] <- rbind(end_scen, extend_scen)
      }
    }
    n_scen_corr <- sapply(1:length(ord.scen), function(i) nrow(ord.scen[[i]]))
    ScenAdding_corr <- matrix(nrow=sum(n_scen_corr), ncol=L) 
    u <- 0
    for (l in 1:L) {
      ScenAdding_corr[(u+1):(u+n_scen_corr[l]),] <- ord.scen[[l]]
      u <- u + n_scen_corr[l]
    }
    Order_corr <- rbind(order.selected.mat, extend_ord)
    out$n.scen.symmetry <- n_scen_corr
    out$order.symmetry <- Order_corr
    out$scen.symmetry <- ScenAdding_corr
    out$n.symmetry <- sum(n_scen_corr)
  }
  return(out)
}



sim.scen <- function(sim, MTC=NULL, dim, TTL=NULL){
  # Inputs:
  #--------------------------------------------------------
  # sim: integer, index of simulation
  # MTC: vector, the location of the MTC.
  # dim: vector, the dimension of the grid.
  # TTL: target toxicity level, value in [0,1].
  # -------------------------------------------------------
  dim1 <- dim[1]; dim2 <- dim[2]; dim3 <- dim[3]
  y <- array(dim=c(dim[1], dim[2], dim[3]))
  if(is.null(MTC)) {
    fix <- c(1,1,1)
    TTL <- 0.001
  } else {
    fix <- MTC
  }
  y[fix[1], fix[2], fix[3]] <- TTL
  if(fix[1]>1) y[1:(fix[1]-1),fix[2],fix[3]] <-sort(runif(fix[1]-1, max=TTL))           # MTC row (start : MTC-1)
  if(fix[1]<dim1) y[(fix[1]+1):dim1,fix[2],fix[3]] <-sort(runif(dim1-fix[1], min=TTL))  # MTC row (MTC+1 : end)
  if(fix[2]>1) y[fix[1],1:(fix[2]-1),fix[3]] <-sort(runif(fix[2]-1, max=TTL))           # MTC column (start : MTC-1)
  if(fix[2]<dim2) y[fix[1],(fix[2]+1):dim2,fix[3]] <-sort(runif(dim2-fix[2], min=TTL))  # MTC column (MTC+1 : end)
  if(fix[3]>1) y[fix[1],fix[2],1:(fix[3]-1)] <-sort(runif(fix[3]-1, max=TTL))           # MTC dim3 (start : MTC-1)
  if(fix[3]<dim3) y[fix[1],fix[2],(fix[3]+1):dim3] <-sort(runif(dim3-fix[3], min=TTL))  # MTC dim3 (MTC+1 : end)
  if(fix[1]>1) {
    for (i in rev(1:(fix[1]-1))) {
      if(fix[2]>1){
        for (j in rev(1:(fix[2]-1))) {
          y[i,j,fix[3]] <- runif(1, max=min(y[i, j+1, fix[3]], y[i+1, j, fix[3]]))    # (left, front) MTC plane
        }
      }
    }
    for (i in rev(1:(fix[1]-1))) {
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[i,j,fix[3]] <- runif(1, min=y[i,j-1,fix[3]], max=y[i+1,j,fix[3]])         # (left, back) MTC plane
        }
      }
    }
  }
  if(fix[1]<dim1) {
    for (i in (fix[1]+1):dim1) {
      if(fix[2]>1){
        for (j in rev(1:(fix[2]-1))) {
          y[i,j,fix[3]] <- runif(1, min=y[i-1,j,fix[3]], max=y[i, j+1, fix[3]])      # (right, front) MTC plane
        }
      }
    }
    for (i in (fix[1]+1):dim1) {
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[i,j,fix[3]] <- runif(1, min=max(y[i, j-1, fix[3]], y[i-1, j, fix[3]]))   # (right, back) MTC plane
        }
      }
    }
  }
  if(fix[3]<dim3) {
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          y[i,fix[2],index3] <- runif(1, min=y[i,fix[2],index3-1], max=y[i+1,fix[2],index3])  # (left, top) MTC plane
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          y[i,fix[2],index3] <- runif(1, min=max(y[i,fix[2],index3-1], y[i-1,fix[2],index3])) # (right, top) MTC plane
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[2]>1) {
        for (j in rev(1:(fix[2]-1))) {
          y[fix[1],j,index3] <- runif(1, min=y[fix[1],j,index3-1], max=y[fix[1],j+1,index3])  # (front, top) MTC plane
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[fix[1],j,index3] <- runif(1, min=max(y[fix[1],j,index3-1], y[fix[1],j-1,index3])) # (back, top) MTC plane
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=y[i,j,index3-1], max=min(y[i, j+1, index3], y[i+1, j, index3])) # (left, front, top)
            }
          }
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i,j-1,index3], y[i,j,index3-1]), max=y[i+1,j,index3])     # (left, back, top)
            }
          }
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=max(y[i-1,j,index3], y[i,j,index3-1]), max=y[i, j+1, index3])  # (right, front, top)
            }
          }
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3], y[i,j,index3-1]))    # (right, back, top)
            }
          }
        }
      }
    }
  }
  if(fix[3]>1) {
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          y[i,fix[2],index3] <- runif(1, max=min(y[i+1,fix[2],index3], y[i,fix[2],index3+1]))  # (left, bottom) MTC plane
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          y[i,fix[2],index3] <- runif(1, min=y[i-1,fix[2],index3], max=y[i,fix[2],index3+1])   # (right, bottom) MTC plane
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[2]>1) {
        for (j in rev(1:(fix[2]-1))) {
          y[fix[1],j,index3] <- runif(1, max=min(y[fix[1],j+1,index3], y[fix[1],j,index3+1]))  # (front, bottom) MTC plane
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[fix[1],j,index3] <- runif(1, min=y[fix[1],j-1,index3], max=y[fix[1],j,index3+1])   # (back, bottom) MTC plane
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, max=min(y[i, j+1, index3], y[i+1, j, index3], y[i,j,index3+1]))  # (left, front, bottom)
            }
          }
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=y[i,j-1,index3], max=min(y[i+1,j,index3], y[i,j,index3+1]))  # (left, back, bottom)
            }
          }
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=y[i-1,j,index3], max=min(y[i, j+1, index3], y[i,j,index3+1]))  # (right, front, bottom)
            }
          }
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3]), max=y[i,j,index3+1])  # (right, back, bottom)
            }
          }
        }
      }
    }
  }
  
  y <- as.vector(y)
  y
}
