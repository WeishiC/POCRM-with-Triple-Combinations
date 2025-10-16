library(pocrm); library(dfcrm)  # for CRM-related functions
library(nnet)
library(combinat)               # for listing orderings
library(dplyr); library(tidyr)
library(ggplot2)

#
#  Simulate scenarios for three dimensional grids of any dimension
#
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
  if(fix[1]>1) y[1:(fix[1]-1),fix[2],fix[3]] <-sort(runif(fix[1]-1, max=TTL))
  if(fix[1]<dim1) y[(fix[1]+1):dim1,fix[2],fix[3]] <-sort(runif(dim1-fix[1], min=TTL))
  if(fix[2]>1) y[fix[1],1:(fix[2]-1),fix[3]] <-sort(runif(fix[2]-1, max=TTL))
  if(fix[2]<dim2) y[fix[1],(fix[2]+1):dim2,fix[3]] <-sort(runif(dim2-fix[2], min=TTL))
  if(fix[3]>1) y[fix[1],1:fix[2],1:(fix[3]-1)] <-sort(runif(fix[3]-1, max=TTL))
  if(fix[3]<dim3) y[fix[1],fix[2],(fix[3]+1):dim3] <-sort(runif(dim3-fix[3], min=TTL))
  if(fix[3]<dim3) {
    if(fix[1]>1) {
      for (i in rev(1:(fix[1]-1))) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,fix[3]] <- runif(1, max=min(y[i, j+1, fix[3]], y[i+1, j, fix[3]]))
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,fix[3]] <- runif(1, min=y[i,j-1,fix[3]], max=y[i+1,j,fix[3]])
          }
        }
      }
    }
    if(fix[1]<dim1) {
      for (i in (fix[1]+1):dim1) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,fix[3]] <- runif(1, min=y[i-1,j,fix[3]], max=y[i, j+1, fix[3]])
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,fix[3]] <- runif(1, min=max(y[i, j-1, fix[3]], y[i-1, j, fix[3]]))
          }
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          y[i,fix[2],index3] <- runif(1, min=y[i,fix[2],index3-1], max=y[i+1,fix[2],index3])
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          y[i,fix[2],index3] <- runif(1, min=max(y[i,fix[2],index3-1], y[i-1,fix[2],index3]))
        }
      }
      if(fix[2]>1) {
        for (j in rev(1:(fix[2]-1))) {
          y[fix[1],j,index3] <- runif(1, min=y[fix[1],j,index3-1], max=y[fix[1],j+1,index3])
        }
      }
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[fix[1],j,index3] <- runif(1, min=max(y[fix[1],j,index3-1], y[fix[1],j-1,index3]))
        }
      }
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=y[i,j,index3-1], max=min(y[i, j+1, index3], y[i+1, j, index3]))
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i,j-1,index3], y[i,j,index3-1]), max=y[i+1,j,index3])
            }
          }
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=max(y[i-1,j,index3], y[i,j,index3-1]), max=y[i, j+1, index3])
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3], y[i,j,index3-1]))
            }
          }
        }
      }
    }
  }
  if(fix[3]>1) {
    if(fix[1]>1) {
      for (i in rev(1:(fix[1]-1))) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,fix[3]] <- runif(1, max=min(y[i, j+1, fix[3]], y[i+1, j, fix[3]]))
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,fix[3]] <- runif(1, min=y[i,j-1,fix[3]], max=y[i+1,j,fix[3]])
          }
        }
      }
    }
    if(fix[1]<dim1) {
      for (i in (fix[1]+1):dim1) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,fix[3]] <- runif(1, min=y[i-1,j,fix[3]], max=y[i, j+1, fix[3]])
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,fix[3]] <- runif(1, min=max(y[i, j-1, fix[3]], y[i-1, j, fix[3]]))
          }
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          y[i,fix[2],index3] <- runif(1, max=min(y[i+1,fix[2],index3], y[i,fix[2],index3+1]))
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          y[i,fix[2],index3] <- runif(1, min=y[i-1,fix[2],index3], max=y[i,fix[2],index3+1])
        }
      }
      if(fix[2]>1) {
        for (j in rev(1:(fix[2]-1))) {
          y[fix[1],j,index3] <- runif(1, max=min(y[fix[1],j+1,index3], y[fix[1],j,index3+1]))
        }
      }
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[fix[1],j,index3] <- runif(1, min=y[fix[1],j-1,index3], max=y[fix[1],j,index3+1])
        }
      }
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, max=min(y[i, j+1, index3], y[i+1, j, index3], y[i,j,index3+1]))
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=y[i,j-1,index3], max=min(y[i+1,j,index3], y[i,j,index3+1]))
            }
          }
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=y[i-1,j,index3], max=min(y[i, j+1, index3], y[i,j,index3+1]))
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3]), max=y[i,j,index3+1])
            }
          }
        }
      }
    }
  }
  
  
  y <- as.vector(y)
  #y <- y[c(1, 4, 5, 8, 9, 11, 12, 17, 20, 21, 23, 24)]
  y
}


#  Motivating example 1 (used in the paper) with 12 combinations
# Listing all 148 orderings
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
