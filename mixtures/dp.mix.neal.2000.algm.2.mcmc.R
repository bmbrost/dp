dp.mix.neal.2000.algm.2.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  
  # browser()
  a0 <- start$a0
  z <- start$z
  n <- length(y)  # number of observations
  P0 <- range(y)+c(-10,10)  # redefine P0 to contain support of y plus extra
  F.P0 <- 1/(max(P0)-min(P0))

  z.save <- matrix(0,n,n.mcmc)
  
  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample z (cluster assignments for observations)
    ###
    
    # Following Neal (2000), algorithm 2

    delta.j <- table(z)  # number of observations per cluster
    phi.c <- sort(unique(z))  # unique clusters; same ordering as delta.j
    
    for(i in 1:length(phi.c)){
      idx <- which(z==phi.c[i])
      y.tmp <- y[idx]
      if(delta.j[i]==1){  # update singleton clusters
        p.old <- delta.j[-i]/(n-1+a0)*dnorm(y.tmp,phi.c[-i],1) # prob of existing phi
        p.new <- a0/(n-1+a0) * (a0*F.P0) # prob of new phi
        z[idx] <- sample(c(phi.c[-1],rnorm(1,y.tmp,1)),1,prob=c(p.old,p.new))  # sample phi  
      }
      if(delta.j[i]>1){  # group update for all other clusters
        A <- delta.j[i]/1
        b <- sum(y.tmp)
        z[idx] <- rnorm(1,1/A*b,1/A)
      }
    }  

    ### Method for keeping track of cluster idx
    
    #     cluster.idx <- dense_rank(z)  # rank of clusters smallest to largest
    #     n.cluster <- max(cluster.idx)
    #     for(i in 1:n.cluster){
    #       idx <- which(cluster.idx==i)
    #       y.tmp <- y[idx]
    #       if(delta.j[i]==1){  # update singleton clusters
    #         p.old <- delta.j[-i]/(n-1+a0)*dnorm(y.tmp,phi.c[-i],1) # prob of existing phi
    #         p.new <- a0/(n-1+a0) * (a0*F.P0) # prob of new phi
    #         z[idx] <- sample(c(phi.c[-1],rnorm(1,y.tmp,1)),1,prob=c(p.old,p.new))  # sample phi  
    #       }
    #       if(delta.j[i]>1){  # group update for all other clusters
    #         A <- delta.j[i]/1
    #         b <- sum(y.tmp)
    #         z[idx] <- rnorm(1,1/A*b,1/A)
    #       }
    #     }  
    
    ###
    ###  Save samples 
    ###
    
    z.save[,k] <- z
    
  }
  
  ###
  ### Write output
  ###
  
  list(z=z.save,n.mcmc=n.mcmc)
 
}