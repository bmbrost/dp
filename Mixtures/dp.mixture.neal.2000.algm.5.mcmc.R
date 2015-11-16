dpmixture.neal.2000.algm.5.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # browser()
  a0 <- start$a0
  z <- start$z
  n <- length(y)  # number of observations
  P0 <- range(y)+c(-10,10)  # redefine P0 to contain support of y plus extra
  F.P0 <- 1/(max(P0)-min(P0))

  z.save <- matrix(0,n,n.mcmc)
  keep <- list(z=0)
  
  
  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample z (cluster assignments for observations)
    ###
    
    # Following Neal (2000), algorithm 5

# browser()    
      
    for(i in 1:n){
      y.tmp <- y[i]
      delta.j <- table(z[-i])  # number of observations per cluster
      phi.c <- sort(unique(z[-i]))  # unique clusters; same ordering as delta.j
      p.old <- delta.j/(n-1+a0)  # prob of existing phi
      p.new <- a0/(n-1+a0)  # prob of new phi
#       z.star <- sample(c(phi.c,rnorm(1,y.tmp,1)),1,prob=c(p.old,p.new))  # sample phi  
      z.star <- sample(c(phi.c,rnorm(1,z[i],tune$z)),1,prob=c(p.old,p.new))  # sample phi  
      mh.star.z <- dnorm(y.tmp,z.star,1,log=TRUE)
      mh.0.z <- dnorm(y.tmp,z[i],1,log=TRUE)
      if(exp(mh.star.z-mh.0.z)>runif(1)){
        z[i] <- z.star
        keep$z <- keep$z+1
      } 
    }    
    
    delta.j <- table(z)
    phi.c <- sort(unique(z))
    for(i in length(phi.c)){
      idx <- which(z==phi.c[i])
      y.tmp <- y[idx]
      A <- delta.j[i]/1
      b <- sum(y.tmp)
      z[idx] <- rnorm(1,1/A*b,1/A)
    }

#     for(i in 1:length(phi.c)){
#       idx <- which(z==phi.c[i])
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
  
  keep$z <- keep$z/(n.mcmc*n)
  cat(paste("\nz acceptance rate:",keep$z)) 
  list(z=z.save,n.mcmc=n.mcmc)
 
}