dp.mix.neal.2000.algm.8.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # browser()
  a0 <- start$a0
  z <- start$z
  m <- priors$m
  n <- length(y)  # number of observations
  P0.tmp <- P0
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
    
    # Following Neal (2000), Algorithm 8

# browser()  
    
    phi.c <- sort(unique(z))
    for(i in 1:n){
      delta.j <- table(z[-i])
      phi.c.tmp <- sort(unique(z[-i]))
      singleton <- length(phi.c)!=length(phi.c.tmp)
# see duplicated() in data.table package    
      if(singleton){
        h <- c(z[i],runif(m-1,min(P0),max(P0)))
        p.old <- (delta.j/(n-1+a0))*dnorm(y[i],phi.c.tmp,1)
        p.new <- ((a0/m)/(n-1+a0))*dnorm(y[i],h,1)
        z[i] <- sample(c(phi.c.tmp,h),1,prob=c(p.old,p.new))  
      } 
      else{
        h <- runif(m,min(P0),max(P0))    
        p.old <- (delta.j/(n-1+a0))*dnorm(y[i],phi.c.tmp,1)
        p.new <- ((a0/m)/(n-1+a0))*dnorm(y[i],h,1)
        z[i] <- sample(c(phi.c.tmp,h),1,prob=c(p.old,p.new))  
      }
    }
    
      
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