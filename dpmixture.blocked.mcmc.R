dpmixture.blocked.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # browser()
  a0 <- start$a0
  z <- start$z
  sigma <- 1
  N <- priors$N
  n <- length(y)  # number of observations
  P0 <- range(y)+c(-10,10)  # redefine P0 to contain support of y plus extra
  F.P0 <- 1/(max(P0)-min(P0))

  
  z.save <- matrix(0,n,n.mcmc)
  a0.save <- numeric(n.mcmc)
  keep <- list(z=0)
  
  
  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {

    # Following Ishwaran and James (2001); Also Gelman et al. (2014), Section 23.3
# browser()  
    
    ###
    ### Sample theta (cluster parameter)
    ###
    
    A <- table(z)/sigma
    b <- tapply(y,z,sum)/sigma
    n.theta <- length(A)
    theta <- c(rnorm(n.theta,b/A,1/A),runif(N-n.theta,min(P0),max(P0)))

    ###
    ### Sample z (cluster assignments)
    ###
    
    for(i in 1:n){
      p.tmp <- pie*dnorm(y[i],theta,sigma)
      idx <- sample(1:N,1,prob=p.tmp/sum(p.tmp))
      z[i] <- theta[idx]
    }
    
    
    ###
    ### Sample p (stick-breaking process)
    ###

    z.tab <- table(z)
    n.theta <- length(z.tab)
    z.tab <- c(z.tab,rep(0,N-n.theta))
    v <- c(rbeta(N-1,1+z.tab[-N],a0+n-cumsum(z.tab[-N])),1)
    pie <- v*c(1,cumprod((1-v[-N])))  

    ###
    ### Sample a0 (concentration parameter); See Gelman section 23.3
    ###

    a0 <- rgamma(1,priors$r+N-1,priors$q-sum(log(1-v[-N])))  
    
    
    ###
    ###  Save samples 
    ###
    
    z.save[,k] <- z
    a0.save[k] <- a0    
  }
  
  ###
  ### Write output
  ###
  
  keep$z <- keep$z/(n.mcmc*n)
  cat(paste("\nz acceptance rate:",keep$z)) 
  list(z=z.save,a0=a0.save,n.mcmc=n.mcmc)
 
}