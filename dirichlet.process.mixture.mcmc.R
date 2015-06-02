dirichlet.process.mixture.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
# browser()
  a0 <- start$a0
  z <- start$z
  N <- priors$N  # maximum number of clusters
  q.P0 <- length(P0)
  F <- 1/q.P0
  n <- length(y)
  
#   browser()
#   tab <- table(y)
#   idx <- as.numeric(names(tmp))
#   y.tab <- numeric(q.P0)
#   y.tab[idx] <- tmp
#   p <- rep(a0*F,q.P0) + y.tab
  
  
  P.save <- matrix(0,q.P0,n.mcmc)
  a0.save <- numeric(n.mcmc)
  z.save <- matrix(0,n,n.mcmc)

  keep <- list(a0=0)

  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample a0
    ###
    
#     a0.star <- rnorm(1,a0,tune$a0)
#     if (a0.star >= 0) {
#       p.star <- rep(a0.star*F,q.P0) + y.tab
#       mh.star.a0 <- log(ddirichlet(y.tab/100,p.star))+dgamma(a0.star,priors$a,priors$b,log=TRUE)
#       mh.0.a0 <- log(ddirichlet(y.tab/100,p))+dgamma(a0,priors$a,priors$b,log=TRUE)
#       if (exp(mh.star.a0-mh.0.a0)>runif(1)) {
#         a0 <- a0.star
#         p <- p.star
#         keep$a0 <- keep$a0+1
#       }  
#     }
    
    ###
    ### Sample z (cluster assignments for observations)
    ###
    
    for(i in 1:n){
      y.tmp <- y[i]
      theta.tmp <- setdiff(z,z[i])
      denom <- a0+sum(dnorm(y.tmp,theta.tmp,1))
      p <- c(dnorm(y[idx],theta.tmp,1)/denom,a0/denom)
      z[i] <- sample(c(theta.tmp,rnorm(1,y.tmp,1)),1,prob=p)  
    }
    
    
    
    
    ###
    ### Sample v
    ###

#     v <- c(rbeta(N-1,1,a0),1)
#     pie <- v*c(1,cumprod((1-v[-N])))
#     
    
    
    ###
    ### Sample P
    ###
    
#     P <- rdirichlet(1,p)
    
      
    ###
    ###  Save samples 
    ###

    
#     P.save[, k] <- P
      z.save[,k] <- z
      a0.save[k] <- a0
  }
  
  ###
  ### Write output
  ###
  
  keep$a0 <- keep$a0/n.mcmc
  list(P=P.save,a0=a0.save,z=z.save,keep=keep,n.mcmc=n.mcmc)
  
}