dirichlet.process.prior.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
  a0 <- start$a0
  F <- 1/(max(P0)-min(P0))
  n <- length(y)
  
  y.tab <- table(y)
  k <- as.numeric(names(y.tab))
  q.k <- length(k)
  
#   browser() 
  p <- rep((a0*F)/(a0+n-1),q.k)+y.tab
      
  P.save <- matrix(0,q.k,n.mcmc)
  a0.save <- numeric(n.mcmc)
  
  keep <- list(a0=0)

  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample a0
    ###
    
    a0.star <- rnorm(1,a0,tune$a0)
    if (a0.star >= 0) {
      p.star <- rep((a0.star*F)/(a0.star+n-1),q.k)+y.tab
      mh.star.a0 <- log(ddirichlet(y.tab/100,p.star))+dgamma(a0.star,priors$a,priors$b,log=TRUE)
      mh.0.a0 <- log(ddirichlet(y.tab/100,p))+dgamma(a0,priors$a,priors$b,log=TRUE)
      if (exp(mh.star.a0-mh.0.a0)>runif(1)) {
        a0 <- a0.star
        p <- p.star
        keep$a0 <- keep$a0+1
      }  
    }
    
    
    ###
    ### Sample P
    ###
    
    P <- rdirichlet(1,p)
    
      
    ###
    ###  Save samples 
    ###

    P.save[, k] <- P
    a0.save[k] <- a0
  }
  
  ###
  ### Write output
  ###
  
  keep$a0 <- keep$a0/n.mcmc
  cat(paste("\na0 acceptance rate:",round(keep$a0,2)))  
  list(P=P.save,a0=a0.save,keep=keep,n.mcmc=n.mcmc)
  
}