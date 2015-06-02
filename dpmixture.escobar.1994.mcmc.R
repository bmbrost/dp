dpmixture.escobar.1994.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  a0 <- start$a0
  z <- start$z
  n <- length(y)
  
  z.save <- matrix(0,n,n.mcmc)

  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample z (cluster assignments for observations)
    ###
    
    # Following Escobar (1994) Eq. 3
    
    for(i in 1:n){
      y.tmp <- y[i]
      theta.tmp <- setdiff(z,z[i])
      denom <- a0+sum(dnorm(y.tmp,theta.tmp,1))  # denominator of eq. 3 
      p <- c(dnorm(y.tmp,theta.tmp,1)/denom,a0/denom)  # probability
      z[i] <- sample(c(theta.tmp,rnorm(1,y.tmp,1)),1,prob=p)  # 'new' cluster
    }
    
    
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