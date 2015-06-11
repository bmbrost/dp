dpmixture.blocked.mcmc <- function(y,P0,priors,tune,start,n.mcmc,n.cores=NULL){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  library(data.table)  # for tabulating and summing
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  library(doParallel) #For parallel processing
  library(foreach) #For parallel processing	
  
  ###
  ###  Create cluster for parallel processing
  ###
  
#   if(is.null(n.cores)) n.cores <- detectCores() - 1
#   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
#   mcoptions <- list(preschedule=TRUE)
#   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  # browser()
  a0 <- start$a0
  z <- start$z
  sigma <- 1
  pie <- start$pie
  H <- priors$H
  n <- length(y)  # number of observations
  P0 <- range(y)+c(-10,10)  # redefine P0 to contain support of y plus extra
#   F.P0 <- 1/(max(P0)-min(P0))
  
  z.save <- matrix(0,n,n.mcmc)
  a0.save <- numeric(n.mcmc)
  sigma.save <- numeric(n.mcmc)  
  keep <- list(sigma=0)

  dt <- as.data.table(cbind(y,z))
  tab <- dt[,.(.N,y.sum=sum(y)),by="z"]
  setkey(tab,N)
  n.theta <- tab[,.N]
  z.tab <- rev(tab[,N])


  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {

    # Following Ishwaran and James (2001); Also Gelman et al. (2014), Section 23.3
    
    ###
    ### Sample theta (cluster parameter)
    ###

# browser()
    A <- z.tab/sigma
    y.sum <- rev(tab[,y.sum])
    b <- y.sum/sigma
#     A <- table(z)/sigma
#     b <- tapply(y,z,sum)/sigma
    theta <- c(rnorm(n.theta,b/A,1/A),runif(H-n.theta,min(P0),max(P0)))

    ###
    ### Sample z (cluster assignments)
    ###
    
#     browser()
    idx <- sapply(dt[,y],function(x) sample(1:H,1,prob=pie*dnorm(x,theta,sigma)))
    dt[,z := theta[idx]] 

#     for(i in 1:n){
#       p.tmp <- pie*dnorm(y[i],theta,sigma)
#       idx <- sample(1:N,1,prob=p.tmp)  # normalizing probability is unneseccary
#       z[i] <- theta[idx]
#     }

#     idx <- sapply(y,function(x) sample(1:H,1,prob=pie*dnorm(x,theta,sigma)))
#     z <- theta[idx]

    
    ###
    ### Sample p (stick-breaking process)
    ###
    
    tab <- dt[,.(.N,y.sum=sum(y)),by="z"]
    setkey(tab,N)
    n.theta <- tab[,.N]
    
    # Number of members per cluster in decreasing order
    z.tab <- rev(tab[,N])
# z.tab <- tab[,N]
    z.tab.tmp <- c(z.tab,rep(0,H-n.theta-1))

    # Update stick-breaking weights
    v <- c(rbeta(H-1,1+z.tab.tmp,a0+n-cumsum(z.tab.tmp)),1)
    pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

#     z.tab <- table(z)
#     n.theta <- length(z.tab)
#     z.tab <- c(z.tab,rep(0,N-n.theta))
#     v <- c(rbeta(N-1,1+z.tab[-N],a0+n-cumsum(z.tab[-N])),1)
#     pie <- v*c(1,cumprod((1-v[-N])))  
    
    
    ###
    ### Sample a0 (concentration parameter); See Gelman section 23.3
    ###

    a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
    
    
    ###
    ### Sample sigma (observation error)
    ###

# browser()
    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      mh.star.sigma <- sum(dnorm(dt[,y],dt[,z],sigma.star,log=TRUE))
      mh.0.sigma <- sum(dnorm(dt[,y],dt[,z],sigma,log=TRUE))
      if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        sigma <- sigma.star
        keep$sigma <- keep$sigma+1
      } 
    }


    ###
    ###  Save samples 
    ###
    
    z.save[,k] <- dt[,z]
    a0.save[k] <- a0    
    sigma.save[k] <- sigma
  }
  
  ###
  ### Write output
  ###
  
  keep$sigma <- keep$sigma/n.mcmc
  cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  list(z=z.save,a0=a0.save,sigma=sigma.save,keep=keep,n.mcmc=n.mcmc)
 
}