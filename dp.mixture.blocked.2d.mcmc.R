dpmixture.blocked.2d.mcmc <- function(y,P0,priors,tune,start,n.mcmc,n.cores=NULL){
  
  t.start <- Sys.time()
  
  ###
  ### Libraries and Subroutines
  ###
  
  library(MCMCpack)  # for Dirichlet distribution functions
  library(data.table)  # for tabulating and summing
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  library(doParallel) #For parallel processing
  library(foreach) #For parallel processing  
  
  get.mu.0 <- function(x,dt,tab,sigma){
    y.tmp <- dt[z1==tab[x,z1]&z2==tab[x,z2],.(y1,y2)]
    y.tmp <- as.matrix(y.tmp)
    A <- solve(sigma^2*diag(2))*tab[x,N]
	A.inv <- solve(A)
    b <- colSums(y.tmp%*%solve(sigma^2*diag(2)))
    t(A.inv%*%b)
    rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))
  }
  
  
  ###
  ###  Create cluster for parallel processing
  ###
  
  #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  #   mcoptions <- list(preschedule=TRUE)
  #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  ###
  ### Starting values and priors
  ###
  
  a0 <- start$a0
  z <- start$z
  sigma <- start$sigma
  pie <- start$pie
  H <- priors$H
  
  
  ###
  ###  Setup Variables 
  ###
  
  #   browser()  
  n <- nrow(y)  # number of observations
  dt <- as.data.table(cbind(y1=y[,1],y2=y[,2],z1=z[,1],z2=z[,2]))
  tab <- dt[,.N,by=.(z1,z2)]  
  setkey(tab ,N)
  n.cluster <- tab[,.N]  
  ord <- n.cluster:1
  
  z.save <- array(0,dim=c(n,2,n.mcmc))
  a0.save <- numeric(n.mcmc)
  sigma.save <- numeric(n.mcmc)  
  keep <- list(sigma=0)

    
  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    if(k%%1000==0) cat(k,"");flush.console()
    
    # Following Ishwaran and James (2001); Also Gelman et al. (2014), Section 23.3
    
    ###
    ### Sample theta (cluster parameter)
    ###
    
    # Sampler currently disregards support P0
# browser()
    mu.0 <- t(sapply(1:n.cluster,function(x) get.mu.0(x,dt,tab,sigma)))
    # mu.0 <- cbind(rnorm(n.cluster,mu.0[ord,1],sigma/sqrt(tab[ord,N])),
      # rnorm(n.cluster,mu.0[ord,2],sigma/sqrt(tab[ord,N])))  
    n.cluster.new <- H-n.cluster
    mu.0 <- rbind(mu.0[ord,],cbind(runif(n.cluster.new,P0[1,1],P0[2,1]),
      runif(n.cluster.new,P0[1,2],P0[3,2])))
       
  
    ###
    ### Sample h.t (cluster assignments)
    ###
    
    # browser()
    h.t <- sapply(1:n,function(x) sample(1:H,1,
      prob=pie*dnorm(y[x,1],mu.0[,1],sigma)*dnorm(y[x,2],mu.0[,2],sigma)))
    dt[,z1:=mu.0[h.t,1]]
    dt[,z2:=mu.0[h.t,2]]    


    ###
    ### Sample sigma (observation error)
    ###

    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      mh.star.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma.star,log=TRUE)+
        dnorm(dt[,y2],dt[,z2],sigma.star,log=TRUE))
      mh.0.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma,log=TRUE)+
        dnorm(dt[,y2],dt[,z2],sigma,log=TRUE))
      if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        sigma <- sigma.star
        keep$sigma <- keep$sigma+1
      } 
    }


    ###
    ### Sample p (stick-breaking process)
    ###

    tab <- dt[,.N,by=.(z1,z2)]  
    setkey(tab ,N)
    n.cluster <- tab[,.N]  
    ord <- n.cluster:1
    tab.tmp <- c(tab[ord,N],rep(0,H-n.cluster-1))

    # Update stick-breaking weights
    v <- c(rbeta(H-1,1+tab.tmp,a0+n-cumsum(tab.tmp)),1)
    pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

    
    ###
    ### Sample a0 (concentration parameter); See Gelman section 23.3
    ###
    
    a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
    

    ###
    ###  Save samples 
    ###

    z.save[,1,k] <- dt[,z1]
    z.save[,2,k] <- dt[,z2]
    a0.save[k] <- a0    
    sigma.save[k] <- sigma
  }
  
  ###
  ### Write output
  ###
  
  keep$sigma <- keep$sigma/n.mcmc
  cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
  list(z=z.save,a0=a0.save,sigma=sigma.save,keep=keep,n.mcmc=n.mcmc)
  
}