dpmixture.blocked.2d.mcmc <- function(y,P0,priors,tune,start,n.mcmc,n.cores=NULL){
  
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
  
  
#   browser()
  a0 <- start$a0
  z <- start$z
  sigma <- start$sigma
  pie <- start$pie
  H <- priors$H
  n <- nrow(y)  # number of observations
    
  z.save <- array(0,dim=c(n,2,n.mcmc))
  a0.save <- numeric(n.mcmc)
  sigma.save <- numeric(n.mcmc)  
  keep <- list(sigma=0)
  
  uz <- paste(z[,1],z[,2])
  A <- solve(sigma^2*diag(2))
  tapply(y,uz,function(x),)

  A <- solve(sigma^2*diag(2))
  b <- y%*%A
  dt <- as.data.table(cbind(y1=y[,1],y2=y[,2],z1=z[,1],z2=z[,2],b%*%solve(A)))
  tab <- dt[,.(.N,y1.sum=sum(V5),y2.sum=sum(V6)),by=.(z1,z2)]
  tab[, c("y1.sum","y2.sum"):=list(y1.sum/N,y2.sum/N)]
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
    
    # Sampler currently disregards support P0
   
    theta.star <- cbind(rnorm(n.theta,rev(tab[,y1.sum]),1/sqrt(A[1,1])),
      rnorm(n.theta,rev(tab[,y2.sum]),1/sqrt(A[2,2])))  
    n.theta.new <- H-n.theta
    theta <- rbind(theta.star,cbind(runif(n.theta.new,P0[1,1],P0[2,1]),
      runif(n.theta.new,P0[1,2],P0[3,2])))
    
#     mu.0.star <- rnorm(2,mu.0,tune$mu.0)
#     if(mu.0.star[1]>S.tilde[1,1]&mu.0.star[1]<S.tilde[2,1]& #Reject proposals for mu.0 not in S.tilde
#          mu.0.star[2]>S.tilde[1,2]&mu.0.star[2]<S.tilde[3,2]){
#       idx <- which(z==0) #Update using z==0 only
#       mh.star.mu.0 <- sum(dtnorm(mu[idx,1],mu.0.star[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
#          dtnorm(mu[idx,2],mu.0.star[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
#       mh.0.mu.0 <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
#          dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
#       if(exp(mh.star.mu.0-mh.0.mu.0)>runif(1)){
#         mu.0 <- mu.0.star
#         keep$mu.0 <- keep$mu.0+1
#       } 
#     }    
    
    # browser()
#     A <- z.tab/sigma
#     y.sum <- rev(tab[,y.sum])
#     b <- y.sum/sigma
#     #     A <- table(z)/sigma
#     #     b <- tapply(y,z,sum)/sigma
#     theta <- c(rnorm(n.theta,b/A,1/A),runif(H-n.theta,min(P0),max(P0)))
    
    ###
    ### Sample z (cluster assignments)
    ###
    
#         browser()

    idx <- sapply(1:n,function(x) sample(1:H,1,
      prob=pie*dnorm(dt[x,y1],theta[,1],sigma)*dnorm(dt[x,y2],theta[,2],sigma)))
    dt[,z1:=theta[idx,1]]
    dt[,z2:=theta[idx,2]]

    # dt[,.(z1=theta[idx,1],z2=theta[idx,2])]        
    # DT[, c("V1","V2") := list (round(exp(V1),2), LETTERS [4:6])]

# plot(y)
# points(theta,col=2)
# segments(y[,1],y[,2],theta[idx,1],theta[idx,2])    
# 
# points(z,col=3)
# points(tab[,.(y1.sum,y2.sum)],col=4)

#     idx <- sapply(dt[,y],function(x) sample(1:H,1,prob=pie*dnorm(x,theta,sigma)))
#     dt[,z := theta[idx]] 
    
    #     for(i in 1:n){
    #       p.tmp <- pie*dnorm(y[i],theta,sigma)
    #       idx <- sample(1:N,1,prob=p.tmp)  # normalizing probability is unneseccary
    #       z[i] <- theta[idx]
    #     }
    
    #     idx <- sapply(y,function(x) sample(1:H,1,prob=pie*dnorm(x,theta,sigma)))
    #     z <- theta[idx]
    

    ###
    ### Sample sigma (observation error)
    ###

    # browser()
    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      mh.star.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma.star,log=TRUE)+
        dnorm(dt[,y2],dt[,z2],sigma.star,log=TRUE))
      mh.0.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma,log=TRUE)+
        dnorm(dt[,y2],dt[,z2],sigma,log=TRUE))
      if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        sigma <- sigma.star
        keep$sigma <- keep$sigma+1
        A <- solve(sigma^2*diag(2))
        b <- y%*%(A)
        tmp <- b%*%solve(A)
        dt[,V5:=tmp[,1]]
        dt[,V6:=tmp[,2]]
#         dt[,c("V5","V6"):=list(tmp[,1],tmp[,2])]
      
      } 
    }
    


    ###
    ### Sample p (stick-breaking process)
    ###
      
    tab <- dt[,.(.N,y1.sum=sum(V5),y2.sum=sum(V6)),by=.(z1,z2)]
    tab[, c("y1.sum","y2.sum"):=list(y1.sum/N,y2.sum/N)]
    setkey(tab,N)
    n.theta <- tab[,.N]
    z.tab <- rev(tab[,N])

#     tab <- dt[,.(.N,y.sum=sum(y)),by="z"]
#     setkey(tab,N)
#     n.theta <- tab[,.N]
    
    # Number of members per cluster in decreasing order
#     z.tab <- rev(tab[,N])
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
    ###  Save samples 
    ###

#     z.save[,,k] <- dt[,.(z1,z2)]
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
  list(z=z.save,a0=a0.save,sigma=sigma.save,keep=keep,n.mcmc=n.mcmc)
  
}