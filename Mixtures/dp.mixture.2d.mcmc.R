# Dirichlety process mixture model for 2-dimensional data

dpmixture.2d.mcmc <- function(s,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
  
	t.start <- Sys.time()
  
	###
	### Libraries and Subroutines
	###
  
	# library(MCMCpack)  # for Dirichlet distribution functions
	# library(data.table)  # for tabulating and summing
	# library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# library(doParallel) #For parallel processing
	# library(foreach) #For parallel processing  
  
  	get.mu.0 <- function(x,s,h,sigma){
		# browser()
		idx <- which(h==x)
		n <- length(idx)
		Sigma.inv <- solve(sigma^2*diag(2))
		b <- colSums(s[idx,]%*%Sigma.inv)
		A.inv <- solve(n*Sigma.inv)
		rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
	}

	###
	### Starting values and priors
	###
  
	theta <- start$theta
	mu.0 <- start$mu.0
	sigma <- start$sigma
	pie <- start$pie
	H <- priors$H
 
  
	###
	###  Setup Variables 
	###
  
    # browser()  
	T <- nrow(s)  # number of observations

	# Setup up Dirichlet process mixture variables
	mu.0 <- unique(mu.0) #  unique clusters
	m <- nrow(mu.0)  # number of clusters
	ht <- c(1:m)[match(start$mu.0[,1],mu.0[,1])]  # idx assigning each obs to record in mu.0
	tab <- table(ht)  # tabulate cluster membership
	ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
	samp <- as.numeric(names(tab))[ord]

	mu.0 <- rbind(mu.0,matrix(0,H-m,2))  # Add unoccupied 'dummy' clusters


	###
	### Create receptacles for output
	###

	mu.0.save <- array(0,dim=c(T,2,n.mcmc))
	theta.save <- numeric(n.mcmc)
	sigma.save <- numeric(n.mcmc)  
	m.save <- numeric(n.mcmc)  
	keep <- list(sigma=0)

    
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) cat(k,"");flush.console()
# print(k)    
		###
		### Dirichlet process parameters
		###	

		# Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# and Gelman et al. (2014), Section 23.3
	
		# Sample 'occupied' mu.0 (true location of clusters with non-zero membership)		 
	    # Sampler currently disregards support P0
# browser()
		mu.0[samp,] <- t(sapply(samp,function(x) get.mu.0(x,s,ht,sigma)))
		# mu.0 <- t(sapply(samp,function(x) get.mu.0(x,s,ht,sigma)))
		
		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior 
	    mu.0[-samp,] <- cbind(runif(H-m,min(S.tilde[,1]),max(S.tilde[,1])),
			runif(H-m,min(S.tilde[,2]),max(S.tilde[,2])))

	    # mu.0 <- rbind(mu.0,cbind(runif(H-m,min(S.tilde[,1]),max(S.tilde[,1])),
			# runif(H-m,min(S.tilde[,2]),max(S.tilde[,2]))))

	    # Note: sampling order matters for the remaining DP updates. Cluster parameters 
	    # must be sampled in same order as pie, i.e., sorted by decreasing membership
    
	    # Sample cluster assignment indicator
		samp <- c(samp,c(1:H)[-samp])
		ht <- sapply(1:T,function(x) sample(samp,1,
			prob=exp(log(pie)+dnorm(s[x,1],mu.0[samp,1],sigma,log=TRUE)+
			dnorm(s[x,2],mu.0[samp,2],sigma,log=TRUE))))
		mu.0.tmp <- mu.0[ht,]

		# Tabulate cluster membership with base functions
		tab <- table(ht) #  tabulate cluster membership
		m <- length(tab) #  number of clusters
		ord <- order(tab,decreasing=TRUE) #  idx of tab ordered by membership
		samp <- as.numeric(names(tab))[ord]
		
		# Stick-breaking process
		pad <- ifelse(m<H,H-m-1,0)
		tab.tmp <- c(tab[ord],rep(0,pad))  # membership in decreasing order
		eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities
   
		# Sample theta (concentration parameter)
	    # See Gelman section 23.3, Ishwaran and Zarepour (2000)
       	
       	# theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  

		tmp <- rbeta(1,theta+1,T)
		c <- priors$r
		d <- priors$q
		p.tmp <- (c+m-1)/(c+m-1+T*(d-log(tmp)))
		p.tmp <- rbinom(1,1,p.tmp)
		theta <- ifelse(p.tmp==1,rgamma(1,c+m,d-log(tmp)),rgamma(1,c+m-1,d-log(tmp)))
		
		###
		### Sample sigma (observation error)
	    ###

    	sigma.star <- rnorm(1,sigma,tune$sigma)
	    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
			mh.star.sigma <- sum(dnorm(s[,1],mu.0.tmp[,1],sigma.star,log=TRUE)+
	        	dnorm(s[,2],mu.0.tmp[,2],sigma.star,log=TRUE))
    		mh.0.sigma <- sum(dnorm(s[,1],mu.0.tmp[,1],sigma,log=TRUE)+
		        dnorm(s[,2],mu.0.tmp[,2],sigma,log=TRUE))
		    if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        		sigma <- sigma.star
		        keep$sigma <- keep$sigma+1
    		} 
	    }
 

	    ###
    	###  Save samples 
	    ###

    	mu.0.save[,,k] <- mu.0.tmp
		m.save[k] <- m
    	theta.save[k] <- theta    
	    sigma.save[k] <- sigma
	}
  
	###
	### Write output
	###
  
	keep$sigma <- keep$sigma/n.mcmc
	cat(paste("\nsigma acceptance rate:",keep$sigma)) 
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	list(mu.0=mu.0.save,theta=theta.save,sigma=sigma.save,m=m.save,keep=keep,n.mcmc=n.mcmc)
 }


####
#### This version uses data.table for tabluation. Reference simulation code that's 
#### commented out at the bottom of dp.mixture.2d.sim.R
####

# dpmixture.blocked.2d.mcmc <- function(y,P0,priors,tune,start,n.mcmc,n.cores=NULL){
  
  # t.start <- Sys.time()
  
  # ###
  # ### Libraries and Subroutines
  # ###
  
  # library(MCMCpack)  # for Dirichlet distribution functions
  # library(data.table)  # for tabulating and summing
  # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # library(doParallel) #For parallel processing
  # library(foreach) #For parallel processing  
  
  # get.mu.0 <- function(x,dt,tab,sigma){
    # y.tmp <- dt[z1==tab[x,z1]&z2==tab[x,z2],.(y1,y2)]
    # y.tmp <- as.matrix(y.tmp)
    # A <- solve(sigma^2*diag(2))*tab[x,N]
	# A.inv <- solve(A)
    # b <- colSums(y.tmp%*%solve(sigma^2*diag(2)))
    # t(A.inv%*%b)
    # rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))
  # }
  
  
  # ###
  # ###  Create cluster for parallel processing
  # ###
  
  # #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  # #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  # #   mcoptions <- list(preschedule=TRUE)
  # #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  # ###
  # ### Starting values and priors
  # ###
  
  # a0 <- start$a0
  # z <- start$z
  # sigma <- start$sigma
  # pie <- start$pie
  # H <- priors$H
  
  
  # ###
  # ###  Setup Variables 
  # ###
  
  # #   browser()  
  # n <- nrow(y)  # number of observations
  # dt <- as.data.table(cbind(y1=y[,1],y2=y[,2],z1=z[,1],z2=z[,2]))
  # tab <- dt[,.N,by=.(z1,z2)]  
  # setkey(tab ,N)
  # n.cluster <- tab[,.N]  
  # ord <- n.cluster:1
  
  # z.save <- array(0,dim=c(n,2,n.mcmc))
  # a0.save <- numeric(n.mcmc)
  # sigma.save <- numeric(n.mcmc)  
  # keep <- list(sigma=0)

    
  # ###
  # ### Begin MCMC loop
  # ###
  
  # for (k in 1:n.mcmc) {
    # if(k%%1000==0) cat(k,"");flush.console()
    
    # # Following Ishwaran and James (2001); Also Gelman et al. (2014), Section 23.3
    
    # ###
    # ### Sample theta (cluster parameter)
    # ###
    
    # # Sampler currently disregards support P0
# # browser()
    # mu.0 <- t(sapply(1:n.cluster,function(x) get.mu.0(x,dt,tab,sigma)))
    # # mu.0 <- cbind(rnorm(n.cluster,mu.0[ord,1],sigma/sqrt(tab[ord,N])),
      # # rnorm(n.cluster,mu.0[ord,2],sigma/sqrt(tab[ord,N])))  
    # n.cluster.new <- H-n.cluster
    # mu.0 <- rbind(mu.0[ord,],cbind(runif(n.cluster.new,P0[1,1],P0[2,1]),
      # runif(n.cluster.new,P0[1,2],P0[3,2])))
       
  
    # ###
    # ### Sample h.t (cluster assignments)
    # ###
    
    # # browser()
    # h.t <- sapply(1:n,function(x) sample(1:H,1,
      # prob=pie*dnorm(y[x,1],mu.0[,1],sigma)*dnorm(y[x,2],mu.0[,2],sigma)))
    # dt[,z1:=mu.0[h.t,1]]
    # dt[,z2:=mu.0[h.t,2]]    


    # ###
    # ### Sample sigma (observation error)
    # ###

    # sigma.star <- rnorm(1,sigma,tune$sigma)
    # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      # mh.star.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma.star,log=TRUE)+
        # dnorm(dt[,y2],dt[,z2],sigma.star,log=TRUE))
      # mh.0.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma,log=TRUE)+
        # dnorm(dt[,y2],dt[,z2],sigma,log=TRUE))
      # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        # sigma <- sigma.star
        # keep$sigma <- keep$sigma+1
      # } 
    # }


    # ###
    # ### Sample p (stick-breaking process)
    # ###

    # tab <- dt[,.N,by=.(z1,z2)]  
    # setkey(tab ,N)
    # n.cluster <- tab[,.N]  
    # ord <- n.cluster:1
    # tab.tmp <- c(tab[ord,N],rep(0,H-n.cluster-1))

    # # Update stick-breaking weights
    # v <- c(rbeta(H-1,1+tab.tmp,a0+n-cumsum(tab.tmp)),1)
    # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

    
    # ###
    # ### Sample a0 (concentration parameter); See Gelman section 23.3
    # ###
    
    # a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
    

    # ###
    # ###  Save samples 
    # ###

    # z.save[,1,k] <- dt[,z1]
    # z.save[,2,k] <- dt[,z2]
    # a0.save[k] <- a0    
    # sigma.save[k] <- sigma
  # }
  
  # ###
  # ### Write output
  # ###
  
  # keep$sigma <- keep$sigma/n.mcmc
  # cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
  # list(z=z.save,a0=a0.save,sigma=sigma.save,keep=keep,n.mcmc=n.mcmc)
  
# }