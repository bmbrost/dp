# Hierarchical Dirichlety process mixture model for 2-dimensional data

hdpmixture.2d.mcmc <- function(s,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
  
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
		mu.0 <- t(sapply(samp,function(x) get.mu.0(x,s,ht,sigma)))
		
		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior 
	    mu.0 <- rbind(mu.0,cbind(runif(H-m,min(S.tilde[,1]),max(S.tilde[,1])),
			runif(H-m,min(S.tilde[,2]),max(S.tilde[,2]))))

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
       	
       	theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  

		
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