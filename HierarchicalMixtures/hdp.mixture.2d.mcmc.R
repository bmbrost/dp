# Hierarchical Dirichlety process mixture model for 2-dimensional data

hdpmixture.2d.mcmc <- function(s,j,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
  
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

	mu <- start$mu
	H <- priors$H
	sigma <- start$sigma
	pie.j <- start$pie.j
	pie.0 <- start$pie.0
	theta.j <- start$theta.j
	theta.0 <- start$theta.0
 
  
	###
	###  Setup Variables 
	###
  
    # browser()  
	T <- nrow(s)  # total number of observations
	T.j <- table(j)  # number of observation per group
	J <- length(T.j)  # number of groups
	
	# Setup Dirichlet process mixture variables
	mu.0 <- unique(mu) #  unique clusters
	m0 <- nrow(mu.0)  # number of clusters
	h <- c(1:m0)[match(mu[,1],mu.0[,1])]  # idx assigning each obs in s to record in mu.0

	# Summarize classification over all groups
	n0 <- table(h)  # tabulate overall cluster membership
	ord.0 <- order(n0,decreasing=TRUE)  # clusters ordered by membership
	samp.0 <- as.numeric(names(n0))[ord.0]

	# Summarize classification by group
	n.j <- tapply(h,j,table)  # number of clusters by group
	ord.j <- lapply(n.j,order,decreasing=TRUE) #idx of tab.j ordered by membership
	samp.j <- sapply(1:J,function(x) as.numeric(names(n.j[[x]]))[ord.j[[x]]])
	m.j <- sapply(1:J,function(x) length(n.j[[x]]))  # number of clusters by group

	mu.0 <- rbind(mu.0,matrix(0,H-m0,2))  # Add unoccupied 'dummy' clusters
	
	###
	### Create receptacles for output
	###

	mu.save <- array(0,dim=c(T,2,n.mcmc))
	theta.j.save <- matrix(0,n.mcmc,J)
	theta.0.save <- numeric(n.mcmc)
	sigma.save <- numeric(n.mcmc)  
	m0.save <- numeric(n.mcmc)  
	m.j.save <- matrix(0,n.mcmc,J)  
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

		# Global stick-breaking process
		# sample in same order as pie.0, i.e., in order of decreasing membership n0

		pad <- ifelse(m0<H,H-m0-1,0)
		n.tmp <- c(n0[ord.0],rep(0,pad))  # membership in decreasing order
		eta.0 <- c(rbeta(H-1,1+n.tmp,theta.0+T-cumsum(n.tmp)),1)  # stick-breaking weights
		pie.tmp <- eta.0*c(1,cumprod((1-eta.0[-H])))  # mixture component probabilities	

		# Catalogue in same order as mu.0
		pie.0[samp.0] <- pie.tmp[1:m0]
		pie.0[-samp.0] <- pie.tmp[(m0+1):H]

		# Group-level stick-breaking process
		# sample in same order as pie.j, i.e., in order of membership n.j
		
		pad <- sapply(m.j,function(x) ifelse(x<H,H-x-1,0))	
		for(i in 1:J){
# if(k==100) browser()
			theta.tmp <- theta.j[i]
			samp.tmp <- lapply(samp.j,function(x) c(x,c(1:H)[-x]))
			samp.tmp <- samp.tmp[[i]]
			
			# within group membership in decreasing order
			# n.tmp <- c(n.j[[i]][ord.j[[i]]],rep(0,pad[i]))  	
			# eta.tmp <- c(rbeta(H-1,1+n.tmp,theta.tmp+T-cumsum(n.tmp)),1)  # stick-breaking
				 # weights
			
			# global membership in order of decreasing group membership
			idx <- sapply(samp.j[[i]],function(x) which(samp.0==x))
			n.tmp <- c(n0[ord.0][idx],n0[ord.0][-idx],rep(0,H-m0))  
			eta.tmp <- c(rbeta(H-1,theta.tmp*pie.0+n.tmp,
				theta.tmp*(1-cumsum(c(pie.0[samp.tmp],pie.0[-samp.tmp])))+T-cumsum(n.tmp)),1)
# plot(eta.tmp,type="b")
		    pie.tmp <- eta.tmp*c(1,cumprod((1-eta.tmp[-H])))  # mixture component probabilities
# plot(pie.tmp,type="b")
		
			# eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
		    # pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities
			pie.j[samp.j[[i]],i] <- pie.tmp[1:m.j[i]]
			pie.j[-samp.j[[i]],i] <- pie.tmp[(m.j[i]+1):H]
		}
# browser()
# matplot(pie.j,type="l")
# lines(pie.0)

		# Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# and Gelman et al. (2014), Section 23.3
	
		# Sample 'occupied' mu.0 (true location of clusters with non-zero membership)		 
	    # Sampler currently disregards support P0
# browser()
		
# plot(mu.0)
		mu.0[samp.0,] <- t(sapply(samp.0,function(x) get.mu.0(x,s,h,sigma)))
# points(mu.0,col=2)				
		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior 
	    mu.0[-samp.0,] <- cbind(runif(H-m0,min(S.tilde[,1]),max(S.tilde[,1])),
			runif(H-m0,min(S.tilde[,2]),max(S.tilde[,2])))

	    # Note: sampling order matters for the remaining DP updates. Cluster parameters 
	    # must be sampled in same order as pie, i.e., sorted by decreasing membership
    
	    # Sample cluster assignment indicator
		samp.tmp <- lapply(samp.j,function(x) c(x,c(1:H)[-x]))
		for(i in 1:J){
			idx <- which(j==i)
			s.tmp <- s[idx,]
			T.tmp <- T.j[i]
			samp <- samp.tmp[[i]]
				
			h[idx] <- apply(s.tmp,1,function(x) sample(1:H,1,
				prob=exp(log(pie.j[,i])+dnorm(x[1],mu.0[,1],sigma,log=TRUE)+
				dnorm(x[2],mu.0[,2],sigma,log=TRUE))))
			# h[idx] <- sapply(1:T.tmp,function(x) sample(samp,1,
				# prob=exp(log(pie[,i])+dnorm(s.tmp[x,1],mu.0[samp,1],sigma,log=TRUE)+
				# dnorm(s.tmp[x,2],mu.0[samp,2],sigma,log=TRUE))))
		}		

# browser()
# b <- 3*c(-sigma,sigma) # Plot buffer for errors
# plot(0,0,xlim=range(S.tilde[,1])+b,ylim=range(S.tilde[,2])+b,
	# pch="",yaxt="n",xaxt="n",xlab="",ylab="")
# polygon(x=S.tilde[,1],y=S.tilde[,2],col="gray85")
# segments(s[,1],s[,2],mu.0[h,1],mu.0[h,2],col="gray60")
# points(s,pch=19,cex=0.25,col=rep(2:(J+1),each=T))
# points(mu.0,pch=19,cex=0.5,col=1)

		# ht <- sapply(1:T,function(x) sample(samp,1,
			# prob=exp(log(pie)+dnorm(s[x,1],mu.0[samp,1],sigma,log=TRUE)+
			# dnorm(s[x,2],mu.0[samp,2],sigma,log=TRUE))))
		mu.tmp <- mu.0[h,]

	n0 <- table(h)  # tabulate overall cluster membership
	m0 <- length(n0)  # number of clusters

	# Summarize classification over all groups
	ord.0 <- order(n0,decreasing=TRUE)  # clusters ordered by membership
	samp.0 <- as.numeric(names(n0))[ord.0]

	# Summarize classification by group
	n.j <- tapply(h,j,table)  # number of clusters by group
	ord.j <- lapply(n.j,order,decreasing=TRUE) #idx of tab.j ordered by membership
	samp.j <- sapply(1:J,function(x) as.numeric(names(n.j[[x]]))[ord.j[[x]]])
	m.j <- sapply(1:J,function(x) length(n.j[[x]]))  # number of clusters by group


		# Tabulate cluster membership with base functions
		# tab <- table(ht) #  tabulate cluster membership
		# m <- length(tab) #  number of clusters
		# ord <- order(tab,decreasing=TRUE) #  idx of tab ordered by membership
		# samp <- as.numeric(names(tab))[ord]
		


		# Stick-breaking process
		# pad <- ifelse(m<H,H-m-1,0)
		# tab.tmp <- c(tab[ord],rep(0,pad))  # membership in decreasing order
		# eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    # pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities
  
		# Sample theta (concentration parameter)
	    # See Gelman section 23.3, Ishwaran and Zarepour (2000)
       	
       	# theta.0 <- rgamma(1,priors$r.0+H-1,priors$q.0-sum(log(1-eta.0[-H])))  

		# Escobar and West (1995) and West (1997?) white paper
		tmp <- rbeta(1,theta.0+1,T)
		c <- priors$r.0
		d <- priors$q.0
		p.tmp <- (c+m0-1)/(c+m0-1+T*(d-log(tmp)))
		p.tmp <- rbinom(1,1,p.tmp)
		theta.0 <- ifelse(p.tmp==1,rgamma(1,c+m0,d-log(tmp)),rgamma(1,c+m0-1,d-log(tmp)))

		
		
		###
		### Sample sigma (observation error)
	    ###

    	# sigma.star <- rnorm(1,sigma,tune$sigma)
	    # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
			# mh.star.sigma <- sum(dnorm(s[,1],mu.0.tmp[,1],sigma.star,log=TRUE)+
	        	# dnorm(s[,2],mu.0.tmp[,2],sigma.star,log=TRUE))
    		# mh.0.sigma <- sum(dnorm(s[,1],mu.0.tmp[,1],sigma,log=TRUE)+
		        # dnorm(s[,2],mu.0.tmp[,2],sigma,log=TRUE))
		    # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        		# sigma <- sigma.star
		        # keep$sigma <- keep$sigma+1
    		# } 
	    # }
 

	    ###
    	###  Save samples 
	    ###

    	mu.save[,,k] <- mu.tmp
		m0.save[k] <- m0
		m.j.save[k,] <- m.j
    	theta.0.save[k] <- theta.0    
    	theta.j.save[k,] <- theta.0    
	    sigma.save[k] <- sigma
	}
  
	###
	### Write output
	###
  
	keep$sigma <- keep$sigma/n.mcmc
	cat(paste("\nsigma acceptance rate:",keep$sigma)) 
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	list(mu.0=mu.save,theta.0=theta.0.save,theta.j=theta.j.save,m0=m0.save,m.j=m.j.save,
		sigma=sigma.save,keep=keep,n.mcmc=n.mcmc)
 }