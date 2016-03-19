# Simulate 2-dimensional data for a hierarchical Dirichlet process
# mixture model. Use hdp.mix.2d.mcmc.R for model fitting.

rm(list=ls())

# library(cluster)
# library(msm)
# library(pscl)
library(MCMCpack)  # for rdirichlet(...)
library(DPpackage)  # for gamma prior elicitation

# source("/Users/brost/Documents/git/DP/mixtures/dp.utils.R")  # sim functions

#############################################################################
### Define support of cluster centroids
#############################################################################

S.tilde <- cbind(c(-2,2,2,-2,-2),c(-2,-2,2,2,-2))  # 2D uniform base probability measure


#############################################################################
### Simulate clusters and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
#############################################################################

J <- 5  # numer of groups or individuals to simulate  
T <- 100  # number of observations to simulate per group
H <- 100  # maximum number of parent clusters for truncation approximation

# Simulate possible clusters from S.tilde
mu.0.tmp <- cbind(runif(H,min(S.tilde[,1]),max(S.tilde[,1])),
	runif(H,min(S.tilde[,2]),max(S.tilde[,2])))
plot(mu.0.tmp) #  possible cluster locations


###
### Simulate parent clusters (level 1)
###

theta.0 <- 10.0  # Dirichlet process mixture concentration parameter

# Prior elicitation for theta.0
# theta.0*log(T)
# E.m.0 <- theta.0*(digamma(theta.0+T*J)-digamma(theta.0))  # expected number of clusters
# theta.0.priors <- DPelicit(T,mean=E.m.0,std=10,method="JGL")$inp  # Gamma(r,q) prior for theta

# Stick-breaking process of parent clusters
eta.0 <- c(rbeta(H-1,1,theta.0),1)  # stick-breaking weights
pie.0 <- eta.0*c(1,cumprod(1-eta.0[-H]))  # probability mass
plot(pie.0,type="b")

# # Simulate realized clusters locations
# h.0 <- sample(1:H,T,replace=TRUE,prob=pie.0)  # latent cluster assignments
# m.0 <- length(unique(h.0))  # number of clusters
# tab.0 <- table(h.0)  # tabulate cluster membership
# mu.0 <- mu.0.tmp[h.0,] #  occupied clusters

# # Plot realized cluster locations
# points(mu.0,pch=19,cex=0.5,col=2)


###
### Simulate child clusters (level 2)
###

theta.j <- rep(1,J)

# Prior elicitation for theta.0
# E.m <- theta*(digamma(theta+J)-digamma(theta))  # expected number of clusters
# theta.priors <- DPelicit(T,mean=E.m[1],std=3,method="JGL")$inp  # Gamma(r,q) prior for theta

# Group-level stick-breaking process
eta.j <- sapply(1:J,function(x) rbeta(H,theta.j[x]*pie.0,theta.j[x]*(1-cumsum(pie.0))))  
pie.j <- apply(eta.j,2,function(x) x*c(1,cumprod(1-x[-H])))  # probability mass
matplot(pie.j,type="l",lty=1)


###
### Simulate realized clusters
###

# Observation-level values
h <- lapply(1:J,function(x) sample(1:H,T,replace=TRUE,prob=pie.j[,x])) # latent cluster
	# assignments for j=1,...,J and t=1,...,T
mu <- lapply(h,function(x) mu.0.tmp[x,])  # clusters for j=1,...,J and t=1,...,T
mu.mat <- do.call(rbind,lapply(mu,matrix,ncol=2))  # convert mu.t to matrix
plot(mu.mat)

# Group-level summaries
h.j <- lapply(h,function(x) unique(x))  # unique cluster assignments by group
m.j <- sapply(1:J,function(x) length(unique(h[[x]])))  # number of clusters by group
n.j <- lapply(h,function(x) table(x))  # tabulate cluster membership by group
ord.j <- lapply(n.j,order,decreasing=TRUE)  # idx of tab.j ordered by membership

# Overall summaries
h0 <- unique(unlist(h))  # unique cluster assignments 
mu.0 <- mu.0.tmp[h0,]  # unique clusters
m0 <- length(h0)  # total number of clusters

# Plot realized cluster locations
plot(mu.0,pch=19,cex=0.75,col=2)



#############################################################################
### Simulate observations conditional on cluster assignments
#############################################################################

sigma <- 0.25
s <- lapply(mu,function(x) matrix(rnorm(T*2,x,sigma),,2))
s.mat <- do.call(rbind,lapply(s,matrix,ncol=2))

b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=range(S.tilde[,1])+b,ylim=range(S.tilde[,2])+b,
	pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.tilde[,1],y=S.tilde[,2],col="gray85")
segments(s.mat[,1],s.mat[,2],mu.mat[,1],mu.mat[,2],col="gray60")
points(s.mat,pch=19,cex=0.25,col=rep(2:(J+1),each=T))
points(mu.0,pch=19,cex=0.5,col=1)


#############################################################################
### Fit models
#############################################################################
E.m <- theta.0*(digamma(theta.0+T*J)-digamma(theta.0))  # expected number of clusters
theta.priors <- DPelicit(T*J,mean=10,std=3,method="JGL")$inp  # Gamma(r,q) prior for theta


start <- list(theta.j=theta.j,theta.0=theta.0,mu=mu.mat,pie.j=pie.j,pie.0=pie.0,sigma=sigma)
# priors <- list(H=H,r=theta.priors[1],q=theta.priors[2],sigma.l=0,sigma.u=5)
hist(rgamma(1000,1,1))
priors <- list(H=H,sigma.l=0,sigma.u=5,
	r.0=1.1,q.0=0.1,  # priors for theta.0; see Teh et al. 2004, HPD UC-Berkeley report
	r=1,q=1)  # priors for theta; see Teh et al. 2004, HPD UC-Berkeley report
tune <- list(sigma=0.025)
source("/Users/brost/Documents/git/DP/hierarchical/hdp.mix.2d.mcmc.R")
out1 <- hdp.mix.2d.mcmc(s.mat,j=rep(1:J,each=T),S.tilde,
	priors=priors,tune=tune,start=start,n.mcmc=1000)

mod <- out1
idx <- 1:1000

# True clusters
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=range(S.tilde[,1])+b,ylim=range(S.tilde[,2])+b,
	pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.tilde[,1],y=S.tilde[,2],col="gray85")
points(s.mat,pch=19,cex=0.25,col=rep(2:(J+1),each=T))
points(mod$mu[,1,idx],mod$mu[,2,idx],cex=0.01,col=3)
points(mu.0,pch=19,cex=0.75,col=2)

pt.idx <- 205
pt.idx <- 1:100+00
plot(0,0,xlim=range(S.tilde[,1])+b,ylim=range(S.tilde[,2])+b,
	pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.tilde[,1],y=S.tilde[,2],col="gray85")
points(s.mat[,1],s.mat[,2],pch=19,cex=0.1,col="grey")
points(s.mat[pt.idx,1],s.mat[pt.idx,2],pch=19,cex=0.25,col=1)
points(mod$mu[pt.idx,1,idx],mod$mu[pt.idx,2,idx],cex=0.01,col=3)
points(mu.0,pch=19,cex=0.75,col=2)

# Inference on theta: DP concentration parameter
hist(mod$theta.0[idx],breaks=100);abline(v=theta.0,col=2,lty=2) 
mean(mod$theta.0[idx])*log(T)

# Inference on m0: number of global clusters
plot(mod$m0,type="l");abline(h=m0,col=2,lty=2)  # true number of clusters  
barplot(table(mod$m0)) 

# Inference on m.j: number of clusters per group
# library(lattice)
idx.tmp <- 1
plot(mod$m.j[,idx.tmp],type="l");abline(h=m.j[idx.tmp],col=2,lty=2)  # true number of clusters  
barplot(table(mod$m.j[,idx.tmp])) 


# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)

# Examine individual observation
pt.idx <- 3
plot(mod$mu[pt.idx,1,idx],mod$mu[pt.idx,2,idx],pch=19,col=rgb(0,0,0,0.25),cex=0.25,
	ylim=range(c(mod$mu[pt.idx,2,idx],s.mat[pt.idx,2])),
	xlim=range(c(mod$mu[pt.idx,1,idx],s.mat[pt.idx,1])),asp=1)
points(s.mat[pt.idx,1],s.mat[pt.idx,2],col=2,pch=19)
points(s.mat,col=3,cex=0.5,pch=19)
identify(s.mat)


