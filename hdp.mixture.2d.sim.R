# Simulate 2-dimensional data for a hierarchical Dirichlet process
# mixture model. Use hdp.mixture.2d.mcmc.R for model fitting.

rm(list=ls())

library(cluster)
# library(msm)
# library(pscl)
library(MCMCpack)  # for rdirichlet(...)
library(DPpackage)

# source("/Users/brost/Documents/git/DPMixtures/dp.utils.R")  # sim functions

#############################################################################
### Define support of cluster centroids
#############################################################################

S.tilde <- cbind(c(-2,2,2,-2,-2),c(-2,-2,2,2,-2))  # 2D uniform base probability measure


#############################################################################
### Simulate clusters and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
#############################################################################

T <- 500  # number of observations to simulate
theta <- 1.5  # Dirichlet process mixture concentration parameter
H <- 50  # maximum number of clusters for truncation approximation

# Prior elicitation for theta
E.m <- theta*(digamma(theta+T)-digamma(theta))  # expected number of clusters
theta.priors <- DPelicit(T,mean=E.m,std=3,method="JGL")$inp  # Gamma(r,q) prior for theta

# Stick-breaking process
eta <- c(rbeta(H-1,1,theta),1)  # stick-breaking weights
pie <- eta*c(1,cumprod((1-eta[-H])))  # probability mass
plot(pie,type="b")

# Simulate possible clusters from S.tilde
mu.0.tmp <- cbind(runif(H,min(S.tilde[,1]),max(S.tilde[,1])),
	runif(H,min(S.tilde[,2]),max(S.tilde[,2])))
plot(mu.0.tmp) #  possible cluster locations

# Simulate realized clusters locations
ht <- sample(1:H,T,replace=TRUE,prob=pie)  # latent cluster assignments
m <- length(unique(ht))  # number of clusters
tab <- table(ht)  # tabulate cluster membership
mu.0 <- mu.0.tmp[ht,] #  occupied clusters

# Plot realized cluster locations
points(mu.0,pch=19,cex=0.5,col=2)


#############################################################################
### Simulate observations conditional on cluster assignments
#############################################################################

sigma <- 0.25
s <- matrix(rnorm(T*2,mu.0,sigma),,2)

b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=range(S.tilde[,1])+b,ylim=range(S.tilde[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.tilde[,1],y=S.tilde[,2],col="gray85")
points(s,pch=19,cex=0.5,col=rgb(0,0,0,0.25))
points(mu.0,pch=19,cex=0.75,col=2)


#############################################################################
### Fit models
#############################################################################

start <- list(theta=theta,mu.0=mu.0,pie=pie,sigma=sigma)
priors <- list(H=H,r=theta.priors[1],q=theta.priors[2],sigma.l=0,sigma.u=5)
tune <- list(sigma=0.025)
source("/Users/brost/Documents/git/DPMixtures/hdp.mixture.2d.mcmc.R")
out1 <- hdpmixture.2d.mcmc(s,S.tilde,priors=priors,tune=tune,start=start,n.mcmc=5000)

mod <- out1
idx <- 1:5000

# True clusters
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=range(S.tilde[,1])+b,ylim=range(S.tilde[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.tilde[,1],y=S.tilde[,2],col="gray85")
points(s,pch=19,cex=0.5,col=rgb(0,0,0,0.25))
points(mod$mu.0[,1,idx],mod$mu.0[,2,idx],cex=0.01,col=3)
points(mu.0,pch=19,cex=0.75,col=2)

# Inference on theta: DP concentration parameter
hist(mod$theta[idx],breaks=100);abline(v=theta,col=2,lty=2) 
mean(mod$theta[idx])*log(T)

# Inference on m: modeled number of clusters
plot(mod$m,type="l");abline(h=m,col=2,lty=2)  # true number of clusters  
barplot(table(mod$m)) 

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)

# Examine individual observation
pt.idx <- 2
plot(mod$mu.0[pt.idx,1,idx],mod$mu.0[pt.idx,2,idx],pch=19,col=rgb(0,0,0,0.25),cex=0.25,
	ylim=range(c(mod$mu.0[pt.idx,2,idx],s[pt.idx,2])),
	xlim=range(c(mod$mu.0[pt.idx,1,idx],s[pt.idx,1])),asp=1)
points(s[pt.idx,1],s[pt.idx,2],col=2,pch=19)
points(s,col=3,cex=0.5,pch=19)
identify(s)


