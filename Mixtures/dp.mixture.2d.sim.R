# Simulate 2-dimensional data for a Dirichlet process mixture.
# Use dp.mixture.2d.mcmc.R for model fitting.

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

T <- 10000  # number of observations to simulate
theta <- 5  # Dirichlet process mixture concentration parameter
H <- 100  # maximum number of clusters for truncation approximation

# Prior elicitation for theta
# theta*log(1+T/theta)  # expected number of clusters (Escobar and West 1995)
E.m <- theta*(digamma(theta+T)-digamma(theta))  # expected number of clusters
theta.priors <- DPelicit(T,mean=E.m,std=10,method="JGL")$inp  # Gamma(r,q) prior for theta

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
# priors <- list(H=H,r=1,q=0.1,sigma.l=0,sigma.u=5)
tune <- list(sigma=0.025)
source("/Users/brost/Documents/git/DPMixtures/dp.mixture.2d.mcmc.R")
out1 <- dpmixture.2d.mcmc(s,S.tilde,priors=priors,tune=tune,start=start,n.mcmc=1000)

mod <- out1
idx <- 1:1000
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
pt.idx <- 2000
plot(mod$mu.0[pt.idx,1,idx],mod$mu.0[pt.idx,2,idx],pch=19,col=rgb(0,0,0,0.25),cex=0.25,
	ylim=range(c(mod$mu.0[pt.idx,2,idx],s[pt.idx,2])),
	xlim=range(c(mod$mu.0[pt.idx,1,idx],s[pt.idx,1])),asp=1)
points(s[pt.idx,1],s[pt.idx,2],col=2,pch=19)
points(s,col=3,cex=0.5,pch=19)
identify(s)




####
#### Use the simulation code below for the version of the MCMC algorithm
#### that's commented out at the bottom of dp.mixture.2d.mcmc.R
####

# rm(list=ls())

# library(cluster)
# # library(msm)
# # library(pscl)
# library(MCMCpack)  # for rdirichlet(...)

# source("/Users/brost/Documents/git/DPMixtures/dp.utils.R")  # sim functions

# ###
# ### Simulation 2-dimensional Dirichlet process mixture
# ###

# # Using stick-breaking process (Gelman et al. 2014, BDA Section 23.2)
# n <- 500  # number of observations to simulate
# a0 <- 1.5  # concentration parameter
# P0 <- cbind(c(-2,2,2,-2,-2),c(-2,-2,2,2,-2))  # 2-D uniform base probability measure
# H <- 50  # maximum number of clusters for truncation approximation
# sim1 <- stick.2d(n,P0,a0,H)
# z <- sim1$z  # cluster assignments
# z.tab <- table(z)

# # Generate observations conditional on cluster assignments
# sigma <- 0.1
# y <- matrix(rnorm(n*2,z,sigma),,2)

# b <- 3*c(-sigma,sigma) # Plot buffer for errors
# plot(0,0,xlim=range(P0[,1])+b,ylim=range(P0[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
# polygon(x=P0[,1],y=P0[,2],col="gray85")
# points(y,pch=19,cex=0.5)
# points(z,pch=19,cex=1,col=rgb(1,0,0,0.15))


# ###
# ### Fit models
# ###

# # Fit model using blocked Gibbs sampler 
# source("/Users/brost/Documents/git/DPMixtures/dp.mixture.blocked.2d.mcmc.R")
# # hist(rgamma(1000,2,2),breaks=100)
# # hist(rgamma(1000,1,1),breaks=100)
# start <- list(a0=a0,z=z,#z=fitted(kmeans(y,rpois(1,10))),
  # sigma=sigma,pie=rdirichlet(1,rep(1/H,H))) #sim1$pie)  # 
# start <- list(a0=a0,z=z,#z=fitted(kmeans(y,rpois(1,10))),
  # sigma=sigma,pie=sim1$pie)  # 
# out1 <- dpmixture.blocked.2d.mcmc(y,P0,
  # priors=list(H=H,r=2,q=2,sigma.l=0,sigma.u=5),
  # tune=list(z=0.5,sigma=0.01),start=start,n.mcmc=2000)

# mod <- out1
# idx <- 1:100
# idx <- 1:2000
# idx <- 1000:2000
# idx <- 1:2500
# idx <- 9000:10000

# # True clusters
# b <- 3*c(-sigma,sigma) # Plot buffer for errors
# plot(0,0,xlim=range(P0[,1])+b,ylim=range(P0[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
# polygon(x=P0[,1],y=P0[,2],col="gray85")
# points(mod$z[,1,idx],mod$z[,2,idx],pch=19,cex=0.5,col=rgb(0,0,0,0.0025))
# points(y,pch=19,cex=0.25,col=3)
# points(z,pch=19,cex=0.5,col=rgb(1,0,0,1))

# cl <- kmeans(apply(mod$z,2,I),11)
# points(cl$centers,col=4,pch=19)

# # Concentration parameter
# hist(mod$a0[idx],breaks=100,col="gray90");abline(v=a0,col=2,lty=2) 
# mean(mod$a0[idx])*log(n)

# # Observation error
# hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)

# # Modeled number of clusters
# nclust <- apply(mod$z[,,idx],c(3),function(x) nrow(unique(x)))
# plot(nclust,type="l")
# abline(h=nrow(unique(z)),col=2,lty=2)  # true number of clusters  
# barplot(table(nclust))



# plot(apply(mod$z[,idx],2,max),type="l")
# cl.ranks <- apply(mod$z[,idx],2,dense_rank)
# plot(c(mod$z[,idx])[c(cl.ranks)==8],type="l")

# hist(c(mod$z[,idx])[c(test)==1],breaks=500,xlim=range(mod$z[idx]),ylim=c(0,5),prob=TRUE)  
# hist(c(mod$z[,idx])[c(test)==5],col=5,breaks=500,add=TRUE,prob=TRUE)  

# pt.idx <- 94
# plot(mod$z[pt.idx,idx],type="l");abline(h=z[pt.idx],col="red",lty=2)
# hist(mod$z[,idx],breaks=5000,xlim=c(range(mod$z[pt.idx,])+c(-10,10)),prob=TRUE)
# hist(mod$z[pt.idx,idx],breaks=50,col="red",add=TRUE,border="red",prob=TRUE);abline(v=z[pt.idx],col="red",lty=2)
# points(y[pt.idx],-0.010,pch=19)

# which.max(apply(mod$z,1,var))
# abline(v=29.75)
# which.min(sapply(y,function(x) dist(c(x,29.75))))