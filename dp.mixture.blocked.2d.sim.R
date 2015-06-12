rm(list=ls())

library(cluster)
library(msm)
library(pscl)

source("/Users/brost/Documents/git/DPMixtures/dp.utils.R")  # sim functions

###
### Simulation 2-dimensional Dirichlet process mixture
###

# Using stick-breaking process (Gelman et al. 2014, BDA Section 23.2)
n <- 100  # number of observations to simulate
a0 <- 2  # concentration parameter
P0 <- cbind(c(-2,2,2,-2,-2),c(-2,-2,2,2,-2))  # 2-D uniform base probability measure
H <- 50  # maximum number of clusters for truncation approximation
sim1 <- stick.2d(n,P0,a0,H)
z <- sim1$z  # cluster assignments

# Generate observations conditional on cluster assignments
sigma <- 0.1
y <- matrix(rnorm(n*2,z,sigma),,2)

b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=range(P0[,1])+b,ylim=range(P0[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=P0[,1],y=P0[,2],col="gray85")
points(z,pch=19,cex=1,col=rgb(1,0,0,0.15))
points(y,pch=19,cex=0.5)


###
### Fit models
###

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/DPMixtures/dp.mixture.blocked.2d.mcmc.R")
# hist(rgamma(1000,2,2),breaks=100)
# hist(rgamma(1000,1,1),breaks=100)
start <- list(a0=a0,z=z,#z=fitted(kmeans(y,rpois(1,10))),
  sigma=sigma,pie=rdirichlet(1,rep(1/H,H)))
out1 <- dpmixture.blocked.2d.mcmc(y,P0,
  priors=list(H=H,r=20,q=10,sigma.l=0,sigma.u=5),
  tune=list(z=0.5,sigma=0.1),start=start,n.mcmc=100)

mod <- out1
idx <- 1:1000
idx <- 1:5000
idx <- 9000:10000

# True clusters
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=range(P0[,1])+b,ylim=range(P0[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=P0[,1],y=P0[,2],col="gray85")
points(z,pch=19,cex=1,col=rgb(1,0,0,0.15))
points(mod$z[,1,],mod$z[,2,],pch=19,cex=0.1)
points(y,pch=19,cex=0.5,col=3)



# Concentration parameter
hist(mod$a0[idx],breaks=100);abline(v=a0,col=2,lty=2) 
mean(mod$a0[idx])*log(n)

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)

# Modeled number of clusters
plot(apply(mod$z[,idx],2,function(x) length(unique(x))),type="l")  
abline(h=length(z.tab),col=2,lty=2)  # true number of clusters


plot(apply(mod$z[,idx],2,max),type="l")
cl.ranks <- apply(mod$z[,idx],2,dense_rank)
plot(c(mod$z[,idx])[c(cl.ranks)==8],type="l")

hist(c(mod$z[,idx])[c(test)==1],breaks=500,xlim=range(mod$z[idx]),ylim=c(0,5),prob=TRUE)  
hist(c(mod$z[,idx])[c(test)==5],col=5,breaks=500,add=TRUE,prob=TRUE)  

pt.idx <- 94
plot(mod$z[pt.idx,idx],type="l");abline(h=z[pt.idx],col="red",lty=2)
hist(mod$z[,idx],breaks=5000,xlim=c(range(mod$z[pt.idx,])+c(-10,10)),prob=TRUE)
hist(mod$z[pt.idx,idx],breaks=50,col="red",add=TRUE,border="red",prob=TRUE);abline(v=z[pt.idx],col="red",lty=2)
points(y[pt.idx],-0.010,pch=19)

which.max(apply(mod$z,1,var))
abline(v=29.75)
which.min(sapply(y,function(x) dist(c(x,29.75))))
