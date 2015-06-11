rm(list=ls())

library(cluster)

###
### Simulation 1-dimensional Dirichlet process mixture
###

# Using stick-breaking process
a0 <- 2  # concentration parameter
P0 <- c(0,100)  # upper and lower limits to base probability measure
n <- 100  # number of observations

a0*log(n)  # approximates expected number of components

N <- 50  # maximum number of clusters for truncation approximation
  # to DPM (Gelman et al. 2014, section 23.2)
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))
plot(pie)

theta <- runif(N,min(P0),max(P0))  # clusters randomly drawn from P0
z <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments for observations
hist(z,breaks=1000)
z.tab <- table(z)
z.tab
length(z.tab)

# Generate observations condition on theta, i.e., [y_i|theta_i]
sigma <- 1
y <- rnorm(n,z,sigma)
plot(y,rep(1,n),col=rgb(0,0,0,0.25),pch=19)
points(z,rep(1,n),col="red",pch=19)


###
### Fit model
###

# Fit model according to Escobar (1994), Eq. 3: Gibbs update
# Same as Neal (2000), Algorithm 2
source("/Users/brost/Documents/git/DPMixtures/dpmixture.escobar.1994.mcmc.R")
out1 <- dpmixture.escobar.1994.mcmc(y,P0,start=list(a0=a0,z=y),n.mcmc=1000)

# Fit model according to Neal (2000), Algorithm 2: Gibbs group update of z
source("/Users/brost/Documents/git/DPMixtures/dpmixture.neal.2000.algm.2.mcmc.R")
out2 <- dpmixture.neal.2000.algm.2.mcmc(y,P0,start=list(a0=a0,z=y),n.mcmc=1000)

# Fit model according to Neal (2000), Algorithm 5: MH update
source("/Users/brost/Documents/git/DPMixtures/dpmixture.neal.2000.algm.5.mcmc.R")
out3 <- dpmixture.neal.2000.algm.5.mcmc(y,P0,tune=list(z=0.5),
  start=list(a0=a0,z=y),n.mcmc=1000)

# Fit model according to Neal (2000), Algorithm 7: MH update
source("/Users/brost/Documents/git/DPMixtures/dpmixture.neal.2000.algm.7.mcmc.R")
out4 <- dpmixture.neal.2000.algm.7.mcmc(y,P0,tune=list(z=0.5),
  start=list(a0=a0,z=y),n.mcmc=1000)

# Fit model according to Neal (2000), Algorithm 8: auxilliary variables
source("/Users/brost/Documents/git/DPMixtures/dpmixture.neal.2000.algm.8.mcmc.R")
out5 <- dpmixture.neal.2000.algm.8.mcmc(y,P0,priors=list(m=3),tune=list(z=0.5),
  start=list(a0=a0,z=z),n.mcmc=1000)

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/DPMixtures/dpmixture.blocked.mcmc.R")
# hist(rgamma(1000,2,2),breaks=100)
# hist(rgamma(1000,1,1),breaks=100)
out6 <- dpmixture.blocked.mcmc(y,P0,
    priors=list(H=50,r=20,q=10,sigma.l=0,sigma.u=5),
    tune=list(z=0.5,sigma=0.1),
    start=list(a0=a0,z=fitted(kmeans(y,rpois(1,10))),pie=pie,sigma=sigma),
    n.mcmc=1000)

mod <- out1
mod <- out2
mod <- out3
mod <- out4
mod <- out5
mod <- out6
idx <- 1:1000
idx <- 1:5000
idx <- 9000:10000
hist(z,breaks=1000,prob=TRUE,col="red",border="red",ylim=c(0,1))
hist(mod$z[,idx],breaks=1000,prob=TRUE,add=TRUE)
points(y,rep(-0.01,n),col=rgb(0,0,0,0.25),pch="|",cex=0.75)
points(z,rep(-0.01,n),col="red",pch="|",cex=0.75)
hist(mod$a0[idx],breaks=100);abline(v=a0,col=2,lty=2)
mean(mod$a0[idx])*log(n)
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)

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

plot(apply(mod$z[,idx],2,function(x) length(unique(x))),type="l")
abline(h=length(z.tab),col=2,lty=2)

which.max(apply(mod$z,1,var))
abline(v=29.75)
which.min(sapply(y,function(x) dist(c(x,29.75))))










###
### Fit model
###

source("/Users/brost/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
out1 <- dirichlet.process.prior.mcmc(k,P0,priors=list(a=0.01,b=0.1),tune=list(a0=0.1),
                                     start=list(a0=a0),n.mcmc=1000)
tab <- hist(k,breaks=0:n,plot=FALSE)
boxplot(c(out1$P)*n~rep(1:n,out1$n.mcmc),pch=19,cex=0.25)
points(1:n,tab$counts,col="red")
points(1:n,apply(out1$P*n,1,mean),col="blue")





# # Generate data according to Chinese restaurant process (Gelman et al. 2014, section 23.2)
# a0 <- 0.1  # concentration parameter
# K <- 100  # total number of clusters 
# P0 <- seq(1,K,1)  # base probability measure; discrete uniform distribution 
# N <- 20  # maximum number of clusters
# 
# n <- 100
# k <- numeric(n)
# 
# k[1] <- sample(P0,1)
# for (i in 2:n) {
#   tab <- table(k)
#   idx <- as.numeric(names(tab))>0
#   occ <- as.numeric(names(tab)[idx])
#   p <- rep(a0/(a0+n-1),100)
#   p[occ] <- tab[idx]/(a0+n-1)
#   k[i] <- sample(P0,1,prob=p)
# }
# 
# hist(k,breaks=0:n)
# 
# 
# ###
# ### Fit model
# ###
# 
# source("/Users/brost/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
# out1 <- dirichlet.process.prior.mcmc(k,P0,priors=list(a=0.01,b=0.1),tune=list(a0=0.1),
#    start=list(a0=a0),n.mcmc=1000)
# tab <- hist(k,breaks=0:n,plot=FALSE)
# boxplot(c(out1$P)*n~rep(1:n,out1$n.mcmc),pch=19,cex=0.25)
# points(1:n,tab$counts,col="red")
# points(1:n,apply(out1$P*n,1,mean),col="blue")
# 
