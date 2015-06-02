rm(list=ls())

###
### Simulation 1-dimensional Dirichlet process mixture
###


# Following Escobar (1994)

a0 <- 2  # concentration parameter
P0 <- c(0,100)  # upper and lower limits to base probability measure
n <- 1000  # number of observations

N <- 50  # maximum number of clusters for truncation approximation to DPM
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))
plot(pie)

theta <- runif(N,min(P0),max(P0))  # clusters randomly drawn from P0
# z <- sample(1:length(theta),n,replace=TRUE,prob=pie)  # cluster assignments for observations
z <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments for observations
# hist(theta[z],breaks=1000)
hist(z,breaks=1000)

z.tab <- table(z)
z.tab

# Generate observations condition on theta, i.e., [y_i|theta_i]
y <- rnorm(n,z,1)
plot(y,rep(1,n),col=rgb(0,0,0,0.25),pch=19)
points(z,rep(1,n),col="red",pch=19)

# Make sure P0 contains support of y; redefine P0 if necessary
range(y)
# P0 <- range(y) + c(-5,5)

idx <- 1
y.tmp <- y[idx]
theta.tmp <- setdiff(z,z[idx])
denom <- a0+sum(dnorm(y[idx],theta.tmp,1))
p <- c(dnorm(y[idx],theta.tmp,1)/denom,a0/denom)
z[idx] <- sample(c(theta.tmp,rnorm(1,y[idx],1)),1,prob=p)


source("/Users/brost/Documents/git/Haulouts/dirichlet.process.mixture.mcmc.R")
out1 <- dirichlet.process.mixture.mcmc(y,P0,priors=list(a=0.01,b=0.1),tune=list(a0=0.1),
  start=list(a0=a0,z=y),n.mcmc=1000)
hist(out1$z,breaks=1000,ylim=c(0,1000))
points(y,rep(-15,n),col=rgb(0,0,0,0.25),pch="|",cex=0.75)
points(z,rep(-15,n),col="red",pch="|",cex=0.75)












# Using blocked Gibbs sampler per Gelman et al. 2014, section 23.2
a0 <- 1  # concentration parameter
P0 <- seq(0,100,1)  # base probability measure; discrete uniform distribution 
n <- 100  # number of observations

# Generate data according to stick-breaking process

# Discrete uniform base distribution
N <- 20  # maximum number of clusters for truncation approximation to DPM
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))
theta <- sample(P0,N,replace=FALSE)  # clusters randomly drawn from P0
k <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments of observations
hist(k,breaks=P0)
tab.k <- table(k)

# Continuous uniform base distribution
N <- 20  # maximum number of clusters for truncation approximation to DPM
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))

theta <- runif(N,min(P0),max(P0))  # clusters randomly drawn from P0
k <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments of observations
hist(k,breaks=1000)
# hist(k,breaks=P0)  # for comparison to discrete base distribution
k.tab <- table(k)



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
