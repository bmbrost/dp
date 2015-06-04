rm(list=ls())

###
### Simulation 1-dimensional Dirichlet process mixture
###

# Using stick-breaking process
a0 <- 2  # concentration parameter
P0 <- c(0,100)  # upper and lower limits to base probability measure
n <- 100  # number of observations

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

# Generate observations condition on theta, i.e., [y_i|theta_i]
y <- rnorm(n,z,1)
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


mod <- out1
mod <- out2
mod <- out3
mod <- out4
mod <- out5
idx <- 1:1000
idx <- 1:5000
idx <- 9000:10000
hist(z,breaks=1000,prob=TRUE,col="red",border="red")
hist(mod$z[,idx],breaks=1000,prob=TRUE,add=TRUE)
points(y,rep(-0.05,n),col=rgb(0,0,0,0.25),pch="|",cex=0.75)
points(z,rep(-0.05,n),col="red",pch="|",cex=0.75)

plot(apply(mod$z[,idx],2,min),type="l")
  
pt.idx <- 70
plot(mod$z[pt.idx,idx],type="l");abline(h=z[pt.idx],col="red",lty=2)
hist(mod$z[,idx],breaks=5000,xlim=c(range(mod$z[pt.idx,])+c(-10,10)),prob=TRUE)
hist(mod$z[pt.idx,idx],breaks=50,col="red",add=TRUE,border="red",prob=TRUE);abline(v=z[pt.idx],col="red",lty=2)
points(y[pt.idx],-0.010,pch=19)

apply(mod$z[,idx],2,function(x) length(unique(x)))
z.tab

which.max(apply(mod$z,1,var))
abline(v=31.5)
which.min(sapply(y,function(x) dist(c(x,70.0))))


table(mod$z[,1000])
z.tab









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
