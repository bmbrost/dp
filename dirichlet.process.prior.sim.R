###
### Simulation 1-dimensional Dirichlet process prior (Gelman et al. 2014, section 23.2)
###

a0 <- 20  # concentration parameter
P0 <- seq(0,6,1)  # base probability measure
n <- 100  # number of observations

# Generate data according to stick-breaking process

# Discrete uniform base distribution
N <- length(P0)  # maximum number of clusters
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))

theta <- sample(P0,N,replace=FALSE)  # clusters randomly drawn from P0
k <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments of observations
hist(k,breaks=P0)
k.tab <- table(k)

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

source("/Users/bmb/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
source("/Users/brost/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
hist(rgamma(10000,2,.2),breaks=100)
out1 <- dirichlet.process.prior.mcmc(k,P0,priors=list(a=2,b=2),tune=list(a0=1.5),
   start=list(a0=a0),n.mcmc=10000)
# boxplot(c(out1$P)*n~rep(as.numeric(names(k.tab)),out1$n.mcmc),pch=19,cex=0.25)
plot(rep(as.numeric(names(k.tab)),out1$n.mcmc),c(out1$P)*n,pch=19,cex=0.25,col=rgb(0,0,0,0.05))
points(as.numeric(names(k.tab)),k.tab,col="red")
points(as.numeric(names(k.tab)),apply(out1$P*n,1,mean),col="blue")

hist(out1$a0,breaks=100);abline(v=a0,col=2,lty=2)
mean(out1$a0)
