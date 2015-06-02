###
### Simulation 1-dimensional Dirichlet process
###

a.0 <- 0.1  # strength parameter
K <- 100  # total number of clusters 
G.0 <- seq(1,K,1)  # discrete base distribution 

# Generate data according to Chinese restaurant process
n <- 100
k <- numeric(n)

k[1] <- sample(G.0,1)
for (i in 2:n) {
  tab <- table(k)
  idx <- as.numeric(names(tab))>0
  occ <- as.numeric(names(tab)[idx])
  p <- rep(a.0/(a.0+n-1),100)
  p[occ] <- tab[idx]/(a.0+n-1)
  k[i] <- sample(G.0,1,prob=p)
}

k
hist(k,breaks=100)
