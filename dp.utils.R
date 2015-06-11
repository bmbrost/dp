stick <- function(n,P0,a0,H){
  
  # Simulate Dirichlet process data using stick-breaking construction
  # See Gelman et al. (2014) BDA, Section 23.2
  
  # Arguments: n = number of observations to simulate,
  # P0 = upper and lower limits of uniform base probability measure,
  # a0 = concentration parameter, H = maximum number of clusters
  # for truncation approximation.
  
  v <- c(rbeta(H-1,1,a0),1)
  pie <- v*c(1,cumprod((1-v[-N])))
  theta <- runif(H,min(P0),max(P0))  # clusters randomly drawn from P0
  z <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments for observations  
  list(z=z,theta=theta,pie=pie,v=v,n=n,P0=P0,a0=a0,H=H)  
}


crp <- function(n,P0,a0){
  
  # Simulate Dirichlet process data using Chinese restaurant process
    
  # Arguments: n = number of observations to simulate,
  # P0 = upper and lower limits of discrete uniform base probability measure,
  # a0 = concentration parameter
  
  library(data.table)
    
  z <- numeric(n)
  z <- as.data.table(z)
  z[1] <- runif(1,P0[1],P0[2])
 
  for (i in 2:n) {
    tab <- z[1:(i-1),.N,by="z"]
    p.old <- tab[,N]/(a0+n-1)
    p.new <- a0/(a0+n-1)
    z[i] <- sample(c(tab[,z],runif(1,P0[1],P0[2])),1,prob=c(p.old,p.new))
  }
  list(z=z[,z],n=n,P0=P0,a0=a0)  
}


polya <- function(n,P0,a0){
  
  # Simulate Dirichlet process data using Polya urn scheme
  
  # Arguments: n = number of observations to simulate,
  # P0 = upper and lower limits of discrete uniform base probability measure,
  # a0 = concentration parameter
  
}
