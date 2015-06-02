dpmixture.neal.2000.eq.3.6.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
  # browser()
  a0 <- start$a0
  z <- start$z
  n <- length(y)  # number of observations
  P0 <- range(y)+c(-10,10)  # redefine P0 to contain support of y plus extra
  F.P0 <- 1/(max(P0)-min(P0))
  z.save <- matrix(0,n,n.mcmc)
  
  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample z (cluster assignments for observations)
    ###
    
    # Following Escobar (1994) Eq. 3 
    # Same Gibbs sampler as Neal (2000) Eq. 3.2
# browser()  
    delta.j <- table(z)  # number of observations per cluster
    phi.c <- sort(unique(z))  # unique clusters locations; same ordering as delta.j
    idx <- which(delta.j==1)
    
#     idx <- which(z%in%phi.c[delta.j==1])  # index of observations associated with no 
      # other observations
    for(i in idx){
#       i <- idx[1]
      phi.tmp <- phi.c[i]
      y.tmp <- y[which(z==phi.tmp)]
      p.old <- delta.j[-i]/(n-1+a0)*dnorm(y.tmp,phi.c[-i],1) # prob of existing phi
      p.new <- a0/(n-1+a0) * (a0*F.P0) # prob of new phi
      z[i] <- sample(c(phi.c[-1],rnorm(1,y.tmp,1)),1,prob=c(p.old,p.new))  # sample phi
    }  
    for(i in unique(z)){  # sample new phi from it's full conditional
#       i <- unique(z)[1]
      idx <- which(z==i)
      A <- length(idx)/1
      b <- sum(z[idx])
      z[idx] <- rnorm(1,1/A*b,1/A)
    }  
      
#     for(i in 1:n){ 
#       y.tmp <- y[i]
#       idx <- which(phi.c==z[i])  
#       phi.tmp <- phi.c[-idx]
# # Note: unsure if normalizing constant c.norm and p.new are being calculated
# # due to integral over normal and uniform densities
# # fxn <- function(x) {dnorm(x,y.tmp)*dunif(x,min(P0),max(P0))}
# # integrate(fxn, min(P0), max(P0))
# #       c.norm <- a0*F.P0+sum(dnorm(y.tmp,phi.tmp,1))  # normalizing constant, i.e., A(Y)
#       p.old <- dnorm(y.tmp,phi.tmp,1) #/c.norm  # prob of existing phi
#       p.new <- (a0*F.P0) #/c.norm  # prob of new phi
#       z[i] <- sample(c(phi.tmp,rnorm(1,y.tmp,1)),1,prob=c(p.old,p.new))  # sample phi
#     }
     
    ###
    ###  Save samples 
    ###
    
    z.save[,k] <- z
    
  }
  
  ###
  ### Write output
  ###
  
  list(z=z.save,n.mcmc=n.mcmc)
  
}