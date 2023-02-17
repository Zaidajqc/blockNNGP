## Description:
##
# R implementation of block-NNGP latent effect with rgeneric

# A variable theta is defined by INLA in the code to store theta=(theta1, theta2)  to provide an internal representation of the hyperparameters (sigma2 and  phi, respectively) to make numerical optimization easier.
# In order to define the block-NNGP latent effect in INLA, we need to define the next functions:
# The mean of the latent effects: mu
# The precision of the latent effects: Q(theta)
# A ‘graph’, with a binary representation of the precision matrix: W
# The initial values of the parameters.
# A log-normalizing constant.
# The log-prior of theta.

## Arguments:
##
#blockNNGP.model <- inla.rgeneric.define(inla.rgeneric.blockNNGP.model, W = W, n= n, n.blocks= n.blocks,nb =nb,ind_obs1=ind_obs1,num1=num1,indb=indb,coords.D=coords.D)

'inla.rgeneric.blockNNGPST.model' <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
    "log.prior", "quit"),
  theta = NULL) {

#  interpret.theta function will take the parameters in the internal scale and return the marginal variance and phi parameters:
  interpret.theta <- function() {
    return(
      list(sigmasq = exp(-theta[1L]),
      phi = 0.2 - (0.195)/ (1 + exp(theta[2L])),  # b - (b-a)/(1+exp(theta2))
      a = (2*exp(theta[3L])/(1+exp(theta[3L]))) -1))
  }
 # the graph function represents the entries of the precision matrix that are non-zero. W  must be passed as a sparse matrix (as defined in package Matrix) and the returned matrix must be sparse too.
  graph <- function(){
    require(Matrix)

    return(Diagonal(nrow(W), x = 1) + W)
  }


  # meancov_nn function computes the matrices Bbk and Fbk for each block k.
  meancov_nn <- function(Sigma, ind_obs, ind_neigblocks, indnum){
    invC_nbi <- chol2inv(chol( Sigma[ind_neigblocks, ind_neigblocks] ))
    B_bi     <- (Sigma[ind_obs, ind_neigblocks])%*%invC_nbi
    F_bi     <- Sigma[ind_obs, ind_obs] -(B_bi%*%Sigma[ind_neigblocks, ind_obs])
    invFbi   <- chol2inv(chol(F_bi))
    Bstar_bi <- matrix(0,length(ind_obs), dim(Sigma)[1])
    
    Bstar_bi[,ind_neigblocks] <- -B_bi
    Bstar_bi[indnum]       <- 1
    
    val     <- list( invFbi = invFbi, Bstar_bi = Bstar_bi)
    return(val)
  }

  # PrecNNGP function computes the precision matrix of blockNNGP
  PrecblockNNGP  <- function(nloc, n.blocks, Sigma, nb, ind_obs1, num1, indb){
    
    Fs_1         <- matrix(0,nloc,nloc)
    Bb           <- matrix(0,nloc,nloc)
    
    Fs_1[1:nb[1],1:nb[1]]   <- chol2inv(chol(Sigma[ind_obs1,ind_obs1]))
    
    Bstar_bi                 <- matrix(0,length(ind_obs1),nloc)
    diag(Bstar_bi[num1,num1]) <- 1
    Bb[1:nb[1],]             <- Bstar_bi
    
    #    system.time(
    for (j in 2:n.blocks){
      ress <- meancov_nn(Sigma,indb[[j-1]][[1]],indb[[j-1]][[2]],indb[[j-1]][[4]])
      Fs_1[(nb[j-1]+1):(nb[j]),(nb[j-1]+1):(nb[j])] <- ress$invFbi
      Bb[(nb[j-1]+1):(nb[j]),] <- ress$Bstar_bi
    }
    #    )
    
    Bbb <- as(Bb , "dgCMatrix")
    Fs_11 <- as(Fs_1 , "dgCMatrix")
    
    invCs        <- crossprod(Bbb,Fs_11)%*%Bbb
    
    return(invCs)
  }

  
  precision.ar1 = function(N, rho){
    Qar = matrix(0, N, N)
     diag(Qar) = 1+rho^2
    for (l in 1:(N-1)) {
      Qar[l, l+1] = -rho
      Qar[l+1, l] = -rho
    }
    Qar[1,1] = 1
    Qar[N,N] = 1
   Qar =as(Qar , "dgCMatrix")
    return(Qar)
  }

  
 # Q function defines the precision matrix which is defined in a similar way of W. Here we define the precision matrix of the blockNNGP latent effect.
  Q <- function() {
    require(Matrix)

    param <- interpret.theta()
    nu<- 0.5
    R <- (coords.D*param$phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=coords.D*param$phi, nu=nu)
    diag(R) <- 1
    C <- param$sigmasq*R
 
    Qs <- PrecblockNNGP(n, n.blocks,C,nb,ind_obs1,num1,indb)
    Qt <- precision.ar1(T, param$a) # T=3
    Q <- kronecker(Qt,Qs)
    
       return( Q )
  }

# mu function is the mean of the blockNNGP latent effect which is zero.
  mu <- function()
  {
    return(numeric(0))
  }

# log.norm.const function computes the normalizing constant, INLA computes it if numeric is zero.
  log.norm.const <- function() {
    return(numeric(0))

  }

# log.prior function computes the pdf of prior distributions for sigmasq and phi. In particular, for the marginal variance we set a gamma distribution  with parameters 1 and 0.00005,and for phi=2/range we set a uniform distribution on (a,b), where a and b are associated to the minimum and maximum distance between locations, in this case a=1 and b=30. extra terms that appear in the definition of the log-density of the prior are due to the change of variable involved. INLA works with (theta1, theta2)  internally, but the prior is set on (sigmasq, phi).
  log.prior <- function() {
    param   <- interpret.theta()
    #    a = 1
    #    b = 30
    res <- dgamma(param$sigmasq, 1, 5e-05, log = TRUE) + log(param$sigmasq) +
      (2*log(1/(0.195))) + log(param$phi-0.005) + log(0.2 - param$phi) +
      log((1+param$a)/(1-param$a)) - 2*log(2/(1-param$a))
    return(res)
  }

  # function to set the initial values of the parameters in the internal scale must be provided. This implies that the initial values of  sigmsq and phi are 1  and   # b - (b-a)/2, respectively.
  initial <- function() {
    return(c(0.1, 0.1, 0.1))
  }

  # A quit() function is called when all computations are finished before exiting the C code. In this case, we will simply return nothing.
  quit <- function() {
    return(invisible())
  }

  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

