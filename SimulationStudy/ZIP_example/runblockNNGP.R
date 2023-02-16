############################################################################
## Author: Z. Quiroz, M. O. Prates, D. Dey and H. Rue
## Date: 17.09.2021
##
## Description:
##
##    The code performs a full Bayesian analysis of NNGP and block-NNGP models using
##    Integrated Nested Laplace approximation (INLA). The BlockNNGP latent effect is implemented
##    using the rgeneric funcion.

## Arguments:
## case:  'simNNGP' for NNGP models, 'regular' for blockNNGP models with regular blocks, 
## and 'irregular' for blockNNGP models with irregular blocks. 
## formula:  NNGP  runs  ‘inla’ formula like
##  ‘y ~ 1 + x + f(idx, model = NNGP.model)’ 
## formula: blockNNGP_reg and blockNNGP  run  ‘inla’ formula like
##  ‘y ~ 1 + x + f(idx, model = blockNNGP.model)’ 
## family: A string indicating the likelihood family. For a list of possible
##          alternatives and use ‘inla.doc’ for detailed docs for
##          individual families.
## loc:  locations
## y: observed data
## X: covariates
## n.blocks: number of  blocks  (regular or irregular)
## num.nb: number of neighbors or neighbor blocks. 
##  Value:
##      ‘inla’ returns an object of class ‘"inla"’. 
##
#############################################################################

rm(list=ls())


setwd("//Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/blockNNGP/review3STCO/block-NNGP-SIM3/")


#######################
## functions required
####################### 

source("blockNNGPfunctionREGULAR.R")
source("blockNNGPfunctionIRREGULAR.R")
source("NNGPfunction.R")
source("blockNNGPrgeneric.R")
source("NNGPrgeneric.R")
source("Irregblock.R")
source("utils.R")

#######################
## libraries required
####################### 
library(INLA)
library(fields)
library(lattice)
library(akima) 
library(Matrix)
library(slam)
library(igraph)
library(coda)
library(MBA)
library(mvtnorm)
library(ggforce)
library(Rcpp)
library(tidyverse)
library(raster)


dir.save = getwd()


## simulate data
n 	<- 2500



loc 	<- cbind(runif(n,0,1), runif(n,0,1))


##Matern##
sigma.sq <- 1
nu 	 <- 0.5
phi 	 <- 12

range    <- sqrt(8 * nu) / phi
D 	 <- as.matrix(dist(loc))
R 	 <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
diag(R)  <- 1


nloc 	 <- dim(loc)[1]

C 	 <- sigma.sq*R

D 	 	<- chol(C)
rnorm_n.obs 	<- rnorm(n)
w 		<- t(matrix(rnorm_n.obs, ncol=(n))%*%D )

X 		<- as.matrix(cbind(1, rnorm(nloc))) ## X = intercept + covariate

B 		<- as.matrix(c(1,5))
p 		<- length(B)

eta         <- X%*%B + w

prob = 0.2

lambda = exp(eta)
y = rpois(n, lambda = lambda)
y[ rbinom(n, size=1, prob=prob) == 1 ] = 0

# for more information of this Zero inflated Poisson model visit
# https://inla.r-inla-download.org/r-inla.org/doc/likelihood/zeroinflated.pdf


########################################
#### Run regular block-NNGP models 
#######################################

case='regular'


#number of blocks (n.partition=8 ---> n.blocks= 8^2 = 64 blocks)
n.partition 	<- 8
n.blocks    	<- n.partition^2
#number of neigbhor blocks
num.nb 	<- 2
res1 	 	<- blockNNGP_reg(case, loc,  y, X, w,  dir.save, n.blocks, num.nb)
#summary(res1)


##########################################
#### Run irregular block-NNGP models 
#### for now you can only set 
#### nexp=2,3,4,5,6,7
#### n.blocks = 2^2, 2^3,2^4,2^5,2^6,2^7 
#########################################

case='irregular'


#number of blocks (nexp=7 ---> n.blocks=2^7)
nexp 	 <- 7
n.blocks <- 2^nexp
#number of neigbhor blocks
num.nb   <- 2
res2	 <- blockNNGP(case, loc,  y, X, w, dir.save, n.blocks, num.nb)
#summary(res2)

#number of blocks (nexp=6 ---> n.blocks=2^6)
nexp 	 <- 6
n.blocks <- 2^nexp
#number of neigbhor blocks
num.nb   <- 2
res3 	 <- blockNNGP(case, loc,  y, X, w,  dir.save, n.blocks, num.nb)
#summary(res3)


########################################
#### Run NNGP models 
#######################################

case='NNGP'


#number of neigbhors
num.nb   <- 20
res4 	 <- NNGP(case, loc,  y, X, w,   dir.save,  num.nb)
#summary(res4)



