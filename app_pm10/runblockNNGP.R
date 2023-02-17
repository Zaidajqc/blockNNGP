############################################################################
## Author: Z. Quiroz, M. O. Prates, D. Dey and H. Rue
## Date: 17.09.2021
##
## Description:
##
##    The code performs a full Bayesian analysis of block-NNGP models using
##    Integrated Nested Laplace approximation (INLA). The BlockNNGP latent effect is implemented
##    using the rgeneric funcion.

## Arguments:
## case: 'irregular' for blockNNGP models with irregular blocks.
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
library(abind)

dir.save = getwd()



##--- Path for the Covariates directory
pm10.path = "/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/blockNNGP/review3STCO/block-NNGP-PM10/Cameletti2012/"


##--- Source the function for selecting the covariates for a given day
source(paste(pm10.path,"Covariates/Covariates_selector.R",sep=""))

setwd("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/blockNNGP/review3STCO/block-NNGP-PM10/Cameletti2012/")

## ################################
## Load the data
## ################################
##--- for the 24 stations and 182 days
Piemonte_data <-read.table("Piemonte_data_byday.csv",header=TRUE,sep=",")
coordinates <-read.table("coordinates.csv",header=TRUE,sep=",")


rownames(coordinates) = coordinates[,"Station.ID"]


##--- Borders of Piemonte (in km)
borders <-read.table("Piemonte_borders.csv",header=TRUE,sep=",")

## ################################
## Work out how many stations and days there are
## ################################
n_stations <- length(coordinates$Station.ID)
n_data <- length(Piemonte_data$Station.ID)
n_days <- as.integer(n_data/n_stations)

##--- Check that the data is OK
if (n_data %% n_stations != 0) {
    print("The number of data points needs to be an integer multiple of the number of stations!")
    return
}


## ################################
##Standardize covariates and calculate LOG of PM10
## ################################
##--- The covariates are standardised using the mean and std.dev. for
##--- the 24 data sites.
mean_covariates = apply(Piemonte_data[,3:10],2,mean)
sd_covariates = apply(Piemonte_data[,3:10],2,sd)

Piemonte_data[,3:10] =
    scale(Piemonte_data[,3:10],
          mean_covariates, sd_covariates)

Piemonte_data$logPM10 = log(Piemonte_data$PM10)

Piemonte_data$time = rep(1:n_days,each = n_stations)



##  final data

loc 	<- cbind(coordinates$UTMX, coordinates$UTMY)
X 		<- as.matrix(cbind(1, Piemonte_data[,3:10])) ## X = intercept + covariate
y 		<-  Piemonte_data$logPM10
T       <-  n_days
n        <- n_stations


#######################
## Estimation
#######################

setwd("/Users/zquiroz/Documents/ZAIDA_MAC/PAPERS/blockNNGP/review3STCO/block-NNGP-PM10/")
dir.save = getwd()


#######################
## functions required
#######################

source("blockNNGPfunctionIRREGULAR.R")
source("blockNNGPrgenericST.R")
source("Irregblock.R")
source("utils.R")


##########################################
#### Run irregular block-NNGP models 
#### for now you can only set 
#### nexp=2,3,4,5,6,7
#### n.blocks = 2^2, 2^3,2^4,2^5,2^6,2^7 
#########################################


case='irregular'


#number of blocks (nexp=7 ---> n.blocks=2^2)
nexp 	 <- 2
n.blocks <- 2^nexp
#number of neigbhor blocks
num.nb   <- 2
res1	 <- blockNNGP(case, loc,  y, X,  dir.save, n.blocks, num.nb,T, borders)
#summary(res1)


#number of blocks (nexp=7 ---> n.blocks=2^2)
nexp      <- 2
n.blocks <- 2^nexp
#number of neigbhor blocks
num.nb   <- 1
res2     <- blockNNGP(case, loc,  y, X,  dir.save, n.blocks, num.nb,T, borders)
#summary(res2)

#number of blocks (nexp=7 ---> n.blocks=2^3)
nexp      <- 3
n.blocks <- 2^nexp
#number of neigbhor blocks
num.nb   <- 2
res3     <- blockNNGP(case, loc,  y, X,  dir.save, n.blocks, num.nb,T, borders)
#summary(res3)


#number of blocks (nexp=7 ---> n.blocks=2^3)
nexp      <- 3
n.blocks <- 2^nexp
#number of neigbhor blocks
num.nb   <- 1
res4     <- blockNNGP(case, loc,  y, X,  dir.save, n.blocks, num.nb,T, borders)
#summary(res4)

