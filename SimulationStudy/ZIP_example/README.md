## Fast Bayesian inference of Block Nearest Neighbor Gaussian models for ZIP data

We provide the R code needed to run a simulation study of NNGP and block-NNGP models for Zero-inflated Poisson (ZIP) data ( using Integrated Nested Laplace approximation (INLA). 

- [main]SimFit/runblockNNGP.R: simulate   data from a gamma distribution and fit the NNGP and blockNNGP models to the simulated data. 
- SimFit/blockNNGPfunctionREGULAR.R: auxiliary code that contains blockNNGP functions for regular blocks. We set family="gamma" in inla() function. 

  resf <- inla(f.blockNNGP, data = as.data.frame(data1), family = "zeroinflatedpoisson1")

- SimFit/blockNNGPfunctionIRREGULAR.R: auxiliary code that contains blockNNGP functions for irregular blocks. We set family="zeroinflatedpoisson1" in inla() function. Please see ZIP type 1 in https://inla.r-inla-download.org/r-inla.org/doc/likelihood/zeroinflated.pdf

  resf <- inla(f.blockNNGP, data = as.data.frame(data1), family = "zeroinflatedpoisson1")

- SimFit/Irregblock.R:  auxiliary functions to build irregular blocks.
- SimFit/blockNNGPrgeneric.R: INLA-rgeneric code for blockNNGP.
- SimFit/NNGPfunction.R: auxiliary code that contains NNGP functions. We set family="zeroinflatedpoisson1" in inla() function. 

  resf <- inla(f.NNGP, data = as.data.frame(data11), family = "zeroinflatedpoisson1")

- SimFit/NNGPrgeneric.R: INLA-rgeneric code for NNGP.
- SimFit/utils.R:  auxiliary functions.
