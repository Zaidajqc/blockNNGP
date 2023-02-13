
## Fast Bayesian inference of Block Nearest Neighbor Gaussian models for GAMMA data

We provide the R code needed to run a simulation study of NNGP and block-NNGP models for gaussian data using Integrated Nested Laplace approximation (INLA).

- [main]SimFit/runblockNNGP.R: simulate   data from a gamma distribution and fit the NNGP and blockNNGP models to the simulated data. 
- SimFit/blockNNGPfunctionREGULAR.R: auxiliary code that contains blockNNGP functions for regular blocks. We set family="gamma" in inla() function. 
- SimFit/blockNNGPfunctionIRREGULAR.R: auxiliary code that contains blockNNGP functions for irregular blocks. We set family="gamma" in inla() function. 
- SimFit/Irregblock.R:  auxiliary functions to build irregular blocks.
- SimFit/blockNNGPrgeneric.R: INLA-rgeneric code for blockNNGP.
- SimFit/NNGPfunction.R: auxiliary code that contains NNGP functions. We set family="gamma" in inla() function. 
- SimFit/NNGPrgeneric.R: INLA-rgeneric code for NNGP.
- SimFit/utils.R:  auxiliary functions.
