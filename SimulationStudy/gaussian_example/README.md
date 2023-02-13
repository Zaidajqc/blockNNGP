# Fast Bayesian inference of Block Nearest Neighbor Gaussian models for gaussian data

- [main]SimFit/runblockNNGP.R: simulate  gaussian data and fit the NNGP and blockNNGP models to the simulated data. 
- SimFit/blockNNGPfunctionREGULAR.R: auxiliary code that contains blockNNGP functions for regular blocks. We set family="gaussian" in inla() function. 
- SimFit/blockNNGPfunctionIRREGULAR.R: auxiliary code that contains blockNNGP functions for irregular blocks. We set family="gaussian" in inla() function. 
- SimFit/Irregblock.R:  auxiliary functions to build irregular blocks.
- SimFit/blockNNGPrgeneric.R: INLA-rgeneric code for blockNNGP.
- SimFit/NNGPfunction.R: auxiliary code that contains NNGP functions. We set family="gaussian" in inla() function. 
- SimFit/NNGPrgeneric.R: INLA-rgeneric code for NNGP.
- SimFit/utils.R:  auxiliary functions.
