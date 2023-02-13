# Fast Bayesian inference of Block Nearest Neighbor Gaussian models for large data

## Zaida Quiroz, Marcos Prates, Dipak Dey and Håvard Rue.

We provide the R code needed to run a simulation study of NNGP and block-NNGP models using Integrated Nested Laplace approximation (INLA). 

All results are obtained using R and the packages "fields", "lattice", "akima", "Matrix", "slam", "igraph", "coda", "MBA", , "mvtnorm", "ggforce", "Rcpp", "tidyverse" and "raster" all available on CRAN. You have to add "INLA"  from https://www.r-inla.org/download-install.  Main code (with example usage) is denoted by [main]. 


- [main]SimFit/runblockNNGP.R: simulate  data and fit the NNGP and blockNNGP models to the simulated data. The user can simulate data from gaussian or other distributions available in INLA, please see https://www.r-inla.org/documentation#h.ldt45kpbu659 for some examples. And use the following command in R: inla.list.models()
- SimFit/blockNNGPfunctionREGULAR.R: auxiliary code that contains blockNNGP functions for regular blocks. For simulated gaussian data set family="gaussian" in inla() function and for simulated non-gaussian data change the family in inla() function. 
- SimFit/blockNNGPfunctionIRREGULAR.R: auxiliary code that contains blockNNGP functions for irregular blocks.  For simulated gaussian data set family="gaussian" in inla() function and for simulated non-gaussian data change the family in inla() function.
- SimFit/Irregblock.R:  auxiliary functions to build irregular blocks.
- SimFit/blockNNGPrgeneric.R: INLA-rgeneric code for blockNNGP.
- SimFit/NNGPfunction.R: auxiliary code that contains NNGP functions.  For simulated gaussian data set family="gaussian" in inla() function and for simulated non-gaussian data change the family in inla() function.
- SimFit/NNGPrgeneric.R: INLA-rgeneric code for NNGP.
- SimFit/utils.R:  auxiliary functions.


```
Zaida Quiroz
Science Department
School of Mathematics
Pontificia Universidad Católica del Perú
Lima, Perú
E-mail: zquiroz@pucp.edu.pe
```
