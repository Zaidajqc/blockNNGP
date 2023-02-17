# Fast Bayesian inference of Block Nearest Neighbor Gaussian models for large data

## Zaida Quiroz, Marcos Prates, Dipak Dey and Håvard Rue.

We provide the R code needed to run spatio-temporal block-NNGP models using Integrated Nested Laplace approximation (INLA). 

All results are obtained using R and the packages "fields", "lattice", "akima", "Matrix", "slam", "igraph", "coda", "MBA", , "mvtnorm", "ggforce", "Rcpp", "tidyverse", "raster" and "abind" all available on CRAN. You have to add "INLA"  from https://www.r-inla.org/download-install.  Main code (with example usage) is denoted by [main]. 


- [main]SimFit/runblockNNGP.R: call pm10 data and fit spatio-temporal blockNNGP model. 
- SimFit/blockNNGPfunctionIRREGULAR.R: auxiliary code that contains blockNNGP functions for irregular blocks. 
- SimFit/Irregblock.R:  auxiliary functions to build irregular blocks.
- SimFit/blockNNGPSTrgeneric.R: INLA-rgeneric code for spatio-temporal blockNNGP random effect.
- SimFit/utils.R:  auxiliary functions.


```
Zaida Quiroz
Science Department
School of Mathematics
Pontificia Universidad Católica del Perú
Lima, Perú
E-mail: zquiroz@pucp.edu.pe
```
