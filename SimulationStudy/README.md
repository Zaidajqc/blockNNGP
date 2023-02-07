# Fast Bayesian inference of Block Nearest Neighbor Gaussian models for large data

## Zaida Quiroz, Marcos Prates, Dipak Dey and Håvard Rue.

We provide the R code needed to run a simulation study of NNGP and block-NNGP models using Integrated Nested Laplace approximation (INLA). 

All results are obtained using R and the packages "fields", "lattice", "akima", "Matrix", "slam", "igraph", "coda", "MBA", , "mvtnorm", "ggforce", "Rcpp", "tidyverse" and "raster" all available on CRAN. You have to add "INLA"  from https://www.r-inla.org/download-install.  Main code (with example usage) are denoted by [main]. 


- SimFit/runblockNNGP.R: fit the NNGP and blockNNGP models to the simulated data. 
- SimFit/Irregblock.R:
"utils.R",
"blockNNGPrgeneric.R", "blockNNGPfunction.R",  "blockNNGPfunction1.R", "NNGPrgeneric.R",
"NNGPfunction.R", "runblockNNGP.R".

```
Zaida Quiroz
Science Department
School of Mathematics
Pontificia Universidad Católica del Perú
Lima, Perú
E-mail: zquiroz@pucp.edu.pe
```
