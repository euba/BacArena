Title:    BacArena: Agent based modeling framework for cellular communities  
Author:   Eugen Bauer and Johannes Zimmermann

# BacArena <img src="man/bacarena.png" align="right" height = 150/>
[![DOI:10.1371/journal.pcbi.1005544](https://zenodo.org/badge/DOI/10.1371/journal.pcbi.1005544.svg)](https://doi.org/10.1371/journal.pcbi.1005544)  
[![CRAN version](http://www.r-pkg.org/badges/version/BacArena?color=blue)](https://cran.r-project.org/package=BacArena)
![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/BacArena)
![CRAN downloads/month](https://cranlogs.r-pkg.org/badges/BacArena)
[![twitter](https://img.shields.io/twitter/follow/_bacarena?style=social&logo=twitter)](https://twitter.com/_bacarena)

BacArena is an agent based modeling framework for cellular communities.

Existing genome-scale metabolic models can be used to simulate growth and interactions in time and space.
In BacArena two well established methods are merged: a) Flux balance analysis to predict the activity of
metabolic reactions and b) Agent based modelling in order to provide an environment.

It has been  efficiently implemented in R language (with some C++ routines) and is freely available [CRAN](https://cran.r-project.org/package=BacArena).

Features:
- Each organism is represented individually
- Simulation of >10 different species and thousands of organisms on your desktop computer
- Diffusion of substances
- Screening of phenotypes
- Detection of crossfeeding
- Chemotaxis
- Kinetics of reactions
- Separation of simulation and evaluation
- Rich evaluation methods (data mining)
- Reproducible simulations
- Object oriented implementation
- Easily expandable due to rule based approach


## Installation

- Install the latest CRAN release: 
```r
install.packages("BacArena")
```

- Install the development version:
```r
library(devtools)
install_github("euba/bacarena")
```

- Special hints for linux user:
  - glpk header files needed, e.g. for debian install package: ``libglpk-dev``

- Special hints for windows user:
  - Besides [R](https://cran.r-project.org/bin/windows/base/) you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)


## Quick start
```r
library("BacArena")
openArena()
```

## matlab models
- cobra matlab model files can be imported
- [tutorial](https://gist.github.com/jotech/2ec33f33fb86a400fb40b816277f4147)
```r
readMATmod("model.mat")
```

## SBML support
- sybilSBML is needed for SBML input
- currently sybilSBML is not available on CRAN because the CRAN test servers do not have the latest version of libsbml yet
- manual installation of sybilSBML (for linux):
```
wget https://www.cs.hhu.de/fileadmin/redaktion/Oeffentliche_Medien/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Informatik/Bioinformatik/sybilSBML_3.0.5.tar.gz
R CMD INSTALL  sybilSBML_3.0.5.tar.gz
```


## Documentation

- A tutorial is available: [Introduction](https://CRAN.R-project.org/package=BacArena/vignettes/BacArena-Introduction.pdf) 
- Short tutorials are on [github gist](https://gist.github.com/jotech)


## Issues

Please report bugs, disorders or features you would like to see: [Issues](https://github.com/euba/BacArena/issues)
