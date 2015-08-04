# BacArena

BacArena is an agent based modeling framework for cellular communities.
Existing genome-scale metabolic models can be used to simulate growth and interactions in time and space.
In BacArena two well established methods are merged: a) Flux balance analysis to predict the activity of
metabolic reactions and b) Agenda based modelling 
It has been  efficiently implemented in R language (with some C++ routines) and is freely available (soon on cran).

Features:
- Each organism is represented individually
- Simulation of >10 different species and thousands of organisms on your desktop computer
- Diffusion of substances
- Detection of different phenotypes
- Chemotaxis
- Kinetics of reactions
- Separation of simulation and evaluation
- Rich evaluation methods (data mining)
- Reproducible simulations
- Object oriented implementation
- Easily expandable due to rule based approach


## Installation

- To install the latest release:
  ``R CMD install https://github.com/euba/BacArena/releases/download/0.9/BacArena_0.9.tar.gz``

- To install the development version:
  ```R
library(devtools)
install_github("euba/bacarena")
```

## Documentation

A tutorial is available: https://github.com/euba/BacArena/blob/master/vignettes/BacArena-Introduction.pdf


## Issues

Please report bugs, disorders or features you would like to see: (https://github.com/euba/BacArena/issues)
