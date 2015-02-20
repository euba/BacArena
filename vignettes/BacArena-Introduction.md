---
title: "Introduction to BacArena"
author: "Eugen Bauer and Johannes Zimmermann"
date: "2015-02-20"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

Microbial communities are essential for global ecosystems and human health. Computational modeling of microbial consortia is thus a major goal in systems biology and microbial ecology. BacArena is a project to simulate bacterial behaviour in communities. A lot of progress is done in the last years to gain genome wide metabolic reconstructions of certain organisms, which open a wide field of mathematical analysis. One of this new methods is  fux balanced analysis (fba) to estimate optimal metabolic  fluxes under certain constraints. With this work, advanced models are build which are available in an exchangeable format (SBML). Some of these models can be found in^[ http://bigg.ucsd.edu/]. The idea of this project is to use this existing reconstructions and put them in a spatial and temporal environment to study their possible interactions. This is achieved by the combination of agent based modeling with fba. Each bacterium is considered as an agent with individual states, own properties and rules to act. Agents are located on a grid where they can move and interact via metabolic exchanges computed by fba. The starting point for our project is curiosity of what could be done with this huge models. We just throw those models into an arena to see what kind of actions will evolve.

## Installing

BacArena is avaible as a R package. The developmental version can be installed from GitHub with the devtools package:

```r
library("devtools")
```

```
## Warning: package 'devtools' was built under R version 3.1.2
```

```r
install_github("euba/BacArena", ref="rpkg")
```

```
## Downloading github repo euba/BacArena@rpkg
## Installing BacArena
## "C:/PROGRA~1/R/R-31~1.0/bin/x64/R" --vanilla CMD INSTALL  \
##   "C:/Users/eugen.bauer/AppData/Local/Temp/RtmpsXKCg9/devtools10a0188d1f49/euba-BacArena-99ccdb5"  \
##   --library="C:/Users/eugen.bauer/Documents/R/win-library/3.1"  \
##   --install-tests
```

```
## Error: Command failed (1)
```

## Getting started

First we have to load the installed package in the workspace with

```r
library("BacArena")
```

```
## Loading required package: sybil
## Loading required package: Matrix
## Loading required package: lattice
## Loading required package: RcppArmadillo
```

```
## Warning: package 'RcppArmadillo' was built under R version 3.1.2
```

```
## Loading required package: glpkAPI
## using GLPK version 4.47
```
Next we set-up the *Escherichia coli* core metabolic model

```r
ecore <- model
```
The model is integrated in BacArena by default. Alternatively you can also load your own model of interest with commands of the sybil package and libSBML. After we loaded the model, we convert it into an object of class Bac by calling the constructor

```r
bac <- Bac(ecore,deathrate=0.1,duplirate=1,
           growthlimit=0.05,growtype="exponential")
```
Now we have to set up an environment in which the organisms can interact

```r
arena <- Arena(100,100)
```
Here, we chose 100 times 100 as the grid size of the environment. Next we want to put our created organism in its environment by

```r
addOrg(arena,bac,amount=1,x=50,y=50)
```
With this command we added one individual of our bacterium in the middle of the environment (by its x and y position). Next we can add the substances to the environment

```r
addSubs(arena,40)
```
Now we added all possible substances to the environment arena with a concentration 40 mmol per gridcell. Finally we can start the simulation with

```r
eval <- simEnv(arena,20)
```

```
## iter: 1 bacs: 1 
##    user  system elapsed 
##    0.03    0.00    0.05 
##    user  system elapsed 
##    0.13    0.00    0.12 
## iter: 2 bacs: 2 
##    user  system elapsed 
##    0.03    0.00    0.03 
##    user  system elapsed 
##    0.11    0.03    0.14 
## iter: 3 bacs: 4 
##    user  system elapsed 
##    0.03    0.00    0.03 
##    user  system elapsed 
##    0.12    0.00    0.12 
## iter: 4 bacs: 8 
##    user  system elapsed 
##    0.08    0.00    0.08 
##    user  system elapsed 
##    0.13    0.00    0.12 
## iter: 5 bacs: 16 
##    user  system elapsed 
##    0.13    0.00    0.12 
##    user  system elapsed 
##    0.22    0.00    0.22 
## iter: 6 bacs: 16 
##    user  system elapsed 
##    0.14    0.00    0.14 
##    user  system elapsed 
##    0.12    0.00    0.13 
## iter: 7 bacs: 32 
##    user  system elapsed 
##    0.25    0.00    0.25 
##    user  system elapsed 
##    0.12    0.00    0.13 
## iter: 8 bacs: 32 
##    user  system elapsed 
##    0.33    0.00    0.32 
##    user  system elapsed 
##    0.13    0.00    0.12 
## iter: 9 bacs: 64 
##    user  system elapsed 
##    0.52    0.00    0.51 
##    user  system elapsed 
##    0.13    0.00    0.13 
## iter: 10 bacs: 119 
##    user  system elapsed 
##    0.81    0.00    0.81 
##    user  system elapsed 
##    0.14    0.00    0.14 
## iter: 11 bacs: 199 
##    user  system elapsed 
##    1.18    0.00    1.18 
##    user  system elapsed 
##    0.14    0.00    0.14 
## iter: 12 bacs: 203 
##    user  system elapsed 
##    1.28    0.00    1.28 
##    user  system elapsed 
##    0.14    0.00    0.14 
## iter: 13 bacs: 322 
##    user  system elapsed 
##    1.45    0.00    1.45 
##    user  system elapsed 
##    0.15    0.00    0.16 
## iter: 14 bacs: 436 
##    user  system elapsed 
##    2.05    0.00    2.04 
##    user  system elapsed 
##    0.17    0.00    0.17 
## iter: 15 bacs: 448 
##    user  system elapsed 
##    1.95    0.00    1.96 
##    user  system elapsed 
##    0.19    0.00    0.19 
## iter: 16 bacs: 602 
##    user  system elapsed 
##    2.30    0.00    2.32 
##    user  system elapsed 
##    0.21    0.00    0.20 
## iter: 17 bacs: 746 
##    user  system elapsed 
##    2.59    0.00    2.60 
##    user  system elapsed 
##    0.24    0.00    0.23 
## iter: 18 bacs: 889 
##    user  system elapsed 
##    3.41    0.00    3.42 
##    user  system elapsed 
##    0.24    0.00    0.23 
## iter: 19 bacs: 915 
##    user  system elapsed 
##    3.61    0.00    3.60 
##    user  system elapsed 
##    0.25    0.00    0.25 
## iter: 20 bacs: 1142 
##    user  system elapsed 
##    3.96    0.00    3.98 
##    user  system elapsed 
##    0.28    0.00    0.28
```
The object eval stores all 20 simulation steps, that we performed. After we retrieve the eval object we can plot now the results of the simulation

```r
plotCurves(eval)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 
This will plot the growth curve and curves of substance concentration changes over the 20 simulation steps. If we are interested in the spatial and temporal changes of our constructed population we can use 

```r
evalArena(eval,sims=c(1,20))
```

```
## Error in evalArena(eval, sims = c(1, 20)): unused argument (sims = c(1, 20))
```
This will produce multiple plots one by one for each simulation step with the spatial structure of the population (black dots represent individuals). We can also investigate the spatial change of the population together with the main subtrate glucose

```r
evalArena(eval,c("population","EX_glc(e)"),sims=c(1,20))
```

```
## Error in evalArena(eval, c("population", "EX_glc(e)"), sims = c(1, 20)): unused argument (sims = c(1, 20))
```
Here we only plot the first and the last result of the simulation steps given by the parameter sims. At the same time we can also integrate the visualization of different phenotypes into the population 

```r
evalArena(eval,c("population","EX_glc(e)"),phencol=T,sims=c(1,20))
```

```
## Error in evalArena(eval, c("population", "EX_glc(e)"), phencol = T, sims = c(1, : unused argument (sims = c(1, 20))
```
Now we can see that the periphery of the population has a different color than the individuals in the center. This indicates that individuals on the outside of the population use a different metabolism (respiration of glucose) than the center (fermentation of glucose and acetate). To visualize the differences of the apparent phenotypes we can use

```r
minePheno(eval)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 
This will create a PCA plot with the similarity of the different phenotypes. If we are interested in the definition of the phenotypes we can retreive the original phenotype matrix with

```r
pmat <- getPhenoMat(eval)
```
The object pmat carries now the different phenotypes which are defined by used exchange reactions within individuals on the population. A value of 1 means secretion, 2 means uptake and 0 means no usage of the substance of interest.
