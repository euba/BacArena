setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
SYBIL_SETTINGS("SOLVER", "glpkAPI")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

seed <- sample(1:1000,1)
cat("seed:", seed)
set.seed(seed)

data(Ec_core)

bac = Bac(model=Ec_core, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
            speed=1, type="ecore")

arena = Arena(n=10, m=10, stir=F)
addOrg(arena, bac, amount=1)
addSubs(arena, smax=20, difunc="cpp", difspeed=1)
res <- simEnv(arena, time=50)

plotCurves(res)
plotCurves2(res)

plotTotFlux(res)

  