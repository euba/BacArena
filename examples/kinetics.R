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

bac <- Bac(model=Ec_core, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=1, type="ecore")

setKinetics(bac, exchangeR="EX_glc(e)", Km=5, vmax=10)

arena = Arena(n=10, m=10, stir=F)
addOrg(arena, bac, amount=1)
#addSubs(arena, smax=20, difunc="cpp", difspeed=1)
addSubs(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)", "EX_h2o(e)", "EX_h2o(e)", "EX_h(e)", "EX_nh4(e)"), difunc="cpp", difspeed=1)
res <- simEnv(arena, time=50)

plotCurves2(res)
