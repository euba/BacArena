setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
SYBIL_SETTINGS("SOLVER", "glpkAPI")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")


data(Ec_core)
bac = Bac(model=Ec_core, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential", speed=1, type="ecoli")

bac@lbnd[["EX_o2(e)"]] <- -18.5
setKinetics(bac, "EX_fru(e)", Km=39, vmax=5.58)
#setKinetics(bac, "EX_succ(e)", Km=15.5, vmax=3.18)


arena = Arena(n=10, m=10, stir=F)
addOrg(arena, bac, amount=5)
addSubs(arena,20,c("EX_fru(e)", "EX_o2(e)","EX_pi(e)", "EX_h2o(e)", "EX_h(e)", "EX_nh4(e)"), difunc="cpp", difspeed=1)

sim <- simEnv(arena, time=100)

