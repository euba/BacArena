##Loading the required packages## 
library(ReacTran)
library(sybil)
library(sybilSBML)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(igraph)
setwd('P:/GitRep/BacArena')

source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")

load("P:/GitRep/Vis_Simulation_Paper/010416/All_orgs_renamed_mk_ready.RData")
orglist = orglist[c(1,2,3,5,7,8,9,11,13,14)]

#########################################################################
sim_cplex = list()
library(cplexAPI)
SYBIL_SETTINGS("SOLVER","cplexAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_cplex[[i]] = simEnv(arena, time=10, sec_obj="none")
}

sim_gurobi= list()
library(sybilGUROBI)
SYBIL_SETTINGS("SOLVER","sybilGUROBI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_gurobi[[i]] = simEnv(arena, time=10, sec_obj="none")
}

sim_glpk = list()
library(glpkAPI)
SYBIL_SETTINGS("SOLVER","glpkAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_glpk[[i]] = simEnv(arena, time=10, sec_obj="none")
}
#################################################################
sim_cplex_mtf = list()
library(cplexAPI)
SYBIL_SETTINGS("SOLVER","cplexAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_cplex_mtf[[i]] = simEnv(arena, time=10, sec_obj="mtf")
}

sim_gurobi_mtf= list()
library(sybilGUROBI)
SYBIL_SETTINGS("SOLVER","sybilGUROBI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_gurobi_mtf[[i]] = simEnv(arena, time=10, sec_obj="mtf")
}

sim_glpk_mtf = list()
library(glpkAPI)
SYBIL_SETTINGS("SOLVER","glpkAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_glpk_mtf[[i]] = simEnv(arena, time=10, sec_obj="mtf")
}

###########################################################
#########################################################################
###########################################################
sim_cplex_opt_rxn = list()
library(cplexAPI)
SYBIL_SETTINGS("SOLVER","cplexAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_cplex_opt_rxn[[i]] = simEnv(arena, time=10, sec_obj="opt_rxn")
}

sim_gurobi_opt_rxn= list()
library(sybilGUROBI)
SYBIL_SETTINGS("SOLVER","sybilGUROBI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_gurobi_opt_rxn[[i]] = simEnv(arena, time=10, sec_obj="opt_rxn")
}

sim_glpk_opt_rxn = list()
library(glpkAPI)
SYBIL_SETTINGS("SOLVER","glpkAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_glpk_opt_rxn[[i]] = simEnv(arena, time=10, sec_obj="opt_rxn")
}
#################################################################
sim_cplex_opt_ex = list()
library(cplexAPI)
SYBIL_SETTINGS("SOLVER","cplexAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_cplex_opt_ex[[i]] = simEnv(arena, time=10, sec_obj="opt_ex")
}

sim_gurobi_opt_ex= list()
library(sybilGUROBI)
SYBIL_SETTINGS("SOLVER","sybilGUROBI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_gurobi_opt_ex[[i]] = simEnv(arena, time=10, sec_obj="opt_ex")
}

sim_glpk_opt_ex = list()
library(glpkAPI)
SYBIL_SETTINGS("SOLVER","glpkAPI")
for(i in 1:10){
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1, seed=i)
  for(i in 1:length(orglist)){
    model = orglist[[i]]
    bac = Bac(model=model, growtype="exponential", speed=10)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=10, difspeed=6.7e-6, unit='mmol/cm2')
  sim_glpk_opt_ex[[i]] = simEnv(arena, time=10, sec_obj="opt_ex")
}
