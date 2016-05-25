##Loading the required packages## 
library(ReacTran)
library(sybil)
library(sybilSBML)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(igraph)
library(parallel)
setwd('P:/GitRep/BacArena')

source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")

load("P:/GitRep/Vis_Simulation_Paper/010416/All_orgs_renamed_mk_ready.RData")
orglist = orglist[c(1,2,3,5,7,8,9,11,13,14)]

orglist = orglist[c(1,6)]

#########################################################################
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="none")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_gurobi <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","sybilGUROBI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="none")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_glpk <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="none")
}) ))
stopCluster(cl)
#################################################################
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex_mtf <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="mtf")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_gurobi_mtf <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","sybilGUROBI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="mtf")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_glpk_mtf <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="mtf")
}) ))
stopCluster(cl)
###########################################################
#########################################################################
###########################################################
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex_opt_rxn <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="opt_rxn")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_gurobi_opt_rxn <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","sybilGUROBI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="opt_rxn")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_glpk_opt_rxn <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="opt_rxn")
}) ))
stopCluster(cl)
#################################################################
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex_opt_ex <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="opt_ex")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_gurobi_opt_ex <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","sybilGUROBI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="opt_ex")
}) ))
stopCluster(cl)

cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_glpk_opt_ex <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="opt_ex")
}) ))
stopCluster(cl)


save(sim_cplex, file="P:/BACARENA/AOS/sim_cplex.RData")
save(sim_gurobi, file="P:/BACARENA/AOS/sim_gurobi.RData")
save(sim_glpk, file="P:/BACARENA/AOS/sim_glpk.RData")
save(sim_cplex_mtf, file="P:/BACARENA/AOS/sim_cplex_mtf.RData")
save(sim_gurobi_mtf, file="P:/BACARENA/AOS/sim_gurobi_mtf.RData")
save(sim_glpk_mtf, file="P:/BACARENA/AOS/sim_glpk_mtf.RData")

save(sim_cplex_opt_rxn, file="P:/BACARENA/AOS/sim_cplex_opt_rxn.RData")
save(sim_gurobi_opt_rxn, file="P:/BACARENA/AOS/sim_gurobi_opt_rxn.RData")
save(sim_glpk_opt_rxn, file="P:/BACARENA/AOS/sim_glpk_opt_rxn.RData")
save(sim_cplex_opt_ex, file="P:/BACARENA/AOS/sim_cplex_opt_ex.RData")
save(sim_gurobi_opt_ex, file="P:/BACARENA/AOS/sim_gurobi_opt_ex.RData")
save(sim_glpk_opt_ex, file="P:/BACARENA/AOS/sim_glpk_opt_ex.RData")

#########################################################################
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex_bar <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  sybil::SYBIL_SETTINGS("METHOD","baropt")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="none")
}) ))
stopCluster(cl)
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex_prim <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  sybil::SYBIL_SETTINGS("METHOD","primopt")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="none")
}) ))
stopCluster(cl)
cl <- makeCluster(5, type="PSOCK",outfile="cluster_sihumi.log") # PSOCK works with win/mac/lin
clusterExport(cl, "orglist")
print(system.time(sim_cplex_hbar <- parLapply(cl, 1:10, function(i){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  setwd('P:/GitRep/BacArena')
  source(file="R/Arena.R")
  source(file="R/Stuff.R")
  source(file="R/Substance.R")
  source(file="R/Organism.R")
  Rcpp::sourceCpp("src/diff.cpp")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  sybil::SYBIL_SETTINGS("METHOD","hybbaropt")
  arena <- Arena(50, 50, stir=F, Lx=0.01, Ly=0.01, tstep=1)
  for(j in 1:length(orglist)){
    model = orglist[[j]]
    bac = Bac(model=model, growtype="exponential", speed=5)
    arena = addOrg(arena, bac, amount=5)
  }
  arena=addSubs(arena, smax=5*10^-5, difspeed=6.7e-6, unit='mmol/cm2')
  sim = simEnv(arena, time=10, sec_obj="none")
}) ))
stopCluster(cl)

#######################################################################################
################################################### Analysis
#######################################################################################







