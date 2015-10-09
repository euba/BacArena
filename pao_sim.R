setwd("P:/BACARENA/Comparison/MatNet/P_aeruginosa/")

##Loading the required packages## 
library(sybil)
#library(sybilSBML)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(compiler)
setwd("/Users/euba/GitRep/BacArena/")
#setwd('P:/GitRep/BacArena')
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

SYBIL_SETTINGS("SOLVER","sybilGUROBI") #setting solver to GUROBI

library(sybilSBML)
#model = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/modelPOA.xml")
load('/Users/euba/GitRep/BacArena/poa_model.RData')

set.seed(5000)
# paramters from MatNet:
# if cellMass >= 2
# cellDeathThreshold = 0
model@lowbnd[grep("EX",model@react_id)]
modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-20)


bace1 = Bac(model=modelP, deathrate=0.05, duplirate=1, growthlimit=0.05, growtype="exponential",
            speed=3, type="PAO", lyse=F)
arena = Arena(n=100, m=100, stir=F, tstep=0.5)
addOrg(arena, bace1, amount=1, x=50, y=50)
addSubs(arena, smax=50, difunc="cpp", difspeed=1,
        mediac=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)])

print(system.time(evalsim <- simEnv(arena, time=50)))
# 
# library(animation)
# oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
# saveVideo({
#   ani.options(interval = 0.5)
#   evalArena(evalsim, phencol=T, plot_items=c('population','EX_EC0007'))
# },video.name = "PAO_pop_phen_bio2.avi", other.opts = "-b 600k")
# 
# library(animation)
# oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
# saveVideo({
#   ani.options(interval = 0.5)
#   evalArena(evalsim, phencol=T)
# },video.name = "PAO_pop5.avi", other.opts = "-b 600k")
