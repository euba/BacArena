setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(compiler) # byte code 
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")

#
# Simulation
#
data(Ec_core)
bac = Bac(model=Ec_core, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=1, type="ecore", lyse=T)
arena = Arena(n=50, m=50, stir=F, seed=8904)
addOrg(arena, bac, amount=11, x=c((arena@n/2-5):(arena@n/2+5)), y=c((arena@m/2-5):(arena@m/2+5)))
addSubs(arena, smax=20, difunc="pde", difspeed=1)
sim <- simEnv(arena, time=1000)

#
# Evaluation
#
plotCurves(sim)
evalArena(sim, plot_items = c("Population", "EX_glc(e)", "EX_o2(e)", "EX_for(e)"), time=10)
evalArena(sim, plot_items = c("Population", "EX_glc(e)", "EX_ac(e)", "EX_o2(e)"), time=10) 
evalArena(sim, plot_items = "Population")
plotCurves2(sim)

corr <- getCorrM(sim, reactions=False)
checkCorr(sim,tocheck = "o2")
plotTotFlux(sim, legendpos = "topleft")


library(animation)
#oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
saveVideo({
  ani.options(interval = 0.5)
  evalArena(sim, phencol=T)
},video.name = "ecoli.avi", other.opts = "-b 600k")