setwd("~/uni/bacarena")
setwd('/Users/euba/GitRep/BacArena')
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(ReacTran)
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")

SYBIL_SETTINGS("SOLVER","sybilGUROBI") #setting solver to GUROBI

#
# Simulation
#
data(Ec_core)
Ec_core = changeBounds(Ec_core,react_id(findExchReact(Ec_core)),lb=-1000)
Ec_core = changeBounds(Ec_core,'EX_glc(e)',lb=-10)
bac = Bac(model=Ec_core, growtype="exponential", cellarea=4.42, lyse=F)
setKinetics(bac, exchangeR="EX_glc(e)", Km=0.01, vmax=7.56)
arena = Arena(n=200, m=200, stir=F, seed=8904, Lx=0.05, Ly=0.05, tstep=0.2)
#addOrg(arena, bac, amount=11, x=c((arena@n/2-5):(arena@n/2+5)), y=c((arena@m/2-5):(arena@m/2+5)))
addOrg(arena, bac, amount=1, x=arena@n/2, y=arena@m/2,growth = 0.9)
#addOrg(arena, bac, amount=arena@n, x=1:arena@n, y=arena@m/2)
addSubs(arena, smax=0.05, difspeed=6.7e-6, unit='mM') #0.02412#0.072
#addSubs(arena, smax=10, mediac=c("EX_o2(e)","EX_h(e)","EX_co2(e)","EX_o2(e)","EX_pi(e)"), difunc="pde", difspeed=rep(0.072,5))
#createGradient(arena,smax=20,mediac="EX_o2(e)",position='left',steep=0.5)
#createGradient(arena,smax=20,mediac=arena@mediac,position='left',steep=0.5)
sim <- simEnv(arena, time=10, lrw=26937744)

SYBIL_SETTINGS("TOLERANCE",1E-6) #set tolerance of FBA value higher
SYBIL_SETTINGS("MAXIMUM",1) #set tolerance of FBA value higher
mod = changeBounds(Ec_core,react_id(findExchReact(Ec_core)),lb=-10000)
optimizeProb(mod)

optimizeProb(Ec_core)
#1x1 -> 1
#10x10 -> 1870
#20x20 -> 7420
#30x30 -> 16670
#40x40 -> 29620
#50x50 -> 46270
#60x60 -> 66620
#70x70 -> 90670
#80x80 -> 118420
#90x90 -> 149870
#100x100 -> 185020
#200x200 -> 740020
#300x300 -> 1665020
#400x400 -> 2960020
#500x500 -> 4625020
#600x600 -> 6660020
#700x700 -> 9065020
#800x800 -> 11840020
#900x900 -> 14985020
#1000x1000 -> 18500020

are = c(10*10,20*20,30*30,40*40,50*50,60*60,70*70,80*80,90*90,100*100)#,200*200,300*300,400*400,500*500,600*600,700*700,800*800,900*900,1000*1000)
ns = c(10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000)
lrw = c(1870,7420,16670,29620,46270,66620,90670,118420,149870,185020)#,740020,1665020,2960020,4625020,6660020,9065020,11840020,14985020,18500020)
lm(lrw ~ are)$coefficients

are*18.5 + 20

plot(ns,are,type="b")
plot(ns,lrw,type="b")
plot(are,lrw,type="b")
#
# Evaluation
#
plotCurves(sim)
evalArena(sim, plot_items = c("Population", "EX_glc(e)", "EX_o2(e)", "EX_ac(e)"),phencol=T)#, time=10)
evalArena(sim, plot_items = c("Population", "EX_glc(e)", "EX_ac(e)", "EX_o2(e)"),phencol=T, time=120) 
evalArena(sim, plot_items = "Population")
plotCurves2(sim)

corr <- getCorrM(sim, reactions=False)
checkCorr(sim,tocheck = "o2")
plotTotFlux(sim, legendpos = "topleft")


library(animation)
oopts = ani.options(ffmpeg = "/Users/euba/bin/ffmpeg/ffmpeg")
saveVideo({
  ani.options(interval = 0.5)
  videoPhen(evalsim)
},video.name = "pao_sim.mp4", other.opts = "-b 600k")
