setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(ReacTran)
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
arena = Arena(n=90, m=90, stir=F, seed=8904)
#addOrg(arena, bac, amount=11, x=c((arena@n/2-5):(arena@n/2+5)), y=c((arena@m/2-5):(arena@m/2+5)))
addOrg(arena, bac, amount=1, x=round(n(arena)/2), y=round(m(arena)/2))
addSubs(arena, smax=10,difspeed=0.1)
#addSubs(arena, smax=c(20,50), mediac=c("EX_glc(e)","EX_o2(e)"), difunc="pde", difspeed=c(0.2,5))
#createGradient(arena,smax=20,mediac="EX_o2(e)",position='left',steep=0.5)
#createGradient(arena,smax=20,mediac=arena@mediac,position='left',steep=0.5)
sim <- simEnv(arena, time=30)

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
evalArena(sim, plot_items = c("Population", "EX_glc(e)", "EX_o2(e)", "EX_for(e)"),phencol=T)#, time=10)
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