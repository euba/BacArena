# load libraries and other R files to have everything in place
setwd('E:/GitRep/BacArena')
setwd("/Users/euba/GitRep/BacArena/")
setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(compiler) # byte code 
SYBIL_SETTINGS("SOLVER", "sybilGUROBI")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

set.seed(5000)
data(Ec_core)
ecore = Ec_core
bace = Bac(model=ecore, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=10, type="ec.ore1", lyse=F)
bace2 = Bac(model=ecore, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=10, type="ec.ore2", lyse=F)
arena = Arena(n=100, m=100, stir=F)
addOrg(arena, bace, amount=100, x=1:100, y=rep(1,100))
addOrg(arena, bace2, amount=1, x=50, y=50)
addSubs(arena, smax=20, difunc="cpp", difspeed=1)

setwd('benchmark')
Rprof(filename = "benchmark_2.out")
evalsim <- simEnv(arena, time=10)
Rprof(NULL)

bench = summaryRprof("benchmark_2.out")
rtime = c((bench$by.total["\"simEnv\"",1]-bench$by.total["\"simBac\"",1]),(bench$by.total["\"simBac\"",1]-bench$by.total["\"optimizeLP\"",1]),bench$by.total["\"optimizeLP\"",1])
names(rtime) = c("Background","Individual","FBA")
barplot(rtime)
bench$by.total[1,1]


bench$by.total





set.seed(5000)
data(Ec_core)
ecore = Ec_core
ecore1 = changeBounds(ecore, names(bace1@lbnd[bace1@medium][which(bace1@lbnd[bace1@medium]==0)]), -1000)
ecore2 = changeBounds(ecore, c('EX_o2(e)'), 0, 0)

bace1 = Bac(model=ecore, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=1, type="ecore1", lyse=T)
bace2 = Bac(model=ecore2, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=0, type="ecore2")
arena = Arena(n=30, m=30, stir=F)
addOrg(arena, bace1, amount=10)
addOrg(arena, bace2, amount=10,x=1:10,y=1:10)
addSubs(arena, smax=20, difunc="cpp", difspeed=1)

print(system.time(evalsim <- simEnv(arena, time=10)))
format(object.size(evalsim), units='Mb')
evalArena(evalsim)

evalArena(evalsim, plot_items=c('population','EX_o2(e)'), phencol=F, retdata=F)
plotCurves(evalsim, remove=T, retdata=F)
minePheno(evalsim, time=8)

selPheno(evalsim,time=1,type='ecore1')

library(animation)
saveVideo(evalArena(evalsim, plot_items=c('population'), phencol=T, retdata=F), video.name="Pop4_phen.mp4", other.opts="-b 300k")  # higher bitrate, better quality


setwd("P:/BACARENA/benchmark")
set.seed(5000)
data(Ec_core)
ecore = Ec_core
pao = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/modelPOA.xml")

bace = Bac(model=ecore, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
            speed=1, type="ecore1", lyse=F)
arena = Arena(n=100, m=100, stir=F)
addOrg(arena, bace, amount=100, x=1:100, y=rep(1,100))
addSubs(arena, smax=20, difunc="cpp", difspeed=1)

Rprof()
evalsim <- simEnv(arena, time=10)
Rprof(NULL)

bench = summaryRprof("Rprof.out")

rtime = c((bench$by.total["\"simEnv\"",1]-bench$by.total["\"simBac\"",1]),(bench$by.total["\"simBac\"",1]-bench$by.total["\"optimizeLP\"",1]),bench$by.total["\"optimizeLP\"",1])
names(rtime) = c("Background","Individual","FBA")
barplot(rtime)

