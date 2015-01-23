# load libraries and other R files to have everything in place
setwd("~/BacArena")
library(Rcpp)
library(sybil)
library(glpkAPI) 
#library(microbenchmark)
#library(ggplot2)
library(compiler) # byte code 
SYBIL_SETTINGS("SOLVER", "glpkAPI")
source(file="cpp_source.R")
source(file="class/class_baggage.R")
source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Bac.R")
source(file="class/Organism.R")
source(file="class/Eval.R")
source(file="class/Human.R")
Rcpp::sourceCpp("cpp/diff.cpp")
#Rcpp::sourceCpp("cpp/addBac.cpp")

set.seed(5000)
load("data/Bcoli_model.R")
ecore = model
ecore1 = changeBounds(ecore, names(bace1@lbnd[bace1@medium][which(bace1@lbnd[bace1@medium]==0)]), -1000)
ecore2 = changeBounds(ecore, c('EX_o2(e)'), 0, 0)

bace1 = Bac(model=ecore, deathrate=0.1, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=2, type="ecore1", lyse=T)
bace2 = Bac(model=ecore2, deathrate=0.1, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=2, type="ecore2")
arena = Arena(n=100, m=100, tstep=1)
addOrg(arena, bace1, amount=10)
addOrg(arena, bace2, amount=10)
addSubs(arena, smax=2000)

print(system.time(evalsim <- simulate(arena, time=20)))
format(object.size(evalsim), units='Mb')

evalArena(evalsim, plot_items=c('population'), phencol=F, retdata=F)
plotCurves(evalsim, remove=F, retdata = F)
minePheno(evalsim)


library(animation)
saveVideo(evalArena(evalsim, plot_items=c('population'), phencol=T, retdata=F), video.name="Pop4_phen.mp4", other.opts="-b 300k")  # higher bitrate, better quality
