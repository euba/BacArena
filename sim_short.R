# load libraries and other R files to have everything in place
setwd('C:/Users/eugen.bauer/Documents/GitHub/BacArena')
setwd("~/BacArena")
library(Rcpp)
library(sybil)
library(glpkAPI) 
#library(microbenchmark)
#library(ggplot2)
library(compiler) # byte code 
SYBIL_SETTINGS("SOLVER", "glpkAPI")
#source(file="cpp_source.R")
#source(file="R/class_baggage.R")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Bac.R")
source(file="R/Organism.R")
source(file="R/Eval.R")
#source(file="data/Human.R")
Rcpp::sourceCpp("cpp/diff.cpp")
#Rcpp::sourceCpp("cpp/addBac.cpp")

set.seed(5000)
load("models/ecore_model.R")
ecore = model
ecore1 = changeBounds(ecore, names(bace1@lbnd[bace1@medium][which(bace1@lbnd[bace1@medium]==0)]), -1000)
ecore2 = changeBounds(ecore, c('EX_o2(e)'), 0, 0)

bace1 = Bac(model=ecore, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=2, type="ecore1", lyse=T)
bace2 = Bac(model=ecore2, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
           speed=2, type="ecore2")
arena = Arena(n=1000, m=1000, tstep=0.5)
addOrg(arena, bace1, amount=10)
addOrg(arena, bace2, amount=10)
addSubs(arena, smax=100)

print(system.time(evalsim <- simulate(arena, time=20)))
format(object.size(evalsim), units='Mb')

evalArena(evalsim, plot_items=c('population','EX_o2(e)'), phencol=T, retdata=F)
plotCurves(evalsim, remove=T, retdata=F)
minePheno(evalsim)


library(animation)
saveVideo(evalArena(evalsim, plot_items=c('population'), phencol=T, retdata=F), video.name="Pop4_phen.mp4", other.opts="-b 300k")  # higher bitrate, better quality
