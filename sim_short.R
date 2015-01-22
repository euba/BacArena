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
load("data/ecore_model.R")
ecore = model
ecore = changeBounds(ecore, c('EX_ac(e)','EX_akg(e)','EX_etoh(e)','EX_for(e)',
                      'EX_fum(e)','EX_lac_D(e)','EX_pyr(e)','EX_succ(e)'), -20)

bace1 = Bac(model=ecore, deathrate=0.2, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=2, type="ecore1")
bace2 = Bac(model=ecore, deathrate=0.2, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=2, type="ecore2")
arena = Arena(n=100, m=100, tstep=1)
addOrg(arena, bace1, amount=10)
addOrg(arena, bace2, amount=10)
addSubs(arena, smax=20)

print(system.time(evalsim <- simulate(arena, time=30)))
format(object.size(evalsim), units='Mb')

evalArena(evalsim, plot_items=c('population'), phencol=T, retdata=F)
plotCurves(evalsim, remove=T, retdata = F)


library(animation)
saveVideo(evalArena(evalsim, plot_items=c('population'), phencol=F, retdata=F), video.name = "Pop4.mp4", other.opts = "-b 300k")  # higher bitrate, better quality
