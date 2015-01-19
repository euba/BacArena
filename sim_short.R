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
Rcpp::sourceCpp("cpp/diff.cpp")
#Rcpp::sourceCpp("cpp/addBac.cpp")

set.seed(5000)
load("data/ecore_model.R")
ecore = model
ecore = changeBounds(ecore, c('EX_ac(e)','EX_akg(e)','EX_etoh(e)','EX_for(e)',
                      'EX_fum(e)','EX_lac_D(e)','EX_pyr(e)','EX_succ(e)'), -20)

bace = Bac(model=ecore, deathrate=0.2, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=2, budge=F, lyse=F, chem="EX_glc(e)")
arena = Arena(n=100, m=100, tstep=1, stir=T)
#addOrg(arena, bace, amount=1, x=50, y=50)
addOrg(arena, bace, amount=1)
#print(system.time(addOrg(arena, bace, amount=50000)))
addSubs(arena, smax=1000, mediac=c("EX_glc(e)","EX_pi(e)"))

simlist <- simulate(arena, time=30)

plot(arena@orgdat[,c('x','y')])
image(as.matrix(arena@occmat))
stirEnv(arena,sublb = 0)


changeSub(arena, "EX_o2(e)", 1000)
format(object.size(arena), units='b')

#arena@media[["EX_o2(e)"]]@diffmat <- Matrix(0, nrow=1000, ncol=1000)
#arena@media[["EX_glc(e)"]]@diffmat <- Matrix(0, nrow=100, ncol=100)
#arena@media[["EX_glc(e)"]]@diffmat[50,50] <- 5000
simlist <- simulate(arena, time=1)
addEval(simlist, arena, replace=T)
for(i in 1:40){
  print(i)
  simlist <- simulate(simlist, time=1)
  arena <- getArena(simlist)
  changeSub(arena, "EX_o2(e)", 1000)
  addEval(simlist, arena, replace=T)
}

#print(system.time(simlist <- simulate(arena, time=30)))
format(object.size(simlist), units='Mb')

evalArena(simlist, plot_items=c('population'), phencol=T, retdata=F)
plotCurves(simlist, remove=T, retdata = F)
