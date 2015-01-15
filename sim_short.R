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
Rcpp::sourceCpp("cpp/diff.cpp")
#Rcpp::sourceCpp("cpp/addBac.cpp")

set.seed(5000)
load("data/ecore_model.R")
ecore = model

bace = Bac(model=ecore, deathrate=0.5, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=0, lyse=T, chem="EX_glc(e)")
arena = Arena(n=100, m=100, tstep=0.5)
addOrg(arena, bace, amount=1, x=50, y=50)
#addOrg(arena, bace, amount=500)
addSubs(arena, smax=30)
format(object.size(arena), units='b')
#arena@media[["EX_glc(e)"]]@diffmat <- Matrix(0, nrow=100, ncol=100)
#arena@media[["EX_glc(e)"]]@diffmat[50,50] <- 5000
print(system.time(simlist <- simulate(arena, time=30)))

for(i in 1:length(simlist)){
  print(i)
  popana <- as.matrix(simlist[[i]]@occmat)
  #dat <- simlist[[i]]@orgdat
  #for(j in 1:nrow(dat)){
  #  popana[dat[j,]$x,dat[j,]$y] <- dat[j,]$phenotype + 1
  #}
  par(mfrow = c(1, 2))
  image(popana, axes=F, main='Population')
  image(as.matrix(simlist[[i]]@media[["EX_glc(e)"]]@diffmat), axes=F,main='Glucose')
}
