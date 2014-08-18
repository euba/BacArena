# load libraries and other R files to have everything in place
setwd("~/BacArena")
library(snow) # parallel computing
library(Rcpp)
library(inline)
library(sybil)
#SYBIL_SETTINGS("SOLVER", "glpkAPI")
SYBIL_SETTINGS("SOLVER", "clpAPI")
#load class definitions
source(file="cpp_source.R")
source(file="class/class_baggage.R")
source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Bac.R")
source(file="class/Organism.R")
source(file="class/Population.R")

Rcpp::sourceCpp("diff.cpp")

#load ecoli core model to play around
load("data/ecore_model.R")
mod <- model
specs=list()
specs[[1]]=mod

pop = Population(specs, specn=rep(20, length(specs)), n=10, m=5)
for(i in 1:100){
  #print(system.time(moveRand(pop)))
  #lapply(pop@media, diffuseNaive)
  print(system.time(for(j in seq_along(pop@media)){
    #diffuseNaiveR(pop@media[[j]])
    diffuseNaiveCpp(pop@media[[j]]@diffmat)
  }))
  diffuseNaive(pop@media)
  #print(system.time(lapply(pop@media, diffuseNaiveR)))
  max <- seq_along(pop@orglist)
  print(system.time(for(j in max){
    move(pop@orglist[[j]],pop)
    medcon = getmed(pop,pop@orglist[[j]]@x,pop@orglist[[j]]@y)
    constrain(pop@orglist[[j]], names(medcon), lb=-medcon)
    optimizeLP(pop@orglist[[j]])
    pop@media = consume(pop@orglist[[j]],pop@media)
    #print(pop@orglist[[j]]@fbasol$fluxes[which(react_id(pop@orglist[[j]]@model) == "EX_glc(e)")])
    #print(pop@orglist[[j]]@growth)
    #pop@orglist=growth(pop@orglist[[j]], pop, j)
    growth(pop@orglist[[j]], pop, j)
  }))
  
  #dmat <- pop@media[["EX_glc(e)"]]@diffmat
  #image(dmat)
  #image(pop2mat(pop))
  print(paste("iter:", i, "bacs:",length(pop@orglist)))
  #print(pop@media[["EX_glc(e)"]]@diffmat)
  #print(sum(unlist(lapply(pop@orglist, function(x){print(x@growth)}))))
  #print(pop@orgn)
}

o1 = pop@orglist[[1]]
optimizeLP(o1)

