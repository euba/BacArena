# load libraries and other R files to have everything in place
setwd("~/BacArena")
library(snow) # parallel computing
library(Rcpp)
library(inline)
library(sybil)
SYBIL_SETTINGS("SOLVER", "clpAPI")
#load class definitions
source(file="cpp_source.R")
source(file="class/class_baggage.R")
source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Bac.R")
source(file="class/Organism.R")
source(file="class/Population.R")

#load ecoli core model to play around
load("data/ecore_model.R")
mod <- model
specs=list()
specs[[1]]=mod

pop = Population(specs, specn=rep(1000, length(specs)), n=100, m=100)
for(i in 1:20){
  print(system.time(moveRand(pop)))
  #lapply(pop@media, diffuseNaive)
  print(system.time(for(j in seq_along(pop@media)){
    diffuseNaive(pop@media[[j]])
  }))
  #print(system.time(lapply(pop@media, diffuseNaiveR)))
  media <- pop@media
  print(system.time(for(j in seq_along(pop@orglist)){
    medcon = getmed(pop,pop@orglist[[j]]@x,pop@orglist[[j]]@y)
    constrain(pop@orglist[[j]], names(medcon), lb=-medcon)
    #for(k in seq_along(media)){
    #  pop@media[[k]] <- consume(pop@orglist[[j]], media[[k]])
    #}
    optimizeLP(pop@orglist[[j]])
    #print(pop@orglist[[j]]@fbasol$fluxes[which(react_id(pop@orglist[[j]]@model) == "EX_glc(e)")])
    #print(pop@orglist[[j]]@growth)
    growth(pop@orglist[[j]], j, pop)
  }))
  
  #dmat <- pop@media[["EX_glc(e)"]]@diffmat
  #image(dmat)
  #image(pop2mat(pop))
  print(paste("iter:", i, "bacs:",length(pop@orglist)))
  #print(pop@media[["EX_glc(e)"]]@diffmat)
  #print(sum(unlist(lapply(pop@orglist, function(x){print(x@growth)}))))
  #print(pop@orgn)
}




