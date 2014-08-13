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

pop = Population(specs, specn=rep(1, length(specs)), n=2, m=2)

for(i in 1:2){
  moveRand(pop)
  #lapply(pop@media, diffuseNaive)
  for(j in seq_along(pop@media)){
    diffuseNaive(pop@media[[j]])
  }
  bacs <- pop@orglist
  media <- pop@media
  for(j in seq_along(bacs)){
    bac <- bacs[[j]]
    for(k in seq_along(media)){
      constrain(bac, media[[k]]@name, lb=-media[[k]]@diffmat[bac@x,bac@y])
      consume(bac, pop@media[[k]])
    }
    print(bac@)
    optimizeLP(bac)
  }
  
  #dmat <- pop@media[["EX_glc(e)"]]@diffmat
  #image(dmat)
  image(pop2mat(pop))
  print(i)
  print(pop@media[["EX_glc(e)"]]@diffmat)
  #print(sum(unlist(lapply(pop@orglist, function(x){print(x@growth)}))))
  #print(pop@orgn)
}




