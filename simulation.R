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

pop = Population(specs, specn=rep(50, length(specs)), n=40, m=40)

for(i in 1:50){
  moveRand(pop)
  #lapply(pop@media, diffuseNaive)
  for(j in seq_along(pop@media)){
    diffuseNaive(pop@media[[j]])
  }
  bacs <- pop@orglist
  media <- pop@media
  for(i in seq_along(bacs)){
    bac <- bacs[[i]]
    for(k in seq_along(media)){
      constrain(bac, media[[k]]@name, lb=-media[[k]]@diffmat[bac@x,bac@y])
      consume(bac)
    }
    optimizeLP(bac)
    consume()
  }
  
  dmat <- pop@media[["EX_glc(e)"]]@diffmat
  image(dmat)
  image(pop2mat(pop))
  #print(i)
  #print(sum(unlist(lapply(pop@orglist, function(x){print(x@growth)}))))
  #print(pop@orgn)
}




