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

simlist <- list()

pop = Population(specs, specn=rep(10, length(specs)), n=100, m=100)
addSub(pop, "EX_glc(e)", 20)
addSub(pop, "EX_h2o(e)", 20)
addSub(pop, "EX_o2(e)", 20)
addSub(pop, "EX_pi(e)", 20)

for(i in 1:100){
  simlist[[i]] <- pop
  print(system.time(for(j in seq_along(pop@media)){
    #diffuseNaiveR(pop@media[[j]])
    diffuseNaiveCpp(pop@media[[j]]@diffmat, donut=FALSE)
  }))
  j = 0
  print(system.time(while(j+1 <= length(pop@orglist)){
    j<-j+1
    move(pop@orglist[[j]],pop)
    medcon = getmed(pop,pop@orglist[[j]]@x,pop@orglist[[j]]@y)
    constrain(pop@orglist[[j]], names(medcon), lb=-medcon)
    optimizeLP(pop@orglist[[j]])
    pop@media = consume(pop@orglist[[j]],pop@media)
    growth(pop@orglist[[j]], pop, j,lifecosts=0.6)
  }))
  if(length(pop@orglist)==0){
    print("All bacs dead!")
    break
  } 
  cat("iter:", i, "bacs:",length(pop@orglist),"\n\n")
}
