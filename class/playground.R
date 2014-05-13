# this file is a playground to test how the current classes work in action
# the main goal of this file is to construct a basic framework for BacArena, which can then be merged with diffbac
# it is actually a little bit like diffbac.R, but for the current oop version of BacArena

setwd("~/BacArena")

optimizeProb(bac1@model, solverParm = list(PRESOLVE = "GLP_ON"))


# load libraries and other R files to have everything in place
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

#testing constructor
bac1 = Bac(x=1, y=1, model=mod, growth=1, n=1, m=1, type="test")
#bac1 = Bac(x=1, y=1, model=mod, growth=1, fobj="ATPM")

findUpt(bac1, flag=F)

constrain(bac1, findUpt(bac1), lb=0, ub=1000)

optimizeLP(bac1)

optimizeProb(bac1@model, solverParm = list(PRESOLVE = GLP_ON))
optimizeProb(bac1@model, solverParm = list(warmUpGLPK = GLP_ON))
bac1@lpobj
optimizer(bac1@model)

checkOptSol(bac1@lpobj)

constrain(bac1, "EX_o2(e)", lb=-1000, ub=0)
optimizeLP(bac1)
bac1@lpobj

bac1@model = changeFobj(bac1@model, "ATPM")
org1 <- Organism(x=1, y=2, model=mod, n=1, m=1)
#org1@model = org1@changeObj(org1@model, "ATPM")

#bac1 = Bac(org1, growth=1) #this does not work, but would be cool if...
#testing constructor
org1 <- Organism(x=1, y=2, model=mod, n=1, m=1)
#checking the organisms functions
org1@lpobj
constrain(org1, "EX_o2(e)", lb=0, ub=0)
optimizeLP(org1)
org1@lpobj
changeFobj(org1, "ATPM")
org1@lpobj

#testing linear growth of bacteria
growLin(bac1)
bac1@growth
for(i in 1:100){
  growLin(bac1)
  print(bac1@growth)
}

#testing constructor
n=100
m=100
smax=50
diffmat = matrix(smax, nrow=n, ncol=m)
diffmat[(n/2-n/4):(n/2+n/4), (m/2-m/4):(m/2+m/4)] = 0
sub1 <- Substance(n=n, mm, smax=smax, name="test", diffmat=diffmat)
jpeg(paste("plot", formatC(1, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
image(sub1@diffmat, axes=FALSE, col=rev(heat.colors(200)))
dev.off()
diffuseNaive(sub1)
diffuseNaiveR(sub1)
substrate=list()
dmat = sub1@diffmat
substrate[[1]] <- dmat
for(i in 2:100){
  diffuseNaive(sub1)
  #dmat = sub1@diffmat
  #jpeg(paste("plot", formatC(i, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
  #image(dmat, axes=FALSE, col=rev(heat.colors(200)))
  #jpeg(paste("plot", formatC(i, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
  image(sub1@diffmat, axes=FALSE, col=rev(heat.colors(200)))
  #dev.off()
  print(i)
}
diffusion(list(dmat))

substrate <- list()


specs3 = list(specs[[1]], specs[[2]], specs[[3]])
specs=list()
specs[[1]]=mod

pop = Population(specs3, specn=rep(5, length(specs3)), n=20, m=20)
repliDie(pop)
#jpeg(paste("plot", formatC(1, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
image(pop2mat(pop))
#dev.off()
for(i in 2:200){
  moveRand(pop)
  #jpeg(paste("plot", formatC(i, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
  image(pop2mat(pop))
  #dev.off()
  print(i)
}


pop = Population(specs, specn=10, n=10, m=10)
for(i in 1:50){
  moveRand(pop)
  repliDie(pop)
  for(j in seq_along(pop@media)){
    diffuseNaive(pop@media[[j]])
  }
  dmat <- pop@media[["EX_glc(e)"]]@diffmat
  image(dmat)
  image(pop2mat(pop))
  print(i)
  print(sum(unlist(lapply(pop@orglist, function(x){print(x@growth)}))))
  print(pop@orgn)
}
