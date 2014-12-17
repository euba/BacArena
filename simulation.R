# load libraries and other R files to have everything in place
#setwd("~/BacArena")
library(snow) # parallel computing
library(Rcpp)
library(inline)
library(sybil)
#SYBIL_SETTINGS("SOLVER", "glpkAPI")
SYBIL_SETTINGS("SOLVER", "clpAPI")
source(file="cpp_source.R")
source(file="class/class_baggage.R")
#source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Bac.R")
source(file="class/Organism.R")
#source(file="class/Population.R")
Rcpp::sourceCpp("diff.cpp")

setwd("C:/Users/eugen.bauer/Documents/GitHub/BacArena/") #have to update orgn after deletion of bacteria?
source(file="class/Arena.R")
#load("/home/eugen/specsm.RData")

load("data/Bcoli_model.R")
ecoli = model
load("data/ecore_model.R")
ecore = model
load("data/clos_model.R")
clos = model
load("data/barkeri_model.R")
meth = model

bacm = Bac(model=meth, deathrate=0.3, duplirate=1.5, growthlimit=0.05, growtype="exponential")
bacc = Bac(model=clos, deathrate=0.3, duplirate=1.5, growthlimit=0.05, growtype="exponential", ex="ex_")
bace = Bac(model=ecore, deathrate=1, duplirate=1.5, growthlimit=0.05, growtype="exponential")
#bace = Bac(model=ecoli, deathrate=1, duplirate=1.5, growthlimit=0.05, growtype="exponential")

arena = Arena(n=200, m=200)
#addOrg(arena, bacm, amount=1, x=25, y=25)
#addOrg(arena, bacc, amount=1, x=24, y=25)
addOrg(arena, bace, amount=10)
addSubs(arena, smax=30)
format(object.size(arena), units='Mb')

print(system.time(simlist <- simulate(arena, time=20)))

format(object.size(simlist), units='Mb')

for(i in seq_along(simlist)){
  #popana <- as.matrix(simlist[[i]]@occmat)
  popana <- as.matrix(simlist[[i]]@media[['EX_for(e)']]@diffmat)
  image(popana)
}









object.size(arena)
#object.size(arena@occmat)
#object.size(as.matrix(arena@occmat))
addSubs(arena, smax=20)#, mediac=names(arena@phenotypes[[1]][[1]]))
object.size(arena)
format(object.size(arena), units='Mb')

print(system.time(simlist <- simulate(arena, time=1, reduce=T)))
format(object.size(simlist), units='Mb')

print(system.time(for(i in 1:10000){optimizeLP(bace)}))
format(object.size(simlist[[1]]), units='Mb')

bac <- Bac(specs, deathrate=0.2, duplirate=2, growthlimit=0.1, growtype='exponential')
org = Organism(model=specs)
print(system.time(for(i in 1:100){optimizeLP(bac)}))
print(system.time(for(i in 1:1000000){optimizeLP(org)}))

object.size(arena@media[[1]]@diffmat)
object.size(Matrix(arena@media[[1]]@diffmat, sparse = TRUE))

system.time(
addBac(arena, baclist=list(specs), amount=100*100)
)
format(object.size(arena),units='Mb')
object.size(arena)
#object.size(arena@orglist)
#object.size(arena@orglist[[1]])

#object.size(arena@orglist[[1]]@ubnd)
#object.size(arena@orglist[[1]]@lbnd)
addSubs(arena, smax=20)
arena@orglist[[1]]
simlist <- simulate(arena, time=2)

#3.569.236.768 bytes
#2177396768 bytes
#785556768 bytes
#397896128 bytes
#267256128 bytes
#70216336 bytes
#46.456.336 bytes

#3568800040 bytes
#2176960040 bytes
#785120040 bytes
#397440040 bytes
#266800040 bytes
#69760040 bytes
#46.000.040 bytes

#for(i in seq_along(specs)){
#  print(i)
#  addListBac(arena, bacmod=specs[[i]], amount=1)
#}
#addBac(arena, amount=10)

arena
arena@orglist[[1]]
#changeSub(arena, "EX_glc(e)", 20)
#changeSub(arena, "EX_h2o(e)", 20)
#changeSub(arena, "EX_o2(e)", 20)
#changeSub(arena, "EX_pi(e)", 20)
simlist <- simulate(arena, time=5)

lapply(popana@media, function(x){print(x@name)})

#growthc <- data.frame(time=seq_along(simlist))
#snames <- unique(unlist(lapply(simlist[[1]]@orglist,function(x){return(x@type)})))
for(i in seq_along(simlist)){
  popana <- as.matrix(simlist[[i]]@occmat)
  #popana <- as.matrix(simlist[[i]]@media[['EX_for(e)']]@diffmat)
  image(popana)
}

A <- as(regMat, "sparseMatrix")       # see also `vignette("Intro2Matrix")`
B <- Matrix(regMat, sparse = TRUE)  

object.size(arena@occmat)
object.size(Matrix(arena@occmat, sparse=TRUE))

object.size(simlist)
object.size(simlist[[1]])
object.size(simlist[[1]]@orglist)
object.size(simlist[[1]]@orglist[[1]])
object.size(simlist[[1]]@orglist[[1]]@lpobj)
object.size(simlist[[1]]@orglist[[1]]@fbasol)
object.size(simlist[[1]]@orglist[[1]]@lbnd)
object.size(simlist[[1]]@orglist[[1]]@ubnd)



