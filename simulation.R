setwd("~/BacArena") #have to update orgn after deletion of bacteria?
source(file="class/Arena.R")
#load("/home/eugen/specsm.RData")

load("data/Bcoli_model.R")
specs = model

#org = Organism(1,2,specs) #69352 bytes -> 17576 bytes

arena = Arena(n=100, m=100)
addBac(arena, baclist=list(specs), amount=100*100)
#for(i in seq_along(specs)){
#  print(i)
#  addListBac(arena, bacmod=specs[[i]], amount=1)
#}
#addBac(arena, amount=10)
addSubs(arena, smax=20)

arena
arena@orglist[[1]]
#changeSub(arena, "EX_glc(e)", 20)
#changeSub(arena, "EX_h2o(e)", 20)
#changeSub(arena, "EX_o2(e)", 20)
#changeSub(arena, "EX_pi(e)", 20)
simlist <- simulate(arena, time=5)

lapply(popana@media, function(x){print(x@name)})

growthc <- data.frame(time=seq_along(simlist))
snames <- unique(unlist(lapply(simlist[[1]]@orglist,function(x){return(x@type)})))
for(i in seq_along(simlist)){
  popana <- simlist[[i]]
  image(popana@occmat)
  #image(popana@media[[7]]@diffmat)
  popocc <- table(popana@occmat)[-1]
  for(j in seq_along(popocc)){
    growthc[i,snames[j]] <- popocc[j]
  }
}

A <- as(regMat, "sparseMatrix")       # see also `vignette("Intro2Matrix")`
B <- Matrix(regMat, sparse = TRUE)  

object.size(arena@occmat)
object.size(Matrix(arena@occmat, sparse=TRUE))
