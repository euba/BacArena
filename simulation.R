setwd("C:/Users/eugen.bauer/Documents/GitHub/BacArena/") #have to update orgn after deletion of bacteria?
source(file="class/Arena.R")
#load("/home/eugen/specsm.RData")

load("data/Bcoli_model.R")
specs = model

org = Organism(1,2,specs) #69352 bytes -> 17576 bytes

arena = Arena(n=10, m=10)
system.time(
addBac(arena, baclist=list(specs), amount=10*10)
)
format(object.size(arena),units='Mb')
object.size(arena)
#object.size(arena@orglist)
#object.size(arena@orglist[[1]])

#object.size(arena@orglist[[1]]@ubnd)
#object.size(arena@orglist[[1]]@lbnd)
addSubs(arena, smax=20)
arena@orglist[[1]]
simlist <- simulate(arena, time=5)

#3569236768 bytes
#2177396768 bytes
#785556768 bytes
#397896128 bytes
#267256128 bytes
#70216336 bytes
#46456336 bytes

#3568800040 bytes
#2176960040 bytes
#785120040 bytes
#397440040 bytes
#266800040 bytes
#69760040 bytes
#46000040 bytes

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

object.size(simlist)
object.size(simlist[[1]])
object.size(simlist[[1]]@orglist)
object.size(simlist[[1]]@orglist[[1]])
object.size(simlist[[1]]@orglist[[1]]@lpobj)
object.size(simlist[[1]]@orglist[[1]]@fbasol)
object.size(simlist[[1]]@orglist[[1]]@lbnd)
object.size(simlist[[1]]@orglist[[1]]@ubnd)



