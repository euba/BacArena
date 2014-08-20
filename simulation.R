setwd("~/BacArena")
source(file="class/Arena.R")

arena = Arena(n=10, m=10)
addBac(arena, amount=50)

#load ecoli core model to play around
load("data/ecore_model.R")
mod <- model
specs=list()
specs[[1]]=mod

simlist <- list()

addSub(pop, "EX_glc(e)", 20)
addSub(pop, "EX_h2o(e)", 20)
addSub(pop, "EX_o2(e)", 20)
addSub(pop, "EX_pi(e)", 20)

#pop = Population(specs, specn=rep(10, length(specs)), n=100, m=100)

for(i in 1:10){
  simlist[[i]] <- arena
  print(system.time(for(j in seq_along(arena@media)){
    #diffuseNaiveR(arena@media[[j]])
    diffuseNaiveCpp(arena@media[[j]]@diffmat, donut=FALSE)
  }))
  j = 0
  print(system.time(while(j+1 <= length(arena@orglist)){
    j<-j+1
    move(arena@orglist[[j]],arena)
    medcon = getmed(arena,arena@orglist[[j]]@x,arena@orglist[[j]]@y)
    constrain(arena@orglist[[j]], names(medcon), lb=-medcon)
    optimizeLP(arena@orglist[[j]])
    arena@media = consume(arena@orglist[[j]],arena@media)
    growth(arena@orglist[[j]], arena, j,lifecosts=0.6)
  }))
  if(length(arena@orglist)==0){
    print("All bacs dead!")
    break
  } 
  cat("iter:", i, "bacs:",length(arena@orglist),"\n\n")
}
