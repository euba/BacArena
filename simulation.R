setwd("~/BacArena")
source(file="class/Arena.R")
load("data/ecore_model.R")
mod <- model

arena = Arena(n=100, m=100)
addBac(arena, amount=10)
addSubs(arena, 0)
changeSub(arena, "EX_glc(e)", 20)
changeSub(arena, "EX_h2o(e)", 20)
changeSub(arena, "EX_o2(e)", 20)
changeSub(arena, "EX_pi(e)", 20)


simlist <- list()
for(i in 1:100){
  simlist[[i]] <- arena
  print(system.time(for(j in seq_along(arena@media)){
    #diffuseNaiveR(arena@media[[j]])
    diffuseNaiveCpp(arena@media[[j]]@diffmat, donut=FALSE)
  }))
  j = 0
  orgl <- arena@orglist
  print(system.time(while(j+1 <= length(orgl) && j+1 <= length(arena@orglist)){
    j<-j+1
    move(arena@orglist[[j]],arena)
    medcon = getmed(arena,arena@orglist[[j]]@x,arena@orglist[[j]]@y)
    constrain(arena@orglist[[j]], names(medcon), lb=-medcon)
    optimizeLP(arena@orglist[[j]])
    arena@media = consume(arena@orglist[[j]],arena@media)
    growth(arena@orglist[[j]], arena, j,lifecosts=0.1)
  }))
  if(length(arena@orglist)==0){
    print("All bacs dead!")
    break
  } 
  cat("iter:", i, "bacs:",length(arena@orglist),"\n\n")
}
