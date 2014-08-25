setwd("~/BacArena")
source(file="class/Arena.R")
load("/home/eugen/specsm.RData")

arena = Arena(n=50, m=50)
addListBac(arena, baclist=specs, amount=1)
#for(i in seq_along(specs)){
#  print(i)
#  addListBac(arena, bacmod=specs[[i]], amount=1)
#}
#addBac(arena, amount=10)
addSubs(arena, smax=20)
#changeSub(arena, "EX_glc(e)", 20)
#changeSub(arena, "EX_h2o(e)", 20)
#changeSub(arena, "EX_o2(e)", 20)
#changeSub(arena, "EX_pi(e)", 20)
simlist <- simulate(arena, time=30)
