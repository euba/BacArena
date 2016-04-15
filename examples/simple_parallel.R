library(devtools)
install_github(repo="euba/bacarena")
#install_local("~/uni/bacarena")
library(BacArena)
library(parallel)

# create cluster with all available cores
#nr_cores <- detectCores()
nr_cores <- 4
print(paste("cores available:", nr_cores))
cl <- makeCluster(nr_cores, type="PSOCK") # PSOCK works with win/mac/lin

# variables which will be used in cluster have to be exported
bcoli <- readRDS(file = "data/bcoli_orth.RDS")
clusterExport(cl, "bcoli")

print(system.time(simlist <- parLapply(cl, 1:nr_cores, function(i){
  bac   <- BacArena::Bac(model=bcoli)
  arena <- BacArena::Arena(n=50, m=50, stir=F, Lx=0.0125, Ly=0.0125)
  arena <- BacArena::addOrg(arena, bac, amount=50)
  arena <- BacArena::addMinMed(arena, bac)
  arena <- BacArena::addSubs(arena, mediac="EX_o2(e)", smax=3, unit="fmol/cell", add = F)
  arena <- BacArena::rmSubs(arena, mediac="EX_co2(e)")
  #sim   <- BacArena::simEnv(arena, time=10)   # 44.845 
  sim   <- BacArena::simEnv(arena, time=10, diff_par = TRUE, cl_size = 2)   # 104.373 
}) ))
stopCluster(cl)

subs=c("EX_glc(e)", "EX_ac(e)", "EX_o2(e)", "EX_etoh(e)", "EX_co2(e)")
lapply(simlist, BacArena::plotCurves2, growthCurve=FALSE, subs=subs)
