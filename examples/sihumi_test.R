library(BacArena)

sihumi <- readRDS("~/uni/dat/mod/r/sihumi.RDS")

arena <- Arena(n=20, m=20)

for(i in seq_along(sihumi)){
  bac <- Bac(model=sihumi[[i]])
  arena <- addOrg(arena, specI=bac, amount=5)
}
arena <- addSubs(arena, smax=5, unit="mM")

sim <- simEnv(arena, time=10)

save(sim, file = "~/uni/bacarena/data/sihumi_test.rda")
