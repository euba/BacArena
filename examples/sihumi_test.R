library(BacArena)
SYBIL_SETTINGS("SOLVER", "cplexAPI")

#sihumi <- readRDS("~/uni/dat/mod/r/sihumi.RDS")
#latest_agora <- readRDS("~/uni/dat/mod/r/agora1.01.western.SWcorr.RDS")
#sihumi.idx <- match(sapply(sihumi, mod_desc), sapply(latest_agora, mod_desc))
#sihumi_new <- latest_agora[sihumi.idx]
#saveRDS(sihumi_new, "~/uni/dat/mod/r/sihumi_agora1.01.western.SWcorr.RDS")

sihumi <- readRDS("~/uni/dat/mod/r/sihumi_agora1.01.western.SWcorr.RDS")

arena <- Arena(n=20, m=20)
for(i in seq_along(sihumi)){
  bac <- Bac(model=sybil::upgradeModelorg(sihumi[[i]]))
  arena <- addOrg(arena, specI=bac, amount=5)
  arena <- addDefaultMed(arena, org = bac, unit = "mmol/arena")
}
sihumi_test <- simEnv(arena, time=7, sec_obj = "mtf", with_shadow = T)

#save(sihumi_test, file = "~/uni/bacarena/data/sihumi_test.rda", compress="xz")
#load("~/uni/bacarena/data/sihumi_test.rda")

plotGrowthCurve(sihumi_test)[[1]]
plotSpecActivity(sihumi_test)[[2]]

findFeeding3(sihumi_test, time=5, mets=names(head(getVarSubs(sihumi_test),30)) )


BacArena::plotShadowCost(sihumi_test, spec_nr=7)[[2]]
sihumi_test@shadowlist[[4]]$Bacteroides_thetaiotaomicron_VPI_5482
