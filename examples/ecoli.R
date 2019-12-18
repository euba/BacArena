library(BacArena)
data("Ec_core")

SYBIL_SETTINGS("SOLVER", "cplexAPI")
#SYBIL_SETTINGS("SOLVER", "glpkAPI")

bac <- Bac(Ec_core)
arena <- Arena(n=20, m=20)
arena <- addOrg(arena,bac,amount=20)
arena <- addSubs(arena, smax=0.005, mediac="EX_glc(e)", unit="mM")
arena <- addSubs(arena, smax=1, mediac=c("EX_pi(e)", "EX_h2o(e)",
                                         "EX_o2(e)", "EX_nh4(e)"), unit="mM") 
arena@removeM[c(1,20),1:20] <- 1
arena@removeM[1:20,c(1,20)] <- 1
image(arena@removeM)
sim   <- simEnv(arena, time=5, with_shadow = T, sec_obj = F)



getSubHist(sim, "EX_glc(e)")

sim@shadowlist[[3]]
plotShadowCost(sim, 1)[[2]] # fru, glc -2

plotGrowthCurve(sim.old)

sim <- simEnv(arena, time=2)
arena2 <- getArena(sim,2)
arena2@orgdat <- arena2@orgdat[-sample(nrow(arena2@orgdat), round(nrow(arena2@orgdat) / 10)), ] #remove ~10% of individuals randomly 
simEnv(arena2, time=2)

arena2 <- getArena(sim,2)
neworgdat <- arena2@orgdat #get the current orgdat
neworgdat <- neworgdat[-sample(nrow(neworgdat), round(nrow(neworgdat) / 10)), ] #remove ~10% of individuals randomly 
arena2 <- changeOrg(arena2,neworgdat)
