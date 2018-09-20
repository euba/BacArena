library(BacArena)
data("Ec_core")

bac <- Bac(Ec_core)
arena <- Arena(n=20, m=20)
arena <- addOrg(arena,bac,amount=20)
arena <- addSubs(arena, smax=0.005, mediac="EX_glc(e)", unit="mM")
arena <- addSubs(arena, smax=1, mediac=c("EX_pi(e)", "EX_h2o(e)",
                                         "EX_o2(e)", "EX_nh4(e)"), unit="mM") 
sim   <- simEnv(arena, time=5, with_shadow = T, sec_obj = T)

getSubHist(sim, "EX_glc(e)")

sim@shadowlist[[3]]
plotShadowCost(sim, 1)[[2]]
