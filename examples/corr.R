setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
SYBIL_SETTINGS("SOLVER", "glpkAPI")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

seed <- sample(1:1000,1)
cat("seed:", seed)
set.seed(seed)

data(Ec_core)

ecore_anaerob = changeBounds(Ec_core, c('EX_o2(e)'), 0, 0)

bac1 = Bac(model=Ec_core, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
            speed=1, type="ec_aerob")
bac2 = Bac(model=ecore_anaerob, deathrate=0.05, duplirate=0.5, growthlimit=0.05, growtype="exponential",
            speed=0, type="ec_aenaerob")
arena = Arena(n=30, m=30, stir=F)
addOrg(arena, bac2, amount=10,x=1:10,y=1:10)
addOrg(arena, bac1, amount=10)
addSubs(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)", "EX_h2o(e)", "EX_h2o(e)", "EX_h(e)", "EX_nh4(e)"), difunc="cpp", difspeed=1)

sim <- simEnv(arena, time=100)

evalArena(sim)
plotCurves(sim)
plotCurves2(sim)


corr <- getCorrM(sim)
corrplot(corr, order = "hclust", tl.pos="n", type="upper")

checkCorr(sim, corr, c("co2"))
checkCorr(sim, corr, c("ec_aerob"))

library(corrplot)
corr2 <- getCorrM(sim, reactions=F)
corrplot(corr2)
corrplot(corr2, order = "hclust")
corrplot(corr2, method = "square", order = "hclust", type="upper")
corrplot(corr2, method = "number",tl.cex=0.1)

cex.before <- par("cex")
par(cex = 0.5)
corrplot.mixed(corr2, insig = "blank", method = "color",
         addCoef.col="grey", 
         order = "AOE", tl.cex = 1/par("cex"),
         cl.cex = 1/par("cex"), addCoefasPercent = TRUE)
corrplot.mixed(corr2, tl.cex=0.5/par("cex"), cl.cex = 1/par("cex"))
par(cex = cex.before)


#library(caret) # very big! just for findCorrelation??
#corr_high <- corr[,findCorrelation(corr, 0.999)]
#corrplot(corr_high, type="upper", tl.pos="n")
#corrplot.mixed(corr_high)
