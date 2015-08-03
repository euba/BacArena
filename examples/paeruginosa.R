setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(sybilSBML)
SYBIL_SETTINGS("SOLVER", "glpkAPI")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

model = readSBMLmod("modelPOA_named_EX.xml")

#model@lowbnd[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] == 0)] <- -1000 # set not defined exchange reactions to possible output


# paramters from MatNet:
# if cellMass >= 2
# cellDeathThreshold = 0

#model@lowbnd[grep("EX",model@react_id)]
#modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-10)

bac = Bac(model=model, deathrate=0.05, duplirate=1, growthlimit=0.05, growtype="exponential",
          speed=2, type="PAO", lyse=F)

#arena = Arena(n=100, m=100, stir=F, tstep=1)
#addOrg(arena, bac, amount=22, x=c(10:20,80:90), y=rep(1,22))
#addOrg(arena, bac, amount=20, x=c(1:10,91:100), y=c(1:10,91:100))

arena = Arena(n=10, m=10, stir=F, tstep=1)
#addOrg(arena, bac, amount=22, x=c(10:20,80:90), y=rep(1,22))
addOrg(arena, bac, amount=10, x=c(1:5,6:10), y=c(1:5,6:10))


#arena = Arena(n=10, m=10, stir=F, tstep=0.5)
#addOrg(arena, bac, amount=6, x=c(1:3,8:10), y=rep(1,6))
addSubs(arena, smax=50, difunc="cpp", difspeed=1,
        mediac=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)])

sim <- simEnv(arena, time=500)

plotCurves2(sim, ignore=c("EX_H2O","EX_H"), num=6)
evalArena(sim, time=1, phencol = T)



library(animation)
oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
saveVideo({
  ani.options(interval = 0.5)
  evalArena(evalsim, phencol=T, plot_items=c('population','EX_EC0007'))
},video.name = "PAO_pop_phen_bio2.avi", other.opts = "-b 600k")

library(animation)
oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
saveVideo({
  ani.options(interval = 0.5)
  evalArena(evalsim, phencol=T)
},video.name = "PAO_pop5.avi", other.opts = "-b 600k")
