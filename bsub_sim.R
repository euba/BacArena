# load libraries and other R files to have everything in place
setwd('P:/GitRep/BacArena')
setwd("/Users/euba/GitRep/BacArena/")
setwd("~/uni/bacarena")
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(compiler) # byte code 
SYBIL_SETTINGS("SOLVER", "sybilGUROBI")
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

setwd("E:/BACARENA/Comparison/MatNet/P_aeruginosa/")
setwd("/Volumes/PHD/BACARENA/B_subtilis/")
setwd("P:/BACARENA/B_subtilis/")
setwd("/Users/euba/GitRep/BacArena/B_subtilis/")

library(sybilSBML)
model = readSBMLmod("Bs_iYO844_flux1.xml")
msgg = c("EX_k(e)","EX_mops(e)","EX_mg2(e)","EX_ca2(e)","EX_mn2(e)","EX_fe3(e)","EX_zn2(e)","EX_thym(e)",
         "EX_glyc(e)","EX_glu_L(e)","EX_co2(e)","EX_o2(e)","EX_pi(e)","EX_h2o(e)","EX_h(e)",
         "EX_so4(e)") #added as a sulfur source
modelB = changeBounds(model,model@react_id[grep("EX_",model@react_id)],lb=0)
#modelB = changeBounds(modelB,msgg,lb=-c(5,100,2,0.7,0.05,0.1,0.001,0.002,68.4,29.6,1000,1000,1000,1000,1000,10))
modelB = changeBounds(model,model@react_id[grep("EX_",model@react_id)],lb=-20)
modelB = changeBounds(modelB,c("EX_co2(e)","EX_o2(e)","EX_pi(e)","EX_h2o(e)","EX_h(e)","EX_no2(e)","EX_no3(e)"),lb=-1000)

optimizeProb(modelB)
optimizeProb(model)
#EX_so4(e),EX_pi(e),EX_o2(e),EX_na1(e),EX_nh4(e),EX_mg2(e),EX_k(e),EX_h(e),EX_h2o(e),EX_glc(e),EX_fe3(e),EX_co2(e),EX_ca2(e)

#minmed = model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)]<0)]
#modelB = changeBounds(model,model@react_id[grep("EX_",model@react_id)],lb=-10)
#modelB = changeBounds(model,setdiff(model@react_id[grep("EX_",model@react_id)],minmed),lb=-50)

bace1 = Bac(model=modelB, deathrate=0.1, duplirate=1, growthlimit=0.01, growtype="exponential",
            speed=1, type="Bsubtilis", lyse=F)
arena = Arena(n=100, m=100, stir=F, tstep=0.1)
addOrg(arena, bace1, amount=1, x=arena@n/2, y=arena@m/2, growth = 0.5)
#addSubs(arena, smax=10000, difunc="cpp", difspeed=1, mediac = msgg)
addSubs(arena, difunc="cpp", difspeed=1, mediac = c(msgg),
        smax = c(5,100,2,0.7,0.05,0.1,0.001,0.002,68.4,29.6,1000,1000,1000,1000,1000,10))

print(system.time(evalsim <- simEnv(arena, time=90)))

evalArena(evalsim,phencol=T)
plotCurves(evalsim)
plotCurves(evalsim,medplot=c('EX_nh4(e)','EX_glu_L(e)'))
plotCurves(evalsim,medplot=c('EX_nh4(e)'))
plotCurves(evalsim,medplot=c('EX_glu_L(e)'))
evalArena(evalsim,plot_items = c('Population'))
evalArena(evalsim,plot_items = c('EX_nh4(e)'))
evalArena(evalsim,plot_items = c('EX_glu_L(e)'))
evalArena(evalsim,phencol=T,plot_items = c('Population',"EX_glu_L(e)","EX_gln_L(e)","EX_nh4(e)"))
