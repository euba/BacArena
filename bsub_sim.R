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

library(sybilSBML)
model = readSBMLmod("Bs_iYO844_flux1.xml")
msgg = c("EX_k(e)","EX_mops(e)","EX_mg2(e)","EX_ca2(e)","EX_mn2(e)","EX_fe3(e)","EX_zn2(e)","EX_thym(e)",
         "EX_glyc(e)","EX_glu_L(e)","EX_co2(e)","EX_o2(e)","EX_pi(e)","EX_h2o(e)","EX_h(e)",
         "EX_so4(e)") #added as a sulfur source
modelB = changeBounds(model,model@react_id[grep("EX_",model@react_id)],lb=0)
modelB = changeBounds(modelB,msgg,lb=-c(5,100,2,0.7,0.05,0.1,0.001,0.002,68.4,29.6,1000,1000,1000,1000,1000,10))

optimizeProb(modelB)
optimizeProb(model)
#EX_so4(e),EX_pi(e),EX_o2(e),EX_na1(e),EX_nh4(e),EX_mg2(e),EX_k(e),EX_h(e),EX_h2o(e),EX_glc(e),EX_fe3(e),EX_co2(e),EX_ca2(e)

#minmed = model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)]<0)]
#modelB = changeBounds(model,model@react_id[grep("EX_",model@react_id)],lb=-10)
#modelB = changeBounds(model,setdiff(model@react_id[grep("EX_",model@react_id)],minmed),lb=-50)

bace1 = Bac(model=modelB, deathrate=0, duplirate=1, growthlimit=0.01, growtype="exponential",
            speed=2, type="Bsubtilis", lyse=F)
arena = Arena(n=100, m=100, stir=F, tstep=1)
addOrg(arena, bace1, amount=1, x=arena@n/2, y=arena@m/2, growth = 0.5)
addSubs(arena, smax=100000, difunc="cpp", difspeed=1)

print(system.time(evalsim <- simEnv(arena, time=30)))

evalArena(evalsim,phencol=T)

# find out what external metabolites (medium) the model has
metrans = model@met_name
names(metrans) = model@met_id
metrans[paste(gsub('EX_','',model@react_id[grep("EX",model@react_id)]),'[None]',sep='')]

minmed = model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)]
metrans[paste(gsub('EX_','',minmed),'[None]',sep='')]

set.seed(5000)
# paramters from MatNet:
# if cellMass >= 2
# cellDeathThreshold = 0
model@lowbnd[grep("EX",model@react_id)]
modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-10)
modelP = changeBounds(model,c("EX_EC0007", #Oxygen
                              "EX_EC0007", #Molybdate
                              "EX_EC0144", #Cobalt
                              "EX_EC0073", #Nitrite
                              "EX_EC0001", #H20
                              "EX_EC0236", #Nickel
                              "EX_EC0518", #Nitrogen
                              "EX_EC0034", #Zinc
                              "EX_EC0034", #Sodium
                              "EX_EC0994", #Cadmium
                              "EX_EC0056", #Copper
                              "EX_EC0957", #Amonium
                              "EX_EC0011", #CO2
                              "EX_EC0201", #Nitrate
                              "EX_EC9512", #Fe3
                              "EX_EC0048", #Sulfate
                              "EX_EC0030", #Manganese
                              "EX_EC0021", #Iron
                              "EX_EC0065", #Proton
                              "EX_EC0197" #Potassium
),lb=-Inf)

bace1 = Bac(model=modelP, deathrate=0.05, duplirate=1, growthlimit=0.05, growtype="exponential",
            speed=5, type="PAO", lyse=F)
arena = Arena(n=100, m=100, stir=F, tstep=0.5)
addOrg(arena, bace1, amount=1, x=arena@n/2, y=arena@m/2)
addSubs(arena, smax=50, difunc="cpp", difspeed=1, mediac=minmed)

print(system.time(evalsim <- simEnv(arena, time=100)))
