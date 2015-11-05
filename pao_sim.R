##Loading the required packages## 
library(sybil)
#library(sybilSBML)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(compiler)
#setwd("/Users/euba/GitRep/BacArena/")
setwd("uni/bacarena")
#setwd('P:/GitRep/BacArena')
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")

#SYBIL_SETTINGS("SOLVER","sybilGUROBI") #setting solver to GUROBI
SYBIL_SETTINGS("SOLVER") #setting solver to GUROBI

library(sybilSBML)
#model = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/modelPOA.xml")
model = readSBMLmod("modelPOA.xml")
load("poa_model.RData")

set.seed(100)
# paramters from MatNet:
# if cellMass >= 2
# cellDeathThreshold = 0
#model@lowbnd[grep("EX",model@react_id)]
#model = changeBounds(model,"EX_EC0027",lb=-10)
#model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)]==-10)]
modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-1000)
#modelP = model

modelP = changeBounds(modelP,"EX_EC0011",lb=0)
modelP = changeBounds(modelP,"EX_EC0065",lb=0)
modelP = changeBounds(modelP,"EX_EC0027",lb=-10)

# modelPm = changeBounds(modelP,"EX_EC0036",lb=0)
# #modelPm = changeBounds(modelPm,"EX_EC0065",lb=0)
# modelPm = changeBounds(modelPm,"EX_EC0957",ub=0)
# modelPm = changeBounds(modelPm,"EX_EC0115",lb=0)
# #modelPm = changeBounds(modelPm,"EX_EC0035",lb=0)

bace = Bac(model=modelP, deathrate=0.05, duplirate=1, growthlimit=0.05, growtype="exponential",
            speed=2, type="PAO", lyse=F)
arena = Arena(n=150, m=150, stir=F, tstep=0.5)
addOrg(arena, bace, amount=1, x=75, y=75)
addSubs(arena, smax=50, difunc="cpp", difspeed=1,
        mediac=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)])

#print(system.time(evalsim <- simEnv(arena, time=50)))

preps = list()
for(i in 1:20){
  print(i)
  preps[[i]] <- simEnv(arena, time=60)
}

evalArena(evalsim, phencol=T)

par(mfrow=c(2,3))
evalArena(evalsim, phencol=T, time = c(20,40,50))
evalArena(evalsim_m, phencol=T, time = c(20,40,50))
par(mfrow=c(1,1))
# 
# library(animation)
# oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
# saveVideo({
#   ani.options(interval = 0.5)
#   evalArena(evalsim, phencol=T, plot_items=c('population','EX_EC0007'))
# },video.name = "PAO_pop_phen_bio2.avi", other.opts = "-b 600k")
# 
# library(animation)
# oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
# saveVideo({
#   ani.options(interval = 0.5)
#   evalArena(evalsim, phencol=T)
# },video.name = "PAO_pop5.avi", other.opts = "-b 600k")

#################################################################
######################## Phenotype analysis
#################################################################
cmp = 'EX_EC0036'
#evalArena(evalsim,plot_items = cmp)
i = 35
#max(matrix(evalsim@medlist[[i]][[cmp]],100,100))

#interesting compounds: EX_EC0065, EX_EC0011,EX_EC0036, EX_EC0021

#evalsim@simlist[[i]]$phenotype
#getPhenoMat(evalsim,time=i)[,cmp]
#table(getPhenoMat(evalsim,time=i)[,cmp])

#plot ammonia uptake phenotypes
nuse = getPhenoMat(evalsim)[,cmp] #other interesting compounds: 'EX_h(e)','EX_pi(e)','EX_man1p(e)','EX_man6p(e)','EX_chor(e)','EX_succ(e)','EX_fum(e)','EX_for(e)','EX_cit(e)','EX_6pgc(e)','EX_acac(e)','EX_pep(e)','EX_btd_RR(e)','EX_ac(e)','EX_ppa(e)','EX_dha(e)','EX_lac_L(e)','EX_pyr(e)','EX_tyr_L(e)','EX_thym(e)','EX_glyclt(e)'
pop = evalsim@simlist[[i]]
pop$phenotype_n = 1
pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1

#pop$phenotype_n[which(pop$phenotype!=0)]
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=1, col=pop$phenotype_n, pch=19)

phenmat <- getPhenoMat(evalsim)
trans = modelP@met_name
names(trans) = paste('EX',gsub('\\[None\\]','',modelP@met_id),sep='_')
phensel <- phenmat[which(phenmat[,cmp]==1),]
colnames(phensel) <- trans[colnames(phensel)]
phensel[,which(apply(phensel,2,sum)!=0)]

length(which(sort(evalsim@subchange)!=0))

par(mfrow=c(1,1))

par(mfrow=c(3,5))
for(i in 1:length(evalsim@simlist)){
  for(sub in names(which(sort(evalsim@subchange)!=0))){
    nuse = getPhenoMat(evalsim)[,sub]
    pop = evalsim@simlist[[i]]
    pop$phenotype_n = 1
    pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
    plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
         axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19, main=sub)
  }
}

par(mfrow=c(4,5))
for(sub in names(which(sort(evalsim@subchange)!=0))){
  nuse = getPhenoMat(evalsim)[,sub]
  pop = evalsim@simlist[[i]]
  pop$phenotype_n = 1
  pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
  plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
       axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19, main=sub)
}

modelP@met_name[which(modelP@met_id=='EC0115[None]')]
modelP@met_name[which(modelP@met_id=='EC0036[None]')]

modelP@met_name[which(modelP@met_id=='EC0029[None]')]
modelP@met_name[which(modelP@met_id=='EC0065[None]')]

modelP@met_name[which(modelP@met_id=='EC0071[None]')]
modelP@met_name[which(modelP@met_id=='EC0957[None]')]
modelP@met_name[which(modelP@met_id=='EC0115[None]')]
modelP@met_name[which(modelP@met_id=='EC0035[None]')]
modelP@met_name[which(modelP@met_id=='EC0029[None]')]

par(mfrow=c(4,4))
for(sub in names(which(sort(evalsim_m@subchange)!=0))){
  nuse = getPhenoMat(evalsim_m)[,sub]
  pop = evalsim_m@simlist[[i]]
  pop$phenotype_n = 1
  pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
  plot(pop[,c('x','y')],xlim=c(0,evalsim_m@n),ylim=c(0,evalsim_m@m),xlab='',ylab='',
       axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19, main=sub)
}

#################################################################
######################## Analysis
#################################################################

par(mfrow=c(3,3))
par(mfrow=c(1,1))
#c(20,35,50)
#cmp = 'EX_EC0036'
cmp = 'EX_EC0029'
i = 50
nuse = getPhenoMat(evalsim)[,cmp] #other interesting compounds: 'EX_h(e)','EX_pi(e)','EX_man1p(e)','EX_man6p(e)','EX_chor(e)','EX_succ(e)','EX_fum(e)','EX_for(e)','EX_cit(e)','EX_6pgc(e)','EX_acac(e)','EX_pep(e)','EX_btd_RR(e)','EX_ac(e)','EX_ppa(e)','EX_dha(e)','EX_lac_L(e)','EX_pyr(e)','EX_tyr_L(e)','EX_thym(e)','EX_glyclt(e)'
pop = evalsim@simlist[[i]]
pop$phenotype_n = 1
pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19)

par(mfrow=c(3,3))
par(mfrow=c(1,1))
#c(20,35,50)
#cmp = 'EX_EC0036'
cmp = 'EX_EC0029'
i = 50
nuse = getPhenoMat(evalsim_m)[,cmp] #other interesting compounds: 'EX_h(e)','EX_pi(e)','EX_man1p(e)','EX_man6p(e)','EX_chor(e)','EX_succ(e)','EX_fum(e)','EX_for(e)','EX_cit(e)','EX_6pgc(e)','EX_acac(e)','EX_pep(e)','EX_btd_RR(e)','EX_ac(e)','EX_ppa(e)','EX_dha(e)','EX_lac_L(e)','EX_pyr(e)','EX_tyr_L(e)','EX_thym(e)','EX_glyclt(e)'
pop = evalsim_m@simlist[[i]]
pop$phenotype_n = 1
pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
plot(pop[,c('x','y')],xlim=c(0,evalsim_m@n),ylim=c(0,evalsim_m@m),xlab='',ylab='',
     axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19)
