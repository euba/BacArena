##Loading the required packages## 
library(ReacTran)
library(sybil)
library(sybilSBML)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
#setwd("/Users/euba/GitRep/BacArena/")
setwd("~/uni/bacarena")
#setwd('P:/GitRep/BacArena')
#setwd('P:/GitRep/BacArena')

source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")

#library(sybilSBML)
#model = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/modelPOA.xml")
#model = readSBMLmod("modelPOA.xml")
#load("poa_model.RData")
library(cplexAPI)
SYBIL_SETTINGS("SOLVER","cplexAPI") #setting solver to GUROBI
#model = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/modelPOA.xml")
load('P:/GitRep/BacArena/model/poa_model_new.RData')
#load('/Users/euba/GitRep/BacArena/poa_model.RData')

medium = read.csv('P:/GitRep/BacArena/model/poa_minimal_medium.csv')
modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-1000)
#modelP = changeBounds(modelP,"EX_EC0027",lb=-10)
#modelP@met_name[which(modelP@met_id=='EC0027[None]')]

# bac = Bac(model=modelP, growtype="exponential", lyse=F)
# setKinetics(bac, exchangeR="EX_EC0027", Km=0.01, vmax=7.56)
# arena = Arena(n=200, m=200, stir=F, seed=8904, Lx=0.05, Ly=0.05, tstep=0.5)
# addOrg(arena, bac, amount=1, x=arena@n/2, y=arena@m/2,growth = 0.9)
# addSubs(arena, smax=0, difspeed=1.675e-6, unit='mM')
# #addSubs(arena, smax=0.05, difspeed=1.675e-6, unit='mM',
# #        mediac=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)])
# addSubs(arena, smax=medium$Concentration, difspeed=medium$diffusion.constant, unit='mM', 
#         mediac=as.character(medium$Reaction.ID))
# #addSubs(arena, smax=10, mediac=c("EX_o2(e)","EX_h(e)","EX_co2(e)","EX_o2(e)","EX_pi(e)"), difunc="pde", difspeed=rep(0.072,5))
# #createGradient(arena,smax=20,mediac="EX_o2(e)",position='left',steep=0.5)
# #createGradient(arena,smax=20,mediac=arena@mediac,position='left',steep=0.5)
# print(system.time(sim <- simEnv(arena, time=96)))

simlist = list()
for(i in 1:5){
  print(i)
  bac = Bac(model=modelP, growtype="exponential", deathrate=0)
  setKinetics(bac, exchangeR="EX_EC0027", Km=0.01, vmax=7.56)
  arena = Arena(n=200, m=200, stir=F, seed=i*100, Lx=0.05, Ly=0.05, tstep=0.5)
  addOrg(arena, bac, amount=1, x=arena@n/2, y=arena@m/2,growth = 0.9)
  addSubs(arena, smax=0, difspeed=1.675e-6, unit='mM')
  addSubs(arena, smax=0.05, difspeed=medium$diffusion.constant, unit='mM', 
          mediac=as.character(medium$Reaction.ID))
  print(system.time(simlist[[i]] <- simEnv(arena, time=96)))
}
sim = simlist[[1]]
evalsim = sim
#save(simlist, file='simlist.RData')

modelP@met_name[which(modelP@met_id=='EC0029[None]')]


write.csv(cbind(rxn=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)],
gsub('\\[e\\]','',sapply(paste(gsub('EX_','',rxn),'[None]',sep=''),function(x,mod){mod@met_name[which(mod@met_id == x)]},mod=modelP))),
file='Minimal medium.csv')
# modelPm = changeBounds(modelP,"EX_EC0036",lb=0)
# #modelPm = changeBounds(modelPm,"EX_EC0065",lb=0)
# modelPm = changeBounds(modelPm,"EX_EC0957",ub=0)
# modelPm = changeBounds(modelPm,"EX_EC0115",lb=0)
# #modelPm = changeBounds(modelPm,"EX_EC0035",lb=0)

bace = Bac(model=modelP, deathrate=0.05, duplirate=1, growthlimit=0.05, growtype="exponential",
            speed=2, type="PAO", lyse=F)
arena = Arena(n=100, m=100, stir=F, tstep=0.5)
addOrg(arena, bace, amount=1, x=50, y=50)
addSubs(arena, smax=50, difunc="pde", difspeed=1,
        mediac=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)])

print(system.time(evalsim <- simEnv(arena, time=50)))

preps = list()
for(i in 1:20){
  print(i)
  preps[[i]] <- simEnv(arena, time=60)
}

plotCurves(sim)
evalArena(sim, phencol=T,plot_items = c('Population','EX_EC0029',"EX_EC0027","EX_EC0007"),time=96)
max(sim@simlist[[97]]$growth)

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
    pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('modelPOA.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
    plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
         axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19, main=sub)
  }
}

par(mfrow=c(4,4))
for(sub in names(which(sort(evalsim@subchange)!=0))){
  nuse = getPhenoMat(evalsim)[,sub]
  pop = evalsim@simlist[[i]]
  pop$phenotype_n = 1
  pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('modelPOA.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
  plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
       axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19, main=sub)
}

modelP@met_name[which(modelP@met_id=='EC0115[None]')]
modelP@met_name[which(modelP@met_id=='EC0036[None]')]

modelP@met_name[which(modelP@met_id=='EC0029[None]')]
modelP@met_name[which(modelP@met_id=='EC0065[None]')]

modelP@met_name[which(modelP@met_id=='EC0021[None]')]
modelP@met_name[which(modelP@met_id=='EC0957[None]')]
modelP@met_name[which(modelP@met_id=='EC0115[None]')]
modelP@met_name[which(modelP@met_id=='EC0035[None]')]
modelP@met_name[which(modelP@met_id=='EC0007[None]')]

par(mfrow=c(4,4))
for(sub in names(which(sort(evalsim_m@subchange)!=0))){
  nuse = getPhenoMat(evalsim_m)[,sub]
  pop = evalsim_m@simlist[[i]]
  pop$phenotype_n = 1
  pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('modelPOA.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
  plot(pop[,c('x','y')],xlim=c(0,evalsim_m@n),ylim=c(0,evalsim_m@m),xlab='',ylab='',
       axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19, main=sub)
}

#################################################################
######################## Analysis
#################################################################

#evalsim=simlist[[q]] #4
par(mfrow=c(3,3))
par(mfrow=c(1,1))
#c(20,35,50)
#cmp = 'EX_EC0214'
#cmp = 'EX_EC0098'
#cmp = 'EX_EC0359'
cmp = 'EX_EC0029'
cmp = 'EX_EC0036'
cmp = 'EX_EC0027'
cmp = 'EX_EC0007'
#cmp = 'EX_EC0007'
#cmp = 'EX_EC0021'
modelP@met_name[which(modelP@met_id=='EC0957[None]')]
modelP@met_name[which(modelP@met_id=='EC0065[None]')]
modelP@met_name[which(modelP@met_id=='EC0029[None]')]
modelP@met_name[which(modelP@met_id=='EC0036[None]')]
modelP@met_name[which(modelP@met_id=='EC0011[None]')]
modelP@met_name[which(modelP@met_id=='EC0009[None]')]
i = 72
nuse = getPhenoMat(evalsim)[,cmp] #other interesting compounds: 'EX_h(e)','EX_pi(e)','EX_man1p(e)','EX_man6p(e)','EX_chor(e)','EX_succ(e)','EX_fum(e)','EX_for(e)','EX_cit(e)','EX_6pgc(e)','EX_acac(e)','EX_pep(e)','EX_btd_RR(e)','EX_ac(e)','EX_ppa(e)','EX_dha(e)','EX_lac_L(e)','EX_pyr(e)','EX_tyr_L(e)','EX_thym(e)','EX_glyclt(e)'
pop = evalsim@simlist[[i]]
pop$phenotype_n = 1
pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('modelPOA.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=(pop$growth/max(pop$growth)), col=pop$phenotype_n, pch=19)
hist(pop$growth,breaks=100)
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=(pop$growth/max(pop$growth)), col=pop$phenotype+1, pch=19)

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


##############################################################

load('P:/BACARENA/RESULTS/POA/simlist_poa.RData')

library(scales)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

modelP
pmat = getPhenoMat(simlist[[3]])

pop = matrix(0,nrow=length(simlist),ncol=length(simlist[[3]]@simlist))
phenlist = lapply(c('dead',apply(getPhenoMat(simlist[[3]]),1,function(x){paste(x,collapse='')})),function(x){
  return(matrix(0,nrow=length(simlist),ncol=length(simlist[[3]]@simlist)))
})
names(phenlist) = c('dead',apply(getPhenoMat(simlist[[3]]),1,function(x){paste(x,collapse='')}))
subs = lapply(which(simlist[[3]]@subchange!=0),function(x){
  return(matrix(0,nrow=length(simlist),ncol=length(simlist[[3]]@simlist)))
})
for(i in 1:length(simlist)){
  dat = plotCurves(simlist[[i]],retdata=T) 
  pop[i,] = dat$Population
  pheno = c('dead',apply(getPhenoMat(simlist[[i]]),1,function(x){paste(x,collapse='')}))
  names(pheno)[1] = 'modelPOA.0'
  for(j in names(which(simlist[[1]]@subchange!=0))){
    subs[[j]][i,] = (dat$Substance[j,])/(0.01*6.25e-08)
  }
  for(j in 1:length(pheno)){
    for(k in 2:length(simlist[[i]]@simlist)){
      ptab = table(simlist[[i]]@simlist[[k]]$phenotype)
      if(gsub("modelPOA.",'',names(pheno)[j]) %in% names(ptab)){
        #print(phenlist[[pheno[j]]])
        phenlist[[pheno[j]]][i,k] = ptab[gsub("modelPOA.",'',names(pheno)[j])]
      }
    }
  }
}
names(phenlist) = 0:(length(phenlist)-1)
phenlist[['7']] = phenlist[['7']]+phenlist[['9']]
phenlist[['14']] = phenlist[['1']]+phenlist[['3']]+phenlist[['4']]+phenlist[['10']]+
  phenlist[['11']]+phenlist[['12']]+phenlist[['13']]+phenlist[['14']]
phenlist = phenlist[-which(names(phenlist) %in% c('1','3','4','9','10','11','12','13'))]
names(phenlist) = c('Metabolically inactive','Glucose respiration','Succinate, Acetate, CO2 production',
                    'Acetate production / Succinate consumption','Succinate, Acetate consumption',
                    'Succinate, Acetate production','Other phenotypes')
pcol=c('black',"#96CED5","#76D73C","#D4477D","#D19EC4","#C17E73","#CBC594")

pgrowth = apply(pop,2,mean)
psd = apply(pop,2,sd)
sconc = do.call(rbind,lapply(subs,function(x){apply(x,2,mean)}))
ssd = do.call(rbind,lapply(subs,function(x){apply(x,2,sd)}))

do.call(rbind,lapply(subs,function(x){apply(x,2,mean)}))

#phendat = melt(rbind(do.call(rbind,lapply(phenlist,function(x){apply(x,2,mean)})),apply(pop,2,mean)))
phendat = melt(rbind(do.call(rbind,lapply(phenlist,function(x){apply(x,2,mean)}))))
colnames(phendat) = c('Phenotypes','time','mean')
#phendat$sd = melt(rbind(do.call(rbind,lapply(phenlist,function(x){apply(x,2,sd)})),apply(pop,2,mean)))[,3]
phendat$sd = melt(do.call(rbind,lapply(phenlist,function(x){apply(x,2,function(x){sd(x)/sqrt(length(x))})})))[,3]
phendat$time = (phendat$time-1)/2
#phendat = phendat[-which(phendat$phens %in% c(3,4,10:14)),]
ggplot(phendat, aes(x=time, y=mean, colour=Phenotypes, group=Phenotypes)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.3,size=1) +
  geom_line(size=1.2) +
  #geom_point(size=5, shape=20) + # 21 is filled circle
  xlab("Time in h") +
  ylab("Number of bacterial cells") +
  geom_vline(xintercept = c(17,20,30,48), colour="black", linetype = "longdash",size=1.5) +
  scale_y_continuous(label=scientific_format()) +
  scale_color_manual(values=pcol) +
  theme_bw(base_size = 30) +
  theme(#legend.position='none',
    legend.text=element_text(size=14),
    #legend.key=element_blank(),
    legend.title = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20),
    axis.title.y = element_text(size=30,vjust=0.5),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour='black',size=2),
    axis.ticks = element_line(size=1,color='black'),
    plot.title = element_text(size=20)) #1200x400           # Position legend in bottom right

evalsim = simlist[[4]]
phen = getPhenoMat(evalsim)[-1,]
i = 48*2 #30,40,45,60,96
pop = evalsim@simlist[[i+1]]
pop$phenotype[which(pop$phenotype==9)] = 7
pop$phenotype[which(pop$phenotype %in% c(3,4,10:14))] = 14
pop$phenotype[which(pop$phenotype==0)] = 1
pop$phenotype[which(pop$phenotype==2)] = 2
pop$phenotype[which(pop$phenotype==5)] = 3
pop$phenotype[which(pop$phenotype==6)] = 4
pop$phenotype[which(pop$phenotype==7)] = 5
pop$phenotype[which(pop$phenotype==8)] = 6
pop$phenotype[which(pop$phenotype==14)] = 7
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xaxt='n',yaxt='n',ann=FALSE,
     axes=T,cex=0.5, col=pcol[pop$phenotype], pch=19)#1000

popdat = data.frame('mean'=apply(pop,2,mean),'sd'=apply(pop,2,sd),'time'=0:96/2)
ggplot(popdat, aes(x=time, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.3,size=1) +
  geom_line(size=1.2) +
  #geom_point(size=5, shape=20) + # 21 is filled circle
  geom_vline(xintercept = c(17,20,30,48), colour="black", linetype = "longdash",size=1.5) +
  xlab("Time in h") +
  ylab("Number of bacterial cells") +
  scale_y_continuous(label=scientific_format()) +
  theme_bw(base_size = 30) +
  theme(legend.position='none',
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=20)) #12x6           # Position legend in bottom right

subcol = setdiff(colpal2,pcol)
names(subcol) = c('Acetate','CO2','D-Glucose','Glycolaldehyde','NH4+','Orthophosphate','Oxygen','Succinate','Sulfate')
trans = modelP@met_name
names(trans) = paste('EX',gsub('\\[None\\]','',modelP@met_id),sep='_')
trans = gsub('\\[e\\]','',trans)
subdat = melt(do.call(rbind,lapply(subs,function(x){apply(x,2,mean)})))
colnames(subdat) = c('Substances','time','mean')
subdat$sd = melt(do.call(rbind,lapply(subs,function(x){apply(x,2,function(x){sd(x)/sqrt(length(x))})})))[,3]
subdat$time = (subdat$time-1)/2
subdat$Substances = trans[gsub('_E','_',subdat$Substances)]
subdat = subdat[-which(subdat$Substances %in% c('Biomass','H2O','H+')),]
ggplot(subdat, aes(x=time, y=mean, colour=Substances, group=Substances)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.3,size=1) +
  geom_line(size=1.2) +
  #geom_point(size=5, shape=20) + # 21 is filled circle
  geom_vline(xintercept = c(17,20,30,48), colour="black", linetype = "longdash",size=1.5) +
  xlab("Time in h") +
  ylab("Concentration in mM") +
  scale_y_continuous(label=scientific_format()) +
  scale_color_manual(values=subcol) +
  theme_bw(base_size = 30) +
  theme(#legend.position='none',
        legend.text=element_text(size=14),
        #legend.key=element_blank(),
        legend.title = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=20)) #12x6           # Position legend in bottom right

library(NMF)
library(vegan)
library(RColorBrewer)

evalsim = simlist[[1]]
i = 48*2 #17,20,30,48
sub = matrix(evalsim@medlist[[i+1]][['EX_EC0007']],ncol=200,nrow=200)
image(sub)
aheatmap(sub/(0.01*6.25e-08),scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,
         color=colorRampPalette(c("seashell", subcol['Oxygen']))(100))#580x533

i = 48*2 #17,20,30,48
#sub = matrix(0,ncol=200,nrow=200)
sub = matrix(evalsim@medlist[[i+1]][['EX_EC0029']],ncol=200,nrow=200)
aheatmap(sub/(0.01*6.25e-08),scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,breaks=0.025,
         color=colorRampPalette(c("seashell", subcol['Acetate']))(100))#580x533
i = 48*2 #17,20,30,48
sub = matrix(evalsim@medlist[[i+1]][['EX_EC0036']],ncol=200,nrow=200)
aheatmap(sub/(0.01*6.25e-08),scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,breaks=0.025,
         color=colorRampPalette(c("seashell", subcol['Succinate']))(100))#580x533



evalsim = simlist[[1]]
cols = colpal2[c(1,3,7,8)]
evalsim = simlist[[4]]
phen = getPhenoMat(evalsim)[-1,]
phen = phen[,which(apply(phen,2,sum)!=0)]
colnames(phen) = trans[gsub('_E','_',colnames(phen))]
cmp = 'EX_EC0029'
cmp = 'EX_EC0036'
cmp = 'EX_EC0027'
cmp = 'EX_EC0007'
i = 96 #30,40,45,60,96
pop = evalsim@simlist[[i]]
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xaxt='n',yaxt='n',ann=FALSE,
     axes=T,cex=(pop$growth/max(pop$growth)), col=pcol[pop$phenotype+1], pch=19)#700x750
videoPhen = function(evalsim){
  for(i in 1:length(evalsim@simlist)){
    pop = evalsim@simlist[[i]]
    plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xaxt='n',yaxt='n',ann=FALSE,
         axes=T,cex=(pop$growth/max(pop$growth)), col=pop$phenotype+1, pch=19)#700x750
  }
}

2 = 'Glycolaldehyde production'
#3 = 'Glycolaldehyde, Acetate production'
#4 = 'Glycolaldehyde production / Acetate consumption'
5 = 'Succinate, Acetate, CO2 production'
6 = 'Acetate production / Succinate consumption'
7 = 'Succinate, Acetate consumption'
8 = 'Succinate, Acetate production'
#9 = 'Succinate, Acetate consumption'= 7


10 = 'Glycolaldehyde, Acetate production / Succinate consumption' = 6
11 = 'Succinate, Acetate production'
12 = 'Succinate, Acetate production' = 11
13 = 'Succinate, Acetate production' = 11
14 = 'Glycolaldehyde, Succinate, Acetate production' = 5
phen[,c(2,5,9)]
c(8,9,10,12,14,15)