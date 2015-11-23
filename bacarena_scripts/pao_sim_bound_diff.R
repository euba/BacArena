setwd("~/uni/bacarena")

library(ReacTran)
library(sybil)
source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")

par(mfrow=c(1,1))


load("20151119-pao_sim.RData")
evalArena(sim, phencol=T, time = c(34, 40, 60, 96))
m<-getPhenoMat(sim, sparse=T)
plotCurves2(sim)
plotCurves2(sim, phencol = T)

heatmap(m, )

require(gplots)
heatmap.2(m, scale="none",Rowv=NA,Colv=NA,col = rev(brewer.pal(11,"RdBu")),margins=c(5,5),cexRow=0.5, cexCol=1.0,key=TRUE,keysize=1.5, trace="none")
image(m, axes = FALSE, col = grey(seq(0, 1, length = 256)))

#hm <- heatmap.2(m, scale="none", Rowv=NA, Colv=NA,
hm <- heatmap.2(t(m), scale="none",
                col = c("white", "blue", "darkblue"), ## using your colors
                breaks = c(0,0.9, 1.1,2), ## using your breaks
                #dendrogram = "none",  ## to suppress warnings
                margins=c(5,5), cexRow=0.5, cexCol=1.0, key=TRUE, keysize=1.5,
                trace="none")


load("poa_model.RData")
modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-1000)
medium = read.csv('Minimal_medium.csv')
bac = Bac(model=modelP, growtype="exponential")
setKinetics(bac, exchangeR="EX_EC0027", Km=0.01, vmax=7.56)
arena = Arena(n=20, m=20, stir=F, Lx=0.005, Ly=0.005, tstep=0.5)
#arena = Arena(n=200, m=200, stir=F, seed=i*100, Lx=0.05, Ly=0.05, tstep=0.5)
addOrg(arena, bac, amount=10)
addSubs(arena, smax=0, difspeed=1.675e-6, unit='mM')
addSubs(arena, smax=0.05, difspeed=medium$diffusion.constant, unit='mM', 
        mediac=as.character(medium$Reaction.ID))
changeSub(arena, 0, "EX_EC0007") # oxygen
changeSub(arena, 0, "EX_EC0027") # glucose
arena@media$EX_EC0007@pde="InfluxBoundDiff2d"
arena@media$EX_EC0007@boundS = 3.125e-15
arena@media$EX_EC0027@pde="InfluxBoundDiff2d"
arena@media$EX_EC0027@boundS = 3.125e-16
sim <- simEnv(arena, time=150)


evalArena(sim, time=c(10, 20, 40,60, 80, 100),phencol=T)
table(sim@simlist[[50]]$phenotype)
evalArena(sim, plot_items="EX_EC0007", time=c(10, 20, 40,60, 80, 100))
sim@medlist[[20]]$EX_EC0007
image(matrix(sim@medlist[[2]]$EX_EC0007, nrow = arena@n, ncol=arena@m))

modelP@met_name[which(modelP@met_id=='EC0211[None]')]

rxn=model@react_id[grep("EX",model@react_id)]
dic <- gsub('\\[e\\]','', sapply(paste(gsub('EX_','',rxn),'[None]',sep=''),function(x,mod){mod@met_name[which(mod@met_id == x)]},mod=modelP))
m <- getPhenoMat(sim)
colnames(m) <- dic
m <- m[-1,] # first phenotyp is NA
m <- m[,which(colSums(abs(m))>0)]
hm <- heatmap.2(m, scale="none",
                col = c("white", "darkgreen", "darkred"), ## using your colors
                breaks = c(0,0.9, 1.1,2), ## using your breaks
                #dendrogram = "none",  ## to suppress warnings
                margins=c(5,5), cexRow=0.5, cexCol=1.0, key=TRUE, keysize=1.5,
                trace="none")

selPheno(sim, time=0, type="modelPOA", reduce=T)





# EC0007 oxygen
# EC0027 glucose
# EC0957 bh4
# EC0009 Pi
# EC0048 Sulfate
# EC0001 h2o
# EC9324 biomass
# EC0064 H+
# EC0029 acetate
# EC0036 succ
# EC0011 co2

evalArena(sim, phencol=T)



simlist = list()
for(i in 1:1){
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