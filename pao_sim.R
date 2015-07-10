setwd("P:/BACARENA/Comparison/MatNet/P_aeruginosa/")

library(sybilSBML)
model = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/modelPOA.xml")


set.seed(5000)
# paramters from MatNet:
# if cellMass >= 2
# cellDeathThreshold = 0
model@lowbnd[grep("EX",model@react_id)]
modelP = changeBounds(model,model@react_id[grep("EX",model@react_id)],lb=-10)


bace1 = Bac(model=modelP, deathrate=0.05, duplirate=1, growthlimit=0.05, growtype="exponential",
            speed=3, type="PAO", lyse=F,chem="EX_EC0007")
arena = Arena(n=100, m=100, stir=F, tstep=0.5)
addOrg(arena, bace1, amount=22, x=c(10:20,80:90), y=rep(1,22))
addSubs(arena, smax=50, difunc="cpp", difspeed=1,
        mediac=model@react_id[grep("EX",model@react_id)][which(model@lowbnd[grep("EX",model@react_id)] < -1)])

print(system.time(evalsim <- simEnv(arena, time=50)))

library(animation)
oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
saveVideo({
  ani.options(interval = 0.5)
  evalArena(evalsim, phencol=T, plot_items=c('population','EX_EC0007'))
},video.name = "PAO_pop_phen_bio.avi", other.opts = "-b 600k")

library(animation)
oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
saveVideo({
  ani.options(interval = 0.5)
  evalArena(evalsim, phencol=T)
},video.name = "PAO_pop4.avi", other.opts = "-b 600k")
