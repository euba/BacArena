model = readSBMLmod("P:/BACARENA/Comparison/MatNet/P_aeruginosa/POA_test.xml")


set.seed(5000)
# paramters from MatNet:
# if cellMass >= 2
# cellDeathThreshold = 0

bace1 = Bac(model=model, deathrate=0, duplirate=2, growthlimit=0.05, growtype="exponential",
            speed=1, type="PAO", lyse=T)
arena = Arena(n=30, m=30, stir=F, tstep=1/12)
addOrg(arena, bace1, amount=100, x=1:100, y=rep(1,100))
addSubs(arena, smax=20, difunc="cpp", difspeed=1)

print(system.time(evalsim <- simEnv(arena, time=200)))


library(animation)
oopts = ani.options(ffmpeg = "C:/ffmpeg.exe")
saveVideo({
  ani.options(interval = 0.5)
  evalArena(evalsim, phencol=T, plot_items=c('population','EX_EC0007'))
},video.name = "PAO_pop_phen_bio.avi", other.opts = "-b 600k")
