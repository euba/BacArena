library(RNetLogo)

nlpath <- "/home/eugen/netlogo-5.0.4"
modelpath <- "/home/eugen/BacArena/Quorum_Sensing.nlogo"

NLStart("/home/eugen/netlogo-5.0.4", gui=TRUE, obj.name=NULL, nl.version=5, is3d=FALSE)
NLLoadModel(modelpath)

#re <- vector("numeric")
NLCommand("setup")
for(i in 1:10){
  NLCommand("go")
  #re[i] <- NLReport("EnergyG")
  pos <- NLGetAgentSet(c("xcor","ycor"),"turtles")
  #plot(pos)
  patch <- NLGetPatches("pcolor","patches")
  plot(patch)
}
#NLDoCommand(100,"go")
#NLGetPatches("pcolor","patches")
NLDoReport(100, "go", "EnergyG")
