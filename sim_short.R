# load libraries and other R files to have everything in place
setwd("~/BacArena")
library(Rcpp)
library(sybil)
library(glpkAPI) 
#library(microbenchmark)
#library(ggplot2)
library(compiler) # byte code 
SYBIL_SETTINGS("SOLVER", "glpkAPI")
source(file="cpp_source.R")
source(file="class/class_baggage.R")
source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Bac.R")
source(file="class/Organism.R")
source(file="class/Eval.R")
Rcpp::sourceCpp("cpp/diff.cpp")
#Rcpp::sourceCpp("cpp/addBac.cpp")

set.seed(5000)
load("data/ecore_model.R")
ecore = model
ecore = changeBounds(ecore, c('EX_ac(e)','EX_akg(e)','EX_etoh(e)','EX_for(e)',
                      'EX_fum(e)','EX_lac_D(e)','EX_pyr(e)','EX_succ(e)'), -20)

bace = Bac(model=ecore, deathrate=0.4, duplirate=1.5, growthlimit=0.05, growtype="exponential",
           speed=2, budge=F, lyse=F, chem="EX_glc(e)")
arena = Arena(n=100, m=100, tstep=0.5)
addOrg(arena, bace, amount=1, x=50, y=50)
#addOrg(arena, bace, amount=500)
addSubs(arena, smax=30)
format(object.size(arena), units='b')

arena@media[["EX_o2(e)"]]@diffmat <- Matrix(0, nrow=100, ncol=100)
#arena@media[["EX_glc(e)"]]@diffmat <- Matrix(0, nrow=100, ncol=100)
#arena@media[["EX_glc(e)"]]@diffmat[50,50] <- 5000
print(system.time(simlist <- simulate(arena, time=30)))
format(object.size(simlist), units='Mb')

# simlist <- list()
# simlist[[1]] <- simulate(arena, time=1)[[2]]
# for(i in 2:50){
#   print(i)
#   simlist[[i-1]]@media[["EX_o2(e)"]]@diffmat <- Matrix(1000, nrow=100, ncol=100)
#   simlist[[i]] <- simulate(simlist[[i-1]], time=1)[[2]]
# }

for(i in 1:length(simlist)){
  print(i)
  popana <- as.matrix(simlist[[i]]@occmat)
  #dat <- simlist[[i]]@orgdat
  #for(j in 1:nrow(dat)){
  #  popana[dat[j,]$x,dat[j,]$y] <- dat[j,]$phenotype + 1
  #}
  par(mfrow = c(2, 2),bg = 'white')
  plot(simlist[[i]]@orgdat[,c('x','y')],xlim=c(0,100),ylim=c(0,100),xlab='',ylab='',
       pch=20,axes=FALSE,cex=0.4,main='Population', col=simlist[[i]]@orgdat$phenotype+1)
  #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  #image(popana, axes=F, main='Population')
  image(as.matrix(simlist[[i]]@media[["EX_glc(e)"]]@diffmat), axes=F,main='Glucose')
  image(as.matrix(simlist[[i]]@media[["EX_o2(e)"]]@diffmat), axes=F,main='Oxygen')
  image(as.matrix(simlist[[i]]@media[["EX_ac(e)"]]@diffmat), axes=F,main='Acetate')
}
