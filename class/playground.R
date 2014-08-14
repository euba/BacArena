# this file is a playground to test how the current classes work in action
# the main goal of this file is to construct a basic framework for BacArena, which can then be merged with diffbac
# it is actually a little bit like diffbac.R, but for the current oop version of BacArena

medcon = getmed(pop,bac@x,bac@y)
constrain(bac, names(medcon), lb=-medcon)


ltest <- list()
for(i in 1:1000){
  ltest[[i]] <- sysBiolAlg(mod, algorithm = "fba")
  mod = changeBounds(mod, "EX_o2(e)", lb=runif(1, min = -100, max = 0))
  ltest[[i]] <- mod
  #ltest[[i]] <- org1
}

cl <- makeCluster(4, type = "SOCK")
clusterExport(cl, list("ltest", "optimizeProb"))
#clusterApply(cl, library(sybil))
#clusterCall(cl, optimizeLP, org1)
#clusterCall(cl, optimizeProb, org1@model)
#system.time(
  #parLapply(cl, ltest, function(x){print(optimizeProb(x))})
  test <- parLapply(cl, ltest, optimizeProb)
  #parLapply(cl, ltest, optimizeLP)
#)
stopCluster(cl)

system.time(
  lapply(ltest, optimizeProb)
)

setwd("~/BacArena")

optimizeProb(bac1@model, solverParm = list(PRESOLVE = "GLP_ON"))

fba <- sysBiolAlg(mod, algorithm = "fba")
system.time(
for(i in 1:1000){
  #lb_new = mod@lowbnd
  #lb_new[which(react_id(mod)=="EX_o2(e)")] = 0
  #sol = optimizeProb(fba, lb=lb_new)
  sol = optimizeProb(fba)
})
system.time(
for(i in 1:1000){
  #lb_new = mod@lowbnd
  #lb_new[which(react_id(mod)=="EX_o2(e)")] = 0
  #sol = optimizeProb(mod, lb=lb_new)
  sol = optimizeProb(mod)
  #print(sol@lp_obj)
})

source(file="class/Organism.R")

sf = vector()
so = vector()
sb = vector()
org1 <- Organism(x=1, y=2, model=mod, n=1, m=1)
bac1 <- Bac(x=1, y=2, model=mod, n=1, m=1)
fba <- sysBiolAlg(mod, algorithm = "fba")
for(j in 1:20){
  sf[j] = system.time(for(i in 1:1000){sol = optimizeProb(fba)})
  so[j] = system.time(for(i in 1:1000){optimizeLP(org1)})[1]
  sb[j] = system.time(for(i in 1:1000){optimizeLP(bac1)})[1]
}

boxplot(cbind(sf,so,sb))

test <- optimizeProb(fba)
optimizeProb(test)

promptSysBiolAlg()

# load libraries and other R files to have everything in place
library(snow)
library(Rcpp)
library(inline)
library(sybil)
SYBIL_SETTINGS("SOLVER", "clpAPI")
#load class definitions
source(file="cpp_source.R")
source(file="class/class_baggage.R")
source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Bac.R")
source(file="class/Organism.R")
source(file="class/Population.R")

#load ecoli core model to play around
load("data/ecore_model.R")
mod <- model

org1 = Organism(x=1, y=1, model=mod, n=1, m=1)
org1  = optimizeLP(org1)
org1@fbasol$obj

#testing constructor
bac1 = Bac(x=1, y=1, model=mod, growth=1, n=1, m=1, type="test")
#bac1 = Bac(x=1, y=1, model=mod, growth=1, fobj="ATPM")

findUpt(bac1, flag=F)

constrain(bac1, findUpt(bac1), lb=0, ub=1000)

optimizeLP(bac1)

optimizeProb(bac1@model, solverParm = list(PRESOLVE = GLP_ON))
optimizeProb(bac1@model, solverParm = list(warmUpGLPK = GLP_ON))
bac1@lpobj
optimizer(bac1@model)

checkOptSol(bac1@lpobj)

constrain(bac1, "EX_o2(e)", lb=-1000, ub=0)
optimizeLP(bac1)
bac1@lpobj

bac1@model = changeFobj(bac1@model, "ATPM")
org1 <- Organism(x=1, y=2, model=mod, n=1, m=1)
#org1@model = org1@changeObj(org1@model, "ATPM")

#bac1 = Bac(org1, growth=1) #this does not work, but would be cool if...
#testing constructor
org1 <- Organism(x=1, y=2, model=mod, n=1, m=1)
#checking the organisms functions
org1@lpobj
constrain(org1, "EX_o2(e)", lb=0, ub=0)
optimizeLP(org1)
org1@lpobj
changeFobj(org1, "ATPM")
org1@lpobj

#testing linear growth of bacteria
growLin(bac1)
bac1@growth
for(i in 1:100){
  growLin(bac1)
  print(bac1@growth)
}

#testing constructor
n=100
m=100
smax=50
diffmat = matrix(smax, nrow=n, ncol=m)
diffmat[(n/2-n/4):(n/2+n/4), (m/2-m/4):(m/2+m/4)] = 0
sub1 <- Substance(n=n, mm, smax=smax, name="test", diffmat=diffmat)
#jpeg(paste("plot", formatC(1, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
#image(sub1@diffmat, axes=FALSE, col=rev(heat.colors(200)))
#dev.off()
diffuseNaive(sub1)
diffuseNaiveR(sub1)
substrate=list()
dmat = sub1@diffmat
substrate[[1]] <- dmat
for(i in 2:100){
  diffuseNaiveR(sub1)
  #diffuseNaive(sub1)
  #dmat = sub1@diffmat
  #jpeg(paste("plot", formatC(i, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
  #image(dmat, axes=FALSE, col=rev(heat.colors(200)))
  #jpeg(paste("plot", formatC(i, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
  #image(sub1@diffmat, axes=FALSE, col=rev(heat.colors(200)))
  #dev.off()
  print(i)
}
diffusion(list(dmat))

substrate <- list()


specs3 = list(specs[[1]], specs[[2]], specs[[3]])
specs=list()
specs[[1]]=mod

system.time(Population(specs, specn=rep(100, length(specs)), n=100, m=100))
repliDie(pop)
#jpeg(paste("plot", formatC(1, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
image(pop2mat(pop))
image(pop@media$EX_for@diffmat)
#dev.off()
for(i in 2:200){
  moveRand(pop)
  #jpeg(paste("plot", formatC(i, width = 4, format = "d", flag = "0"), ".jpg"  ,sep=""), width = 800, height = 800)
  image(pop2mat(pop))
  #dev.off()
  print(i)
}


pop = Population(specs, specn=10, n=10, m=10)
for(i in 1:50){
  moveRand(pop)
  repliDie(pop)
  for(j in seq_along(pop@media)){
    diffuseNaive(pop@media[[j]])
  }
  dmat <- pop@media[["EX_glc(e)"]]@diffmat
  image(dmat)
  image(pop2mat(pop))
  print(i)
  print(sum(unlist(lapply(pop@orglist, function(x){print(x@growth)}))))
  print(pop@orgn)
}


####################################
View(cbind(mod@lowbnd,pop@orglist[[1]]@model@lowbnd))


#source(file="class/Organism.R")
#source(file="class/Bac.R")
#source(file="class/Population.R")
pop = Population(list(mod), specn=10, n=20, m=20, smax=50)

cl <- makeCluster(6, type = "SOCK")
clusterExport(cl, list("optimizeProb"))
#system.time(
for(i in 1:50){
  moveRand(pop)
  pop@media <- lapply(pop@media, function(x){
    diffuseNaiveR(x)
    return(x)})
  #dmat <- pop@media[["EX_glc(e)"]]@diffmat
  #image(dmat)
  #image(pop2mat(pop))

  fbares = parLapply(cl, pop@orglist, function(x){
    return(optimizeProb(x@model, retOptSol=F))
    })
  n <- pop@n
  m <- pop@m
  bmat <- pop2imat(pop)
  bmatn <- matrix(NA, nrow=n+2, ncol=m+2) #define environment with boundary conditions
  bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
  spes <- pop@orglist
  sprm <- rep(0, length(spes))
  for(j in seq_along(fbares)){
    spec <- spes[[j]]
    if(is.na(fbares[[j]]$obj)){
      fbares[[j]]$obj = 0
      spec@fbasol = fbares[[j]]
    }else{
      spec@fbasol = fbares[[j]]
    }
    ic = spec@x
    jc = spec@y
    upts <- findUpt(spec)
    constrain(spec, upts, lb=0) #define the medium in the next step
    mediaspec <- pop@media[upts]
    for(k in seq_along(mediaspec)){
      constrain(spec, mediaspec[[k]]@name, lb=-mediaspec[[k]]@diffmat[spec@x,spec@y])
    }
    
    ######################################################################################
    growLin(spec, 0.1) #linear growth
    ######################################################################################
    
    mediaspec <- lapply(mediaspec, function(x, bac){consume(bac, x)}, bac=spec) #account for the consumption of metbaolites
    pop@media[names(mediaspec)] <- mediaspec #update media composition in the original object
    spes[[j]] <- spec #update the species list
    ## now let them die
    if(spec@growth < 0.2){
      sprm[j] = 1
      bmat[ic, jc] <- 0
      next
    }
    ## now let them replicate
    if(spec@growth > 2){ #test if they are able to replicate (enough accumulated biomass)
      neighbours <- c(bmatn[ic,jc], 
                      bmatn[ic+1,jc], 
                      bmatn[ic+2,jc], 
                      bmatn[ic+2,jc+1],
                      bmatn[ic+2,jc+2], 
                      bmatn[ic+1,jc+2],
                      bmatn[ic,jc+2],
                      bmatn[ic,jc+1])
      pos <- which(neighbours==0)
      if(length(pos) > 1){
        pos = sample(pos, 1)
      }else{
        if(length(pos) == 0){
          repli(spec, 0, 0, bd=T)
          next
        }
      }
      switch(pos,
            {bmat[ic-1,jc-1] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic-1, jc-1)},
            {bmat[ic,jc-1] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic, jc-1)},
            {bmat[ic+1,jc-1] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic+1, jc-1)},
            {bmat[ic+1,jc] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic+1, jc)},
            {bmat[ic+1,jc+1] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic+1, jc+1)},
            {bmat[ic,jc+1] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic, jc+1)},
            {bmat[ic-1,jc+1] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic-1, jc+1)},
            {bmat[ic-1,jc] <- bmat[ic,jc]; spes[[length(spes)+1]] <- repli(spec, ic-1, jc)})
    } 
  }
  if(sum(sprm)!=0){
    pop@orglist <- spes[-which(sprm==1)]
  }else{
    pop@orglist <- spes
  }
  pop@orgn <- length(spes)
  if(length(pop@orglist)==0){
    print("All are Dead")
    break
  }
  print(length(pop@orglist))
  #save()
}
#)
stopCluster(cl)
