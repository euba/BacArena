

# Arena is the class, that contains a list of of all organisms in the Arena

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

setClass("Arena",
         representation(
           orglist= "list", # list of organisms objects in the Arena
           specs=   "list", #
           orgn= "numeric", # number of organisms in orglist
           media= "list", # media composition of mixed organisms
           mediac= "vector",
           feat= "data.frame", # data frame with features of the organisms, first column should contain the class (e.g Bac, Organism)
           occmat="matrix",    # occupacy matrix (showing which cells have bacs 
           n        = "numeric",  # grid size
           m        = "numeric"  # grid size
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n, m,  ex="EX", ...){
  
  # load libraries and other R files to have everything in place
  #setwd("~/BacArena")
  library(snow) # parallel computing
  library(Rcpp)
  library(inline)
  library(sybil)
  #SYBIL_SETTINGS("SOLVER", "glpkAPI")
  SYBIL_SETTINGS("SOLVER", "clpAPI")
  source(file="cpp_source.R")
  source(file="class/class_baggage.R")
  #source(file="class/Arena.R")
  source(file="class/Substance.R")
  source(file="class/Bac.R")
  source(file="class/Organism.R")
  #source(file="class/Population.R")
  Rcpp::sourceCpp("diff.cpp")
  
  specs=list()

  new("Arena", n=n, m=m, specs=list(), orglist=list(), orgn=0, media=list(), mediac=vector(), feat=data.frame(), occmat=matrix(0, nrow=n, ncol=m))
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#setGeneric("show", function(object){standardGeneric("show")})
#setMethod("show", "Arena", function(object){
#  cat("Object of class Arena with", object@orgn, "Individuals of types", levels(object@feat[,1]))
#})


# Add all substances defined by exchange reactions of the available bacs

setGeneric("addSubs", function(object, smax=20, ex="EX"){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax=20, ex="EX"){
  newmedia = object@media
  specs = object@specs
  mediac = object@mediac
  n = object@n
  m = object@m
  
  if(length(mediac)==0){
    mediac = unique(unlist(lapply(specs, function(x,ex){
      allreact <- react_id(x) #uses flags to find exchange reactions
      upt <- allreact[grep(ex, allreact)]
      return(upt)
    }, ex=ex)))
  }
  sapply(mediac, function(x, smax, n, m){
    newmedia[[x]] <<- Substance(n, m, smax, name=x)
  }, smax=smax, n=n, m=m)
  
  eval.parent(substitute(object@media <- newmedia))
})


setGeneric("changeSub", function(object, subname, value){standardGeneric("changeSub")})
setMethod("changeSub", "Arena", function(object, subname, value){
  if (subname %in% names(object@media)) eval.parent(substitute(object@media[subname] <- Substance(object@n, object@m, smax=value, name=subname)))
  #return(Substance(object@n, object@m, smax=value, name=subname))
  else stop("Substance does not exist in medium")
})


setGeneric("addBac", function(object, mediac={}, bac=model, amount, ex="EX",feat=data.frame("Type"=rep("Bac", length(newspecs)), 
                                                                       "Motility"=rep("random", length(newspecs)))){standardGeneric("addBac")})
setMethod("addBac", "Arena", function(object, mediac={}, bac=model, amount, ex="EX", feat=data.frame("Type"=rep("Bac", length(newspecs)), 
                                                                               "Motility"=rep("random", length(newspecs)))){
  if(amount+object@orgn > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  
  mod <- bac
  newspecs = unique(c(object@specs, mod))
  specn=rep(amount, length(newspecs))
  neworglist = object@orglist
  newoccmat = object@occmat
  media = object@media
  mediac= object@mediac
  n = object@n
  m = object@m
  
  if(length(mediac)==0){
    mediac = unique(unlist(lapply(newspecs, function(x,ex){
      allreact <- react_id(x) #uses flags to find exchange reactions
      upt <- allreact[grep(ex, allreact)]
      return(upt)
    }, ex=ex)))
  }
  for(i in seq_along(newspecs)){
    switch(as.character(feat[i,1]),
           "Bac"= {specI <- Bac(x=sample(1:n, 1), y=sample(1:m, 1), model=newspecs[[i]], growth=1)},
           "Organism"= {specI <- Organism(x=sample(1:n, 1), y=sample(1:m, 1), model=newspecs[[i]])},
           stop("Your Organism class is not defined yet."))
    #exs <- findExchReact(specI@model) #sybil function looks for bounds
    #upt <- uptReact(exs)
    upt <- findUpt(specI)
    constrain(specI, upt, lb=0) #constrain the uptake reaction to zero to define medium in next step
    #constrain(specI, mediac, lb=-smax)
    for(j in 1:specn[i]){
      while(newoccmat[specI@x,specI@y] != 0){
        specI@x <- sample(1:n, 1)
        specI@y <- sample(1:m, 1)
      }
      newoccmat[specI@x,specI@y]=i
      neworglist=c(neworglist, specI)
    }
  }
  sapply(mediac, function(x, smax, n, m){
    #media[[x]] <<- Substance(n, m, smax, name=x)
    #media[[x]] <<- Substance(n, m, smax, name=x)
    media[[x]] <<- Substance(n, m, 0, name=x)
  }, smax=smax, n=n, m=m)
  
  eval.parent(substitute(object@orglist <- neworglist))
  eval.parent(substitute(object@occmat <- newoccmat))
  eval.parent(substitute(object@orgn <- length(neworglist)))
  eval.parent(substitute(object@specs <- newspecs))
})


setGeneric("addListBac", function(object, mediac={}, baclist, amount, ex="EX",feat=data.frame("Type"=rep("Bac", length(newspecs)), 
                                                                                                               "Motility"=rep("random", length(newspecs)))){standardGeneric("addListBac")})
setMethod("addListBac", "Arena", function(object, mediac={}, baclist, amount, ex="EX", feat=data.frame("Type"=rep("Bac", length(newspecs)), 
                                                                                                                        "Motility"=rep("random", length(newspecs)))){
  if(amount+object@orgn > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  
  newspecs = baclist
  specn=rep(amount, length(baclist))
  neworglist = object@orglist
  newoccmat = object@occmat
  media = object@media
  mediac= object@mediac
  n = object@n
  m = object@m
  
  if(length(mediac)==0){
    mediac = unique(unlist(lapply(newspecs, function(x,ex){
      allreact <- react_id(x) #uses flags to find exchange reactions
      upt <- allreact[grep(ex, allreact)]
      return(upt)
    }, ex=ex)))
  }
  for(i in seq_along(newspecs)){
    cat("Adding Organism type", i, "\n")
    switch(as.character(feat[i,1]),
           "Bac"= {specI <- Bac(x=sample(1:n, 1), y=sample(1:m, 1), model=newspecs[[i]], growth=1)},
           "Organism"= {specI <- Organism(x=sample(1:n, 1), y=sample(1:m, 1), model=newspecs[[i]])},
           stop("Your Organism class is not defined yet."))
    #exs <- findExchReact(specI@model) #sybil function looks for bounds
    #upt <- uptReact(exs)
    upt <- findUpt(specI)
    constrain(specI, upt, lb=0) #constrain the uptake reaction to zero to define medium in next step
    #constrain(specI, mediac, lb=-smax)
    for(j in 1:specn[i]){
      while(newoccmat[specI@x,specI@y] != 0){
        specI@x <- sample(1:n, 1)
        specI@y <- sample(1:m, 1)
      }
      newoccmat[specI@x,specI@y]=i
      neworglist=c(neworglist, specI)
    }
  }
  sapply(mediac, function(x, smax, n, m){
    #media[[x]] <<- Substance(n, m, smax, name=x)
    #media[[x]] <<- Substance(n, m, smax, name=x)
    media[[x]] <<- Substance(n, m, 0, name=x)
  }, smax=smax, n=n, m=m)
  
  eval.parent(substitute(object@orglist <- neworglist))
  eval.parent(substitute(object@occmat <- newoccmat))
  eval.parent(substitute(object@orgn <- length(neworglist)))
  eval.parent(substitute(object@specs <- newspecs))
})

#function for getting a vector of media concentrations for a specific position

setGeneric("getmed", function(object, xp, yp){standardGeneric("getmed")})
setMethod("getmed", "Arena", function(object, xp, yp){
  return(unlist(lapply(object@media, function(x, xp, yp){return(x@diffmat[xp,yp])}, xp=xp, yp=yp)))
})

setGeneric("simulate", function(object, time){standardGeneric("simulate")})
setMethod("simulate", "Arena", function(object, time){
  simlist <- list()
  arena <- object
  for(i in 1:time){
    simlist[[i]] <- arena
    for(j in seq_along(arena@media)){
      #diffuseNaiveR(arena@media[[j]])
      diffuseNaiveCpp(arena@media[[j]]@diffmat, donut=FALSE)
    }
    j = 0
    orgl <- arena@orglist
    print(system.time(while(j+1 <= length(orgl) && j+1 <= length(arena@orglist)){ #time: 21
      j<-j+1
      move(arena@orglist[[j]],arena) #time: 2
      medcon = getmed(arena,arena@orglist[[j]]@x,arena@orglist[[j]]@y) #time: 1
      constrain(arena@orglist[[j]], names(medcon), lb=-medcon) #time: 5
      optimizeLP(arena@orglist[[j]]) #time: 2
      arena@media = consume(arena@orglist[[j]],arena@media) #time: 1
      growth(arena@orglist[[j]], arena, j,lifecosts=0.1) #time: 15 -> Problem: overwriting of orglist (is too big)
    })) #background time: 0
    
    if(length(arena@orglist)==0){
      print("All bacs dead!")
      break
    } 
    #cat("iter:", i, "bacs:",length(arena@orglist),"\n\n")
  }
  return(simlist)
})

#show function for class Arena

removeMethod(show, "Arena")
setMethod(show, "Arena", function(object){
  print(paste('Arena of size ',object@n,'x',object@m,' with ',sum(object@occmat),
              ' organisms of ',length(object@specs),' species.',sep=''))
})


