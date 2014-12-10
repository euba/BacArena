

# Arena is the class, that contains a list of of all organisms in the Arena

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

setClass("Arena",
         representation(
           orglist= "list", #list of organisms objects in the Arena
           #specs=   "list", #
           orgn= "numeric", #number of organisms in orglist
           media= "list", #media composition of mixed organisms
           mediac= "character",
           feat= "list", #list of lists with features of the organisms such as the type, lpobj and contraints
           occmat= "matrix",    #occupacy matrix (showing which cells have bacs 
           n        = "numeric",  #grid size
           m        = "numeric"  #grid size
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n, m,  ex="EX_", ...){
  
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

  new("Arena", n=n, m=m, orglist=list(), orgn=0, media=list(), mediac=character(),
      feat=list(), occmat=matrix(0, nrow=n, ncol=m))
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#setGeneric("show", function(object){standardGeneric("show")})
#setMethod("show", "Arena", function(object){
#  cat("Object of class Arena with", object@orgn, "Individuals of types", levels(object@feat[,1]))
#})


# Add all substances defined by exchange reactions of the available bacs

setGeneric("addSubs", function(object, smax=20, ex="EX_"){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax=20, ex="EX_"){
  newmedia = object@media
  mediac = object@mediac
  n = object@n
  m = object@m

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


# Add bacteria to the arena

setGeneric("addBac", function(object, baclist, amount, ex="EX_", type='Bac'){standardGeneric("addBac")})
setMethod("addBac", "Arena", function(object, baclist, amount, ex="EX_", type='Bac'){
  if(amount+object@orgn > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  
  newspecs = baclist
  specn=rep(amount, length(baclist))
  neworglist = object@orglist
  newoccmat = object@occmat
  media = object@media
  medias = object@mediac
  n = object@n
  m = object@m
  feats = object@feat
  
  if(length(feats)==0){
    for(i in seq_along(newspecs)){
      feats[[mod_desc(newspecs[[i]])]] <- list()
      feats[[mod_desc(newspecs[[i]])]]$type <- type
    }
  }
  
  if(length(medias)==0){
    medias = as.character(unique(unlist(lapply(newspecs, function(x,ex){
      allreact <- react_id(x) #uses flags to find exchange reactions
      upt <- allreact[grep(ex, allreact)]
      return(upt)
    }, ex=ex))))
  }
  for(i in seq_along(newspecs)){
    cat("Adding Organism type", i, "\n")
    switch(feats[[i]]$type,
           "Bac"= {specI <- Bac(x=sample(1:n, 1), y=sample(1:m, 1), model=newspecs[[i]], growth=1)},
           "Organism"= {specI <- Organism(x=sample(1:n, 1), y=sample(1:m, 1), model=newspecs[[i]])},
           stop("Your Organism class is not defined yet."))
    #exs <- findExchReact(specI@model) #sybil function looks for bounds
    #upt <- uptReact(exs)
    upt <- findUpt(specI)
    #constrain(specI, upt, lb=0) #constrain the uptake reaction to zero to define medium in next step
    #constrain(specI, mediac, lb=-smax)
    feats[[i]]$lbnd <- specI@lbnd
    feats[[i]]$ubnd <- specI@ubnd
    feats[[i]]$lpobj <- specI@lpobj
    #feats[[i]]$upts <- upt
    specI@fbasol$fluxes <- specI@fbasol$fluxes[which(names(specI@lbnd) %in% upt)]
    specI@lbnd <- unname(specI@lbnd[upt])
    #specI@ubnd <- 0
    #specI@lpobj <- list()
    specI@lbnd <- NULL
    specI@ubnd <- NULL
    specI@lpobj <- NULL

    for(j in 1:specn[i]){
      while(newoccmat[specI@x,specI@y] != 0){
        specI@x <- sample(1:n, 1)
        specI@y <- sample(1:m, 1)
      }
      newoccmat[specI@x,specI@y]=i
      neworglist=c(neworglist, specI)
    }
  }
  
  eval.parent(substitute(object@mediac <- medias))
  eval.parent(substitute(object@feat <- feats))
  eval.parent(substitute(object@orglist <- neworglist))
  eval.parent(substitute(object@occmat <- newoccmat))
  eval.parent(substitute(object@orgn <- length(neworglist)))
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
    print(system.time(while(j+1 <= length(orgl) && j+1 <= length(arena@orglist)){ #time: 11.2
      j <- j+1
      orgfeat = arena@feat[[arena@orglist[[j]]@type]]
      move(arena@orglist[[j]],arena) #time: 1
      medcon = getmed(arena, arena@orglist[[j]]@x, arena@orglist[[j]]@y) #time: 0
      
      #constrain(arena@orglist[[j]], names(medcon), lb=-medcon) #time: 3
      
      orgfeat$lbnd[names(medcon)] = -medcon
      optimizeLP(arena@orglist[[j]], orgfeat$lpobj, inds=which(names(orgfeat$lbnd) %in% arena@mediac),
                 lb=orgfeat$lbnd, ub=orgfeat$ubnd) #time: 5
      arena@media = consume(arena@orglist[[j]], arena@media, fname=arena@mediac) #time: 2
      growth(arena@orglist[[j]], arena, j, lifecosts=0.1) #time: 7 -> Problem: overwriting of orglist (is too big)
    })) #background time: 0
    
    if(length(arena@orglist)==0){
      print("All organisms died!")
      break
    } 
    cat("iter:", i, "bacs:",length(arena@orglist),"\n\n")
  }
  return(simlist)
})

#show function for class Arena

removeMethod(show, "Arena")
setMethod(show, "Arena", function(object){
  print(paste('Arena of size ',object@n,'x',object@m,' with ',sum(object@occmat),
              ' organisms of ',length(unique(unlist(lapply(object@orglist,function(x){return(x@type)})))),' species.',sep=''))
})

