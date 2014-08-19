

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
  setwd("~/BacArena")
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

setGeneric("addSubs", function(object, smax=20){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax=20){
  media = object@media
  specs = object@specs
  
  if(length(mediac)==0){
    mediac = unique(unlist(lapply(specs, function(x,ex){
      allreact <- react_id(x) #uses flags to find exchange reactions
      upt <- allreact[grep(ex, allreact)]
      return(upt)
    }, ex=ex)))
  }
  sapply(mediac, function(x, smax, n, m){
    media[[x]] <<- Substance(n, m, 0, name=x)
  }, smax=smax, n=n, m=m)
  
  eval.parent(substitute(object@media <- media))
})


setGeneric("changeSub", function(object, subname, value){standardGeneric("changeSub")})
setMethod("changeSub", "Arena", function(object, subname, value){
  if (subname %in% names(object@media)) eval.parent(substitute(object@media[subname] <- Substance(object@n, object@m, smax=value, name=subname)))
  #return(Substance(object@n, object@m, smax=value, name=subname))
  else stop("Substance does not exist in medium")
})


setGeneric("addBac", function(object, mediac={}, bacfile="data/ecore_model.R", amount, feat=data.frame("Type"=rep("Bac", length(specs)), 
                                                                       "Motility"=rep("random", length(specs)))){standardGeneric("addBac")})
setMethod("addBac", "Arena", function(object, mediac={}, bacfile="data/ecore_model.R", amount, feat=data.frame("Type"=rep("Bac", length(specs)), 
                                                                               "Motility"=rep("random", length(specs)))){
  if(amount > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  
  load(bacfile)
  mod <- model
  specs = c(object@specs, mod)
  specn=rep(amount, length(specs))
  orglist = object@orglist
  occmat = object@occmat
  media = object@media
  mediac= object@mediac
  n = object@n
  m = object@m
  
  if(length(mediac)==0){
    mediac = unique(unlist(lapply(specs, function(x,ex){
      allreact <- react_id(x) #uses flags to find exchange reactions
      upt <- allreact[grep(ex, allreact)]
      return(upt)
    }, ex=ex)))
  }
  for(i in seq_along(specs)){
    switch(as.character(feat[i,1]),
           "Bac"= {specI <- Bac(x=sample(1:n, 1), y=sample(1:m, 1), model=specs[[i]], growth=1, n=n, m=m)},
           "Organism"= {specI <- Organism(x=sample(1:n, 1), y=sample(1:m, 1), model=specs[[i]], n=n, m=m)},
           stop("Your Organism class is not defined yet."))
    exs <- findExchReact(specI@model) #sybil function looks for bounds
    upt <- uptReact(exs)
    constrain(specI, upt, lb=0) #constrain the uptake reaction to zero to define medium in next step
    constrain(specI, mediac, lb=-smax)
    for(j in 1:specn[i]){
      while(occmat[specI@x,specI@y] != 0){
        specI@x <- sample(1:n, 1)
        specI@y <- sample(1:m, 1)
      }
      occmat[specI@x,specI@y]=i
      orglist=c(orglist, specI)
    }
  }
  sapply(mediac, function(x, smax, n, m){
    #media[[x]] <<- Substance(n, m, smax, name=x)
    #media[[x]] <<- Substance(n, m, smax, name=x)
    media[[x]] <<- Substance(n, m, 0, name=x)
  }, smax=smax, n=n, m=m)
  
  eval.parent(substitute(object@orglist <- orglist))
  eval.parent(substitute(object@occmat <- occmat))
  eval.parent(substitute(object@orgn <- length(orglist)))
})


#function for getting a vector of media concentrations for a specific position

setGeneric("getmed", function(object, xp, yp){standardGeneric("getmed")})
setMethod("getmed", "Arena", function(object, xp, yp){
  return(unlist(lapply(object@media, function(x, xp, yp){return(x@diffmat[xp,yp])}, xp=xp, yp=yp)))
})





