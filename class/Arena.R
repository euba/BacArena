source(file="class/Grid.R")

# Arena is the class, that contains a list of of all organisms in the Arena

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

setClass("Arena",
         contains="Grid",
         representation(
           orglist= "list", # list of organisms objects in the Arena
           orgn= "numeric", # number of organisms in orglist
           media= "list", # media composition of mixed organisms
           feat= "data.frame", # data frame with features of the organisms, first column should contain the class (e.g Bac, Organism)
           occmat="matrix"    # occupacy matrix (showing which cells have bacs 
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(specs, specn, n, m, mediac={}, feat=data.frame("Type"=rep("Bac", length(specs)), 
                                            "Motility"=rep("random", length(specs))), smax=20, ex="EX", ...){
  if(sum(specn) > n*m){
    stop("More individuals than space on the grid")
  }
  orglist = list()
  occmat = matrix(0, nrow=n, ncol=m) #presence absence matrix with Organism position
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
   media = list()
   sapply(mediac, function(x, smax, n, m){
     #media[[x]] <<- Substance(n, m, smax, name=x)
     #media[[x]] <<- Substance(n, m, smax, name=x)
     media[[x]] <<- Substance(n, m, 0, name=x)
     }, smax=smax, n=n, m=m)
  new("Arena", Grid(n=n, m=m), orglist=orglist, orgn=length(orglist), media=media, feat=feat, occmat=occmat,...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#setGeneric("show", function(object){standardGeneric("show")})
#setMethod("show", "Arena", function(object){
#  cat("Object of class Arena with", object@orgn, "Individuals of types", levels(object@feat[,1]))
#})

setGeneric("addSub", function(object, subname, value){standardGeneric("addSub")})
setMethod("addSub", "Arena", function(object, subname, value){
  if (subname %in% names(object@media)) eval.parent(substitute(object@media[subname] <- Substance(object@n, object@m, smax=value, name=subname)))
    #return(Substance(object@n, object@m, smax=value, name=subname))
  else stop("Substance does not exist in medium")
})


#function for getting a vector of media concentrations for a specific position

setGeneric("getmed", function(object, xp, yp){standardGeneric("getmed")})
setMethod("getmed", "Arena", function(object, xp, yp){
  return(unlist(lapply(object@media, function(x, xp, yp){return(x@diffmat[xp,yp])}, xp=xp, yp=yp)))
})





