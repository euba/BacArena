

# Arena is the class, that contains a list of of all organisms in the Arena

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

setClass("Arena",
         representation(
           orgdat="data.frame", # data frame of individuals in the Arena
           specs="list", # list of organism types in the Arena
           media="list", # media composition of mixed organisms
           phenotypes="list", # list of unique phenotypes of the individuals
           mediac="character",
           occmat="Matrix", # occupacy matrix (showing which cells have bacs) -> sparse Matrix
           n="numeric",  # grid size
           m="numeric"  # grid size
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n, m){
  
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
  
  new("Arena", n=n, m=m, orgdat=data.frame(), specs=list(), media=list(), mediac=character(),
      phenotypes=list(), occmat=Matrix(as.integer(0), nrow=n, ncol=m, sparse=T))
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

# Add Individuals to the arena

setGeneric("addOrg", function(object, specI, amount, x=NULL, y=NULL, growth=1, ...){standardGeneric("addOrg")})
setMethod("addOrg", "Arena", function(object, specI, amount, x=NULL, y=NULL, growth=1, ...){
  if(amount+sum(object@occmat) > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  n <- object@n
  m <- object@m
  spectype <- specI@type
  newoccmat <- object@occmat
  neworgdat <- object@orgdat
  newspecs <- object@specs
  newphens <- object@phenotypes[[spectype]]
  
  newspecs[[spectype]] <- specI
  if(length(newphens)!=0){
    ptype <- checkPhen(object, specI)
    newphens <- object@phenotypes[[spectype]]
  }else{
    newphens[[1]] <- getPhenotype(specI)
    ptype=1
  }
  lastind <- nrow(object@orgdat)
  if(length(x*y)==0){
    for(i in 1:amount){
      neworgdat[lastind+i,'growth']=growth
      neworgdat[lastind+i,'type']=as.integer(which(names(newspecs)==spectype))
      neworgdat[lastind+i,'phenotype']=as.integer(ptype)
      x <- as.numeric(sample(1:n, 1))
      y <- as.numeric(sample(1:m, 1))
      while(newoccmat[x,y] != 0){
        x <- as.numeric(sample(1:n, 1))
        y <- as.numeric(sample(1:m, 1))
      }
      neworgdat[lastind+i,'x']=x
      neworgdat[lastind+i,'y']=y
      newoccmat[x,y] <- as.numeric(which(names(newspecs)==spectype))
    }
  }else{
    for(i in 1:amount){
      neworgdat[lastind+i,'growth']=growth
      neworgdat[lastind+i,'type']=as.integer(which(names(newspecs)==spectype))
      neworgdat[lastind+i,'phenotype']=as.integer(ptype)
      neworgdat[lastind+i,'x']=x[i]
      neworgdat[lastind+i,'y']=y[i]
      newoccmat[x[i],y[i]] <- as.numeric(which(names(newspecs)==spectype))
    }
  }

  eval.parent(substitute(object@occmat <- newoccmat))
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@specs <- newspecs))
  eval.parent(substitute(object@phenotypes[[spectype]] <- newphens))
  eval.parent(substitute(object@mediac <- union(object@mediac, specI@medium)))
})

# Add all substances defined by exchange reactions of the available bacs

setGeneric("addSubs", function(object, smax, mediac=object@mediac){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax, mediac=object@mediac){
  newmedia = object@media
  n = object@n
  m = object@m
  mediac = setdiff(mediac, names(arena@media))
  if(length(mediac)==0){
    print("No novel compound can be included to the medium.")
  }else{
    sapply(mediac, function(x, smax, n, m){
      newmedia[[x]] <<- Substance(n, m, smax, name=x)
    }, smax=smax, n=n, m=m)
    
    eval.parent(substitute(object@media <- c(object@media,newmedia)))
  }
})


setGeneric("changeSub", function(object, subname, value){standardGeneric("changeSub")})
setMethod("changeSub", "Arena", function(object, subname, value){
  if (subname %in% names(object@media)) eval.parent(substitute(object@media[subname] <- Substance(object@n, object@m, smax=value, name=subname)))
  #return(Substance(object@n, object@m, smax=value, name=subname))
  else stop("Substance does not exist in medium")
})

#function for getting a vector of media concentrations for a specific position

setGeneric("getmed", function(object, xp, yp){standardGeneric("getmed")})
setMethod("getmed", "Arena", function(object, xp, yp){
  return(unlist(lapply(object@media, function(x, xp, yp){return(x@diffmat[xp,yp])}, xp=xp, yp=yp)))
})

#function for checking if a phenotype is emergent

setGeneric("checkPhen", function(object, org){standardGeneric("checkPhen")})
setMethod("checkPhen", "Arena", function(object, org){
  ptype <- 0
  phenotypes <- object@phenotypes[[org@type]]
  phenspec <- getPhenotype(org)
  if(length(phenspec) != 0){
    for(i in 1:length(phenotypes)){
      inlist <- intersect(names(phenotypes[[i]]),names(phenspec))
      if(sum(phenotypes[[i]][inlist]==phenspec[inlist])==length(inlist)){
        ptype=i
        break
      }
    }
    if(ptype==0){
      ptype = length(phenotypes)+1
      phenotypes[[ptype]] <- phenspec
      object2 <- object
      object2@phenotypes[[org@type]] <- phenotypes
      eval.parent(substitute(object <- object2)) #has to be like this, otherwise there is a problem with the slot name!
      #eval.parent(substitute(object@phenotypes[[org@type]] <- phenotypes))
    }
  }
  return(ptype)
})

#main function for simulation of the whole arena

setGeneric("simulate", function(object, time, reduce=F){standardGeneric("simulate")})
setMethod("simulate", "Arena", function(object, time, reduce=F){
  simlist <- list()
  arena <- object
  for(i in 1:time){
    simlist[[i]] <- arena
    if(reduce){ #drop some items to reduce the overall size of simlist
      simlist[[i]]@specs <- list()
      simlist[[i]]@mediac <- character()
      simlist[[i]]@occmat <- Matrix()
    }
    for(j in seq_along(arena@media)){
      #diffuseNaiveR(arena@media[[j]])
      submat <- as.matrix(arena@media[[j]]@diffmat)
      diffuseNaiveCpp(submat, donut=FALSE)
      arena@media[[j]]@diffmat <- Matrix(submat, sparse=T)
    }
    j = 0
    orgl <- arena@orgdat
    print(system.time(while(j+1 <= nrow(orgl) && j+1 <= nrow(arena@orgdat)){ #time: 11.2
      j <- j+1
      org <- arena@specs[[arena@orgdat[j,'type']]]
      switch(class(org),
             "Bac"= {arena <- simBac(org, arena, j)},
             stop("Simulation function for Organism object not defined yet."))
    })) #background time: 0
    
    if(nrow(arena@orgdat)==0){
      print("All organisms died!")
      break
    } 
    cat("iter:", i, "bacs:",nrow(arena@orgdat),"\n\n")
  }
  return(simlist)
})

#show function for class Arena

removeMethod(show, "Arena")
setMethod(show, "Arena", function(object){
  print(paste('Arena of size ',object@n,'x',object@m,' with ',sum(object@occmat),
              ' organisms of ',length(object@specs),' species.',sep=''))
})

