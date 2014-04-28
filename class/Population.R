source(file="class/Arena.R")

# Population is the class, that contains a list of of all organisms in the Arena

########################################################################################################
###################################### Population CLASS ################################################
########################################################################################################

setClass("Population",
         contains="Arena",
         representation(
           orglist= "list", # list of organisms objects in the Arena
           orgn= "numeric", # number of organisms in orglist
           media= "character", # media composition of mixed organisms
           feat= "data.frame" # data frame with features of the organisms, first column should contain the class (e.g Bac, Organism)
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Population <- function(specs, specn, n, m, media={}, feat=data.frame("Type"=rep("Bac", length(specs)), 
                                                     "Motility"=rep("random", length(specs))), ...){
  if(sum(specn) > n*m){
    stop("More individual than space on the grid")
  }
  orglist = list()
  pamat = matrix(0, nrow=n, ncol=m)
  for(i in seq_along(specs)){
    switch(feat[i,1],
           "Bac"= {specI <- Bac(x=sample(1:n, 1), y=sample(1:m, 1), model=specs[[i]], growth=1, n=n, m=m)},
           "Organism"= {specI <- Organism(x=sample(1:n, 1), y=sample(1:m, 1), model=specs[[i]], n=n, m=m)},
           stop("Your Organism class is not defined yet."))
    for(j in 1:specn[i]){
      while(pamat[specI@x,specI@y] != 0){
        specI@x <- sample(1:n, 1)
        specI@y <- sample(1:m, 1)
      }
      pamat[specI@x,specI@y]=1
      orglist=c(orglist, specI)
    }
  }
  media <- unique(unlist(lapply(orglist, function(x){findUpt(x, ...)}, ...))) # defining media as the union of Exchange reactions 
                                                                              # (can also be changed by altering the flag argument)
  new("Population", Arena(n=n, m=m), orglist=orglist, orgn=length(orglist), media=media, feat=feat, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#setGeneric("show", function(object){standardGeneric("show")})
#setMethod("show", "Population", function(object){
#  cat("Object of class Population with", object@orgn, "Individuals of types", levels(object@feat[,1]))
#})

#function for converting Population object into a position/type data.frame of all organisms involved

setGeneric("pop2dat", function(object){standardGeneric("pop2dat")})
setMethod("pop2dat", "Population", function(object){
  tmp <- lapply(object@orglist, function(x){return(c(x@x, x@y, x@type))})
  tmp <- t(as.data.frame(tmp))
  popdat <- data.frame("x"=as.numeric(tmp[,1]), "y"=as.numeric(tmp[,2]), "type"=as.factor(tmp[,3]), row.names=1:nrow(tmp))
  return(popdat)
})

#function for converting Population object into a presence/absence matrix of all organisms involved

setGeneric("pop2mat", function(object){standardGeneric("pop2mat")})
setMethod("pop2mat", "Population", function(object){
  popdat <- pop2dat(object)
  popmat <- matrix(0, object@n, object@m)
  apply(popdat, 1, function(x, types){
    popmat[as.numeric(x[1]), as.numeric(x[2])] <<- which(types==x[3])
  }, types=levels(popdat$type))
  return(popmat)
})

#function for converting Population object into a presence/absence matrix of all individuals involved

setGeneric("pop2imat", function(object){standardGeneric("pop2imat")})
setMethod("pop2imat", "Population", function(object){
  popdat <- pop2dat(object)
  popmat <- matrix(0, object@n, object@m)
  for(i in 1:nrow(popdat)){
    popmat[popdat[i,]$x, popdat[i,]$y] <- i
  }
  return(popmat)
})

#function for random walk (movement) of the whole population through the grid space

setGeneric("moveRand", function(object){standardGeneric("moveRand")})
setMethod("moveRand", "Population", function(object){
  n <- object@n
  m <- object@m
  bmat <- pop2imat(object)
  bmatn <- matrix(NA, nrow=n+2, ncol=m+2) #define environment with boundary conditions
  bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
  bdat <- object@orglist
  for(i in seq_along(bdat)){
    bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
    ic = bdat[[i]]@x
    jc = bdat[[i]]@y
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
        next
      }
    }
    switch(pos,
          {bmat[ic-1,jc-1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic,jc-1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic+1,jc-1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic+1,jc] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic+1,jc+1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic,jc+1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic-1,jc+1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
          {bmat[ic-1,jc] <- bmat[ic,jc]; bmat[ic,jc] <- 0})
    newpos <- which(bmat==i, arr.ind=T)
    bdat[[i]]@x <- newpos[1]
    bdat[[i]]@y <- newpos[2]
  }
  eval.parent(substitute(object@orglist <- bdat))
})

