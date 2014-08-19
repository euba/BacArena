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


#function for converting Arena object into a position/type data.frame of all organisms involved

setGeneric("pop2dat", function(object){standardGeneric("pop2dat")})
setMethod("pop2dat", "Arena", function(object){
  tmp <- lapply(object@orglist, function(x){return(c(x@x, x@y, x@type))})
  tmp <- t(as.data.frame(tmp))
  popdat <- data.frame("x"=as.numeric(tmp[,1]), "y"=as.numeric(tmp[,2]), "type"=as.factor(tmp[,3]), row.names=1:nrow(tmp))
  return(popdat)
})

#function for converting Arena object into a position/object data.frame of all organisms involved

setGeneric("pop2dat2", function(object){standardGeneric("pop2dat2")})
setMethod("pop2dat2", "Arena", function(object){
  tmp <- lapply(object@orglist, function(x){return(c(x@x, x@y, x))})
  tmp <- t(as.data.frame(tmp))
  popdat <- data.frame("x"=as.numeric(tmp[,1]), "y"=as.numeric(tmp[,2]), "type"=as.factor(tmp[,3]), row.names=1:nrow(tmp))
  return(popdat)
})

#function for converting Arena object into a presence/absence matrix of all organisms involved

setGeneric("pop2mat", function(object){standardGeneric("pop2mat")})
setMethod("pop2mat", "Arena", function(object){
  popdat <- pop2dat(object)
  popmat <- matrix(0, object@n, object@m)
  apply(popdat, 1, function(x, types){
    popmat[as.numeric(x[1]), as.numeric(x[2])] <<- which(types==x[3])
  }, types=levels(popdat$type))
  return(popmat)
})


#function for converting Arena object into a presence/absence matrix of all individuals involved

setGeneric("pop2imat", function(object){standardGeneric("pop2imat")})
setMethod("pop2imat", "Arena", function(object){
  popdat <- pop2dat(object)
  popmat <- matrix(0, object@n, object@m)
  for(i in 1:nrow(popdat)){
    popmat[popdat[i,]$x, popdat[i,]$y] <- i
  }
  return(popmat)
})

#function for getting a vector of media concentrations for a specific position

setGeneric("getmed", function(object, xp, yp){standardGeneric("getmed")})
setMethod("getmed", "Arena", function(object, xp, yp){
  return(unlist(lapply(object@media, function(x, xp, yp){return(x@diffmat[xp,yp])}, xp=xp, yp=yp)))
})


#function for random walk (movement) of the whole Arena through the grid space

setGeneric("moveRand", function(object){standardGeneric("moveRand")})
setMethod("moveRand", "Arena", function(object){
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
    neighbours = ifelse(is.na(neighbours),1,neighbours)
    pos <- which(neighbours==0)
    if(length(pos) >= 1){
      if(length(pos)!=1){
        pos = sample(pos, 1)
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
    }else{
      if(length(pos) == 0){
        next
      }
    }
    newpos <- which(bmat==i, arr.ind=T)
    tmp <- object@occmat[[bdat[[i]]@x,bdat[[i]]@y]]
    object@occmat[[bdat[[i]]@x,bdat[[i]]@y]] <- 0
    bdat[[i]]@x <- newpos[1]
    bdat[[i]]@y <- newpos[2]
    object@occmat[[bdat[[i]]@x,bdat[[i]]@y]] <- tmp
  }
  eval.parent(substitute(object@orglist <- bdat))
  eval.parent(substitute(object@occmat <- object@occmat))
})

#function for replication and death for the individuals in the Arena
# 
# setGeneric("repliDie", function(object, ...){standardGeneric("repliDie")})
# setMethod("repliDie", "Arena", function(object, ...){
#   object2 <- object
#   n <- object2@n
#   m <- object2@m
#   bmat <- pop2imat(object2)
#   bmatn <- matrix(NA, nrow=n+2, ncol=m+2) #define environment with boundary conditions
#   bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
#   specs <- object2@orglist
#   for(i in seq_along(object2@orglist)){
#     bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
#     spec <- specs[[i]] # isolate the specific Bac object
#     ic = spec@x
#     jc = spec@y
#     
#     ## first let them grow and eat
#     upts <- findUpt(spec)#, ...)
#     constrain(spec, upts, lb=0) #define the medium in the next step
#     mediaspec <- object2@media[upts]
#     tmp <- lapply(mediaspec, function(x){ #constrain according to media composition
#       spec@model <<- constrain(spec, x@name, lb=-x@diffmat[spec@x,spec@y])
#       })
#     optimizeLP(spec) #run fba
#     growLin(spec) #linear growth
#     mediaspec <- lapply(mediaspec, function(x, bac){consume(bac, x)}, bac=spec) #account for the consumption of metbaolites
#     object2@media[names(mediaspec)] <- mediaspec #update media composition in the original object
#     specs[[i]] <- spec #update the species list
#     
#     ## now let them die
#     if(spec@growth < 0.1){
#       specs <- specs[-i]
#       bmat[ic, jc] <- 0
#       next
#     }
#     
#     ## now let them replicate
#     if(spec@growth > 2){ #test if they are able to replicate (enough accumulated biomass)
#       neighbours <- c(bmatn[ic,jc], 
#                       bmatn[ic+1,jc], 
#                       bmatn[ic+2,jc], 
#                       bmatn[ic+2,jc+1],
#                       bmatn[ic+2,jc+2], 
#                       bmatn[ic+1,jc+2],
#                       bmatn[ic,jc+2],
#                       bmatn[ic,jc+1])
#       pos <- which(neighbours==0)
#       if(length(pos) > 1){
#         pos = sample(pos, 1)
#       }else{
#         if(length(pos) == 0){
#           repli(specs[[i]], 0, 0, bd=T)
#           next
#         }
#       }
#       switch(pos,
#             {bmat[ic-1,jc-1] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic-1, jc-1)},
#             {bmat[ic,jc-1] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic, jc-1)},
#             {bmat[ic+1,jc-1] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic+1, jc-1)},
#             {bmat[ic+1,jc] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic+1, jc)},
#             {bmat[ic+1,jc+1] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic+1, jc+1)},
#             {bmat[ic,jc+1] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic, jc+1)},
#             {bmat[ic-1,jc+1] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic-1, jc+1)},
#             {bmat[ic-1,jc] <- bmat[ic,jc]; specs[[length(specs)+1]] <- repli(specs[[i]], ic-1, jc)})
#     }
#   }
#   object2@orglist <- specs
#   object2@orgn <- length(specs)
#   eval.parent(substitute(object <- object2)) #subsitute the original object with the new one
# })
# 

# cpp movement
movement_cpp <- cxxfunction(signature(input_frame = "data.frame", seed = "integer", length_n = "integer", length_m = "integer"), body = src_movement, plugin="Rcpp")
setGeneric("movement", function(object){standardGeneric("movement")})
setMethod("movement", "Arena", function(object){
  pop <- object
  bac <- pop2dat2(pop)
  
  new_bac <- movement_cpp(bac, 123, 10, 10)
  lapply(new_bac, function(x){return(x$type)})
  
  lapply()
  eval.parent(substitute(object <- pop))
})


