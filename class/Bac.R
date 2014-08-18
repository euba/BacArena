source(file="class/Organism.R")

# Bac is a subclass of Organism containing bacteria specific features

########################################################################################################
###################################### BAC CLASS #######################################################
########################################################################################################

setClass("Bac",
         contains="Organism",
         representation(
           growth="numeric" # growth (biomass) of the individual
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Bac <- function(x, y, model, growth=1, ...){
  new("Bac", Organism(x=x, y=y, model=model, ...), growth=growth)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

# setMethod("show", "Bac",
#           function(object){
#             cat("Object of class", class(object), "with attributes:\n")
#             print(names(attributes(object)))
#           }
# )

#function for letting bacteria grow by adding the calculated growthrate to the already present growth value -> linear growth
#requires as input: organism object

setGeneric("growLin", function(object, dfactor){standardGeneric("growLin")})
setMethod("growLin", "Bac", function(object, dfactor){
  grow_accum <- object@fbasol$obj + object@growth - dfactor
  return(grow_accum)
})

#function for letting bacteria grow by adding the calculated growthrate multiplied with the current growth plus to the already present growth value -> exp growth
#requires as input: organism object

setGeneric("growExp", function(object, dfactor){standardGeneric("growExp")})
setMethod("growExp", "Bac", function(object, dfactor){
  grow_accum <- (object@fbasol$obj * object@growth + object@growth) - dfactor
  return(grow_accum)
})


#function to get moore-neighbourhood of a bac together with its relative position

setGeneric("getHood", function(object, population){standardGeneric("getHood")})
setMethod("getHood", "Bac", function(object, population){
  x=object@x
  y=object@y
  if(x-1==0) dx=0 else dx=1
  if(x+1>pop@n) dx2=0 else dx2=1
  if(y-1==0) dy=0 else dy=1
  if(y+1>pop@m) dy2=0 else dy2=1
  return(list(as.matrix(population@occmat[,(y-dy):(y+dy2)])[(x-dx):(x+dx2),], c(1+dx,1+dy)))
})


#function to check if the there is some free place in the neighbourhood

setGeneric("emptyHood", function(object, population){standardGeneric("emptyHood")})
setMethod("emptyHood", "Bac", function(object, population){
  hood <- getHood(object, population)
  free <- which(hood[[1]]==0, arr.ind = T)
  if(length(free) == 0) return(NULL)
  else {
    abs <- free[sample(length(free[,1]),1),]
    abs[1] <- abs[1] - hood[[2]][1] + object@x
    abs[2] <- abs[2] - hood[[2]][2] + object@y
    return(abs)
  }
})


# function with the growth model of a bac (biomass growth, replication, death)

setGeneric("growth", function(object, population, j, exp=T, lifecosts=0.1, repli=2){standardGeneric("growth")})
setMethod("growth", "Bac", function(object, population, j, exp=T, lifecosts=0.1, repli=2){
  newpoplist <- population@orglist
  if(exp) newpoplist[[j]]@growth <-growExp(newpoplist[[j]], 0.1) 
  else newpoplist[[j]]@growth <- growLin(newpoplist[[j]], 0.1)
  if(newpoplist[[j]]@growth > 2){
    hood <- emptyHood(newpoplist[[j]], population)
    if(length(hood) != 0){
      doughter <- newpoplist[[j]]
      newg <- newpoplist[[j]]@growth/2.0
      doughter@growth <- newg
      doughter@x <- hood[1]
      doughter@y <- hood[2]
      newpoplist[[j]]@growth = newg
      newpoplist <- c(newpoplist, doughter)
      type <- population@occmat[newpoplist[[j]]@x,newpoplist[[j]]@y]
      eval.parent(substitute(population@occmat[doughter@x,doughter@y] <- type))
      print(paste("bac", newpoplist[[j]]@x, newpoplist[[j]]@y, " replicates:", doughter@x, doughter@y))
      print(population@occmat)
    }
    eval.parent(substitute(population@orglist <- newpoplist))
    #eval.parent(substitute(population@orglist[[j]]@growth <- newg))
  }
  else if(population@orglist[[j]]@growth < 0.1){
    print("bac dies")
    newpoplist <- population@orglist[-j]
    eval.parent(substitute(population@orglist <- newpoplist))
    eval.parent(substitute(population@occmat[newpoplist[[j]]@x,newpoplist[[j]]@y] <- 0))
  }
})


setGeneric("move", function(object, population){standardGeneric("move")})
setMethod("move", "Bac", function(object, population){
  hood <- emptyHood(object, population)
  if(length(hood) != 0){
    xp = hood[1]
    yp = hood[2]
    type <- population@occmat[x,y]
    eval.parrent(substitute(population@occmat[x,y] <- 0))
    eval.parrent(substitute(population@occmat[xp,yp] <- type))
    eval.parent(substitute(object@x <- xp))
    eval.parent(substitute(object@y <- yp))
  }
})

