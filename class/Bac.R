source(file="class/Organism.R")

# Bac is a subclass of Organism containing bacteria specific features

########################################################################################################
###################################### BAC CLASS #######################################################
########################################################################################################

setClass("Bac",
         contains="Organism",
         representation(
           #growth="numeric", # growth (biomass) of the individual
           deathrate="numeric", # factor by which growth is reduced
           duplirate="numeric", # grow cut-off for test of duplication
           speed="integer", # speed by which bacterium is moving (given by cell per iteration)
           growthlimit="numeric",
           growtype="character" # functional type for growth (linear or exponential)
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Bac <- function(model, deathrate, duplirate, speed=2, growthlimit, growtype, ...){
  new("Bac", Organism(model=model, ...), speed=as.integer(speed), deathrate=deathrate, duplirate=duplirate,
      growthlimit=growthlimit, growtype=growtype)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for letting bacteria grow by adding the calculated growthrate to the already present growth value -> linear growth
#requires as input: organism object

setGeneric("growLin", function(object, growth){standardGeneric("growLin")})
setMethod("growLin", "Bac", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- object@fbasol$obj + growth
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})

#function for letting bacteria grow by adding the calculated growthrate multiplied with the current growth plus to the already present growth value -> exp growth
#requires as input: organism object

setGeneric("growExp", function(object, growth){standardGeneric("growExp")})
setMethod("growExp", "Bac", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- (object@fbasol$obj * growth + growth)
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})


#function to get moore-neighbourhood of a bac together with its relative position

setGeneric("getHood", function(object, occmat, x, y){standardGeneric("getHood")})
setMethod("getHood", "Bac", function(object, occmat, x, y){
  occmat <- as.matrix(occmat) #dangerous!
  if(x-1==0) dx=0 else dx=1
  if(x+1>nrow(occmat)) dx2=0 else dx2=1
  if(y-1==0) dy=0 else dy=1
  if(y+1>ncol(occmat)) dy2=0 else dy2=1
  return(list(as.matrix(occmat[,(y-dy):(y+dy2)])[(x-dx):(x+dx2),], c(1+dx,1+dy)))
})


#function to check if the there is some free place in the neighbourhood

setGeneric("emptyHood", function(object, occmat, x, y){standardGeneric("emptyHood")})
setMethod("emptyHood", "Bac", function(object, occmat, x, y){
  hood <- getHood(object, occmat, x, y)
  free <- which(hood[[1]]==0, arr.ind = T)
  if(nrow(free) == 0) return(NULL)
  else {
    #print(class(free))
    #print(paste(x,' ',y))
    abs <- free[sample(nrow(free),1),]
    #abs <- free[sample(1:nrow(free),1),]
    #abs <- free[round(runif(1,min=1,max=nrow(free)),0),]
    #abs <- free[1,]
    abs[1] <- abs[1] - hood[[2]][1] + x
    abs[2] <- abs[2] - hood[[2]][2] + y
    return(abs)
  }
})


# function with the growth model of a bac (biomass growth, replication, death)

setGeneric("growth", function(object, population, j){standardGeneric("growth")})
setMethod("growth", "Bac", function(object, population, j){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$growth <- growLin(object, popvec$growth)},
         "exponential"= {popvec$growth <- growExp(object, popvec$growth)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'growth'] <- popvec$growth
  if(popvec$growth > object@duplirate){
    hood <- emptyHood(object, population@occmat, popvec$x, popvec$y)
    if(length(hood) != 0){
      daughter <- popvec
      daughter$growth <- popvec$growth/2
      daughter$x <- hood[1]
      daughter$y <- hood[2]
      popvec$growth = popvec$growth/2
      neworgdat[nrow(neworgdat)+1,] <- daughter
      neworgdat[j,] <- popvec
      eval.parent(substitute(population@occmat[daughter$x,daughter$y] <- as.numeric(daughter$type)))
    }
  }
  else if(popvec$growth < object@growthlimit){
    #print("bac dies")
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    neworgdat[j,'growth'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

# function for random movement

setGeneric("move", function(object, population, j){standardGeneric("move")})
setMethod("move", "Bac", function(object, population, j){
  popvec <- population@orgdat[j,]
  hood <- emptyHood(object, population@occmat, popvec$x, popvec$y)
  if(length(hood) != 0){
    xp = hood[1]
    yp = hood[2]
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    eval.parent(substitute(population@occmat[xp,yp] <- as.numeric(popvec$type)))
    eval.parent(substitute(population@orgdat[j,]$x <- xp))
    eval.parent(substitute(population@orgdat[j,]$y <- yp))
  }
})

#function for one iteration for Bac class

setGeneric("simBac", function(object, arena, j, sublb){standardGeneric("simBac")})
setMethod("simBac", "Bac", function(object, arena, j, sublb){
  constrain(object, object@medium, lb=-sublb[j,object@medium])
  optimizeLP(object)
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,])))
  dead <- growth(object, arena, j)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object))
  if(!dead && object@speed != 0){
    sapply(1:object@speed,function(x){
      move(object, arena, j)
      arena <<- arena})
  }
  return(arena)
})

#show function for class Bac

removeMethod(show, signature(object="Bac"))
setMethod(show, signature(object="Bac"), function(object){
  print(paste('Bacterium ',object@type,' of class Bac.',sep=''))
})

