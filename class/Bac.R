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
           growthlimit="numeric",
           growtype="character" # functional type for growth (linear or exponential)
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Bac <- function(model, deathrate, duplirate, growthlimit, growtype, ...){
  new("Bac", Organism(model=model, ...), deathrate=deathrate, duplirate=duplirate,
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
  if(length(free) == 0) return(NULL)
  else {
    abs <- free[sample(length(free[,1]),1),]
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
         "exponential"= {popvec$growth <- growLin(object, popvec$growth)},
         stop("Growth type must be either linear or exponential"))

  if(popvec$growth > object@duplirate){
    hood <- emptyHood(object, population@occmat, popvec$x, popvec$y)
    if(length(hood) != 0){
      doughter <- popvec
      doughter$growth <- popvec$growth/2
      doughter$x <- hood[1]
      doughter$y <- hood[2]
      popvec$growth = popvec$growth/2
      neworgdat[nrow(neworgdat)+1,] <- doughter
      neworgdat[j,] <- popvec
      eval.parent(substitute(population@occmat[doughter$x,doughter$y] <- as.numeric(doughter$type)))
    }
  }
  
  else if(popvec$growth < object@growthlimit){
    #print("bac dies")
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    neworgdat <- neworgdat[-j,]
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
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

setGeneric("simBac", function(object, arena, j){standardGeneric("simBac")})
setMethod("simBac", "Bac", function(object, arena, j){
  #arena@occmat <- as.matrix(arena@occmat)
  #print('move')
  move(object, arena, j) #time: 1
  #print('medcon')
  medcon = getmed(arena, arena@orgdat[j,'x'], arena@orgdat[j,'y']) #time: 0
  #print('constrain')
  constrain(object, names(medcon), lb=-medcon) #time: 3
  #print('optimizeLP')
  optimizeLP(object) #time: 5
  #print('consume')
  arena@media = consume(object, arena@media, fname=object@medium, arena@orgdat[j,'x'], arena@orgdat[j,'y']) #time problems!
  #print('growth')
  growth(object, arena, j) #time: 7 -> Problem: overwriting of orglist (is too big)
  #arena@occmat <- Matrix(arena@occmat, sparse=T)
  return(arena)
})

#show function for class Bac

removeMethod(show, signature(object="Bac"))
setMethod(show, signature(object="Bac"), function(object){
  print(paste('Bacterium ',object@type,' of class Bac.',sep=''))
})

