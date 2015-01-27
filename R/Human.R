source(file="Organism.R")

# Human is a subclass of Organism containing human specific features

########################################################################################################
###################################### HUMAN CLASS #####################################################
########################################################################################################

setClass("Human",
         contains="Organism",
         representation(
           objective="character" # name of the objective function, which can be changed
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Human <- function(model, deathrate, duplirate, growthlimit, growtype,
                  objective=model@react_id[which(model@obj_coef==1)], ...){
  model <- changeObjFunc(model, objective)
  new("Human", Organism(model=model, deathrate=deathrate, duplirate=duplirate, growtype=growtype,
                      growthlimit=growthlimit, ...), objective=objective)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("objective", function(object){standardGeneric("objective")})
setMethod("objective", "Human", function(object){return(object@objective)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the objective function of the model -> might be interesting for dynamic changes in varying environments
#requires as input: organism object -> changes the model, fobj slot and lpobj of the object

setGeneric("changeFobj", function(object, new_fobj, model, alg="fba"){standardGeneric("changeFobj")})
setMethod("changeFobj", "Human", function(object, new_fobj, model, alg="fba"){
  eval.parent(substitute(object@objective <- new_fobj)) #(pseudo) call by reference implementation
  model <- changeObjFunc(object@model, new_fobj)
  eval.parent(substitute(object@lpobj <- sysBiolAlg(model, algorithm=alg))) #the lp object has to be updated according to the new objective
})

# function with the growth model of a human cell (biomass growth, replication, death)

setGeneric("cellgrowth", function(object, population, j){standardGeneric("cellgrowth")})
setMethod("cellgrowth", "Human", function(object, population, j){
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
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    neworgdat[j,'growth'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

#function for one iteration for Human class

setGeneric("simHum", function(object, arena, j, sublb){standardGeneric("simHum")})
setMethod("simHum", "Human", function(object, arena, j, sublb){
  lobnd <- constrain(object, object@medium, lb=-sublb[j,object@medium],
                     dryweight=arena@orgdat[j,"growth"], time=arena@tstep)
  optimizeLP(object, lb=lobnd)
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,])))
  dead <- cellgrowth(object, arena, j)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object))
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, names(arena@media), sublb[j,])))
  }
  return(arena)
})

#show function for class Human

setMethod(show, signature(object="Human"), function(object){
  print(paste('Cell culture ',object@type,' of class Human.',sep=''))
})

