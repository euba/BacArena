source(file="class/Arena.R")

# Organism is the class which contains the metabolic model and other features of the organisms in the arena

########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

setClass("Organism",
         contains="Arena",
         representation(
           x="numeric", # x position on grid
           y="numeric", # y position on grid
           model="modelorg", # sybil model object
           type="character", # description of the organism
           fobj="character", # name of the objective fundtion to be optimized
           lpobj="optsol" # sybil optimization object
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(x, y, model, type=mod_desc(model), fobj={}, algorithm = "fba", n, m, ...){
  if(length(fobj) != 0){ # test if the objective was changed, or still default value
    model <- changeObj(model, fobj)
  }
  lpobj <- optimizeProb(model, ...)
  new("Organism", Arena(n=n, m=m), x=x, y=y, model=model, type=mod_desc(model), fobj=model@react_id[which(model@obj_coef==1)], lpobj=lpobj, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("constrain", function(object, reacts, lb, ub){standardGeneric("constrain")})
setMethod("constrain", "Organism", function(object, reacts, lb, ub){
  eval.parent(substitute(object@model <- changeBounds(object@model, reacts, lb=lb, ub=ub))) #(pseudo) call by reference implementation
})

#function for computing the linear programming according to the model structure -> this can be changed for comp. speed (warmstart)
#requires as input: organism object -> changes the optimization object slot

setGeneric("optimizeLP", function(object){standardGeneric("optimizeLP")})
setMethod("optimizeLP", "Organism", function(object){
  eval.parent(substitute(object@lpobj <- optimizeProb(object@model))) #(pseudo) call by reference implementation
})

#function for changing the objective function of the model -> might be interesting for dynamic changes in varying environments
#requires as input: organism object -> changes the model, fobj slot and lpobj of the object

setGeneric("changeFobj", function(object, new_fobj){standardGeneric("changeFobj")})
setMethod("changeFobj", "Organism", function(object, new_fobj){
  eval.parent(substitute(object@fobj <- new_fobj)) #(pseudo) call by reference implementation
  eval.parent(substitute(object@model <- changeObj(object@model, new_fobj)))
  eval.parent(substitute(object@lpobj <- optimizeProb(object@model))) #the lp object has to be updated according to the new objective
})

#function for finding uptake reactions of a model for getting the media conditions

setGeneric("findUpt", function(object, flag=F, ex="EX"){standardGeneric("findUpt")})
setMethod("findUpt", "Organism", function(object, flag=F, ex="EX"){
  if(flag){
    ex <- findExchReact(object@model)
    upt <- uptReact(ex)
  }else{
    allreact <- react_id(object@model)
    upt <- allreact[grep(ex, allreact)]
  }
  return(upt)
})

#function to account for the consumption and production of Substances

setGeneric("consume", function(object, sub){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, sub){
  flux <- as.numeric(object@lpobj@fluxdist@fluxes[which(react_id(object@model) == sub@name),])
  sub@diffmat[object@x, object@y] <- sub@diffmat[object@x, object@y] + flux
  return(sub)
})

