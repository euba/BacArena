source(file="class/Grid.R")

# Organism is the class which contains the metabolic model and other features of the organisms in the arena

########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

setClass("Organism",
         contains="Grid",
         representation(
           x="numeric", # x position on grid
           y="numeric", # y position on grid
           model="modelorg", # sybil model object
           type="character", # description of the organism
           fobj="character", # name of the objective function to be optimized
           lpobj="sysBiolAlg_fba", # sybil optimization object
           #lbounds="numeric", # lower bounds of the model, which can change dynamically
           fbasol="list" # fba solution
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(x, y, model, type=mod_desc(model), fobj={}, algorithm = "fba", n, m, ...){
  if(length(fobj) != 0){ # test if the objective was changed, or still default value
    model <- changeObj(model, fobj)
  }
  lpobj <- sysBiolAlg(mod, algorithm = "fba")
  fbasol = optimizeProb(lpobj)
  #fbasol = list()
  new("Organism", n=n, m=m, x=x, y=y, model=model, type=mod_desc(model),
      fobj=model@react_id[which(model@obj_coef==1)], lpobj=lpobj, #lbounds=mod@lowbnd,
      fbasol=fbasol, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("constrain", function(object, reacts, lb){standardGeneric("constrain")})
setMethod("constrain", "Organism", function(object, reacts, lb){
  eval.parent(substitute(object@model <- changeBounds(object@model , reacts, lb))) #(pseudo) call by reference implementation
})


setGeneric("getExch", function(object){standardGeneric("getExch")})
setMethod("getExch", "Organism", function(object){
  a <- object@fbasol$fluxes
  names(a) <- react_id(object@model)
  #print(a[which(a, <0)])
  ex <- findExchReact(object@model)
  #ex_nr <- met_pos(ex) sybil bug ....
  ex_nr <- seq(17,30)
  b <- a[ex_nr]
  print(c("input", round(b[which(b <0)],1) ))
  print(c("output",round(b[which(b >0)],1) ))
})


#function for computing the linear programming according to the model structure -> this can be changed for comp. speed (warmstart)
#requires as input: organism object -> changes the optimization object slot

setGeneric("optimizeLP", function(object){standardGeneric("optimizeLP")})
setMethod("optimizeLP", "Organism", function(object){ #this function has to be extended to contain also additional solvers
  switch(problem(object@lpobj)@solver,
         glpkAPI=setColsBndsGLPK(problem(object@lpobj)@oobj, 1:object@model@react_num, #specific for GLPK!
                                    lb=object@model@lowbnd, ub=object@model@uppbnd),
         clpAPI=chgColLowerCLP(problem(object@lpobj)@oobj, lb=object@model@lowbnd))
  eval.parent(substitute(object@fbasol <- optimizeProb(object@lpobj)))
  #fbasolnew <- optimizeProb(object@lpobj)
  #fbasolnew$fluxes <- round(fbasolnew$fluxes,2)
  #fbasolnew$obj <- round(fbasolnew$obj,2)
  #eval.parent(substitute(object@fbasol <- fbasolnew))
  
  #eval.parent(substitute(object@fbasol <- optimizeProb(sysBiolAlg(changeBounds(object@model, object@model@react_id,
  #                                                                             lb=object@model@lowbnd),algorithm="fba"))))
  #eval.parent(substitute(object@fbasol <- solfba)) #(pseudo) call by reference implementation
})

#function for changing the objective function of the model -> might be interesting for dynamic changes in varying environments
#requires as input: organism object -> changes the model, fobj slot and lpobj of the object

setGeneric("changeFobj", function(object, new_fobj){standardGeneric("changeFobj")})
setMethod("changeFobj", "Organism", function(object, new_fobj){
  eval.parent(substitute(object@fobj <- new_fobj)) #(pseudo) call by reference implementation
  eval.parent(substitute(object@model <- changeObj(object@model, new_fobj)))
  eval.parent(substitute(object@lpobj <- sysBiolAlg(mod, algorithm = "fba"))) #the lp object has to be updated according to the new objective
})

#function for finding uptake reactions of a model for getting the media conditions

setGeneric("findUpt", function(object, ex="EX"){standardGeneric("findUpt")})
setMethod("findUpt", "Organism", function(object, ex="EX"){
  #if(flag){
  #  upt <- uptReact(findExchReact(object@model)) #sybil function looks for bounds
  #}else{
  #  upt <- react_id(object@model)[grep(ex, allreact)] #uses flags to find exchange reactions
  #}
  return(react_id(object@model)[grep(ex, allreact)]) #uses flags to find exchange reactions
})

#function to account for the consumption and production of Substances

setGeneric("consume", function(object, subs){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, subs){
  if(object@fbasol$obj>=0){
    subs = lapply(subs, function(sub,org){
      sub@diffmat[org@x,org@y]<-sub@diffmat[org@x,org@y]+org@fbasol$fluxes[which(react_id(org@model)==sub@name)]
      return(sub)
    },org=object)
    #sub@diffmat[object@x, object@y]  <- sub@diffmat[object@x,object@y] + object@fbasol$fluxes[which(react_id(object@model) == sub@name)]
    #return(sub)
  }
  return(subs)
})

