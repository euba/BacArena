#source(file="class/Grid.R")

# Organism is the class which contains the metabolic model and other features of the organisms in the arena

########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

setClass("Organism",
         #contains="Grid",
         representation(
           x="numeric", # x position on grid
           y="numeric", # y position on grid
           #model="modelorg", # sybil model object -> this is superbig: takes too much space and time to process
           lbnd="numeric", #lower bounds, which can change dynamically
           ubnd="numeric", #upper bounds, which can change dynamically
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

Organism <- function(x, y, model, type=mod_desc(model), fobj={}, algorithm = "fba", ...){ #the constructor requires the model, after that it is not stored anymore
  if(length(fobj) != 0){ # test if the objective was changed, or still default value
    model <- changeObj(model, fobj)
  }
  lpobj <- sysBiolAlg(model, algorithm = "fba")
  fbasol = optimizeProb(lpobj)
  lbnd = lowbnd(model)
  names(lbnd) = react_id(model)
  ubnd = uppbnd(model)
  names(ubnd) = react_id(model)
  new("Organism", x=x, y=y, lbnd=lbnd, ubnd=ubnd, type=mod_desc(model),
      fobj=model@react_id[which(model@obj_coef==1)], lpobj=lpobj,
      fbasol=fbasol, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("constrain", function(object, reacts, lb){standardGeneric("constrain")})
setMethod("constrain", "Organism", function(object, reacts, lb){
  #inds = which(reacts %in% names(object@lbnd))
  eval.parent(substitute(object@lbnd[reacts] <- lb)) #(pseudo) call by reference implementation
})


setGeneric("getExch", function(object, exn){standardGeneric("getExch")}) ###
setMethod("getExch", "Organism", function(object, exn='EX_'){
  a <- object@fbasol$fluxes
  rnam <- names(object@lbnd)
  names(a) <- rnam
  #print(a[which(a, <0)])
  ex <- rnam[grep(exn, rnam)]
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
         glpkAPI=setColsBndsGLPK(problem(object@lpobj)@oobj, 1:length(object@lbnd), #specific for GLPK!
                                    lb=object@lbnd, ub=object@ubnd),
         clpAPI=chgColLowerCLP(problem(object@lpobj)@oobj, lb=object@lbnd))
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
# -> this function hast to be implemented using the Bac constructor with another objective function
#setGeneric("changeFobj", function(object, new_fobj){standardGeneric("changeFobj")})
#setMethod("changeFobj", "Organism", function(object, new_fobj){
#  eval.parent(substitute(object@fobj <- new_fobj)) #(pseudo) call by reference implementation
#  eval.parent(substitute(object@model <- changeObj(object@model, new_fobj)))
#  eval.parent(substitute(object@lpobj <- sysBiolAlg(mod, algorithm = "fba"))) #the lp object has to be updated according to the new objective
#})

#function for finding uptake reactions of a model for getting the media conditions

setGeneric("findUpt", function(object, ex="EX_"){standardGeneric("findUpt")})
setMethod("findUpt", "Organism", function(object, ex="EX_"){
  #if(flag){
  #  upt <- uptReact(findExchReact(object@model)) #sybil function looks for bounds
  #}else{
  #  upt <- react_id(object@model)[grep(ex, allreact)] #uses flags to find exchange reactions
  #}
  return(names(object@lbnd)[grep(ex, names(object@lbnd))]) #uses flags to find exchange reactions
})

#function to account for the consumption and production of Substances

setGeneric("consume", function(object, subs){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, subs){
  if(object@fbasol$obj>=0){
    subs = subs[intersect(names(subs), names(object@lbnd))]
    subs = lapply(subs, function(sub,org){
      sub@diffmat[org@x,org@y]<-sub@diffmat[org@x,org@y]+org@fbasol$fluxes[which(names(org@lbnd)==sub@name)]
      return(sub)
    },org=object)
    #sub@diffmat[object@x, object@y]  <- sub@diffmat[object@x,object@y] + object@fbasol$fluxes[which(react_id(object@model) == sub@name)]
    #return(sub)
  }
  return(subs)
})

#show function for class Organism

#removeMethod(show, signature(object="Organism"))
setMethod(show, signature(object="Organism"), function(object){
  print(paste('Organism ',object@type,' of class Organism on position x: ',object@x,' and y: ',object@y,'.',sep=''))
})
