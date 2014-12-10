#source(file="class/Grid.R")

setClassUnion("numericORNULL", c("numeric","NULL"))
setClassUnion("LPORNULL", c("sysBiolAlg_fba","NULL"))

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
           lbnd="numericORNULL", #lower bounds, which can change dynamically
           ubnd="numericORNULL", #upper bounds, which can change dynamically
           type="character", # description of the organism
           #fobj="character", # name of the objective function to be optimized
           lpobj="LPORNULL", # sybil optimization object
           #lbounds="numeric", # lower bounds of the model, which can change dynamically
           fbasol="list" # fba solution
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(x, y, model, type=mod_desc(model), algorithm = "fba", ...){ #the constructor requires the model, after that it is not stored anymore
  lpobj <- sysBiolAlg(model, algorithm = "fba")
  fbasol = optimizeProb(lpobj)[-c(1,3,5,6)]
  lbnd = lowbnd(model)
  names(lbnd) = react_id(model)
  ubnd = uppbnd(model)
  names(ubnd) = react_id(model)
  new("Organism", x=x, y=y, lbnd=lbnd, ubnd=ubnd, type=mod_desc(model),
      lpobj=lpobj, fbasol=fbasol, ...)
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

setGeneric("optimizeLP", function(object, lpob=object@lpobj, inds=c(1:object@lbnd),
                                  lb=object@lbnd, ub=object@ubnd){standardGeneric("optimizeLP")})
setMethod("optimizeLP", "Organism", function(object, lpob=object@lpobj, inds=c(1:object@lbnd),
                                             lb=object@lbnd, ub=object@ubnd){ #this function has to be extended to contain also additional solvers
  switch(problem(lpob)@solver,
         glpkAPI=setColsBndsGLPK(problem(lpob)@oobj, 1:length(lb), #specific for GLPK!
                                    lb=object@lbnd, ub=ub),
         clpAPI=chgColLowerCLP(problem(lpob)@oobj, lb=lb))
  fbasl <- optimizeProb(lpob)[-c(1,3,5,6)]
  fbasl$fluxes <- fbasl$fluxes[inds]
  eval.parent(substitute(object@fbasol <- fbasl))
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
setGeneric("consume", function(object, subs, fname=names(object@lbnd)){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, subs, fname=names(object@lbnd)){
  if(object@fbasol$obj>=0){
    med = intersect(names(subs), fname)
    subs = subs[med]
    flux = object@fbasol$fluxes
    names(flux) = med
    subs = lapply(subs, function(sub,xcoord,ycoord,flux){
      sub@diffmat[xcoord,ycoord]<-sub@diffmat[xcoord,ycoord]+flux[sub@name]
      return(sub)
    },xcoord=object@x,ycoord=object@y,flux=flux)
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
