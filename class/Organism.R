#source(file="class/Grid.R")

# Organism is the class which contains the metabolic model and other features of the organisms in the arena

########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

setClass("Organism",
         representation(
           #x="numeric", # x position on grid
           #y="numeric", # y position on grid
           #model="modelorg", # sybil model object -> this is superbig: takes too much space and time to process
           lbnd="numeric", #lower bounds, which can change dynamically
           ubnd="numeric", #upper bounds, which can change dynamically
           type="character", # description of the organism
           medium="character", #character vector with exchange reactions of Organism
           fobj="character", # name of the objective function to be optimized
           lpobj="sysBiolAlg_fba", # sybil optimization object
           fbasol="list" # fba solution
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(model, typename=mod_desc(model), algo="fba", ex="EX_",
                     objective=model@react_id[which(model@obj_coef==1)], ...){ #the constructor requires the model, after that it is not stored anymore
  model <- changeObjFunc(model, objective)
  rxname = react_id(model)
  lpobject <- sysBiolAlg(model, algorithm=algo)
  fbasol <- optimizeProb(lpobject)
  names(fbasol$fluxes) = rxname
  
  lobnd = lowbnd(model)
  names(lobnd) = rxname
  upbnd = uppbnd(model)
  names(upbnd) = rxname
  
  new("Organism", lbnd=lobnd, ubnd=upbnd, type=typename, medium=rxname[grep(ex, rxname)],
      lpobj=lpobject, fbasol=fbasol, fobj=objective, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the objective function of the model -> might be interesting for dynamic changes in varying environments
#requires as input: organism object -> changes the model, fobj slot and lpobj of the object

#setGeneric("changeFobj", function(object, new_fobj){standardGeneric("changeFobj")})
#setMethod("changeFobj", "Organism", function(object, new_fobj){
#  eval.parent(substitute(object@fobj <- new_fobj)) #(pseudo) call by reference implementation
#  eval.parent(substitute(object@model <- changeObjFunc(object@model, new_fobj)))
#  eval.parent(substitute(object@lpobj <- sysBiolAlg(object@model, algorithm = "fba"))) #the lp object has to be updated according to the new objective
#})

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("constrain", function(object, reacts, lb){standardGeneric("constrain")})
setMethod("constrain", "Organism", function(object, reacts, lb){
  eval.parent(substitute(object@lbnd[reacts] <- lb)) #(pseudo) call by reference implementation
})

#function for computing the linear programming according to the model structure -> this can be changed for comp. speed (warmstart)
#requires as input: organism object -> changes the optimization object slot

setGeneric("optimizeLP", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd){standardGeneric("optimizeLP")})
setMethod("optimizeLP", "Organism", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd){ #this function has to be extended to contain also additional solvers
  switch(problem(lpob)@solver,
         glpkAPI=setColsBndsGLPK(problem(lpob)@oobj, 1:length(lb), #specific for GLPK!
                                    lb=lb, ub=ub),
         clpAPI=chgColLowerCLP(problem(lpob)@oobj, lb=lb))
  fbasl <- optimizeProb(lpob)
  names(fbasl$fluxes) <- names(object@lbnd)
  eval.parent(substitute(object@fbasol <- fbasl))
})

#function to account for the consumption and production of Substances

# setGeneric("consume", function(object, subs, fname=names(object@lbnd), x, y, cutoff=1e-6){standardGeneric("consume")})
# setMethod("consume", "Organism", function(object, subs, fname=names(object@lbnd), x, y, cutoff=1e-6){
#   if(object@fbasol$obj>=cutoff){
#     med = intersect(names(subs), fname)
#     flux = object@fbasol$fluxes[med]
#     flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux))
#     med = names(flux)
#     subsel = subs[med]
#     subsel = lapply(subsel, function(sub,xcoord,ycoord,flux){
#       sub@diffmat[xcoord,ycoord]<-as.matrix(sub@diffmat)[xcoord,ycoord]+flux[sub@name]
#       return(sub)
#     },xcoord=x,ycoord=y,flux=flux)
#     subs[med] = subsel
#   }
#   return(subs)
# })
setGeneric("consume", function(object, sublb, cutoff=1e-6){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, sublb, cutoff=1e-6){
  if(object@fbasol$obj>=cutoff){
    flux = object@fbasol$fluxes[object@medium]
    flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux))
    sublb[names(flux)] = sublb[names(flux)] + flux
  }
  return(sublb)
})

#function to extract the phenotype (what is consumed and what produced from the medium)

setGeneric("getPhenotype", function(object, cutoff=1e-6){standardGeneric("getPhenotype")})
setMethod("getPhenotype", "Organism", function(object, cutoff=1e-6){
  exflux=object@fbasol$fluxes[object@medium]
  exflux=ifelse(abs(exflux)<cutoff,0,1)*exflux
  exflux=ifelse(exflux<0,-1,exflux)
  exflux=ifelse(exflux>0,1,exflux)
  return(exflux[which(exflux!=0)])
})

#show function for class Organism

#removeMethod(show, signature(object="Organism"))
setMethod(show, signature(object="Organism"), function(object){
  print(paste('Organism ',object@type,' of class Organism.',sep=''))
})
