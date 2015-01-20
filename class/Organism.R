#source(file="class/Grid.R")

# Organism is the class which contains the metabolic model and other features of the organisms in the arena

########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

setClass("Organism",
         representation(
           lbnd="numeric", #lower bounds, which can change dynamically
           ubnd="numeric", #upper bounds, which can change dynamically
           type="character", # description of the organism
           medium="character", #character vector with exchange reactions of Organism
           lpobj="sysBiolAlg_fba", # sybil optimization object
           fbasol="list", # fba solution
           lyse="logical", # boolean variable for lysis
           feat="list", #list containing conditional features for the object (contains at the momement only biomass components for lysis)
           deathrate="numeric", # factor by which growth is reduced
           duplirate="numeric", # grow cut-off for test of duplication
           growthlimit="numeric",
           growtype="character" # functional type for growth (linear or exponential)
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(model, typename=mod_desc(model), algo="fba", ex="EX_", deathrate, duplirate,
                     growthlimit, growtype="exponential", lyse=F, feat=list(), ...){ #the constructor requires the model, after that it is not stored anymore
  rxname = react_id(model)
  lpobject <- sysBiolAlg(model, algorithm=algo)
  fbasol <- optimizeProb(lpobject)
  names(fbasol$fluxes) = rxname
  lobnd = lowbnd(model)
  names(lobnd) = rxname
  upbnd = uppbnd(model)
  names(upbnd) = rxname
  if(lyse){
    stoch <- S(model)[,which(model@obj_coef==1)] #find stochiometry of biomass components
    names(stoch) <- met_id(model)
    feat[["biomass"]] <- stoch[-which(stoch==0)]
  }
  new("Organism", lbnd=lobnd, ubnd=upbnd, type=typename, medium=rxname[grep(ex, rxname)],
      lpobj=lpobject, fbasol=fbasol, lyse=lyse, feat=feat, deathrate=deathrate, duplirate=duplirate,
      growthlimit=growthlimit, growtype=growtype, ...)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("lbnd", function(object){standardGeneric("lbnd")})
setMethod("lbnd", "Organism", function(object){return(object@lbnd)})
setGeneric("ubnd", function(object){standardGeneric("ubnd")})
setMethod("ubnd", "Organism", function(object){return(object@ubnd)})
setGeneric("type", function(object){standardGeneric("type")})
setMethod("type", "Organism", function(object){return(object@type)})
setGeneric("medium", function(object){standardGeneric("medium")})
setMethod("medium", "Organism", function(object){return(object@medium)})
setGeneric("lpobj", function(object){standardGeneric("lpobj")})
setMethod("lpobj", "Organism", function(object){return(object@lpobj)})
setGeneric("fbasol", function(object){standardGeneric("fbasol")})
setMethod("fbasol", "Organism", function(object){return(object@fbasol)})
setGeneric("lyse", function(object){standardGeneric("lyse")})
setMethod("lyse", "Organism", function(object){return(object@lyse)})
setGeneric("feat", function(object){standardGeneric("feat")})
setMethod("feat", "Organism", function(object){return(object@feat)})
setGeneric("deathrate", function(object){standardGeneric("deathrate")})
setMethod("deathrate", "Organism", function(object){return(object@deathrate)})
setGeneric("duplirate", function(object){standardGeneric("duplirate")})
setMethod("duplirate", "Organism", function(object){return(object@duplirate)})
setGeneric("growthlimit", function(object){standardGeneric("growthlimit")})
setMethod("growthlimit", "Organism", function(object){return(object@growthlimit)})
setGeneric("growtype", function(object){standardGeneric("growtype")})
setMethod("growtype", "Organism", function(object){return(object@feat)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("constrain", function(object, reacts, lb, dryweight, time){standardGeneric("constrain")})
setMethod("constrain", "Organism", function(object, reacts, lb, dryweight, time){
  lobnd <- object@lbnd*dryweight*time #costrain according to flux definition: mmol/(gDW*hr)
  lobnd[reacts] <- ifelse(lb<=lobnd[reacts], lobnd[reacts], lb) #check if lower bounds in biological relevant range
  return(lobnd)
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

setGeneric("consume", function(object, sublb, cutoff=1e-6){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, sublb, cutoff=1e-6){
  if(object@fbasol$obj>=cutoff){
    flux = object@fbasol$fluxes[object@medium]
    flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux))
    sublb[names(flux)] = round(sublb[names(flux)]+flux, 6)
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

#function for letting bacteria grow by adding the calculated growthrate to the already present growth value -> linear growth
#requires as input: organism object

setGeneric("growLin", function(object, growth){standardGeneric("growLin")})
setMethod("growLin", "Organism", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- object@fbasol$obj + growth
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})

#function for letting bacteria grow by adding the calculated growthrate multiplied with the current growth plus to the already present growth value -> exp growth
#requires as input: organism object

setGeneric("growExp", function(object, growth){standardGeneric("growExp")})
setMethod("growExp", "Organism", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- (object@fbasol$obj * growth + growth)
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})

#function for lysis of bacterial cells by adding biomass_compounds * growth to the medium

setGeneric("lysis", function(object, medium, sublb, factor=object@growthlimit, csuffix="\\[c\\]", esuffix="(e)", ex="EX_"){standardGeneric("lysis")})
setMethod("lysis", "Organism", function(object, medium, sublb, factor=object@growthlimit, csuffix="\\[c\\]", esuffix="(e)", ex="EX_"){
  stoch = object@feat[["biomass"]]
  names(stoch) <- paste(ex,gsub(csuffix, esuffix, names(stoch)),sep='')
  lysate = round(abs(na.omit(stoch[medium]))*factor, 6)
  sublb[names(lysate)] = sublb[names(lysate)] + lysate
  return(sublb)
})

#function to get moore-neighbourhood of a bac together with its relative position

setGeneric("getHood", function(object, occmat, x, y){standardGeneric("getHood")})
setMethod("getHood", "Organism", function(object, occmat, x, y){
  occmat <- as.matrix(occmat) #dangerous!
  if(x-1==0) dx=0 else dx=1
  if(x+1>nrow(occmat)) dx2=0 else dx2=1
  if(y-1==0) dy=0 else dy=1
  if(y+1>ncol(occmat)) dy2=0 else dy2=1
  return(list(as.matrix(occmat[,(y-dy):(y+dy2)])[(x-dx):(x+dx2),], c(1+dx,1+dy)))
})


#function to check if the there is some free place in the neighbourhood

setGeneric("emptyHood", function(object, occmat, x, y){standardGeneric("emptyHood")})
setMethod("emptyHood", "Organism", function(object, occmat, x, y){
  hood <- getHood(object, occmat, x, y)
  free <- which(hood[[1]]==0, arr.ind = T)
  if(nrow(free) == 0) return(NULL)
  else {
    abs <- free[sample(nrow(free),1),]
    abs[1] <- abs[1] - hood[[2]][1] + x
    abs[2] <- abs[2] - hood[[2]][2] + y
    return(abs)
  }
})

#show function for class Organism

removeMethod(show, signature(object="Organism"))
setMethod(show, signature(object="Organism"), function(object){
  print(paste('Organism ',object@type,' of class Organism.',sep=''))
})
