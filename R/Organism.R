########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

#' Structure of the S4 class "Organism"
#' 
#' Structure of the S4 class \code{Organism} representing the organisms present in the environment.
#' @export Organism
#' @exportClass Organism
#' @import sybil
#' @importFrom stats na.omit
#' @rdname Organism
#'
#' @slot lbnd A numeric vector containing the lower bounds of the model structure.
#' @slot ubnd A numeric vector containing the upper bounds of the model structure.
#' @slot type A character vector containing the description of the organism.
#' @slot medium A character vector containing all exchange reactions of the organism.
#' @slot lpobj A sybil optimization object containing the linear programing problem.
#' @slot fbasol A list with the solutions of the flux balance analysis.
#' @slot lyse A boolean variable indicating if the organism should lyse after death.
#' @slot feat A list containing conditional features for the object (contains at the momement only biomass components for lysis).
#' @slot deathrate A numeric value giving the factor by which the biomass should be reduced in every iteration if no growth is possible (default (E.coli): 0.21 pg)
#' @slot minweight A numeric value giving the growth limit at which the organism dies. (default (E.coli): 0.083 pg)
#' @slot growtype A character vector giving the functional type for growth (linear or exponential).
#' @slot kinetics A List containing Km and v_max values for each reactions.
#' @slot speed A integer vector representing the speed by which bacterium is moving (given by cell per iteration).
#' @slot cellarea A numeric value indicating the surface that one organism occupies (default (E.coli): 4.42 mu_m^2)
#' @slot maxweight A numeric value giving the maximal dry weight of single organism (default (E.coli): 1.172 pg)
#' @slot cellweight_mean A numeric giving the mean of starting biomass (default (E.coli): 0.489 pg)
#' @slot cellweight_sd A numeric giving the standard derivation of starting biomass (default (E.coli): 0.132 pg)
#' @slot model Object of class sybil::modelorg containging the genome sclae metabolic model
#' @slot algo Algorithm to be used during optimization (default fba)
#' @slot rbiomass Name of biomass reactions which is used for growth model (set automatically but needs input if objective is not biomass optimization)
#' @slot limit_growth If true then a upper bound on growth will be set, see maxweight (default: True).
#' @slot coupling_constraints List with coupling parameters.
#' @slot predator Name of organism which can kill this one. 
setClass("Organism",
         representation(
           lbnd="numeric",
           ubnd="numeric",
           type="character",
           medium="character",
           lpobj="sysBiolAlg",
           fbasol="list",
           lyse="logical",
           feat="list",
           deathrate="numeric",
           minweight="numeric",
           growtype="character",
           kinetics="list",
           cellarea="numeric",
           maxweight="numeric",
           cellweight_mean = "numeric",
           cellweight_sd = "numeric",
           speed="numeric",
           model="modelorg",
           algo="character",
           rbiomass="character",
           limit_growth="logical",
           coupling_constraints="list",
           predator="character"
         ),
         prototype(
           deathrate = 0.21,
           minweight = 0.083,
           growtype = "exponential",
           kinetics = list(),
           cellarea = 4.42,
           maxweight = 1.172,
           cellweight_mean = 0.489,
           cellweight_sd = 0.132,
           speed = 2,
           limit_growth = TRUE
         )
)


#' Find biomass reaction in model
#' @description Helper function to search for biomass reaction in available reactions of a model
#' 
#' @param model Object of class sybil::modelorg containging the genome sclae metabolic model
#' @param keys Vector with strings which are used to find biomass reaction in model
#' @return Vector with reaction ids for biomass reaction(s)
findrBiomass <- function(model, keys=c("biom", "cpd11416")){
  ex     <- sybil::findExchReact(model)
  ex_pos <- ex@react_pos
  ex_biom<- c(grep(paste0(keys, collapse = "|"), ex@met_id, ignore.case = TRUE),
              grep(paste0(keys, collapse = "|"), met_name(model)[ex@met_pos], ignore.case = TRUE))
  rbio <- vector()
  for(k in keys){
    idx <- grep(k, sybil::react_id(model), ignore.case = TRUE)
    if(length(idx)==0) idx <- grep(k, sybil::react_name(model), ignore.case = TRUE)
    if(length(idx)>0) rbio <- c(rbio, idx)
  }
  if( length(ex_biom) > 0 ){ # take care of biomass metabolite
    rbio   <- c(rbio, ex@react_pos[ex_biom])
    ex_pos <- setdiff(ex_pos, ex@react_pos[ex_biom]) 
  } 
  if(length(rbio)==0) return(NULL)
  
  rbio <- setdiff(rbio, ex_pos) # exclude exchange reactions
  return(sybil::react_id(model)[rbio])
}



########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Organism-class}}
#' 
#' @export
#' @name Organism-constructor
#' 
#' @param model Object of class sybil::modelorg containging the genome sclae metabolic model
#' @param algo A single character string giving the name of the algorithm to use. See \link[sybil]{SYBIL_SETTINGS}
#' @param ex Identifier for exchange reactions
#' @param ex_comp Defining exchange reactions whose compounds should be added to the medium of the arena (default: all)
#' @param csuffix suffix for intern metabolites used by lysis function.
#' @param esuffix suffix for external metabolites used by lysis function.
#' @param feat A list containing conditional features for the object (contains at the momement only biomass components for lysis).
#' @param typename A string defining the name (set to model name in default case)
#' @param lyse A boolean variable indicating if the organism should lyse after death.
#' @param setExInf Enable if all lower bounds of exchange reaction which are set to zero (i.e. no uptake possible!) should be set to -infitity (default: true)
#' @param setAllExInf Enable if all lower bounds of exchange reaction should be set to -infitity (default: false)
#' @param coupling_constraints List with coupling parameters.
#' @param predator Name of organism which can kill this one.
#' @param ... Arguments of \code{\link{Organism-class}}
#' @return Object of class Organism
Organism <- function(model, algo="fba", ex="EX_", ex_comp=NA, csuffix="\\[c\\]", esuffix="\\[e\\]", lyse=FALSE, feat=list(), 
                     typename=NA, setExInf=TRUE, setAllExInf=FALSE, coupling_constraints=list(), predator="", ...){
  pot_biomass <- findrBiomass(model)
  if(all(model@obj_coef==0)){
    if(length(pot_biomass)==0) stop("No objection function set in model")
    print("No objective function set, try set automatically...")
    print(paste("found:", pot_biomass))
    print(paste("set new biomass function:", pot_biomass[1]))
    sybil::obj_coef(model)[which(sybil::react_id(model)==pot_biomass[1])] <- 1
    rbiomass <- pot_biomass[1]
  }else{
    idx_obj <- which(sybil::obj_coef(model)!=0)
    idx_bio <- which(sybil::react_id(model) %in% pot_biomass)
    if(any(idx_obj %in% idx_bio)) rbiomass <- sybil::react_id(model)[idx_obj[which(idx_obj %in% idx_bio)]] # easy case: objective is (also) optimiziation of biomass
    else if(length(idx_bio)>0){ # if objective is not optimization of biomass
      cat("Optimization of biomass seems to be not the objective. Even if not optimized, a biomass reactions is needed for the growth model.\n")
      cat("Objective functions:", sybil::react_id(model)[idx_obj], "\n")
      cat("Available biomass reactions:", sybil::react_id(model)[idx_bio], "\n")
      cat("Biomass reaction used for growth model:", sybil::react_id(model)[idx_bio][1], "\n")
      rbiomass <- sybil::react_id(model)[idx_bio][1]
    }else{ # if no biomass reaction is found
      cat("Optimization of biomass seems to be not the objective. Even if not optimized, a biomass reactions is needed for the growth model.\n")
      cat("Objective functions ID:", sybil::react_id(model)[idx_obj], "\t", "Objective functions name:", sybil::react_name(model)[idx_obj], "\n")
      if(length(idx_obj) == 0){
        stop("No objective found for the model")  
      }
      rbiomass <- sybil::react_id(model)[idx_obj]
      warning("No biomass objective set. Please check that current objective makes sense.")
    }
  }
  if(algo=="coupling" & length(coupling_constraints)==0) stop("Please provide couling constraints for coupling!")
  
  if(is.na(ex)){
    medc <- sybil::react_id(sybil::findExchReact(model))
    names(medc) <- model@met_name[sybil::findExchReact(model)@met_pos]
  }else{
    exf <- sybil::findExchReact(model)
    medc <- exf@react_id[grep(ex, exf@react_id)]
    names(medc) <- model@met_name[exf@met_pos[grep(ex, exf@react_id)]]
  }
  if(!all(duplicated(medc)==FALSE)){
    warning("Model file contains duplicated reaction IDs. Attention, duplicated exchange reaction will be removed.")
    print(medc[which(medc==medc[duplicated(medc)])])
    medc <- unique(medc)}
  if(!is.na(ex_comp)){
    medc <- medc[grep(ex_comp, medc)]
  }
  
  rxname<- sybil::react_id(model)
  lobnd <- sybil::lowbnd(model)
  upbnd <- sybil::uppbnd(model) 
  names(lobnd) = rxname
  if( setAllExInf ){
    lobnd[ which( names(lobnd) %in% medc ) ] <- -1000
  }
  if( setExInf ){ # if setExInf is true then set lower bound of all exchange reactions which have zero values to -INF
    lobnd[ which(names(lobnd) %in% medc & lobnd==0) ] <- -1000
  }
  lobnd.ex <- lobnd[match(medc, rxname)]
  lobnd.ex.med <- stats::median(lobnd.ex[ which( lobnd.ex < 0 & lobnd.ex > -1000 ) ])
  if( !is.na(lobnd.ex.med) ){
    print(paste0("Median lower bound for non-zero and non-Inf exchanges is:", round(lobnd.ex.med), 6))
    if( lobnd.ex.med > -10 ){
      warning("Many lower bounds of of the model seems to be set to non -infinity. Please be aware that they will be used as maximal uptake rates even when the available medium is more abundant! (set setAllExInf=TRUE to reset all exchanges to -INF)")
      #print( lobnd.ex[ which( lobnd.ex < 0 & lobnd.ex >  lobnd.ex.med ) ] )
    }    
  }
  
  if(is.na(typename)) typename <- ifelse( length(sybil::mod_desc(model)) > 0, sybil::mod_desc(model), 
                                          ifelse( length(sybil::mod_name(model)) > 0, sybil::mod_name(model), 
                                                  ifelse( length(sybil::mod_id(model)) > 0, sybil::mod_id(model), "test" )))
  if(algo=="coupling"){
    #lpobject <- sybil::sysBiolAlg(model, algorithm="mtfEasyConstraint2", easyConstraint=coupling_constraints, pFBAcoeff = 1e-5)
    lpobject <- sybil::sysBiolAlg(model, algorithm="mtfEasyConstraint", easyConstraint=coupling_constraints)
  }else lpobject <- sybil::sysBiolAlg(model, algorithm=algo)
  fbasol <- sybil::optimizeProb(lpobject, react=1:length(lobnd), ub=upbnd, lb=lobnd)
  names(fbasol$fluxes) = rxname
  upbnd = sybil::uppbnd(model)
  names(upbnd) = rxname
  
  #browser()
  if(lyse | predator != ""){
    stochmat <- as.matrix(sybil::S(model))
    colnames(stochmat) <- rxname
    rownames(stochmat) <- sybil::met_id(model)
    stoch <- stochmat[,rbiomass]
    biomets <- stoch[-which(stoch==0)]
    exs <- sybil::findExchReact(model)
    extrans <- sybil::react_id(exs)
    names(extrans) <- sybil::met_id(exs)
    extrans <- extrans[which(extrans %in% medc)]
    names(biomets) <- gsub(csuffix,esuffix,names(biomets))
    biomets <- biomets[which(names(biomets) %in% names(extrans))]
    names(biomets) <- extrans[names(biomets)]
    feat[["biomass"]] <- biomets
  }

  new("Organism", lbnd=lobnd, ubnd=upbnd, type=typename, medium=medc, lpobj=lpobject,
      fbasol=fbasol, feat=feat, lyse=lyse, model=model, algo=algo,rbiomass=rbiomass, 
      coupling_constraints=coupling_constraints, predator=predator, ...)
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
setGeneric("minweight", function(object){standardGeneric("minweight")})
setMethod("minweight", "Organism", function(object){return(object@minweight)})
setGeneric("growtype", function(object){standardGeneric("growtype")})
setMethod("growtype", "Organism", function(object){return(object@feat)})
setGeneric("kinetics", function(object){standardGeneric("kinetics")})
setMethod("kinetics", "Organism", function(object){return(object@kinetics)})
setGeneric("speed", function(object){standardGeneric("speed")})
setMethod("speed", "Organism", function(object){return(object@speed)})
setGeneric("model", function(object){standardGeneric("model")})
setMethod("model", "Organism", function(object){return(object@model)})
setGeneric("maxweight", function(object){standardGeneric("maxweight")})
setMethod("maxweight", "Organism", function(object){return(object@maxweight)})
setGeneric("cellweight_mean", function(object){standardGeneric("cellweight_mean")})
setMethod("cellweight_mean", "Organism", function(object){return(object@cellweight_mean)})
setGeneric("cellarea", function(object){standardGeneric("cellarea")})
setMethod("cellarea", "Organism", function(object){return(object@cellarea)})
setGeneric("algo", function(object){standardGeneric("algo")})
setMethod("algo", "Organism", function(object){return(object@algo)})
setGeneric("coupling_constraints", function(object){standardGeneric("coupling_constraints")})
setMethod("coupling_constraints", "Organism", function(object){return(object@coupling_constraints)})
setGeneric("predator", function(object){standardGeneric("predator")})
setMethod("predator", "Organism", function(object){return(object@predator)})


########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for constraining the models based on metabolite concentrations
#'
#' @description The generic function \code{constrain} changes the constraints of the model representation of an organism.
#' @export
#' @rdname constrain
#'
#' @param object An object of class Organisms.
#' @param reacts A character vector giving the names of reactions which should be constrained.
#' @param lb A numeric vector giving the constraint values of lower bounds (e.g. avaible metabolite concentrations
#' @param dryweight A number giving the current dryweight of the organism.
#' @param tstep A number giving the time intervals for each simulation step.
#' @param scale A numeric defining the scaling (units for linear programming has to be in certain range)
#' @param j debuging index to track cell
#' @param cutoff value used to define numeric accuracy while interpreting optimization results
#' @return Returns the lower bounds, which carry the constraints and names of relevant reactions.
#' @details The constraints are calculated according to the flux definition as mmol/(gDW*hr) with the parameters \code{dryweight} and \code{tstep}.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize an organism
#' lobnds <- constrain(org,org@medium,org@lbnd[org@medium],1,1,1,1,1)
setGeneric("constrain", function(object, reacts, lb, dryweight, tstep, scale, j, cutoff=1e-6){standardGeneric("constrain")})
#' @export
#' @rdname constrain
setMethod("constrain", "Organism", function(object, reacts, lb, dryweight, tstep, scale, j, cutoff=1e-6){
  reacts = unique(reacts)
  lb = lb[reacts]
  lobnd <- object@lbnd
  upbnd <- object@ubnd
  if(dryweight<Inf){lobnd[reacts] <- object@lbnd[reacts]*(dryweight/object@cellweight_mean)*tstep} #costrain according to flux definition: mmol/(gDW*hr)
  #lobnd[reacts] <- ifelse(lb<=lobnd[reacts], ifelse(lobnd[reacts]==0, lb, lobnd[reacts]), lb) #check if lower bounds in biological relevant range
  lobnd[reacts] <- ifelse(lb<=lobnd[reacts], lobnd[reacts], lb) #check if lower bounds in biological relevant range
  #if(j==1) browser()
  if(length(object@kinetics) != 0){
    lobnd[names(object@kinetics)] <- unlist(lapply(names(object@kinetics), function(name){
      Km  <- (object@kinetics[[name]][["Km"]]*0.01*scale)*10^12 #scale mM to fmol/gridcell
      vmax <- object@kinetics[[name]][["vmax"]]*(dryweight/object@cellweight_mean)*tstep #scale to fmol/pgDW*hr
      s   <- -lb[name] # change sign to get concentrations
      lnew = -(vmax*s/(Km + s))*tstep
      if(abs(lnew)>s){if(-s>=lobnd[name]){lnew=-s}else{lnew=lobnd[name]}}
      return(lnew)
    }))
  }
  if( object@limit_growth ){ # set upper bound for growth
    growth_limit <- (object@maxweight*1.5) - dryweight # 1.5 is factor to allow some room for the biomass to accumulate
    # if growth is too fast, then growth_limit could become negative (hight dryweight). Use low cutoff value so that dryweight can be decreased by cell division.
    upbnd[which(object@model@react_id == object@rbiomass)] <- ifelse( growth_limit>0, growth_limit, cutoff )
  }
  
  return(list(lobnd, upbnd))
})


#' @title Function to set Michaelis-Menten kinetics for uptake of a substance
#'
#' @description The generic function \code{setKinetics} provides kinetics for exchange reactions.
#' @export
#' @rdname setKinetics
#'
#' @param object An object of class Organisms.
#' @param exchangeR Name of an exchange reaction
#' @param Km Parameter Michaelis-Menten-Kinetics (in mM)
#' @param vmax Parameter Michaelis-Menten-Kinetics (in mmol/(g*h))
#'
setGeneric("setKinetics", function(object, exchangeR, Km, vmax){standardGeneric("setKinetics")})
#' @export
#' @rdname setKinetics
setMethod("setKinetics", "Organism", function(object, exchangeR, Km, vmax){
  if( !(exchangeR %in% names(object@lbnd)) ){
    stop("Incorrect exchange reaction setup for kinetics")
  }
  nKinetics <- object@kinetics
  if(exchangeR %in% names(nKinetics)){
    warning("Kinetics for this exchange reaction already given, set to new value")
    nKinetics[[exchangeR]] <- list()
  }
  nKinetics[[exchangeR]] <- c(nKinetics[[exchangeR]], "Km"=Km, "vmax"=vmax)
  object@kinetics <- nKinetics
  return(object)
})


#' @title Function for computing the linear programming according to the model structure 
#'
#' @description The generic function \code{optimizeLP} implements a linear programming based on the problem structure and refined constraints.
#' @export
#' @rdname optimizeLP
#'
#' @param object An object of class Organisms.
#' @param lpob A linear programing object encoding the problem to solve.
#' @param lb A numeric vector giving the constraint values of lower bounds.
#' @param ub A numeric vector giving the constraint values of upper bounds.
#' @param cutoff value used to define numeric accuracy while interpreting optimization results
#' @param j debuging index to track cell
#' @param sec_obj character giving the secondary objective for a bi-level LP if wanted. Use "mtf" for minimizing total flux, "opt_rxn" for optimizing a random reaction, "opt_ex" for optimizing a random exchange reaction, and "sumex" for optimizing the sum of all exchange fluxes.
#' @param with_shadow True if shadow costs should be stored (default false).
#' @return Modified problem object according to the constraints and then solved with \code{optimizeProb}.
#' @details The parameter for sec_obj can be used to optimize a bi-level LP with a secondary objective if wanted. This can be helpful to subselect the solution space and create less alternative optimal solution. The secondary objective can be set to "mtf" to minimize the total flux, to simulate minimal enzyme usage of an organisms. If set to "opt_rxn" or "opt_ex", the secondary objective is picked as a random reaction or exchange reaction respectively everytime a fba is performed. This means that every individual of a population will select a different secondary reaction to optimize. The "sumex" option maximizes the secretion of products.
#' @seealso \code{\link{Organism-class}}, \code{\link{optimizeProb}} and \code{\link{sysBiolAlg}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' org <- Organism(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a organism
#' org@fbasol <- optimizeLP(org)
setGeneric("optimizeLP", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd, cutoff=1e-6, j, sec_obj="none", with_shadow=FALSE){standardGeneric("optimizeLP")})
#' @export
#' @rdname optimizeLP
setMethod("optimizeLP", "Organism", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd, cutoff=1e-6, j, sec_obj="none", with_shadow=FALSE){ 
  fbasl <- sybil::optimizeProb(lpob, react=1:length(lb), ub=ub, lb=lb, resetChanges = FALSE) # resetChanges needed for cplex reduced/shadow costs
  #browser()
  #cat("stat:",fbasl$stat, "\tobj:",fbasl$obj,"\n")
  #summary(getRedCosts(lpob@problem)[ex@react_pos])
  solve_ok <- getLPstat(opt_sol=fbasl, solver=lpob@problem@solver)
  if(is.na(solve_ok)){solve_ok = FALSE}
  if(!solve_ok | is.na(fbasl$obj) | fbasl$obj<cutoff){
    fbasl$obj <- 0
    fbasl$fluxes <- rep(0, length(fbasl$fluxes))
  }

  if(sec_obj!="none" && fbasl$obj!=0){
    switch(sec_obj,
           mtf =   {mod = sybil::changeBounds(object@model, react=1:length(lb), lb=lb, ub=ub);
                    mtf_sol <- sybil::optimizeProb(mod, algorithm="mtf", mtfobj=fbasl$obj, retOptSol=F);
                    mtf_ok  <- getLPstat(opt_sol=mtf_sol, solver=lpob@problem@solver)
                    if( mtf_ok ) fbasl$fluxes = mtf_sol$fluxes},
           opt_rxn={mod = sybil::changeBounds(object@model, object@model@react_id, lb=lb, ub=ub);
                    mod = sybil::changeBounds(mod, mod@react_id[which(obj_coef(mod)==1)], lb=fbasl$obj, ub=fbasl$obj);
                    mod <- sybil::changeObjFunc(mod, sample(mod@react_id,1));
                    solrand <- sybil::optimizeProb(mod, lpdir=sample(c("max","min"),1));
                    solrand@lp_obj
                    fbasl$fluxes = sybil::getFluxDist(solrand)},
           opt_ex ={mod = sybil::changeBounds(object@model, object@model@react_id, lb=lb, ub=ub);
                    mod = sybil::changeBounds(mod, mod@react_id[which(obj_coef(mod)==1)], lb=fbasl$obj, ub=fbasl$obj);
                    mod <- sybil::changeObjFunc(mod, sample(object@medium,1));
                    solrand <- sybil::optimizeProb(mod, lpdir=sample(c("max","min"),1));
                    solrand@lp_obj
                    fbasl$fluxes = sybil::getFluxDist(solrand)},
           sumex = {mod <- sybil::changeBounds(object@model, object@model@react_id, lb=lb, ub=ub);
                    old_obj_coef <- which(obj_coef(mod)==1)
                    mod <- sybil::changeBounds(mod, mod@react_id[old_obj_coef], lb=fbasl$obj*0.05);
                    mod <- sybil::changeObjFunc(mod, unname(object@medium));
                    solrand <- sybil::optimizeProb(mod);
                    fbasl$fluxes = sybil::getFluxDist(solrand);
                    fbasl$obj <- fbasl$fluxes[old_obj_coef]},
           stop("Secondary objective not suported!"))
  }
  # get shadow cost from dual problem (only for glpk and cplex)
  if(!with_shadow | !solve_ok){
    shadow=NULL
  } else{
    ex <- findExchReact(object@model)
    idx <- c(ex@react_pos, grep(object@rbiomass, object@model@react_id))
    switch(lpob@problem@solver, # use reduced costs, shadow costs are not supported by sybil and direct access is causing problems with cplex
           glpkAPI =  {shadow <- sybil::getRedCosts(lpob@problem)[idx]},
           cplexAPI = {shadow <- sybil::getRedCosts(lpob@problem)[idx]},
           shadow=NULL)
    names(shadow) <- object@model@react_id[idx]
    #if(all(is.na(shadow))) browser()
  }
  
  names(fbasl$fluxes) <- names(object@lbnd)
  return(list(fbasl, shadow))
})


#' @title Function to account for the consumption and production of substances
#'
#' @description The generic function \code{consume} implements the consumption and production of substances based on the flux of exchange reactions of organisms
#' @export
#' @rdname consume
#'
#' @param object An object of class Organisms.
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param cutoff A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @param bacnum Integer indicating the number of bacteria individuals per gridcell
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @return Returns the updated vector containing the substance concentrations in the current position of the individual of interest.
#' @details The consumption is implemented by adding the flux of the exchange reactions to the current substance concentrations.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' NULL
setGeneric("consume", function(object, sublb, cutoff=1e-6, bacnum, fbasol){standardGeneric("consume")})
#' @export
#' @rdname consume
setMethod("consume", "Organism", function(object, sublb, cutoff=1e-6, bacnum, fbasol){
  if(fbasol$obj>=cutoff && !is.na(fbasol$obj)){
    flux = fbasol$fluxes[object@medium]*bacnum #scale flux to whole population size
    #flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux)) ?
    sublb[names(flux)] = round(sublb[names(flux)]+flux, round(-log10(cutoff))) # use cutoff also in this case
  }
  return(sublb)
})


#' @title Function to extract the phenotype of an organism object
#'
#' @description The generic function \code{getPhenotype} implements an identification of organism phenotypes.
#' @export
#' @rdname getPhenotype
#'
#' @param object An object of class Organisms.
#' @param cutoff A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @param par A boolean indicating if running in parallel mode.
#' @return Returns the phenotype of the organisms where the uptake of substances is indicated by a negative and production of substances by a positive number
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages. Uptake of substances is indicated by a negative and production of substances by a positive number.
#' @seealso \code{\link{Organism-class}}, \code{\link{checkPhen}} and \code{\link{minePheno}}
setGeneric("getPhenotype", function(object, cutoff=1e-6, fbasol, par=FALSE){standardGeneric("getPhenotype")})
#' @export
#' @rdname getPhenotype
setMethod("getPhenotype", "Organism", function(object, cutoff=1e-6, fbasol, par=FALSE){
  exflux=fbasol$fluxes[object@medium]
  exflux=ifelse(abs(exflux)<cutoff,0,1)*exflux
  exflux=ifelse(exflux>0,1,exflux)
  exflux=ifelse(exflux<0,2,exflux)
  if(!par) return(exflux[which(exflux!=0)]) else return(exflux)
})



#' @title Function for letting organisms grow linearly
#'
#' @description The generic function \code{growLin} implements a growth model of organisms in their environment.
#' @export
#' @rdname growLin
#'
#' @param object An object of class Organisms.
#' @param biomass A number indicating the current biomass, which has to be updated. 
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @param tstep A number giving the time intervals for each simulation step.
#' @return Returns the updated biomass of the organisms of interest.
#' @details Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value.
#' @seealso \code{\link{Organism-class}} and \code{\link{optimizeLP}}
setGeneric("growLin", function(object, biomass, fbasol, tstep){standardGeneric("growLin")})
#' @export
#' @rdname growLin
setMethod("growLin", "Organism", function(object, biomass, fbasol, tstep){
  growth <- fbasol$fluxes[object@rbiomass]
  if(growth > 0){
    grow_accum <- growth + biomass
  } else grow_accum <- biomass - object@deathrate*tstep
  #cat("\t", growth, biomass, grow_accum, "\n")
  return(grow_accum)
})



#' @title Function for letting organisms grow exponentially
#'
#' @description The generic function \code{growExp} implements a growth model of organisms in their environment.
#' @export
#' @rdname growExp
#'
#' @param object An object of class Organisms.
#' @param biomass A number indicating the current biomass, which has to be updated. 
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @param tstep A number giving the time intervals for each simulation step.
#' @return Returns the updated biomass of the organisms of interest.
#' @details Exponential growth of organisms is implemented by adding the calculated growthrate multiplied with the current growth calculated by \code{optimizeLP} plus to the already present growth value
#' @seealso \code{\link{Organism-class}} and \code{\link{optimizeLP}}
setGeneric("growExp", function(object, biomass, fbasol, tstep){standardGeneric("growExp")})
#' @export
#' @rdname growExp
setMethod("growExp", "Organism", function(object, biomass, fbasol, tstep){
  growth <- fbasol$fluxes[object@rbiomass]
  if(growth > 0){
    grow_accum <- (growth * biomass + biomass)
  } else grow_accum <- biomass - object@deathrate*biomass*tstep
  #cat("\t", growth, biomass, grow_accum, "\n")
  return(grow_accum)
})



#' @title Lysis function of organismal cells by adding biomass_compounds to the medium
#'
#' @description The generic function \code{lysis} implements cell lysis by the stochiometric concentration of the biomass compounds of organisms to the concentration of substances in the environment
#' @export
#' @rdname lysis
#'
#' @param object An object of class Organisms.
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param factor A number given the factor with which the biomass compound concentrations are multiplied to achieve the final concentration which is added to the environment
#' @return Returns the updated vector containing the substance concentrations in the current position of the dead individual of interest.
#' @details Lysis is implemented by taking the intersect between biomass compounds and the substances in the environment and adding the normalized stochiometric concentrations of the biomass compounds to the medium.
#' @seealso \code{\link{Organism-class}} and \code{\link{optimizeLP}}
#' @examples
#' NULL
setGeneric("lysis", function(object, sublb, factor=object@minweight){standardGeneric("lysis")})
#' @export
#' @rdname lysis
setMethod("lysis", "Organism", function(object, sublb, factor=object@minweight){
  stoch = object@feat[["biomass"]]
  lysate = round(abs(stoch)*factor, 6)
  sublb[names(lysate)] = sublb[names(lysate)] + lysate
  return(sublb)
})

#' @title Function to check if the there is a free place in the Moore neighbourhood
#'
#' @description The generic function \code{emptyHood} gives a free space which is present in the Moore neighbourhood of an individual of interest.
#' @export
#' @rdname emptyHood
#'
#' @param object An object of class Organisms.
#' @param x A number giving the x position of the individual of interest in its environment.
#' @param y A number giving the y position of the individual of interest in its environment.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param pos A dataframe with all occupied x and y positions 
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @param inverse Return occupied positions instead
#' @return Returns the free position in the Moore neighbourhood, which is not occupied by other individuals. If there is no free space \code{NULL} is returned.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' NULL
setGeneric("emptyHood", function(object, pos, n, m, x, y, occupyM, inverse=FALSE){standardGeneric("emptyHood")})
#' @export
#' @rdname emptyHood
setMethod("emptyHood", "Organism", function(object, pos, n, m, x, y, occupyM, inverse=FALSE){
  occ <- which(occupyM!=0, arr.ind = T)[,c(2,1)]
  colnames(occ) <- c("x","y")
  pos <- rbind(pos, occ) # block also occupied areas
  xp = c(x-1,x,x+1)
  yp = c(y-1,y,y+1)
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  #xp = xp[xp>0 & xp<=n]
  #xp = xp[yp>0 & yp<=m]
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  if( !inverse ){
    freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))  
  }else{
    freenb=paste(pos$x,pos$y,sep='_')
  }
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
})

#' @title Function to check if the there is a free place in the Moore neighbourhood
#'
#' @description The generic function \code{NemptyHood} gives a free space which is present in the Moore neighbourhood of an individual of interest.
#' @export
#' @rdname NemptyHood
#'
#' @param object An object of class Organisms.
#' @param x A number giving the x position of the individual of interest in its environment.
#' @param y A number giving the y position of the individual of interest in its environment.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param pos A dataframe with all occupied x and y positions 
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @param inverse Return occupied positions instead
#' @return Returns the free position in the Moore neighbourhood, which is not occupied by other individuals. If there is no free space \code{NULL} is returned.
#' @seealso \code{\link{Organism-class}}
#' @examples
#' NULL
setGeneric("NemptyHood", function(object, pos, n, m, x, y, occupyM, inverse=FALSE){standardGeneric("NemptyHood")})
#' @export
#' @rdname NemptyHood
setMethod("NemptyHood", "Organism", function(object, pos, n, m, x, y, occupyM, inverse=FALSE){
  occ <- which(occupyM!=0, arr.ind = T)[,c(2,1)]
  colnames(occ) <- c("x","y")
  pos <- rbind(pos, occ) # block also occupied areas
  xp = c(x-1,x,x+1)
  yp = c(y-1,y,y+1)
  for(i in 2:object@speed){
    xp = c(xp,x-i,x+i)
    yp = c(yp,y-i,y+i)
  }
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  if( !inverse ){
    freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))  
  }else{
    freenb=paste(pos$x,pos$y,sep='_')
  }
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
})

#' @title Function for random movement of organisms
#'
#' @description The generic function \code{move} implements a random movement in the Moore neighbourhood of an individual.
#' @export
#' @rdname move
#'
#' @param object An object of class Organism.
#' @param j The number of the iteration of interest.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param pos A dataframe with all occupied x and y positions 
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @details Organisms move in a random position the Moore neighbourhood, which is not occupied by other individuals. If there is no free space the individuals stays in the same position.
#' @seealso \code{\link{Organism-class}}, \code{\link{emptyHood}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' move(bac,n=20,m=20,j=1,pos=arena@orgdat[,c('x','y')], occupyM=arena@occupyM)
setGeneric("move", function(object, pos, n, m, j, occupyM){standardGeneric("move")})
#' @export
#' @rdname move
setMethod("move", "Organism", function(object, pos, n, m, j, occupyM){
  if(object@speed == 1){
    freenb <- emptyHood(object, pos, n, m, pos[j,1], pos[j,2], occupyM)
  }else{
    freenb <- NemptyHood(object, pos, n, m, pos[j,1], pos[j,2], occupyM)
  }
  if(length(freenb) != 0){
    npos = freenb[sample(length(freenb),1)]
    npos = as.numeric(unlist(strsplit(npos,'_')))
    if(occupyM[npos[2], npos[1]] == 0){ # check if there is no obstacle
      pos[j,] = npos
    }
  }
  return(pos)
})

#show function for class Organism

setMethod(show, signature(object="Organism"), function(object){
  print(paste('Organism ',object@type,' of class Organism.',sep=''))
})

# Bac is a subclass of Organism containing bacteria specific features

########################################################################################################
###################################### BAC CLASS #######################################################
########################################################################################################

#' Structure of the S4 class "Bac"
#' 
#' Structure of the S4 class \code{Bac} inheriting from class \code{\link{Organism-class}} representing bacterial cells.
#' @export Bac
#' @exportClass Bac
#' @rdname Bac
#'
#' @slot deathrate A numeric value giving the factor by which the biomass should be reduced in every iteration if no growth is possible (default (E.coli): 0.21 pg)
#' @slot minweight A numeric value giving the growth limit at which the organism dies. (default (E.coli): 0.083 pg)
#' @slot cellarea A numeric value indicating the surface that one organism occupies (default (E.coli): 4.42 mu_m^2)
#' @slot maxweight A numeric value giving the maximal dry weight of single organism (default (E.coli): 1.172 pg)
#' @slot chem A character vector indicating name of substance which is the chemotaxis attractant. Empty character vector if no chemotaxis.
setClass("Bac",
         contains="Organism",
         representation(
           chem="character" # name of substance which is the chemotaxis attractant
         ),
         prototype(
           deathrate = 0.21,
           minweight = 0.083,
           cellarea = 4.42,
           maxweight = 1.172
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Bac-class}}
#'
#' @name Bac-Constructor
#' @export
#' 
#' @return Object of class \code{\link{Bac-class}}
#' @param model model
#' @param ... Arguments of \code{\link{Organism}}
#' @param chem A character vector indicating name of substance which is the chemotaxis attractant. Empty character vector if no chemotaxis.
Bac <- function(model, chem='', ...){
  new("Bac", Organism(model=model, ...), chem=chem)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("chem", function(object){standardGeneric("chem")})
setMethod("chem", "Bac", function(object){return(object@chem)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function implementing a growth model of a bacterium
#'
#' @description The generic function \code{growth} implements different growth models for an object of class Bac.
#' @export
#' @rdname growth
#'
#' @param object An object of class Bac.
#' @param population An object of class Arena.
#' @param j The index of the organism of interest in orgdat.
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @param tstep A number giving the time intervals for each simulation step.
#' @return Boolean variable of the \code{j}th individual indicating if individual died.
#' @details Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value. Exponential growth of organisms is implemented by adding the calculated growthrate multiplied with the current growth calculated by \code{optimizeLP} plus to the already present growth value
#' @seealso \code{\link{Bac-class}}, \code{\link{growLin}} and \code{\link{growExp}}
setGeneric("growth", function(object, population, j, occupyM, fbasol, tstep){standardGeneric("growth")})
#' @export
#' @rdname growth
setMethod("growth", "Bac", function(object, population, j, occupyM, fbasol, tstep){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$biomass <- growLin(object, popvec$biomass, fbasol, tstep)},
         "exponential"= {popvec$biomass <- growExp(object, popvec$biomass, fbasol, tstep)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'biomass'] <- popvec$biomass
  while( popvec$biomass > object@maxweight ){
  #if(popvec$biomass > object@maxweight){ 
    freenb <- emptyHood(object, neworgdat[,c('x','y')],
              population@n, population@m, popvec$x, popvec$y, occupyM)
    if(length(freenb) != 0){
      npos = freenb[sample(length(freenb),1)]
      npos = as.numeric(unlist(strsplit(npos,'_')))
      if(occupyM[npos[2], npos[1]] == 0){ # check if there is no obstacle
        popvec$biomass <- popvec$biomass/2
        daughter <- popvec
        daughter$biomass <- popvec$biomass
        daughter$x <- npos[1]
        daughter$y <- npos[2]
        neworgdat[nrow(neworgdat)+1,] <- daughter
        neworgdat[j,'biomass'] <- popvec$biomass
      }
    }else
      break
  }
  if(popvec$biomass < object@minweight){
    neworgdat[j,'biomass'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})


#' @title Function implementing a growth model of a bacterium
#'
#' @description The generic function \code{growth_par} implements different growth models for an object of class Bac.
#' @export
#' @rdname growth_par
#' 
#' @param object An object of class Bac.
#' @param population An object of class Arena.
#' @param j The index of the organism of interest in orgdat.
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @param tstep A number giving the time intervals for each simulation step.
#' @return A list 
setGeneric("growth_par", function(object, population, j, fbasol, tstep){standardGeneric("growth_par")})
#' @export
#' @rdname growth_par
setMethod("growth_par", "Bac", function(object, population, j, fbasol, tstep){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$biomass <- growLin(object, popvec$biomass, fbasol, tstep)},
         "exponential"= {popvec$biomass <- growExp(object, popvec$biomass, fbasol, tstep)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'biomass'] <- popvec$biomass
  if(popvec$biomass < object@minweight){
    neworgdat[j,'biomass'] <- NA
    dead <- T
  }
  return(list(dead, neworgdat))
})


# chemotaxis still needs to be edited

#' @title Function for chemotaxis of bacteria to their prefered substrate
#'
#' @description The generic function \code{chemotaxis} implements a bacterial movement in the Moore neighbourhood to the highest substrate concentration.
#' @export
#' @rdname chemotaxis
#'
#' @param object An object of class Bac.
#' @param population An object of class Arena.
#' @param j The number of the iteration of interest.
#' @param chemo The vector that contains the prefered substrate.
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @details Bacteria move to a position in the Moore neighbourhood which has the highest concentration of the prefered substrate, which is not occupied by other individuals. The prefered substance is given by slot \code{chem} in the \code{Bac} object. If there is no free space the individuals stays in the same position. If the concentration in the Moore neighbourhood has the same concentration in every position, then random movement is implemented.
#' @seealso \code{\link{Bac-class}} and \code{\link{emptyHood}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05, chem = "EX_o2(e)",
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' chemotaxis(bac,arena,1, "EX_glc(e)", arena@occupyM)
setGeneric("chemotaxis", function(object, population, j, chemo, occupyM){standardGeneric("chemotaxis")})
#' @export
#' @rdname chemotaxis
setMethod("chemotaxis", "Bac", function(object, population, j, chemo, occupyM){
  popvec <- population@orgdat[j,]
  attract <- population@media[[chemo]]@diffmat
  freenb <- emptyHood(object, population@orgdat[,c('x','y')],
                      population@n, population@m, popvec$x, popvec$y, occupyM)
  if(length(freenb) != 0){
    conc <- sapply(freenb, function(x, attract){
      npos = as.numeric(unlist(strsplit(x,'_')))
      return(attract[npos[1],npos[2]])
    }, attract=attract)
    abs <- freenb[which(conc==max(conc))]
    if(length(abs)!=1){
      abs <- abs[sample(1:length(abs),1)]
    }
    npos = as.numeric(unlist(strsplit(abs,'_')))
    #eval.parent(substitute(population@orgdat[j,]$x <- npos[1]))
    #eval.parent(substitute(population@orgdat[j,]$y <- npos[2]))
    return(npos)
  }
})

#function for one iteration for Bac class
#' @title Function for one simulation iteration for objects of Bac class
#'
#' @description The generic function \code{simBac} implements all neccessary functions for the individuals to update the complete environment. 
#' @export
#' @rdname simBac
#'
#' @param object An object of class Bac.
#' @param arena An object of class Arena defining the environment.
#' @param j The index of the organism of interest in orgdat.
#' @param bacnum integer indicating the number of bacteria individuals per gridcell
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param sec_obj character giving the secondary objective for a bi-level LP if wanted.
#' @param cutoff value used to define numeric accuracy.
#' @param pcut A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @param with_shadow True if shadow cost should be stores (default off).
#' @return Returns the updated enivironment of the \code{population} parameter with all new positions of individuals on the grid and all new substrate concentrations.
#' @details Bacterial individuals undergo step by step the following procedures: First the individuals are constrained with \code{constrain} to the substrate environment, then flux balance analysis is computed with \code{optimizeLP}, after this the substrate concentrations are updated with \code{consume}, then the bacterial growth is implemented with \code{growth}, the potential new phenotypes are added with \code{checkPhen}, finally the additional and conditional functions \code{lysis}, \code{move} or \code{chemotaxis} are performed. In case of many compounds in the vector of \code{chemotaxis}, the change of the position takes place by the order of the compounds in the vector of \code{chemotaxis}. Can be used as a wrapper for all important bacterial functions in a function similar to \code{simEnv}.
#' @seealso \code{\link{Bac-class}}, \code{\link{Arena-class}}, \code{\link{simEnv}}, \code{constrain}, \code{optimizeLP}, \code{consume}, \code{growth}, \code{checkPhen}, \code{lysis}, \code{move} and \code{chemotaxis}
#' @examples
#' NULL
setGeneric("simBac", function(object, arena, j, sublb, bacnum, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE){standardGeneric("simBac")})
#' @export
#' @rdname simBac
setMethod("simBac", "Bac", function(object, arena, j, sublb, bacnum, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE){
  predator_found <- FALSE
  if( object@predator != ""){
    pos  <- arena@orgdat[,c('x','y')]
    nb <- emptyHood(object, pos, arena@n, arena@m, pos[j,1], pos[j,2], arena@occupyM, inverse=T)  
    unlist(strsplit(nb,'_'))
    nb_types <- sapply(strsplit(nb,'_'), function(coord){
      arena@orgdat[arena@orgdat$x==coord[1] & arena@orgdat$y==coord[2], "type"]
    })
    nb_names <- unique(names(arena@specs)[unlist(nb_types)])
    if( object@predator %in% nb_names){
      predator_found <- TRUE
    }
  }
  if( predator_found ){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,], factor=arena@orgdat[j,"biomass"])))
    dead <- TRUE
    arena@orgdat[j, "biomass"] <- NA
    #print("died because of predator")
  }else{
    const <- constrain(object, object@medium, lb=-sublb[j,object@medium]/bacnum, #scale to population size
                       dryweight=arena@orgdat[j,"biomass"], tstep=arena@tstep, scale=arena@scale, j)
    lobnd <- const[[1]]; upbnd <- const[[2]]
    optimization <- optimizeLP(object, lb=lobnd, ub=upbnd, j=j, sec_obj=sec_obj, cutoff=cutoff, with_shadow=with_shadow)
    fbasol <- optimization[[1]]
    
    eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,], bacnum=bacnum, fbasol=fbasol, cutoff) )) #scale consumption to the number of cells?
    
    dead <- growth(object, arena, j, arena@occupyM, fbasol=fbasol, tstep=arena@tstep)
    arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, org=object, fbasol=fbasol, cutoff=pcut))
    
    type <- object@type
    arena@mflux[[type]]  <- arena@mflux[[type]] + fbasol$fluxes # remember active fluxes
    arena@shadow[[type]] <- arena@shadow[[type]]+ optimization[[2]]
    idx <- match(arena@mediac, names(fbasol$fluxes))
    exchanges <- data.frame(type,t(fbasol$fluxes[idx]))
    colnames(exchanges) <- c("species", unname(arena@mediac))
    arena@exchanges <- rbind(arena@exchanges, exchanges) # remember exchanges
  }
  
  
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,])))
  }
  if(!dead && !arena@stir && object@speed != 0){
    if(object@chem[1] == ''){
      pos <- arena@orgdat[,c('x','y')]
      mov_pos <- move(object, pos, arena@n, arena@m, j, arena@occupyM)
      arena@orgdat[,c('x','y')] <- mov_pos
    }else{
      for (v in seq_along(object@chem)){
      chemo <- object@chem[[v]]
      chemo_pos <- chemotaxis(object, arena, j, chemo, arena@occupyM)
      if(!is.null(chemo_pos)){arena@orgdat[j,c('x','y')] <- chemo_pos}
      }
    }
  }
  return(arena)
})

#' @title Function for one simulation iteration for objects of Bac class
#'
#' @description The generic function \code{simBac_par} implements all neccessary functions for the individuals to update the complete environment. 
#' @export
#' @rdname simBac_par
#'
#' @param object An object of class Bac.
#' @param arena An object of class Arena defining the environment.
#' @param j The index of the organism of interest in orgdat.
#' @param bacnum integer indicating the number of bacteria individuals per gridcell
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param lpobject linar programming object (copy of organism@lpobj) that have to be a deep copy in parallel due to pointer use in sybil.
#' @param sec_obj character giving the secondary objective for a bi-level LP if wanted.
#' @param cutoff value used to define numeric accuracy
#' @param with_shadow True if shadow cost should be stores (default off).
#' @return Returns the updated enivironment of the \code{population} parameter with all new positions of individuals on the grid and all new substrate concentrations.
#'
setGeneric("simBac_par", function(object, arena, j, sublb, bacnum, lpobject, sec_obj="none", cutoff=1e-6, with_shadow=FALSE){standardGeneric("simBac_par")})
#' @export
#' @rdname simBac_par
setMethod("simBac_par", "Bac", function(object, arena, j, sublb, bacnum, lpobject, sec_obj="none", cutoff=1e-6, with_shadow=FALSE){
  const <- constrain(object, object@medium, lb=-sublb[j,object@medium]/bacnum, #scale to population size
                     dryweight=arena@orgdat[j,"biomass"], tstep=arena@tstep, scale=arena@scale)
  lobnd <- const[[1]]; upbnd <- const[[2]]
  fbasol <- optimizeLP(object, lb=lobnd, ub=upbnd, j=j, sec_obj=sec_obj, cutoff=cutoff, with_shadow=with_shadow)

  sublb[j,] <- consume(object, sublb=sublb[j,], bacnum=bacnum, fbasol=fbasol, cutoff=1e-6) #scale consumption to the number of cells?
  
  growth <- growth_par(object, arena, j, fbasol, tstep=arena@tstep)
  dead <- growth[[1]]
  neworgdat <- growth[[2]]
  
  if(dead && object@lyse){
    sublb[j,] <- lysis(object, sublb[j,])
  }
  return(list(neworgdat[j,], sublb[j,], fbasol))
})




#show function for class Bac

setMethod(show, signature(object="Bac"), function(object){
  # print status of exchange reactions
  ex_lb <- object@lbnd[which(names(object@lbnd) %in% object@medium)]
  group_ex_lb <- split(ex_lb, factor(unlist(unname(ex_lb))))
  lapply(seq_along(group_ex_lb), function(i){
    if(as.numeric(names(group_ex_lb)[i])==0){
      print("Exchange reaction with _NO_ uptake set:")
    } else print(paste("Exchange reaction with uptake set of", names(group_ex_lb)[i]))
    print(names(group_ex_lb[[i]]))
    cat("\n")
  })
  
  print(paste('Bacterium ',object@type,' of class Bac.',sep=''))
})

# Human is a subclass of Organism containing human specific features

########################################################################################################
###################################### HUMAN CLASS #####################################################
########################################################################################################

#' Structure of the S4 class "Human"
#' 
#' Structure of the S4 class \code{Human} inheriting from class \code{\link{Organism-class}} representing human cells.
#' @export Human
#' @exportClass Human
#' @rdname Human
#'
#' @slot objective A character vector representing the current reaction which should be used as an objective function for the flux balance analysis.
setClass("Human",
         contains="Organism",
         representation(
           objective="character"
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Human-class}}
#' 
#' @name Human-constructor
#' @export
#' 
#' @param model model
#' @param objective A character vector representing the current reaction which should be used as an objective function for the flux balance analysis.
#' @param speed A integer vector representing the speed by which bacterium is moving (given by cell per iteration).
#' @param ... Arguments of \code{\link{Organism}}
#' @return Object of class \code{\link{Human-class}}
Human <- function(model, objective=model@react_id[which(model@obj_coef!=0)], speed=0, ...){
  model <- sybil::changeObjFunc(model, objective, obj_coef=model@obj_coef[which(model@obj_coef!=0)])
  new("Human", Organism(model=model, speed=speed, ...), objective=objective)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("objective", function(object){standardGeneric("objective")})
setMethod("objective", "Human", function(object){return(object@objective)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for changing the objective function of the model
#'
#' @description The generic function \code{changeFobj} changes the objective function, which is used for the linear programming in \code{optimizeLP}.
#' @export
#' @rdname changeFobj
#'
#' @param object An object of class Human.
#' @param new_fobj A character vector giving the reaction name of the new objective function.
#' @param model The original model structure which is converted into a problem object used for the next optimization.
#' @param alg A character vector giving the algorithm which should be used for the optimization (default is flux balance analysis).
#' @details To avoid the bias to just one particular objective function, the objective can be changed dynamically in this function. 
#' @seealso \code{\link{Human-class}} and \code{\link{optimizeLP}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' human <- Human(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' changeFobj(human,'EX_glc(e)',Ec_core)
setGeneric("changeFobj", function(object, new_fobj, model, alg="fba"){standardGeneric("changeFobj")})
#' @export
#' @rdname changeFobj
setMethod("changeFobj", "Human", function(object, new_fobj, model, alg="fba"){
  eval.parent(substitute(object@objective <- new_fobj)) #(pseudo) call by reference implementation
  model <- sybil::changeObjFunc(model, new_fobj)
  eval.parent(substitute(object@lpobj <- sysBiolAlg(model, algorithm=alg))) #the lp object has to be updated according to the new objective
})

#' @title Function implementing a growth model of a human cell
#'
#' @description The generic function \code{cellgrowth} implements different growth models for an object of class Human.
#' @export
#' @rdname cellgrowth
#'
#' @param object An object of class Human.
#' @param population An object of class Arena.
#' @param j The number of the iteration of interest.
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @param tstep A number giving the time intervals for each simulation step.
#' @return Boolean variable of the \code{j}th individual indicating if individual died.
#' @details Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value. Exponential growth of organisms is implemented by adding the calculated growthrate multiplied with the current growth calculated by \code{optimizeLP} plus to the already present growth value.
#' @seealso \code{\link{Human-class}}, \code{\link{growLin}} and \code{\link{growExp}}
setGeneric("cellgrowth", function(object, population, j, occupyM, fbasol, tstep){standardGeneric("cellgrowth")})
#' @export
#' @rdname cellgrowth
setMethod("cellgrowth", "Human", function(object, population, j, occupyM, fbasol, tstep){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$biomass <- growLin(object, popvec$biomass, fbasol, tstep)},
         "exponential"= {popvec$biomass <- growExp(object, popvec$biomass, fbasol, tstep)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'biomass'] <- popvec$biomass
  if(popvec$biomass > object@maxweight){
    freenb <- emptyHood(object, population@orgdat[,c('x','y')],
                        population@n, population@m, popvec$x, popvec$y, occupyM)
    if(length(freenb) != 0){
      npos = freenb[sample(length(freenb),1)]
      npos = as.numeric(unlist(strsplit(npos,'_')))
      if(occupyM[npos[2], npos[1]] == 0){ # check if there is no obstacle
        daughter <- popvec
        daughter$biomass <- popvec$biomass/2
        daughter$x <- npos[1]
        daughter$y <- npos[2]
        neworgdat[nrow(neworgdat)+1,] <- daughter
        neworgdat[j,'biomass'] <- popvec$biomass/2
      }
    }
  }else if(popvec$biomass < object@minweight){
    neworgdat[j,'biomass'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

#' @title Function for one simulation iteration for objects of Human class
#'
#' @description The generic function \code{simHum} implements all neccessary functions for the individuals to update the complete environment. 
#' @export
#' @rdname simHum
#'
#' @param object An object of class Human.
#' @param j The number of the iteration of interest.
#' @param arena An object of class Arena defining the environment.
#' @param bacnum integer indicating the number of bacteria individuals per gridcell
#' @param sublb A vector containing the substance concentrations in the current position of the individual of interest.
#' @param cutoff value used to define numeric accuracy.
#' @param pcut A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @param sec_obj character giving the secondary objective for a bi-level LP if wanted.
#' @param with_shadow True if shadow cost should be stores (default off).
#' @return Returns the updated enivironment of the \code{arena} parameter with all new positions of individuals on the grid and all new substrate concentrations.
#' @details Human cell individuals undergo the step by step the following procedures: First the individuals are constrained with \code{constrain} to the substrate environment, then flux balance analysis is computed with \code{optimizeLP}, after this the substrate concentrations are updated with \code{consume}, then the cell growth is implemented with \code{cellgrowth}, the potential new phenotypes are added with \code{checkPhen}, finally the conditional function \code{lysis} is performed. Can be used as a wrapper for all important cell functions in a function similar to \code{simEnv}.
#' @seealso \code{\link{Human-class}}, \code{\link{Arena-class}}, \code{\link{simEnv}}, \code{constrain}, \code{optimizeLP}, \code{consume}, \code{cellgrowth}, \code{checkPhen} and \code{lysis}
#' @examples
#' NULL
setGeneric("simHum", function(object, arena, j, sublb, bacnum, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE){standardGeneric("simHum")})
#' @export
#' @rdname simHum
setMethod("simHum", "Human", function(object, arena, j, sublb, bacnum, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE){
  const <- constrain(object, object@medium, lb=-sublb[j,object@medium]/bacnum, #scale to population size
                     dryweight=arena@orgdat[j,"biomass"], tstep=arena@tstep, scale=arena@scale, j)
  lobnd <- const[[1]]; upbnd <- const[[2]]
  optimization <- optimizeLP(object, lb=lobnd, ub=upbnd, j=j, sec_obj=sec_obj, cutoff=cutoff, with_shadow=with_shadow)
  fbasol <- optimization[[1]]
  
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,], bacnum=bacnum, fbasol=fbasol, cutoff) )) #scale consumption to the number of cells?
  
  dead <- cellgrowth(object, arena, j, arena@occupyM, fbasol=fbasol, tstep=arena@tstep)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object, fbasol=fbasol, cutoff=pcut))
  type <- object@type
  arena@mflux[[type]]  <- arena@mflux[[type]] + fbasol$fluxes # remember active fluxes
  arena@shadow[[type]] <- arena@shadow[[type]]+ optimization[[2]]

  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, names(arena@media), sublb[j,])))
  }
  pos <- arena@orgdat[,c('x','y')]
  if(!dead && !arena@stir && object@speed != 0){
    if(object@chem == ''){
      mov_pos <- move(object, pos, arena@n, arena@m, j, arena@occupyM)
      arena@orgdat[,c('x','y')] <- mov_pos
    }else{
      chemo_pos <- chemotaxis(object, arena, j, arena@occupyM)
      arena@orgdat[j,c('x','y')] <- chemo_pos
    }
  }
  return(arena)
})

#show function for class Human

setMethod(show, signature(object="Human"), function(object){
  print(paste('Cell culture ',object@type,' of class Human.',sep=''))
})
