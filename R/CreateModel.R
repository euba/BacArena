########################################################################################################
############################### Create new species models in BacArena ##################################
########################################################################################################
### NOTE: The names of the reactions for the metabolic models should be in BiGG nomenclature!
### EB 14/12/13: I separated functions for model creation and constraint definition, because sometimes
### contraints need to be redefined and it would be better to not load the sbml every time

create_model <- function(name, sbml, biomass, maintenance="R_ATPM",
                         mets, translate, ngam, gam){ #stores sbml data in BacArena specific data type
  #reading model in
  print("Reading SBML model...")
  model <- read.sbml(sbml, "_b$")
  
  #setting biomass function
  biomassf <- biomass
  #setting ATP maintenance
  maintenancef <- maintenance
  #setting reaction names
  sub_ex <- mets
  sname <- translate
  names(sub_ex) <- sname
  
  #setting non-growth and growth associated maintenance
  model_ngam <- ngam
  model_gam <- gam
  
  print("Model successfully created!")
  return(list(sbml=model, biomassf=biomassf, maintenancef=maintenancef,
              ngam=model_ngam, gam=model_gam))
}

constrain <- function(name, model, minimal_medium, contraints, substrat){ #add model contraints (from literature or own)
  #storing stochiometric matrix
  stochmatrix <- model$sbml$stoch
  
  #lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  lower_bound <- model$sbml$lb
  upper_bound <- model$sbml$ub
  lower_bound[which(colnames(stochmatrix)==get_maintenancef(name))] <- model$ngam # milne paper 2011
  upper_bound[which(colnames(stochmatrix)==get_maintenancef(name))] <- model$ngam
  lower_bound[grep("R_EX", colnames(stochmatrix))] <- 0 # define growth media
  
  #set lower bounds based on given minimal medium
  apply(minimal_medium, 1, function(x, lower_bound, stochmatrix){
    lower_bound[which(colnames(stochmatrix)==x[1])] <<- x[2]
  }, lower_bound=lower_bound, stochmatrix=stochmatrix)
  
  #set lower bound based on current substrat concentration -> not a good solution, maybe put in different function?
  sapply(names(sapply(substrat, names)), function(x, beijerinckii_stochmatrix, beijerinckii_sub_ex, beijerinckii_lower_bound){
    if(x %in% names(beijerinckii_sub_ex)) beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)==beijerinckii_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
  },beijerinckii_stochmatrix=beijerinckii_stochmatrix, beijerinckii_sub_ex=beijerinckii_sub_ex, beijerinckii_lower_bound=beijerinckii_lower_bound)
  
  #set lower bounds based on given additional constraints -> what about conditional contraints (if o2, then...)
  apply(constraints, 1, function(x, lower_bound, stochmatrix){
    lower_bound[which(colnames(stochmatrix)==x[1])] <<- x[2]
  }, lower_bound=lower_bound, stochmatrix=stochmatrix)
  
  return(lower_bound)
}