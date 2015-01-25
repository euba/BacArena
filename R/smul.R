if(!exists("smul_sbml", mode="list")) smul_sbml <- read.sbml("/home/johannes/sulfuro/src/sbml-sim/sm.xml", "_boundary$") #_b$  regexpr for filtering external species (special handling necessary for equilibrium condition in fba)

smul_biomassf <- "R_Biomass"

smul_maintenancef <- "R_ATPM"

smul_stochmatrix <- smul_sbml$stoch

#
# associate each substrate with an exchange reaction (sbml specific!!)
#
sname <- c("acetate", "co2", "h2o", "iphosphate", "pyruvate", "h2", "o2")
smul_sub_ex <- c("R_EX_ACET_p", "R_EX_CARBONDIOXIDE_c", "R_EX_WATER_c", "R_EX_Pi_p", "R_EX_PYRUVATE_c", "R_EX_HYDROGENMOLECULE_p", "R_EX_OXYGENMOLECULE_c")
names(smul_sub_ex) <- sname

smul_lower_bound <- smul_sbml$lb
smul_upper_bound <- smul_sbml$ub

smul_ngam <- 8.39 #
smul_gam <- 45.73 #from biomass function

smul_set_lower_bound <- function(substrat){
  #
  # lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  #
  smul_lower_bound <- smul_sbml$lb
  smul_upper_bound <- smul_sbml$ub
  smul_lower_bound[which(colnames(smul_stochmatrix)==get_maintenancef("smul"))] <- smul_ngam
  smul_upper_bound[which(colnames(smul_stochmatrix)==get_maintenancef("smul"))] <- smul_ngam
  smul_lower_bound[grep("R_EX", colnames(smul_stochmatrix))] <- 0 # define growth media
  
  #set minimal medium
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_ANTHRANILATE_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="EX_WATER_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_Pi_p")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_PANTOTHENATE_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_THIAMINEP_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_NA_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_FE2_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_FE3_c")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_SULFATE_p")] = -100
  smul_lower_bound[which(colnames(smul_stochmatrix)=="R_EX_AMMONIUM_c")] = -100
  
  
  #print(substrat)
  sapply(names(sapply(substrat, names)), function(x, smul_stochmatrix, smul_sub_ex, smul_lower_bound){
    if(x %in% names(smul_sub_ex)) {
      #print(paste(x, -substrat[[x]]))
      smul_lower_bound[which(colnames(smul_stochmatrix)==smul_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
    }
  },smul_stochmatrix=smul_stochmatrix, smul_sub_ex=smul_sub_ex, smul_lower_bound=smul_lower_bound)
  
  #max uptake rates:
  #if(substrat[["pyruvate"]] > 20) smul_lower_bound[which(colnames(smul_stochmatrix)==smul_sub_ex[["pyruvate"]])] <- -20
  #if(substrat[["o2"]] > 2.5) smul_lower_bound[which(colnames(smul_stochmatrix)==smul_sub_ex[["o2"]])] <- -2.5
    
  #print(smul_lower_bound)
  
  #slb <- smul_lower_bound
  #names(slb) <- colnames(smul_stochmatrix)
  #print(slb[which(colnames(smul_stochmatrix) %in% smul_sub_ex)])
  return(smul_lower_bound)
}
