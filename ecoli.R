ecoli_sbml <- read.sbml("data/ecoli_core.xml", "_b$")
# last delivered variable is regexpr for filtering external reactions (special handling necessary for equilibrium condition in fba)

ecoli_biomassf <- "R_Biomass_Ecoli_core_N__w_GAM_"

ecoli_maintenancef <- "R_ATPM"

ecoli_stochmatrix <- ecoli_sbml$stoch

#
# associate each substrate with an exchange reaction (sbml specific!!)
#
#ecoli_sub_ex <- character(length(s))
#names(ecoli_sub_ex) <- s
#ecoli_sub_ex[["acetate"]]         <- "R_EX_ac_e_"
#ecoli_sub_ex[["aketoglutarate"]]  <- "R_EX_akg_e_"
#ecoli_sub_ex[["co2"]]             <- "R_EX_co2_e_"
#ecoli_sub_ex[["ethanol"]]          <- "R_EX_etoh_e_"
#ecoli_sub_ex[["fumarate"]]        <- "R_EX_fum_e_"
#ecoli_sub_ex[["formiate"]]        <- "R_EX_for_e_"
#ecoli_sub_ex[["glucose"]]         <- "R_EX_glc_e_"
#ecoli_sub_ex[["h2o"]]             <- "R_EX_h2o_e_"
#ecoli_sub_ex[["proton"]]          <- "R_EX_h_e_"
#ecoli_sub_ex[["lactate"]]         <- "R_EX_lac_D_e_"
#ecoli_sub_ex[["o2"]]              <- "R_EX_o2_e_"
#ecoli_sub_ex[["iphosphate"]]      <- "R_EX_pi_e_"
#ecoli_sub_ex[["pyruvate"]]        <- "R_EX_pyr_e_"
#ecoli_sub_ex[["succinate"]]       <- "R_EX_succ_e_"

sname <- c("acetate","aketoglutarate", "co2", "ethanol", "formiate", "fumarate","glucose", "h2o", "proton", "lactate","o2", "iphosphate", "pyruvate", "succinate")
ecoli_sub_ex <- c("R_EX_ac_e_", "R_EX_akg_e_", "R_EX_co2_e_", "R_EX_etoh_e_", "R_EX_for_e_", "R_EX_fum_e_",  "R_EX_glc_e_", "R_EX_h2o_e_", "R_EX_h_e_", "R_EX_lac_D_e_", "R_EX_o2_e_", "R_EX_pi_e_", "R_EX_pyr_e_", "R_EX_succ_e_")
names(ecoli_sub_ex) <- sname

ecoli_lower_bound <- ecoli_sbml$lb
ecoli_upper_bound <- ecoli_sbml$ub

#
# set lower bounds to current substrat concentration in cell
#
#Varma and Palsson 1994:
#maximum oxygen utilization rate (15 mmol of 02 per g [dry weight] per h)
#the maximum aerobic glucose utilization rate (10.5 mmol of Glc per g [dryweight] per h), the maximum anaerobic glucose utilization rate (18.5 mmol of Glc per g [dry weight] per h),  
#the non-growth-associated maintenance requirements (7.6 mmol of ATP per g [dry weight] per h), and the growth-associated maintenance requirements (13 mmol of ATP per g of biomass).

ecoli_set_lower_bound <- function(substrat){
  #
  # lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  #
  ecoli_lower_bound <- ecoli_sbml$lb
  ecoli_upper_bound <- ecoli_sbml$ub
  ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==get_maintenancef("ecoli"))] <- 7.6
  ecoli_upper_bound[which(colnames(ecoli_stochmatrix)==get_maintenancef("ecoli"))] <- 7.6
  ecoli_lower_bound[grep("R_EX", colnames(ecoli_stochmatrix))] <- 0 # define growth media
  
  
  sapply(names(sapply(substrat, names)), function(x, ecoli_stochmatrix, ecoli_sub_ex, lb){
    if(x %in% names(ecoli_sub_ex)) lb[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
  },ecoli_stochmatrix=ecoli_stochmatrix, ecoli_sub_ex=ecoli_sub_ex, lb=lb)
  
  
  #ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <-  -substrat[["glucose"]]
  #if anaerobic conditions overwrite maximal Glucose uptake:
  
  if(substrat[["o2"]] < 0.01){
    #if(substrat[["glucose"]] < 18.5) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -substrat[["glucose"]]
    #else ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -18.5
    if(substrat[["glucose"]] > 18.5) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -18.5
  }else{
    #if(substrat[["glucose"]] < 10.5) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -substrat[["glucose"]]
    #else ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -10.5
    if(substrat[["glucose"]] > 10.5) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -10.5
  }
  #ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -18.5
  ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["h2o"]])] <-  -substrat[["h2o"]]
  ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["proton"]])]   <-  -substrat[["proton"]]
  ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["o2"]])]  <-  -substrat[["o2"]]
  #if(substrat[["o2"]] < 15) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["o2"]])] <- -substrat[["o2"]]
  #  else ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["o2"]])] <- -15
  ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["iphosphate"]])]  <-  -substrat[["iphosphate"]]
  return(ecoli_lower_bound)
}
