if(!exists("beijerinckii_sbml", mode="list")) beijerinckii_sbml <- read.sbml("data/Cbeijerinckii_iCM925.xml", "_b$") #_b$  regexpr for filtering external species (special handling necessary for equilibrium condition in fba)

beijerinckii_biomassf <- "R_biomass"

beijerinckii_maintenancef <- "R_ATPM"

beijerinckii_stochmatrix <- beijerinckii_sbml$stoch

#
# associate each substrate with an exchange reaction (sbml specific!!)
#
sname <- c("acetate", "co2", "ethanol", "glucose", "h2o", "proton", "lactate", "iphosphate", "succinate", "h2", "acetone", "butyrate", "butanol")
beijerinckii_sub_ex <- c("R_ex_ac_e", "R_ex_co2_e", "R_ex_etoh_e", "R_ex_glc_D_e", "R_ex_h2o_e", "R_ex_h_e", "R_ex_lac_D_e", "R_ex_pi_e", "R_ex_succ_e", "R_ex_h2_e", "R_ex_acetone_e", "R_ex_but_e", "R_ex_butoh_e")
names(beijerinckii_sub_ex) <- sname

beijerinckii_lower_bound <- beijerinckii_sbml$lb
beijerinckii_upper_bound <- beijerinckii_sbml$ub

beijerinckii_ngam <- 8.5 #Price Paper 2006
beijerinckii_gam <- 40 #from biomass function


beijerinckii_set_lower_bound <- function(substrat){
  #
  # lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  #
  beijerinckii_lower_bound <- beijerinckii_sbml$lb
  beijerinckii_upper_bound <- beijerinckii_sbml$ub
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)==get_maintenancef("beijerinckii"))] <- 8.5 # milne paper 2011
  beijerinckii_upper_bound[which(colnames(beijerinckii_stochmatrix)==get_maintenancef("beijerinckii"))] <- 8.5
  beijerinckii_lower_bound[grep("R_EX", colnames(beijerinckii_stochmatrix))] <- 0 # define growth media
  
  #set minimal medium milne paper 2011
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_4abz_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_ac_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_acetone_e")] = -100
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_btn_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_but_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_butoh_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_co2_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_etoh_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_glc-D_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_h_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_h2_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_h2o_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_lac-D_e")] = -100
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_lac-L_e")] = -100
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_nh4_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_pi_e")] = -100
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_so4_e")] = -100
  beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_thm_e")] = -100
  #beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)=="R_ex_succ_e")] = -100
  
      
  #print(substrat)
  sapply(names(sapply(substrat, names)), function(x, beijerinckii_stochmatrix, beijerinckii_sub_ex, beijerinckii_lower_bound){
    if(x %in% names(beijerinckii_sub_ex)) beijerinckii_lower_bound[which(colnames(beijerinckii_stochmatrix)==beijerinckii_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
  },beijerinckii_stochmatrix=beijerinckii_stochmatrix, beijerinckii_sub_ex=beijerinckii_sub_ex, beijerinckii_lower_bound=beijerinckii_lower_bound)
  #print(beijerinckii_lower_bound)
  return(beijerinckii_lower_bound)
}