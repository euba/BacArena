Bcoli_sbml <- read.sbml("data/ecoli_iAF1260.xml", "_b$")
# last delivered variable is regexpr for filtering external reactions (special handling necessary for equilibrium condition in fba)

Bcoli_biomassf <- "R_Ec_biomass_iAF1260_core_59p81M"

Bcoli_maintenancef <- "R_ATPM"

Bcoli_stochmatrix <- Bcoli_sbml$stoch

#colnames(Bcoli_stochmatrix)[grep("biomass", colnames(Bcoli_stochmatrix))] #for finding specific reactions
colnames(Bcoli_stochmatrix)[grep("R_ATP", colnames(Bcoli_stochmatrix))] #for finding specific reactions

#
# associate each substrate with an exchange reaction (sbml specific!!)
#
sname <- c("acetate","aketoglutarate", "co2", "ethanol", "formiate", "fumarate","glucose", "h2o", "proton", "lactate","o2", "iphosphate", "pyruvate", "succinate")
Bcoli_sub_ex <- c("R_EX_ac_LPAREN_e_RPAREN_", "R_EX_akg_LPAREN_e_RPAREN_", "R_EX_co2_LPAREN_e_RPAREN_", "R_EX_etoh_LPAREN_e_RPAREN_", "R_EX_for_LPAREN_e_RPAREN_", "R_EX_fum_LPAREN_e_RPAREN_",  "R_EX_glc_LPAREN_e_RPAREN_", "R_EX_h2o_LPAREN_e_RPAREN_", "R_EX_h_LPAREN_e_RPAREN_", "R_EX_lac_D_LPAREN_e_RPAREN_", "R_EX_o2_LPAREN_e_RPAREN_", "R_EX_pi_LPAREN_e_RPAREN_", "R_EX_pyr_LPAREN_e_RPAREN_", "R_EX_succ_LPAREN_e_RPAREN_")

names(Bcoli_sub_ex) <- sname

Bcoli_lower_bound <- Bcoli_sbml$lb
Bcoli_upper_bound <- Bcoli_sbml$ub

#
# set lower bounds to current substrat concentration in cell
#
#Varma and Palsson 1994:
#maximum oxygen utilization rate (15 mmol of 02 per g [dry weight] per h)
#the maximum aerobic glucose utilization rate (10.5 mmol of Glc per g [dryweight] per h), the maximum anaerobic glucose utilization rate (18.5 mmol of Glc per g [dry weight] per h),  
#the non-growth-associated maintenance requirements (7.6 mmol of ATP per g [dry weight] per h), and the growth-associated maintenance requirements (13 mmol of ATP per g of biomass).

Bcoli_set_lower_bound <- function(substrat){
  #
  # lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  #
  Bcoli_lower_bound <- Bcoli_sbml$lb
  Bcoli_upper_bound <- Bcoli_sbml$ub
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)==get_maintenancef("Bcoli"))] <- 7.6
  Bcoli_upper_bound[which(colnames(Bcoli_stochmatrix)==get_maintenancef("Bcoli"))] <- 7.6
  Bcoli_lower_bound[grep("R_EX", colnames(Bcoli_stochmatrix))] <- 0 # define growth media
  
  #set minimal medium -> see suppl. Table2 material Feist et al 2007
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_ca2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_cl_LPAREN_e_RPAREN_")] = -100
  #Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_co2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_cobalt2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_cu2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_fe2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_fe3_LPAREN_e_RPAREN_")] = -100
  #Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_h_LPAREN_e_RPAREN_")] = -100
  #Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_h2o_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_k_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_mg2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_mn2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_mobd_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_na1_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_nh4_LPAREN_e_RPAREN_")] = -100
  #Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_pi_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_so4_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_tungs_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_zn2_LPAREN_e_RPAREN_")] = -100
  Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)=="R_EX_cbl1_LPAREN_e_RPAREN_")] = -100
  
  #print(substrat)
  sapply(names(sapply(substrat, names)), function(x, Bcoli_stochmatrix, Bcoli_sub_ex, Bcoli_lower_bound){
    #print(x)
    if(x %in% names(Bcoli_sub_ex)) Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)==Bcoli_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
  },Bcoli_stochmatrix=Bcoli_stochmatrix, Bcoli_sub_ex=Bcoli_sub_ex, Bcoli_lower_bound=Bcoli_lower_bound)
  
  
  #limit flux exchange rate:
  if(substrat[["o2"]] < 0.01){
    if(substrat[["glucose"]] > 18.5) Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)==Bcoli_sub_ex[["glucose"]])] <- -18.5
  }else{
    if(substrat[["glucose"]] > 10.5) Bcoli_lower_bound[which(colnames(Bcoli_stochmatrix)==Bcoli_sub_ex[["glucose"]])] <- -10.5
  }
  return(Bcoli_lower_bound)
}
