barkeri_sbml <- read.sbml("data/barkeri_iAF692.xml", "_b$") #_b$  regexpr for filtering external species (special handling necessary for equilibrium condition in fba)

barkeri_biomassf <- "R_Mb_biomass_30"

barkeri_maintenancef <- "R_ATPM"

barkeri_stochmatrix <- barkeri_sbml$stoch

#
# associate each substrate with an exchange reaction (sbml specific!!)
#
sname <- c("acetate", "co2", "h2o", "proton", "iphosphate", "pyruvate", "h2", "methanol", "methane")
barkeri_sub_ex <- c("R_EX_ac_LPAREN_e_RPAREN_", "R_EX_co2_LPAREN_e_RPAREN_", "R_EX_h2o_LPAREN_e_RPAREN_", "R_EX_h_LPAREN_e_RPAREN_", "R_EX_pi_LPAREN_e_RPAREN_", "R_EX_pyr_LPAREN_e_RPAREN_", "R_EX_h2_LPAREN_e_RPAREN_", "R_EX_meoh_LPAREN_e_RPAREN_", "R_EX_ch4_LPAREN_e_RPAREN_")
names(barkeri_sub_ex) <- sname

barkeri_lower_bound <- barkeri_sbml$lb
barkeri_upper_bound <- barkeri_sbml$ub


barkeri_set_lower_bound <- function(substrat){
  #
  # lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  #
  barkeri_lower_bound <- barkeri_sbml$lb
  barkeri_upper_bound <- barkeri_sbml$ub
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)==get_maintenancef("barkeri"))] <- 1.75 # feist paper 2013
  barkeri_upper_bound[which(colnames(barkeri_stochmatrix)==get_maintenancef("barkeri"))] <- 1.75
  barkeri_lower_bound[grep("R_EX", colnames(barkeri_stochmatrix))] <- 0 # define growth media

  #set minimal medium
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_4abz_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_cobalt2_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_h2s_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_n2_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_na1_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_nac_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_nh4_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_ni2_LPAREN_e_RPAREN_")] = -100
  barkeri_lower_bound[which(colnames(barkeri_stochmatrix)=="R_EX_so3_LPAREN_e_RPAREN_")] = -100
    
  #print(substrat)
  sapply(names(sapply(substrat, names)), function(x, barkeri_stochmatrix, barkeri_sub_ex, barkeri_lower_bound){
    if(x %in% names(barkeri_sub_ex)) barkeri_lower_bound[which(colnames(barkeri_stochmatrix)==barkeri_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
  },barkeri_stochmatrix=barkeri_stochmatrix, barkeri_sub_ex=barkeri_sub_ex, barkeri_lower_bound=barkeri_lower_bound)
  #print(barkeri_lower_bound)
  return(barkeri_lower_bound)
}
