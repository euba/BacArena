mod <- readSBMLmod("data/ecoli_core.xml", bndCond = FALSE)

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
ecoli_sub_ex <- react_id(findExchReact(mod))
names(ecoli_sub_ex) <- sname

#
# set lower bounds to current substrat concentration in cell
#
#Varma and Palsson 1994:
#maximum oxygen utilization rate (15 mmol of 02 per g [dry weight] per h)
#the maximum aerobic glucose utilization rate (10.5 mmol of Glc per g [dryweight] per h), the maximum anaerobic glucose utilization rate (18.5 mmol of Glc per g [dry weight] per h),  
#the non-growth-associated maintenance requirements (7.6 mmol of ATP per g [dry weight] per h), and the growth-associated maintenance requirements (13 mmol of ATP per g of biomass).

ecoli_set_lower_bound <- function(substrat){ 
  #lowbnd(mod) = rep(0, length(lowbnd(mod))) #define growth media
  sapply(names(substrat), function(x, y, z, substrat, ecoli_sub_ex){
    rxn <- ecoli_sub_ex[x]
    if(!is.na(rxn)){
      mod <<- changeBounds(mod, rxn, lb=-substrat[[x]][y,z])
    }
  }, substrat=substrat, y=1, z=1, ecoli_sub_ex=ecoli_sub_ex)
  mod <- changeBounds(mod, c(ecoli_sub_ex["glucose"], ecoli_sub_ex["o2"]), lb=-c(11, 18.2))
  
  #
  # lower/upper bound is init in read.sbml because of irreversible reactions => lb=0
  #
  #ecoli_lower_bound <- ecoli_sbml$lb
  #ecoli_upper_bound <- ecoli_sbml$ub
  #ecoli_lower_bound[grep("R_EX", colnames(ecoli_stochmatrix))] <- 0 # define growth media
  
  #print(substrat)
  #sapply(names(sapply(substrat, names)), function(x, ecoli_stochmatrix, ecoli_sub_ex, ecoli_lower_bound){
    #print(x)
  #  if(x %in% names(ecoli_sub_ex)) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[[x]])] <<- - substrat[[x]] # "<<-" is necessary for extern variable modification
  #},ecoli_stochmatrix=ecoli_stochmatrix, ecoli_sub_ex=ecoli_sub_ex, ecoli_lower_bound=ecoli_lower_bound)
  
  
  #limit flux exchange rate:
  #if(substrat[["o2"]] < 0.01){
  #  if(substrat[["glucose"]] > 18.5) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -18.5
  #}else{
  #  if(substrat[["glucose"]] > 10.5) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -10.5
  #}
  
  #Feist et al 2007:
  #if(substrat[["glucose"]] > 11) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["glucose"]])] <- -11
  #if(substrat[["o2"]] > 18.2) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["o2"]])] <- -18.2
  #Orth et al 2011:
  #if(substrat[["lactate"]] > 16) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["lactate"]])] <- -16
  #if(substrat[["succinate"]] > 16) ecoli_lower_bound[which(colnames(ecoli_stochmatrix)==ecoli_sub_ex[["succinate"]])] <- -16
  #return(ecoli_lower_bound)
}

optimizeProb(mod, algorithm = "fba", retOptSol = F, solver = "clpAPI")$obj #test objective
