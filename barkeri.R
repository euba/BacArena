barkeri_sbml <- read.sbml("data/barkeri_iAF692.xml", "_b$")

barkeri_biomassf <- "R_Mb_biomass_30"

barkeri_maintenancef <- "R_ATPM"

barkeri_stochmatrix <- barkeri_sbml$stoch

barkeri_ex_pattern <- "_b$"  # regexpr for filtering external reactions (special handling necessary for equilibrium condition in fba)

#
# associate each substrate with an exchange reaction (sbml specific!!)
#
barkeri_sub_ex <- character(length(s))
names(barkeri_sub_ex) <- s
barkeri_sub_ex[["acetate"]]         <- "R_EX_ac_LPAREN_e_RPAREN_"
#ecoli_sub_ex[["aketoglutarate"]]  <- "R_EX_akg_e_"
barkeri_sub_ex[["co2"]]             <- "R_EX_co2_LPAREN_e_RPAREN_"
#ecoli_sub_ex[["etanol"]]          <- "R_EX_etoh_e_"
#ecoli_sub_ex[["fumarate"]]        <- "R_EX_fum_e_"
#ecoli_sub_ex[["formiate"]]        <- "R_EX_for_e_"
#ecoli_sub_ex[["glucose"]]         <- "R_EX_glc_e_"
barkeri_sub_ex[["h2o"]]             <- "R_EX_h2o_LPAREN_e_RPAREN_"
barkeri_sub_ex[["proton"]]          <- "R_EX_h_LPAREN_e_RPAREN_"
#ecoli_sub_ex[["lactate"]]         <- "R_EX_lac_D_e_"
#ecoli_sub_ex[["o2"]]              <- "R_EX_o2_e_"
barkeri_sub_ex[["iphosphate"]]      <- "R_EX_pi_LPAREN_e_RPAREN_"
barkeri_sub_ex[["pyruvate"]]        <- "R_EX_pyr_LPAREN_e_RPAREN_"
#ecoli_sub_ex[["succinate"]]       <- "R_EX_succ_e_"

barkeri_sub_ex[["h2"]]              <- "R_EX_h2_LPAREN_e_RPAREN_"
barkeri_sub_ex[["methanol"]]        <- "R_EX_meoh_LPAREN_e_RPAREN_"
barkeri_sub_ex[["methane"]]         <-  "R_EX_ch4_LPAREN_e_RPAREN_"