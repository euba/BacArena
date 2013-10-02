ecoli_big_sbml <- read.sbml("data/ecoli_iAF1260.xml", "_b$")
# last delivered variable is regexpr for filtering external reactions (special handling necessary for equilibrium condition in fba)

ecoli_big_biomassf <- "R_Biomass_Ecoli_core_N__w_GAM_"

ecoli_big_maintenancef <- "R_ATPM"

ecoli_big_stochmatrix <- ecoli_big_sbml$stoch
