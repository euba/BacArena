ecoli_sbml <- read.sbml("data/ecoli_iAF1260.xml", "_b$")
# last delivered variable is regexpr for filtering external reactions (special handling necessary for equilibrium condition in fba)

ecoli_biomassf <- "R_Biomass_Ecoli_core_N__w_GAM_"

ecoli_maintenancef <- "R_ATPM"

ecoli_stochmatrix <- ecoli_sbml$stoch
