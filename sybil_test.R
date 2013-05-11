library(sybil)
library(lpSolveAPI)

SYBIL_SETTINGS("SOLVER", "lpSolveAPI")
SYBIL_SETTINGS("METHOD", "lp_solve")
SYBIL_SETTINGS()

source("http://bioconductor.org/biocLite.R")
biocLite("rsbml")

test <- optObj(solver = "lpSolveAPI") 

data(Ec_core)
optL <- optimizeProb(Ec_core, algorithm = "fba", retOptSol = FALSE)
optL <- optimizeProb(Ec_core)

ex <- findExchReact(Ec_core)
upt <- uptReact(ex)
ex[upt]
mod <- changeBounds(Ec_core, ex[c("EX_glc(e)", "EX_lac_D(e)")], lb = c(0, -10))
findExchReact(mod)
optL <- optimizeProb(Ec_core, algorithm = "fba", retOptSol = FALSE)

optimizeProb(Ec_core)
