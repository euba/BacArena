Ec_core@obj_coef

library(MASS)
IC = ginv(as.matrix(t(Ec_core@S)))


sc = IC %*% ( Ec_core@obj_coef -  )

getRedCosts(fbasl)
redC <- getRedCosts(lpob@problem)

names(sc) = Ec_core@met_id

