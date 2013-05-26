library(lpSolveAPI)
m <- matrix(c(1,-1,-1,0,0,0,0,0,0, 0,1,0,-1,0,-2,0,0,0, 0,0,1,0,0,1,-1,0,0, 0,0,0,2,0,0,1,-1,0, 0,0,0,-1,1,0,0,0,0, 0,0,0,0,0,1,0,0,-1), nrow=6, ncol=9, byrow=T)
model <- make.lp(0, 9)
lp.control(model,sense='max')
c <- rep(0,9)
c[8] <- 1
set.objfn(model, c)
for(i in 1:nrow(m)){
  add.constraint(model, m[i,], "=", 0)
}
set.bounds(model, lower=rep(0,9))
set.bounds(model, upper=10, columns=1)
#ColNames <- c("r1","r2","r3","r4","r5","r6","r7","r8","r9")
#RowNames <- c("a","b","c","d","e","f","g")
#dimnames(model) <- list(RowNames, ColNames)
solve(model)
get.objective(model)

write.lp(model,"test", NULL)




lps.model <- make.lp(0, 3)
xt <- c(6,2,4)
add.constraint(lps.model, xt, "<=", 150)
xt <- c(1,1,6)
add.constraint(lps.model, xt, ">=", 0)
xt <- c(4,5,4)
add.constraint(lps.model, xt, "=", 40)
set.objfn(lps.model, c(-3,-4,-3))

solve(lps.model)
get.objective(lps.model)






set.objfn(linp, stoch[,which(colnames(stoch)=="R_Biomass_Ecoli_core_N__w_GAM_")])
add.constraint(linp, stoch[1,])
set.bounds(linp,lower=lb)
set.bounds(linp,upper=ub)
dimnames(linp) <- 
  solve(linp)
get.objective(linp)  

setpwd()



f.con <- t(matrix(c(1,-1,-1,0,0,0,0,0,0, 0,1,0,-1,0,-2,0,0,0, 0,0,1,0,0,1,-1,0,0, 0,0,0,2,0,0,1,-1,0, 0,0,0,-1,1,0,0,0,0, 0,0,0,0,0,1,0,0,-1, 1,-1,-1,0,0,0,0,0,0), nrow=7, ncol=9, byrow=T))
c <- rep(0,9)
c[8] <- 1
f.obj <- c
f.dir <- c(rep(">=", 9), "<=")
f.rhs <- c(rep(0,9),10)
lp ("max", f.obj, f.con, f.dir, f.rhs)


f.obj <- c(1, 9, 3)
f.con <- matrix (c(1, 2, 3, 3, 2, 2), nrow=2, byrow=TRUE)
f.dir <- c("<=", "<=")
f.rhs <- c(9, 15)
#
# Now run.
#
lp ("max", f.obj, f.con, f.dir, f.rhs)








# linear programming
c <- rep(0, dim(stoch)[2])
c[which(colnames(stoch)=="R_Biomass_Ecoli_core_N__w_GAM_")] <- 1 
#f.obj <- stoch[,which(colnames(stoch)=="R_Biomass_Ecoli_core_N__w_GAM_")]
f.obj <- c
f.con <- rbind(stoch,stoch)
f.dir <- c(rep(">=", dim(stoch)[1]), rep("<=", dim(stoch)[1]))
f.rhs <- c(lb, ub)
lp ('max', f.obj, f.con, f.dir, f.rhs)