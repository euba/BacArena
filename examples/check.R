# script to check basic functionality of bacarena
library(sybil)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(compiler)
library(ReacTran)
setwd('P:/GitRep/BacArena')
source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")

SYBIL_SETTINGS("SOLVER","cplexAPI") #setting solver to GUROBI


#
# 1) basic growth model ecoli core
#
openArena()


#
# 2a) Diffusion test
#
n=100; m=100; D=2; tsteps=8; L=10; init=1
arena <- Arena(n=n, m=m, Lx=L, Ly=L)
arena@media[[1]]  <- Substance(arena@n, arena@m, smax=0, gridgeometry=arena@gridgeometry, id="test", name="test substance", difspeed=D)
names(arena@media) <- "test"
arena@media[[1]]@diffmat[round(arena@n/2), 1] <- init
sim_diff <- simEnv(arena, time=tsteps, continue = TRUE)

par(mfrow=c(3,3))
lapply(sim_diff@medlist, function(x){image(matrix(x$test, nrow=arena@n, byrow=FALSE))})
par(mfrow=c(1,1))

outa <- lapply(seq_along(sim_diff@medlist), function(i){c(i, sim_diff@medlist[i])})
outa <- matrix(unlist(outa), byrow=TRUE, nrow=tsteps+1) # t_0
colnames(outa) <- c("time", paste(1:(n*m)))
attributes(outa)$class <- c("deSolve", "matrix")
attributes(outa)$dimens <- c(n,m)
attributes(outa)$nspec <- 1
image(outa, ask = FALSE, mfrow = c(3, 3), main = paste("time", times))

# 2b) Comparison to pde only
library(ReacTran)
x.grid  <- setup.grid.1D(x.up = 0, L = L, N = n)
y.grid  <- setup.grid.1D(x.up = 0, L = L, N = m)
grid2D <- setup.grid.2D(x.grid, y.grid)
D.grid    <- setup.prop.2D(value = D, y.value = D, grid = grid2D)
y0 <- matrix(nrow = n, ncol = m, data = 0)
y0[round(n/2), 1] <- 1
Diff2d_2 <- function (t, y, pars)  {
  with (as.list(pars), {
    CONC  <- matrix(nrow = n, ncol = m, data = y)
    dCONC <- tran.2D(CONC, grid = grid2D, D.grid = D.grid)$dC
    return (list(dCONC))
  })
}
outb <- ode.2D (y = y0, func = Diff2d_2, t = 1:tsteps+1, parms=NULL, # t_0
                dim = c(n, m), lrw = 1706236)
image(outb, ask = FALSE, mfrow = c(3, 3), main = paste("time", times))
