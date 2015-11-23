library(microbenchmark)

setwd("~/uni/bacarena")
#library(Rcpp)
#library(RcppArmadillo)
#Rcpp::sourceCpp("src/diff.cpp")
library(Matrix)
library(sybil)
library(ReacTran)
library(deSolve)
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")

data(Ec_core)
arena = Arena(n=11, m=11, stir=F, Lx=5, Ly=5)
bac = Bac(model=Ec_core, growtype="exponential",
          speed=0, deathrate=0, type="ecore", lyse=F)
bac = Bac(model=Ec_core, growtype="exponential", cellarea=4.42, lyse=F)
addOrg(arena, bac, amount=1, x=1, y=1)
addSubs(arena, smax=5, unit="mM", difunc="pde", difspeed=1)
changeSub(arena, smax=0, "EX_co2(e)")
#changeSub(arena,20,c("EX_glc(e)","EX_pi(e)", "EX_h2o(e)", "EX_h(e)", "EX_nh4(e"))
arena@media$`EX_co2(e)`@difspeed=0.5
arena@media$`EX_co2(e)`@boundS=0.1
arena@media$`EX_co2(e)`@pde="InfluxBoundDiff2d"
sim <- simEnv(arena, time=15, continue=T)
out <- lapply(sim@medlist, function(x){matrix(x$`EX_co2(e)`, nrow=arena@n)})
par(mfrow = c(ceiling(sqrt(length(out))), ceiling(sqrt(length(out)))))
out_min <- min(sapply(out, min))
out_max <- max(sapply(out, max))
for(i in 1:length(out)){image(out[[i]], main=paste("time",i), zlim=c(out_min, out_max))}
t(round(out[[3]],2))



min <- min(sapply(out, min))
max <- max(sapply(out, max))


par(mfrow=c(1,1))

arena@media$`EX_o2(e)`@diffmat
plotCurves2(sim)





arena = Arena(n=11, m=11, stir=F, Lx=10, Ly=5)
#show(arena)
data(Ec_core)
bac = Bac(model=Ec_core, growtype="exponential",
          speed=0, deathrate=0, type="ecore", lyse=F)
bac = Bac(model=Ec_core, growtype="exponential", cellarea=4.42, lyse=F, deathrate=0)
addOrg(arena, bac, amount=1, x=1, y=1)
addSubs(arena, smax=0, difunc="pde", difspeed=1)
#arena@media$`EX_co2(e)`@diffmat[ceiling(arena@n/2),ceiling(arena@m/2)] <- 100
arena@media$`EX_co2(e)`@difspeed=0.1
arena@media$`EX_co2(e)`@pde="BoundDiff2d"
#arena@media$`EX_co2(e)`@pde="Diff2d"
sim <- simEnv(arena, time=14, continue=T)
#sim@medlist[[2]]$`EX_co2(e)`[15]
out <- lapply(sim@medlist, function(x){matrix(x$`EX_co2(e)`, nrow=arena@n)})
par(mfrow = c(ceiling(sqrt(length(out))), ceiling(sqrt(length(out)))))
for(i in 1:length(out)){image(out[[i]], main=paste("time",i))}

par(mfrow=c(1,1))
conservation <- lapply(sim@medlist, function(x){sum(matrix(x$`EX_co2(e)`, nrow=arena@n))})
plot(seq(1, length(conservation)), unlist(conservation), type="b")

t(round(out[[2]],2))

par(mfrow=c(1,1))
plotCurves2(sim)



























x=c(10*10, 25*25, 51*51, 61*61, 71*71, 81*81, 91*91, 101*101)
y=c(3901, 29911, 160000, 230000, 330000, 430000, 580000, 710000)
lm <- lm(y~x)
summary(lm)
plot(x,y)
(z <- line(x,y))
abline(coef(lm))
abline(coef=c(0, lm$coefficients[2]))


co2_dat <- c(evalArena())
evalArena(sim, plot_items = "EX_co2(e)", time=2)
evalArena(sim, plot_items = c("Population", "EX_co2(e)", "EX_o2(e)", "EX_for(e)"), time=5)
image(outb, ask = FALSE, mfrow = c(3, 3), main = paste("time", times))



setwd("~/uni/bacarena")
Rcpp::sourceCpp("src/diff.cpp")


library(ReacTran)
n=50
m=50
x.grid  <- setup.grid.1D(x.up = 0, L = 10, N = n)
y.grid  <- setup.grid.1D(x.up = 0, L = 10, N = m)
grid2D <- setup.grid.2D(x.grid, y.grid)
C.x.up   <- rep(1, times = n)
C.x.down <- rep(0, times = n)
C.y.up   <- rep(1, times = m)
C.y.down <- rep(0, times = m)
D.grid    <- setup.prop.2D(value = 1, y.value = 1, grid = grid2D)
Diff2d <- function (t, y, parms)  {
    CONC  <- matrix(nrow = n, ncol = m, data = y)
    dCONC <- tran.2D(CONC, grid = grid2D, D.grid = D.grid)$dC
    return (list(dCONC))
}
ma <- matrix(data = rep(0, n*m), nrow = n)
ma[ceiling(n)/2,ceiling(m)/2] <- 100
image(ma)
print(ma)
#outb <- ode.2D (y = ma, func = Diff2d, t = 1:2, parms=NULL,
#                dim = c(n, m), lrw = 950000)
outb <- ode.2D (y = ma, func = Diff2d, t = 1:2, parms=NULL,
                dim = c(n, m), lrw=16000)

ma <- matrix(outb[2,][-1], ncol=m, nrow=n)
print(ma)
image(ma)
image(outb, ask = FALSE, mfrow = c(3, 3), main = paste("time", times))



m <- matrix(data = rep(0, 10*10), nrow = 10)
m[5,5] <- 100
image(m)
print (m)
diffuseSteveCpp(y = m, donut = T, D = 3, h = 1, tstep=0.1)
image(m)
print (m)





m <- matrix(data = rep(0, 100*100), nrow = 100)
m[50,50] <- 100
for(t in seq(100)){
  image(m)
  diffuseSteveCpp(y = m, donut = T, D = 1, h = 1, tstep=0.1)
}





matrix(data = runif(100*100), nrow = 100)

b <- microbenchmark(
  cpp <- diffuseNaiveCpp(y=matrix(data = runif(100*100), nrow = 100), donut = T),
  grajdeanu <- diffuseGrajdeanuCpp(y=matrix(data = runif(100*100), nrow = 100), mu = 1, donut = T),
  steve <- diffuseSteveCpp(y = matrix(data = runif(100*100), nrow = 100), D = 1, h = 1, tstep=0.1),
  pde <- ode.2D(y=matrix(data = runif(100*100), nrow = 100), func = Diff2d, t = 1:2, parms=NULL, dim = c(100, 100),  method="rk4"),
  times=100L)
b
boxplot(b)


m = matrix(data = runif(100*100), nrow = 100)
sum(m)
diffuseSteveCpp(m, donut=T, D = 1, h = 1, tstep = 0.1)
sum(m)
image(m)


