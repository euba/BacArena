library(microbenchmark)

setwd("~/uni/bacarena")
Rcpp::sourceCpp("src/diff.cpp")


library(ReacTran)
n=100
m=100
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
                dim = c(n, m), method="rk4")

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


