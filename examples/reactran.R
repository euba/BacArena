library(ReacTran)

n     <- 10
m     <- 10
D     <- 1           # diffusion coeff, X- and Y-direction
ini   <- 1           # initial value at x=0

x.grid  <- setup.grid.1D(x.up = 0, L = 10, N = n)
y.grid  <- setup.grid.1D(x.up = 0, L = 10, N = m)
grid2D <- setup.grid.2D(x.grid, y.grid)




# Boundary conditions: fixed concentration  
C.x.up   <- rep(1, times = n)
C.x.down <- rep(0, times = n)
C.y.up   <- rep(1, times = m)
C.y.down <- rep(0, times = m)

Diff2d <- function (t, y, pars)  {
  with (as.list(pars), {
    print(test)
    CONC  <- matrix(nrow = n, ncol = m, data = y)
  #  dCONC <- tran.2D(CONC, grid = grid2D, D.grid = D.grid,
  #                   C.x.up = C.x.up, C.x.down = C.x.down,
  #                   C.y.up = C.y.up, C.y.down = C.y.down)$dC
    dCONC <- tran.2D(CONC, grid = grid2D, D.grid = D.grid)$dC
    return (list(dCONC))
  })
}

# initial condition: 0 everywhere, except in central point
y0 <- matrix(nrow = n, ncol = m, data = 0)
#y0[ceiling(n/2),ceiling(m/2)] <- ini  # initial concentration in the central point...
y0[1,1] <- ini

# solve for 8 time units
times <- 1:8
outb <- ode.2D (y = y0, func = Diff2d, t = times, parms=c(D=4, test="asd"),
                dim = c(n, m), lrw = 160000)

image(outb, ask = FALSE, mfrow = c(3, 3), main = paste("time", times))

matrix(data=outb[2,][-1], ncol=m, nrow=n)
sum(matrix(data=outb[2,][-1], ncol=m, nrow=n))
