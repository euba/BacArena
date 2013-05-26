library("simecol")

n <- 100
m <- 100

test <- new("gridModel",
            main = function(time, init, parms) {
              x <- init
              y <- x
              for (i in 0:(parms$n-1)){
                for(j in 0:(parms$m-1)){
                  if(y[i+1,j+1] != 0){
                    a <- (i + round(runif(1,-1,1))) %% parms$n 
                    b <- (j + round(runif(1,-1,1))) %% parms$m 
                    if(x[a+1,b+1] == 0){
                      x[a+1,b+1] <- 1
                      x[i+1,j+1] <- 0
                      #print(c(i,j,a,b))
                    }
                  }
                }
              }
              print(sum(apply(x, 1, sum)))
              #print(x)
              return(x)
            },
            parms = list(n=n, m=m),
            times = c(from=1, to=100, by=1),
            init = matrix(round(runif(n*m)), nrow=n, ncol=m),
            #init = matrix(c(rep(0,24),1), nrow=n, ncol=m),
            solver = "iteration"
)

## and to run this example:

plot(sim(test))
