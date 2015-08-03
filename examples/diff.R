library(microbenchmark)

setwd("~/uni/bacarena")
Rcpp::sourceCpp("src/diff.cpp")

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
  new <- diffuseSteveCpp(y = matrix(data = runif(100*100), nrow = 100), donut = T, D = 1, h = 1, tstep=0.1),
  times=1000L)
b
boxplot(b)


m = matrix(data = runif(100*100), nrow = 100)
sum(m)
diffuseSteveCpp(m, donut=T, D = 1, h = 1, tstep = 0.1)
sum(m)
image(m)


