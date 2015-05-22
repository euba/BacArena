setwd("~/uni/bacarena")
Rcpp::sourceCpp("src/diff.cpp")

# Benchmark
library(microbenchmark)
library(ggplot2)
sub <- Substance(n=20,m=20,smax=40,name='test',difunc='r')
b <- microbenchmark(
  diffuseGrajdeanuCpp(matrix(data=runif(10000), ncol=100, nrow=100), mu=10, donut=F),
  diffuseNaiveCpp(matrix(data=runif(10000), ncol=100, nrow=100), donut=F),
  diffuseR(sub))
autoplot(b)

# random matrix
m = matrix(data=runif(25), ncol=5, nrow=5)
print(m)
sum(m)
diffuseGrajdeanuCpp(m, mu=10000, donut=T)
print(m)
sum(m)

# R function
library(Matrix)
sub <- Substance(n=20,m=20,smax=40,name='test',difunc='r')
diffuseR(sub)


m = matrix(data=runif(4), ncol=2, nrow=2)
print(m)
sum(m)
diffuseNaiveCpp(m, donut=F)
print(m)
sum(m)


# zero matrix with a center peak
m = matrix(data=c(0,0,0,0,1,0,0,0,0), nrow=3, ncol=3, byrow = TRUE)
diffuseGrajdeanuCpp(m, donut=T)
sum(m)

# vector
m = matrix(data=c(0,1,0), nrow=1, ncol=3, byrow = TRUE)
diffuseGrajdeanuCpp(m, donut=T)
