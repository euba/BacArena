setwd("~/uni/bacarena")
Rcpp::sourceCpp("src/diff.cpp")

# random matrix
m = matrix(data=runif(25), ncol=5, nrow=5)
print(m)
sum(m)
diffuseGrajdeanuCpp(m, donut=T)
print(m)
sum(m)

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
