#Diffusion and movement

library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)

source(file="fba.R")
source(file="cpp_source.R")
sbml <- read.sbml("data/ecoli_core.xml")


#Variable Declaration

n <- 16
m <- 16
iter <- 300

bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)


movement <- cxxfunction(signature(A = "numeric"), body = src_movement, plugin="Rcpp")
diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")

movement2 <- function(bac, n, m){
  y <- bac
  for (i in 0:(n-1)){
    for(j in 0:(m-1)){
      if(y[i+1,j+1] != 0){
        a <- (i + round(runif(1,-1,1))) %% n 
        b <- (j + round(runif(1,-1,1))) %% m 
        if(bac[a+1,b+1] == 0){
          bac[a+1,b+1] <- 1
          bac[i+1,j+1] <- 0
        }
      }
    }
  }
  return(bac)
}


#enum <- as.vector(sapply(letters, function(x){
#  for(i in 1:9){
#    a[i] <- paste(x, i, sep="")
#  }
#  return(a)
#}))

#Initiation of agents

#
# intial Substrate distribution
#
s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  #matrix(sample(1:100, n*m, replace=T), nrow=n, ncol=m) # integer  in [0:100]
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # doubles in [0:100]
  
  #matrix(c(runif(n), rep(0, n*m-n)), nrow=n, ncol=m)
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m)
  matrix(c(rep(0,(n*m-2*n)/2), rep(100,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s

#
# intial bacterial distribution
#
#bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
#bac <- matrix(c(rep(1,2*n), rep(0,n*m-2*n)), nrow=n, ncol=m)
bac <- matrix(c(rep(0,(n*m-2*n)/2), rep(1,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)


#Iteration with rules to apply for each agent
for(time in 1:iter){      
  #diffusion of substrates by mean over neighbourhood
  #substrat <- lapply(substrat, function(x){
  #  anb <- eightneighbours(x)
  #  nb <- neighbours(x)
  #  mat <- (anb+x)/(nb+1)
  #  return(mat)
  #})
  diffusion(substrat)
  
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  image(substrat$M_glc_b, zlim=c(0,100), col=colorRampPalette(c("white", "black", "red"))(40))
  #par(mfrow=c(3,3))
  #for(i in 1:9){
  #  image(x[[i]], zlim=c(0,1), 
  #    col=colorRampPalette(c("white", "black", rainbow(14)[i]))(40),
  #    main=paste(names(x)[i], paste("step:", time)))
  #}
  #dev.off()
  print(paste("Sum of glucose:", sum(apply(substrat$M_glc_b, 1, sum))))
  
  #random movement of bacteria:
  #bac <- movement(bac)  
  #bac <- movement2(bac,n,m)
  bac <- t(apply(bac, 1, sample))  # movement by using random permutation matrix
  
  #fba
  bacnum <- sum(apply(bac, 1, sum))
  print(paste("Bacs:", bacnum))
  gvec <- 1:bacnum
  k <- 0
  for (i in 1:n){
    for(j in 1:m){
      if(bac[i,j] == 1){ # if there is a Bacterial
        spos <- lapply(substrat, function(x, i, j){
          return(x[i,j])
        },i=i, j=j)
        growth <- fba(spos, sbml$stoch, sbml$lb, sbml$ub, sbml$ex)
        gvec[k]=growth
        k <- k + 1
      }
    }
  }
  print(sd(gvec))
    
  
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  image(bac, col=c("white", "black"))
  #dev.off()
}
