#Diffusion and movement

library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)

source(file="fba.R")
source(file="cpp_source.R")
sbml <- read.sbml("data/ecoli_core.xml")


#Variable Declaration

n <- 4
m <- 4
iter <- 1

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
        if(bac[a+1,b+1] == 1){
          bac[a+1,b+1] <- 1
          bac[i+1,j+1] <- 0
        }
      }
    }
  }
  return(bac)
}

a <- vector(mode="character")
enum <- as.vector(sapply(letters, function(x){
  for(i in 1:9){
    a[i] <- paste(x, i, sep="")
  }
  return(a)
}))

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
  matrix(c(rep(0,(n*m-2*n)/2), rep(1000,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s
# associate each substrate with an exchange reaction (sbml specific!!)
sub_ex <- character(length(s))
names(sub_ex) <- s
sub_ex[["M_ac_b"]]    <- "R_EX_ac_e_"
sub_ex[["M_akg_b"]]   <- "R_EX_akg_e_"
sub_ex[["M_co2_b"]]   <- "R_EX_co2_e_"
sub_ex[["M_etoh_b"]]  <- "R_EX_etoh_e_"
sub_ex[["M_fum_b"]]   <- "R_EX_fum_e_" 
sub_ex[["M_for_b"]]   <- "R_EX_for_e_"
sub_ex[["M_glc_b"]]   <- "R_EX_glc_e_"
sub_ex[["M_h2o_b"]]   <- "R_EX_h2o_e_"
sub_ex[["M_h_b"]]     <- "R_EX_h_e_"
sub_ex[["M_lac_D_b"]] <- "R_EX_lac_D_e_"
sub_ex[["M_o2_b"]]    <- "R_EX_o2_e_"
sub_ex[["M_pi_b"]]    <- "R_EX_pi_e_"
sub_ex[["M_pyr_b"]]   <- "R_EX_pyr_e_"
sub_ex[["M_succ_b"]]  <- "R_EX_succ_e_"
                    

#
# intial bacterial distribution
#
#bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
#bac <- matrix(c(rep(1,2*n), rep(0,n*m-2*n)), nrow=n, ncol=m)
bac <- matrix(c(rep(0,(n*m-2*n)/2), rep(1,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)


#Iteration with rules to apply for each agent
mgvec <- vector("numeric")
sgvec <- vector("numeric")
for(time in 1:iter){      
  #diffusion of substrates by mean over neighbourhood
  #substrat <- lapply(substrat, function(x){
  #  anb <- eightneighbours(x)
  #  nb <- neighbours(x)
  #  mat <- (anb+x)/(nb+1)
  #  return(mat)
  #})
  #
  #diffusion(substrat)
  
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=1000, height=1000)
  par(mfrow=c(2,2))
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
        spos <- lapply(substrat, function(x, i, j){ # get current substrat vector
          return(x[i,j])
        },i=i, j=j)
        growth <- fba(spos, sbml$stoch, sbml$lb, sbml$ub, sbml$ex, sbml$reac)
        #
        # new substrat vector after metabolism
        print("")
        print(paste("glucose before: ", substrat[["M_glc_b"]][[i,j]]))
        sapply(names(sapply(substrat, names)),function(x,i,j,substrat){
          growth[[sub_ex[[x]]]]
          #substrat[[x]][i,j] <<- substrat[[x]][i,j] + growth[[sub_ex[[x]]]] # "<<-" is necessary for extern variable modification
        },i=i,j=j,substrat=substrat)
        print(paste("glucose uptake by bac: (", i, ",", j, ")=", growth[["R_EX_glc_e_"]]))
        print(paste("growth: ",growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]))
        print(paste("glucose after: ", substrat[["M_glc_b"]][[i,j]]))
        #
        #print(growth)
        gvec[k]=growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]
        k <- k + 1
      }
    }
  }
  #
  # fba test (test123 is artificial substrate vector)
  #
  #test123 <- c(0,0,0,0,0,0,1000,1000,1000,0,1000,1000,0,0)
  #names(test123) <- s
  #fba(test123, sbml$stoch, sbml$lb, sbml$ub, sbml$ex, sbml$reac)
  #
  
  image(bac, col=c("white", "black"))
  mgvec[time] <- mean(gvec)
  sgvec[time] <- sd(gvec)
  #plot(1:time, mgvec, type="b", xlab="Iteration", ylab="mean replication rate")
  #plot(1:time, sgvec, type="b", xlab="Iteration", ylab="standard deviation replication rate")
  
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  #dev.off()
}
