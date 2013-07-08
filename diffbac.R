#Diffusion and movement

library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)

setwd("~/BacArena")
source(file="fba.R")
source(file="cpp_source.R")
#sbml <- read.sbml("data/ecoli_core.xml")
sbml <- read.sbml("/home/eugen/BacArena/data/ecoli_core.xml")

movement <- cxxfunction(signature(input_matrix = "matrix", input_frame = "data.frame"), body = src_movement, plugin="Rcpp")
diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")

#Variable Declaration

n <- 5
m <- 5
iter <- 5
bacs <- 5

#bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
bac <- data.frame(x=round(runif(bacs, min=1, max=m)), y=round(runif(bacs, min=1, max=n)), 
                  type=rep("ecoli", bacs), growth=rep(1, bacs))
bac <- bac[!duplicated(bac[,1:2]),]
rownames(bac) <- 1:nrow(bac) #change indices in data.frame

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
  matrix(c(rep(0,(n*m-2*n)/2), rep(10,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)
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
#bac <- matrix(c(rep(0,(n*m-2*n)/2), rep(1,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)

#Initializatio of reporter variables
mgvec <- vector("numeric")
sgvec <- vector("numeric")


bacnum <- dim(bac)[1]
gvec <- 1:bacnum

#Varma and Palsson 1994:
#the non-growth-associated maintenance requirements (7.6 mmol of ATP per g [dry weight] per h), and the growth-associated maintenance requirements (13 mmol of ATP per g of biomass).

#Iteration with rules to apply for each agent
for(time in 1:iter){        
  #
  #plotting functions:
  #
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=1000, height=1000)
  par(mfrow=c(1,2))
  image(substrat$M_glc_b, zlim=c(0,50), col=colorRampPalette(c("white", "black", "red"))(40))
  #par(mfrow=c(3,3))
  #for(i in 1:9){
  #  image(x[[i]], zlim=c(0,1), 
  #    col=colorRampPalette(c("white", "black", rainbow(14)[i]))(40),
  #    main=paste(names(x)[i], paste("step:", time)))
  #}
  #dev.off()
  bacnum <- dim(bac)[1]
  #print(paste("Bacs:", bacnum))
  #gvec <- 1:bacnum
  
  #mgvec[time] <- mean(gvec)
  #sgvec[time] <- sd(gvec)
  #plot(1:time, mgvec, type="b", xlab="Iteration", ylab="mean replication rate")
  #plot(1:time, sgvec, type="b", xlab="Iteration", ylab="standard deviation replication rate")
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  #dev.off()
  
  #Reporter:
  #print(paste("Sum of glucose:", sum(apply(substrat$M_glc_b, 1, sum))))
  #print(paste("Sum of O2:", sum(apply(substrat$M_o2_b, 1, sum))))
  #single element in matrix:
  #print(paste("Glucose:", substrat$M_glc_b[n/2,m/2]))
  #print(paste("O2:", substrat$M_o2_b[n/2,m/2]))
  #print(paste("Acetate:", substrat$M_ac_b[n/2,m/2]))
  #print(paste("CO2:", substrat$M_co2_b[n/2,m/2]))
  #print(paste("Ethanol:", substrat$M_etoh_b[n/2,m/2]))
  #print(paste("Succinate:", substrat$M_succ_b[n/2,m/2]))
  #print(paste("Akg:", substrat$M_akg_b[n/2,m/2]))
  #print(paste("Fumarate:", substrat$M_fum_b[n/2,m/2]))
  #print(paste("Formate:", substrat$M_for_b[n/2,m/2]))
  #print("")
  #diffusion of substrates by mean over neighbourhood
  #substrat <- lapply(substrat, function(x){
  #  anb <- eightneighbours(x)
  #  nb <- neighbours(x)
  #  mat <- (anb+x)/(nb+1)
  #  return(mat)
  #})
  #
  #diffusion(substrat)
  
  #print(paste("Glucose:", substrat$M_glc_b[n/2,m/2]))
  #print(paste("O2:", substrat$M_o2_b[n/2,m/2]))
  #print(paste("Acetate:", substrat$M_ac_b[n/2,m/2]))
  #print(paste("CO2:", substrat$M_co2_b[n/2,m/2]))
  #print(paste("Ethanol:", substrat$M_etoh_b[n/2,m/2]))
  #print(paste("Succinate:", substrat$M_succ_b[n/2,m/2]))
  #print(paste("Akg:", substrat$M_akg_b[n/2,m/2]))
  #print(paste("Fumarate:", substrat$M_fum_b[n/2,m/2]))
  #print(paste("Formate:", substrat$M_for_b[n/2,m/2]))
  #print("")
  #random movement of bacteria:
  
  #bac <- movement(bac)  
  #print(bac)
  bac_img <- movement(matrix(0,n,m), bac) # move bacs and return matrix for printing
  image(bac_img, col=c("white", "black"))
  #print("")
  #print(bac)
  
  #bac <- movement2(bac,n,m)
  #bac <- t(apply(bac, 1, sample))  # movement by using random permutation matrix
  
  #fba
  #bacnum <- sum(apply(bac, 1, sum))
  bacnum <- dim(bac)[1]
  #print(paste("Bacs:", bacnum))
  gvec <- 1:bacnum
  
  xr <- round(runif(bacnum, min = -1, max = 1))
  yr <- round(runif(bacnum, min = -1, max = 1))
  for(l in 1:bacnum){
    i <- bac[l,][1,1]
    j <- bac[l,][1,2]
    spos <- lapply(substrat, function(x, i, j){ # get current substrat vector
      return(x[i,j])
    },i=i, j=j)
    growth <- fba(spos, sbml$stoch, sbml$lb, sbml$ub, sbml$ex, sbml$reac)
    bac[l,][1,4] <- bac[l,][1,4] + growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]
    #
    # new substrat vector after metabolism
    #print("")
    #print(paste("glucose before: ", substrat[["M_glc_b"]][[i,j]]))
    sapply(names(sapply(substrat, names)),function(x,i,j,substrat){
      #growth[[sub_ex[[x]]]]
      substrat[[x]][i,j] <<- substrat[[x]][i,j] + growth[[sub_ex[[x]]]] # "<<-" is necessary for extern variable modification
    },i=i,j=j,substrat=substrat)
    #print(paste("glucose uptake by bac: (", i, ",", j, ")=", growth[["R_EX_glc_e_"]]))
    #print(paste("growth: ",growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]))
    #print(paste("glucose after: ", substrat[["M_glc_b"]][[i,j]]))
    #
    #print(growth)
    gvec[l]=growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]
    
    #live and die
    if(bac[l,]$growth>1)  bac[-l,]
    #movement
    a <- (i + xr[l])
    b <- (j + yr[l])
    if(a == -1){a = n-1}
    if(b == -1){b = m-1}
    if(a == n){a = 0}
    if(b == m){b = 0}
    #if(!(bac[,1:2]==c(a,b))){ # if empty go for it!
    bac[l,1:2] <- c(a,b)
    #}
  }
  
  #print(paste("Glucose:", substrat$M_glc_b[n/2,m/2]))
  #print(paste("O2:", substrat$M_o2_b[n/2,m/2]))
  #print(paste("Acetate:", substrat$M_ac_b[n/2,m/2]))
  #print(paste("CO2:", substrat$M_co2_b[n/2,m/2]))
  #print(paste("Ethanol:", substrat$M_etoh_b[n/2,m/2]))
  #print(paste("Succinate:", substrat$M_succ_b[n/2,m/2]))
  #print(paste("Akg:", substrat$M_akg_b[n/2,m/2]))
  #print(paste("Fumarate:", substrat$M_fum_b[n/2,m/2]))
  #print(paste("Formate:", substrat$M_for_b[n/2,m/2]))
  #print(paste("Bac:", bac[n/2,m/2]))
  #print("")
  #
  # fba test (test123 is artificial substrate vector)
  #
  #("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
  #test123 <- c(0,0,0,0,1000,1000,1000,1000,1000,0,2.5,1000,0,0)
  #test123 <- c(1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000)
  #names(test123) <- s
  #print(fba(test123, sbml$stoch, sbml$lb, sbml$ub, sbml$ex, sbml$reac))
  #
}
