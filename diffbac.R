#Diffusion and movement

library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)

setwd("~/BacArena")
source(file="fba.R")
source(file="cpp_source.R")
sbml <- read.sbml("data/ecoli_core.xml")

movement <- cxxfunction(signature(input_matrix = "matrix", input_frame = "data.frame"), body = src_movement, plugin="Rcpp")
diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")

#
# Variable Declaration
#

n <- 20
m <- 20
iter <- 50
bacs <- 8

#
# Initiation of agents
#

bac <- data.frame(x=round(runif(bacs, min=1, max=n)), y=round(runif(bacs, min=1, max=m)), 
                  type=rep("ecoli", bacs), growth=rep(1, bacs))
bac <- bac[!duplicated(bac[,1:2]),]
rownames(bac) <- 1:nrow(bac) #change indices in data.frame


#
# intial Substrate distribution
#
s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # random substrate
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m) # downstairs substrate
  #matrix(c(rep(0,(n*m-2*n)/2), rep(10,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m) # substrate in the middle of our street ohooo
  matrix(20,n,m) # homogen substrate distribution
}, n=n, m=m)
names(substrat) <- s
#
# associate each substrate with an exchange reaction (sbml specific!!)
#
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
#Iteration with rules to apply for each agent
#
bac_history <- vector(mode="numeric")
for(time in 1:iter){        
  #
  #plotting functions
   par(mfrow=c(2,2))
   image(substrat$M_glc_b, zlim=c(0,50), col=colorRampPalette(c("white", "black", "red"))(40), main="glucose concentration")
   image(substrat$M_o2_b, zlim=c(0,50), col=colorRampPalette(c("white", "black", "red"))(40), main="oxygen concentration")
   bacnum <- dim(bac)[1]
   
   bac_history[time] <- bacnum
   plot(1:time, bac_history, type="b", main="growth curve")
  
  #
  # Model of Diffusion
  #
  diffusion(substrat)
  
  #
  # Model of Movement
  #
  #print(bac)
  #tmp <- movement(matrix(0,n,m), bac) # move bacs and return matrix for printing
  #bac_img <- tmp$matrix
  #bac <- tmp$df
  #image(bac_img, col=c("white", "black"))
  
  
  mat = matrix(0,n,m) # conversion of data frame into bac matrix
  apply(bac[,1:2], 1, function(x){
    mat[as.numeric(x[1]), as.numeric(x[2])] <<- 1
  })
   image(mat, col=c("white", "black"), main="bacterial movement")
  
  
  #
  # FBA
  #
  gvec <- 1:bacnum
  #print(bacnum)
  
  xr <- round(runif(bacnum, min = -1, max = 1))
  yr <- round(runif(bacnum, min = -1, max = 1))
  for(l in 1:bacnum){
    i <- bac[l,][1,1]
    j <- bac[l,][1,2]
    spos <- lapply(substrat, function(x, i, j){ # get current substrat vector
      return(x[i,j])
    },i=i, j=j)
    growth <- fba(spos, sbml$stoch, sbml$lb, sbml$ub, sbml$ex, sbml$reac, bac[l,][1,4], sub_ex)
    bac[l,][1,4] <- bac[l,][1,4] + growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]

    sapply(names(sapply(substrat, names)),function(x,i,j,substrat){
      substrat[[x]][i,j] <<- substrat[[x]][i,j] + growth[[sub_ex[[x]]]] # "<<-" is necessary for extern variable modification
    },i=i,j=j,substrat=substrat)
    
    gvec[l]=growth[["R_Biomass_Ecoli_core_N__w_GAM_"]]
    
    #
    # Movement in R 
    #
    a <- (i + xr[l])
    b <- (j + yr[l])
    if(a == 0){a = n}
    if(b == 0){b = m}
    if(a == n+1){a = 1}
    if(b == m+1){b = 1}
    test <- apply(bac[,1:2], 1, function(x, p){
      if(sum(x==p)==2){
        return(T)
      }else{
        return(F)
      }
    }, p=c(a,b))
    if(bac[l,]$growth>2){ # test for duplication
      bac[l,]$growth <- bac[l,]$growth/2
      bac <- rbind(bac, bac[l,])
      dupli <- T
    }
    if(!(sum(test)>=1)){ # if empty go for it!
      bac[l,1:2] <- c(a,b)
    }else{
      if(dupli){ # if neighbour not empty and cell duplicated, kill doughter cell
        bac <- bac[-(bacnum+1),]
      }
    }
  }
  
  #
  # Live and die
  #
  bac$growth <- bac$growth-0.2 #the cost of living
  bac <- bac[!(bac$growth<0.5),] #death
  if(dim(bac)[1]==0){
    print("ALL BACTERIA DIED")
    break
  }
}
