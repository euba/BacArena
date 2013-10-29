library(Rcpp)
library(inline)
library(rbenchmark)
library(digest) # hashes

setwd("~/BacArena")
source(file="fba.R")
source(file="cpp_source.R")
source(file="baggage.R")
source(file="config.R")

set.seed(seed)

if(!exists("movement", mode="function")){ #test if function already exists -> saves time for testing source code
movement <- cxxfunction(signature(input_matrix = "matrix", input_frame = "data.frame", seed = "integer"), body = src_movement, plugin="Rcpp")
diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")
}


# variables for diagnostic plot
growth_vec_history <- list()
substrat_history <- matrix(data=0, nrow=length(substrat), ncol=iter)
bac_history <- sapply(levels(bac[,3]), function(x){list()[[x]]}) # init list with entry for each bac type

time_history <- list(mode="numeric") # time comsumption
fba_hash <-new.env() # hash fba results to improve performance

BacArena_data<- list()

########################################################################################################
###################################### MAIN LOOP #######################################################
########################################################################################################

for(time in 1:iter){
  #speed test
  time_fba    <- 0
  time_diff   <- 0
  time_mov    <- 0
  time_total  <- 0
  time_plot   <- 0
  time_tmp    <- proc.time()
  time_unk    <- 0
  
  #iterration status
  no_fba_found <- 0
  hash_uses <- 0
  
  growth_vec <- vector(mode="numeric") # for diagnostic plot
  substrat_history[,time] <- unlist(lapply(substrat,FUN=mean))
  rownames(substrat_history) <- names(substrat)

  bacnum <- dim(bac)[1]
  
  #time_tmp4 <- proc.time()
  #time_unk <- time_unk + proc.time() - time_tmp4
  
  time_tmp3 <- proc.time()
    min_tmp <- min(unlist(lapply(substrat, min)))
    diffusion(substrat)
    min_diff <- min(unlist(lapply(substrat, min)))
    if(min_diff < 0) {
      print(min_tmp)
      print(min_diff)
      stop("diffusion -> negative substrate")
    }
  time_diff <- proc.time() - time_tmp3
  
  xr <- round(runif(bacnum, min = -1, max = 1))
  yr <- round(runif(bacnum, min = -1, max = 1))
  
  ########################################################################################################
  ###################################### BACTERIA LOOP ###################################################
  ########################################################################################################
  
  for(l in 1:bacnum){     
    time_tmp4 <- proc.time()
    #
    # get variables according to bac type
    #
    sbml <- get_sbml(bac[l,]$type) # get sbml according to bac type
    biomassf <- get_biomassf(bac[l,]$type)
    sub_ex <- get_sub_ex(bac[l,]$type)
        
    i <- bac[l,][1,1]
    j <- bac[l,][1,2]
        
    spos <- lapply(substrat, function(x, i, j){ # get current substrat vector
      return(x[i,j])
    },i=i, j=j)
    time_unk <- time_unk + proc.time() - time_tmp4
    
    ########################################################################################################
    ########################################### FBA ########################################################
    ########################################################################################################
    
    
    time_tmp3 <- proc.time()

    #
    # Hashing to reuse old fba results
    #
    
      hash_spos <- digest(floor(unlist(spos))) # rounding!!
      if(exists(hash_spos, envir=fba_hash)) {
        growth <- fba_hash[[hash_spos]]
        hash_uses <- hash_uses + 1
      #  stop("ATTENTION: TAKING HASH!!!!")
      }
      else {
        growth <- fba(spos, sbml$stoch, sbml$ex, sbml$reac, bac[l,][1,4], sub_ex, bac[l,]$type)
        assign(hash_spos, growth, envir=fba_hash)
      }
    #growth <- fba(spos, sbml$stoch, sbml$ex, sbml$reac, bac[l,][1,4], sub_ex, bac[l,]$type)
    
    time_fba <- time_fba + proc.time() - time_tmp3
    time_tmp4 <- proc.time()
    #
    # check prog solutions first!
    #
    if(length(growth)==1){ # <=> if(growth == "DEAD"){
      no_fba_found <- no_fba_found+1

       # new growth rate according to advanced magic calculation *muah*
       starving_growth <- -(get_ngam(bac[l,]$type)/get_gam(bac[l,]$type)) 
       bac[l,]$growth <- starving_growth + bac[l,]$growth
       growth_vec[l] <- starving_growth  # for diagnostic plot
    }
    #
    # if there is a feasable fba solution continue:
    #
    else{
      growth_vec[l] <- growth[[biomassf]] # for diagnostic plot
      if(growth[[biomassf]] != 0){ # continue only if there is growth !
        bac[l,][1,4] <- bac[l,][1,4] + growth[[biomassf]]
        
        sapply(names(sapply(substrat, names)),function(x,i,j,substrat){
          if(x %in% names(sub_ex)) { # only update substrate which are metabolic relevant for current organism
            #print(substrat[[x]][i,j])
            #print(growth[[sub_ex[[x]]]])
            
            #
            # rounding to avoid gray bug (different variable accuracy leads to artefacts -> very small negative differences from zero which produces a gray plot)
            #
            #substrat[[x]][i,j] <<- substrat[[x]][i,j] + growth[[sub_ex[[x]]]] # "<<-" is necessary for extern variable modification
            substrat[[x]][i,j] <<- round(substrat[[x]][i,j] + growth[[sub_ex[[x]]]], digits=-log10(epsilon)) # "<<-" is necessary for extern variable modification
            
            if(substrat[[x]][i,j] < 0){
              print (growth[[sub_ex[[x]]]])
              print (substrat[[x]][i,j])
              stop("negative stubstrate!")
            } 
          }
        },i=i,j=j,substrat=substrat)
      }
    }
    time_unk <- time_unk + proc.time() - time_tmp4
########################################################################################################
###################################### MOVEMENT & DOUBLING #############################################
########################################################################################################
    
    time_tmp2 <- proc.time()
    dupli <- F # boolean variable to test for duplication
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
        bac <- bac[-(dim(bac)[1]),]
      }
    }
    time_mov <- time_mov + proc.time() - time_tmp2
  }
  
  
########################################################################################################
##################################### TO BE OR NOT TO BE ###############################################
########################################################################################################
  
  time_tmp4 <- proc.time()
  #
  bac <- bac[!(bac$growth<0),] #death
  #
  if(dim(bac)[1]==0){
    print("ALL BACTERIA DIED")
    break
  }
  time_unk <- time_unk + proc.time() - time_tmp4
  
  growth_vec_history[[time]] <- growth_vec # for diagnostic plot
  
  time_tot <- proc.time() - time_tmp
  time_cur <- cbind(time_tot[3], time_diff[3], time_mov[3], time_fba[3], time_unk[3])
  colnames(time_cur) <- c("total", "diffusion", "movement", "fba", "unknown")
  #print(time_cur)
  time_history[[time]] <- time_cur
  
  ########################################################################################################
  ##################################### DATA GENERATION ##################################################
  ########################################################################################################
  
  BacArena_data[[time]] <- list(substrat=substrat, bac=bac)
  
  #
  # diagnostic plots 
  #
  #diag_plot(BacArena_data, time)
    
  #print interation status
  baccounter <- table(bac$type)
  iter_print <- c(time, no_fba_found, time_tot[3], hash_uses, seed, unname(baccounter))
  names(iter_print) <- c("iteration", "no_fba/growth", "time_elapsed", "#hashing", "seed", names(baccounter))
  print(iter_print)
}

#plot time consumption #save as 6 × 6 inch
m <- do.call(rbind, time_history)
plot(1:dim(m)[1], m[,1], type="l", col=1, pch=1, , ylab="computation time", xlab="time") #set max y-value to highest product conentration
for(i in 2:(dim(m)[2])){
  lines(1:dim(m)[1], m[,i], col=i, pch=i, type="l")
}
legend("topleft", colnames(time_cur), pch=1, col=c(1:dim(m)[2]), cex=0.64, bty="n")
plot(1:dim(m)[1], m[,1], ylim=c(0,1), type="n", col=1, pch=1, , ylab="rel. computation time", xlab="time") #set max y-value to highest product conentration
for(i in 2:(dim(m)[2])){
  lines(1:dim(m)[1], m[,i]/m[,1], col=i, pch=i, type="l")
}
legend("left", colnames(time_cur), pch=1, col=c(1:dim(m)[2]), cex=0.64, bty="n")

#save data
save(BacArena_data, file = "BacArena_data.RData")
