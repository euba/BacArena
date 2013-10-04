library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)
#just for plotting
#library(ggplot2)
#library(reshape2)
#library(reshape)
#library(scales)
#library(gridExtra)

setwd("~/BacArena")
source(file="fba.R")
source(file="cpp_source.R")
source(file="baggage.R")
source(file="config.R")

if(!exists("movement", mode="function")){ #test if function already exists -> saves time for testing source code
movement <- cxxfunction(signature(input_matrix = "matrix", input_frame = "data.frame"), body = src_movement, plugin="Rcpp")
diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")
}


plot_list <- list()
bac_history <- vector(mode="numeric")
substrat_history <- matrix(data=0, nrow=length(substrat), ncol=iter)
max_glucose = max(substrat$glucose)
max_acetate = 100
growth_vec <- vector(mode="numeric") # all current growth rate for statistics
growth_vec_history <- list(mode="numeric") # all  growth rate for statistics

########################################################################################################
###################################### MAIN LOOP #######################################################
########################################################################################################

for(time in 1:iter){
  bacnum <- dim(bac)[1]
  bac_history[time] <- bacnum
  substrat_history[,time] <- unlist(lapply(substrat,FUN=mean))
  rownames(substrat_history) <- names(substrat)
  
  diffusion(substrat)
  
  gvec <- 1:bacnum
  #print(bacnum)
    
  xr <- round(runif(bacnum, min = -1, max = 1))
  yr <- round(runif(bacnum, min = -1, max = 1))
  
  
  ########################################################################################################
  ###################################### BACTERIA LOOP ###################################################
  ########################################################################################################
  
  for(l in 1:bacnum){
  #if(dim(bac)[1]==0) break # well a little bit ugly but our sins will be forgiven, every time you use a break god kills another kitten
    #print(l)
    #print(bac)
        
    #
    # get variables according to bac type
    #
    sbml <- get_sbml(bac[l,]$type) # get sbml according to bac type
    biomassf <- get_biomassf(bac[l,]$type)
    sub_ex <- get_sub_ex(bac[l,]$type)
        
    i <- bac[l,][1,1]
    j <- bac[l,][1,2]
    
    substrat <- lapply(substrat, function(x){round(x, digits=2)}) # ROUNDING!!!
    
    spos <- lapply(substrat, function(x, i, j){ # get current substrat vector
      return(x[i,j])
    },i=i, j=j)

    
    ########################################################################################################
    ########################################### FBA ########################################################
    ########################################################################################################
    
    growth <- fba(spos, sbml$stoch, sbml$ex, sbml$reac, bac[l,][1,4], sub_ex, bac[l,]$type)
    
    #
    # check prog solutions first!
    #
    if(growth == "DEAD"){
      print("----")
      cat("no fba solution found for: ")
      print(bac[l,])
      #print(t(spos))
      print("----")
      #bac <- bac[-l,]
    }
    #
    # if there is a feasable fba solution continue:
    #
    else{
      growth <- lapply(growth, function(x){round(x, digits=2)}) # ROUNDING!!!
      growth_vec[l] <- growth[[biomassf]] # save current growth rate to plot it    
      if(growth[[biomassf]] != 0){ # continue only if there is growth !
        bac[l,][1,4] <- bac[l,][1,4] + growth[[biomassf]]
        
        sapply(names(sapply(substrat, names)),function(x,i,j,substrat){
          if(x %in% names(sub_ex)) { # only update substrate which are metabolic relevant for current organism
            #print(substrat[[x]][i,j])
            #print(growth[[sub_ex[[x]]]])
            
            substrat[[x]][i,j] <<- substrat[[x]][i,j] + growth[[sub_ex[[x]]]] # "<<-" is necessary for extern variable modification
          }
        },i=i,j=j,substrat=substrat)
        
        gvec[l]=growth[[biomassf]]
        #print(growth)
      }
 
########################################################################################################
###################################### MOVEMENT & DOUBLING #############################################
########################################################################################################

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
      if(bac[l,]$growth>1){ # test for duplication
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
  }

  
########################################################################################################
##################################### TO BE OR NOT TO BE ###############################################
########################################################################################################
  
  bac$growth <- bac$growth-0.08 #the cost of living
  bac <- bac[!(bac$growth<0.1),] #death
  #
  if(dim(bac)[1]==0){
    print("ALL BACTERIA DIED")
    break
  }
  growth_vec_history[[time]] <- growth_vec
  plot_list[[time]] <- plot.bacs(time=time, bac=bac, growth_vec_history=growth_vec_history, subnam1="glucose", subnam2="co2", subnam3="h2", prodnam="methane")
}


###############plotting

#lapply(plot_list[-1], function(x){
#  grid.newpage() # Open a new page on grid device
#  pushViewport(viewport(layout = grid.layout(3, 2)))
#  print(x$sub, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
#  print(x$prod, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
#  print(x$subs, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
#  print(x$bacpos, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
#  print(x$growth, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
#})
