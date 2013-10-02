#Diffusion and movement

library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)

setwd("~/BacArena")
source(file="fba.R")
source(file="cpp_source.R")
source(file="baggage.R")

#movement <- cxxfunction(signature(input_matrix = "matrix", input_frame = "data.frame"), body = src_movement, plugin="Rcpp")
diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")

#
# Variable Declaration
#

n <- 5
m <- 5
iter <- 2
bacs <- 1
smax <- 70
#ecoli
#s <- c("acetate","aketoglutarate", "co2", "ethanol", "formiate", "fumarate", "glucose", "h2o", "proton", "lactate","o2", "iphosphate", "pyruvate", "succinate")
#mbarkeri
s <- c("acetate","aketoglutarate", "co2", "ethanol", "formiate", "fumarate", "glucose", "h2o", "proton", "lactate","o2", "iphosphate", "pyruvate", "succinate", "h2", "methanol", "methane")

#
# loading bacteria
#
source(file="ecoli.R")
source(file="barkeri.R")


#
# Initiation of agents
#

#bac <- data.frame(x=round(runif(bacs, min=1, max=n)), y=round(runif(bacs, min=1, max=m)),
#                  type=rep("ecoli", bacs), growth=rep(1, bacs))
#bac <- bac[!duplicated(bac[,1:2]),]
#rownames(bac) <- 1:nrow(bac) #change indices in data.frame


#bac <- data.frame(x=round(n/2), y=round(m/2),type="ecoli", growth=1) # one cell in the centre
bac <- data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1) # one cell in the centre
#bac <- data.frame(x=n, y=m,type="ecoli", growth=1) # one cell in the centre

#
# intial Substrate distribution
#
substrat <- lapply(s, function(x, n, m){
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # random substrate
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m) # downstairs substrate
  #matrix(c(rep(0,(n*m-2*n)/2), rep(10,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m) # substrate in the middle of our street ohooo
  matrix(smax,n,m) # homogen substrate distribution
}, n=n, m=m)
names(substrat) <- s
#substrat[["o2"]] <- matrix(0,n,m)
#substrat[["acetate"]] <- matrix(0,n,m)
#substrat[["formiate"]] <- matrix(0,n,m)
#substrat[["ethanol"]] <- matrix(0,n,m)
#substrat[["fumarate"]] <- matrix(0,n,m)
#substrat[["h2o"]] <- matrix(99999,n,m)



#
#Iteration with rules to apply for each agent
#
bac_history <- vector(mode="numeric")
substrat_history <- matrix(data=0, nrow=length(substrat), ncol=iter)
max_glucose = max(substrat$glucose)
max_acetate = 100
for(time in 1:iter){
  #
  #plotting functions
  par(mfrow=c(3,2))
  image(substrat$glucose, zlim=c(0,max_glucose), col=colorRampPalette(c("white", "green"))(40), main="glucose concentration")
  #image(substrat$glucose, zlim=c(0,70), col=colorRampPalette(c("white", "green"))(40), main="glucose concentration")

  #image(substrat$acetate, zlim=c(0,200), col=colorRampPalette(c("white", "orange"))(40), main="acetate concentration")
  image(substrat$acetate, zlim=c(0,max_acetate), col=colorRampPalette(c("white", "orange"))(40), main="acetate concentration")
  
  substrat_history[,time] <- unlist(lapply(substrat,FUN=mean))
  rownames(substrat_history) <- names(substrat)
  plot(1:time, substrat_history["h2o",1:time], col="blue", ylim=c(0,max(substrat_history[,1:time])), ylab="concentration", xlab="time") #set max y-value to highest product conentration
  lines(1:time, substrat_history["glucose",1:time], col="green", type="b")
  lines(1:time, substrat_history["o2",1:time], col="cyan", type="b")
  lines(1:time, substrat_history["co2",1:time], col="magenta", type="b")
  legend("left", c("water", "glucose", "o2", "co2"), pch=1,col=c("blue", "green", "cyan", "magenta"))
  
  plot(1:time, substrat_history["acetate",1:time], col="orange", pch=1, ylim=c(0,max_acetate), ylab="concentration", xlab="time")
  lines(1:time, substrat_history["fumarate",1:time], col="gray", pch=2, type="b")
  lines(1:time, substrat_history["formiate",1:time], col="red", pch=3, type="b")
  lines(1:time, substrat_history["ethanol",1:time], col="brown", pch=4, type="b")
  legend("left", c("acetate", "fumarate", "formiate", "ethanol"), pch=c(1,2,3,4),col=c("orange", "gray", "red", "brown"))
  
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
    #growth <- fba(spos, sbml$stoch, lb, ub, sbml$ex, sbml$reac, bac[l,][1,4], sub_ex, bac[l,]$type)
    #print(spos)
    print(bac)
    growth <- fba(spos, sbml$stoch, sbml$ex, sbml$reac, bac[l,][1,4], sub_ex, bac[l,]$type)
    
    lb <- get_lower_bound(bac[l,]$type)
    ub <- get_upper_bound(bac[l,]$type)
    
    # check for feasable lin prog solutions first!
    if(growth == "DEAD"){
      print("----")
      cat("no fba solution found for: ")
      print(bac[l,])
      print(t(spos))
      print("----")
      bac[-l,]
    }else{
      growth <- lapply(growth, function(x){round(x, digits=2)}) # ROUNDING!!!
      if(growth[[biomassf]] != 0){ # continue only if there is growth !
        bac[l,][1,4] <- bac[l,][1,4] + growth[[biomassf]]
        
        sapply(names(sapply(substrat, names)),function(x,i,j,substrat){
          if(x %in% names(sub_ex)) { # only update substrate which are metabolic relevant for current organism
            #print(substrat[[x]][i,j])
            #print(growth[[sub_ex[[x]]]])
            
            substrat[[x]][i,j] <<- substrat[[x]][i,j] + growth[[sub_ex[[x]]]] # "<<-" is necessary for extern variable modification
          }
        },i=i,j=j,substrat=substrat)
        
        # check for errorlike substrate uptake/concentration
        if(max(sapply(substrat,max))>1e+20){
          print(t(growth))
          print("")
          print("lower bound")
          # get lower bound with names
          lbound <- sapply(names(sapply(substrat, names)), function(x, stoch, sub_ex, lb){
            if(x %in% names(sub_ex)) lb[which(colnames(stoch)==sub_ex[[x]])]
          },stoch=sbml$stoch, sub_ex=sub_ex, lb=lb)
          print(t(lbound))
          print("")
          print("upper bound")
          rbound <- sapply(names(sapply(substrat, names)), function(x, stoch, sub_ex, ub){
            if(x %in% names(sub_ex)) ub[which(colnames(stoch)==sub_ex[[x]])]
          },stoch=sbml$stoch, sub_ex=sub_ex, ub=ub)
          print(t(rbound))
          print("")
          print(t(spos))
          print(bac[l,])
          print(c(x, substrat[[x]][i,j], " uptake: ", -growth[[sub_ex[[x]]]]))
          stop("FBA ERROR: Unusual high substrate concentrations!! (~infinity)")
        }
        
        gvec[l]=growth[[biomassf]]
        #print(growth)
      }
      
      #
      # Movement in R
      #
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
  
  #
  # Live and die
  #
  bac$growth <- bac$growth-0.1 #the cost of living
  bac <- bac[!(bac$growth<0.1),] #death
  
  if(dim(bac)[1]==0){
    print("ALL BACTERIA DIED")
    break
  }
}