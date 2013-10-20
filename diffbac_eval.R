#just for plotting
library(ggplot2)
library(reshape2)
library(reshape)
library(scales)
library(gridExtra)

setwd("~/BacArena")
source(file="baggage.R")

#load("BacArena_data.RData")
#load("results/barkeri_50x50_.RData")
#load("results/ecoli_20x20_aerob_seed55.RData")
#load("results/beijerinckii_20x20_seed943.RData")
#load("results/barkeri_20x20_seed9659.RData")
#load("results/barkeri_h2_20x20_seed6083.RData")
load("results/barkeri_beijerinckii_20x20_seed6764.RData")


#fav_subs <-c("acetate","co2","ethanol","formiate","glucose","o2") #for ecoli
#fav_subs <-c("acetate","co2","glucose","h2","butyrate", "succinate") #for beijerinckii
#fav_subs <-c("methanol","co2","methane","h2o") #for barkeri
#fav_subs <-c("h2o","co2","methane","h2") #for barkeri on h2/co2
fav_subs <-c("h2o","co2","methane","h2","acetate","glucose","butyrate", "succinate") #for barkeri,beijerinckii
#fav_subs <- names(BacArena_data[[1]]$substrat)
plot_list <- list()
growth_vec_history <- list()
substrat <- BacArena_data[[1]]$substrat
substrat_history <- matrix(data=0, nrow=length(substrat[fav_subs]), ncol=length(BacArena_data))
bac_history <- sapply(levels(bac[,3]), function(x){list()[[x]]}) # init list with entry for each bac type

for(time in 1:length(BacArena_data)){
  substrat <- BacArena_data[[time]]$substrat
  bac <- BacArena_data[[time]]$bac
  growth_vec <- bac$growth
  growth_vec_history[[time]] <- growth_vec
  if (dim(bac)[1] >= 1) {
    bac_color <- as.numeric(as.factor(levels(bac[,3])))
    names(bac_color) <- levels(bac[,3])
  }
  
  substrat_history[,time] <- unlist(lapply(substrat,FUN=mean))[fav_subs]
  rownames(substrat_history) <- fav_subs
  
  sapply(levels(bac[,3]), function(x,time,tab,bac_history){
    bac_history[[x]][time] <<- tab[x]
  },time=time,tab=table(bac$type),bac_history=bac_history)
  
  
  #plot.bacs(time=time, bac=bac, substrate=substrat, growth_vec_history=growth_vec_history, sub_his=substrat_history, bac_his=bac_history, subnam1="glucose", subnam2="co2", subnam3="o2", prodnam="acetate", bac_color=bac_color)
  plot_list[[time]] <- plot.bacs.cool(bac=bac, time=time, substrat=substrat, sub="glucose", sub2="h2", sub3="co2", prod="methane", bac_col=bac_color)
}

#save as 8.49 × 7.72 inch
plot_list[[60]]$bac
plot_list[[100]]$bac
plot_list[[160]]$bac

#save as: 6x10
plot_list[[time]]$subs
plot_list[[time]]$growth
###############plotting for ggplot
# # 
#   lapply(plot_list[-1], function(x){
#    grid.newpage() # Open a new page on grid device
#    pushViewport(viewport(layout = grid.layout(3, 2)))
#    print(x$sub, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
#    print(x$prod, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
#    print(x$subs, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
#    print(x$bacpos, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
#    print(x$growth, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
#   })
