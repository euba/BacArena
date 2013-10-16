#just for plotting
library(ggplot2)
library(reshape2)
library(reshape)
library(scales)
library(gridExtra)

source(file="baggage.R")

setwd("~/BacArena")
load("BacArena_data.RData")

plot_list <- list()
growth_vec_history <- list()
substrat_history <- matrix(data=0, nrow=length(substrat), ncol=iter)
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
  
  substrat_history[,time] <- unlist(lapply(substrat,FUN=mean))
  rownames(substrat_history) <- names(substrat)
  
  sapply(levels(bac[,3]), function(x,time,tab,bac_history){
    bac_history[[x]][time] <<- tab[x]
  },time=time,tab=table(bac$type),bac_history=bac_history)
  
  
  plot.bacs(time=time, bac=bac, substrate=substrat, growth_vec_history=growth_vec_history, sub_his=substrat_history, bac_his=bac_history, subnam1="glucose", subnam2="co2", subnam3="o2", prodnam="acetate", bac_color=bac_color)
  #plot_list[[time]] <- plot.bacs.cool(bac=bac, time=time, substrat=substrat, sub="glucose", sub2="o2", sub3="acetate", prod="co2", bac_col=bac_color)
}

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
