library(ggplot2)

growthc <- data.frame(time=seq_along(simlist))
snames <- unique(unlist(lapply(simlist[[1]]@orglist,function(x){return(x@type)})))
for(i in seq_along(simlist)){
  popana <- simlist[[i]]
  image(popana@occmat)
  popocc <- table(popana@occmat)[-1]
  for(j in seq_along(popocc)){
    growthc[i,snames[j]] <- popocc[j]
  }
}

ggplot(growthc, aes(x=time, y=growthc[,2:ncol(growthc)])) +
  geom_line() +
  
  theme_bw(base_size=20)
