library(ggplot2)
library(reshape2)
library(animation)

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

evalSub <- function(simlist, sub, lc="gray90", hc="darkred"){
  for(i in seq_along(simlist)){
    ymax = max(simlist[[i]]@media[[sub]]@diffmat)
    subHeat <- ggplot(melt(simlist[[i]]@media[[sub]]@diffmat), aes(Var1, Var2)) +
      geom_tile(aes(fill = value), colour = "white") +
      scale_fill_gradient(low=lc, high=hc, limits=c(0, ymax)) +
      scale_x_continuous(labels = function(x){round(x)}) +
      scale_y_continuous(labels = function(x){round(x)}) +
      ggtitle(paste(sub,"concentration in mM")) +
      ylab("") +
      xlab("") +
      theme(plot.title = element_text(size=rel(2), vjust=0.2),
            legend.position="top",
            legend.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.background = element_blank())
    plot(subHeat)
  }
  #return(subHeat)
}

saveGIF({
  ani.options(nmax = 30)
  evalSub(simlist,"EX_ac(e)")
}, interval = 0.05, movie.name = "bm_demo.gif", ani.width = 600, ani.height = 600)

ggplot(growthc, aes(x=time, y=growthc[,2:ncol(growthc)])) +
  geom_line() +
  theme_bw(base_size=20)

subHis <- t(data.frame(lapply(simlist, function(x){return(unlist(lapply(x@media,function(y){return(mean(y@diffmat))})))})))
rownames(subHis) = 1:nrow(subHis)
subs <- ggplot(melt(as.matrix(sub_his[,1:time])), aes(x=X2, y=value, group=X1, colour=X1)) +
  geom_line(size=1) +
  scale_x_continuous(labels = function(x){floor(x)}) +
  scale_y_continuous(labels = function(x){floor(x)}) +
  ylab("average concentration in mM") +
  xlab("time") +
  theme(axis.line = element_line(colour="black", size=1),
        axis.ticks = element_line(colour="black", size=1),
        axis.text = element_text(size=rel(1.5), colour="black"),
        axis.title.y = element_text(size=rel(1.8), vjust=0.2),
        axis.title.x = element_text(size=rel(1.8), vjust=0.2),
        plot.title = element_text(size=rel(2), vjust=0.2),
        legend.text = element_text(size=rel(1.3), colour="black"),
        legend.title = element_blank(),
        legend.key = element_rect(colour="white", fill="white"),
        panel.grid.major = element_line(colour="grey", size=0.2),
        panel.grid.minor = element_line(colour="grey", size=0.1),
        #panel.border = element_rect(colour="black", size=1),
        panel.background = element_blank())