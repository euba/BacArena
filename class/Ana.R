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

evalSub <- function(simlist, sub, lc="gray90", hc="darkred", file=){
  for(i in seq_along(simlist)){
    ymax = max(simlist[[i]]@media[[sub]]@diffmat)
    subHeat[[i]] <- ggplot(melt(simlist[[i]]@media[[sub]]@diffmat), aes(Var1, Var2)) +
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
    plot(subHeat[[i]])
  }
  #return(subHeat)
}

saveGIF({
  ani.options(nmax = 30)
  evalSub(simlist,"EX_glc(e)")
}, interval = 0.05, movie.name = "bm_demo.gif", ani.width = 600, ani.height = 600)

ggplot(growthc, aes(x=time, y=growthc[,2:ncol(growthc)])) +
  geom_line() +
  theme_bw(base_size=20)
