library(ggplot2)
library(reshape2)
library(animation)
library(RColorBrewer)

evalBac <- function(simlist){
  growthc <- data.frame(time=seq_along(simlist))
  snames <- unique(unlist(lapply(simlist[[1]]@orglist,function(x){return(x@type)})))
  for(i in seq_along(simlist)){
    popana <- simlist[[i]]
    #image(popana@occmat, col=terrain.colors(length(levels(as.factor(popana@occmat)))), xaxt='n',yaxt='n')
    image(popana@occmat, xaxt='n',yaxt='n')
    popocc <- table(popana@occmat)
    if("0" %in% names(popocc)){
      popocc = popocc[-which(names(popocc)=="0")]
    }
    for(j in seq_along(popocc)){
      growthc[i,snames[j]] <- popocc[j]
    }
  }
}
evalSubs <- function(simlist){
  for(i in seq_along(simlist)){
    popana <- simlist[[i]]
    #image(popana@occmat, col=c("white",rainbow(simlist[[1]]@orgn)), xaxt='n',yaxt='n')
    fac = length(popana@media)
    par(mfrow=c(10,20), mar=c(0,0,0,0), oma=c(0,0,0,0))
    #for(j in seq_along(popana@media)){
    for(j in 1:fac){
      image(popana@media[[j]]@diffmat, xaxt='n',yaxt='n')
    }
  }
}

saveVideo({
  ani.options(interval=0.35,nmax=500,ani.width=800,ani.height=800)
  evalSubs(simlist)
}, video.name = "mult_sub.mp4", other.opts = "-b 600k")

saveVideo({
  ani.options(interval=0.1,nmax=500,ani.width=1000,ani.height=1000)
  evalBac(simlist, cols="black")
}, video.name = "ecoli_move.mp4", other.opts = "-b 1500k")

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

lapply(simlist[[1]]@media,function(x){print(x@name)})

saveVideo({
  ani.options(interval=0.1,nmax=800,ani.width=1000,ani.height=1000)
  evalSub(simlist,"EX_cpd00039_e0")
}, video.name = "bm_demo3.mp4", other.opts = "-b 900k")

ggplot(growthc, aes(x=time, y=growthc[,2:ncol(growthc)])) +
  geom_line() +
  theme_bw(base_size=20)

subHis <- t(data.frame(lapply(simlist, function(x){return(unlist(lapply(x@media,function(y){return(mean(y@diffmat))})))})))
rownames(subHis) = 1:nrow(subHis)
subs <- ggplot(melt(as.matrix(subHis)), aes(x=Var1, y=value, group=Var2, colour=Var2)) +
  geom_line(size=1) +
  scale_x_continuous(labels = function(x){floor(x)}) +
  scale_y_continuous(labels = function(x){floor(x)}) +
  ylab("average concentration in mM") +
  xlab("time") +
  theme_bw(base_size=20)
#   theme(axis.line = element_line(colour="black", size=1),
#         axis.ticks = element_line(colour="black", size=1),
#         axis.text = element_text(size=rel(1.5), colour="black"),
#         axis.title.y = element_text(size=rel(1.8), vjust=0.2),
#         axis.title.x = element_text(size=rel(1.8), vjust=0.2),
#         plot.title = element_text(size=rel(2), vjust=0.2),
#         legend.text = element_text(size=rel(1.3), colour="black"),
#         legend.title = element_blank(),
#         legend.key = element_rect(colour="white", fill="white"),
#         panel.grid.major = element_line(colour="grey", size=0.2),
#         panel.grid.minor = element_line(colour="grey", size=0.1),
#         #panel.border = element_rect(colour="black", size=1),
#         panel.background = element_blank())

subHis <- t(data.frame(lapply(simlist, function(x){return(unlist(lapply(x@media,function(y){return(mean(y@diffmat))})))})))
rownames(subHis) = 1:nrow(subHis)
msub = melt(as.matrix(subHis))
colnames(msub) = c("Time","Substances","Conc")
subs <- ggplot(msub, aes(x=Time, y=Conc, group=Substances, colour=Substances)) +
  geom_line(size=0.1) +
  scale_x_continuous(labels = function(x){floor(x)}) +
  scale_y_continuous(labels = function(x){floor(x)}) +
  ylab("Average concentration in mM") +
  xlab("Time") +
  theme_bw(base_size=20)
plot(subs)#12x6

bacHis <- matrix(0, nrow=length(simlist), ncol=length(levels(as.factor((simlist[[1]]@occmat))))-1)
for(i in seq_along(simlist)){
  bacocc = table(as.vector(simlist[[i]]@occmat))
  if("0" %in% names(bacocc)){
    bacocc = bacocc[-which(names(bacocc)=="0")]
  }
  bacHis[i,as.numeric(names(bacocc))] = as.vector(bacocc)
}
rownames(bacHis) = 1:nrow(bacHis)
mbac=melt(as.matrix(bacHis))
colnames(mbac)=c("Time","Species","Abundance")
#types="Ecoli"
#mbac$Species=types[mbac$Species]
bacs <- ggplot(mbac, aes(x=Time, y=Abundance, group=as.factor(Species), colour=as.factor(Species))) +
  geom_line(size=0.01) +
  ylab("Number of Bacteria") +
  xlab("Time") +
  theme_bw(base_size=20)
plot(bacs)#12x6

test <- sort(table(simlist[[length(simlist)]]@occmat))#131
trans <- read.csv("P:/MULTISPEC/model_stats.csv", header=T, row.names=1)
names(test) = trans[sapply(as.numeric(names(test)),function(x,specs){return(specs[[x]]@mod_desc)},specs=specs),6]


names(attributes(simlist[[1]]))


evalBac <- function(simlist, cols){
  for(i in seq_along(simlist)){
    dat=melt(simlist[[i]]@occmat)
    if(min(dat$value)==0){
      dat=dat[-which(dat$value==0),]
    }else{
      dat=dat
    }
    bp <- ggplot(data=dat, aes(x=Var1, y=Var2)) +
      geom_point(size=3, aes(colour = factor(value))) +
      scale_x_continuous(limits = c(0, simlist[[i]]@n)) +
      scale_y_continuous(limits = c(0, simlist[[i]]@m)) +
      #scale_colour_manual(values=cols) +
      theme_bw()+
      theme(
        axis.text = element_blank()
        ,axis.ticks = element_blank()
        ,axis.title = element_blank()
        ,plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,legend.position = "none"
      ) 
    jpeg(paste("bacs",i,".jpeg", sep=""))
    plot(bp)
    dev.off()
  }
}
evalBac(simlist)

evalSub <- function(simlist, sub, lc="gray90", hc="darkred"){
  for(i in seq_along(simlist)){
    subconc = ifelse(simlist[[i]]@media[[sub]]@diffmat<0,0,simlist[[i]]@media[[sub]]@diffmat)
    ymax = max(subconc)
    subHeat <- ggplot(melt(subconc), aes(Var1, Var2)) +
      geom_tile(aes(fill = value)) +
      scale_fill_gradient(low=lc, high=hc, limits=c(0, ymax)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw()+
      theme(
        axis.text = element_blank()
        ,axis.ticks = element_blank()
        ,axis.title = element_blank()
        ,plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,legend.position = "none"
      ) 
    plot(subHeat)
  }
  #return(subHeat)
}

evalSubs <- function(simlist, lc=rep("gray90",length(simlist[[1]]@media)),
                     hc=rep("darkred",length(simlist[[1]]@media))){
  #ymax=max(unlist(lapply(simlist, function(x){return(unlist(lapply(x@media,function(y){max(y@diffmat)})))})))
  for(i in seq_along(simlist)){
    popana <- simlist[[i]]
    #image(popana@occmat, col=c("white",rainbow(simlist[[1]]@orgn)), xaxt='n',yaxt='n')
    fac = length(popana@media)
    #par(mfrow=c(11,20), mar=c(0,0,1.5,0), oma=c(0,0,0,0))
    par(mfrow=c(11,20), mar=c(0,0,0,0), oma=c(0,0,0,0))
    #for(j in seq_along(popana@media)){
    for(j in 1:fac){
      ymax=max(unlist(lapply(simlist, function(x,ind){return(max(x@media[[ind]]@diffmat))},ind=j)))
      getPalette = colorRampPalette(c(lc[j],hc[j]))
      subconc = ifelse(popana@media[[j]]@diffmat<0,0,popana@media[[j]]@diffmat)
      #ymax = max(subconc)
      #image(subconc, xaxt='n',yaxt='n', zlim=c(0,ymax), col=getPalette(100), main=popana@media[[j]]@name)
      image(subconc, xaxt='n',yaxt='n', zlim=c(0,ymax), col=getPalette(100))
    }
  }
}

evalSubs(simlist[1:10])

saveVideo({
  ani.options(interval=0.15,nmax=500,ani.width=1200,ani.height=1000)
  evalSubs(simlist,hc=colorRampPalette(brewer.pal(8, "Accent"))(length(simlist[[1]]@media)))
}, video.name = "mult_subs2.mp4", other.opts = "-b 1500k")

saveVideo({
  ani.options(interval=0.15,nmax=500,ani.width=800,ani.height=800)
  evalBac(simlist)
}, video.name = "mult_bac2.mp4", other.opts = "-b 600k")

