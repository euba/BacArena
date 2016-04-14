load('/Users/euba/GitRep/pao_reps.RData')

preps = simlist

popsall = list()
for(i in 1:length(preps[[1]]@specs)){
  popsall[[names(preps[[1]]@specs)[i]]] = matrix(0,nrow=length(preps),
                                                      ncol=length(preps[[1]]@simlist))
}
for(i in 1:length(preps)){
  pdat = plotCurves(preps[[i]],retdata=T)$Population
  for(j in 1:nrow(pdat)){
    popsall[[rownames(pdat)[j]]][i,] = pdat[j,]
  }
}
means <- lapply(popsall, function(x){
  return(rbind(apply(x,2,mean),apply(x,2,sd)))
})
popdat = data.frame("means"=means[[1]][1,], "sd"=means[[1]][2,], "type"=rep(names(means)[1],ncol(means[[1]])),
                    "time"=(0:(ncol(means[[1]])-1)))
for(i in 2:length(means)){
  app = data.frame("means"=means[[i]][1,], "sd"=means[[i]][2,], "type"=rep(names(means)[i],ncol(means[[i]])),
                   "time"=(0:(ncol(means[[i]])-1)))
  popdat = rbind(popdat,app)
}
ggplot(popdat, aes(x=time, y=means, col=type)) + 
  geom_errorbar(aes(ymin=means-sd, ymax=means+sd, col=type), width=0.4,size=1.2) +
  geom_line(size=1.2) +
  geom_point(size=5, shape=20) + # 21 is filled circle
  xlab("Time in h") +
  ylab("Number of individuals") +
  scale_color_manual(values=bcol) +
  theme_bw(base_size = 30) +
  theme(#legend.position='none',
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=20)) #15x5           # Position legend in bottom right









for(i in 1:length(preps)){
  pdat = plotCurves(preps[[i]],retdata=T)$Population
  for(j in 1:nrow(pdat)){
    popsall[[rownames(pdat)[j]]][i,] = pdat[j,]
  }
}
popsall = popsall[[1]]
popdat = data.frame('mean'=apply(popsall,2,mean),
                    'se'=apply(popsall,2,sd),
                    'time'=1:ncol(popsall))

ggplot(popdat, aes(x=time, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=1,size=1.2) +
  geom_line(size=1.2) +
  geom_point(size=5, shape=20) + # 21 is filled circle
  xlab("Time in h") +
  ylab("Number of individuals") +
  theme_bw(base_size = 30) +
  theme(legend.position='none',
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=20)) #12x6           # Position legend in bottom right


load('/Users/euba/GitRep/BacArena/poa_model.RData')
trans = model@met_name
names(trans) = paste('EX',gsub('\\[None\\]','',model@met_id),sep='_')

popsall = list()
for(i in 1:length(preps[[1]]@media)){
  popsall[[names(preps[[1]]@media)[i]]] = matrix(0,nrow=length(preps),
                                                 ncol=length(preps[[1]]@simlist))
}
for(i in 1:length(preps)){
  pdat = plotCurves(preps[[i]],retdata=T)$Substances
  for(j in 1:nrow(pdat)){
    popsall[[rownames(pdat)[j]]][i,] = pdat[j,]
  }
}
popsall_sel = popsall[names(which(preps[[1]]@subchange > 100))]
print(length(popsall_sel))
popdat = data.frame('mean'=unlist(lapply(popsall_sel, function(x){return(apply(x,2,mean))})),
                    'se'=unlist(lapply(popsall_sel, function(x){return(apply(x,2,sd))})),
                    'org'=as.vector(sapply(names(popsall_sel),function(x){rep(x,61)})),
                    'time'=rep(1:61,length(popsall_sel)))
popdat$trans = trans[popdat$org]
ggplot(popdat, aes(x=time, y=mean, colour=trans, group=trans)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=1,size=1.2) +
  geom_line(size=1.2) +
  geom_point(size=5, shape=20) + # 21 is filled circle
  xlab("Time in h") +
  ylab("Concentration in mmol per gridcell") +
  theme_bw(base_size = 30) +
  theme(#legend.position='none',
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=20)) #12x6           # Position legend in bottom right


evalsim = preps[[1]]

library(RColorBrewer)
library(NMF)

i = 51
pop = evalsim@simlist[[i]]
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=0.5, col=pop$phenotype+1, pch=19)
#8x8
#evalArena(evalsim,plot_items=c('EX_EC0027'),time=c(20,40,50))
#image(matrix(evalsim@medlist[[i]]['EX_EC0027'][[1]],nrow=evalsim@n,ncol=evalsim@m),axes=FALSE, col=colorRampPalette(c("seashell2","red4"))(100))
aheatmap(matrix(evalsim@medlist[[i]]['EX_EC0027'][[1]],nrow=evalsim@n,ncol=evalsim@m),
         scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,
         color=colorRampPalette(c("seashell2","red4"))(100))#9.5x8

i=51
cmp = 'EX_EC0036'#Succinate
nuse = getPhenoMat(evalsim)[,cmp] #other interesting compounds: 'EX_h(e)','EX_pi(e)','EX_man1p(e)','EX_man6p(e)','EX_chor(e)','EX_succ(e)','EX_fum(e)','EX_for(e)','EX_cit(e)','EX_6pgc(e)','EX_acac(e)','EX_pep(e)','EX_btd_RR(e)','EX_ac(e)','EX_ppa(e)','EX_dha(e)','EX_lac_L(e)','EX_pyr(e)','EX_tyr_L(e)','EX_thym(e)','EX_glyclt(e)'
pop = evalsim@simlist[[i]]
pop$phenotype_n = 1
pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19)

aheatmap(matrix(evalsim@medlist[[i]]['EX_EC0036'][[1]],nrow=evalsim@n,ncol=evalsim@m),
         scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,
         color=colorRampPalette(c("seashell2","seagreen4"))(100))#9.5x8

i=51
cmp = 'EX_EC0029'#Acetate
nuse = getPhenoMat(evalsim)[,cmp] #other interesting compounds: 'EX_h(e)','EX_pi(e)','EX_man1p(e)','EX_man6p(e)','EX_chor(e)','EX_succ(e)','EX_fum(e)','EX_for(e)','EX_cit(e)','EX_6pgc(e)','EX_acac(e)','EX_pep(e)','EX_btd_RR(e)','EX_ac(e)','EX_ppa(e)','EX_dha(e)','EX_lac_L(e)','EX_pyr(e)','EX_tyr_L(e)','EX_thym(e)','EX_glyclt(e)'
pop = evalsim@simlist[[i]]
pop$phenotype_n = 1
pop$phenotype_n[which(pop$phenotype!=0)] = nuse[paste('PAO.',pop$phenotype[which(pop$phenotype!=0)],sep='')] + 1
plot(pop[,c('x','y')],xlim=c(0,evalsim@n),ylim=c(0,evalsim@m),xlab='',ylab='',
     axes=FALSE,cex=0.5, col=pop$phenotype_n, pch=19)

aheatmap(matrix(evalsim@medlist[[i]]['EX_EC0029'][[1]],nrow=evalsim@n,ncol=evalsim@m),
         scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,
         color=colorRampPalette(c("seashell2","royalblue4"))(100))#9.5x8