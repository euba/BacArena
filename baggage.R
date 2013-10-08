get_sbml <- function(type){
  if (type=="Bcoli") return(Bcoli_sbml)
  if (type=="ecoli") return(ecoli_sbml)
  if (type=="barkeri") return(barkeri_sbml)
  if(type=="beijerinckii") return(beijerinckii_sbml)
}

get_biomassf <- function(type){
  if (type=="Bcoli") return(Bcoli_biomassf)
  if (type=="ecoli") return(ecoli_biomassf)
  if (type=="barkeri") return(barkeri_biomassf)
  if(type=="beijerinckii") return(beijerinckii_biomassf)
}

get_sub_ex <- function(type){
  if (type=="Bcoli") return(Bcoli_sub_ex)
  if (type=="ecoli") return(ecoli_sub_ex)
  if (type=="barkeri") return(barkeri_sub_ex)
  if(type=="beijerinckii") return(beijerinckii_sub_ex)
}

get_maintenancef <- function(type){
  if (type=="Bcoli") return(Bcoli_maintenancef)
  if (type=="ecoli") return(ecoli_maintenancef)
  if (type=="barkeri") return(barkeri_maintenancef)
  if(type=="beijerinckii") return(beijerinckii_maintenancef)
}

get_lower_bound <- function(type){
  if (type=="Bcoli") return(Bcoli_lower_bound)
  if (type=="ecoli") return(ecoli_lower_bound)
  if (type=="barkeri") return(barkeri_lower_bound)
  if(type=="beijerinckii") return(beijerinckii_lower_bound)
}

get_upper_bound <- function(type){
  if (type=="Bcoli") return(Bcoli_upper_bound)
  if (type=="ecoli") return(ecoli_upper_bound)
  if (type=="barkeri") return(barkeri_upper_bound)
  if(type=="beijerinckii") return(beijerinckii_upper_bound)
}

set_lower_bound <- function(type, substrat){
  if (type=="Bcoli") return(Bcoli_set_lower_bound(substrat))
  if (type=="ecoli") return(ecoli_set_lower_bound(substrat))
  if (type=="barkeri") return(barkeri_set_lower_bound(substrat))
  if(type=="beijerinckii") return(beijerinckii_set_lower_bound(substrat))
}

get_ngam <- function(type, substrat){
  if (type=="Bcoli") return(Bcoli_ngam)
  if (type=="ecoli") return(ecoli_ngam)
  if (type=="barkeri") return(barkeri_ngam)
  if(type=="beijerinckii") return(beijerinckii_ngam)
}

get_gam <- function(type, substrat){
  if (type=="Bcoli") return(Bcoli_gam)
  if (type=="ecoli") return(ecoli_gam)
  if (type=="barkeri") return(barkeri_gam)
  if(type=="beijerinckii") return(beijerinckii_gam)
}

########################################################################################################
###################################### PLOT FUNCTIONS ##################################################
########################################################################################################

plot.bacs.cool <- function(substrate=substrat[[7]], product=substrat[[1]], sub_his=substrat_history,
                      bac_his=bac_history, bac, time,
                      subnam=names(substrat)[7], prodnam=names(substrat)[1]){ #substrate and product as matrices
  data.m <- melt(substrate)
  sub <- ggplot(data.m, aes(X1, X2)) + geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "gray90", high = "darkgreen") +
    scale_x_continuous(labels = function(x){round(x)}) +
    scale_y_continuous(labels = function(x){round(x)}) +
    ggtitle(paste(subnam,"concentration in mM")) +
    ylab("") +
    xlab("") +
    theme(plot.title = element_text(size=rel(2), vjust=0.2),
          legend.position="top",
          legend.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank())
  
  data.m <- melt(product)
  prod <- ggplot(data.m, aes(X1, X2)) + geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "gray90", high = "darkred") +
    scale_x_continuous(labels = function(x){round(x)}) +
    scale_y_continuous(labels = function(x){round(x)}) +
    ggtitle(paste(prodnam,"concentration in mM")) +
    ylab("") +
    xlab("") +
    theme(plot.title = element_text(size=rel(2), vjust=0.2),
          legend.position="top",
          legend.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank())
  print(prod)
  subs <- ggplot(melt(as.matrix(sub_his[,1:time])), aes(x=X2, y=value, group=X1, colour=X1)) +
    geom_line(size=1) +
    scale_x_continuous(labels = function(x){floor(x)}) +
    scale_y_continuous(labels = function(x){floor(x)}) +
    ylab("concentration in mM") +
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
  growth <- ggplot(as.data.frame(cbind(1:length(bac_his), bac_his)), aes(x=V1, y=bac_his)) +
    geom_line(size=1) +
    scale_x_continuous(labels = function(x){floor(x)}) +
    scale_y_continuous(labels = function(x){floor(x)}) +
    ylab("number of bacteria") +
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

  mat = matrix(0,n,m) # conversion of data frame into bac matrix
  apply(bac[,1:2], 1, function(x){
    mat[as.numeric(x[1]), as.numeric(x[2])] <<- 1
  })
  data.m <- melt(mat)
  bacpos <- ggplot(data.m, aes(X1, X2)) + geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "gray90", high = "black") +
    scale_x_continuous(labels = function(x){round(x)}) +
    scale_y_continuous(labels = function(x){round(x)}) +
    ggtitle("bacterial movement") +
    ylab("") +
    xlab("") +
    theme(plot.title = element_text(size=rel(2), vjust=0.2),
          legend.position="none",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank())
  return(list(sub=sub, prod=prod, bacpos=bacpos, growth=growth, subs=subs))
}

plot.bacs <- function(substrate=substrat, growth_vec_history=growth_vec_history,
                      subnam1, subnam2, subnam3, prodnam, sub_his=substrat_history,
                      bac_his=bac_history, bac, time, bac_color=bac_color){ #substrate and product as matrices
  par(mfrow=c(3,4))
  
  image(substrate[[subnam1]], zlim=c(0,max(substrate[[subnam1]])), col=colorRampPalette(c("white", "green"))(40), main=subnam1)
  image(substrate[[subnam2]], zlim=c(0,max(substrate[[subnam2]])), col=colorRampPalette(c("white", "red"))(40), main=subnam2)
  image(substrate[[subnam3]], zlim=c(0,max(substrate[[subnam3]])), col=colorRampPalette(c("white", "steelblue"))(40), main=subnam3)
  image(substrate[[prodnam]], zlim=c(0,max(substrate[[prodnam]])), col=colorRampPalette(c("white", "orange"))(40), main=prodnam)
  
  plot(1:time, sub_his[1,1:time], col=1, pch=1, ylim=c(0,max(sub_his[,1:time])), ylab="concentration", xlab="time") #set max y-value to highest product conentration
  for(i in 2:(dim(sub_his)[1])){
    lines(1:time, sub_his[i,1:time], col=i, pch=i, type="b")
  }
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("top", row.names(sub_his), pch=1:(dim(sub_his)[1]), col=1:dim(sub_his)[1], cex=0.5, pt.cex=0.5, bty="n", y.intersp=0.5)
  
  plot(unlist(bac_his),type="n",xlim=c(1,max(sapply(bac_his,length))), main="growth curve", xlab="time", ylab="number of bacteria")
  mapply(lines,bac_his,col=bac_color, type="b")
  
  mat = matrix(length(bac_color)+1,n,m) # conversion of data frame into bac matrix
  apply(bac[,1:3], 1, function(x){
    mat[as.numeric(x[1]), as.numeric(x[2])] <<- bac_color[x[3]]
  })
  #image(mat, col=c("white", "black"), main="bacterial movement")
  image(mat, col=terrain.colors(length(bac_color)+1), main="bacterial movement")
  boxplot(growth_vec_history, main="growth rates")
}
