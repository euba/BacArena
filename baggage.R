get_sbml <- function(type){
  if (type=="ecoli") return(ecoli_sbml)
  if (type=="barkeri") return(barkeri_sbml)
}

get_biomassf <- function(type){
  if (type=="ecoli") return(ecoli_biomassf)
  if (type=="barkeri") return(barkeri_biomassf)
}

get_sub_ex <- function(type){
  if (type=="ecoli") return(ecoli_sub_ex)
  if (type=="barkeri") return(barkeri_sub_ex)
}

get_maintenancef <- function(type){
  if (type=="ecoli") return(ecoli_maintenancef)
  if (type=="barkeri") return(barkeri_maintenancef)
}

get_lower_bound <- function(type){
  if (type=="ecoli") return(ecoli_lower_bound)
  if (type=="barkeri") return(barkeri_lower_bound)
}

get_upper_bound <- function(type){
  if (type=="ecoli") return(ecoli_upper_bound)
  if (type=="barkeri") return(barkeri_upper_bound)
}

set_lower_bound <- function(type, substrat){
  if (type=="ecoli") return(ecoli_set_lower_bound(substrat))
  if (type=="barkeri") return(barkeri_set_lower_bound(substrat))
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

plot.bacs <- function(substrate=substrat$glucose, product=substrat$acetate, sub_his=substrat_history,
                      bac_his=bac_history, bac, time){ #substrate and product as matrices
  par(mfrow=c(3,2))
  
  image(substrate, zlim=c(0,max(substrate)), col=colorRampPalette(c("white", "green"))(40), main="substrate concentration")
  image(product, zlim=c(0,max(product)), col=colorRampPalette(c("white", "red"))(40), main="product concentration")
  
  plot(1:time, sub_his[1,1:time], col=1, pch=1, ylim=c(0,max(sub_his[,1:time])), ylab="concentration", xlab="time") #set max y-value to highest product conentration
  for(i in 2:(dim(sub_his)[1])){
    lines(1:time, sub_his[i,1:time], col=i, pch=i, type="b")
  }
  plot(0,0, col="white")
  legend("top", row.names(sub_his), pch=1:(dim(sub_his)[1]), col=1:dim(sub_his)[1], cex=0.64)
  
  
  plot(1:length(bac_his), bac_his, type="b", main="growth curve", xlab="time", ylab="number of bacteria")
  
  mat = matrix(10,n,m) # conversion of data frame into bac matrix
  apply(bac[,1:3], 1, function(x){
    if(x[3] == "barkeri") mat[as.numeric(x[1]), as.numeric(x[2])] <<- 1
    if(x[3] == "ecoli") mat[as.numeric(x[1]), as.numeric(x[2])] <<- 5
  })
  #image(mat, col=c("white", "black"), main="bacterial movement")
  image(mat, col=terrain.colors(10), main="bacterial movement")
  print(mat)
}
