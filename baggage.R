get_sbml <- function(type){
  switch(type,
         ecoli = return(ecoli_sbml),
         barkeri = return(barkeri_sbml))
}

get_biomassf <- function(type){
  switch(type,
         ecoli = return(ecoli_biomassf),
         barkeri = return(barkeri_biomassf))
}

get_sub_ex <- function(type){
  switch(type,
         ecoli = return(ecoli_sub_ex),
         barkeri = return(barkeri_sub_ex))
}

get_maintenancef <- function(type){
  switch(type,
         ecoli = return(ecoli_maintenancef),
         barkeri = return(barkeri_maintenancef))
}

get_lower_bound <- function(type){
  switch(type,
         ecoli = return(ecoli_lower_bound),
         barkeri = return(barkeri_lower_bound))
}

get_upper_bound <- function(type){
  switch(type,
         ecoli = return(ecoli_upper_bound),
         barkeri = return(barkeri_upper_bound))
}

set_lower_bound <- function(type, substrat){
  switch(type,
         ecoli = return(ecoli_set_lower_bound(substrat)),
         barkeri = return(barkeri_set_lower_bound(substrat)))
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
  legend("top", row.names(sub_his), pch=1:(dim(sub_his)[1]), col=1:dim(sub_his)[1])
  
  plot(1:length(bac_his), bac_his, type="b", main="growth curve", xlab="time", ylab="number of bacteria")
  
  mat = matrix(0,n,m) # conversion of data frame into bac matrix
  apply(bac[,1:2], 1, function(x){
    mat[as.numeric(x[1]), as.numeric(x[2])] <<- 1
  })
  image(mat, col=c("white", "black"), main="bacterial movement")
}

