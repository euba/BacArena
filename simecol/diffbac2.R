#Diffusion and movement

library(simecol)

#Variable Declaration

n <- 20
m <- 20
iter <- 10
enum <- as.vector(sapply(letters, function(x){
  for(i in 1:9){
    a[i] <- paste(x, i, sep="")
  }
  return(a)
}))

#Initiation of agents

s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_D_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  matrix(runif(n*m), nrow=n, ncol=m)
  #matrix(c(runif(n), rep(0, n*m-n)), nrow=n, ncol=m)
  #matrix(c(rep(1, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s

#bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
#bac <- matrix(c(rep(1,2*n), rep(0,n*m-2*n)), nrow=n, ncol=m)
bac <- matrix(c(rep(0,(n*m-2*n)/2), rep(1,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)

#Iteration with rules to apply for each agent

for(time in 1:iter){      
  #diffusion of substrates by mean over neighbourhood
  substrat <- lapply(substrat, function(x){
    anb <- eightneighbours(x)
    nb <- neighbours(x)
    mat <- (anb+x)/(nb+1)
    return(mat)
  })
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  #image(substrat$M_co2_b, zlim=c(0,1), col=colorRampPalette(c("white", "black", "red"))(40))
  #par(mfrow=c(3,3))
  #for(i in 1:9){
  #  image(x[[i]], zlim=c(0,1), 
  #    col=colorRampPalette(c("white", "black", rainbow(14)[i]))(40),
  #    main=paste(names(x)[i], paste("step:", time)))
  #}
  #print(x$M_co2_b)
  #dev.off()
  
  #random movement of bacteria:
  y <- bac
  for (i in 0:(n-1)){
    for(j in 0:(m-1)){
      if(y[i+1,j+1] != 0){
        a <- (i + round(runif(1,-1,1))) %% n 
        b <- (j + round(runif(1,-1,1))) %% m 
        if(bac[a+1,b+1] == 0){
          bac[a+1,b+1] <- 1
          bac[i+1,j+1] <- 0
        }
      }
    }
  }
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  image(bac, col=c("white", "black"))
  #dev.off()
}
