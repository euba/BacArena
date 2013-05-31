#Diffusion and movement

library(simecol)
library(Rcpp)
library(inline)
library(rbenchmark)

source(file="fba.R")
sbml <- read.sbml("/home/eugen/BacArena/data/ecoli_core.xml")


#Variable Declaration

n <- 100
m <- 100
iter <- 100

bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)

src <- '
  const Rcpp::NumericMatrix  source(A);
  Rcpp::NumericMatrix tmp = Rcpp::clone(source);

   /* initialize random seed: */
  srand (time(NULL));
  /* generate secret number between 1 and 10: */
  //iSecret = rand() % 3;
  
  int n = tmp.nrow();
  int m = tmp.ncol();

  for (int i = 0; i < n; i++){
    for (int j = 1; j < m; j++){
      if(source(i,j) != 0){
        int a = (i + rand() % 3 - 1) % n; // get an integer between [-1:1]
        int b = (j + rand() % 3 - 1) % m; // get an integer between [-1:1]
        if (a == -1) a = n -1; //ugly. we are not satisfied with this...
        if (b == -1) b = m -1;
        if(tmp(a,b) == 0){ // if empty go for it!
          tmp(a,b) = 1;
          tmp(i,j) = 0;
        }
      }
    }
  }
  return tmp;
'
movement <- cxxfunction(signature(A = "numeric"), body = src, plugin="Rcpp")

movement2 <- function(bac, n, m){
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
  return(bac)
}

a <- vector(mode="character")
enum <- as.vector(sapply(letters, function(x){
  for(i in 1:9){
    a[i] <- paste(x, i, sep="")
  }
  return(a)
}))

#Initiation of agents

#
# intial Substrate distribution
#
s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  #matrix(sample(1:100, n*m, replace=T), nrow=n, ncol=m) # integer  in [0:100]
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # doubles in [0:100]
  
  #matrix(c(runif(n), rep(0, n*m-n)), nrow=n, ncol=m)
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m)
  matrix(c(rep(0,(n*m-2*n)/2), rep(100,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s

#
# intial bacterial distribution
#
#bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
#bac <- matrix(c(rep(1,2*n), rep(0,n*m-2*n)), nrow=n, ncol=m)
bac <- matrix(c(rep(0,(n*m-2*n)/2), rep(1,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m)


#Iteration with rules to apply for each agent
mgvec <- vector("numeric")
sgvec <- vector("numeric")
for(time in 1:iter){      
  #diffusion of substrates by mean over neighbourhood
  substrat <- lapply(substrat, function(x){
    anb <- eightneighbours(x)
    nb <- neighbours(x)
    mat <- (anb+x)/(nb+1)
    return(mat)
  })
  #plotting functions:
  jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=1000, height=1000)
  par(mfrow=c(2,2))
  image(substrat$M_glc_b, zlim=c(0,100), col=colorRampPalette(c("white", "black", "red"))(40))
  #par(mfrow=c(3,3))
  #for(i in 1:9){
  #  image(x[[i]], zlim=c(0,1), 
  #    col=colorRampPalette(c("white", "black", rainbow(14)[i]))(40),
  #    main=paste(names(x)[i], paste("step:", time)))
  #}
  #dev.off()
  
  #random movement of bacteria:
  #bac <- movement(bac)  
  #bac <- movement2(bac,n,m)
  bac <- t(apply(bac, 1, sample))  # movement by using random permutation matrix
  
  bacnum <- sum(apply(bac, 1, sum))
  print(bacnum)
  gvec <- 1:bacnum
  k <- 0
  for (i in 1:n){
    for(j in 1:m){
      if(bac[i,j] == 1){ # if there is a Bacterial
        spos <- lapply(substrat, function(x, i, j){
          return(x[i,j])
        },i=i, j=j)
        growth <- fba(spos, sbml$stoch, sbml$lb, sbml$ub, sbml$ex)
        gvec[k]=growth
        k <- k + 1
      }
    }
  }
  image(bac, col=c("white", "black"))
  mgvec[time] <- mean(gvec)
  sgvec[time] <- sd(gvec)
  plot(1:time, mgvec, type="b", xlab="Iteration", ylab="mean replication rate")
  plot(1:time, sgvec, type="b", xlab="Iteration", ylab="standard deviation replication rate")
  
  #plotting functions:
  #jpeg(paste("~/BacArena/plot", paste(enum[time], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
  dev.off()
}
