library(Rcpp)
library(inline)
library(rbenchmark)


n <- 2000
m <- 2000
bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)



## benchmark diffusion

s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_D_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  matrix(runif(n*m), nrow=n, ncol=m)
  #matrix(c(runif(n), rep(0, n*m-n)), nrow=n, ncol=m)
  #matrix(c(rep(1, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s

#diffusion of substrates by mean over neighbourhood
substrat <- lapply(substrat, function(x){
  anb <- eightneighbours(x)
  nb <- neighbours(x)
  mat <- (anb+x)/(nb+1)
  return(mat)
})




## benchmark bug movement (cpp vs. R loops)

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

res <- benchmark(movement2(bac,n,m),
                 movement(bac),
                 columns=c("test", "replications", "elapsed","relative", "user.self", "sys.self"),
                 order="relative",
                 replications=1)
print(res) ## show result



## benchmark matrix operations

add <- function(y,n,m){
  for (i in 0:(n-1)){
    for(j in 0:(m-1)){
      y [i,j] <- y[i,j] + y[i,j]
    }
  }
  return(y)
}

add2 <- function(y){
  return(y+y)
}



src <- '
  Rcpp::NumericMatrix tmp(A);
  int n = tmp.nrow();
  int m = tmp.ncol();
  for (int i = 0; i < n; i++){
    for (int j = 1; j < m; j++){
      //if(m[i+1][j+1] != 0){
        tmp(i,j) += tmp(i,j);
      //}
    }
  }
  return tmp;
'
fun <- cxxfunction(signature(A = "numeric"), body = src, plugin="Rcpp")

bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
library(rbenchmark)
res <- benchmark(add(bac,n,m),
                 add2(bac),
                 fun(bac),
                 columns=c("test", "replications", "elapsed","relative", "user.self", "sys.self"),
                 order="relative",
                 replications=1)
print(res) ## show result
