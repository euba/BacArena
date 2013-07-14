library(Rcpp)
library(inline)
library(rbenchmark)


n <- 100
m <- 100
bac <- matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)
bacs <- round(n*m*0.1)
iter <- 20
s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # random substrate
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m) # downstairs substrate
  matrix(c(rep(0,(n*m-2*n)/2), rep(10,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m) # substrate in the middle of our street ohooo
}, n=n, m=m)
names(substrat) <- s


#
# benchmark movement (c++ vs. R)
#

bac <- data.frame(x=round(runif(bacs, min=1, max=n)), y=round(runif(bacs, min=1, max=m)), 
                  type=rep("ecoli", bacs), growth=rep(1, bacs))
bac <- bac[!duplicated(bac[,1:2]),]
rownames(bac) <- 1:nrow(bac) #change indices in data.frame


movement <- cxxfunction(signature(input_matrix = "matrix", input_frame = "data.frame"), body = src_movement, plugin="Rcpp")


R_mov = function(bac, substrat, n, m, iter){
  for(time in 1:iter){        
  #
  #plotting functions
  bacnum <- dim(bac)[1]
  
  mat = matrix(0,n,m) # conversion of data frame into bac matrix
  apply(bac[,1:2], 1, function(x){
    mat[as.numeric(x[1]), as.numeric(x[2])] <<- 1
  })

  
  #
  # FBA
  #
  bacnum <- dim(bac)[1]
  gvec <- 1:bacnum
  
  xr <- round(runif(bacnum, min = -1, max = 1))
  yr <- round(runif(bacnum, min = -1, max = 1))
  for(l in 1:bacnum){
    i <- bac[l,][1,1]
    j <- bac[l,][1,2]
      
    #
    # Movement in R 
    #
    a <- (i + xr[l])
    b <- (j + yr[l])
    if(a == 0){a = n}
    if(b == 0){b = m}
    if(a == n+1){a = 1}
    if(b == m+1){b = 1}
    test <- apply(bac[,1:2], 1, function(x, p){
      if(sum(x==p)==2){
        return(T)
      }else{
        return(F)
      }
    }, p=c(a,b))
    if(!(sum(test)>=1)){ # if empty go for it!
      bac[l,1:2] <- c(a,b)
    }
  }
}
}

cpp_mov <- function(bac, substrat, n, m, iter){
  for(time in 1:iter){        
  #
  #plotting functions
  bacnum <- dim(bac)[1]
  
  #
  # Model of Movement
  #
  tmp <- movement(matrix(0,n,m), bac) # move bacs and return matrix for printing
  bac_img <- tmp$matrix
  bac <- tmp$df
  
  #
  # FBA
  #
  bacnum <- dim(bac)[1]
  gvec <- 1:bacnum
  
  for(l in 1:bacnum){
    i <- bac[l,][1,1]
    j <- bac[l,][1,2]
  }
}
}

res <- benchmark(           
  R_mov(bac, substrat, n, m, iter),
  cpp_mov(bac, substrat, n, m, iter),
  columns=c("test", "replications", "elapsed","relative", "user.self", "sys.self"),
  order="relative",
  replications=1)
print(res) ## show result



##benchmark matrix generation

src_matrix <- '
  int n = as<int>(N);
  int m = as<int>(M);
  Rcpp::NumericMatrix tmp(n,m);
  return tmp;
'
zero_matrix <- cxxfunction(signature(N = "integer", M= "integer"), body = src_matrix, plugin="Rcpp")
res <- benchmark(           
  diffusion2(substrat),
  diffusion(substrat),
  columns=c("test", "replications", "elapsed","relative", "user.self", "sys.self"),
  order="relative",
  replications=1)
print(res) ## show result



## benchmark diffusion

s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_D_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")
substrat <- lapply(s, function(x, n, m){
  matrix(runif(n*m), nrow=n, ncol=m)
  #matrix(c(runif(n), rep(0, n*m-n)), nrow=n, ncol=m)
  #matrix(c(rep(1, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s

#diffusion of substrates by mean over neighbourhood
diffusion2 <- function(substrat){ 
    blob <- lapply(substrat, function(x){
    anb <- eightneighbours(x)
    nb <- neighbours(x)
    mat <- (anb+x)/(nb+1)
    #return(mat)
})}

diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")

res <- benchmark(           
                 diffusion2(substrat),
                 diffusion(substrat),
                 columns=c("test", "replications", "elapsed","relative", "user.self", "sys.self"),
                 order="relative",
                 replications=1)
print(res) ## show result



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
