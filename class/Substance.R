source(file="class/Arena.R")

# Substance inherits from Arena and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

setClass("Substance",
         contains="Arena",
         representation(
           smax     = "numeric",  # substrate start concentration
           diffmat  = "matrix",   # matrix containing concentrations
           name  = "character"   # String describing object
           #diffconst= "numeric",  # diffusion constant
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Substance <- function(n, m, smax, diffmat={}, name, ...){
  if(length(diffmat)==0){
    diffmat = matrix(smax, nrow=n, ncol=m)
  }
  new("Substance", smax=smax, diffmat=diffmat, name=name, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the mediacomposition and gradients

setGeneric("changeSub", function(object, diffmat){standardGeneric("changeSub")})
setMethod("changeSub", "Substance", function(object, diffmat){
  eval.parent(substitute(object@diffmat <- diffmat)) #(pseudo) call by reference implementation
})

#R function for naive diffusion (neighbourhood) of the Substance matrix

setGeneric("diffuseNaiveR", function(object){standardGeneric("diffuseNaiveR")})
setMethod("diffuseNaiveR", "Substance", function(object){
  smat <- object@diffmat
  smatn <- matrix(NA, nrow=dim(smat)[1]+2, ncol=dim(smat)[2]+2) #define environment with boundary conditions
  smatn[2:(dim(smat)[1]+1), 2:(dim(smat)[2]+1)] <- smat #put the values into the environment
  i <- sample(1:dim(smat)[1], dim(smat)[1])
  j <- sample(1:dim(smat)[2], dim(smat)[2])
  for(ic in seq_along(i)){
    for(jc in seq_along(j)){
      neighbours <- c(smatn[ic,jc], 
                      smatn[ic+1,jc], 
                      smatn[ic+2,jc], 
                      smatn[ic+2,jc+1],
                      smatn[ic+2,jc+2], 
                      smatn[ic+1,jc+2],
                      smatn[ic,jc+2],
                      smatn[ic,jc+1])
      minc <- min(neighbours, na.rm = T)
      if(smat[ic,jc] > minc){
        nmin <- which(neighbours == minc)
        if(length(nmin)!=1){
          nmin <- sample(nmin, 1)
        }
        switch(nmin,
              {mn <- mean(c(smat[ic,jc], smat[ic-1,jc-1])); smat[ic,jc] <- mn; smat[ic-1,jc-1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic,jc-1])); smat[ic,jc] <- mn; smat[ic,jc-1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic+1,jc-1])); smat[ic,jc] <- mn; smat[ic+1,jc-1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic+1,jc])); smat[ic,jc] <- mn; smat[ic+1,jc] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic+1,jc+1])); smat[ic,jc] <- mn; smat[ic+1,jc+1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic,jc+1])); smat[ic,jc] <- mn; smat[ic,jc+1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic-1,jc+1])); smat[ic,jc] <- mn; smat[ic-1,jc+1] <- mn},
              {mn <- mean(c(smat[ic,jc], smat[ic-1,jc])); smat[ic,jc] <- mn; smat[ic-1,jc] <- mn})
      }
    }
  }
  eval.parent(substitute(object@diffmat <- smat))
})

#R function for naive diffusion (neighbourhood) of the Substance matrix

diffusion <- cxxfunction(signature(A = "numeric"), body = src_diffusion, plugin="Rcpp")
setGeneric("diffuseNaive", function(object){standardGeneric("diffuseNaive")})
setMethod("diffuseNaive", "Substance", function(object){
  smat <- object@diffmat
  diffusion(list(smat))
  eval.parent(substitute(object@diffmat <- smat))
})