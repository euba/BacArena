source(file="class/Arena.R")

# Substance inherits from Arena and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

setClass("Substance",
         contains="Arena",
         representation(
           smax     = "numeric",  # substrate start concentration
           diffconst= "numeric",  # diffusion constant
           diffmat  = "matrix"   # matrix containing concentrations
           #gradient = "matrix"    # gradient matrix containing infos about how substrate conc is distributed
           # on Arena. Different objects (substrates) may contain different gradients.
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Substance <- function(diffconst, n, m, smax, gradient=matrix(0,0,0), ...){
  diffmat = matrix(smax, nrow=n, ncol=m)
  diffmat[(n/2-n/4):(n/2+n/4), (m/2-m/4):(m/2+m/4)] = 0
  #if(sum(dim(gradient)) != 0){
  #  apply(gradient, function(x, diffmat){
  #    diffmat[x[1],x[2]] <<- x[3]
  #  })
  #}
  new("Substance", smax=smax, diffconst=diffconst, diffmat=diffmat, ...)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("diffuseNaive", function(object){standardGeneric("diffuseNaive")})
setMethod("diffuseNaive", "Substance", function(object){
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

