# Substance inherits from Grid and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

setClass("Substance",
         representation(
           smax = "numeric", # substrate start concentration
           diffmat = "Matrix", # matrix containing concentrations -> sparse Matrix
           name = "character", # String describing object
           difunc = "character", # String describing the function for diffusion
           difspeed = "numeric" # String describing the function for diffusion
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Substance <- function(n, m, smax, diffmat={}, name, difunc="cpp", difspeed=1, ...){
  if(length(diffmat)==0){
    diffmat = Matrix(smax, nrow=n, ncol=m, sparse=T)
  }
  new("Substance", smax=smax, diffmat=diffmat, name=name, difunc=difunc, difspeed=round(difspeed), ...)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("smax", function(object){standardGeneric("smax")})
setMethod("smax", "Substance", function(object){return(object@smax)})
setGeneric("diffmat", function(object){standardGeneric("diffmat")})
setMethod("diffmat", "Substance", function(object){return(object@diffmat)})
setGeneric("name", function(object){standardGeneric("name")})
setMethod("name", "Substance", function(object){return(object@name)})
setGeneric("difunc", function(object){standardGeneric("difunc")})
setMethod("difunc", "Substance", function(object){return(object@difunc)})
setGeneric("difspeed", function(object){standardGeneric("difspeed")})
setMethod("difspeed", "Substance", function(object){return(object@difspeed)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################


#R function for naive diffusion (neighbourhood) of the Substance matrix

setGeneric("diffuseR", function(object){standardGeneric("diffuseR")})
setMethod("diffuseR", "Substance", function(object){
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
  #eval.parent(substitute(object@diffmat <- smat))
})

#C++ function for naive diffusion (neighbourhood) of the Substance matrix

setGeneric("diffuseCpp", function(object){standardGeneric("diffuseCpp")})
setMethod("diffuseCpp", "Substance", function(object){
  submat <- as.matrix(object@diffmat)
  diffuseNaiveCpp(submat, donut=FALSE)
  eval.parent(substitute(object@diffmat <- Matrix(submat, sparse=T)))
})

#show function for class Substance

setMethod(show, signature(object="Substance"), function(object){
  print(paste('Compound ',object@name,' of class Substance with a total concentration of ',
              round(sum(object@diffmat)/length(c(object@diffmat)),4),' mmol per gridcell.',sep=''))
})
