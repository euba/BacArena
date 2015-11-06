# Substance inherits from Grid and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

#' Structure of the S4 class "Substance"
#' 
#' Structure of the S4 class \code{Substance} representing substances in the environment which can be produced or consumed.
#'
#' @slot smax A number representing the start concentration of the substance for each grid cell in the environment. 
#' @slot diffmat A sparse matrix containing all concentrations of the substance in the environment.
#' @slot name A character vector representing the name of the substance.
#' @slot difunc A character vector ("pde","cpp" or "r") describing the function for diffusion.
#' @slot difspeed A number indicating the diffusion speed (given by number of cells per iteration).
setClass("Substance",
         representation(
           smax = "numeric",
           diffmat = "Matrix",
           name = "character",
           difunc = "character",
           difspeed = "numeric"
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Substance <- function(n, m, smax, diffmat={}, name, difunc="pde", difspeed=1, ...){
  if(length(diffmat)==0){
    diffmat = Matrix(smax, nrow=n, ncol=m, sparse=T)
  }
  new("Substance", smax=smax, diffmat=diffmat, name=name, difunc=difunc, difspeed=difspeed, ...)
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

#' @title Function for naive diffusion (neighbourhood) of the Substance matrix
#'
#' @description The generic function \code{diff} implements the diffusion in the Moore neighbourhood in \code{R}.
#'
#' @param object An object of class Substance.
#' @details The diffusion is implemented by iterating through each cell in the grid and taking the cell with the lowest concentration in the Moore neighbourhood to update the concentration of both by their mean.
#' @seealso \code{\link{Substance-class}} and \code{\link{diffuseCpp}}
#' @examples
#' \dontrun{
#' sub <- Substance(n=20,m=20,smax=40,name='test',difunc='r') #initialize test substance
#' diffuseR(sub)
#' }
setGeneric("diffuseR", function(object){standardGeneric("diffuseR")})
setMethod("diffuseR", "Substance", function(object){
  smat <- as(object@diffmat, "matrix")
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
  eval.parent(substitute(object@diffmat <- as(smat, "sparseMatrix")))
})


setGeneric("diffusePDE", function(object, init_mat, geometry){standardGeneric("diffusePDE")})
setMethod("diffusePDE", "Substance", function(object, init_mat, geometry){
  #init_mat <- as.matrix(object@diffmat)
  solution <- ode.2D(y = init_mat, func = Diff2d, t = 1:2, parms = c(geometry=geometry, D=object@difspeed),
                   dim = c(geometry$grid2D$x.N, geometry$grid2D$y.N), method="lsodes", lrw=geometry$lrw)
  diff_mat <- matrix(data=solution[2,][-1], ncol=ncol(init_mat), nrow=nrow(init_mat))
  return(diff_mat)
})



#show function for class Substance

setMethod(show, signature(object="Substance"), function(object){
  print(paste('Compound ',object@name,' of class Substance with a total concentration of ',
              round(sum(object@diffmat)/length(c(object@diffmat)),4),' mmol per gridcell.',sep=''))
})
