# Substance inherits from Grid and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

#' Structure of the S4 class "Substance"
#' 
#' Structure of the S4 class \code{Substance} representing substances in the environment which can be produced or consumed.
#' @import ReacTran deSolve
#' @export Substance
#' @exportClass Substance
#' @rdname Substance
#'
#' @slot smax A number representing the start concentration of the substance for each grid cell in the environment. 
#' @slot diffmat A sparse matrix containing all concentrations of the substance in the environment.
#' @slot name A character vector representing the name of the substance.
#' @slot id A character vector representing the identifier of the substance.
#' @slot difunc A character vector ("pde","cpp" or "r") describing the function for diffusion.
#' @slot difspeed A number indicating the diffusion rate (given by cm^2/h). Default is set to glucose diffusion in a aqueous solution (6.7e-6 cm^2/s * 3600 s/h = 0.02412 cm^2/h ).
#' @slot advspeed A number indicating the advection rate in x direction (given by cm/h).
#' @slot diffgeometry Diffusion coefficient defined on all grid cells (initially set by constructor).
#' @slot pde Choose diffusion transport reaction to be used (default is diffusion only)
#' @slot boundS A number defining the attached amount of substance at the boundary (Warning: boundary-function must be set in pde!)
setClass("Substance",
         representation(
           smax = "numeric",
           diffmat = "Matrix",
           name = "character",
           id = "character",
           difunc = "character",
           difspeed = "numeric",
           advspeed = "numeric",
           diffgeometry = "list",
           pde = "character",
           boundS = "numeric"
         ),
         prototype(
           difunc = "pde",
           pde="Diff2d",
           boundS = 0
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{Substance}
#' 
#' The constructor to get a new object of class \code{Substance}
#' @export
#' @name Substance-constructor
#' 
#' @param smax A number representing the start concentration of the substance for each grid cell in the environment. 
#' @param difspeed A number indicating the diffusion speed in x and y direction (given by cm^2/h). For more complex setup define Dgrid.
#' @param advspeed A number indicating the advection speed in x direction (given by cm/h). For more complex setup define Vgrid.
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param gridgeometry A list containing grid geometry parameter 
#' @param occupyM A matrix indicating grid cells that are obstacles
#' @param diffmat A matrix with spatial distributed initial concentrations (unit in fmol) (if not set, a homogenous matrix using smax is created)
#' @param template True if diffmat matrix should be used as tempalte only (will be multiplied with smax to obtain cocentrations)
#' @param Dgrid A matrix indicating the diffusion speed in x and y direction (given by cm^2/h).
#' @param Vgrid A number indicating the advection speed in x direction (given by cm/h).
#' @param ... Arguments of \code{\link{Substance-class}}
#' @return Object of class \code{Substance}
Substance <- function(n, m, smax, gridgeometry, difspeed=0.02412, advspeed=0, occupyM, Dgrid=NULL, Vgrid=NULL, diffmat=NULL, template=FALSE, ...){
  if(length(diffmat)==0){
    diffmat <- Matrix::Matrix(smax, ncol=n, nrow=m, sparse=TRUE)
  }else{
      if(template){
        diffmat <- Matrix::Matrix(smax * diffmat, sparse=T)
      }else{
        diffmat <- Matrix::Matrix(diffmat, sparse=T)
      }
  }
  if(ncol(diffmat)!=n && nrow(diffmat)!=m){
    print(paste("arena dimensions  :", n, m))
    print(paste("diffmat dimension:", ncol(diffmat), nrow(diffmat)))
    stop("Dimension of diffmat is invalid")
  } 
  
  new_occupyM <- apply(occupyM, 2, function(x)ifelse(x==0, 1, 0)) # switch 0 and zero for better processing
  if(is.vector(new_occupyM)){new_occupyM=t(as.matrix(new_occupyM))} #conversion to matrix is important if n=1, because otherwise previous apply will output a vector
  if(length(Dgrid)==0) {
    Dgrid <- ReacTran::setup.prop.2D(value=difspeed, grid = gridgeometry$grid2D)
    Dgrid <- lapply(Dgrid, # set obstacle cells to zero
                    function(x) {
                      x[1:nrow(new_occupyM),1:ncol(new_occupyM)] <- x[1:nrow(new_occupyM),1:ncol(new_occupyM)]*new_occupyM; x} )}
  if(length(Vgrid)==0) {
    Vgrid <- ReacTran::setup.prop.2D(value=0, y.value=advspeed, grid = gridgeometry$grid2D)
    Vgrid <- lapply(Vgrid, # set obstacle cells to zero
                    function(x) {x[1:nrow(new_occupyM),1:ncol(new_occupyM)] <- x[1:nrow(new_occupyM),1:ncol(new_occupyM)]*new_occupyM; x} )}
  diffgeometry <- list(Dgrid=Dgrid, Vgrid=Vgrid)
  
  new("Substance", smax=smax, diffmat=diffmat, difspeed=difspeed, diffgeometry=diffgeometry, ...)
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
#' @description The generic function \code{diffuseR} implements the diffusion in the Moore neighbourhood in \code{R}.
#' @export
#' @rdname diffuseR
#'
#' @param object An object of class Substance.
#' @details The diffusion is implemented by iterating through each cell in the grid and taking the cell with the lowest concentration in the Moore neighbourhood to update the concentration of both by their mean.
#' @seealso \code{\link{Substance-class}} and \code{\link{diffusePDE}}
setGeneric("diffuseR", function(object){standardGeneric("diffuseR")})
#' @export
#' @rdname diffuseR
setMethod("diffuseR", "Substance", function(object){
  smat <- matrix(object@diffmat, nrow=nrow(object@diffmat), ncol=ncol(object@diffmat))
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

#' @title Function for diffusion of the Substance matrix
#'
#' @description The generic function \code{diffusePDE} implements the diffusion by the solving diffusion equation.
#' @export
#' @rdname diffusePDE
#'
#' @param object An object of class Substance.
#' @param init_mat A matrix with values to be used by the diffusion.
#' @param gridgeometry A list specifying the geometry of the Arena
#' @param tstep A numeric value giving the time step of integration
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in simEnv() by the function estimate_lrw())
#' @details Partial differential equation is solved to model 2d diffusion process in the arena.
#' @seealso \code{\link{Substance-class}} and \code{\link{diffuseR}}
setGeneric("diffusePDE", function(object, init_mat, gridgeometry, lrw=NULL, tstep){standardGeneric("diffusePDE")})
#' @export
#' @rdname diffusePDE
setMethod("diffusePDE", "Substance", function(object, init_mat, gridgeometry, lrw=NULL, tstep){
  if(is.null(lrw)){
    lrw=estimate_lrw(gridgeometry$grid2D$x.N, gridgeometry$grid2D$y.N)}
  solution <- deSolve::ode.2D(y = init_mat, func = get(object@pde), times=c(1,1+tstep), parms = c(gridgeometry=gridgeometry, diffgeometry=object@diffgeometry, boundS=object@boundS),
                     dimens = c(gridgeometry$grid2D$x.N, gridgeometry$grid2D$y.N), method="lsodes", lrw=lrw)#160000
  diff_mat <- matrix(data=solution[2,][-1], ncol=ncol(init_mat), nrow=nrow(init_mat))
  return(diff_mat)
})



#show function for class Substance

setMethod(show, signature(object="Substance"), function(object){
  print(paste('Compound ',object@name,' of class Substance with a total concentration of ',
              sum(object@diffmat)/length(object@diffmat),' mmol per gridcell.',sep=''))
})
