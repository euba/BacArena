source(file="class/Arena.R")

# Substance inherits from Arena and contains the matrices with concentrations

########################################################################################################
###################################### SUBSTANCE CLASS #################################################
########################################################################################################

setClass("Substance",
         contains="Arena",
         representation(
           diffconst= "numeric",  # diffusion constant
           diffmat  = "matrix",   # matrix containing concentrations
           gradient = "matrix"    # gradient matrix containing infos about how substrate conc is distributed
           # on Arena. Different objects (substrates) may contain different gradients.
         )
)


########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Substance <- function(diffconst, n, m, smax, gradient=matrix(0,0,0)){
  diffmat = matrix(smax, nrow=n, ncol=m)
  if(sum(dim(gradient)) != 0){
    apply(gradient, function(x, diffmat){
      diffmat[x[1],x[2]] <<- x[3]
    })
  }
  new("Substance", diffconst=diffconst, diffmat=diffmat, gradient=gradient)
}

############################ just some tests

sub1 <- Substance(diffconst=5, n=5, m=5, smax=5)
arena1 <- Arena(n=5, m=5, smax=5, seed=5, epsilon=5, iter=5)
