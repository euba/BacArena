source(file="class/Organism.R")

# Bac is a subclass of Organism containing bacteria specific features

########################################################################################################
###################################### BAC CLASS #######################################################
########################################################################################################

setClass("Bac",
         contains="Organism",
         representation(
           growth="numeric" # growth (biomass) of the individual
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Bac <- function(x, y, model, growth=1, ...){
  new("Bac", Organism(x=x, y=y, model=model, ...), growth=growth)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

# setMethod("show", "Bac",
#           function(object){
#             cat("Object of class", class(object), "with attributes:\n")
#             print(names(attributes(object)))
#           }
# )

#function for letting bacteria grow by adding the calculated growthrate to the already present growth value -> linear growth
#requires as input: organism object

setGeneric("growLin", function(object, dfactor){standardGeneric("growLin")})
setMethod("growLin", "Bac", function(object, dfactor){
  grow_accum <- spec@fbasol$obj + object@growth
  grow_accum <- grow_accum - dfactor
  eval.parent(substitute(object@growth <- grow_accum)) #(pseudo) call by reference implementation
})

#function for letting bacteria replicate to the specified position by returning doughter cell

setGeneric("repli", function(object, x, y, bd=F){standardGeneric("repli")})
setMethod("repli", "Bac", function(object, x, y, bd=F){
  doughter <- object
  newg <- object@growth/2
  doughter@growth <- newg
  doughter@x <- x
  doughter@y <- y
  eval.parent(substitute(object@growth <- newg))
  if(!bd){
    return(doughter)
  }
})