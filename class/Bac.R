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

setGeneric("growLin", function(object){standardGeneric("growLin")})
setMethod("growLin", "Bac", function(object){
  #optimizeLP(object) #this does not work, because you can't change an object twice by the call by reference in one function
  grow_val <- object@lpobj@lp_obj
  grow_accum <- grow_val + object@growth
  eval.parent(substitute(object@growth <- grow_accum)) #(pseudo) call by reference implementation
})
