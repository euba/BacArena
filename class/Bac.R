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

Bac <- function(x, y, model, growth, ...){
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

#testing constructor
bac1 = Bac(x=1, y=1, model=mod, growth=1)
#bac1@model = changeObj(bac1@model, "ATPM")
#org1 <- Organism(x=1, y=2, model=mod)
#org1@model = org1@changeObj(org1@model, "ATPM")

#bac1 = Bac(org1, growth=1) #this does not work, but would be cool if...
