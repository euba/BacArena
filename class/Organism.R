source(file="class/Arena.R")

# Organism is the class which contains the metabolic model and other features of the organisms in the arena

########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

setClass("Organism",
         contains="Arena",
         representation(
           x="numeric", # x position on grid
           y="numeric", # y position on grid
           model="modelorg", # sybil model object
           type="character", # description of the organism
           objlp="character", # name of the objective fundtion to be optimized
           lpobj="optsol" # sybil optimization object
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(x, y, model, type=mod_desc(model), objlp="", algorithm = "fba", ...){
  if(objlp != ""){ # test if the objective was changed, or still default value
    model <- changeObj(model, objlp)
  }
  lpobj <- optimizeProb(model, ...)
  new("Organism", x=x, y=y, model=model, type=mod_desc(model), objlp=model@react_id[which(model@obj_coef==1)], lpobj=lpobj)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

setGeneric("changeObj", function(mod, newobj) standardGeneric("changeObj"))
setMethod("changeObj", "Organism",
          function(mod, newobj){ # changes the objective function of a model according to the given reaction id
            mod@obj_coef <- rep(0, length(mod@obj_coef))
            mod@obj_coef[which(mod@react_id==newobj)] <- 1
            return(mod)
})

#testing constructor
org1 <- Organism(x=1, y=2, model=mod)
