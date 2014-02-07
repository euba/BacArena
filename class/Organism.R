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



