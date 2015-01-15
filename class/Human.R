source(file="class/Organism.R")

# Human is a subclass of Organism containing human specific features

########################################################################################################
###################################### HUMAN CLASS #####################################################
########################################################################################################

setClass("Human",
         contains="Organism",
         representation(
           deathrate="numeric", # factor by which growth is reduced
           duplirate="numeric", # grow cut-off for test of duplication
           speed="integer", # speed by which bacterium is moving (given by cell per iteration)
           growthlimit="numeric",
           lyse="logical", #flag indicating, if bacterial lysis should be implemented
           feat="list", #list containing conditional features for the object (contains at the momement only biomass components for lysis)
           growtype="character", # functional type for growth (linear or exponential)
           chem='character' # name of substance which is the chemotaxis attractant
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Human <- function(model, deathrate, duplirate, speed=2, growthlimit, growtype, lyse=F, chem='', ...){
  feat <- list()
  if(lyse){
    stoch <- S(model)[,which(model@obj_coef==1)] #find stochiometry of biomass components
    names(stoch) <- met_id(model)
    feat[["biomass"]] <- stoch[-which(stoch==0)]
  }
  new("Bac", Organism(model=model, ...), speed=as.integer(speed), deathrate=deathrate, duplirate=duplirate,
      growthlimit=growthlimit, lyse=lyse, feat=feat, growtype=growtype, chem=chem)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################
