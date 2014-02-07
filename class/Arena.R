# Arena is the main class from which all others inherit

########################################################################################################
###################################### ARENA CLASS #####################################################
########################################################################################################

setClass("Arena",
         representation(
           n        = "numeric",  # grid size
           m        = "numeric",  # grid size
           iter     = "numeric",  # iterations
           seed     = "numeric",  # reproducible random numbers
           epsilon  = "numeric",  # accuracy in substrate representation
           #orglist  = "list",     # list with Organism objects
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n, m, iter, seed, epsilon){
  new("Arena", n=n, m=m, iter=iter, seed=seed, epsilon=epsilon)
}

