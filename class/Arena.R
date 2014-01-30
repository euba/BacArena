# Arena is the main class from which all others inherit


########################################################################################################
###################################### ARENA CLASS #####################################################
########################################################################################################

setClass("Arena",
         representation(
           n        = "numeric",  # grid size
           m        = "numeric",  # grid size
           iter     = "numeric",  # iterations
           smax     = "numeric",  # substrate start concentration
           seed     = "numeric",  # reproducible random numbers
           epsilon  = "numeric"   # accuracy in substrate representation
         )
)



########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n, m, iter, smax, seed, epsilon){
  new("Arena", n=n, m=m, iter=iter, smax=smax, seed=seed, epsilon=epsilon)
}