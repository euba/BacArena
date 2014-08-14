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
  grow_accum <- object@fbasol$obj + object@growth
  grow_accum <- grow_accum - dfactor
  eval.parent(substitute(object@growth <- grow_accum)) #(pseudo) call by reference implementation
})

#function for letting bacteria grow by adding the calculated growthrate multiplied with the current growth plus to the already present growth value -> exp growth
#requires as input: organism object

setGeneric("growExp", function(object, dfactor){standardGeneric("growExp")})
setMethod("growExp", "Bac", function(object, dfactor){
  grow_accum <- object@fbasol$obj * object@growth + object@growth
  grow_accum <- grow_accum - dfactor
  eval.parent(substitute(object@growth <- grow_accum)) #(pseudo) call by reference implementation
})


#function to get moore-neighbourhood of a bac

setGeneric("getHood", function(object, population){standardGeneric("getHood")})
setMethod("getHood", "Bac", function(object, population){
  return(pop2mat(population)[,(object@y-1):(object@y+1)][(object@x-1):(object@x+1),])
})


#function to check if the there is some free place in the neighbourhood

setGeneric("emptyHood", function(object, j, population){standardGeneric("emptyHood")})
setMethod("emptyHood", "Bac", function(object, population){
  hood <- getHood(object, population)
  rnd <- sample(1:2, 16, replace=T)-1
  for(i in rnd[1:8]){
    for(j in rnd[9:16])
      if(hood[i,j]==0) return(c(i,j))
  }
  return(-1)
})


# function with the growth model of a bac (biomass growth, replication, death)

setGeneric("growth", function(object, j, population, exp=T){standardGeneric("growth")})
setMethod("growth", "Bac", function(object, j, population, exp=T){
  if(exp) growExp(object, 0.1)
  else growLin(object, 0.1)
  if(object@growth > 2){
    hood <- emptyHood(object, population)
    if(hood != -1){
      doughter <- population@orglist[[j]]
      newg <- population@orglist[[j]]@growth/2
      doughter@growth <- newg
      doughter@x <- hood[1]
      doughter@y <- hood[2]
      eval.parent(substitute(object@growth <- newg))
      eval.parent(substitute(population@orglist <- c(population@orglist, doughter)))
    }
  }
  else if(object@growth < 0.1){
    val.parent(substitute(population@orglist <- population@orglist[[j]]))
  }
})
  


