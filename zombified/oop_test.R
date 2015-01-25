library(sybil)
library(sybilSBML)
load("data/ecore_model.R")
mod <- model


######################## class definition

setClass("Bac",
         representation(
           x="numeric", # x position on grid
           y="numeric", # y position on grid
           type="character", # description of the metabolic model
           growth="numeric", # growth (biomass) of the individual
           model="modelorg", # sybil model object
           lpobj="optsol" # sybil optimization object
        )
)

######################## constructor definition

Bac <- function(x, y, model, growth, ...){
  lpobj <- optimizeProb(model, algorithm = "fba", ...)
  new("Bac", x=x, y=y, type=mod_desc(model), growth=growth, model=model, lpobj=lpobj)
}

######################## definition of class methods to get the values of the slots (alternatively use Bac@model...)

setGeneric("x", function(x) standardGeneric("x"))
setMethod("x", "Bac", function(x) x@x)
 
setGeneric("y", function(x) standardGeneric("y"))
setMethod("y", "Bac", function(x) x@y)

setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "Bac", function(x) x@type)

setGeneric("growth", function(x) standardGeneric("growth"))
setMethod("growth", "Bac", function(x) x@growth)

setGeneric("model", function(x) standardGeneric("model"))
setMethod("model", "Bac", function(x) x@model)

setGeneric("lpobj", function(x) standardGeneric("lpobj"))
setMethod("lpobj", "Bac", function(x) x@lpobj)

######################## class method definitions

#setGeneric("show", function(object, ...) standardGeneric("show"))
setMethod("show", "Bac", # gives errors, because there are also other show methods?
          function(object){
            cat("Object of class", class(object), "with attributes:\n")
            print(names(attributes(object)))
          }
)

###### fragwuerdige funktion....
setGeneric("grow", function(object, ...) standardGeneric("grow"))
setMethod("grow", "Bac",
          function(object, ...){
            opt <- optimizeProb(model(object), algorithm = "fba", ...)
            object@lpobj <- opt
            object@growth <- object@growth + lp_obj(opt)
            return(object)
          }
)


bac1 <- Bac(x=1, y=2, model=mod, growth=1, solver="clpAPI")
bac1 <- grow(bac1)
