source(file="class/Arena.R")

# Eval is a subclass of Arena containing function to reduce the size of simulations and evalution of results

########################################################################################################
###################################### EVAL CLASS ######################################################
########################################################################################################

setClass("Eval",
         contains="Arena",
         representation(
           medlist="list", # list of simulation results with medium concentration
           simlist="list" # list of simulation results with organism features
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Eval <- function(arena){
  new("Eval", n=arena@n, m=arena@m, tstep=arena@tstep, specs=arena@specs, mediac=arena@mediac, occmat=Matrix(),
      phenotypes=arena@phenotypes, media=arena@media, orgdat=arena@orgdat, medlist=list(), simlist=list())
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the Eval object -> adding simulations

setGeneric("addEval", function(object, arena){standardGeneric("addEval")})
setMethod("addEval", "Eval", function(object, arena){
  eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x){
    return(as.vector(x@diffmat))
  })))
  eval.parent(substitute(object@simlist[[length(object@simlist)+1]] <- arena@orgdat))
  eval.parent(substitute(object@phenotypes <- arena@phenotypes))
})


#function for re-creating an Arena object from Evalution -> usefull, if you want start a new simulation from a previous time point

setGeneric("getArena", function(object, time=length(object@medlist)){standardGeneric("getArena")})
setMethod("getArena", "Eval", function(object, time=length(object@medlist)){
  newmedia <- lapply(object@media, function(x, meds, n, m){
    x@diffmat <- Matrix(meds[[x@name]],nrow=n,ncol=m,sparse=T)
    return(x)
  },meds=object@medlist[[time]], n=object@n, m=object@m)
  
  occdat <- object@simlist[[time]]
  newoccmat <- matrix(0, nrow=object@n, ncol=object@m)
  apply(occdat, 1, function(x){
    newoccmat[x[4],x[5]] <<- x[2]
  })
  
  arena <- Arena(n=object@n, m=object@m, tstep=object@tstep, specs=object@specs, mediac=object@mediac,
        phenotypes=object@phenotypes , media=newmedia, orgdat=occdat, occmat=Matrix(newoccmat,sparse=T))
  return(arena)
})

#function for plotting spatial and temporal change of populations and/or concentrations

setGeneric("evalArena", function(object, plot_items, phencol=F){standardGeneric("evalArena")})
setMethod("evalArena", "Eval", function(object, plot_items, phencol=F){
  for(i in 1:length(object@simlist)){
    subnam <- names(object@medlist[[i]])
    inds <- which(plot_items %in% subnam)
    if(length(plot_items)==1){
      if(length(inds)!=0){
        for(j in 1:length(inds)){
          image(matrix(object@medlist[[i]][[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
        }
      }
    }else if(length(plot_items)<=6){
      par(mfrow=c(2,ceiling(length(plot_items)/2)))
      for(j in 1:length(inds)){
        image(matrix(object@medlist[[i]][[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
      }
    }else{
      par(mfrow=c(3,ceiling(length(plot_items)/3)))
      for(j in 1:length(inds)){
        image(matrix(object@medlist[[i]][[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
      }
    }
    if(plot_items[1]=='population'){
      if(phencol){
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=20,axes=FALSE,cex=0.4,main='Population', col=object@simlist[[i]]$phenotype+1)
      }else{
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=20,axes=FALSE,cex=0.4,main='Population')
      }
    }
  }
})

#show function for class Eval

removeMethod(show, signature(object="Eval"))
setMethod(show, signature(object="Eval"), function(object){
  print(paste('Evaluation results of ',length(object@medlist),' simulation steps.',sep=''))
})