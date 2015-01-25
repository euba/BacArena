source(file="Arena.R")

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
      phenotypes=arena@phenotypes, media=arena@media, orgdat=arena@orgdat, medlist=list(), simlist=list(), stir=arena@stir)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("medlist", function(object){standardGeneric("medlist")})
setMethod("medlist", "Eval", function(object){return(object@medlist)})
setGeneric("simlist", function(object){standardGeneric("simlist")})
setMethod("simlist", "Eval", function(object){return(object@simlist)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the Eval object -> adding or replacing simulations

setGeneric("addEval", function(object, arena, replace=F){standardGeneric("addEval")})
setMethod("addEval", "Eval", function(object, arena, replace=F){
  if(!replace){
    eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x){
      return(as.vector(x@diffmat))
    })))
    eval.parent(substitute(object@simlist[[length(object@simlist)+1]] <- arena@orgdat))
    eval.parent(substitute(object@phenotypes <- arena@phenotypes))
  }else{
    eval.parent(substitute(object@medlist[[length(object@medlist)]] <- lapply(arena@media, function(x){
      return(as.vector(x@diffmat))
    })))
    eval.parent(substitute(object@simlist[[length(object@simlist)]] <- arena@orgdat))
    eval.parent(substitute(object@phenotypes <- arena@phenotypes))
  }
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
        phenotypes=object@phenotypes , media=newmedia, orgdat=occdat, occmat=Matrix(newoccmat,sparse=T), stir=object@stir)
  return(arena)
})

#function for plotting spatial and temporal change of populations and/or concentrations

setGeneric("evalArena", function(object, plot_items='population', phencol=F, retdata=F){standardGeneric("evalArena")})
setMethod("evalArena", "Eval", function(object, plot_items='population', phencol=F, retdata=F){
  old.par <- par(no.readonly = TRUE)
  if(retdata){
    retlist = list()
    for(i in 1:length(plot_items)){
      retlist[[i]] = list()
    }
    names(retlist) = plot_items
  }
  for(i in 1:length(object@simlist)){
    subnam <- names(object@medlist[[i]])
    inds <- which(subnam %in% plot_items)
    if(length(plot_items)==1){
      if(length(inds)!=0){
        for(j in 1:length(inds)){
          if(retdata){
            retlist[[subnam[inds[j]]]][[j]] = matrix(object@medlist[[i]][[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
          }
          image(matrix(object@medlist[[i]][[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
        }
      }
    }else if(length(plot_items)<=6){
      par(mfrow=c(2,ceiling(length(plot_items)/2)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[j]] = matrix(object@medlist[[i]][[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(object@medlist[[i]][[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
      }
    }else{
      par(mfrow=c(3,ceiling(length(plot_items)/3)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[j]] = matrix(object@medlist[[i]][[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(object@medlist[[i]][[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
      }
    }
    if(plot_items[1]=='population'){
      if(retdata){
        retlist[['population']][[j]] = object@simlist[[i]]
      }
      if(phencol){
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=20,axes=FALSE,cex=1,main='Population', col=object@simlist[[i]]$phenotype+1)
      }else{
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=20,axes=FALSE,cex=1,main='Population', col=object@simlist[[i]]$type)
      }
    }
  }
  if(retdata){
    return(retlist)
  }
  par(old.par)
})

#function for plotting the overall change as curves

setGeneric("plotCurves", function(object, medplot=object@mediac, retdata=F, remove=F){standardGeneric("plotCurves")})
setMethod("plotCurves", "Eval", function(object, medplot=object@mediac, retdata=F, remove=F){
  old.par <- par(no.readonly = TRUE)
  growths <- matrix(0, nrow=length(object@specs), ncol=length(object@simlist))
  rownames(growths) = names(object@specs)
  subs <- matrix(0, nrow=length(medplot), ncol=length(object@simlist))
  rownames(subs) = medplot
  for(i in 1:length(object@simlist)){
    simdat <- object@simlist[[i]]
    count <- table(simdat[,'type'])
    for(j in 1:length(count)){
      growths[j,i] <- count[j]
    }
    subdat <- object@medlist[[i]]
    for(j in 1:length(medplot)){
      subs[medplot[j],i] <- sum(subdat[[medplot[j]]])/(object@n*object@m)
    }
  }
  if(remove){
    test <- which(apply(subs,1,function(x){return(sum(ifelse(x==x[1],T,F))==length(x))}))
    if(length(test)!=0){
      subs <- subs[-test,]
    }
  }
  par(mfrow=c(2,1))
  times = c(1:length(object@simlist)*object@tstep)
  plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(growths)),
       type='n', xlab='time in h', ylab='number of individuals on grid',
       main='Population')
  for(i in 1:nrow(growths)){
    lines(times, growths[i,], col=i)
  }
  legend('bottom',legend=rownames(growths),col=1:nrow(growths),cex=0.4/log10(nrow(growths)+1),lwd=1)
  plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(subs)),
       type='n', xlab='time in h', ylab='concentration in mmol per gridcell',
       main='Substance concentrations')
  for(i in 1:nrow(subs)){
    lines(times, subs[i,], col=i)
  }
  legend('right',legend=rownames(subs),col=1:nrow(subs),cex=0.4/log10(nrow(subs)+1),lwd=1)
  if(retdata){
    return(list('Population'=growths,'Substances'=subs))
  }
  par(old.par)
})

#function for getting a matrix of phenotypes from the dataset

setGeneric("getPhenoMat", function(object){standardGeneric("getPhenoMat")})
setMethod("getPhenoMat", "Eval", function(object){
  numphens <- unlist(lapply(object@phenotypes,function(x){return(length(x))}))
  phentypes <- vector()
  for(i in seq_along(numphens)){
    phentypes <- c(phentypes, rep(names(numphens)[i],numphens[i]))
  }
  phentypes <- as.factor(phentypes)
  phenmat <- matrix(0, nrow=length(phentypes), ncol=length(object@mediac))
  colnames(phenmat) <- object@mediac
  rownames(phenmat) <- phentypes
  pind <- 0
  for(i in 1:length(object@phenotypes)){
    for(j in 1:numphens[i]){
      pvec <- object@phenotypes[[i]][[j]]
      pind <- pind + 1
      phenmat[pind, names(pvec)] <- pvec
    }
  }
  phenmat <- ifelse(phenmat==-1,2,phenmat)
  return(phenmat)
})

#function for mining/analyzing phenotypes which occured on the arena

setGeneric("minePheno", function(object, plot_type="pca"){standardGeneric("minePheno")})
setMethod("minePheno", "Eval", function(object, plot_type="pca"){
  phenmat <- getPhenoMat(object)
  if(nrow(phenmat)<=1){
    stop('not enough phenotypes to analyze.')
  }
  old.par <- par(no.readonly = TRUE)
  pcount <- as.vector(table(rownames(phenmat)))
  plabs <- vector()
  plabs2 <- vector()
  for(i in 1:length(pcount)){
    plabs <- c(plabs,paste(rep(levels(as.factor(rownames(phenmat)))[i],pcount[i]),1:pcount[i],sep='_'))
    plabs2 <- c(plabs2,paste(i,1:pcount[i],sep='_'))
  }
  par(mfrow=c(1,1))
  if(plot_type=="pca"){
    phenpca <- prcomp(phenmat)
    plot(phenpca$x[,1:2], xlab='PC1', ylab='PC2', pch=20, cex=0.8, col=as.numeric(as.factor(rownames(phenpca$x))))
    text(phenpca$x[,1:2], labels=plabs2, col=as.numeric(as.factor(rownames(phenpca$x))), cex=0.7)
    typ <- names(object@specs)
    for(i in 1:length(object@specs)){
      segments(phenpca$x[which(rownames(phenpca$x)==typ[i]),1],phenpca$x[which(rownames(phenpca$x)==typ[i]),2],
               mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),1]),mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),2]),col=i)
      points(mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),1]),mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),2]),col=i,pch=15,cex=1.5)
    }
    legend('topright',legend=names(object@specs),col=1:length(object@specs),cex=0.9,lwd=4)
  }
  if(plot_type=="hclust"){
    rownames(phenmat) <- plabs
    plot(hclust(dist(phenmat)))
    par(old.par)
  }
})

#show function for class Eval

removeMethod(show, signature(object="Eval"))
setMethod(show, signature(object="Eval"), function(object){
  print(paste('Evaluation results of ',length(object@medlist),' simulation steps.',sep=''))
})