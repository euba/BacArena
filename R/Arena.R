# Arena is the class, that contains a list of of all organisms in the Arena

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

setClass("Arena",
         representation(
           orgdat="data.frame", # data frame of individuals in the Arena
           specs="list", # list of organism types in the Arena
           media="list", # media composition of mixed organisms
           phenotypes="list", # list of unique phenotypes of the individuals
           mediac="character",
           occmat="Matrix", # occupacy matrix (showing which cells have bacs) -> sparse Matrix
           tstep="numeric", # time steps per iteration
           stir="logical", # boolean variable indicating if environment should be stirred
           n="integer",  # grid size
           m="integer"  # grid size
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n,m,tstep=1,orgdat=data.frame(growth=numeric(0),type=integer(0),phenotype=integer(0),x=integer(0),y=integer(0)),
                  specs=list(),media=list(),mediac=character(),phenotypes=list(),occmat=Matrix(0L,nrow=n,ncol=m,sparse=T),stir=F){
  new("Arena", n=as.integer(n), m=as.integer(m), tstep=tstep, orgdat=orgdat, specs=specs,
      media=media, mediac=mediac, phenotypes=phenotypes, occmat=occmat, stir=stir)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("orgdat", function(object){standardGeneric("orgdat")})
setMethod("orgdat", "Arena", function(object){return(object@orgdat)})
setGeneric("specs", function(object){standardGeneric("specs")})
setMethod("specs", "Arena", function(object){return(object@specs)})
setGeneric("media", function(object){standardGeneric("media")})
setMethod("media", "Arena", function(object){return(object@media)})
setGeneric("phenotypes", function(object){standardGeneric("phenotypes")})
setMethod("phenotypes", "Arena", function(object){return(object@phenotypes)})
setGeneric("mediac", function(object){standardGeneric("mediac")})
setMethod("mediac", "Arena", function(object){return(object@mediac)})
setGeneric("occmat", function(object){standardGeneric("occmat")})
setMethod("occmat", "Arena", function(object){return(object@occmat)})
setGeneric("tstep", function(object){standardGeneric("tstep")})
setMethod("tstep", "Arena", function(object){return(object@tstep)})
setGeneric("stir", function(object){standardGeneric("stir")})
setMethod("stir", "Arena", function(object){return(object@stir)})
setGeneric("n", function(object){standardGeneric("n")})
setMethod("n", "Arena", function(object){return(object@n)})
setGeneric("m", function(object){standardGeneric("m")})
setMethod("m", "Arena", function(object){return(object@m)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

# Add Individuals to the arena


#' Function with to add individuals to the arena
#'
#' @param object An object of class Arena
#' @param specI An object of class Organism
#' @param amount Number of organisms to add
#' @param x x positions of individuals on the grid
#' @param y y positions of individuals on the grid
#' @param growth Starting biomass of organisms
#' @examples
#' NULL
setGeneric("addOrg", function(object, specI, amount, x=NULL, y=NULL, growth=1, ...){standardGeneric("addOrg")})
setMethod("addOrg", "Arena", function(object, specI, amount, x=NULL, y=NULL, growth=1, ...){
  if(amount+sum(object@occmat) > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  n <- object@n
  m <- object@m
  spectype <- specI@type
  newoccmat <- object@occmat
  neworgdat <- object@orgdat
  newspecs <- object@specs
  newphens <- object@phenotypes[[spectype]]
  newspecs[[spectype]] <- specI
  type <- which(names(newspecs)==spectype)
  
  if(length(newphens)!=0){
    ptype <- as.integer(checkPhen(object, specI))
    newphens <- object@phenotypes[[spectype]]
  }else{
    newphens[[1]] <- getPhenotype(specI)
    ptype=as.integer(1)
  }
  lastind <- nrow(object@orgdat)
  if(length(x*y)==0){
    cmbs = expand.grid(1:n,1:m)
    rownames(cmbs) = paste(cmbs[,1],cmbs[,2],sep='_')
    taken <- paste(object@orgdat$x,object@orgdat$y,sep='_')
    if(length(taken)!=0){
      cmbs <- cmbs[-which(rownames(cmbs) %in% taken),]
    }
    sel <- sample(1:nrow(cmbs),amount)
    xp = cmbs[sel,1]
    yp = cmbs[sel,2]
    neworgdat[(lastind+1):(amount+lastind),'x']=xp
    neworgdat[(lastind+1):(amount+lastind),'y']=yp
    neworgdat[(lastind+1):(amount+lastind),'growth']=rep(growth, amount)
    neworgdat[(lastind+1):(amount+lastind),'type']=rep(type, amount)
    neworgdat[(lastind+1):(amount+lastind),'phenotype']=rep(ptype, amount)
    newoccmat <- as.matrix(newoccmat)
    for(i in 1:length(xp)){
      newoccmat[xp[i],yp[i]] = type
    }
    newoccmat <- Matrix(newoccmat,sparse=T)
  }else{
    neworgdat[(lastind+1):(amount+lastind),'x']=x
    neworgdat[(lastind+1):(amount+lastind),'y']=y
    neworgdat[(lastind+1):(amount+lastind),'growth']=rep(growth, amount)
    neworgdat[(lastind+1):(amount+lastind),'type']=rep(type, amount)
    neworgdat[(lastind+1):(amount+lastind),'phenotype']=rep(ptype, amount)
    newoccmat <- as.matrix(newoccmat)
    for(i in 1:length(x)){
      newoccmat[x[i],y[i]] = type
    }
    newoccmat <- Matrix(newoccmat,sparse=T)
  }
  eval.parent(substitute(object@occmat <- newoccmat))
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@specs <- newspecs))
  eval.parent(substitute(object@phenotypes[[spectype]] <- newphens))
  eval.parent(substitute(object@mediac <- union(object@mediac, specI@medium)))
})


setGeneric("addOrg2", function(object, specI, amount, x=NULL, y=NULL, growth=1, ...){standardGeneric("addOrg2")})
setMethod("addOrg2", "Arena", function(object, specI, amount, x=NULL, y=NULL, growth=1, ...){
  if(amount+sum(object@occmat) > object@n*object@m){
    stop("More individuals than space on the grid")
  }
  spectype <- specI@type
  newspecs <- object@specs
  newphens <- object@phenotypes[[spectype]]
  newspecs[[spectype]] <- specI
  type <- which(names(newspecs)==spectype)
  
  if(length(newphens)!=0){
    ptype <- as.integer(checkPhen(object, specI))
    newphens <- object@phenotypes[[spectype]]
  }else{
    newphens[[1]] <- getPhenotype(specI)
    ptype=as.integer(1)
  }
  l = addBacCpp(object@occmat, object@orgdat, amount, growth, type, ptype)
  newoccmat = l[["occmat"]]
  neworgdat = l[["orgdat"]]
  
  eval.parent(substitute(object@occmat <- newoccmat))
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@specs <- newspecs))
  eval.parent(substitute(object@phenotypes[[spectype]] <- newphens))
  eval.parent(substitute(object@mediac <- union(object@mediac, specI@medium)))
})


# Add all substances defined by exchange reactions of the available bacs

setGeneric("addSubs", function(object, smax=0, mediac=object@mediac){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax=0, mediac=object@mediac){
  if(sum(mediac %in% object@mediac)==length(mediac)){
    newmedia <- list()
    sapply(object@mediac, function(x, n, m){
      newmedia[[x]] <<- Substance(n, m, 0, name=x)
    }, n=object@n, m=object@m)
    for(i in 1:length(mediac)){
      newmedia[mediac[i]] <- Substance(object@n, object@m, smax=smax, name=mediac[i])
    }
    eval.parent(substitute(object@media <- newmedia))
  }else stop("Substance can't be produced or taken up by the organisms on the grid")
})

#function for changing the substances in the environment

setGeneric("changeSub", function(object, smax, mediac){standardGeneric("changeSub")})
setMethod("changeSub", "Arena", function(object, smax, mediac){
  if(length(sum(mediac %in% names(object@media)))==length(mediac)){
    for(i in 1:length(mediac)){
      eval.parent(substitute(object@media[mediac[i]] <- Substance(object@n, object@m, smax=smax, name=mediac[i])))
    }
  }else stop("Substance does not exist in medium")
})

#function for changing the Organisms in the environment

setGeneric("changeOrg", function(object, neworgdat){standardGeneric("changeOrg")})
setMethod("changeOrg", "Arena", function(object, neworgdat){
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@occmat <- Matrix(dat2mat(object), sparse=T)))
})

#function for checking if a phenotype is emergent

setGeneric("checkPhen", function(object, org, cutoff=1e-6){standardGeneric("checkPhen")})
setMethod("checkPhen", "Arena", function(object, org, cutoff=1e-6){
  ptype <- 0
  if(org@fbasol$obj>=cutoff){
    phenotypes <- object@phenotypes[[org@type]]
    phenspec <- getPhenotype(org, cutoff=0.1)
    if(length(phenspec) != 0){
      for(i in 1:length(phenotypes)){
        inlist <- intersect(names(phenotypes[[i]]),names(phenspec))
        if(sum(phenotypes[[i]][inlist]==phenspec[inlist])==length(inlist)){
          ptype=i
          break
        }
      }
      if(ptype==0){
        ptype = length(phenotypes)+1
        phenotypes[[ptype]] <- phenspec
        object2 <- object
        object2@phenotypes[[org@type]] <- phenotypes
        eval.parent(substitute(object <- object2)) #has to be like this, otherwise there is a problem with the slot name!
        #eval.parent(substitute(object@phenotypes[[org@type]] <- phenotypes))
      }
    }
  }
  return(ptype)
})

#main function for simulation of the whole arena

setGeneric("simulate", function(object, time){standardGeneric("simulate")})
setMethod("simulate", "Arena", function(object, time){
  switch(class(object),
         "Arena"={arena <- object; evaluation <- Eval(arena)},
         "Eval"={arena <- getArena(object); evaluation <- object},
         stop("Please supply an Arena object.")) 
  sublb <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
  for(j in seq_along(arena@media)){
    submat <- as.matrix(arena@media[[j]]@diffmat)
    sublb[,j] <- apply(arena@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
  }
  sublb <- cbind(as.matrix(arena@orgdat[,c(4,5)]),sublb)
  colnames(sublb) <- c('x','y',arena@mediac)
  rm("submat")
  for(i in 1:time){
    cat("iter:", i, "bacs:",nrow(arena@orgdat),"\n")
    print(system.time(for(j in 1:nrow(arena@orgdat)){
      org <- arena@specs[[arena@orgdat[j,'type']]]
      switch(class(org),
             "Bac"= {arena = simBac(org, arena, j, sublb)}, #the sublb matrix will be modified within this function
             "Human"= {arena = simHum(org, arena, j, sublb)}, #the sublb matrix will be modified within this function
             stop("Simulation function for Organism object not defined yet.")) 
    }))
    test <- is.na(arena@orgdat$growth)
    if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
    rm("test")
    if(!arena@stir){
      sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
      sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
      print(system.time(for(j in seq_along(arena@media)){ #get information from sublb matrix to media list
        submat <- as.matrix(arena@media[[j]]@diffmat)
        apply(sublb[,c('x','y',arena@media[[j]]@name)],1,function(x){submat[x[1],x[2]] <<- x[3]})
        switch(arena@media[[j]]@difunc,
               "cpp"={for(k in 1:arena@media[[j]]@difspeed){diffuseNaiveCpp(submat, donut=FALSE)}},
               "r"={for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
               stop("Simulation function for Organism object not defined yet.")) 
        arena@media[[j]]@diffmat <- Matrix(submat, sparse=T)
        sublb_tmp[,j] <- apply(arena@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
      }))
      sublb <- cbind(as.matrix(arena@orgdat[,c(4,5)]),sublb_tmp)
      colnames(sublb) <- c('x','y',arena@mediac)
      rm("sublb_tmp")
      rm("submat")
    }else{
      sublb <- stirEnv(arena, sublb)
    }
    addEval(evaluation, arena)
    if(sum(arena@occmat)==0){
      print("All organisms died!")
      break
    }
  }
  return(evaluation)
})

#function for stirring the complete evironment

setGeneric("stirEnv", function(object, sublb){standardGeneric("stirEnv")})
setMethod("stirEnv", "Arena", function(object, sublb){
  #stir all the bacteria
  neworgdat <- object@orgdat
  cmbs = expand.grid(1:object@n,1:object@m)
  if(nrow(neworgdat) > nrow(cmbs)){ #not so nice -> there is a problem
    selength <- nrow(cmbs)
  }else{
    selength <- nrow(neworgdat)
  }
  sel <- sample(1:nrow(cmbs),selength)
  neworgdat[,'x'] <- cmbs[sel,1]
  neworgdat[,'y'] <- cmbs[sel,2]
  newoccmat <- matrix(0,object@n,object@m)
  for(i in 1:nrow(neworgdat)){
    newoccmat[neworgdat[i,'x'],neworgdat[i,'y']] = neworgdat[i,'type']
  }
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@occmat <- Matrix(newoccmat,sparse=T)))
  #stir all the substrates + modify substrates
  sublb_tmp <- matrix(0,nrow=nrow(object@orgdat),ncol=(length(object@mediac)))
  sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
  for(j in seq_along(object@media)){ #get information from sublb matrix to media list
    sval <- sum(sublb[,object@media[[j]]@name])/nrow(sublb)
    submat <- matrix(sval,object@n,object@m)
    eval.parent(substitute(object@media[[j]]@diffmat <- Matrix(submat, sparse=T)))
    sublb_tmp[,j] <- apply(object@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
  }
  sublb <- cbind(as.matrix(object@orgdat[,c(4,5)]),sublb_tmp)
  colnames(sublb) <- c('x','y',object@mediac)
  return(sublb)
})

#function for converting orgdat in a matrix

setGeneric("dat2mat", function(object){standardGeneric("dat2mat")})
setMethod("dat2mat", "Arena", function(object){
  newoccmat <- matrix(0,object@n,object@m)
  for(i in 1:nrow(object@orgdat)){
    newoccmat[object@orgdat[i,'x'],object@orgdat[i,'y']] = object@orgdat[i,'type']
  }
  return(newoccmat)
})

#show function for class Arena

setMethod(show, "Arena", function(object){
  print(paste('Arena of size ',object@n,'x',object@m,' with ',sum(object@occmat),
              ' organisms of ',length(object@specs),' species.',sep=''))
})

# Eval is a subclass of Arena containing function to reduce the size of simulations and evalution of results

########################################################################################################
###################################### EVAL CLASS ######################################################
########################################################################################################

setClass("Eval",
         contains="Arena",
         representation(
           medlist="list", # list of simulation results with medium concentration
           simlist="list", # list of simulation results with organism features
           subchange="numeric" # vector of all substrates with numbers indicating the degree of change
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Eval <- function(arena){
  subc = rep(0, length(arena@mediac))
  names(subc) <- arena@mediac
  new("Eval", n=arena@n, m=arena@m, tstep=arena@tstep, specs=arena@specs, mediac=arena@mediac, occmat=Matrix(), subchange=subc,
      phenotypes=arena@phenotypes, media=arena@media, orgdat=arena@orgdat, medlist=list(), simlist=list(), stir=arena@stir)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("medlist", function(object){standardGeneric("medlist")})
setMethod("medlist", "Eval", function(object){return(object@medlist)})
setGeneric("simlist", function(object){standardGeneric("simlist")})
setMethod("simlist", "Eval", function(object){return(object@simlist)})
setGeneric("subchange", function(object){standardGeneric("subchange")})
setMethod("subchange", "Eval", function(object){return(object@subchange)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the Eval object -> adding or replacing simulations

setGeneric("addEval", function(object, arena, replace=F){standardGeneric("addEval")})
setMethod("addEval", "Eval", function(object, arena, replace=F){
  if(!replace){
    subch = rep(0, length(arena@mediac))
    if(length(object@medlist)!=0){
      names(subch) <- arena@mediac
      sapply(names(subch), function(x, oldmed, newmed){
        subch[x] <<- subch[x]+sum(abs(oldmed[[x]]-as.vector(newmed[[x]]@diffmat)))
      },oldmed=extractMed(object), newmed=arena@media)
      eval.parent(substitute(object@subchange <- object@subchange + subch))
    }
    if(sum(subch)!=0){
      eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x, subc){
        if(subc[x@name]!=0){
          return(as.vector(x@diffmat))
        }else{return(vector())}
      }, subc=subch)))
    }else{
      eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x){
        return(as.vector(x@diffmat))
      })))
    }
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
  },meds=extractMed(object,time), n=object@n, m=object@m)
  
  occdat <- object@simlist[[time]]
  newoccmat <- matrix(0, nrow=object@n, ncol=object@m)
  apply(occdat, 1, function(x){
    newoccmat[x[4],x[5]] <<- x[2]
  })
  
  arena <- Arena(n=object@n, m=object@m, tstep=object@tstep, specs=object@specs, mediac=object@mediac,
                 phenotypes=object@phenotypes , media=newmedia, orgdat=occdat, occmat=Matrix(newoccmat,sparse=T), stir=object@stir)
  return(arena)
})

#function for re-extracting medlist to original object

setGeneric("extractMed", function(object, ind=length(object@medlist)){standardGeneric("extractMed")})
setMethod("extractMed", "Eval", function(object, ind=length(object@medlist)){
  medl <- object@medlist
  medlind <- medl[[ind]]
  for(i in 1:length(object@mediac)){
    if(length(medl[[ind]][[i]])==0){
      j <- ind
      while(length(medl[[j]][[i]])==0){j <- j-1}
      medlind[[i]] <- medl[[j]][[i]]
    }
  }
  return(medlind)
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
    meds <- extractMed(object, i)
    if(length(plot_items)==1){
      if(length(inds)!=0){
        for(j in 1:length(inds)){
          if(retdata){
            retlist[[subnam[inds[j]]]][[j]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
          }
          image(matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
        }
      }
    }else if(length(plot_items)<=6){
      par(mfrow=c(2,ceiling(length(plot_items)/2)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[j]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
      }
    }else{
      par(mfrow=c(3,ceiling(length(plot_items)/3)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[j]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(meds[[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]])
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
    subdat <- extractMed(object, i)
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

setMethod(show, signature(object="Eval"), function(object){
  print(paste('Evaluation results of ',length(object@medlist),' simulation steps.',sep=''))
})

