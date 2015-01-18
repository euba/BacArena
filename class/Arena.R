

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
           n="integer",  # grid size
           m="integer"  # grid size
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Arena <- function(n,m,tstep=1,orgdat=data.frame(growth=numeric(0),type=integer(0),phenotype=integer(0),x=integer(0),y=integer(0)),
                  specs=list(),media=list(),mediac=character(),phenotypes=list(),occmat=Matrix(0L,nrow=n,ncol=m,sparse=T)){
  new("Arena", n=as.integer(n), m=as.integer(m), tstep=tstep, orgdat=orgdat, specs=specs,
      media=media, mediac=mediac, phenotypes=phenotypes, occmat=occmat)
}

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

# Add Individuals to the arena

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
    for(i in 1:amount){
      neworgdat[lastind+i,'growth']=growth
      neworgdat[lastind+i,'type']=type
      neworgdat[lastind+i,'phenotype']=ptype
      x <- sample(1:n, 1)
      y <- sample(1:m, 1)
      while(newoccmat[x,y] != 0){
        x <- sample(1:n, 1)
        y <- sample(1:m, 1)
      }
      neworgdat[lastind+i,'x']=x
      neworgdat[lastind+i,'y']=y
      newoccmat[x,y] = type
    }
  }else{
    for(i in 1:amount){
      neworgdat[lastind+i,'growth']=growth
      neworgdat[lastind+i,'type']=type
      neworgdat[lastind+i,'phenotype']=ptype
      neworgdat[lastind+i,'x']=x[i]
      neworgdat[lastind+i,'y']=y[i]
      newoccmat[x[i],y[i]] = type
    }
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

setGeneric("addSubs", function(object, smax, mediac=object@mediac){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax, mediac=object@mediac){
  newmedia = object@media
  n = object@n
  m = object@m
  mediac = setdiff(mediac, names(arena@media))
  if(length(mediac)==0){
    print("No novel compound can be included to the medium.")
  }else{
    sapply(mediac, function(x, smax, n, m){
      newmedia[[x]] <<- Substance(n, m, smax, name=x)
    }, smax=smax, n=n, m=m)
    
    eval.parent(substitute(object@media <- c(object@media,newmedia)))
  }
  newspecs <- lapply(object@specs, function(x,mediac){
    x@medium = intersect(x@medium, mediac)
    return(x)
  }, mediac=mediac)
  eval.parent(substitute(object@specs <- newspecs))
})


setGeneric("changeSub", function(object, subname, value){standardGeneric("changeSub")})
setMethod("changeSub", "Arena", function(object, subname, value){
  if (subname %in% names(object@media)) eval.parent(substitute(object@media[subname] <- Substance(object@n, object@m, smax=value, name=subname)))
  #return(Substance(object@n, object@m, smax=value, name=subname))
  else stop("Substance does not exist in medium")
})

#function for getting a vector of media concentrations for a specific position

#setGeneric("getmed", function(object, diffmat, sublb, subname){standardGeneric("getmed")})
#setMethod("getmed", "Arena", function(object, diffmat, sublb, subname){
#  sublb[,subname] <- diag(diffmat[object@orgdat$x,object@orgdat$y])
#  return(sublb)
#})

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
             "Bac"= {arena = simBac(org, arena, j, sublb)},
             stop("Simulation function for Organism object not defined yet.")) 
    }))
    test <- is.na(arena@orgdat$growth)
    if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
    rm("test")
    
    sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
    sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
    print(system.time(for(j in seq_along(arena@media)){
      submat <- as.matrix(arena@media[[j]]@diffmat)
      apply(sublb[,c('x','y',arena@media[[j]]@name)],1,function(x){submat[x[1],x[2]] <<- x[3]})
      switch(arena@media[[j]]@difunc,
             "cpp"={for(k in 1:arena@media[[j]]@difspeed){diffuseNaiveCpp(submat, donut=FALSE)}},
             "r"={for(k in 1:arena@media[[j]]@difspeed){diffuseNaiveR(arena@media[[j]])}},
             stop("Simulation function for Organism object not defined yet.")) 
      arena@media[[j]]@diffmat <- Matrix(submat, sparse=T)
      sublb_tmp[,j] <- apply(arena@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
    }))
    sublb <- cbind(as.matrix(arena@orgdat[,c(4,5)]),sublb_tmp)
    colnames(sublb) <- c('x','y',arena@mediac)
    rm("sublb_tmp")
    rm("submat")
    
    addEval(evaluation, arena)
    if(sum(arena@occmat)==0){
      print("All organisms died!")
      break
    }
  }
  return(evaluation)
})

#show function for class Arena

removeMethod(show, "Arena")
setMethod(show, "Arena", function(object){
  print(paste('Arena of size ',object@n,'x',object@m,' with ',sum(object@occmat),
              ' organisms of ',length(object@specs),' species.',sep=''))
})

