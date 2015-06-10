########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

#' Structure of the S4 class "Arena"
#' 
#' Structure of the S4 class \code{Arena} to represent the environment in which Organisms and Substances interact.
#'
#' @slot orgdat A data frame collecting information about the accumulated growth, type, phenotype, x and y position for each individual in the environment.
#' @slot specs A list of organism types and their associated parameters.
#' @slot media A list of objects of class \code{\link{Substance-class}} for each compound in the environment.
#' @slot phenotypes A list of unique phenotypes (metabolites consumed and produced), which occurred in the environment.
#' @slot mediac A character vector containing the names of all substances in the environment.
#' @slot occmat A sparse matrix showing which cells in the environment are occupied by individuals.
#' @slot tstep A number giving the time (in h) per iteration.
#' @slot stir A boolean variable indicating if environment should be stirred.
#' @slot n A number giving the horizontal size of the environment.
#' @slot m A number giving the vertical size of the environment.
setClass("Arena",
         representation(
           orgdat="data.frame",
           specs="list",
           media="list",
           phenotypes="list",
           mediac="character",
           occmat="Matrix",
           tstep="numeric",
           stir="logical",
           n="integer",
           m="integer"
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

#' @title Add individuals to the environment
#'
#' @description The generic function \code{addOrg} adds individuals to the environment.
#'
#' @param object An object of class Arena.
#' @param specI An object of class Organism.
#' @param amount A numeric number giving the number of individuals to add.
#' @param x A numeric vector giving the x positions of individuals on the grid.
#' @param y A numeric vector giving the y positions of individuals on the grid.
#' @param growth A numeric vector giving the starting biomass of the individuals.
#' @details The arguments \code{x} and \code{y} should be in the same length as the number of organisms added (given by the argument \code{amount}).
#' @seealso \code{\link{Arena-class}} and \code{\link{Bac-class}} 
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' }
setGeneric("addOrg", function(object, specI, amount, x=NULL, y=NULL, growth=1){standardGeneric("addOrg")})
setMethod("addOrg", "Arena", function(object, specI, amount, x=NULL, y=NULL, growth=1){
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
  if(sum(duplicated(paste(neworgdat$x,neworgdat$y,sep="_")))!=0){
    stop("You have multiple individuals in the same position! Make sure that your x an y positions are unique")
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

#' @title Add substances to the environment
#'
#' @description The generic function \code{addSubs} adds specific substances to the environment.
#'
#' @param object An object of class Arena.
#' @param smax A number indicating the maximum substance concentration per grid cell.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param smax A number indicating the maximum substance concentration per grid cell.
#' @slot difunc A character vector ("pde","cpp" or "r") describing the function for diffusion.
#' @slot difspeed A number indicating the diffusion speed (given by number of cells per iteration).
#' @details If nothing but \code{object} is given, then all possible substrates are initilized with a concentration of 0. Afterwards, \code{\link{changeSub} can be used to modify the concentrations of specific substances.} 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addSubs(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) #add substances glucose, oxygen and phosphate
#' }
setGeneric("addSubs", function(object, smax=0, mediac=object@mediac, difunc="pde", difspeed=1){standardGeneric("addSubs")})
setMethod("addSubs", "Arena", function(object, smax=0, mediac=object@mediac, difunc="pde", difspeed=1){
  if(sum(mediac %in% object@mediac)==length(mediac)){
    newmedia <- list()
    sapply(object@mediac, function(x, n, m){
      newmedia[[x]] <<- Substance(n, m, 0, name=x)
    }, n=object@n, m=object@m)
    for(i in 1:length(mediac)){
      newmedia[mediac[i]] <- Substance(object@n, object@m, smax=smax, name=mediac[i], difunc=difunc, difspeed=difspeed)
    }
    eval.parent(substitute(object@media <- newmedia))
  }else stop("Substance can't be produced or taken up by the organisms on the grid")
})

#' @title Change substances in the environment
#'
#' @description The generic function \code{changeSub} changes specific substances in the environment.
#'
#' @param object An object of class Arena.
#' @param smax A number indicating the maximum substance concentration per grid cell.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @details If nothing but \code{object} is given, then all possible substrates are initilized with a concentration of 0. Afterwards, \code{\link{changeSub}} can be used to modify the concentrations of specific substances.
#' @seealso \code{\link{Arena-class}} and \code{\link{addSubs}} 
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addSubs(arena) #add all substances with no concentrations.
#' changeSub(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) #add substances glucose, oxygen and phosphate
#' }
setGeneric("changeSub", function(object, smax, mediac){standardGeneric("changeSub")})
setMethod("changeSub", "Arena", function(object, smax, mediac){
  if(sum(mediac %in% names(object@media))==length(mediac)){
    for(i in 1:length(mediac)){
      eval.parent(substitute(object@media[mediac[i]] <- Substance(object@n, object@m, smax=smax, name=mediac[i],
                                                                  difunc=object@media[[mediac[i]]]@difunc,
                                                                  difspeed=object@media[[mediac[i]]]@difspeed)))
    }
  }else stop("Substance does not exist in medium")
})

#' @title Change substance concentration patterns in the environment
#'
#' @description The generic function \code{changeDiff} changes specific substance concentration patterns in the environment.
#'
#' @param object An object of class Arena.
#' @param newdiffmat A matrix giving the new gradient matrix of the specific substances in the environment.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @details This function can be used to add gradients of specific substances in the environment. The default conditions in \code{changeSubs} assumes an equal concentration in every grid cell of the environment. 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSubs}} 
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,30) #add all substances with no concentrations.
#' gradient <- matrix(1:200,20,20)
#' changeDiff(arena,gradient,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) #add substances glucose, oxygen and phosphate
#' }
setGeneric("changeDiff", function(object, newdiffmat, mediac){standardGeneric("changeDiff")})
setMethod("changeDiff", "Arena", function(object, newdiffmat, mediac){
  if(nrow(newdiffmat)==object@n && ncol(newdiffmat)==object@m){
    for(i in 1:length(mediac)){
      eval.parent(substitute(object@media[[mediac[i]]]@diffmat <- Matrix(newdiffmat, sparse=T)))
    }
  }else stop("Given matrix is not compatible in dimensions with the environment.")
})

#' @title Change substance concentration patterns in the environment according to a gradient
#'
#' @description The generic function \code{createGradient} changes specific substance concentration patterns in the environment.
#'
#' @param object An object of class Arena.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param position A character vector giving the position (top, bottom, right and left) of the gradient.
#' @param smax A number giving the maximum concentration of the substance.
#' @param steep A number between 0 and 1 giving the steepness of the gradient (concentration relative to the arena size).
#' @details This function can be used to add gradients of specific substances in the environment. 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSubs}} 
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,30) #add all substances with no concentrations.
#' createGradient(arena,smax=50,mediac=c("EX_glc(e)","EX_o2(e)","EX_pi(e)"),
#'              position='top',steep=0.5)
#' }
setGeneric("createGradient", function(object, mediac, position, smax, steep){standardGeneric("createGradient")})
setMethod("createGradient", "Arena", function(object, mediac, position, smax, steep){
  if(steep<=0 || steep>=1){stop("Steepness must be in between 0 and 1.")}
  mediac = intersect(mediac,object@mediac)
  newdiffmat <- matrix(0,nrow=object@n,ncol=object@m)
  gradn = floor(object@n*steep)
  gradm = floor(object@m*steep)
  switch(position,
         'top'={for(i in 1:object@m){newdiffmat[0:gradm+1,i]=seq(smax,0,length.out=gradm+1)}},
         'bottom'={for(i in 1:object@m){newdiffmat[gradm:object@m,i]=seq(0,smax,length.out=gradm+1)}},
         'right'={for(i in 1:object@n){newdiffmat[i,gradn:object@n]=seq(0,smax,length.out=gradn+1)}},
         'left'={for(i in 1:object@n){newdiffmat[i,0:gradn+1]=seq(smax,0,length.out=gradn+1)}},
         stop("Positions must be top, bottom, right and left."))
  for(i in 1:length(mediac)){
    eval.parent(substitute(object@media[[mediac[i]]]@diffmat <- Matrix(newdiffmat, sparse=T)))
  }
})

#' @title Change organisms in the environment
#'
#' @description The generic function \code{changeOrg} changes organisms in the environment.
#'
#' @param object An object of class Arena.
#' @param neworgdat A data frame with new information about the accumulated growth, type, phenotype, x and y position for each individual in the environment.
#' @details The argument \code{neworgdat} contains the same information as the \code{orgdat} slot of \code{\link{Arena-class}}. The \code{orgdat} slot of an \code{Arena} object can be used to create \code{neworgdat}.
#' @seealso \code{\link{Arena-class}} and \code{\link{addOrg}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' neworgdat <- orgdat(arena) #get the current orgdat
#' neworgdat <- neworgdat[-1,] #remove the first individual
#' changeOrg(arena,neworgdat)
#' }
setGeneric("changeOrg", function(object, neworgdat){standardGeneric("changeOrg")})
setMethod("changeOrg", "Arena", function(object, neworgdat){
  eval.parent(substitute(object@orgdat <- neworgdat))
  eval.parent(substitute(object@occmat <- Matrix(dat2mat(object), sparse=T)))
})

#' @title Function for checking phenotypes in the environment
#'
#' @description The generic function \code{checkPhen} checks and adds the phenotypes of organisms in the environment.
#'
#' @param object An object of class Arena.
#' @param org An object of class Organism.
#' @param cutoff A number giving the cutoff for values of the objective function and fluxes of exchange reactions.
#' @return Returns a number indicating the number of the phenotype in the phenotype list.
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages. Uptake of substances are indicated by a negative and production of substances by a positive number.
#' @seealso \code{\link{Arena-class}} and \code{\link{getPhenotype}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' checkPhen(arena,bac) #returns 1 as the index of the current phenotype in the list.
#' }
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
      }
    }
  }
  return(ptype)
})

#' @title Main function for simulating all processes in the environment
#'
#' @description The generic function \code{simEnv} for a simple simulation of the environment.
#'
#' @param object An object of class Arena or Eval.
#' @param time A number giving the number of iterations to perform for the simulation
#' @return Returns an object of class \code{Eval} which can be used for subsequent analysis steps.
#' @details The returned object itself can be used for a subsequent simulation, due to the inheritance between \code{Eval} and \code{Arena}.
#' @seealso \code{\link{Arena-class}} and \code{\link{Eval-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' }
setGeneric("simEnv", function(object, time){standardGeneric("simEnv")})
setMethod("simEnv", "Arena", function(object, time){
  switch(class(object),
         "Arena"={arena <- object; evaluation <- Eval(arena)},
         "Eval"={arena <- getArena(object); evaluation <- object},
         stop("Please supply an Arena object.")) 
  sublb <- getSublb(arena)
  for(i in 1:time){
    cat("iter:", i, "Organisms:",nrow(arena@orgdat),"\n")
    for(j in 1:nrow(arena@orgdat)){
      org <- arena@specs[[arena@orgdat[j,'type']]]
      switch(class(org),
             "Bac"= {arena = simBac(org, arena, j, sublb)}, #the sublb matrix will be modified within this function
             "Human"= {arena = simHum(org, arena, j, sublb)}, #the sublb matrix will be modified within this function
             stop("Simulation function for Organism object not defined yet.")) 
    }
    test <- is.na(arena@orgdat$growth)
    if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
    rm("test")
    if(!arena@stir){
      sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
      sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
      for(j in seq_along(arena@media)){ #get information from sublb matrix to media list
        submat <- as.matrix(arena@media[[j]]@diffmat)
        apply(sublb[,c('x','y',arena@media[[j]]@name)],1,function(x){submat[x[1],x[2]] <<- x[3]})
        switch(arena@media[[j]]@difunc,
               "pde"={diffuseGrajdeanuCpp(submat, donut=FALSE, mu=arena@media[[j]]@difspeed)},
               "cpp"={for(k in 1:arena@media[[j]]@difspeed){diffuseNaiveCpp(submat, donut=FALSE)}},
               "r"={for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
               stop("Simulation function for Organism object not defined yet.")) 
        arena@media[[j]]@diffmat <- Matrix(submat, sparse=T)
        sublb_tmp[,j] <- apply(arena@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
      }
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

#' @title Function for calculated the substrate concentration for every organism
#'
#' @description The generic function \code{getSublb} calculates the substrate concentration for every individual in the environment based on their x and y position.
#'
#' @param object An object of class Arena.
#' @return Returns the substrate concentration for every individual in the environment with substrates as well as x and y positions as columns and rows for each organism.
#' @seealso \code{\link{Arena-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' sublb <- getSublb(arena)
#' }
setGeneric("getSublb", function(object){standardGeneric("getSublb")})
setMethod("getSublb", "Arena", function(object){
  sublb <- matrix(0,nrow=nrow(object@orgdat),ncol=(length(object@mediac)))
  for(j in seq_along(object@media)){
    submat <- as.matrix(object@media[[j]]@diffmat)
    sublb[,j] <- apply(object@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
  }
  sublb <- cbind(as.matrix(object@orgdat[,c(4,5)]),sublb)
  colnames(sublb) <- c('x','y',object@mediac)
  return(sublb)
})

#' @title Function for stirring/mixing the complete evironment
#'
#' @description The generic function \code{stirEnv} simulates the event of mixing all substrates and organisms in the environment.
#'
#' @param object An object of class Arena.
#' @param sublb A matrix with the substrate concentration for every individual in the environment based on their x and y position.
#' @return Returns the substrate concentration for every individual in the environment with substrates as well as x and y positions as columns and rows for each organism.
#' @details The stirring is implemented as a random permutation of organism positions and the equalization of of all substrate concentrations.
#' @seealso \code{\link{Arena-class}} and \code{\link{getSublb}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' sublb <- getSublb(arena)
#' stirEnv(arena,sublb)
#' }
setGeneric("stirEnv", function(object, sublb){standardGeneric("stirEnv")})
setMethod("stirEnv", "Arena", function(object, sublb){
  #stir all the organism except the ones which are not moving
  neworgdat <- object@orgdat
  cmbs = expand.grid(1:object@n,1:object@m)
  sit <- which(unlist(lapply(object@specs,function(x)(return(x@speed))))==0)
  if(length(sit) != 0){
    siti <- which(neworgdat$type %in% sit)
    sitorgdat <- neworgdat[siti,]
    neworgdat <- neworgdat[-siti,]
    rownames(cmbs) <- paste(cmbs$Var1,cmbs$Var2,sep="_")
    cmbs <- cmbs[setdiff(rownames(cmbs),c(paste(sitorgdat$x,sitorgdat$y,sep="_"))),]
  }
  if(nrow(neworgdat) > nrow(cmbs)){ #not so nice -> there is a problem
    selength <- nrow(cmbs)
  }else{
    selength <- nrow(neworgdat)
  }
  sel <- sample(1:nrow(cmbs),selength)
  neworgdat[,'x'] <- cmbs[sel,1]
  neworgdat[,'y'] <- cmbs[sel,2]
  newoccmat <- matrix(0,object@n,object@m)
  if(length(sit) != 0){neworgdat <- rbind(neworgdat,sitorgdat)}
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

#' @title Function for transforming the organism data frame to a presence/absence matrix of organisms
#'
#' @description The generic function \code{dat2mat} simulates the event of mixing all substrates and organisms in the environment.
#'
#' @param object An object of class Arena.
#' @return Returns the presence/absence matrix of organisms on the grid based on the \code{orgdat} slot of the \code{Arena} class.
#' @seealso \code{\link{Arena-class}} and \code{\link{getSublb}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' occmat <- dat2mat(arena)
#' image(occmat)
#' }
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
  print(paste('Arena of size ',object@n,'x',object@m,' with ',sum(ifelse(as.matrix(object@occmat)==0,0,1)),
              ' organisms of ',length(object@specs),' species.',sep=''))
})

# Eval is a subclass of Arena containing function to reduce the size of simulations and evalution of results

########################################################################################################
###################################### EVAL CLASS ######################################################
########################################################################################################

#' Structure of the S4 class "Eval"
#' 
#' Structure of the S4 class \code{Eval} inheriting from class \code{\link{Arena-class}} for the analysis of simulations.
#'
#' @slot medlist A list of compressed medium concentrations (only changes of concentrations are stored) per time step.
#' @slot simlist A list of the organism features per time step.
#' @slot subchange A vector of all substrates with numbers indicating the degree of change in the overall simulation.
setClass("Eval",
         contains="Arena",
         representation(
           medlist="list",
           simlist="list",
           subchange="numeric"
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

#' @title Function for adding a simulation step
#'
#' @description The generic function \code{addEval} adds results of a simulation step to an \code{Eval} object.
#'
#' @param object An object of class Eval.
#' @param arena An object of class Arena.
#' @param replace A boolean variable indicating if the last simulation step should be replaced by the new simulation step \code{arena}.
#' @details The function \code{addEval} can be used in iterations to manipulate an \code{Arena} object and store the results in an \code{Eval} object.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' addEval(eval,arena)
#' }
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

#' @title Function for re-constructing an Arena object from a simulation step
#'
#' @description The generic function \code{getArena} re-constructs an \code{Arena} object from a simulation step within an \code{Eval} object.
#'
#' @param object An object of class Eval.
#' @param time A number giving the simulation step of interest.
#' @return Returns an object of class \code{Arena} containing the organisms and substance conditions in simulation step \code{time}.
#' @details The function \code{addEval} can be used to manipulate an \code{Arena} object from a simulation step to modify the subsequent simulation steps.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' arena5 <- getArena(eval,5)
#' }
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

#' @title Function for re-constructing a medium concentrations from simulations
#'
#' @description The generic function \code{extractMed} re-constructs a list of vectors of medium concentrations from a simulation step in an \code{Eval} object.
#'
#' @param object An object of class Eval.
#' @param ind A number giving the simulation step of interest.
#' @return Returns a list containing concentration vectors of all medium substances.
#' @details Medium concentrations in slot \code{medlist} of an object of class \code{Eval} store only the changes of concentrations in the simulation process. The function \code{extractMed} reconstructs the original and uncompressed version of medium concentrations.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' med5 <- extractMed(eval,5)
#' }
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

#' @title Function for plotting spatial and temporal change of populations and/or concentrations
#'
#' @description The generic function \code{evalArena} plots heatmaps from the simulation steps in an \code{Eval} object.
#'
#' @param object An object of class Eval.
#' @param plot_items A character vector giving the items, which should be plotted.
#' @param phencol A boolean variable indicating if the phenotypes of the organisms in the environment should be integrated as different colors in the population plot.
#' @param retdata A boolean variable indicating if the data used to generate the plots should be returned.
#' @param sims A numeric vector giving the simulation steps which should be plotted.
#' @return Returns several plots of the chosen plot items. Optional the data to generate the original plots can be returned.
#' @details If \code{phencol} is \code{TRUE} then different phenotypes of the same organism are visualized by varying colors, otherwise different organism types are represented by varying colors. The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' evalArena(eval)
#' ## if animation package is installed a movie of the simulation can be stored:
#' library(animation)
#' saveVideo({evalArena(eval)},video.name="Ecoli_sim.mp4")
#' }
setGeneric("evalArena", function(object, plot_items='population', phencol=F, retdata=F, sims=1:length(object@simlist)){standardGeneric("evalArena")})
setMethod("evalArena", "Eval", function(object, plot_items='population', phencol=F, retdata=F, sims=1:length(object@simlist)){
  old.par <- par(no.readonly = TRUE)
  if(retdata){
    retlist = list()
    for(i in 1:length(plot_items)){
      retlist[[i]] = list()
    }
    names(retlist) = plot_items
  }
  for(i in sims){
    subnam <- names(object@medlist[[i]])
    inds <- which(subnam %in% plot_items)
    meds <- extractMed(object, i)
    if(length(plot_items)==1){
      if(length(inds)!=0){
        for(j in 1:length(inds)){
          if(retdata){
            retlist[[subnam[inds[j]]]][[j]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
          }
          image(matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]],
                zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
        }
      }
    }else if(length(plot_items)<=6){
      par(mfrow=c(2,ceiling(length(plot_items)/2)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[j]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]],
              zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
      }
    }else{
      par(mfrow=c(3,ceiling(length(plot_items)/3)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[j]] = matrix(meds[[subnam[inds[j]]]],nrow=object@n,ncol=object@m)
        }
        image(matrix(meds[[inds[j]]],nrow=object@n,ncol=object@m),axes=F,main=subnam[inds[j]],
              zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
      }
    }
    if(plot_items[1]=='population'){
      if(retdata){
        retlist[['population']][[j]] = object@simlist[[i]]
      }
      if(phencol){
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=object@simlist[[i]]$type-1,axes=FALSE,cex=1,main='Population', col=object@simlist[[i]]$phenotype+1)
      }else{
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=object@simlist[[i]]$type-1,axes=FALSE,cex=1,main='Population', col=object@simlist[[i]]$type)
      }
    }
  }
  if(retdata){
    return(retlist)
  }
  par(old.par)
})

#' @title Function for plotting the overall change as curves
#'
#' @description The generic function \code{plotCurves} plots the growth curves and concentration changes of substances from simulation steps in an \code{Eval} object.
#'
#' @param object An object of class Eval.
#' @param medplot A character vector giving the name of substances which should be plotted.
#' @param retdata A boolean variable indicating if the data used to generate the plots should be returned.
#' @param remove A boolean variable indicating if substances, which don't change in their concentration should be removed from the plot.
#' @param legend Boolean variable indicating if legend should be plotted
#' @return Returns two graphs in one plot: the growth curves and the curves of concentration changes. Optional the data to generate the original plots can be returned.
#' @details The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' plotCurves(eval)
#' }
setGeneric("plotCurves", function(object, medplot=object@mediac, retdata=F, remove=F, legend=F){standardGeneric("plotCurves")})
setMethod("plotCurves", "Eval", function(object, medplot=object@mediac, retdata=F, remove=F, legend=F){
  old.par <- par(no.readonly = TRUE)
  growths <- matrix(0, nrow=length(object@specs), ncol=length(object@simlist))
  rownames(growths) = names(object@specs)
  subs <- matrix(0, nrow=length(medplot), ncol=length(object@simlist))
  rownames(subs) = medplot
  for(i in 1:length(object@simlist)){
    simdat <- object@simlist[[i]]
    #count <- table(simdat[,'type'])
    #for(j in 1:length(count)){
    #  growths[j,i] <- count[j]
    #}
    abun <- table(simdat$type)
    for(j in 1:length(abun)){
      growths[as.numeric(names(abun[j])),i] = abun[j]
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
    lines(times, growths[i,], col=i, type='b', pch=i-1)
  }
  # and length(object@specs) > 1
  if(legend){legend('topleft',legend=rownames(growths),col=1:nrow(growths), cex=ifelse(length(object@specs)==1,1,0.4/log10(nrow(growths)+1)),pch=(0:nrow(growths)-1),lwd=1, bty="n")}
  #if(legend){legend('topleft',legend=rownames(growths),col=1:nrow(growths), cex=0.4/log10(nrow(growths)+1),pch=(0:nrow(growths)-1),lwd=1, bty="n")}
  plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(subs)),
       type='n', xlab='time in h', ylab='concentration in mmol per gridcell',
       main='Substance concentrations')
  for(i in 1:nrow(subs)){
    lines(times, subs[i,], col=i)
  }
  if(legend){legend('right',legend=rownames(subs),col=1:nrow(subs),cex=0.4/log10(nrow(subs)+1),lwd=1)}
  if(retdata){
    return(list('Population'=growths,'Substances'=subs))
  }
  par(old.par)
})

#' @title Function for getting a matrix of phenotypes from the dataset
#'
#' @description The generic function \code{getPhenoMat} reconstructs a matrix with the usage of exchange reactions of the different organisms in the environment.
#'
#' @param object An object of class Eval.
#' @return Returns a matrix with different phenotypes of the organism as rows and all possible exchange reactions as columns. A value of 1 means secretion, 2 means uptake and 0 means no usage of the substance of interest.
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenotype}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' phenmat <- getPhenoMat(eval)
#' }
setGeneric("getPhenoMat", function(object, step="total"){standardGeneric("getPhenoMat")})
setMethod("getPhenoMat", "Eval", function(object, step="total"){
  if(step=="total"){
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
  }else{
    typestep = object@simlist[[step]]$type
    phenstep = object@simlist[[step]]$phenotype
    typenam = names(object@phenotypes)
    phens = list()
    for(i in levels(as.factor(typestep))){
      tphen = object@phenotypes[[as.numeric(i)]]
      phens[[typenam[as.numeric(i)]]] = tphen[as.numeric(levels(as.factor(phenstep[which(typestep==as.numeric(i))])))]
    }
    
    numphens <- unlist(lapply(phens,function(x){return(length(x))}))
    phentypes <- vector()
    for(i in seq_along(numphens)){
      phentypes <- c(phentypes, rep(names(numphens)[i],numphens[i]))
    }
    phentypes <- as.factor(phentypes)
    phenmat <- matrix(0, nrow=length(phentypes), ncol=length(object@mediac))
    colnames(phenmat) <- object@mediac
    rownames(phenmat) <- phentypes
    pind <- 0
    for(i in 1:length(phens)){
      for(j in 1:numphens[i]){
        pvec <- phens[[i]][[j]]
        pind <- pind + 1
        phenmat[pind, names(pvec)] <- pvec
      }
    }
    phenmat <- ifelse(phenmat==-1,2,phenmat)
    return(phenmat)
  }
})

#' @title Function for mining/analyzing phenotypes which occured on the arena
#'
#' @description The generic function \code{minePheno} mines the similarity and differences of phenotypes reconstructed by \code{getPhenoMat} for each simulation step in an \code{Eval} object.
#'
#' @param object An object of class Eval.
#' @param plot_type A character vector giving the plot which should be returned (either "pca" for a principle coordinate analysis or "hclust" for hierarchical clustering).
#' @param legend Boolean variable indicating if legend should be plotted
#' @return Returns a plot for each simulation step representing the similarity of phenotypes of organisms within the environment. 
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenoMat}}
#' @examples
#' \dontrun{
#' ecore <- model #get Escherichia coli core metabolic model
#' bac <- Bac(ecore,deathrate=0.05,duplirate=0.5,
#'            growthlimit=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(20,20) #initialize the environment
#' addOrg(arena,bac,amount=10) #add 10 organisms
#' addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,10)
#' minePheno(eval)
#' }
setGeneric("minePheno", function(object, plot_type="pca", legend=F, step="total"){standardGeneric("minePheno")})
setMethod("minePheno", "Eval", function(object, plot_type="pca", legend=F, step="total"){
  phenmat <- getPhenoMat(object, step)
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
      points(mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),1]),mean(phenpca$x[which(rownames(phenpca$x)==typ[i]),2]),col=i,pch=i-1,cex=1.5)
    }
    if(legend){legend('topright',legend=names(object@specs),col=1:length(object@specs),pch=(0:length(object@specs)-1),cex=0.9,lwd=4)}
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
