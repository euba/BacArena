globalVariables(c("diffuseNaiveCpp","diffuseSteveCpp"))

########################################################################################################
###################################### Arena CLASS ################################################
########################################################################################################

#' Structure of the S4 class "Arena"
#' 
#' Structure of the S4 class \code{Arena} to represent the environment in which Organisms and Substances interact.
#' @import Matrix Rcpp
#' @export Arena
#' @exportClass Arena
#'
#' @slot orgdat A data frame collecting information about the accumulated biomass, type, phenotype, x and y position for each individual in the environment.
#' @slot specs A list of organism types and their associated parameters.
#' @slot media A list of objects of class \code{\link{Substance-class}} for each compound in the environment.
#' @slot phenotypes A list of unique phenotypes (metabolites consumed and produced), which occurred in the environment.
#' @slot mediac A character vector containing the names of all substances in the environment.
#' @slot tstep A number giving the time (in h) per iteration.
#' @slot stir A boolean variable indicating if environment should be stirred. If true, bacteria move to random positions within the environment and substances have a uniform concentration value.
#' @slot mflux A vector containing highly used metabolic reactions within the arena
#' @slot exchanges A data.frame containing last exchanges of each organism.
#' @slot shadow A vector containing shadow prices of metabolites present in the arena
#' @slot n A number giving the horizontal size of the environment.
#' @slot m A number giving the vertical size of the environment.
#' @slot Lx A number giving the horizontal grid size in cm.
#' @slot Ly A number giving the vertical grid size in cm.
#' @slot gridgeometry A list containing grid geometry parameter 
#' @slot seed An integer refering to the random number seed used to be reproducible
#' @slot scale A numeric defining the scale factor used for intern unit conversion.
#' @slot models A list containing Objects of class sybil::modelorg which represent the genome scale metabolic models
#' @slot occupyM A matrix indicating grid cells that are obstacles
#' @slot sublb A data matrix containing positions with amounts of substance for all organism
#' @slot removeM A matrix indicating grid cells from which organisms are removed (i.e. killed) after each time step
setClass("Arena",
         representation(
           orgdat="data.frame",
           specs="list",
           media="list",
           phenotypes="character",
           mediac="character",
           tstep="numeric",
           stir="logical",
           mflux="list",
           exchanges="data.frame",
           shadow="list",
           n="numeric",
           m="numeric",
           gridgeometry="list",
           Lx="numeric",
           Ly="numeric",
           seed="numeric",
           scale="numeric",
           models="list",
           occupyM="matrix",
           sublb="matrix",
           removeM="matrix"
        ),
        prototype(
          orgdat = data.frame(biomass=numeric(0),type=integer(0),phenotype=integer(0),x=integer(0),y=integer(0)),
          specs = list(),
          media = list(),
          phenotypes = character(),
          mediac = character(),
          tstep = 1,
          stir = F,
          mflux = list(),
          shadow = list(),
          models=list(),
          sublb=matrix()
        )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Arena-class}}
#' 
#' @export
#' @name Arena-constructor
#' @param n A number giving the horizontal size of the environment.
#' @param m A number giving the vertical size of the environment.
#' @param Lx A number giving the horizontal grid size in cm.
#' @param Ly A number giving the vertical grid size in cm.
#' @param seed An integer refering to the random number seed used to be reproducible
#' @param ... Arguments of \code{\link{Arena-class}}
Arena <- function(Lx=NULL, Ly=NULL, n=100, m=100, seed=sample(1:10000,1), ...){
  if(is.null(Lx)) Lx <- 0.025/100 * n
  if(is.null(Ly)) Ly <- 0.025/100 * m
  
  set.seed(seed) # remember random seed
  gridgeometry = list(grid2D=ReacTran::setup.grid.2D(ReacTran::setup.grid.1D(x.up = 0, L = Ly, N = m), 
                                                     ReacTran::setup.grid.1D(x.up = 0, L = Lx, N = n)))
  scale   <- (Lx*Ly)/(n*m)
  occupyM <- matrix(0, ncol=n, nrow=m)
  removeM <- matrix(0, ncol=n, nrow=m)
  new("Arena", Lx=Lx, Ly=Ly, n=n, m=m, scale=scale, gridgeometry=gridgeometry, occupyM=occupyM, removeM=removeM, seed=seed, ...)
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
setGeneric("tstep", function(object){standardGeneric("tstep")})
setMethod("tstep", "Arena", function(object){return(object@tstep)})
setGeneric("stir", function(object){standardGeneric("stir")})
setMethod("stir", "Arena", function(object){return(object@stir)})
setGeneric("mflux", function(object){standardGeneric("mflux")})
setMethod("mflux", "Arena", function(object){return(object@shadow)})
setGeneric("shadow", function(object){standardGeneric("shadow")})
setMethod("shadow", "Arena", function(object){return(object@mflux)})
setGeneric("n", function(object){standardGeneric("n")})
setMethod("n", "Arena", function(object){return(object@n)})
setGeneric("m", function(object){standardGeneric("m")})
setMethod("m", "Arena", function(object){return(object@m)})
setGeneric("gridgeometry", function(object){standardGeneric("gridgeometry")})
setMethod("gridgeometry", "Arena", function(object){return(object@gridgeometry)})
setGeneric("scale", function(object){standardGeneric("scale")})
setMethod("scale", "Arena", function(object){return(object@scale)})
setGeneric("seed", function(object){standardGeneric("seed")})
setMethod("seed", "Arena", function(object){return(object@seed)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################


#' @title Add individuals to the environment
#'
#' @description The generic function \code{addOrg} adds individuals to the environment.
#' @export
#' @rdname addOrg
#'
#' @param object An object of class Arena.
#' @param specI An object of class Organism.
#' @param amount A numeric number giving the number of individuals to add.
#' @param x A numeric vector giving the x positions of individuals on the grid.
#' @param y A numeric vector giving the y positions of individuals on the grid.
#' @param posmat A As an alternative to parameter x, y, a matrix with corrdinates an be specified
#' @param n0 Start column of matrix to take free positions from (default 1)
#' @param m0 Start row of matrix to take free positions from (default 1)
#' @param m End row of matrix to take free positions from (default arena@n)
#' @param n End column of matrix to take free positions from (default arena@m)
#' @param biomass A numeric vector giving the starting biomass of the individuals. (unit: fg)
#' @details The arguments \code{x} and \code{y} should be in the same length as the number of organisms added (given by the argument \code{amount}).
#' @seealso \code{\link{Arena-class}} and \code{\link{Bac-class}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' 
#' # Alternative way: adding organisms by giving matrix with positions
#' arena <- Arena(n=20,m=20)
#' mat <- matrix(sample(c(0,1), 400, replace = TRUE), nrow = 20, ncol = 20)
#' bac <- Bac(Ec_core)
#' arena <- addOrg(arena, specI=bac, posmat = mat)
#' 
setGeneric("addOrg", function(object, specI, amount=1, x=NULL, y=NULL, posmat=NULL, biomass=NA, n0=NULL, n=NULL, m0=NULL, m=NULL){standardGeneric("addOrg")})
#' @export
#' @rdname addOrg
setMethod("addOrg", "Arena", function(object, specI, amount=1, x=NULL, y=NULL, posmat=NULL, biomass=NA, n0=NULL, n=NULL, m0=NULL, m=NULL){
  switch(class(object),"Arena"={object <- object}, "Eval"={object <- getArena(object)}, stop("Please supply an object of class Arena or Eval."))
  if(length(posmat)>0){
    if(nrow(posmat)!=object@m | ncol(posmat)!=object@n){
      stop("Matrix posmat has invalid dimensions (should be equal to arena)")}
      idx <- which(posmat==1, arr.ind=TRUE)
      x <- idx[,2]
      y <- idx[,1]
      amount <- sum(posmat)
  }
  if(length(n)==0) n <- object@n; if(length(m)==0) m <- object@m
  if(length(n0)==0) n0 <- 1; if(length(m0)==0) m0 <- 1
  if(n0>n | m0>m) stop("Infeasible n0>n or m0>m set")
  if(amount+nrow(object@orgdat) > n*m){
    stop("More individuals than space on the grid")
  }
  bacnum <- round(object@scale/(specI@cellarea*10^(-8)))
  if(bacnum<1){
    stop("Physical arena size (Lx, Ly) too small. Maximal amount of cells in one grid cell would be zero.")
  }
  if( specI@maxweight <= specI@minweight*2 )
    stop("Maxweight needs to be bigger than two-times minweight otherwise cells will die immediately after duplication")
  
  # check if name is already used by other organism
  idx.dupl <- match(specI@type, names(object@specs))
  if( !is.na(idx.dupl) ){
    if( !identical( specI@model, object@specs[[idx.dupl]]@model ) ){
      spectype <- paste0(specI@type, "_",length(grep(specI@type, names(object@specs))))  
      specI@type <- spectype # needs to be re-set!
      print(paste0("Organism of the same type but with different model already present, added a new one:", spectype))
    }
    else{
      spectype <- specI@type  
      print("Organism of the same type and model already present, merged new organism with old one")
    } 
  } 
  else{
    spectype <- specI@type
  } 
  
  neworgdat <- object@orgdat
  newspecs <- object@specs
  newspecs[[spectype]] <- specI
  type <- which(names(newspecs)==spectype)
  newmflux <- object@mflux
  newshadow <- object@shadow

  # mflux
  newmflux[[spectype]] <- numeric(length(specI@lbnd))
  names(newmflux[[spectype]]) <- names(specI@lbnd)
  #shadow
  ex=sybil::findExchReact(specI@model)
  newshadow[[spectype]] <- numeric(length(ex)+1)
  names(newshadow[[spectype]]) <- c(ex@met_id, specI@rbiomass)

  type <- which(names(newspecs)==spectype) 
  lastind <- nrow(object@orgdat)
  if(length(x*y)==0){
    cmbs = expand.grid(n0:n,m0:m)
    rownames(cmbs) = paste(cmbs[,1],cmbs[,2],sep='_')
    taken <- paste(object@orgdat$x,object@orgdat$y,sep='_')
    obstacles <- which(object@occupyM>0, arr.ind = TRUE) 
    taken <- c(taken, paste(obstacles[,1], obstacles[,2],sep="_")) # extend taken to contain obstacle grid cells
    not_empty <- which(rownames(cmbs) %in% taken)
    if(length(not_empty) > 0){
      cmbs <- cmbs[-which(rownames(cmbs) %in% taken),]}
    sel <- sample(1:nrow(cmbs),amount)
    xp = cmbs[sel,1]
    yp = cmbs[sel,2]
    neworgdat[(lastind+1):(amount+lastind),'x']=xp
    neworgdat[(lastind+1):(amount+lastind),'y']=yp
    if(is.numeric(biomass)) neworgdat[(lastind+1):(amount+lastind),'biomass'] = rep(biomass, amount)
    else neworgdat[(lastind+1):(amount+lastind),'biomass'] = abs(rnorm(amount, mean=specI@cellweight_mean, sd=specI@cellweight_sd))
    neworgdat[(lastind+1):(amount+lastind),'type']=rep(type, amount)
    neworgdat[(lastind+1):(amount+lastind),'phenotype']=rep(NA, amount)
  }else{
    if( any(x<1) || any(x>n) || any(y<1) || any(y>m) ){stop("The positions of the individuals are beyond the dimensions of the environment.")}
    neworgdat[(lastind+1):(amount+lastind),'x']=x
    neworgdat[(lastind+1):(amount+lastind),'y']=y
    if(is.numeric(biomass)) neworgdat[(lastind+1):(amount+lastind),'biomass'] = rep(biomass, amount)
    else neworgdat[(lastind+1):(amount+lastind),'biomass'] = abs(rnorm(amount, mean=specI@cellweight_mean, sd=specI@cellweight_sd))
    neworgdat[(lastind+1):(amount+lastind),'type']=rep(type, amount)
    neworgdat[(lastind+1):(amount+lastind),'phenotype']=rep(NA, amount)
  }
  if(sum(duplicated(paste(neworgdat$x,neworgdat$y,sep="_")))!=0){
    stop("You have multiple individuals in the same position! Make sure that your x an y positions are unique")
  }
  #add initial medium (without concentration) for each organism
  newmet <- specI@medium[which(!specI@medium %in% object@mediac)]
  if(length(newmet) > 0){
    newmedia = list()
    for(i in 1:length(newmet)){
      newmedia[[unname(newmet[i])]] <- Substance(object@n, object@m, smax=0, id=unname(newmet[i]), name=names(newmet[i]), gridgeometry=object@gridgeometry, occupyM=object@occupyM)
    }
    object@media <- c(object@media,newmedia)
    object@media <- object@media[unique(names(object@media))]    
    newmediac <- c(object@mediac, specI@medium)
    object@mediac <- newmediac[!duplicated(newmediac)]
  }
  object@orgdat <- neworgdat
  object@specs <- newspecs
  object@mflux <- newmflux
  object@shadow <- newshadow
  object@models <- c(object@models, specI@model)
  return(object)
})

#' @title Add substances to the environment
#'
#' @description The generic function \code{addSubs} adds specific substances to the environment.
#' @export
#' @rdname addSubs
#'
#' @param object An object of class Arena.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param smax A numeric vector indicating the maximum substance concentration per grid cell.
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
#' @param difunc A character vector ("pde","cpp" or "r") describing the function for diffusion.
#' @param pde Choose diffusion transport reaction to be used (default is diffusion only)
#' @param difspeed A number indicating the diffusion rate (given by cm^2/h). Default is set to glucose diffusion in a aqueous solution (6.7e-6 cm^2/s * 3600 s/h = 0.02412 cm^2/h )
#' @param add A boolean variable defining whether the amount of substance should be summed or replaced
#' @param diffmat A matrix with spatial distributed initial concentrations (if not set, a homogenous matrix using smax is created)
#' @param template True if diffmat matrix should be used as tempalte only (will be multiplied with smax to obtain cocentrations)
#' @param Dgrid A matrix indicating the diffusion speed in x and y direction (given by cm^2/h).
#' @param Vgrid A number indicating the advection speed in x direction (given by cm/h).
#' @param addAnyway If true substance will be added even if there is no connection (i.e. exchanges) with organisms
#' @details If nothing but \code{object} is given, then all possible substrates are initilized with a concentration of 0. Afterwards, \code{\link{changeSub} can be used to modify the concentrations of specific substances.} 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) #add glucose, o2, pi
setGeneric("addSubs", function(object, smax=0, mediac=object@mediac, difunc="pde", pde="Diff2d", difspeed=0.02412, unit="mmol/cell", add=TRUE, diffmat=NULL, template=FALSE, Dgrid=NULL, Vgrid=NULL, addAnyway=FALSE){standardGeneric("addSubs")})
#' @rdname addSubs
#' @export
setMethod("addSubs", "Arena", function(object, smax=0, mediac=object@mediac, difunc="pde", pde="Diff2d", difspeed=0.02412, unit="mmol/cell", add=TRUE, diffmat=NULL, template=FALSE, Dgrid=NULL, Vgrid=NULL, addAnyway=FALSE){
  #switch(class(object),"Arena"={object <- object}, "Eval"={arenalast <- getArena(object)}, stop("Please supply an object of class Arena or Eval."))
  if(!(class(object)=="Arena" || class(object)=="Eval")){stop("Please supply an object of class Arena or Eval.")}
  if(length(smax) != length(mediac) && length(smax) != 1){
    stop("The parameter smax should be of the same size of mediac or equal to 1.")
  }
  if(class(mediac)=="factor") mediac <- as.character(mediac)
  newmet <- setdiff(mediac, object@mediac); names(newmet) <- gsub("EX_", "", newmet)
  if(addAnyway & length(newmet) > 0){  # add substance even if there is no exchange reaction
    newmedia = list()
    for(i in 1:length(newmet)){
      if(length(Dgrid)  > 1) Dgrid_i <- Dgrid[[i]] else Dgrid_i <- Dgrid
      if(length(Vgrid)  > 1) Vgrid_i <- Vgrid[[i]] else Vgrid_i <- Vgrid
      if(length(difunc) > 1) difunc_i <- difunc[[i]] else difunc_i <- difunc
      if(length(pde) > 1) pde_i <- pde[[i]] else pde_i <- pde
      newmedia[[unname(newmet[i])]] <- Substance(object@n, object@m, smax=0, id=unname(newmet[i]), 
                                                 name=names(newmet[i]), gridgeometry=object@gridgeometry, 
                                                 occupyM=object@occupyM, Dgrid=Dgrid_i, Vgrid=Vgrid_i,
                                                 difunc=difunc_i, pde=pde_i)}
    object@media  <- c(object@media,newmedia)
    object@mediac <- c(object@mediac, newmet)
  }
  if(sum(mediac %in% object@mediac) != length(mediac)){
    print(setdiff(mediac, object@mediac))
    warning("Substance does not exist in exchange reactions. It will not be added")
  }
  if(length(object@media)==0){
    stop("Organisms need to be defined first to determine what substances can be exchanged. In case you want to add substances without links to exchanges (addAnyway=TRUE) then please provide parameter mediac")
  }

  if(length(smax) != length(mediac))    {smax = rep(as.numeric(smax),length(mediac))}
  if(length(difspeed) != length(mediac)){difspeed = rep(difspeed,length(mediac))}
  
  # 1) consider units
  switch(unit,
         'mM'         = { conv <- 10^12 * 0.01 * object@scale}, # conversion of mMol in fmol/grid_cell
         'mmol/cm2'   = { conv <- 10^12 * object@scale}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/arena' = { conv <- 10^12 / (object@n*object@m)}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/cell'  = { conv <- 10^12}, # conversion of mmol/cell in fmol/cell
         'fmol/cell'  = { conv <- 1}, # already in fmol/cell
         stop("Wrong unit for concentration."))
  smax <- conv * smax
  if( is.numeric(diffmat) ){
    if( !template )
      diffmat <- conv * diffmat
    else
      smax    <- smax * (object@n * object@m) / sum(diffmat!=0)
  }
  if( is(diffmat, "list") ){
    if( !template )  
      diffmat <- lapply(diffmat, function(x){conv*x})
    else
      smax    <- sapply(seq_along(smax), function(i){smax[i] * (object@n * object@m) / sum(diffmat[[i]]) })
  }
  

  # 2) create and add substances assuming that organisms are already added
  
  if(class(object)=="Arena"){
    newmedia <- object@media
    for(i in 1:length(mediac)){
      if(mediac[[i]] %in% object@mediac){ # add only if possible
        old_diffmat <- object@media[[mediac[i]]]@diffmat
        new_name <- ifelse( length(names(mediac[i]))==0, object@media[[mediac[i]]]@name, names(mediac[i]) )  # add substance names 
        if(is(Dgrid, "list")) Dgrid_i  <- Dgrid[[i]] else Dgrid_i <- Dgrid
        if(is(Vgrid, "list")) Vgrid_i  <- Vgrid[[i]] else Vgrid_i <- Vgrid
        if(is(difunc, "list")) difunc_i <- difunc[[i]] else difunc_i <- difunc
        if(is(pde, "list"))    pde_i    <- pde[[i]] else pde_i <- pde
        if(is(diffmat, "list"))  diffmat_i<- diffmat[[i]] else diffmat_i <- diffmat
        if(is(template, "list")) template_i<- template[[i]] else template_i <- template
        object@media[[mediac[i]]] <- Substance(object@n, object@m, smax=smax[i], id=unname(mediac[i]), name=new_name, gridgeometry=object@gridgeometry, difunc=difunc_i, pde=pde_i, difspeed = difspeed[i], occupyM=object@occupyM, diffmat=diffmat_i, template=template_i, Dgrid=Dgrid_i, Vgrid=Vgrid_i)
        if(add){
          object@media[[mediac[i]]]@diffmat <- object@media[[mediac[i]]]@diffmat + old_diffmat}}}
  }else if(class(object)=="Eval"){ # if eval class then add substance in last time step of medlist
    time <- length(object@medlist)
    media <- lapply(object@media[names(object@medlist[[time]])], function(x, meds, n, m){
      x@diffmat <- Matrix::Matrix(meds[[x@id]],ncol=n,nrow=m,sparse=TRUE)
      return(x)
    },meds=extractMed(object,time), n=object@n, m=object@m)
    for(i in 1:length(mediac)){
      if(mediac[[i]] %in% object@mediac){ # add only if possible
        old_diffmat <- media[[mediac[i]]]@diffmat
        if(is(Dgrid, "list")) Dgrid_i <- Dgrid[[i]] else Dgrid_i <- Dgrid
        if(is(Vgrid, "list")) Vgrid_i <- Vgrid[[i]] else Vgrid_i <- Vgrid
        if(is(difunc, "list")) difunc_i <- difunc[[i]] else difunc_i <- difunc
        if(is(pde, "list")) pde_i <- pde[[i]] else pde_i <- pde
        if(is(diffmat, "list"))  diffmat_i<- diffmat[[i]] else diffmat_i <- diffmat
        if(is(template, "list")) template_i<- template[[i]] else template_i <- template
        new_name <- ifelse( length(names(mediac[i]))==0, object@media[[mediac[i]]]@name, names(mediac[i]) )  # add substance names 
        new_sub  <- Substance(object@n, object@m, smax=smax[i], id=unname(mediac[i]), name=new_name, gridgeometry=object@gridgeometry, difunc=difunc_i, pde=pde_i, difspeed = difspeed[i], occupyM=object@occupyM, diffmat=diffmat_i, template=template_i, Dgrid=Dgrid_i, Vgrid=Vgrid_i)
        if(add){
          object@medlist[[time]][[mediac[i]]] <- as.vector(new_sub@diffmat + old_diffmat)
        }else{
          object@medlist[[time]][[mediac[i]]] <- as.vector(new_sub@diffmat)
        }
        }}
  }
  
  # 3) return changed arena object
  return(object)
})

#' @title Change substances in the environment
#'
#' @description The generic function \code{changeSub} changes specific substances in the environment.
#' @export
#' @rdname changeSub
#'
#' @param object An object of class Arena.
#' @param smax A number or vector of numbers indicating the maximum substance concentration per grid cell.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
#' @details If nothing but \code{object} is given, then all possible substrates are initilized with a concentration of 0. Afterwards, \code{\link{changeSub}} can be used to modify the concentrations of specific substances.
#' @seealso \code{\link{Arena-class}} and \code{\link{addSubs}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena) #add all substances with no concentrations.
#' arena <- changeSub(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) 
#' #add substances glucose, oxygen and phosphate
setGeneric("changeSub", function(object, smax, mediac, unit="mmol/cell"){standardGeneric("changeSub")})
#' @rdname changeSub
#' @export
setMethod("changeSub", "Arena", function(object, smax, mediac, unit="mmol/cell"){
  warning("DEPRECATED: Please use addSubs()")
  if(length(smax)>1 & length(smax) != length(mediac)){
    stop("Number of substances does not match number of given concentrations")
  }
  if(length(smax) == 1){
    smax = rep(as.numeric(smax),length(mediac))
  }
  if(length(setdiff(mediac, names(object@media))) == 0 ){
    switch(unit,
           'mM'={smax <- 10^12 * smax* 0.01 * object@scale}, # conversion of mMol in fmol/grid_cell
           'mmol/cm2'={smax <- 10^12 * smax * object@scale}, # conversion of mmol/arena in fmol/grid_cell
           'mmol/arena'={smax <- 10^12 * smax / (object@n*object@m)}, # conversion of mmol/arena in fmol/grid_cell
           'mmol/cell'={smax <- 10^12 * smax}, # conversion of mmol/cell in fmol/cell
           stop("Wrong unit for concentration."))
    for(i in which(mediac %in% object@mediac)){
      object@media[mediac[i]] <- Substance(object@n, object@m, smax=smax[i], id=mediac[i], name=object@media[[mediac[i]]]@name,
                                                                  difunc=object@media[[mediac[i]]]@difunc,
                                                                  difspeed=object@media[[mediac[i]]]@difspeed, gridgeometry=object@gridgeometry, occupyM=object@occupyM)
      return(object)
    }
  }else stop("Substance does not exist in medium.")
})


#' @title Add default medium of an organism to arena.
#'
#' @description The generic function \code{addDefaultMed} uses the lower bounds defined in an organism's model file to compose minimal medium.
#' @export
#' @rdname addDefaultMed
#'
#' @param object An object of class Arena.
#' @param org An object of class Organism
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
setGeneric("addDefaultMed", function(object, org, unit="mM"){standardGeneric("addDefaultMed")})
#' @rdname addDefaultMed
#' @export
setMethod("addDefaultMed", "Arena", function(object, org, unit="mM"){
  lb_ex <- org@model@lowbnd[which(org@model@react_id %in% unname(org@medium))]
  lb_ex <- ifelse(lb_ex==-Inf,-1000,lb_ex) # change concentration if it is infinite
  min_id  <-  unname(org@medium[which(lb_ex < 0)])
  min_val <-  -lb_ex[which(lb_ex < 0)]
  switch(unit,
         'mM'         = { conv <- 10^12 * 0.01 * object@scale}, # conversion of mMol in fmol/grid_cell
         'mmol/cm2'   = { conv <- 10^12 * object@scale}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/arena' = { conv <- 10^12 / (object@n*object@m)}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/cell'  = { conv <- 10^12}, # conversion of mmol/cell in fmol/cell
         'fmol/cell'  = { conv <- 1}, # already in fmol/cell
         stop("Wrong unit for concentration."))
  min_val <- conv * min_val
  for(id in min_id){
    object@media[[id]]@diffmat = Matrix::Matrix(min_val[[which(min_id==id)]], nrow=object@m, ncol=object@n, sparse=TRUE)}
  return(object)
})

#' @title Add minimal medium of an organism to arena.
#'
#' @description The generic function \code{addEssentialMed} uses flux variability analysis to determine a essential growth medium  components.
#' @export
#' @rdname addEssentialMed
#'
#' @param object An object of class Arena.
#' @param org An object of class Organism
#' @param only_return Set true if essential metabolites should only be returned but not added to arena
#' @param limit A metabolite is considered as essential if its remove whould decrease biomass growth below limit (between 0,100; default 10\%)
setGeneric("addEssentialMed", function(object, org, only_return=FALSE, limit=10){standardGeneric("addEssentialMed")})
#' @rdname addEssentialMed
#' @export
setMethod("addEssentialMed", "Arena", function(object, org, only_return=FALSE, limit=10){
  model <- org@model
  ex <- sybil::findExchReact(model)
  lowbnd(model)[ex@react_pos] <- -1000 # set all lower bounds to -1000
  
  var_r <- sybil::fluxVar(org@model, react=ex@react_id, percentage=limit)
  ex_max <- sybil::maxSol(var_r, "lp_obj")
  min_id  <- ex@react_id[which(ex_max<0)]
  
  min_val <- -ex@lowbnd[which(ex_max<0)] # get lower bounds also for -INF cases
  min_val[min_val==Inf] <- 1000
  
  predefined <- sapply(object@media, function(x){sum(x@diffmat)}) # only use compounds which are not already in present in the media
  idx <- which(min_id %in% names(predefined)[which(predefined==0)])
  
  if(length(idx)==0) return(object)
  
  if(only_return) return(min_id[idx])
  
  for(id in intersect(min_id[idx], object@mediac)){
    object@media[[id]]@diffmat = Matrix::Matrix(min_val[[which(min_id==id)]], nrow=object@m, ncol=object@n, sparse=TRUE)}
  return(object)
})



#' @title Remove substances
#'
#' @description The generic function \code{rmSubs} removes all amounts of substances available in the arena for given compounds.
#' @export
#' @rdname rmSubs
#'
#' @param object An object of class Arena.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
setGeneric("rmSubs", function(object, mediac){standardGeneric("rmSubs")})
#' @rdname rmSubs
#' @export
setMethod("rmSubs", "Arena", function(object, mediac){
  object <- addSubs(object, smax=0, mediac=mediac, add=FALSE)
  return(object)
})


#' @title Remove all substances in the environment
#'
#' @description The generic function \code{flushSubs} removes specific substances in the environment.
#' @export
#' @rdname flushSubs
#'
#' @param object An object of class Arena.
#' @seealso \code{\link{Arena-class}} and \code{\link{addSubs}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena, smax=40) #add all substances with no concentrations.
#' arena <- changeSub(arena,20,c("EX_glc(e)","EX_o2(e)","EX_pi(e)")) 
#' #add substances glucose, oxygen and phosphate
#' arena <- flushSubs(arena) #remove all created substance concentrations
setGeneric("flushSubs", function(object){standardGeneric("flushSubs")})
#' @export
#' @rdname flushSubs
setMethod("flushSubs", "Arena", function(object){
  object@media <- list()
  return(object)
})

#' @title Change substance concentration patterns in the environment
#'
#' @description The generic function \code{changeDiff} changes specific substance concentration patterns in the environment.
#' @export
#' @rdname changeDiff
#'
#' @param object An object of class Arena.
#' @param newdiffmat A matrix giving the new gradient matrix of the specific substances in the environment.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @details This function can be used to add gradients of specific substances in the environment. The default conditions in \code{changeSubs} assumes an equal concentration in every grid cell of the environment. 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,30) #add all substances with no concentrations.
#' gradient <- matrix(1:200,20,20)
#' arena <- changeDiff(arena,gradient,c("EX_glc(e)","EX_o2(e)","EX_pi(e)"))
#' # add substances glucose, oxygen and phosphate
setGeneric("changeDiff", function(object, newdiffmat, mediac){standardGeneric("changeDiff")})
#' @export
#' @rdname changeDiff
setMethod("changeDiff", "Arena", function(object, newdiffmat, mediac){
  if(nrow(newdiffmat)==object@m && ncol(newdiffmat)==object@n){
    for(i in 1:length(mediac)){
      object@media[[mediac[i]]]@diffmat <- Matrix::Matrix(newdiffmat, sparse=TRUE)
      return(object)
    }
  }else stop("Given matrix is not compatible in dimensions with the environment.")
})

#' @title Change substance concentration patterns in the environment according to a gradient
#'
#' @description The generic function \code{createGradient} changes specific substance concentration patterns in the environment.
#' @export
#' @rdname createGradient
#'
#' @param object An object of class Arena.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @param position A character vector giving the position (top, bottom, right and left) of the gradient.
#' @param smax A number giving the maximum concentration of the substance.
#' @param add A boolean variable defining whether the amount of substance should be summed or replaced
#' @param steep A number between 0 and 1 giving the steepness of the gradient (concentration relative to the arena size).
#' @param unit A character used as chemical unit to set the amount of the substances to be added (valid values are: mmol/cell, mmol/cm2, mmol/arena, mM)
#' @details This function can be used to add gradients of specific substances in the environment. 
#' @seealso \code{\link{Arena-class}} and \code{\link{changeSub}} 
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,30) #add all substances with no concentrations.
#' arena <- createGradient(arena,smax=50,mediac=c("EX_glc(e)","EX_o2(e)","EX_pi(e)"),
#'              position='top',steep=0.5, add=FALSE)
setGeneric("createGradient", function(object, mediac, position, smax, steep, add=FALSE, unit='mmol/cell'){standardGeneric("createGradient")})
#' @export
#' @rdname createGradient
setMethod("createGradient", "Arena", function(object, mediac, position, smax, steep, add=FALSE, unit='mmol/cell'){
  if(steep<=0 || steep>=1){stop("Steepness must be in between 0 and 1.")}
  mediac = intersect(mediac,object@mediac)
  switch(unit,
         'mM'={smax <- 10^12 * smax* 0.01 * object@scale}, # conversion of mMol in fmol/grid_cell
         'mmol/cm2'={smax <- 10^12 * smax * object@scale}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/arena'={smax <- 10^12 * smax / (object@n*object@m)}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/cell'={smax <- 10^12 * smax}, # conversion of mmol/cell in fmol/cell
         stop("Wrong unit for concentration."))
  newdiffmat <- matrix(0,nrow=object@m,ncol=object@n)
  gradn = floor(object@n*steep)
  gradm = floor(object@m*steep)
  switch(position,
         'top'={for(i in 1:object@n){newdiffmat[0:gradm+1,i]=rev(seq(0,smax,length.out=gradm+1))}},
         'bottom'={for(i in 1:object@n){newdiffmat[object@n:(object@n-gradm),i]=rev(seq(0,smax,length.out=gradm+1))}},
         'right'={for(i in 1:object@m){newdiffmat[i,gradn:object@m]=seq(0,smax,length.out=gradn+1)}},
         'left'={for(i in 1:object@m){newdiffmat[i,0:gradn+1]=seq(smax,0,length.out=gradn+1)}},
         stop("Positions must be top, bottom, right, or left."))
  for(i in 1:length(mediac)){
    if(add){
      object@media[[mediac[i]]]@diffmat <- Matrix::Matrix(matrix(object@media[[mediac[i]]]@diffmat, nrow=object@m, ncol=object@n)+newdiffmat, sparse=TRUE)
    }else{
      object@media[[mediac[i]]]@diffmat <- Matrix::Matrix(newdiffmat, sparse=TRUE)
    }
  }
  return(object)
})

#' @title Change organisms in the environment
#'
#' @description The generic function \code{changeOrg} changes organisms in the environment.
#' @export
#' @rdname changeOrg
#'
#' @param object An object of class Arena.
#' @param neworgdat A data frame with new information about the accumulated biomass, type, phenotype, x and y position for each individual in the environment.
#' @details The argument \code{neworgdat} contains the same information as the \code{orgdat} slot of \code{\link{Arena-class}}. The \code{orgdat} slot of an \code{Arena} object can be used to create \code{neworgdat}.
#' @seealso \code{\link{Arena-class}} and \code{\link{addOrg}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' neworgdat <- arena@orgdat #get the current orgdat
#' neworgdat <- neworgdat[-1,] #remove the first individual
#' arena <- changeOrg(arena,neworgdat)
setGeneric("changeOrg", function(object, neworgdat){standardGeneric("changeOrg")})
#' @export
#' @rdname changeOrg
setMethod("changeOrg", "Arena", function(object, neworgdat){
  object@orgdat <- neworgdat
  return(object)
})

#' @title Function for checking phenotypes in the environment
#'
#' @description The generic function \code{checkPhen} checks and adds the phenotypes of organisms in the environment.
#' @export
#' @rdname checkPhen
#'
#' @param object An object of class Arena.
#' @param org An object of class Organism.
#' @param cutoff A number giving the cutoff for values of the objective function and fluxes of exchange reactions.
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' @return Returns a number indicating the number of the phenotype in the phenotype list.
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages. Uptake of substances are indicated by a negative and production of substances by a positive number.
#' @seealso \code{\link{Arena-class}} and \code{\link{getPhenotype}}
setGeneric("checkPhen", function(object, org, cutoff=1e-6, fbasol){standardGeneric("checkPhen")})
#' @export
#' @rdname checkPhen
setMethod("checkPhen", "Arena", function(object, org, cutoff=1e-6, fbasol){
  pind <- 0
  if(fbasol$obj>=cutoff){
    test = getPhenotype(org, cutoff=cutoff, fbasol)
    tspec = org@type
    pvec = rep(0,length(object@mediac))
    names(pvec) = object@mediac
    pvec[names(test)] = test
    pvec <- paste(pvec,collapse='')
    phenc <- object@phenotypes
    phensel <- phenc[which(names(phenc)==tspec)]
    pind <- which(phensel==pvec)
    if(length(pind)==0){
      pind = length(phensel)+1
      names(pvec) = tspec
      eval.parent(substitute(object@phenotypes <- c(phenc,pvec)))
    }
  }
  return(pind)
})

#' @title Function for checking phenotypes in the environment
#'
#' @description The generic function \code{checkPhen_par} checks and adds the phenotypes of organisms in the environment.
#' @export
#' @rdname checkPhen_par
#' 
#' @param object An object of class Arena.
#' @param org An object of class Organism.
#' @param cutoff A number giving the cutoff for values of the objective function and fluxes of exchange reactions.
#' @param fbasol Problem object according to the constraints and then solved with \code{optimizeProb}.
#' 
setGeneric("checkPhen_par", function(object, org, cutoff=1e-6, fbasol){standardGeneric("checkPhen_par")})
#' @export
#' @rdname checkPhen_par
setMethod("checkPhen_par", "Arena", function(object, org, cutoff=1e-6, fbasol){
  pind <- 0
  if(fbasol$obj>=cutoff){
    test = getPhenotype(org, cutoff=cutoff, fbasol)
    tspec = org@type
    pvec = rep(0,length(object@mediac))
    names(pvec) = object@mediac
    pvec[names(test)] = test
    pvec <- paste(pvec,collapse='')
    phenc <- object@phenotypes
    phensel <- phenc[which(names(phenc)==tspec)]
    pind <- which(phensel==pvec)
    if(length(pind)==0){
      pind = length(phensel)+1
      names(pvec) = tspec
      #return(list(pind, c(phenc,pvec)))
      return(list(pvec, TRUE))
    }
  }
  return(list(pind, FALSE))
})


# parallel helper funtion
setGeneric("addPhen", function(object, org, pvec){standardGeneric("addPhen")})
setMethod("addPhen", "Arena", function(object, org, pvec){
  tspec = org@type
  phenc <- object@phenotypes
  phensel <- phenc[which(names(phenc)==tspec)]
  pind = length(phensel)+1
  names(pvec) = tspec
  return(list(c(phenc,pvec), pind))
})


#' @title Function for unit conversion
#'
#' @description The generic function \code{unit_conversion} converts units for e.g. substance concentratins
#' @export
#' @rdname unit_conversion
#'
#' @param object An object of class Arena or Eval.
#' @param unit Unit to be converted to fmol/cell
#' @return Conversion factor
setGeneric("unit_conversion", function(object, unit){standardGeneric("unit_conversion")})
#' @export
#' @rdname unit_conversion
setMethod("unit_conversion", "Arena", function(object, unit){
  switch(unit,
         'mM'         = { conv <- 10^12 * 0.01 * object@scale}, # conversion of mMol in fmol/grid_cell
         'mmol/cm2'   = { conv <- 10^12 * object@scale}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/arena' = { conv <- 10^12 / (object@n*object@m)}, # conversion of mmol/arena in fmol/grid_cell
         'mmol/cell'  = { conv <- 10^12}, # conversion of mmol/cell in fmol/cell
         'fmol/cell'  = { conv <- 1}, # already in fmol/cell
         stop("Wrong unit for concentration."))
  return(conv)
})


#' @title Main function for simulating all processes in the environment
#'
#' @description The generic function \code{simEnv} for a simple simulation of the environment.
#' @export
#' @rdname simEnv
#'
#' @param object An object of class Arena or Eval.
#' @param time A number giving the number of iterations to perform for the simulation
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in the simEnv() by the function estimate_lrw())
#' @param continue A boolean indicating whether the simulation should be continued or restarted.
#' @param reduce A boolean indicating if the resulting \code{Eval} object should be reduced
#' @param diffusion True if diffusion should be done (default on).
#' @param diff_par True if diffusion should be run in parallel (default off).
#' @param cl_size If diff_par is true then cl_size defines the number of cores to be used in parallelized diffusion.
#' @param sec_obj character giving the secondary objective for a bi-level LP if wanted. Use "mtf" for minimizing total flux, "opt_rxn" for optimizing a random reaction, "opt_ex" for optimizing a random exchange reaction, and "sumex" for optimizing the sum of all exchange fluxes.
#' @param cutoff value used to define numeric accuracy
#' @param pcut A number giving the cutoff value by which value of objective function is considered greater than 0.
#' @param with_shadow True if shadow cost should be stored.
#' @param verbose Set to false if no status messages should be printed. 
#' @return Returns an object of class \code{Eval} which can be used for subsequent analysis steps.
#' @details The returned object itself can be used for a subsequent simulation, due to the inheritance between \code{Eval} and \code{Arena}. The parameter for sec_obj can be used to optimize a bi-level LP with a secondary objective if wanted. This can be helpful to subselect the solution space and create less alternative optimal solution. The secondary objective can be set to "mtf" to minimize the total flux, to simulate minimal enzyme usage of an organisms. If set to "opt_rxn" or "opt_ex", the secondary objective is picked as a random reaction or exchange reaction respectively everytime a fba is performed. This means that every individual of a population will select a different secondary reaction to optimize. The "sumex" option maximizes the secretion of products.
#' @seealso \code{\link{Arena-class}} and \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
setGeneric("simEnv", function(object, time, lrw=NULL, continue=FALSE, reduce=FALSE, diffusion=TRUE, diff_par=FALSE, cl_size=2, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=TRUE, verbose=TRUE){standardGeneric("simEnv")})
#' @export
#' @rdname simEnv
setMethod("simEnv", "Arena", function(object, time, lrw=NULL, continue=FALSE, reduce=FALSE, diffusion=TRUE, diff_par=FALSE, cl_size=2, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=TRUE, verbose=TRUE){
  if(length(object@media)==0) stop("No media present in Arena!")
  switch(class(object),
         "Arena"={arena <- object; evaluation <- Eval(arena)},
         "Eval"={arena <- getArena(object); evaluation <- object},
         stop("Please supply an object of class Arena or Eval."))
  if( sec_obj=="none" & any(sapply(object@specs, function(s){s@limit_growth})))
    warning("If growth is limited by maximum weight for an organism (max_weight=TRUE) it is recommended to use minimize total flux (sec_obj='mtf').")
  
  if(is.null(lrw)){lrw=estimate_lrw(arena@n,arena@m)}
  for(i in names(arena@specs)){
    phensel <- arena@phenotypes[which(names(arena@phenotypes)==i)]
    if(length(phensel)==0){
      test = getPhenotype(arena@specs[[i]], cutoff=pcut, fbasol=arena@specs[[i]]@fbasol)
      pvec = rep(0,length(arena@mediac))
      names(pvec) = arena@mediac
      pvec[names(test)] = test
      pvec <- paste(pvec,collapse='')
      names(pvec) = i
      arena@phenotypes <- c(arena@phenotypes,pvec)
    }
  }
  if(class(object)!="Eval"){addEval(evaluation, arena)}
  arena@sublb <- getSublb(arena)
  arena@exchanges <- data.frame() # remember exchanges
  diff_t=0
  if(arena@stir){ #create all possible positions on arena
    allxy = expand.grid(1:arena@n,1:arena@m)
    colnames(allxy) = c("x","y")
  }
  if(length(arena@specs) > 0) biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
  for(i in 1:time){
    init_t <- proc.time()[3]
    sublb <- arena@sublb
    if(nrow(arena@orgdat) > 1){
      new_ind = sample(1:nrow(arena@orgdat),nrow(arena@orgdat)) #shuffle through all bacteria to increase randomness
      arena@orgdat = arena@orgdat[new_ind,]
      sublb = sublb[new_ind,] #apply shuffeling also to sublb to ensure same index as orgdat
    }
    #if(verbose) cat("\niteration-start:", i, "\t organisms:",nrow(arena@orgdat), "\t biomass:", sum(arena@orgdat$biomass), "pg \n")
    #if(i==1){ 
    #org_stat <- sapply(seq_along(arena@specs), function(x){dim(arena@orgdat[which(arena@orgdat$type==x),])[1]})
    #if(length(arena@specs) > 0){
    #  old_biomass<-biomass_stat; biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
    #  org_stat <- cbind(org_stat, biomass_stat); rownames(org_stat) <- names(arena@specs); colnames(org_stat) <- c("count", "biomass")
    #  if(verbose) print(as.data.frame(org_stat))}}
    #if(verbose & i!=1) print(as.data.frame(org_stat)[,1:2])
    arena@mflux <- lapply(arena@mflux, function(x){numeric(length(x))}) # empty mflux pool
    arena@shadow <-lapply(arena@shadow, function(x){numeric(length(x))}) # empty shadow pool
    if(nrow(arena@orgdat) > 0){ # if there are organisms left
      org.count <- nrow(arena@orgdat)
      for(j in 1:org.count){ # for each organism in arena
        if(verbose) cat("\rOrganims",j,"/",org.count)
        org <- arena@specs[[arena@orgdat[j,'type']]]
        bacnum = round((arena@scale/(org@cellarea*10^(-8)))) #calculate the number of bacteria individuals per gridcell
        switch(class(org),
               "Bac"= {arena = simBac(org, arena, j, sublb, bacnum, sec_obj=sec_obj, cutoff=cutoff, pcut=pcut, with_shadow=with_shadow)}, #the sublb matrix will be modified within this function
               "Human"= {arena = simHum(org, arena, j, sublb, bacnum, sec_obj=sec_obj, cutoff=cutoff, pcut=pcut, with_shadow=with_shadow)}, #the sublb matrix will be modified within this function
               stop("Simulation function for Organism object not defined yet."))
      }
      test <- is.na(arena@orgdat$biomass)
      if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
      rm("test")
    }
    if(verbose) cat("\r")
    if(diffusion && !arena@stir){
      if(diff_par){
        diff_t <- system.time(arena <- diffuse_par(arena, cluster_size=cl_size, lrw=lrw, sublb=sublb) )[3]
      }else diff_t <- system.time(arena <- diffuse(arena, lrw=lrw, sublb=sublb, verbose=verbose) )[3]
    }
    if(!diffusion){
      if(nrow(sublb)>0){
        for(k in 1:length(arena@media)){
          for(l in 1:nrow(sublb)){
            arena@media[[k]]@diffmat[sublb[l,"y"],sublb[l,"x"]] = sublb[l,k+2] # first two columns are coordinates
          }
        }
      }
      arena@sublb <- getSublb(arena)
    }
    if(arena@stir){ #stir environment -> random movement of bacteria + perfect diffusion
      sublb_tmp = arena@orgdat[,c("x","y")]
      for(sub in names(arena@media)){ #go through each metabolite in medium
        sumc = sum(arena@media[[sub]]@diffmat) #sum of all concentrations
        meanc = sumc/(arena@n*arena@m) #mean per grid cell
        conc = ((sumc-(meanc*nrow(sublb)))+sum(sublb[,sub]))/(arena@n*arena@m) #remove concentrations where bacteria are sitting + add the current concentration in their position
        arena@media[[sub]]@diffmat = Matrix::Matrix(conc,nrow=arena@m,ncol=arena@n,sparse=TRUE) #create matrix with homogen concentration
        sublb_tmp[,sub] = conc #create a new sublb matrix
      }
      newpos = allxy[sample(1:nrow(allxy),nrow(arena@orgdat)),]
      arena@orgdat[,c('x','y')] = newpos
      sublb_tmp[,c('x','y')] = newpos
      arena@sublb = as.matrix(sublb_tmp)
    }
    idx.rm <- which(arena@removeM > 0, arr.ind=TRUE) # remove organisms given removal matrix
    if( nrow(idx.rm) > 0 ){
      idx.rm.str <- apply(idx.rm, 1, function(r){paste0(r,collapse=",")})
      idx.orgdat.str <- apply(arena@orgdat[,c('x','y')], 1, function(r){paste0(r,collapse=",")})
      rm.rows <- which(!is.na(match(idx.orgdat.str, idx.rm.str)))
      if( length(rm.rows)> 0 ){
        arena@orgdat <- arena@orgdat[-rm.rows,]
        if(verbose) cat("removed", length(rm.rows), "organisms\n")
      }
    }
    
    addEval(evaluation, arena)
    if(reduce && i<time){evaluation = redEval(evaluation)}
    if(nrow(arena@orgdat)==0 && !continue){
      if(verbose) print("All organisms died!")
      break
    }
    step_t <- proc.time()[3] - init_t
    if(verbose) cat("\niteration:", i, "\t organisms:",nrow(arena@orgdat), "\t biomass:", sum(arena@orgdat$biomass), "pg \n")
    if(verbose) cat("\r")
    org_stat <- sapply(seq_along(arena@specs), function(x){dim(arena@orgdat[which(arena@orgdat$type==x),])[1]})
    if(length(arena@specs) > 0){
      old_biomass<-biomass_stat; biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
      org_stat <- cbind(org_stat, biomass_stat, 100*(biomass_stat-old_biomass)/old_biomass); rownames(org_stat) <- names(arena@specs); colnames(org_stat) <- c("count", "biomass", "%")
      if(verbose) print(as.data.frame(org_stat))}
    if(verbose) cat("\r")
    if(verbose) cat("\ttime total: ", round(step_t,3), "\tdiffusion: ", round(diff_t,3), " (", 100*round(diff_t/step_t,3),"%)\n\n" )
    if(verbose) cat("--------------------------------------------------------------------\n")
    }
  return(evaluation)
})



#' @title Function for diffusion
#'
#' @description The generic function \code{diffuse} computes the media distribution via diffusion
#' @export
#' @rdname diffuse
#' 
#' @param object An object of class Arena.
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in the simEnv() by the function estimate_lrw())
#' @param sublb A matrix with the substrate concentration for every individual in the environment based on their x and y position.
#' @param verbose Set to false if no status messages should be printed. 
setGeneric("diffuse", function(object, lrw, sublb, verbose=TRUE){standardGeneric("diffuse")})
#' @export
#' @rdname diffuse
setMethod("diffuse", "Arena", function(object, lrw, sublb, verbose=TRUE){
  arena <- object
  diff_init_t <- proc.time()[3]
  sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
  #diff_pre_t <- system.time({
  if(!all(is.na(sublb)) & dim(sublb)[1] > 0){ # if there are organisms
    testdiff <- t(sublb[,-c(1,2)]) == unlist(lapply(arena@media,function(x,n,m){return(mean(x@diffmat))})) #check which mets in sublb have been changed by the microbes
    changed_mets <- which(apply(testdiff,1,sum)/nrow(sublb) < 1) #find the metabolites which are changed by at least one microbe
  } else changed_mets <- list()#})[3]
  #diff_pde_t=0; diff_sublb_t=0
  #diff_loop_t <- system.time({for(j in seq_along(arena@media)){
  for(j in seq_along(arena@media)){
    if(verbose) cat("\rSubstances",j,"/",length(arena@media))
    #skip diffusion if already homogenous (attention in case of boundary/source influx in pde!)
    if(length(changed_mets)>0) homogenous = !(j %in% changed_mets) else homogenous = FALSE
    diffspeed  = arena@media[[j]]@difspeed>0
    diff2d     = arena@media[[j]]@pde=="Diff2d"
    if(diff2d&&!homogenous || !diff2d){
      submat <- matrix(arena@media[[j]]@diffmat, nrow=object@m, ncol=object@n)
      if(!all(is.na(sublb)) && dim(sublb)[1] > 0 && (nrow(sublb) != sum(sublb[,j+2]==mean(submat)))){
        submat[sublb[,c("y","x")]] <- sublb[,arena@media[[j]]@id]
      }
      #diff_pde_t <- diff_pde_t + system.time(switch(arena@media[[j]]@difunc,
      if(diffspeed || !diff2d){
        switch(arena@media[[j]]@difunc,
               "pde"  = {submat <- diffusePDE(arena@media[[j]], submat, gridgeometry=arena@gridgeometry, lrw, tstep=object@tstep)},
               "pde2" = {diffuseSteveCpp(submat, D=arena@media[[j]]@difspeed, h=1, tstep=arena@tstep)},
               "naive"= {diffuseNaiveCpp(submat, donut=FALSE)},
               "r"    = {for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
               stop("Diffusion function not defined yet.")
        )
      }
      arena@media[[j]]@diffmat <- Matrix::Matrix(submat, sparse=TRUE)
    }else submat <- arena@media[[j]]@diffmat
    tryCatch({
      sublb_tmp[,j] <- submat[cbind(arena@orgdat$y,arena@orgdat$x)]
    }, error=function(cond){
      print(cond)
      browser()
    })
  }#})[3]
  if(verbose) cat("\r")
  sublb <- cbind(as.matrix(arena@orgdat[,c("y","x")]),sublb_tmp)
  colnames(sublb) <- c('y','x',arena@mediac)
  arena@sublb <- sublb
  
  #diff_t <- proc.time()[3] - diff_init_t
  #print(paste("diffusion time total", round(diff_t,3), "pre", round(diff_pre_t,3), "loop", round(diff_loop_t,3), "pde", round(diff_pde_t,3), "sublb", round(diff_sublb_t,3), "post", round(diff_post_t,3) ))

  
  #return(list(arena, sublb))
  return(arena)
})



#' @title Main function for simulating in parallel all processes in the environment
#'
#' @description The generic function \code{simEnv_par} for a simple in parallel all simulation of the environment.
#' @export
#' @rdname simEnv_par
#'
#' @param object An object of class Arena or Eval.
#' @param time A number giving the number of iterations to perform for the simulation
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in the simEnv() by the function estimate_lrw())
#' @param continue A boolean indicating whether the simulation should be continued or restarted.
#' @param reduce A boolean indicating if the resulting \code{Eval} object should be reduced
#' @param cluster_size Number of cpu cores to be used.
#' @param diffusion True if diffusion should be done (default on).
#' @param sec_obj character giving the secondary objective for a bi-level LP if wanted.
#' @param cutoff value used to define numeric accuracy
#' @param with_shadow True if shadow cost should be stores (default off).
#' @param verbose Set to false if no status messages should be printed. 
#' @return Returns an object of class \code{Eval} which can be used for subsequent analysis steps.
#' @details The returned object itself can be used for a subsequent simulation, due to the inheritance between \code{Eval} and \code{Arena}.
#' @seealso \code{\link{Arena-class}} and \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
setGeneric("simEnv_par", function(object, time, lrw=NULL, continue=FALSE, reduce=FALSE, cluster_size=NULL, diffusion=TRUE, sec_obj="none", cutoff=1e-6, with_shadow=FALSE, verbose=TRUE){standardGeneric("simEnv_par")})
#' @export
#' @rdname simEnv_par
setMethod("simEnv_par", "Arena", function(object, time, lrw=NULL, continue=FALSE, reduce=FALSE, cluster_size=NULL, diffusion=TRUE, sec_obj="none", cutoff=1e-6, with_shadow=FALSE, verbose=TRUE){
  if(length(object@media)==0) stop("No media present in Arena!")
  switch(class(object),
         "Arena"={arena <- object; evaluation <- Eval(arena)},
         "Eval"={arena <- getArena(object); evaluation <- object},
         stop("Please supply an object of class Arena."))
  if(is.null(lrw)){lrw=estimate_lrw(arena@n,arena@m)}
  for(i in names(arena@specs)){
    phensel <- arena@phenotypes[which(names(arena@phenotypes)==i)]
    if(length(phensel)==0){
      test = getPhenotype(arena@specs[[i]], cutoff=1e-6, fbasol=arena@specs[[i]]@fbasol)
      pvec = rep(0,length(arena@mediac))
      names(pvec) = arena@mediac
      pvec[names(test)] = test
      pvec <- paste(pvec,collapse='')
      names(pvec) = i
      arena@phenotypes <- c(arena@phenotypes,pvec)
    }
  }
  if(class(object)!="Eval"){addEval(evaluation, arena)}
  arena@sublb <- getSublb(arena)
  
  if(length(cluster_size)==0){
    cluster_size <- parallel::detectCores()
  }

  for(i in 1:time){
    diff_t=0; par_t=0; par_post_t=0
    
    init_t <- proc.time()[3]
    arena@orgdat["nr"] <- seq_len(dim(arena@orgdat)[1]) # dummy numbering
    if(verbose) cat("\nparallel iteration:", i, "\t organisms:",nrow(arena@orgdat), "\t biomass:", sum(arena@orgdat$biomass), "pg \n")
    org_stat <- table(arena@orgdat$type)
    names(org_stat) <- names(arena@specs)[as.numeric(names(org_stat))]
    if(verbose) print(org_stat)
    arena@mflux <- lapply(arena@mflux, function(x){numeric(length(x))}) # empty mflux pool
    arena@shadow <-lapply(arena@shadow, function(x){numeric(length(x))}) # empty shadow pool
    if(nrow(arena@orgdat) > 0){ # if there are organisms left
      #if(nrow(arena@orgdat) >= arena@n*arena@m) browser()
      test_init_t <- proc.time()[3]
      sublb <- arena@sublb
      #sublb[,arena@mediac] = sublb[,arena@mediac]*(10^12) #convert to fmol per gridcell
      test_t <- proc.time()[3] - test_init_t
      
      # 1) split orgdat into a data.frames for each species 
      split_orgdat <- split(arena@orgdat, as.factor(arena@orgdat$type))
      # 2) iterate over all species (each has a entry in splited data.frame)
      lapply(1:length(split_orgdat), function(spec_nr){
        splited_species <- split_orgdat[[spec_nr]]
        splited_size <- dim(splited_species)[1]
        # 2.1) in case of big splited data frame go for parallel
        if(splited_size >= 1){ # ATTENTION: magic number, to be defined according to benchmark! (treshold from which parallel is faster than seriell)
          # 2.1.1) group task (each core gets one)
          groups <- split(seq_len(splited_size), cut(seq_len(splited_size), cluster_size))
          # 2.1.2) paralel loop
          names(groups) <- NULL
          #parallel_sol <- lapply(groups, function(g){
          #parallel_sol <- parallel::parLapply(parallelCluster, groups, function(g){
          par_t <<- par_t + system.time(parallel_sol <- parallel::mclapply(groups, function(g){
          #parallel_sol <- parallel::mclapply(groups, function(g){
                                      # 2.1.2.1) critical step: create lpobject for each core 
                                      #(otherwise pointer will corrupt in warm-started optimization)
                                      model <- arena@specs[[spec_nr]]@model
                                      lpobject <- sybil::sysBiolAlg(model, algorithm="fba")
                                      # 2.1.2.2) 
                                      test <- lapply(g, function(i){
                                          org <- arena@specs[[spec_nr]]
                                          bacnum = round((arena@scale/(org@cellarea*10^(-8))))
                                          j <- splited_species$nr[i] # unique id in orgdat (necessary due to split in parallel mode)
                                          simbac <- simBac_par(org, arena, j, sublb, bacnum, lpobject, sec_obj=sec_obj, cutoff=cutoff, with_shadow=with_shadow)
                                          neworgdat <- simbac[[1]]
                                          sublb <- simbac[[2]]
                                          fbasol <- simbac[[3]]
                                          phen_res <- checkPhen_par(arena, org, fbasol=fbasol)
                                          #neworgdat$phenotype <- phen_res[[1]]
                                          todo_pheno_nr <- NULL
                                          if(phen_res[[2]]==TRUE) todo_pheno_nr <- j # unique new phenotype id cannot be determined in parallel
                                          todo_pheno <- phen_res[[1]]
                                          list("neworgdat"=neworgdat, "sublb"=sublb, "fbasol"=fbasol, "todo_pheno"=todo_pheno, "todo_pheno_nr"=todo_pheno_nr)
                                        })
                                      list("neworgdat"=sapply(test, with, test$neworgdat), "sublb"=sapply(test, with, test$sublb), "fbasol_flux"=sapply(test, with, fbasol$fluxes), "todo_pheno"=sapply(test, with, test$todo_pheno), "todo_pheno_nr"=sapply(test, with, test$todo_pheno_nr))
          #})
          }, mc.cores=cluster_size))[3]
          
          par_post_init_t <- proc.time()[3]
          tmpnames <- colnames(arena@orgdat)
          orgdat2 <- data.frame(matrix(unlist(sapply(parallel_sol, with, parallel_sol$neworgdat)), ncol=dim(arena@orgdat)[2], byrow=TRUE))
          colnames(orgdat2) <- tmpnames
          if(all(apply(orgdat2, 1, is.numeric)) != TRUE) browser()
          arena@orgdat <<- orgdat2

          tmpnames <- colnames(sublb)
          sublb2 <- matrix(unlist(sapply(parallel_sol, with, parallel_sol$parallel_solsublb)), ncol=dim(sublb)[2], byrow=TRUE)
          colnames(sublb2) <- tmpnames
          sublb <<- sublb2
          
          fba_fluxes <- sapply(parallel_sol, with, parallel_sol$parallel_solfbasol_flux)
          arena@mflux[[names(arena@specs)[[spec_nr]]]] <<- arena@mflux[[names(arena@specs)[[spec_nr]]]] + colSums(matrix(unlist(fba_fluxes), ncol=length(arena@mflux[[names(arena@specs)[[spec_nr]]]]), byrow = TRUE)) # remember active fluxes
          #arena@shadow[[names(arena@specs)[[spec_nr]]]] <<- arena@shadow[[names(arena@specs)[[spec_nr]]]] + colSums(matrix(unlist(fba_fluxes), ncol=length(arena@mflux[[names(arena@specs)[[spec_nr]]]]), byrow = TRUE)) # remember active fluxes
          todo_pheno <- sapply(parallel_sol, with, parallel_sol$parallel_soltodo_pheno)
          todo_pheno <- as.numeric(unname(unlist(todo_pheno)))
          todo_pheno_nr <- sapply(parallel_sol, with, parallel_sol$parallel_soltodo_pheno_nr)
          todo_pheno_nr <- unlist(todo_pheno_nr)
          if(all(is.na(arena@orgdat$phenotype))) arena@orgdat$phenotype <<- todo_pheno # init case
          if(length(todo_pheno_nr) > 0){ # handle new phenotypes
            unique_todo_pheno <- unique(todo_pheno[todo_pheno_nr])
            lapply(unique_todo_pheno, function(pvec){
              res_addPhen <- addPhen(arena, org=arena@specs[[spec_nr]], pvec)
              arena@phenotypes <<- res_addPhen[[1]]
              arena@orgdat$phenotype[which(todo_pheno == pvec)] <<- res_addPhen[[2]]
            })
          }
          par_post_t <<- par_post_t + (proc.time()[3] - par_post_init_t)
        # 2.2) in case of small splited data frame do seriell work
        }else{
          stop("to be done")
        }
      })
      #movdup_t <- system.time({
      arena@orgdat <- arena@orgdat[,-which(colnames(arena@orgdat)=="nr")] # remove dummy numbering
      movementCpp(arena@orgdat, arena@n, arena@m, arena@occupyM) # call by ref
      arena@orgdat <- duplicateCpp(arena@orgdat, arena@n, arena@m, lapply(arena@specs, function(x){x@maxweight}), arena@occupyM) # call by val
      #})[3]
      
      # delete dead organisms
      test <- is.na(arena@orgdat$biomass)
      if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
      #rm("test")
    }
    
    if(diffusion) diff_t <- system.time(arena <- diffuse_par(arena, lrw, cluster_size, sublb) )[3]

    addEval(evaluation, arena)
    if(reduce && i<time){evaluation = redEval(evaluation)}
    if(nrow(arena@orgdat)==0 && !continue){
      if(verbose) print("All organisms died!")
      break
    }
    step_t <- proc.time()[3] - init_t
    #cat("\ttime total: ", round(step_t,3), "\tdiffusion: ", round(diff_t,3), " (", 100*round(diff_t/step_t,3),"%)\n" )
    if(verbose) cat("\ttime total: ", round(step_t,3), "\tdiffusion: ", round(diff_t,3), " (", 100*round(diff_t/step_t,3),"%)", "\tpar_fba: ", round(par_t,3), " (", 100*round(par_t/step_t,3),"%)", "\tpar_fba_post: ", round(par_post_t,3), " (", 100*round(par_post_t/step_t,3),"%)\n")
    if(verbose) cat("\ttest:", round(test_t,3), " (", 100*round(test_t/step_t,3),"%)")
  }
  #parallel::stopCluster(parallelCluster)
  return(evaluation)
})


#' @title Function for parallelzied diffusion
#'
#' @description The generic function \code{diffuse_par} computes the media distribution via diffusion in parallel
#' @export
#' @rdname diffuse_par
#' 
#' @param object An object of class Arena.
#' @param lrw A numeric value needed by solver to estimate array size (by default lwr is estimated in the simEnv() by the function estimate_lrw())
#' @param cluster_size Amount of cores to be used
#' @param sublb A matrix with the substrate concentration for every individual in the environment based on their x and y position.
setGeneric("diffuse_par", function(object, lrw, cluster_size, sublb){standardGeneric("diffuse_par")})
#' @export
#' @rdname diffuse_par
setMethod("diffuse_par", "Arena", function(object, lrw, cluster_size, sublb){
  diff_init_t <- proc.time()[3]
  cl <- parallel::makeCluster(cluster_size, type="PSOCK")
  #parallel::clusterExport(cl, c(""))
  arena <- object
  sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
  #diff_pre_t <- system.time({ if(dim(sublb)[1] > 0){
  if(!all(is.na(sublb)) && dim(sublb)[1] > 0){ # if there are organisms
      testdiff <- t(sublb[,-c(1,2)]) == unlist(lapply(arena@media,function(x,n,m){return(mean(x@diffmat))})) #check which mets in sublb have been changed by the microbes
      changed_mets <- which(apply(testdiff,1,sum)/nrow(sublb) < 1) #find the metabolites which are changed by at least one microbe
  } else changed_mets <- list()#})[3]
  #diff_pde_t=0; diff_sublb_t=0
  #diff_loop_t <- system.time(parallel_diff <- parallel::mclapply(seq_along(arena@media), function(j){
  #parallel_diff <- parallel::mclapply(seq_along(arena@media), function(j){
  parallel_diff  <- parallel::parLapply(cl, seq_along(arena@media), function(j){
  #parallel_diff <- lapply(seq_along(arena@media), function(j){
  #diff_loop_t <- system.time(parallel_diff <-  foreach(j=seq_along(arena@media)) %dopar% {
    #skip diffusion if already homogenous (attention in case of boundary/source influx in pde!)
    if(length(changed_mets)>0) homogenous = !(j %in% changed_mets) else homogenous = FALSE
    diffspeed  = arena@media[[j]]@difspeed>0
    diff2d     = arena@media[[j]]@pde=="Diff2d"
    if(diff2d&&!homogenous || !diff2d){
      submat <- matrix(arena@media[[j]]@diffmat, nrow=arena@m, ncol=arena@n)
      if(!all(is.na(sublb)) && dim(sublb)[1] > 0 && (nrow(sublb) != sum(sublb[,j+2]==mean(submat)))){
        #diff_sublb_t <<- diff_sublb_t + system.time(submat[sublb[,c("y","x")]] <- sublb[,arena@media[[j]]@id])[3]}
        submat[sublb[,c("y","x")]] <- sublb[,arena@media[[j]]@id]}
      #browser()
      #diff_pde_t <<- diff_pde_t + system.time(switch(arena@media[[j]]@difunc,
      if(diffspeed || !diff2d){
        switch(arena@media[[j]]@difunc,
               "pde"  = {submat <- diffusePDE(arena@media[[j]], submat, gridgeometry=arena@gridgeometry, lrw, tstep=object@tstep)},
               "pde2" = {diffuseSteveCpp(submat, D=arena@media[[j]]@difspeed, h=1, tstep=arena@tstep)},
               "naive"= {diffuseNaiveCpp(submat, donut=FALSE)},
               "r"    = {for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
               stop("Diffusion function not defined yet."))#)[3]
      }
        diffmat_tmp <- Matrix::Matrix(submat, sparse=TRUE)
    }else{
      diffmat_tmp <- arena@media[[j]]@diffmat
      submat <- matrix(arena@media[[j]]@diffmat, nrow=arena@m, ncol=arena@n)
    }
    sublb_tmp  <- submat[cbind(arena@orgdat$y,arena@orgdat$x)]
    list("diffmat"=diffmat_tmp, "sublb"=sublb_tmp)
  #})#)[3]
  #}, mc.cores=cluster_size)#)[3]
  })
  parallel::stopCluster(cl)
  
  #diff_post_t <- system.time({ for(j in seq_along(arena@media)){
  for(j in seq_along(arena@media)){
      arena@media[[j]]@diffmat <- parallel_diff[[j]][[1]]
      sublb_tmp[,j] <- parallel_diff[[j]][[2]]
  }#})[3]
  sublb <- cbind(as.matrix(arena@orgdat[,c(4,5)]),sublb_tmp)
  colnames(sublb) <- c('x','y',arena@mediac)
  arena@sublb <- sublb
  #diff_t <- proc.time()[3] - diff_init_t
  #print(paste("diffusion time total", round(diff_t,3), "pre", round(diff_pre_t,3), "loop", round(diff_loop_t,3), "pde", round(diff_pde_t,3), "sublb", round(diff_sublb_t,3), "post", round(diff_post_t,3) ))
  return(arena)
})




#' @title Function for calculated the substrate concentration for every organism
#'
#' @description The generic function \code{getSublb} calculates the substrate concentration for every individual in the environment based on their x and y position.
#' @export
#' @rdname getSublb
#'
#' @param object An object of class Arena.
#' @return Returns the substrate concentration for every individual in the environment with substrates as well as x and y positions as columns and rows for each organism.
#' @seealso \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' sublb <- getSublb(arena)
setGeneric("getSublb", function(object){standardGeneric("getSublb")})
#' @export
#' @rdname getSublb
setMethod("getSublb", "Arena", function(object){
  sublb <- matrix(0,nrow=nrow(object@orgdat),ncol=(length(object@mediac)))
  for(j in seq_along(object@media)){
    submat <- matrix(object@media[[j]]@diffmat, nrow=object@m, ncol=object@n)
    sublb[,j] <- apply(object@orgdat, 1, function(x,sub){
      tryCatch({return(sub[as.numeric(x[5]),as.numeric(x[4])])
      }, error=function(cond){
        print(cond)
        print(x)
        browser()}
      )
    },sub=submat)
  }
  sublb <- cbind(as.matrix(object@orgdat[,c("y","x")]),sublb)
  colnames(sublb) <- c('y','x',object@mediac)
  return(sublb)
})

#' @title Function for stirring/mixing the complete evironment
#'
#' @description The generic function \code{stirEnv} simulates the event of mixing all substrates and organisms in the environment.
#' @export
#' @rdname stirEnv
#'
#' @param object An object of class Arena.
#' @param sublb A matrix with the substrate concentration for every individual in the environment based on their x and y position.
#' @return Returns the substrate concentration for every individual in the environment with substrates as well as x and y positions as columns and rows for each organism.
#' @details The stirring is implemented as a random permutation of organism positions and the equalization of of all substrate concentrations.
#' @seealso \code{\link{Arena-class}} and \code{\link{getSublb}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' sublb <- getSublb(arena)
#' stirEnv(arena,sublb)
setGeneric("stirEnv", function(object, sublb){standardGeneric("stirEnv")})
#' @export
#' @rdname stirEnv
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
  if(length(sit) != 0){neworgdat <- rbind(neworgdat,sitorgdat)}
  eval.parent(substitute(object@orgdat <- neworgdat))
  #stir all the substrates + modify substrates
  sublb_tmp <- matrix(0,nrow=nrow(object@orgdat),ncol=(length(object@mediac)))
  sublb <- as.data.frame(sublb) #convert to data.frame for faster processing in apply
  for(j in seq_along(object@media)){ #get information from sublb matrix to media list
    sval <- sum(sublb[,object@media[[j]]@id])/nrow(sublb)
    submat <- matrix(sval,object@n,object@m)
    eval.parent(substitute(object@media[[j]]@diffmat <- Matrix::Matrix(submat, sparse=TRUE)))
    sublb_tmp[,j] <- apply(object@orgdat, 1, function(x,sub){return(sub[x[4],x[5]])},sub=submat)
  }
  sublb <- cbind(as.matrix(object@orgdat[,c(4,5)]),sublb_tmp)
  colnames(sublb) <- c('x','y',object@mediac)
  return(sublb)
})

#' @title Function for transforming the organism data frame to a presence/absence matrix of organisms
#'
#' @description The generic function \code{dat2mat} simulates the event of mixing all substrates and organisms in the environment.
#' @export
#' @rdname dat2mat
#'
#' @param object An object of class Arena.
#' @return Returns the presence/absence matrix of organisms on the grid based on the \code{orgdat} slot of the \code{Arena} class.
#' @seealso \code{\link{Arena-class}} and \code{\link{getSublb}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' occmat <- dat2mat(arena)
#' image(occmat)
setGeneric("dat2mat", function(object){standardGeneric("dat2mat")})
#' @export
#' @rdname dat2mat
setMethod("dat2mat", "Arena", function(object){
  newoccmat <- matrix(0,object@n,object@m)
  for(i in 1:nrow(object@orgdat)){
    newoccmat[object@orgdat[i,'y'],object@orgdat[i,'x']] = object@orgdat[i,'type']
  }
  return(newoccmat)
})

#' @title Function for searching a keyword in arena organisms and media
#'
#' @description The generic function \code{findInArena} tries to find information (e.g. full names) about a specific keyword
#' @export
#' @rdname findInArena
#'
#' @param object An object of class Arena.
#' @param pattern A pattern for searching
#' @param search_rea Only search for reactions
#' @param search_sub Only search for substances 
#' @examples
#' data(Ec_core)
#' bac <- Bac(Ec_core)
#' arena <- Arena(n=20,m=20)
#' arena <- addOrg(arena,bac,amount=10)
#' findInArena(arena, "acetate")
setGeneric("findInArena", function(object, pattern, search_rea=TRUE, search_sub=TRUE){standardGeneric("findInArena")})
#' @export
#' @rdname findInArena
setMethod("findInArena", "Arena", function(object, pattern, search_rea=TRUE, search_sub=TRUE){
  if(search_sub){
    res_id <- grep(x=object@mediac, pattern=pattern, ignore.case = TRUE)
    if(length(res_id)>0) print(object@mediac[res_id])
    res_name <- grep(x=names(object@mediac), pattern=pattern, ignore.case = TRUE)
    if(length(res_name)>0) print(object@mediac[res_name])
  }
  
  if(search_rea & length(object@models)>0){
    for(i in 1:length(object@models)){
      model = object@models[[i]]
      cat(paste0("\n\n", i, ". ", model@mod_desc, model@mod_name, "\n"))
      
      res_met_id   <- grep(x=model@met_id,   pattern=pattern, ignore.case = TRUE)
      if(length(res_met_id)>0) print(paste(model@met_id[res_met_id], model@met_name[res_met_id]))
      res_met_name   <- grep(x=model@met_name,   pattern=pattern, ignore.case = TRUE)
      if(length(res_met_name)>0) print(paste(model@met_id[res_met_name], model@met_name[res_met_name]))
      
      res_rea_id   <- grep(x=model@react_id,   pattern=pattern, ignore.case = TRUE)
      if(length(res_rea_id)>0) {
        print(paste(model@react_id[res_rea_id], model@react_name[res_rea_id]))
        sybil::printReaction(model, react=model@react_id[res_rea_id])}
      
      res_rea_name <- grep(x=model@react_name, pattern=pattern, ignore.case = TRUE)
      if(length(res_rea_name)>0) {
        print(paste(model@react_id[res_rea_name], model@react_name[res_rea_name]))
        sybil::printReaction(model, react=model@react_id[res_rea_name])}
    }
  }
})




#show function for class Arena
setMethod(show, "Arena", function(object){
  #
  # 1) goup substances according to concentrations
  all_conc<-lapply(object@media, function(m){
    sum(m@diffmat)/length(m@diffmat)
  })
  group_conc <- split(all_conc, factor(unlist(unname(all_conc))))
  lapply(seq_along(group_conc), function(i){
    if(as.numeric(names(group_conc[i])) != 0){ # ignore substances with zero value
      print(paste("substances with", names(group_conc)[i], "fmol per gridcell:"))
      print(names(group_conc[[i]]))
      cat("\n")
    }
  })
  #
  # 2) general arena info
  print(paste("arena grid cells:",object@n,"x",object@m))
  print(paste("arena grid size [cm]:",object@Lx,"x",object@Ly))
  print(paste("dimension of one grid cell [cm]:",object@Lx/object@n,"x",object@Ly/object@m))
  print(paste("area of one grid cell [cm^2]:", (object@Lx*object@Ly)/(object@n*object@m)))  
  print(paste("flux unit:","mmol/(h*g_dw)"))
  print(paste("1 mM in arena correspons to mmol/grid_cell:", 1/100 * (object@Lx*object@Ly)/(object@n*object@m) ))
  #print(paste("1mM in arena correspons to mmol/grid_cell:", 1/1000 * (object@Lx*object@Ly*sqrt(object@Lx*object@Ly))/(object@n*object@m) ))
  print(paste('Arena of size ',object@n,' x ',object@m,' with ',nrow(object@orgdat),
              ' organisms of ',length(object@specs),' species.',sep=''))
})




# Eval is a subclass of Arena containing function to reduce the size of simulations and evalution of results

########################################################################################################
###################################### EVAL CLASS ######################################################
########################################################################################################

#' Structure of the S4 class "Eval"
#' 
#' Structure of the S4 class \code{Eval} inheriting from class \code{\link{Arena-class}} for the analysis of simulations.
#' @export Eval
#' @exportClass Eval
#' @importFrom  graphics barplot legend lines par points segments text
#' @importFrom stats cor dist hclust na.omit prcomp rnorm
#' @importFrom utils combn
#'
#' @slot medlist A list of compressed medium concentrations (only changes of concentrations are stored) per time step.
#' @slot simlist A list of the organism features per time step.
#' @slot mfluxlist A list of containing highly used metabolic reactions per time step. 
#' @slot shadowlist A list of containing shadow prices per time step. 
#' @slot subchange A vector of all substrates with numbers indicating the degree of change in the overall simulation.
#' @slot exchangeslist A list of containing exchanges per time step. 
setClass("Eval",
         contains="Arena",
         representation(
           medlist="list",
           simlist="list",
           mfluxlist="list",
           shadowlist="list",
           exchangeslist="list",
           subchange="numeric"
         ),
         prototype(
           medlist = list(),
           simlist = list(),
           mfluxlist = list(),
           shadowlist = list(),
           exchangeslist = list()
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

#' Constructor of the S4 class \code{\link{Eval-class}}
#' 
#' @name Eval-constructor
#' @export
#' @param arena An object of class Arena.
Eval <- function(arena){
  subc = rep(0, length(arena@mediac))
  names(subc) <- arena@mediac
  new("Eval", n=arena@n, m=arena@m, Lx=arena@Lx, Ly=arena@Ly, tstep=arena@tstep, specs=arena@specs, mediac=arena@mediac, subchange=subc,
      phenotypes=arena@phenotypes, media=arena@media, orgdat=arena@orgdat, medlist=list(), simlist=list(), stir=arena@stir, mfluxlist=list(), shadowlist=list(), exchangeslist=list() )
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("medlist", function(object){standardGeneric("medlist")})
setMethod("medlist", "Eval", function(object){return(object@medlist)})
setGeneric("simlist", function(object){standardGeneric("simlist")})
setMethod("simlist", "Eval", function(object){return(object@simlist)})
setGeneric("mfluxlist", function(object){standardGeneric("mfluxlist")})
setMethod("mfluxlist", "Eval", function(object){return(object@mfluxlist)})
setGeneric("shadowlist", function(object){standardGeneric("shadowlist")})
setMethod("shadowlist", "Eval", function(object){return(object@shadowlist)})
setGeneric("subchange", function(object){standardGeneric("subchange")})
setMethod("subchange", "Eval", function(object){return(object@subchange)})
setGeneric("exchangeslist", function(object){standardGeneric("exchangeslist")})
setMethod("exchangeslist", "Eval", function(object){return(object@exchangeslist)})


########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#' @title Function for adding a simulation step
#'
#' @description The generic function \code{addEval} adds results of a simulation step to an \code{Eval} object.
#' @export
#' @rdname addEval
#'
#' @param object An object of class Eval.
#' @param arena An object of class Arena.
#' @param replace A boolean variable indicating if the last simulation step should be replaced by the new simulation step \code{arena}.
#' @details The function \code{addEval} can be used in iterations to manipulate an \code{Arena} object and store the results in an \code{Eval} object.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' addEval(eval,arena)
setGeneric("addEval", function(object, arena, replace=F){standardGeneric("addEval")})
#' @export
#' @rdname addEval
setMethod("addEval", "Eval", function(object, arena, replace=F){
  if(!replace){
    subch = rep(0, length(arena@mediac))
    if(length(object@medlist)!=0){
      names(subch) <- arena@mediac
      sapply(names(subch), function(x, oldmed, newmed){
        subch[x] <<- subch[x]+sum(abs(oldmed[[x]]-as.vector(newmed[[x]]@diffmat)))
      },oldmed=extractMed(object), newmed=arena@media)
      subch[names(object@subchange)] <- object@subchange + subch[names(object@subchange)]
      eval.parent(substitute(object@subchange <- subch))
    }
    if(sum(subch)!=0){
      eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x, subc){
        if(subc[x@id]!=0){
          return(as.vector(x@diffmat))
        }else{return(vector())}
      }, subc=subch)))
    }else{
      eval.parent(substitute(object@medlist[[length(object@medlist)+1]] <- lapply(arena@media, function(x){
        return(as.vector(x@diffmat))
      })))
    }
    eval.parent(substitute(object@simlist[[length(object@simlist)+1]] <- arena@orgdat))
    eval.parent(substitute(object@mfluxlist[[length(object@mfluxlist)+1]] <- arena@mflux))
    eval.parent(substitute(object@shadowlist[[length(object@shadowlist)+1]] <- arena@shadow))
    eval.parent(substitute(object@exchangeslist[[length(object@exchangeslist)+1]] <- arena@exchanges))
    eval.parent(substitute(object@phenotypes <- arena@phenotypes))
    eval.parent(substitute(object@specs <- arena@specs))
    eval.parent(substitute(object@mflux <- arena@mflux)) 
    eval.parent(substitute(object@mediac <- arena@mediac))
    eval.parent(substitute(object@media <- arena@media))
    eval.parent(substitute(object@seed <- arena@seed))
    eval.parent(substitute(object@occupyM <- arena@occupyM))
    eval.parent(substitute(object@gridgeometry <- arena@gridgeometry))
    eval.parent(substitute(object@models <- arena@models))
    eval.parent(substitute(object@scale <- arena@scale))
    eval.parent(substitute(object@sublb <- arena@sublb))
    eval.parent(substitute(object@exchanges <- arena@exchanges))
    
    
  }else{
    eval.parent(substitute(object@medlist[[length(object@medlist)]] <- lapply(arena@media, function(x){
      return(as.vector(x@diffmat))
    })))
    eval.parent(substitute(object@simlist[[length(object@simlist)]] <- arena@orgdat))
    eval.parent(substitute(object@mfluxlist[[length(object@mfluxlist)]] <- arena@mflux))
    eval.parent(substitute(object@shadowlist[[length(object@shadowlist)]] <- arena@shadow))
    eval.parent(substitute(object@exchangeslist[[length(object@exchangeslist)]] <- arena@exchanges))
    eval.parent(substitute(object@phenotypes <- arena@phenotypes))
    eval.parent(substitute(object@specs <- arena@specs)) 
    eval.parent(substitute(object@mflux <- arena@mflux)) 
    eval.parent(substitute(object@mediac <- arena@mediac))
    eval.parent(substitute(object@media <- arena@media))
    eval.parent(substitute(object@seed <- arena@seed))
    eval.parent(substitute(object@occupyM <- arena@occupyM))
    eval.parent(substitute(object@gridgeometry <- arena@gridgeometry))
    eval.parent(substitute(object@models <- arena@models))
    eval.parent(substitute(object@scale <- arena@scale))
    eval.parent(substitute(object@sublb <- arena@sublb))
    eval.parent(substitute(object@exchanges <- arena@exchanges))
    
  }
})

#' @title Function for re-constructing an Arena object from a simulation step
#'
#' @description The generic function \code{getArena} re-constructs an \code{Arena} object from a simulation step within an \code{Eval} object.
#' @export
#' @rdname getArena
#'
#' @param object An object of class Eval.
#' @param time A number giving the simulation step of interest.
#' @return Returns an object of class \code{Arena} containing the organisms and substance conditions in simulation step \code{time}.
#' @details The function \code{addEval} can be used to manipulate an \code{Arena} object from a simulation step to modify the subsequent simulation steps.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' arena5 <- getArena(eval,5)
setGeneric("getArena", function(object, time=(length(object@medlist)-1)){standardGeneric("getArena")})
#' @export
#' @rdname getArena
setMethod("getArena", "Eval", function(object, time=(length(object@medlist)-1)){ #index in R start at 1, but the first state is 0
  time = time+1 #index in R start at 1, but the first state is 0
  
  newmedia <- lapply(object@media[names(object@medlist[[time]])], function(x, meds, n, m){
    x@diffmat <- Matrix::Matrix(meds[[x@id]],nrow=m,ncol=n,sparse=TRUE)
    return(x)
  },meds=extractMed(object,time), n=object@n, m=object@m)
  occdat <- object@simlist[[time]]
  
  arena <- Arena(n=object@n, m=object@m, Lx=object@Lx, Ly=object@Ly, tstep=object@tstep, 
                 specs=object@specs, mediac=object@mediac, mflux=object@mfluxlist[[time]],
                 phenotypes=object@phenotypes , media=newmedia, orgdat=occdat, stir=object@stir, 
                 shadow=object@shadowlist[[time]], seed=object@seed, exchanges=object@exchangeslist[[time]])
  arena@occupyM <- object@occupyM
  # reinitialize lp objects
  for(i in seq_along(arena@specs)){ 
    algo <- ifelse(.hasSlot(arena@specs[[i]], "algo"), arena@specs[[i]]@algo, "fba")
    arena@specs[[i]]@lpobj <- sybil::sysBiolAlg(arena@specs[[i]]@model, algorithm=algo)}
  return(arena)
})

#' @title Function for reducing the size of an Eval object by collapsing the medium concentrations
#'
#' @description The generic function \code{redEval} reduces the object size of an \code{Eval} object.
#' @export
#' @rdname redEval
#'
#' @param object An object of class Eval.
#' @param time A number giving the simulation step of interest.
#' @return Returns an object of class \code{Arena} containing the organisms and substance conditions in simulation step \code{time}.
#' @details The function \code{redEval} can be used to reduce the size of an \code{Eval} object from a simulation step.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' eval_reduce <- redEval(eval,5)
setGeneric("redEval", function(object, time="all"){standardGeneric("redEval")})
#' @export
#' @rdname redEval
setMethod("redEval", "Eval", function(object, time=1:length(object@medlist)){ #index in R start at 1, but the first state is 0
  for(i in time){
    object@medlist[[i]] <- lapply(extractMed(object,i),sum)
  }
  return(object)
})

#' @title Function for re-constructing a medium concentrations from simulations
#'
#' @description The generic function \code{extractMed} re-constructs a list of vectors of medium concentrations from a simulation step in an \code{Eval} object.
#' @export
#' @rdname extractMed
#'
#' @param object An object of class Eval.
#' @param time A number giving the simulation step of interest.
#' @param mediac A character vector giving the names of substances, which should be added to the environment (the default takes all possible substances).
#' @return Returns a list containing concentration vectors of all medium substances.
#' @details Medium concentrations in slot \code{medlist} of an object of class \code{Eval} store only the changes of concentrations in the simulation process. The function \code{extractMed} reconstructs the original and uncompressed version of medium concentrations.
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' med5 <- extractMed(eval,5)
setGeneric("extractMed", function(object, time=length(object@medlist), mediac=object@mediac){standardGeneric("extractMed")})
#' @export
#' @rdname extractMed
setMethod("extractMed", "Eval", function(object, time=length(object@medlist), mediac=object@mediac){
  medl <- object@medlist
  medlind <- medl[[time]]
  for(i in which(names(medlind) %in% mediac)){
    if(length(medl[[time]][[i]])==0){
      j <- time
      while(length(medl[[j]][[i]])==0){j <- j-1}
      medlind[[i]] <- medl[[j]][[i]]
    }
  }
  return(medlind[which(names(medlind) %in% mediac)])
})

#' @title Function for plotting spatial and temporal change of populations and/or concentrations
#'
#' @description The generic function \code{evalArena} plots heatmaps from the simulation steps in an \code{Eval} object.
#' @export
#' @rdname evalArena
#'
#' @param object An object of class Eval.
#' @param plot_items A character vector giving the name of the items which should be plotted such as the population structure and several metabolites.
#' @param phencol A boolean variable indicating if the phenotypes of the organisms in the environment should be integrated as different colors in the population plot.
#' @param retdata A boolean variable indicating if the data used to generate the plots should be returned.
#' @param time A numeric vector giving the simulation steps which should be plotted.
#' @param show_legend A boolean variable indicating if a legend shuld be shown.
#' @param legend_pos Position of the legend.
#' @return Returns several plots of the chosen plot items. Optional the data to generate the original plots can be returned.
#' @details If \code{phencol} is \code{TRUE} then different phenotypes of the same organism are visualized by varying colors, otherwise different organism types are represented by varying colors. The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' evalArena(eval)
#'\dontrun{
#' ## if animation package is installed a movie of the simulation can be stored:
#' library(animation)
#' saveVideo({evalArena(eval)},video.name="Ecoli_sim.mp4")
#' }
setGeneric("evalArena", function(object, plot_items='Population', phencol=F, retdata=F, time=(seq_along(object@simlist)-1), show_legend=TRUE, legend_pos="left"){standardGeneric("evalArena")})
#' @export
#' @rdname evalArena
setMethod("evalArena", "Eval", function(object, plot_items='Population', phencol=F, retdata=F, time=(seq_along(object@simlist)-1), show_legend=TRUE, legend_pos="left"){ #index in R start at 1, but the first state is 0
  time = time+1
  #old.par <- par(no.readonly = TRUE)
  if(retdata){
    retlist = list()
    for(i in 1:length(plot_items)){
      retlist[[i]] = list()
    }
    names(retlist) = plot_items
  }
  for(i in time){
    subnam <- names(object@medlist[[i]])
    inds <- which(subnam %in% plot_items)
    meds <- extractMed(object, i)
    if(length(plot_items)==1){
      if(length(inds)!=0){
        for(j in 1:length(inds)){
          if(retdata){
            retlist[[subnam[inds[j]]]][[paste0("time",(i-1))]] = matrix(meds[[subnam[inds[j]]]],ncol=object@n,nrow=object@m)
          }
          image(matrix(meds[[subnam[inds[j]]]],ncol=object@n,nrow=object@m),axes=F,main=paste(subnam[inds[j]], ": #", i),
                zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
        }
      }
    }else if(length(plot_items)<=6){
      par(mfrow=c(2,ceiling(length(plot_items)/2)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[paste0("time",(i-1))]] = matrix(meds[[subnam[inds[j]]]],ncol=object@n,nrow=object@m)
        }
        image(matrix(meds[[subnam[inds[j]]]],ncol=object@n,nrow=object@m),axes=F,main=paste(subnam[inds[j]], ": #", i),
              zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
      }
    }else{
      par(mfrow=c(3,ceiling(length(plot_items)/3)))
      for(j in 1:length(inds)){
        if(retdata){
          retlist[[subnam[inds[j]]]][[paste0("time",(i-1))]] = matrix(meds[[subnam[inds[j]]]],ncol=object@n,nrow=object@m)
        }
        image(matrix(meds[[inds[j]]],ncol=object@n,nrow=object@m),axes=F,main=paste(subnam[inds[j]], ": #", i),
              zlim=c(0,max(unlist(lapply(object@medlist,function(x, snam){return(x[[snam]])},snam=subnam[inds[j]])))))
      }
    }
    if(plot_items[1]=='Population'){
      if(retdata){
        retlist[['Population']][[paste0("time",(i-1))]] = object@simlist[[i]]
      }
      if(phencol){
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=object@simlist[[i]]$type-1,axes=FALSE,cex=1,main=paste('Population', ": #", i), col=object@simlist[[i]]$phenotype+1)
        if(show_legend){
          df_legend <- unique(object@simlist[[i]][,c("type", "phenotype")])
          df_legend <- df_legend[order(df_legend$phenotype),]
          par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE) # extra space for legend
          legend("topright", inset=c(-0.4,0), legend=paste(df_legend$type, df_legend$phenotype), col=df_legend$phenotype+1, pch=df_legend$type-1)
        }
      }else{
        plot(object@simlist[[i]][,c('x','y')],xlim=c(0,object@n),ylim=c(0,object@m),xlab='',ylab='',
             pch=object@simlist[[i]]$type-1,axes=FALSE,cex=1,main=paste('Population', ": #", i), col=object@simlist[[i]]$type)
        if(show_legend){
          par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # extra space for legend
          legend(legend_pos, inset=c(-0.8,0),legend=names(object@specs), col=c(1:length(names(object@specs))), pch=c(1:length(names(object@specs)))-1)
        }
      }
    }
  }
  #title(main=paste("step",i), outer=TRUE)
  #par(old.par)
  if(retdata){
    return(retlist)
  }
})

#' @title Function for plotting the overall change as curves
#'
#' @description The generic function \code{plotCurves} plots the growth curves and concentration changes of substances from simulation steps in an \code{Eval} object.
#' @export
#' @rdname plotCurves
#'
#' @param object An object of class Eval.
#' @param medplot A character vector giving the name of substances which should be plotted.
#' @param retdata A boolean variable indicating if the data used to generate the plots should be returned.
#' @param remove A boolean variable indicating if substances, which don't change in their concentration should be removed from the plot.
#' @param legend Boolean variable indicating if legend should be plotted
#' @param graph True if graphic should be plotted.
#' @return Returns two graphs in one plot: the growth curves and the curves of concentration changes. Optional the data to generate the original plots can be returned.
#' @details The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' plotCurves(eval)
setGeneric("plotCurves", function(object, medplot=object@mediac, retdata=F, remove=F, legend=F, graph=T){standardGeneric("plotCurves")})
#' @export
#' @rdname plotCurves
setMethod("plotCurves", "Eval", function(object, medplot=object@mediac, retdata=F, remove=F, legend=F, graph=T){
  old.par <- par(no.readonly = TRUE)
  growths <- matrix(0, nrow=length(object@specs), ncol=length(object@simlist))
  rownames(growths) = names(object@specs)
  subs <- matrix(0, nrow=length(medplot), ncol=length(object@simlist))
  rownames(subs) = medplot
  for(i in seq_along(object@simlist)){
    simdat <- object@simlist[[i]]
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
  if(graph){
    par(mfrow=c(2,1))
    times = c((seq_along(object@simlist)-1)*object@tstep)
    plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(growths)),
         type='n', xlab='time in h', ylab='number of individuals on grid',
         main='Population')
    for(i in 1:nrow(growths)){
      lines(times, growths[i,], col=i, type='b', pch=i-1)
    }
    if(legend){legend('topleft',legend=rownames(growths),col=1:nrow(growths), cex=ifelse(length(object@specs)==1,1,0.4/log10(nrow(growths)+1)),pch=(0:nrow(growths)-1),lwd=1, bty="n")}
    plot(times, times, xlim=c(0,max(times)), ylim=c(0,max(subs)),
         type='n', xlab='time in h', ylab='concentration in mmol per gridcell',
         main='Substance concentrations')
    for(i in 1:nrow(subs)){
      lines(times, subs[i,], col=i)
    }
    par(old.par)
    if(legend){legend('right',legend=rownames(subs),col=1:nrow(subs),cex=0.4/log10(nrow(subs)+1),lwd=1)}
  }
  if(retdata){
    return(list('Population'=growths,'Substances'=subs))
  }
})

#' @title Function to get varying substances
#'
#' @description The generic function \code{getVarSubs} returns ordered list of substances that showed variance during simulation
#' @export
#' @rdname getVarSubs
#' 
#' @param object An object of class Eval.
#' @param show_products A boolean indicating if only products should be shown
#' @param show_substrates A boolean indicating if only substrates should be shown
#' @param size Maximal number of returned substances (default: show all)
#' @param cutoff Value used to define numeric accuracy while interpreting optimization results
setGeneric("getVarSubs", function(object, show_products=TRUE, show_substrates=TRUE, cutoff=1e-6, size=NULL){standardGeneric("getVarSubs")})
#' @export
#' @rdname getVarSubs
setMethod("getVarSubs", "Eval", function(object, show_products=FALSE, show_substrates=FALSE, cutoff=1e-6, size=NULL){
  prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
  list <- lapply(prelist, function(x){lapply(x, sum)})
  # attention round due to numeric accuracy
  mat <- round(matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist)), round(-log10(cutoff)))
  #mat <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
  mediac <- object@mediac
  #rownames(mat) <- gsub("\\(e\\)","", gsub("EX_","",mediac))
  rownames(mat) <- mediac
  mat_var  <- apply(mat, 1, stats::var)
  if(length(mat_var[which(mat_var>0)]) == 0) return() # no substance having a variance > 0
  if(!(show_products || show_substrates)) {
    ret <- sort(mat_var[which(mat_var>0)], decreasing=TRUE)
    len_ret <- length(ret)
    if(is.null(size) || size >= len_ret) n <- len_ret
    else n <- size
    return(ret[1:n])
  }
  #mat <- mat[which(mat_var>0),]
  rowMin <- apply(mat, 1, min)
  rowMax <- apply(mat, 1, max)
  mat_substrates <- mat_var[which(mat[,1] == rowMax & mat_var > 0)]
  mat_products   <- mat_var[which(mat[,1] == rowMin & mat_var > 0)]
  if( show_products) {
    mat_products <- sort(mat_products, decreasing=TRUE)
    len_mat_products = length(mat_products)
    if(is.null(size) || size >= len_mat_products) n <- len_mat_products
    else n <- size
    return(mat_products[1:n])
  }
  if( show_substrates ){
    mat_substrates <- sort(mat_substrates, decreasing=TRUE)
    len_mat_substrates <- length(mat_substrates)
    if(is.null(size) || size >= len_mat_substrates) n <- len_mat_substrates
    else n <- size
    return(mat_substrates[1:n])
  } 
  return()
})


#' @title Function to get timeline of a substance
#'
#' @description The generic function \code{getSubHist} returns list with amount of substance for each timestep
#' @export
#' @rdname getSubHist
#' 
#' @param object An object of class Eval.
#' @param sub Vector with substances.
#' @param unit Unit to be used
setGeneric("getSubHist", function(object, sub, unit="fmol/cell"){standardGeneric("getSubHist")})
#' @export
#' @rdname getSubHist
setMethod("getSubHist", "Eval", function(object, sub, unit="fmol/cell"){
  sub <- ifelse(!(sub %in% names(object@media)), paste0("EX_", sub, "(e)"), sub)
  in_arena <- sub %in% names(object@media)
  if(!all((in_arena))){
    warning(paste(sub[!in_arena], "does not exist in medium"))
  }
  conv <- 1/unit_conversion(object, unit)
  timeline <- sapply(seq_along(sub), function(i){
    s <- sub[i]
    if(in_arena[i]){
      unlist(lapply(seq_along(object@medlist), function(t){conv*sum(unlist(extractMed(object, time=t, mediac=s)))}))  
    }else{
      rep(NA, length(object@medlist))
    }
  })
  rownames(timeline) <- seq_along(object@medlist)
  return(timeline)
})



#' @title Function for plotting the overall change as curves with maximally distinct colors
#'
#' @description The generic function \code{plotCurves2} plots the growth curves and concentration changes of the most changing substances from simulation steps in an \code{Eval} object using maximally distinct colors.
#' @export
#' @rdname plotCurves2
#'
#' @param object An object of class Eval.
#' @param legendpos A character variable declaring the position of the legend
#' @param ignore A list of character variables with substance names that sould be omitted in the plot
#' @param num An integer defining the number of substrates to be plot
#' @param phencol Boolean variable indicating whether phenotypes should be higlighted
#' @param dict List defining new substance names. List entries are intepreted as old names and the list names as the new ones.
#' @param biomcol A boolean indicating if biomass should be included in gowth curve
#' @param subs List of substance names. If empty, substances with highest variance will be used.
#' @param growthCurve True if growth curve should be shown (default TRUE)
#' @param subCurve True if substance curve should be shown (default TRUE)
#' @return Returns two graphs in one plot: the growth curves and the curves of concentration changes
#' @details The parameter \code{retdata} can be used to access the data used for the returned plots to create own custom plots. 
#' @seealso \code{\link{Eval-class}} and \code{\link{Arena-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' plotCurves2(eval)
setGeneric("plotCurves2", function(object, legendpos="topleft", ignore=c("EX_h(e)","EX_pi(e)", "EX_h2o(e)"),
                                   num=10, phencol=FALSE, biomcol=FALSE, dict=NULL, subs=list(), growthCurve=TRUE, subCurve=TRUE){standardGeneric("plotCurves2")})
#' @export
#' @rdname plotCurves2
setMethod("plotCurves2", "Eval", function(object, legendpos="topright", ignore=c("EX_h(e)","EX_pi(e)", "EX_h2o(e)"), 
                                          num=10, phencol=FALSE, biomcol=FALSE, dict=NULL, subs=list(), growthCurve=TRUE, subCurve=TRUE){
  if(num>length(object@mediac) || num<1) stop("Number of substances invalid")
  
  # 1) print substance curve
  if(subCurve){
    # first get the correct (ie. complete) medlist
    prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
    list <- lapply(prelist, function(x){lapply(x, sum)})
    mat <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
    
    if(length(subs)==0){ # CASE1: plot most varying substances
      #remove substances that should be ignored
      ignore_subs <- which(object@mediac %in% ignore |  gsub("\\[e\\]","", gsub("\\(e\\)","", gsub("EX_","",object@mediac))) %in% ignore)
      if(length(ignore_subs) != 0){
        mat <- mat[-ignore_subs,]
        mediac <- object@mediac[-ignore_subs]
      } else mediac <- object@mediac
      rownames(mat) <-  gsub("\\[e\\]","", gsub("\\(e\\)","", gsub("EX_","",mediac)))
      mat_var  <- apply(mat, 1, stats::var)
      num_var <- length(which(mat_var>0))
      if(num_var>0){
        mat_nice <- tail(mat[order(mat_var),], ifelse(num_var>num, num, num_var))
      }else{
        print("All substances have variance of zero.")
        mat_nice <- tail(mat[order(mat_var),], num)
      }
    }else{ # CASE2: plot only substances given by subs
      subs_index <- which(object@mediac %in% subs | gsub("\\[e\\]","", gsub("\\(e\\)","", gsub("EX_","",object@mediac))) %in% subs)
      if(length(subs_index)==1) mat_nice <- matrix(mat[subs_index,], nrow=1) else  mat_nice <- mat[subs_index,]
      rownames(mat_nice) <- gsub("\\[e\\]","", gsub("\\(e\\)","", gsub("EX_","",object@mediac[subs_index])))
    }
    if(num>length(colpal3)) cols <- colpal1[1:num] else cols <- colpal3[1:num]
    graphics::matplot(t(mat_nice), type='l', col=cols, pch=1, lty=1, lwd=5,
            xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of substance in fmol',
            main='Strongly changing substances')
    if(length(dict) > 0){
      new_names = unlist(lapply(rownames(mat_nice), function(x){dict[[x]]}))
      legend(legendpos, new_names, col=cols, cex=0.7, fill=cols)
    } else legend(legendpos, rownames(mat_nice), col=cols, cex=0.7, fill=cols)
  }
  
  # 2) Plotting growth curve
  if(growthCurve){
    # get bacs
    list <- lapply(object@simlist, function(x){
      occ <- table(x$type)
      unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
    })
    mat_bac  <- do.call(cbind, list)
    rownames(mat_bac) <- names(object@specs)
    
    
    # biomass
    list <- lapply(object@simlist, function(x){
      sum(x$biomass)
    })
    mat_biom  <- do.call(cbind, list)
    rownames(mat_biom) <- "biomass [fg]"
    
    # get pheno
    if(phencol){
      pheno_nr <- table(names(object@phenotypes))
      list <- lapply(object@simlist, function(x){ # time step
        unlist(lapply(seq_along(object@specs), function(j){ # bac type
          occ <- table(x[which(x$type==j),]$phenotype)
          #p <- unlist(lapply(seq_along(object@phenotypes[[j]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
          p <- unlist(lapply(seq(0,pheno_nr[[names(object@specs[j])]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
          names(p) <- paste0(names(object@specs)[j], "_pheno", seq(0,pheno_nr[[names(object@specs[j])]]))
          p
        }))})
      mat_phen  <- do.call(cbind, list)
      if(biomcol) mat_with_phen <- rbind(mat_bac, mat_phen, mat_biom) else mat_with_phen <- rbind(mat_bac, mat_phen)
    } else{
      if(biomcol) mat_with_phen <- rbind(mat_bac, mat_biom) else mat_with_phen <- mat_bac
    }
  
    len <- dim(mat_with_phen)[1]
    if(len>length(colpal3)) cols <- colpal1[1:len] else cols <- colpal3[1:len]
    graphics::matplot(t(mat_with_phen), type='b', col=cols, pch=1, lty=1, lwd=5,
            xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of organisms',
            main='Growth curve')
    legend(legendpos, rownames(mat_with_phen), col=cols, cex=0.7, fill=cols)
  }
})


#' @title Function for plotting the overall change in reaction activity
#'
#' @description The generic function \code{plotTotFlux} plots the time course of reactions with high variation in activity for an \code{Eval} object.
#' @export
#' @rdname plotTotFlux
#'
#' @param object An object of class Eval.
#' @param legendpos A character variable declaring the position of the legend
#' @param num An integer defining the number of substrates to be plot
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' plotTotFlux(eval)
setGeneric("plotTotFlux", function(object, legendpos="topright", num=20){standardGeneric("plotTotFlux")})
#' @export
#' @rdname plotTotFlux
setMethod("plotTotFlux", "Eval", function(object, legendpos="topright", num=20){
  if(num<1) stop("Number of reactions invalid")
  list <- lapply(object@mfluxlist, function(x){
    unlist(x)
  })
  mat  <- do.call(cbind, list)
  mat_var  <- rowSums((mat - rowMeans(mat))^2)/(dim(mat)[2] - 1)
  mat_nice <- tail(mat[order(mat_var),], num)
  
  if(num>length(colpal3)) cols <- colpal1[1:num] else cols <- colpal3[1:num]
  graphics::matplot(t(mat_nice), type='l', col=cols, pch=1, lty=1, lwd=3,
          xlab='time in h', ylab='reaction activity in mmol/(h * g_DW)',
          main='Highly active reactions')
  legend(legendpos, rownames(mat_nice), col=cols, cex=0.6, fill=cols)
  return(rownames(mat_nice))
})


#' @title Function for getting a matrix of phenotypes from the dataset
#'
#' @description The generic function \code{getPhenoMat} reconstructs a matrix with the usage of exchange reactions of the different organisms in the environment.
#' @export
#' @rdname getPhenoMat
#'
#' @param object An object of class Eval.
#' @param time An integer indicating the time step to be used (default value is character "total")
#' @param sparse A boolean indicating whether zero entries should be removed from return matrix
#' @return Returns a matrix with different phenotypes of the organism as rows and all possible exchange reactions as columns. A value of 1 means secretion, 2 means uptake and 0 means no usage of the substance of interest.
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenotype}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' phenmat <- getPhenoMat(eval)
setGeneric("getPhenoMat", function(object, time="total", sparse=F){standardGeneric("getPhenoMat")})
#' @export
#' @rdname getPhenoMat
setMethod("getPhenoMat", "Eval", function(object, time="total", sparse=F){
  phenspec = list()
  for(i in names(object@specs)){
    phenspec[[i]] = object@phenotypes[which(names(object@phenotypes)==i)]
    names(phenspec[[i]]) = seq_along(phenspec[[i]])
  }
  phens = unlist(phenspec)
  if(time != "total"){
    time = time+1 #index in R start at 1, but the first state is 0
    tdat = object@simlist[[time]][,c('type','phenotype')]
    #tdat = tdat[-which(is.na(tdat$phenotype)),]
    tphen = factor(paste(names(object@specs)[tdat$type],tdat$phenotype,sep="."))
    pinds = which(names(phens) %in% levels(tphen))
    if(length(pinds)!=0){phens = phens[pinds]}
  }else{
    time = length(object@medlist)
  }
  phenmat <- matrix(0, nrow=length(phens), ncol=length(names(object@medlist[[time]])))
  colnames(phenmat) <- names(object@medlist[[time]])
  rownames(phenmat) <- names(phens)
  for(i in 1:nrow(phenmat)){
    phenmat[i,1:length(as.numeric(unlist(strsplit(phens[i],split={}))))] = as.numeric(unlist(strsplit(phens[i],split={})))
  }
  if(sparse){
    phenmat <- phenmat[,which(colSums(abs(phenmat))!=0)]
  }
  return(phenmat)
})

#' @title Function for mining/analyzing phenotypes which occured on the arena
#'
#' @description The generic function \code{minePheno} mines the similarity and differences of phenotypes reconstructed by \code{getPhenoMat} for each simulation step in an \code{Eval} object.
#' @export
#' @rdname minePheno
#'
#' @param object An object of class Eval.
#' @param plot_type A character vector giving the plot which should be returned (either "pca" for a principle coordinate analysis or "hclust" for hierarchical clustering).
#' @param legend Boolean variable indicating if legend should be plotted
#' @return Returns a plot for each simulation step representing the similarity of phenotypes of organisms within the environment. 
#' @param time An integer indicating the time step to be used (default value is character "total")
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenoMat}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' minePheno(eval)
setGeneric("minePheno", function(object, plot_type="pca", legend=F, time="total"){standardGeneric("minePheno")})
#' @export
#' @rdname minePheno
setMethod("minePheno", "Eval", function(object, plot_type="pca", legend=F, time="total"){
  phenmat <- getPhenoMat(object, time)
  splitr = strsplit(rownames(phenmat),'\\.')
  rownames(phenmat) = unlist(lapply(splitr, function(x){return(paste(x[1:(length(x)-1)],collapse='.'))}))
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

#' @title Function for selecting phenotypes which occured on the arena from specific iterations and species
#'
#' @description The generic function \code{selPheno} selects phenotypes from specific simulation step in an \code{Eval} object.
#' @export
#' @rdname selPheno
#'
#' @param object An object of class Eval.
#' @param time A numeric vector giving the simulation steps which should be plotted. 
#' @param type A names indicating the species of interest in the arena.
#' @param reduce A boolean variable indicating if the resulting matrix should be reduced.
#' @return Returns a matrix with the substrate usage and the number of individuals using the phenotype. 
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}} and \code{\link{getPhenoMat}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' selPheno(eval,time=5,type='ecoli_core_model',reduce=TRUE)
setGeneric("selPheno", function(object, time, type, reduce=F){standardGeneric("selPheno")})
#' @export
#' @rdname selPheno
setMethod("selPheno", "Eval", function(object, time, type, reduce=F){
  arena = getArena(object, time)
  type_num = which(names(arena@specs)==type)
  arena@orgdat[which(is.na(arena@orgdat$phenotype)),] = 0
  pabund = as.matrix(table(arena@orgdat[which(arena@orgdat$type == type_num),'phenotype']))
  rownames(pabund) = paste(type,'phen',rownames(pabund),sep='_')
  rownames(pabund)[which(rownames(pabund)==paste(type,'phen_0',sep='_'))] = 'inactive'
  colnames(pabund) = 'individuals'
  pmat = getPhenoMat(object, time)
  splitr = strsplit(rownames(pmat),'\\.')
  rownames(pmat) = unlist(lapply(splitr, function(x){return(paste(x[1:(length(x)-1)],collapse='.'))}))
  pmatsp = pmat[which(rownames(pmat) == type),]
  if(is.vector(pmatsp)){
    pmatsp = t(as.matrix(pmatsp))
  }else{
    pmatsp = as.matrix(pmatsp)
  }
  if(reduce){
    if(nrow(pmatsp)==1){
      pmatsp = t(as.matrix(pmatsp[,-which(apply(pmatsp,2,sum)==0)]))
    }else{
      pmatsp = pmatsp[,-which(apply(pmatsp,2,sum)==0)]
    }
  }
  if(length(grep('inactive',rownames(pabund)))!=0){
    rownames(pmatsp) = rownames(pabund)[-which(rownames(pabund)=="inactive")]
    pmatsp = as.data.frame(pmatsp)
    pmatsp['inactive',]=rep(0,ncol(pmatsp))
  }else{
    rownames(pmatsp) = rownames(pabund)
    pmatsp = as.data.frame(pmatsp)
  }
  pmatsp[,'individuals']=rep(NA,nrow(pmatsp))
  pmatsp[rownames(pabund),'individuals'] = pabund[,'individuals']
  pmatsp
  return(as.matrix(pmatsp))
})

#show function for class Eval
setMethod(show, signature(object="Eval"), function(object){
  cat('Evaluation results of ',length(object@medlist)-1,' simulation steps.\n')
  cat("\tarena grid cells:",object@n,"x",object@m,"\n")
  cat("\tarena grid size [cm]:",object@Lx,"x",object@Ly,"\n")
  cat("\tdimension of one grid cell [cm]:",object@Lx/object@n,"x",object@Ly/object@m,"\n")
  cat("\tarea of one grid cell [cm^2]:", (object@Lx*object@Ly)/(object@n*object@m),"\n")
  cat("\tflux unit:","mmol/(h*g_dw)","\n")
  cat("\t1 mM in arena correspons to mmol/grid_cell:", 1/100 * (object@Lx*object@Ly)/(object@n*object@m) ,"\n")
  cat('\tArena of size ',object@n,'x',object@m,' with at first ',nrow(object@orgdat),
              ' organisms of ',length(object@specs),' species.',"\n")
  
})



#' @title Function for investigating a specific phenotype of an organism
#'
#' @description The generic function \code{statPheno} provides statistical and visual information about a certain phenotype.
#' @export
#' @rdname statPheno
#'
#' @param object An object of class Eval.
#' @param type_nr A number indicating the Organism type of the phenotype to be investigated (from orgdat)
#' @param phenotype_nr A number indicating the phenotype to be investigated (from orgdat)
#' @param dict A character vector of all substance IDs with names that should be used instead of possibly cryptic IDs
#' @details The phenotypes are defined by flux through exchange reactions, which indicate potential differential substrate usages.
#' @seealso \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' statPheno(eval, type_nr=1, phenotype_nr=2)
setGeneric("statPheno", function(object, type_nr=1, phenotype_nr, dict=NULL){standardGeneric("statPheno")})
#' @export
#' @rdname statPheno
setMethod("statPheno", "Eval", function(object, type_nr=1, phenotype_nr, dict=NULL){
  spec <- names(object@specs[type_nr])
  all  <- object@phenotypes[ names(object@phenotypes) == spec ]
  phen <- all[phenotype_nr]
  mediac <- gsub("\\(e\\)","", gsub("EX_","",object@mediac))
  if(length(dict) > 0) mediac <- unlist(lapply(mediac, function(x){dict[[x]]}))
  
  v <- as.numeric(unlist(strsplit(phen, split={})))
  names(v) <- mediac
  
  cat(paste("\nproduced by", spec,"\n"))
  print(names(v[which(v==1)]))
  cat(paste("\nconsumed by", spec,"\n"))
  print(names(v[which(v==2)]))
  
  # substrates
  prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
  list <- lapply(prelist, function(x){lapply(x, sum)})
  mat_sub <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
  rownames(mat_sub) <- mediac
  
  occ <- unlist(lapply(seq_along(object@simlist), function(t){dim(object@simlist[[t]][which(object@simlist[[t]]$type==type_nr & object@simlist[[t]]$phenotype==phenotype_nr),])[1]}))
  if(sum(occ)==0){
    cat("\n occurence of phenotype over time\n")
    print(occ)
    stop("This phenotype has never lived?!")
  }
  t_lb <- min(which(occ>0))
  t_ub <- max(which(occ>0))
  cat(paste("\n", spec, "occured between time steps"))
  print(paste0(t_lb, "-", t_ub))
  
  # get substances whose dynamic correlate with investigated phenotype
  mat <- rbind(mat_sub, phenotype=occ)[,t_lb:t_ub]
  corr <- cor(t(mat))
  corr[is.na(corr)] <- 0
  sorted <- sort(corr[dim(mat)[1],])
  high_corr <- c(tail(sorted, n=5), head(sorted, n=5))
  cat("\nsubstance with highest correlation\n")
  print(high_corr)
  
  if( length(high_corr) > 0){
    par(mfrow=c(2,1))
    mat_nice <- mat_sub[names(high_corr)[names(high_corr)!="phenotype"],][,t_lb:t_ub]
    num <- dim(mat_nice)[1]
    if(num>length(colpal3)) cols <- colpal1[1:num] else cols <- colpal3[1:num]
    graphics::matplot(x=seq(t_lb, t_ub), y=t(mat_nice), type='l', col=cols, pch=1, lty=1, lwd=5,
            xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of substance in mmol',
            main='correlated substances')
    legend("topleft", rownames(mat_nice), col=cols, cex=0.4, fill=cols)
    
    mat_pheno <- mat["phenotype",]
    plot(y=mat_pheno, x=seq(t_lb, t_ub), type='l', col="black", pch=1, lty=1, lwd=5,
            xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of organisms',
            main=paste('Growth curve', spec, "phenotype", phenotype_nr))
    
    par(mfrow=c(1,1))
  }
  
})


#' @title Function for investigation of feeding between phenotypes
#'
#' @description The generic function \code{findFeeding} 
#' @export
#' @rdname findFeeding
#' @importFrom igraph graph.empty add.edges delete.edges delete.vertices V E degree vcount layout_with_fr
#' 
#' @param object An object of class Eval.
#' @param tcut Integer giving the minimal mutual occurence ot be considered (dismiss very seldom feedings)
#' @param scut substance names which should be ignored
#' @param legendpos A character variable declaring the position of the legend
#' @param dict List defining new substance names. List entries are intepreted as old names and the list names as the new ones.
#' @param lwd Line thickness scale in graph
#' @param org_dict A named list/vector with names that should replace (eg. unreadable) IDs
#' @return Graph (igraph)
#' 
setGeneric("findFeeding", function(object, dict=NULL, tcut=5, scut=NULL, org_dict=NULL, legendpos="topleft", lwd=1){standardGeneric("findFeeding")})
#' @export
#' @rdname findFeeding
setMethod("findFeeding", "Eval", function(object, dict=NULL, tcut=5, scut=NULL, org_dict=NULL, legendpos="topleft", lwd=1){
  
  # possible problem inactive phenotype is not mentioned in object@phenotypes...

  # 1) Time: get occupation matrix for all phenotypes (occ_phen)
  pheno_nr <- table(names(object@phenotypes))
  pheno_index <- vector()
  list <- lapply(object@simlist, function(x){ # time step
    unlist(lapply(seq_along(object@specs), function(j){ # bac type
      occ <- table(x[which(x$type==j),]$phenotype)
      p <- unlist(lapply(seq(0,pheno_nr[[names(object@specs[j])]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
      if(names(object@specs[j]) %in% org_dict){
        org_name <- names(org_dict[which(org_dict==names(object@specs[j]))])
      }else org_name <- names(object@specs)[j]
      names(p) <- paste0(org_name, "_", seq(0,pheno_nr[[names(object@specs[j])]]))
      p
    }))})
  
  mat_phen  <- do.call(cbind, list)
  mat_tmp <- mat_phen
  mat_tmp[mat_tmp>0] <- 1
  cutoff <- which(rowSums(mat_tmp)>tcut)
  if(length(cutoff) < 2) {
    stop("tcut too high, no phenotypes found")
  } else mat_phen <- mat_phen[which(rowSums(mat_tmp)>tcut),] # reduce considered phenotypes (should exists >= tcut time steps)
  occ_phen  <- mat_phen != 0

  # 2) Substances: get matrix of substrates that could consumed and produced by phenotypes in principle
  mediac <- gsub("\\(e\\)","", gsub("EX_","",object@mediac))
  if(length(dict) > 0) mediac <- unlist(lapply(mediac, function(x){dict[[x]]}))
  phens <- object@phenotypes
  phenmat <- matrix(0, nrow=length(phens), ncol=length(object@mediac))
  colnames(phenmat) <- mediac
  counter = vector("numeric", length(object@specs))
  names(counter) <- lapply(names(object@specs), function(org_name){
    if(org_name %in% org_dict){
      org_name <- names(org_dict[which(org_dict==org_name)])
    }
    org_name
  })
  new_names = unlist(lapply(names(phens), function(x){
    if(x %in% org_dict){
      x = names(org_dict[which(org_dict==x)])
    }
    counter[x] <<- counter[x] + 1
    paste0(x, "_", counter[x])
  }))
  rownames(phenmat) <- new_names
  for(i in 1:nrow(phenmat)){
    phenmat[i,] = as.numeric(unlist(strsplit(phens[i],split={})))
  }
  phenmat_bin <- replace(phenmat, phenmat==2, -1)
  phenmat_abs <- abs(phenmat_bin)
  res <- phenmat_bin[,which(abs(colSums(phenmat_bin)) != colSums(phenmat_abs))]
  if(!is.null(scut) && all(scut %in% object@mediac)) scut <- gsub("\\(e\\)","", gsub("EX_","",scut))
  if(length(intersect(scut, mediac)) > 0) {
    res <- res[,-which(colnames(res) %in% scut)] # reduce substrates
  }else if(!is.null(scut)) print("scut should have valid names (as defined in mediac)")
  
  # graph
  pindex <- rownames(mat_phen)# phenotype index
  cindex <- colnames(res) # substance color index
  g <-igraph::graph.empty(n=length(pindex), directed=TRUE)
  igraph::V(g)$name <- gsub("pheno","",pindex)
  

  # 3) Combinatorics: check for all pairs of phenotypes if they 
  #    i) occured in the same time steps and 
  #   ii) exchange at least one substance
  #combi <- combn(rownames(phenmat), 2)
  combi <- combn(rownames(mat_phen), 2)
  for(i in 1:ncol(combi)){
    if(length(grep("_0$", combi[,i])) > 0) next # do not consider cuples with phenotype0 (i.d. metabolic inactive)
    co_occ <- which(occ_phen[combi[,i][1],] & occ_phen[combi[,i][2],]==T)
    if( length(co_occ) > 0 ){
      ex_both     <- res[c(combi[,i][1], combi[,i][2]),]
      feeding_index <- which(colSums(ex_both)==0 & colSums(abs(ex_both))!=0)
      if(length(feeding_index)>0){
        feeding <- ex_both[,feeding_index]
        if(length(feeding_index)==1){ # if only one substance is exchanged some hack to get name of substance into returned data structure of feeding
          feeding <- as.matrix(ex_both[,feeding_index])
          colnames(feeding) <- colnames(ex_both)[feeding_index]}
        lapply(seq(dim(feeding)[2]), function(x){
          if(feeding[1,x] == -1){
            new_edge <- c(which(pindex==combi[,i][2]), which(pindex==combi[,i][1]))
          }else new_edge <- c(which(pindex==combi[,i][1]), which(pindex==combi[,i][2]))
          col <- colpal3[which(cindex == colnames(feeding)[x])]
          g <<- igraph::add.edges(g, new_edge, color=col, weight=length(co_occ)*lwd)
        })
        #cat("\npossible cross feeding at time steps\n")
        #print(co_occ)
        #print(feeding)
      }
    }
  }
  
  g <- igraph::delete.edges(g, which(igraph::E(g)$weight<tcut)) # delete seldom feedings
  g <- igraph::delete.vertices(g, which(igraph::degree(g, mode="all") == 0)) # delete unconnected
  
  if(igraph::vcount(g) >= 2){
    plot(g, layout=igraph::layout_with_fr, vertex.size=5,
         edge.arrow.size=0.3, edge.width=igraph::E(g)$weight/10)
    legend(legendpos,legend=cindex, col=colpal3, pch=19, cex=0.7)
  }
  return(list(g=g, cindex=cindex))
})


#' @title Function for investigation of feeding between phenotypes
#'
#' @description The generic function \code{findFeeding2} 
#' @export
#' @rdname findFeeding2
#' @importFrom igraph V E graph.data.frame layout.circle
#' 
#' @param object An object of class Eval.
#' @param time A numeric vector giving the simulation steps which should be plotted. 
#' @param mets Character vector of substance names which should be considered
#' @param rm_own A boolean flag indicating if interactions within same species should be plotted
#' @param ind_threshold A number indicating the threshold of individuals to be considered as producers/consumers
#' @param collapse A boolean flag indicating if all phenotypes for every species should be collapsed to either producers or consumers
#' @return Graph (igraph)
#' 
setGeneric("findFeeding2", function(object, time, mets, rm_own=T, ind_threshold=0, collapse=F){standardGeneric("findFeeding2")})
#' @export
#' @rdname findFeeding2
setMethod("findFeeding2", "Eval", function(object, time, mets, rm_own=T, ind_threshold=0, collapse=F){
  mets = intersect(object@mediac,as.character(mets))
  pmat = as.matrix(getPhenoMat(object,time=time)[,mets])
  colnames(pmat) = mets
  time = time+1
  flux = lapply(object@mfluxlist[[time]], function(x){return(x[mets])})
  inter = data.frame(sp1=factor(levels=names(object@specs)),sp2=factor(levels=names(object@specs)),met=factor(levels=mets),
                     producer=numeric(),consumer=numeric(),fluxsp1=numeric(),fluxsp2=numeric())
  for(i in names(object@specs)){
    for(j in names(object@specs)){
      for(k in mets){
        pabi = table(object@simlist[[time]][which(object@simlist[[time]]$type==which(names(object@specs)==i)),"phenotype"])
        pabj = table(object@simlist[[time]][which(object@simlist[[time]]$type==which(names(object@specs)==j)),"phenotype"])
        prod1 = sum(pabi[gsub(paste(i,".",sep=""),"",names(which(pmat[grep(i,rownames(pmat)),k]==1)))])
        cons1 = sum(pabi[gsub(paste(i,".",sep=""),"",names(which(pmat[grep(i,rownames(pmat)),k]==2)))])
        prod2 = sum(pabj[gsub(paste(j,".",sep=""),"",names(which(pmat[grep(j,rownames(pmat)),k]==1)))])
        cons2 = sum(pabj[gsub(paste(j,".",sep=""),"",names(which(pmat[grep(j,rownames(pmat)),k]==2)))])
        if(prod1 > 0){
          if(cons1 > 0){inter[nrow(inter)+1,c("sp1","sp2","met","producer","consumer","fluxsp1","fluxsp2")] = c(i,i,k,prod1,cons1,flux[[i]][k],flux[[i]][k])}
          if(cons2 > 0){inter[nrow(inter)+1,c("sp1","sp2","met","producer","consumer","fluxsp1","fluxsp2")] = c(i,j,k,prod1,cons2,flux[[i]][k],flux[[j]][k])}
        }
        if(prod2 > 0){
          if(cons1 > 0){inter[nrow(inter)+1,c("sp1","sp2","met","producer","consumer","fluxsp1","fluxsp2")] = c(j,i,k,prod2,cons1,flux[[j]][k],flux[[i]][k])}
          if(cons2 > 0){inter[nrow(inter)+1,c("sp1","sp2","met","producer","consumer","fluxsp1","fluxsp2")] = c(j,j,k,prod2,cons2,flux[[j]][k],flux[[j]][k])}
        }
      }
    }
  }
  inter = inter[!duplicated(paste(inter[,1],inter[,2],inter[,3],sep="_")),]
  test = which(paste(inter[,1],inter[,2],sep="_") %in% paste(names(object@specs),names(object@specs),sep="_"))
  if(rm_own && length(test)!=0){inter = inter[-test,]}
  if(collapse){
    for(i in unique(as.character(inter$met))){
      mind = which(inter$met==i)
      interm = inter[mind,]
      for(j in unique(c(as.character(interm$sp1),as.character(interm$sp2)))){
        prod = which(interm$sp1 == j)
        cons = which(interm$sp2 == j)
        if(length(prod) !=0 && length(cons) != 0){
          if(as.numeric(interm[prod[1],"fluxsp1"]) > 0){
            interm = interm[-cons,]
          }else{
            interm = interm[-prod,]
          }
        }
      }
      inter = rbind(inter,interm)
      inter = inter[-mind,]
    }
    rminter = unique(which(as.numeric(inter$fluxsp1)<=0),which(as.numeric(inter$fluxsp2)>0))
    if(length(rminter)!=0){inter = inter[-rminter,]}
  }
  inter$rel_prod = as.numeric(inter$producer)/as.vector(table(object@simlist[[time]]$type))[as.numeric(inter$sp1)]
  inter$rel_cons = as.numeric(inter$consumer)/as.vector(table(object@simlist[[time]]$type))[as.numeric(inter$sp2)]
  if(ind_threshold<1){
    rmtr <- unique(c(which(inter$rel_prod<ind_threshold),which(inter$rel_cons<ind_threshold)))
    if(length(rmtr)==nrow(inter)){stop("No significant crossfeeding detected (try to relax ind_threshold).")}
    if(length(rmtr)!=0){inter = inter[-rmtr,]}
  }else{stop("ind_threshold needs to be between 0 and 1.")}
  vertexatt = data.frame(name = names(object@specs),color=1:length(object@specs),weight=as.vector(table(object@simlist[[time]]$type)))
  g <- igraph::graph.data.frame(inter[,1:2], directed=TRUE, vertices=vertexatt)
  l <- igraph::layout.kamada.kawai(g)
  plot(g,vertex.size=vertexatt$weight/max(vertexatt$weight)*20,edge.color=grDevices::rainbow(length(levels(inter$met)))[as.numeric(inter$met)],
       edge.arrow.size=0.5,edge.width=(inter$rel_prod*inter$rel_cons)*5,vertex.color=vertexatt$color+1,layout=l)
  legend("bottomright",legend=levels(inter$met),col=grDevices::rainbow(length(levels(inter$met))), pch=19, cex=0.7)
  return(list(inter,g))
})

#' @title Function for investigation of feeding between phenotypes
#'
#' @description The generic function \code{findFeeding3} 
#' @export
#' @rdname findFeeding3
#' @importFrom igraph V E graph.data.frame layout.circle
#' 
#' @param object An object of class Eval.
#' @param time A numeric vector giving the simulation steps which should be plotted. 
#' @param mets Character vector of substance names which should be considered
#' @param plot Should the graph also be plotted?
#' @param cutoff Accuracy of crossfeeding interaction (minimal flux to be considered)
#' @return Graph (igraph)
#' 
setGeneric("findFeeding3", function(object, time, mets, plot=TRUE, cutoff=1e-6){standardGeneric("findFeeding3")})
#' @export
#' @rdname findFeeding3
setMethod("findFeeding3", "Eval", function(object, time, mets, plot=TRUE, cutoff=1e-6){
  mets = intersect(object@mediac,as.character(mets))
  time = time+1
  mflux = object@mfluxlist[[time]]
  mfluxmat = do.call(cbind,lapply(mflux,function(x){return(ifelse(is.na(x[mets]),0,x[mets]))}))
  rownames(mfluxmat) = mets
  inter = data.frame()
  for( i in seq_along(rownames(mfluxmat)) ){
    x = mfluxmat[i,]
    interact = matrix(0,ncol=2,nrow=1)
    for(j in names(which(x < -cutoff))){
      if(length(which(x > cutoff))!=0){interact = rbind(interact,cbind(names(which(x > cutoff)),j))}
    }
    if( !is.null(interact) & nrow(interact) > 1){
      interact = interact[-1,,drop=FALSE] # remove zero row
      flux <- sapply(1:nrow(interact), function(k){
        idx.flux <- match(interact[k,], colnames(mfluxmat))
        mfluxmat[i, idx.flux]
      })
      if("character" %in% class(interact)){interact = t(as.matrix(interact))}
      inter = rbind(inter,data.frame(prod=interact[,1],cons=interact[,2],met=rownames(mfluxmat)[i], sim_step=time-1, prod.flux=flux[1,], cons.flux=flux[2,]))
    }
  }
  if(any(dim(inter)==0)) {
    warning(paste("sim_step",(time-1),":","No crossfeeding found. Try other metaboites or time points."))
    #g <- igraph::make_empty_graph()
    return(inter)
  }
  inter$met <- factor(inter$met)
  if (plot) {
  g <- igraph::graph.data.frame(inter[,1:2], directed=TRUE)
  l <- igraph::layout.kamada.kawai(g)
  plot(g,edge.color=grDevices::rainbow(length(levels(inter$met)))[as.numeric(inter$met)],
       edge.label=round(apply(abs(inter[,c("prod.flux", "cons.flux")]),1,min),1),
       edge.width=3,edge.arrow.size=0.8,vertex.color=1:length(igraph::V(g)),layout=l)
  legend("bottomright",legend=levels(inter$met),col=grDevices::rainbow(length(levels(inter$met))), pch=19, cex=0.7)
  return(invisible(list(inter,g)))}
  else return(inter)
})
                          



setGeneric("statSpec", function(object, type_nr=1, dict=NULL,
                                legend_show=TRUE, legend_pos="center", legend_cex=0.75){standardGeneric("statSpec")})
setMethod("statSpec", "Eval", function(object, type_nr=1, dict=NULL, 
                                       legend_show=TRUE, legend_pos="center", legend_cex=0.75){
  if(type_nr <= 0 || type_nr > length(object@specs)){
    stop("Invalid type number, should be number indicating a species of arena@specs")
  }
  sname <- names(object@specs[type_nr])
  pheno_nr <- table(names(object@phenotypes))[[sname]]
  print(paste(sname, "with phenotypes:", pheno_nr))
  
  # zero-phenotyp is inactive!
  phens <- rep(0, length(object@phenotypes)+1) # total number of occuring cells for each phenotyp
  ptimes <- rep(0, length(object@phenotypes)+1) # number of times steps a phenotyp occured
  availTimes <- rep(FALSE, length(object@simlist)) # true/false if cell is living for each timestep
  occ <- lapply(seq(from=2, to=length(object@simlist)), function(t){
    res <- table(object@simlist[[t]][which(object@simlist[[t]]$type==type_nr),]$phenotype)
    #print(phens)
    phens[as.numeric(names(res))+1]   <<- phens[as.numeric(names(res))+1] + res # -1 because of zero-phenotyp
    ptimes[as.numeric(names(res))+1]  <<- ptimes[as.numeric(names(res))+1] + 1
    if(length(res) > 0) availTimes[t] <<- TRUE
    p <- unlist(lapply(seq(0,pheno_nr), function(i){ifelse(i %in% names(res),res[paste(i)], 0)}))
    names(p) <- paste0("pheno", seq(0,pheno_nr))
    p
  })
  
  # species is available in the following time steps
  availTimesteps <- which(availTimes==TRUE)
  print(paste("alive between time steps:", min(availTimesteps), "-", max(availTimesteps)))
  
  len <- length(object@phenotypes)+1
  if(len>length(colpal3)) cols <- colpal1[1:len] else cols <- colpal3[1:len]

  # plot growth curve with all phenotypes
  mat_phen  <- do.call(cbind, occ)
  graphics::matplot(t(mat_phen), type='b', col=cols, pch=1, lty=1, lwd=5,
          xlab=paste0('time in ', ifelse(object@tstep==1, "", object@tstep), 'h'), ylab='amount of organisms',
          main='Growth curve')
  legend(legend_pos, rownames(mat_phen), col=cols, cex=legend_cex, fill=cols)  
  
  # analyze phenotypes: plot total cell number vs. time steps of occurence
  plot(ptimes, phens, bg=cols,type="p", pch = 23, log="y",
       xlab="Number of occured time steps ", ylab="Total number of cells", main=paste(sname,"phenotypes: cells vs. occurrence"))
  if(legend_show){
    legend(legend_pos, pch = 23, pt.bg=cols, cex=legend_cex,
           legend=seq(0, length(object@phenotypes)))
  }
})



#' @title Function to compute and return correlation matrix
#'
#' @description The generic function \code{getCorrM} returns the correlation matrix of several objects.
#' @export
#' @rdname getCorrM
#'
#' @param object An object of class Eval.
#' @param reactions A boolean indicating whether reactions should be included in correlation matrix
#' @param bacs A boolean indicating whether bacteria should be included in correlation matrix
#' @param substrates A boolean indicating whether substrates should be included in correlation matrix
#' @return correlation matrix
#' @details Returns correlation matrix which can be used for statistical analysis
#' @seealso \code{\link{Eval-class}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' getCorrM(eval)
setGeneric("getCorrM", function(object, reactions=TRUE, bacs=TRUE, substrates=TRUE){standardGeneric("getCorrM")})
#' @export
#' @rdname getCorrM
setMethod("getCorrM", "Eval", function(object, reactions=TRUE, bacs=TRUE, substrates=TRUE){
  mat <- matrix(0,0,length(object@medlist))
  if(substrates){
    prelist <- lapply(seq_along(object@medlist), function(i){extractMed(object, i)})
    list <- lapply(prelist, function(x){lapply(x, sum)})
    mat_sub <- matrix(unlist(list), nrow=length(object@media), ncol=length(object@medlist))
    rownames(mat_sub) <- gsub("\\(e\\)","", gsub("EX_","",object@mediac))
    mat <- rbind(mat, mat_sub)
  }
  
  if(reactions){
    list <- lapply(object@mfluxlist, function(x){
      unlist(x)
    })
    mat_rea  <- do.call(cbind, list)
    mat <- rbind(mat, mat_rea)
  }
  
  if(bacs){
    list <- lapply(object@simlist, function(x){
      occ <- table(x$type)
      unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
    })
    mat_bac  <- do.call(cbind, list)
    rownames(mat_bac) <- names(object@specs)
    mat <- rbind(mat, mat_bac)
  }
  
  corr <- cor(t(mat))
  corr[is.na(corr)] <- 0
  return(corr)
})
  
#' @title Function to show correlations of a simulated organism or substrate
#'
#' @description The generic function \code{checkCorr} returns the correlation matrix of several objects.
#' @export
#' @rdname checkCorr
#'
#' @param object An object of class Eval.
#' @param corr A correlation matrix (\code{\link{getCorrM}})
#' @param tocheck A list with substrate, reactions or organism names whose correlations should be shown
#' @details Returns correlation matrix which can be used for statistical analysis
#' @seealso \code{\link{Eval-class}} and \code{\link{getCorrM}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' checkCorr(eval, tocheck="o2")
setGeneric("checkCorr", function(object, corr=NULL, tocheck=list()){standardGeneric("checkCorr")})
#' @export
#' @rdname checkCorr
setMethod("checkCorr", "Eval", function(object, corr=NULL, tocheck=list()){
  if(is.null(corr)) corr <- getCorrM(object)

  lapply(tocheck, function(feature){
    dat <- corr[feature,]
    dat <- dat[-which(names(dat)==feature)]
    dat <- dat[order(dat)]
    dat_interest <- c(head(dat), tail(dat))
    barplot(names.arg=names(dat_interest), height=dat_interest, las=2, main=paste("Highest correlations of", feature))
  })
})

#' @title Function to plot substance shadow costs for a specie
#'
#' @description The generic function \code{plotShadowCost} plots substances have the highest impact on further growth (shadow cost < 0)
#' @export
#' @rdname plotShadowCost
#'
#' @param object An object of class Eval.
#' @param spec_nr Number of the specie
#' @param sub_nr Maximal number of substances to be show
#' @param cutoff Shadow costs should be smaller than cutoff
#' @param noplot Do not plot
#' @param useNames Use substance names instead of ids
#' @details Returns ggplot objects
setGeneric("plotShadowCost", function(object, spec_nr=1, sub_nr=10, cutoff=-1, noplot=FALSE, useNames=FALSE){standardGeneric("plotShadowCost")})
#' @export
#' @rdname plotShadowCost
setMethod("plotShadowCost", "Eval", function(object, spec_nr=1, sub_nr=10, cutoff=-1, noplot=FALSE, useNames=FALSE){
  
  if( length(object@shadowlist[[2]][[spec_nr]]) == 0){
    stop("No shadow costs were calculated during simulation. Try to enable it with simEnv(..., with_shadow=TRUE)")
  }
  
  df <- data.frame(spec=as.character(), sub=as.character(), shadow=as.numeric(), time=as.integer())
  
  m <- matrix(0, ncol=length(object@shadowlist[[1]][[spec_nr]]), nrow=length(object@shadowlist))
  for(t in seq_along(object@shadowlist)){
        m[t,] <- object@shadowlist[[t]][[spec_nr]]
  }
  if(all(m==0, na.rm = T)) {
    print("no shadow costs available")
    return()}
  rxn.last <- names(object@shadowlist[[1]][[spec_nr]])[length(object@shadowlist[[1]][[spec_nr]])]
  if(rxn.last == object@specs[[spec_nr]]@rbiomass & all(m[,-ncol(m)]==0, na.rm = T)){ # last column is reduced cost of biomass
    print("Growth seems to be only limited by duplication rate (check organism's maxweight)")
    return()}
  df <- as.data.frame(m)
  colnames(df) <- names(object@shadowlist[[1]][[spec_nr]])
  
  variance <- apply(m,2,stats::var)
  sorted_var <- sort(variance, decreasing=T, index.return=T)
  
  df <- df[,sorted_var$ix[1:sub_nr], drop=FALSE]
  colmin <- apply(df, 2, min)
  if(sum(colmin<cutoff)==0) stop("No shadow costs for substances found. Try different cutoff.")
  df <- df[,which(colmin<cutoff), drop=FALSE]
  df$time=seq_along(object@shadowlist)
  df <- reshape2::melt(df, id.vars="time")
  colnames(df)[2:3] <- c("sub", "shadow")
  
  if( useNames ) df$sub <- object@specs[[spec_nr]]@model@met_name[match(df$sub, object@specs[[spec_nr]]@model@met_id)]
  
  # do not plot shadow costs but print statistics
  if(noplot){
    df2 <- data.frame()
    for(s in unique(df$sub)){
      df2 <- rbind(df2, summary(df[which(df$sub==s),]$shadow))}
    rownames(df2) <- unique(df$sub)
    colnames(df2) <- names(summary(df[which(df$sub==s),]$shadow))
    return(df2)
  }
  
  q1 <- ggplot2::ggplot(df, ggplot2::aes(x=df$time, y=df$shadow)) + ggplot2::geom_line(ggplot2::aes(col=df$sub), size=1) + ggplot2::xlab("")
  
  q2 <- ggplot2::ggplot(df, ggplot2::aes(factor(df$sub), df$shadow)) + ggplot2::geom_boxplot(ggplot2::aes(color=factor(df$sub), fill=factor(df$sub)), alpha=0.2) +  ggplot2::ggtitle(names(object@specs)[spec_nr]) + ggplot2::xlab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  return(list(q1, q2))
})

#' @title Function to compute flux variability analysis on an simulation object to get min/max of substance usage
#'
#' @description The generic function \code{fluxVarSim} returns a list with the minimum and maximum substance usage of all individuals for each simulation step.
#' @export
#' @rdname fluxVarSim
#'
#' @param object An object of class Eval.
#' @param rnd An integer giving the decimal place to which min/max flux should be rounded.
#' @details Returns a list with the minimum and maximum substance usage for each time point.
#' @seealso \code{\link{Eval-class}} and \code{\link{simEnv}}
setGeneric("fluxVarSim", function(object, rnd){standardGeneric("fluxVarSim")})
#' @export
#' @rdname fluxVarSim
setMethod("fluxVarSim", "Eval", function(object, rnd){
  mflist = list()
  for(i in 1:length(object@simlist)){
    print(paste("Computing FVA for time point",i-1,"...",sep=" "))
    arenait = getArena(object, time=i-1)
    mflmat = matrix(0,nrow=length(arenait@mediac),ncol=2,
                    dimnames=list(arenait@mediac,c("min","max")))
    flmat_max = matrix(NA,nrow=length(arenait@mediac),ncol=nrow(arenait@orgdat))
    flmat_min = matrix(NA,nrow=length(arenait@mediac),ncol=nrow(arenait@orgdat))
    rownames(flmat_max) = arenait@mediac
    rownames(flmat_min) = arenait@mediac
    for(j in 1:nrow(arenait@orgdat)){
      org = arenait@orgdat[j,]
      bact = arenait@specs[[org$type]]
      mconc = unlist(lapply(arenait@media,function(med,x,y){med@diffmat[x,y]},x=org$x,y=org$y))
      mconc = mconc[bact@medium]
      bacnum = round((arenait@scale/(bact@cellarea*10^(-8))))
      lbs = constrain(bact,names(mconc),-mconc/bacnum,org$biomass,arenait@tstep,arenait@scale,j)[[1]]
      fbasl <- optimizeProb(bact@lpobj, react=1:length(lbs), ub=bact@ubnd, lb=lbs)
      model = changeBounds(bact@model,names(lbs),lb=lbs)
      model = changeBounds(model,model@react_id[which(bact@model@obj_coef==1)],lb=fbasl$obj,ub=fbasl$obj)
      nil=utils::capture.output(suppressMessages(fv <- fluxVar(model, bact@medium)))
      
      flmat_max[bact@medium,j] = round(minSol(fv,lp_obj),rnd)
      flmat_min[bact@medium,j] = round(maxSol(fv,lp_obj),rnd)
      #mflmat[bact@medium,"max"] = mflmat[bact@medium,"max"] + round(minSol(fv,lp_obj),rnd)
      #mflmat[bact@medium,"min"] = mflmat[bact@medium,"min"] + round(maxSol(fv,lp_obj),rnd)
    }
    mflmat[rownames(flmat_max),"max"] = apply(flmat_max,1,function(x){return(mean(x,na.rm=T))})
    mflmat[rownames(flmat_max),"min"] = apply(flmat_min,1,function(x){return(mean(x,na.rm=T))})
    mflist[[i]] = mflmat
  }
  return(mflist)
})

#' @title Function to get all reactions fluxes that are associated with the metabolite of a given exchange reactions
#'
#' @description The generic function \code{findRxnFlux} returns a matrix with the flux for each organism and the reaction that is using the metabolite of the given exchange reaction
#' @export
#' @rdname findRxnFlux
#'
#' @param object An object of class Eval.
#' @param ex An exchange reaction of which the metabolite should be shared for in all reactions
#' @param time the time point of the simulation which should be considered
#' @param print_reactions If true the detailed definition of each reactions is printed
#' @param drop_unused If true then inactive reactions will be excluded
#' @details Returns a list with the minimum and maximum substance usage for each time point.
#' @seealso \code{\link{Eval-class}} and \code{\link{simEnv}}
#' @examples
#' data(Ec_core, envir = environment()) #get Escherichia coli core metabolic model
#' bac <- Bac(Ec_core,deathrate=0.05,
#'            minweight=0.05,growtype="exponential") #initialize a bacterium
#' arena <- Arena(n=20,m=20) #initialize the environment
#' arena <- addOrg(arena,bac,amount=10) #add 10 organisms
#' arena <- addSubs(arena,40) #add all possible substances
#' eval <- simEnv(arena,5)
#' fluxlist <- findRxnFlux(eval, "EX_h(e)", 5)
setGeneric("findRxnFlux", function(object, ex, time, print_reactions=FALSE, drop_unused=TRUE){standardGeneric("findRxnFlux")})
#' @export
#' @rdname findRxnFlux
setMethod("findRxnFlux", "Eval", function(object, ex, time, print_reactions=FALSE, drop_unused=TRUE){
  sim=object
  aff_rxn = vector(); aff_rxn_prod = vector(); aff_rxn_cons = vector(); aff_rxn_names = vector()
  for(i in 1:length(sim@specs)){
    model = sim@specs[[i]]@model
    if(ex %in% model@react_id){
      reactind = which(model@react_id==ex) #find the index of the reaction
      metind = which(model@S[,reactind]==-1) #find which metabolites this reaction is affecting
      metcomp = gsub("\\[.\\]","",model@met_id) #strip off compartment id from metabolite
      mets = which(metcomp %in% metcomp[metind]) #find which metabolites are there 
      aff_rid = which(apply(abs(model@S[mets,]),2,sum)!=0) #get id of affected reactions
      aff_rid_prod = unlist(apply(model@S[mets,],1,function(r){which(r==1)})) #get id of  affected producing reactions
      aff_rid_cons = unlist(apply(model@S[mets,],1,function(r){which(r==-1)})) #get id of affected consuming reactions
      
      aff_rxn = union(aff_rxn,model@react_id[aff_rid]) #get the names of the afffected reactions
      aff_rxn_prod <- union(aff_rxn_prod, model@react_id[aff_rid_prod]) #get the names of the afffected producing reactions
      aff_rxn_cons <- union(aff_rxn_cons, model@react_id[aff_rid_cons]) #get the names of the afffected consuming reactions
      aff_rxn_names <-union(aff_rxn_names, sybil::printReaction(model, react=aff_rid, printOut=FALSE)) # save full reaction definitions
    } 
  }
  time=time+1
  mtflux = sim@mfluxlist[[time]] #select right timestep from fluxlist
  mtfmat = matrix(0,nrow=length(mtflux),ncol=length(aff_rxn),
                  dimnames=list(names(sim@specs),aff_rxn)) #create matrix to store flux of each species
  for(i in 1:length(mtflux)){ #iterate through each species
    flux <- ifelse(aff_rxn %in% names(mtflux[[i]]), mtflux[[i]][aff_rxn], 0)
    mtfmat[i,aff_rxn] = flux #store flux of key reactions per species in matrix}
  }
  if(print_reactions) cat(aff_rxn_names, "\n",  sep="\n") # print reaction strings if enabled
  # append index to reaction names to identify if a substance is produced or consumed by iy 
  colnames(mtfmat) <- ifelse(colnames(mtfmat) %in% aff_rxn_prod, paste0(colnames(mtfmat),"{R}"), paste0(colnames(mtfmat), "{L}"))
  cand_transp <- intersect(aff_rxn_prod, aff_rxn_cons)
  if(length(cand_transp)>0){ # identify possible transporters
    pos_transp <- which(aff_rxn %in% cand_transp)
    colnames(mtfmat)[pos_transp] <- paste0(aff_rxn[pos_transp],"{T}")}
  cat("Legend:\t {R} Substance is on the right site of the reaction definition (product)\n\t {L} Substance is on the left site (educt) {T} possible transporter reactions\n\n")
  
  if(drop_unused) mtfmat <- mtfmat[,which(colSums(abs(mtfmat))>0), drop=F]
  if(ncol(mtfmat)>0){
    sumflux = apply(abs(mtfmat),2,sum) #make the sum of the absolute flux for the reactions of each species
    return(mtfmat[,names(sort(sumflux,decreasing=T))]) #return the sorted matrix based on the absolute flux    
  } else print("No active reactions found.")
})


#' @title Function to overview the spatial distribution of a substance over time.
#'
#' @description The generic function \code{plotSubDist} returns a plot for every time step which shows the substance concentration in the environment.
#' @export
#' @rdname plotSubDist
#'
#' @param object An object of class Eval.
#' @param sub Name of a substance.
#' @param times Time points to be considered.
#' @details Returns a plot with 

setGeneric("plotSubDist", function(object, sub, times=NULL){standardGeneric("plotSubDist")})
#' @export
#' @rdname plotSubDist
setMethod("plotSubDist", "Eval", function(object, sub, times=NULL){
  if(length(sub) != 1 | !all(sub %in% object@mediac)) stop("Please use exactly one substance.")
  if(length(times)==0) times <- seq_along(object@medlist)
  outa <- t(sapply(times, function(i){c(i, unlist(extractMed(object, time=i, mediac=sub)))}))
  colnames(outa) <- c("time", paste(1:(object@n*object@m)))
  attributes(outa)$class <- c("deSolve", "matrix")
  attributes(outa)$dimens <- c(object@m,object@n)
  attributes(outa)$nspec <- 1
  mfrow <- sqrt(max(times))
  image(outa, ask = FALSE, mfrow = c(floor(mfrow), ceiling(mfrow)), main = paste(sub, times))
})
  

#' @title Function to overview the spatial distribution of a substance over time.
#'
#' @description The generic function \code{plotSubDist2} returns a plot for every time step which shows the substance concentration in the environment.
#' @export
#' @rdname plotSubDist2
#'
#' @param object An object of class Eval.
#' @param sub Name of a substance.
#' @param times Time points to be considered.
#' @details Returns a plot with 

setGeneric("plotSubDist2", function(object, sub, times=NULL){standardGeneric("plotSubDist2")})
#' @export
#' @rdname plotSubDist2
setMethod("plotSubDist2", "Eval", function(object, sub, times=NULL){
  if(length(sub) != 1 | !all(sub %in% object@mediac)) stop("Please use exactly one substance.")
  if(length(times)==0) times <- seq_along(object@medlist)
  if(max(times) > max(seq_along(object@medlist))) stop("Please use another maximum value in 'times' argument. Your input was out of bounds.")
  if(min(times) < min(seq_along(object@medlist))) stop("Please use another minimum value in 'times argument. Your input was out of bounds.")
  all_df <- data.frame()
  for(t in times){
    m <- matrix(unlist(extractMed(object, time=t, mediac=sub)), nrow = object@m, ncol=object@n)
    df <- reshape2::melt(m, varnames = c("x","y"))
    df$time = t
    all_df <- rbind(all_df, df)
  }
  q <- ggplot2::ggplot(all_df, ggplot2::aes_string("x", "y")) + ggplot2::geom_tile(ggplot2::aes_string(fill = "value")) + 
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") + 
    ggplot2::facet_wrap(~time, labeller = "label_both")+ ggplot2::theme_void() + ggplot2::ggtitle(sub) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(q)
})