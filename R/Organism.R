########################################################################################################
###################################### ORGANISM CLASS ##################################################
########################################################################################################

#' Structure of the S4 class "Organism"
#' 
#' Structure of the S4 class \code{Organism} representing the organisms present in the environment.
#'
#' @slot lbnd A numeric vector containing the lower bounds of the model structure.
#' @slot ubnd A numeric vector containing the upper bounds of the model structure.
#' @slot type A character vector containing the description of the organism.
#' @slot medium A character vector containing all exchange reactions of the organism.
#' @slot lpobj A sybil optimization object containing the linear programing problem.
#' @slot fbasol A list with the solutions of the flux balance analysis.
#' @slot lyse A boolean variable indicating if the organism should lyse after death.
#' @slot feat A list containing conditional features for the object (contains at the momement only biomass components for lysis).
#' @slot deathrate A numeric value giving the factor by which the growth should be reduced in every iteration.
#' @slot duplirate A numeric value giving the growth cut off at which the organism is duplicated.
#' @slot growthlimit A numeric value giving the growth limit at which the organism dies.
#' @slot growtype A character vector giving the functional type for growth (linear or exponential).
setClass("Organism",
         representation(
           lbnd="numeric",
           ubnd="numeric",
           type="character",
           medium="character",
           lpobj="sysBiolAlg",
           fbasol="list",
           lyse="logical",
           feat="list",
           deathrate="numeric",
           duplirate="numeric",
           growthlimit="numeric",
           growtype="character"
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Organism <- function(model, typename=mod_desc(model), algo="fba", ex=NA, deathrate, duplirate, growthlimit,
                     growtype="exponential", lyse=F, feat=list(), csuffix="\\[c\\]", esuffix="\\[e\\]", ...){ #the constructor requires the model, after that it is not stored anymore
  rxname = react_id(model)
  lpobject <- sysBiolAlg(model, algorithm=algo)
  fbasol <- optimizeProb(lpobject)
  names(fbasol$fluxes) = rxname
  lobnd = lowbnd(model)
  names(lobnd) = rxname
  upbnd = uppbnd(model)
  names(upbnd) = rxname
  if(lyse){
    stochmat <- as.matrix(S(model))
    colnames(stochmat) <- react_id(model)
    rownames(stochmat) <- met_id(model)
    stoch <- stochmat[,which(model@obj_coef==1)] #find stochiometry of biomass components
    biomets <- stoch[-which(stoch==0)]
    exs <- findExchReact(model)
    extrans <- react_id(exs)
    names(extrans) <- met_id(exs)
    names(biomets) <- gsub(csuffix,esuffix,names(biomets))
    biomets <- biomets[which(names(biomets) %in% names(extrans))]
    names(biomets) <- extrans[names(biomets)]
    feat[["biomass"]] <- biomets
  }
  if(is.na(ex)){
    medc <- react_id(findExchReact(ecore))
  }else{
    medc <- rxname[grep(ex, rxname)]
  }
  new("Organism", lbnd=lobnd, ubnd=upbnd, type=typename, medium=medc, lpobj=lpobject,
      fbasol=fbasol, lyse=lyse, feat=feat, deathrate=deathrate, duplirate=duplirate,
      growthlimit=growthlimit, growtype=growtype, ...)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("lbnd", function(object){standardGeneric("lbnd")})
setMethod("lbnd", "Organism", function(object){return(object@lbnd)})
setGeneric("ubnd", function(object){standardGeneric("ubnd")})
setMethod("ubnd", "Organism", function(object){return(object@ubnd)})
setGeneric("type", function(object){standardGeneric("type")})
setMethod("type", "Organism", function(object){return(object@type)})
setGeneric("medium", function(object){standardGeneric("medium")})
setMethod("medium", "Organism", function(object){return(object@medium)})
setGeneric("lpobj", function(object){standardGeneric("lpobj")})
setMethod("lpobj", "Organism", function(object){return(object@lpobj)})
setGeneric("fbasol", function(object){standardGeneric("fbasol")})
setMethod("fbasol", "Organism", function(object){return(object@fbasol)})
setGeneric("lyse", function(object){standardGeneric("lyse")})
setMethod("lyse", "Organism", function(object){return(object@lyse)})
setGeneric("feat", function(object){standardGeneric("feat")})
setMethod("feat", "Organism", function(object){return(object@feat)})
setGeneric("deathrate", function(object){standardGeneric("deathrate")})
setMethod("deathrate", "Organism", function(object){return(object@deathrate)})
setGeneric("duplirate", function(object){standardGeneric("duplirate")})
setMethod("duplirate", "Organism", function(object){return(object@duplirate)})
setGeneric("growthlimit", function(object){standardGeneric("growthlimit")})
setMethod("growthlimit", "Organism", function(object){return(object@growthlimit)})
setGeneric("growtype", function(object){standardGeneric("growtype")})
setMethod("growtype", "Organism", function(object){return(object@feat)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for constraining the models based on metabolite concentrations (can be given as vectors or single reaction)
#requires as input: organism object, reaction name, lowerbound, upperbound -> either lowerbound or upperbound can be omitted

setGeneric("constrain", function(object, reacts, lb, dryweight, time){standardGeneric("constrain")})
setMethod("constrain", "Organism", function(object, reacts, lb, dryweight, time){
  lobnd <- object@lbnd*dryweight*time #costrain according to flux definition: mmol/(gDW*hr)
  lobnd[reacts] <- ifelse(lb<=lobnd[reacts], lobnd[reacts], lb) #check if lower bounds in biological relevant range
  return(lobnd)
})

#function for computing the linear programming according to the model structure -> this can be changed for comp. speed (warmstart)
#requires as input: organism object -> changes the optimization object slot

setGeneric("optimizeLP", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd){standardGeneric("optimizeLP")})
setMethod("optimizeLP", "Organism", function(object, lpob=object@lpobj, lb=object@lbnd, ub=object@ubnd){ #this function has to be extended to contain also additional solvers
  switch(problem(lpob)@solver,
         glpkAPI=setColsBndsGLPK(problem(lpob)@oobj, 1:length(lb), #specific for GLPK!
                                    lb=lb, ub=ub),
         clpAPI=chgColLowerCLP(problem(lpob)@oobj, lb=lb))
  fbasl <- optimizeProb(lpob)
  names(fbasl$fluxes) <- names(object@lbnd)
  eval.parent(substitute(object@fbasol <- fbasl))
})

#function to account for the consumption and production of Substances

setGeneric("consume", function(object, sublb, cutoff=1e-6){standardGeneric("consume")})
setMethod("consume", "Organism", function(object, sublb, cutoff=1e-6){
  if(object@fbasol$obj>=cutoff){
    flux = object@fbasol$fluxes[object@medium]
    flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux))
    sublb[names(flux)] = round(sublb[names(flux)]+flux, 6)
  }
  return(sublb)
})

#function to extract the phenotype (what is consumed and what produced from the medium)

setGeneric("getPhenotype", function(object, cutoff=1e-6){standardGeneric("getPhenotype")})
setMethod("getPhenotype", "Organism", function(object, cutoff=1e-6){
  exflux=object@fbasol$fluxes[object@medium]
  exflux=ifelse(abs(exflux)<cutoff,0,1)*exflux
  exflux=ifelse(exflux<0,-1,exflux)
  exflux=ifelse(exflux>0,1,exflux)
  return(exflux[which(exflux!=0)])
})

#function for letting bacteria grow by adding the calculated growthrate to the already present growth value -> linear growth
#requires as input: organism object

setGeneric("growLin", function(object, growth){standardGeneric("growLin")})
setMethod("growLin", "Organism", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- object@fbasol$obj + growth
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})

#function for letting bacteria grow by adding the calculated growthrate multiplied with the current growth plus to the already present growth value -> exp growth
#requires as input: organism object

setGeneric("growExp", function(object, growth){standardGeneric("growExp")})
setMethod("growExp", "Organism", function(object, growth){
  if(object@fbasol$obj > 0) grow_accum <- (object@fbasol$obj * growth + growth)
  else grow_accum <- growth - object@deathrate
  return(grow_accum)
})

#function for lysis of bacterial cells by adding biomass_compounds * growth to the medium

setGeneric("lysis", function(object, sublb, factor=object@growthlimit){standardGeneric("lysis")})
setMethod("lysis", "Organism", function(object, sublb, factor=object@growthlimit){
  stoch = object@feat[["biomass"]]
  lysate = round(abs(stoch)*factor, 6)
  sublb[names(lysate)] = sublb[names(lysate)] + lysate
  return(sublb)
})

#function to get moore-neighbourhood of a bac together with its relative position

setGeneric("getHood", function(object, occmat, x, y){standardGeneric("getHood")})
setMethod("getHood", "Organism", function(object, occmat, x, y){
  occmat <- as.matrix(occmat) #dangerous!
  if(x-1==0) dx=0 else dx=1
  if(x+1>nrow(occmat)) dx2=0 else dx2=1
  if(y-1==0) dy=0 else dy=1
  if(y+1>ncol(occmat)) dy2=0 else dy2=1
  return(list(as.matrix(occmat[,(y-dy):(y+dy2)])[(x-dx):(x+dx2),], c(1+dx,1+dy)))
})


#function to check if the there is some free place in the neighbourhood

setGeneric("emptyHood", function(object, occmat, x, y){standardGeneric("emptyHood")})
setMethod("emptyHood", "Organism", function(object, occmat, x, y){
  hood <- getHood(object, occmat, x, y)
  free <- which(hood[[1]]==0, arr.ind = T)
  if(nrow(free) == 0) return(NULL)
  else {
    abs <- free[sample(nrow(free),1),]
    abs[1] <- abs[1] - hood[[2]][1] + x
    abs[2] <- abs[2] - hood[[2]][2] + y
    return(abs)
  }
})

#show function for class Organism

setMethod(show, signature(object="Organism"), function(object){
  print(paste('Organism ',object@type,' of class Organism.',sep=''))
})

# Bac is a subclass of Organism containing bacteria specific features

########################################################################################################
###################################### BAC CLASS #######################################################
########################################################################################################

#' Structure of the S4 class "Bac"
#' 
#' Structure of the S4 class \code{Bac} inheriting from class \code{\link{Organism}} representing bacterial cells.
#'
#' @slot speed A integer vector representing the speed by which bacterium is moving (given by cell per iteration).
#' @slot budge A boolean vector indicating if budging (bacteria in the souronding area are pushed away) should be implemented.
#' @slot chem A character vector indicating name of substance which is the chemotaxis attractant. Empty character vector if no chemotaxis.
setClass("Bac",
         contains="Organism",
         representation(
           speed="integer", # speed by which bacterium is moving (given by cell per iteration)
           budge="logical", #flag indicating, if budging (veruecktes Labyrinth) should be implemented
           chem="character" # name of substance which is the chemotaxis attractant
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Bac <- function(model, deathrate, duplirate, speed=2, growthlimit, growtype,
                budge=F, chem='', ...){
  new("Bac", Organism(model=model, deathrate=deathrate, duplirate=duplirate, growtype=growtype,
                      growthlimit=growthlimit, ...), budge=budge, speed=as.integer(speed), chem=chem)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("speed", function(object){standardGeneric("speed")})
setMethod("speed", "Bac", function(object){return(object@speed)})
setGeneric("budge", function(object){standardGeneric("budge")})
setMethod("budge", "Bac", function(object){return(object@budge)})
setGeneric("chem", function(object){standardGeneric("chem")})
setMethod("chem", "Bac", function(object){return(object@chem)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

# function with the growth model of a bac (biomass growth, replication, death)

#' Function with the growth model of a bac (biomass growth, replication, death)
#'
#' @param object An object of class Bac
#' @param population An object of class Arena
#' @param j The number of the iteration of interest
#' @return Boolean variable of the \code{j}th individual indicating if individual died. 
#' @examples
#' NULL
setGeneric("growth", function(object, population, j){standardGeneric("growth")})
setMethod("growth", "Bac", function(object, population, j){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$growth <- growLin(object, popvec$growth)},
         "exponential"= {popvec$growth <- growExp(object, popvec$growth)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'growth'] <- popvec$growth
  if(popvec$growth > object@duplirate){
    hood <- emptyHood(object, population@occmat, popvec$x, popvec$y)
    if(length(hood) != 0){
      daughter <- popvec
      daughter$growth <- popvec$growth/2
      daughter$x <- hood[1]
      daughter$y <- hood[2]
      popvec$growth = popvec$growth/2
      neworgdat[nrow(neworgdat)+1,] <- daughter
      neworgdat[j,] <- popvec
      eval.parent(substitute(population@occmat[daughter$x,daughter$y] <- as.numeric(daughter$type)))
    }else if(object@budge){
      hood2 = getHood(object, population@occmat, popvec$x, popvec$y)
      eval.parent(substitute(population <- budging(object, population, j, hood2, repli=T)))
    }
  }
  else if(popvec$growth < object@growthlimit){
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    neworgdat[j,'growth'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

# function for random movement

setGeneric("move", function(object, population, j){standardGeneric("move")})
setMethod("move", "Bac", function(object, population, j){
  popvec <- population@orgdat[j,]
  hood <- emptyHood(object, population@occmat, popvec$x, popvec$y)
  if(length(hood) != 0){
    xp = hood[1]
    yp = hood[2]
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    eval.parent(substitute(population@occmat[xp,yp] <- as.numeric(popvec$type)))
    eval.parent(substitute(population@orgdat[j,]$x <- xp))
    eval.parent(substitute(population@orgdat[j,]$y <- yp))
  }else if(object@budge){
    hood2 = getHood(object, population@occmat, popvec$x, popvec$y)
    eval.parent(substitute(population <- budging(object, population, j, hood2)))
  }
})

# function for chemotaxis: go to direction with highest concentration, otherwise random movement

setGeneric("chemotaxis", function(object, population, j){standardGeneric("chemotaxis")})
setMethod("chemotaxis", "Bac", function(object, population, j){
  popvec <- population@orgdat[j,]
  attract <- population@media[[object@chem]]@diffmat
  hood <- getHood(object, population@occmat, popvec$x, popvec$y)
  free <- which(hood[[1]]==0, arr.ind=T)
  if(nrow(free) != 0){
    conc <- apply(free, 1, function(x, attract, popvec, hood){
      xpos <- x[1] - hood[[2]][1] + popvec$x
      ypos <- x[2] - hood[[2]][2] + popvec$y
      return(attract[xpos,ypos])
    }, attract=attract, popvec=popvec, hood=hood)
    abs <- free[which(conc==max(conc)),]
    if(!is.vector(abs)){
      abs <- abs[sample(nrow(abs),1),]
    }
    abs[1] <- abs[1] - hood[[2]][1] + popvec$x
    abs[2] <- abs[2] - hood[[2]][2] + popvec$y
    hood <- abs
    xp = hood[1]
    yp = hood[2]
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    eval.parent(substitute(population@occmat[xp,yp] <- as.numeric(popvec$type)))
    eval.parent(substitute(population@orgdat[j,]$x <- xp))
    eval.parent(substitute(population@orgdat[j,]$y <- yp))
  }else if(object@budge){
    eval.parent(substitute(population <- budging(object, population, j, hood)))
  }
})

# function for budging of fellow bacteria, while one is moving

setGeneric("budging", function(object, population, j, hood, repli=F){standardGeneric("budging")})
setMethod("budging", "Bac", function(object, population, j, hood, repli=F){
  flag <- T
  orgdat <- population@orgdat
  orgxy <- paste(orgdat$x,orgdat$y,sep='_')
  hood[[1]][2,2] <- 0
  pos <- which(hood[[1]]!=0, arr.ind=T)
  inds <- sample(nrow(pos),nrow(pos))
  i <- 0
  while(i < nrow(pos)){
    i <- i+1
    xy <- pos[inds[i],]
    #test if all positions in this direction are blocked
    xp <- 1 #initialize xp and yp
    yp <- 1
    k <- j
    while(xp<population@n && yp<population@m){
      popvec <- orgdat[k,]
      xp <- xy[1] - hood[[2]][1] + popvec$x
      yp <- xy[2] - hood[[2]][2] + popvec$y
      pxy <- paste(xp,yp,sep='_')
      k <- which(orgxy==pxy)
      if(length(k)==0){flag <- F; break}
    }
  }
  while(!flag){
    popvec <- orgdat[j,]
    xp <- xy[1] - hood[[2]][1] + popvec$x
    yp <- xy[2] - hood[[2]][2] + popvec$y
    if(!repli){population@occmat[popvec$x, popvec$y] <- 0} #if not used in replication function
    population@occmat[xp,yp] <- as.numeric(popvec$type)
    population@orgdat[j,]$x <- xp
    population@orgdat[j,]$y <- yp
    pxy <- paste(xp,yp,sep='_')
    j <- which(orgxy==pxy)
    if(length(j)==0){flag <- T}
    repli <- T
  }
  return(population)
})

#function for one iteration for Bac class

setGeneric("simBac", function(object, arena, j, sublb){standardGeneric("simBac")})
setMethod("simBac", "Bac", function(object, arena, j, sublb){
  lobnd <- constrain(object, object@medium, lb=-sublb[j,object@medium],
                     dryweight=arena@orgdat[j,"growth"], time=arena@tstep)
  optimizeLP(object, lb=lobnd)
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,])))
  dead <- growth(object, arena, j)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object))
  
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,])))
  }
  if(!dead && !arena@stir && object@speed != 0){
    sapply(1:object@speed,function(x){
      if(object@chem == ''){
        move(object, arena, j)
      }else{
        chemotaxis(object, arena, j)
      }
      arena <<- arena})
  }
  return(arena)
})

#show function for class Bac

setMethod(show, signature(object="Bac"), function(object){
  print(paste('Bacterium ',object@type,' of class Bac.',sep=''))
})

# Human is a subclass of Organism containing human specific features

########################################################################################################
###################################### HUMAN CLASS #####################################################
########################################################################################################

#' Structure of the S4 class "Human"
#' 
#' Structure of the S4 class \code{Human} inheriting from class \code{\link{Organism}} representing human cells.
#'
#' @slot objective A character vector representing the current reaction which should be used as an objective function for the flux balance analysis.
setClass("Human",
         contains="Organism",
         representation(
           objective="character"
         )
)

########################################################################################################
###################################### CONSTRUCTOR #####################################################
########################################################################################################

Human <- function(model, deathrate, duplirate, growthlimit, growtype,
                  objective=model@react_id[which(model@obj_coef==1)], ...){
  model <- changeObjFunc(model, objective)
  new("Human", Organism(model=model, deathrate=deathrate, duplirate=duplirate, growtype=growtype,
                        growthlimit=growthlimit, ...), objective=objective)
}

########################################################################################################
###################################### GET METHODS FOR ATTRIBUTES ######################################
########################################################################################################

setGeneric("objective", function(object){standardGeneric("objective")})
setMethod("objective", "Human", function(object){return(object@objective)})

########################################################################################################
###################################### METHODS #########################################################
########################################################################################################

#function for changing the objective function of the model -> might be interesting for dynamic changes in varying environments
#requires as input: organism object -> changes the model, fobj slot and lpobj of the object

setGeneric("changeFobj", function(object, new_fobj, model, alg="fba"){standardGeneric("changeFobj")})
setMethod("changeFobj", "Human", function(object, new_fobj, model, alg="fba"){
  eval.parent(substitute(object@objective <- new_fobj)) #(pseudo) call by reference implementation
  model <- changeObjFunc(object@model, new_fobj)
  eval.parent(substitute(object@lpobj <- sysBiolAlg(model, algorithm=alg))) #the lp object has to be updated according to the new objective
})

# function with the growth model of a human cell (biomass growth, replication, death)

setGeneric("cellgrowth", function(object, population, j){standardGeneric("cellgrowth")})
setMethod("cellgrowth", "Human", function(object, population, j){
  neworgdat <- population@orgdat
  popvec <- neworgdat[j,]
  switch(object@growtype,
         "linear"= {popvec$growth <- growLin(object, popvec$growth)},
         "exponential"= {popvec$growth <- growExp(object, popvec$growth)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'growth'] <- popvec$growth
  if(popvec$growth > object@duplirate){
    hood <- emptyHood(object, population@occmat, popvec$x, popvec$y)
    if(length(hood) != 0){
      daughter <- popvec
      daughter$growth <- popvec$growth/2
      daughter$x <- hood[1]
      daughter$y <- hood[2]
      popvec$growth = popvec$growth/2
      neworgdat[nrow(neworgdat)+1,] <- daughter
      neworgdat[j,] <- popvec
      eval.parent(substitute(population@occmat[daughter$x,daughter$y] <- as.numeric(daughter$type)))
    }
  }
  else if(popvec$growth < object@growthlimit){
    eval.parent(substitute(population@occmat[popvec$x, popvec$y] <- 0))
    neworgdat[j,'growth'] <- NA
    dead <- T
  }
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
})

#function for one iteration for Human class

setGeneric("simHum", function(object, arena, j, sublb){standardGeneric("simHum")})
setMethod("simHum", "Human", function(object, arena, j, sublb){
  lobnd <- constrain(object, object@medium, lb=-sublb[j,object@medium],
                     dryweight=arena@orgdat[j,"growth"], time=arena@tstep)
  optimizeLP(object, lb=lobnd)
  eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,])))
  dead <- cellgrowth(object, arena, j)
  arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, object))
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, names(arena@media), sublb[j,])))
  }
  return(arena)
})

#show function for class Human

setMethod(show, signature(object="Human"), function(object){
  print(paste('Cell culture ',object@type,' of class Human.',sep=''))
})