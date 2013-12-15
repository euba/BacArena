library(libSBML)
library(lpSolveAPI)


#
# make stoichiometric matrix
#

read.sbml <- function(sbml_file, ex_pattern){
  doc <- readSBML(sbml_file)
  sbml <- SBMLDocument_getModel(doc)
  Model_getNumSpecies(sbml)
  Model_getNumReactions(sbml)  

  spec <- vector(mode = "character", length = 0)
  for(i in seq_len(Model_getNumSpecies(sbml))){  
    sp  =  Model_getSpecies(sbml,	i-1);	
    spec[i] <- Species_getId(sp);
  }
  ex <- rep(0, Model_getNumSpecies(sbml)) 
  ex <- regexpr(ex_pattern,spec) #mark external species/metabolites

  #spec <- as.factor(spec)
  reac <- vector(mode = "character", length = 0)
  for(i in seq_len(Model_getNumReactions(sbml))){  
    re  =	Model_getReaction(sbml,	i-1);	
    reac[i] <- Reaction_getId(re);
  }
  stoch <- matrix(data = 0, nrow = length(spec), ncol = length(reac))#, rownames=spec, colnames=reac)
  colnames(stoch) <- reac
  rownames(stoch) <- spec

  # initialize
  lb <- rep(-1e+30, Model_getNumReactions(sbml)) # lower bound
  ub <- rep(1e+30, Model_getNumReactions(sbml))  # upper bound
  
  #
  # fill stoichiometric matrix
  #
  #tmp<-matrix(numeric(0), dim(stoch)[1],0)
  for(i in seq_len(Model_getNumReactions(sbml))){  
    re  =  Model_getReaction(sbml,	i-1);	
    for(j in seq_len(Reaction_getNumReactants(re))){
      ed = Reaction_getReactant(re, j-1);
      stoch[SimpleSpeciesReference_getSpecies(ed),i] <- -SpeciesReference_getStoichiometry(ed)
    }
    for(j in seq_len(Reaction_getNumProducts(re))){
      pr = Reaction_getProduct(re, j-1);
      stoch[SimpleSpeciesReference_getSpecies(pr),i] <- stoch[SimpleSpeciesReference_getSpecies(pr),i] + SpeciesReference_getStoichiometry(pr)
    }
    #
    # consider reversible case?? (only lower bounds or add reversible reaction twice?)
    #
    #if(Reaction_getReversible(re)){
    #  tmp <- cbind(tmp, -stoch[,i])
    #colnames(stoch) <- c(colnames(stoch), paste(colnames(stoch[,i]),"rev", sep=""))
    #}
    if(!Reaction_getReversible(re)){
      lb[i] <- 0 
    }
  }
  #stoch <- cbind(stoch, tmp)
  str(stoch)
  return(list(stoch=stoch, lb=lb, ub=ub, ex=ex, reac=reac))
}  



fba<-function(substrat, stoch, ex, reac, growth, sub_ex, type){

#substrat <- lapply(substrat, function(x){round(x, digits=2)}) # ROUNDING!!!

# objective function
cs <- rep(0, dim(stoch)[2])
cs[which(colnames(stoch)==get_biomassf(type))] <- 1 

# set boundaries
lb <- set_lower_bound(type, substrat) # current substrate defines lower bounds
ub <- get_upper_bound(type)

#det boundaries ngam
lb[which(colnames(stoch)==get_maintenancef(type))] <- get_ngam(type)
ub[which(colnames(stoch)==get_maintenancef(type))] <- get_ngam(type)


# linear programming
#linp <- make.lp(0, dim(stoch)[2], verbose = "full")
linp <- make.lp(0, dim(stoch)[2])
#lp.control(linp,sense='max',verbose='full', anti.degen="none")
lp.control(linp,sense='max',verbose='important')
#row.add.mode(linp, "on") # performance improvement!
set.objfn(linp, cs)
for(i in 1:nrow(stoch)){
  # only add constraint if it's not an external metabolite!!
  if(ex[i] == -1) add.constraint(linp, stoch[i,], "=", 0) 
}
#row.add.mode(linp, "off") # performance improvement!
set.bounds(linp,lower=lb)
set.bounds(linp,upper=ub)

# print substrate
#print("substrate")
#print(t(substrat))
#print("")

# print lower bound
#print("lower bound")
# get lower bound with names
#lbound <- sapply(names(sapply(substrat, names)), function(x, stoch, sub_ex, lb){
#  if(x %in% names(sub_ex)) lb[which(colnames(stoch)==sub_ex[[x]])]
#},stoch=stoch, sub_ex=sub_ex, lb=lb)
#print(t(lbound))

#Solve opt problem
status <- solve(linp)
#print(status)
if(status!=0) return("DEAD") # very important!! if no feasible solutions is found there is still some (useless) result
value <- get.objective(linp)
#get opt fluxes
flux <- get.variables(linp)
names(flux) <- reac

# check for conistency (is fba return bounded correctly?)
sapply(names(sapply(substrat, names)),function(x,substrat){
  if(x %in% names(sub_ex)) { # only update substrate which are metabolic relevant for current organism
    #if(round(as.numeric(substrat[x]),2) < -round(flux[[sub_ex[[x]]]],2)){ # test for negative fba return
    # workaround because round is too slow!
    if(as.numeric(substrat[x]) + epsilon < -flux[[sub_ex[[x]]]]){ # test for negative fba return
      print(t(flux))
      print("")
      print("lower bound")
      # get lower bound with names
      lbound <- sapply(names(sapply(substrat, names)), function(x, stoch, sub_ex, lb){
        if(x %in% names(sub_ex)) lb[which(colnames(stoch)==sub_ex[[x]])]
      },stoch=sbml$stoch, sub_ex=sub_ex, lb=lb)
      print(t(lbound))
      print("")
      print("upper bound")
      rbound <- sapply(names(sapply(substrat, names)), function(x, stoch, sub_ex, ub){
        if(x %in% names(sub_ex)) ub[which(colnames(stoch)==sub_ex[[x]])]
      },stoch=sbml$stoch, sub_ex=sub_ex, ub=ub)
      print(t(rbound))
      print("")
      print(c(x, substrat[[x]], " uptake: ", -flux[[sub_ex[[x]]]]))
      stop("FBA ERROR: return flux excesses available substrate!")
    }
  }
},substrat=substrat)

#considering the case that well growing bacs have higher fluxes!!
pos_uptake <- sapply(names(sapply(substrat, names)),function(x,substrat,growth){
  if(x %in% names(sub_ex == T)) return(flux[[sub_ex[[x]]]] * growth + substrat[[x]])
  else return(0)
},substrat=substrat, growth=growth)
if(all(pos_uptake>=0) == T) flux <- flux * growth

# print uptake
# if(type=="barkeri"){
#   uptake <- sapply(names(sapply(substrat, names)),function(x,substrat){
#     if(x %in% names(sub_ex) == T) return(flux[[sub_ex[[x]]]])
#     else return(0)
#   },substrat=substrat)
#   names(uptake) <- names(substrat)
#   print("uptake")
#   print(t(uptake))
#   print("growth rate:")
#   print(flux[[get_biomassf(type)]])
# }

# debug file
write.lp(linp,"test", NULL) # debug 

rm(linp)
return(flux)

}


starvation_fees <-function(type){
  biomassf <- get_biomassf(type)
  maintenancef <- get_maintenancef(type)
  sbml <- get_sbml(type)
  stoch <- sbml$stoch
  ub <- get_upper_bound(type)
  lb <- get_lower_bound(type)
  sub_ex <- get_sub_ex(type)
  
  # add an input for all biomass metabolites
  iseq <- which (stoch[,biomassf] < 0)
  for (j in iseq){
    stoch <- cbind(stoch, c(rep(0, j-1), 1, rep(0, nrow(stoch) - j)))
    stoch <- rbind(stoch, c(rep(0, ncol(stoch)-1), -1))
    ub <- c(ub, -stoch[,biomassf][j])
    #lb <- c(lb, 0)
    lb <- c(lb, -1e+30)
  }
  
  # only biomassf as input/substrat
  for(x in sub_ex){
  #  lb[which(colnames(stoch)==x)] <- 0
    ub[which(colnames(stoch)==x)] <- 0
  #  print(x)
  }
  
  lb[which(colnames(sbml$stoch)==biomassf)] <- 0
  #ub[which(colnames(sbml$stoch)==biomassf)] <- 1
  
  #sbml$stoch[,biomassf]["M_atp_c"] <- 0
  max_ngam <- 0

  # complete successive a new biomassfunction (prevent accumulation of certain biomass components)
  #sbml$stoch[,biomassf] <- c(rep(0, j-1), biomassf_org[j], rep(0, nrow(sbml$stoch) - j))

  
  # objective function
  cs <- rep(0, dim(stoch)[2])
  cs[which(colnames(stoch)==maintenancef)] <- 1 
  
  # linear programming
  linp <- make.lp(0, dim(stoch)[2])
  lp.control(linp,sense='max',verbose='important')
  set.objfn(linp, cs)
  for(i in 1:nrow(sbml$stoch)){
    # only add constraint if it's not an external metabolite!!
    if (sbml$ex[i] == -1 ) add.constraint(linp, stoch[i,], "=", 0)
    #add.constraint(linp, stoch[i,], "=", 0)
    
    #if (i %in% iseq){
    #  add.constraint(linp, sbml$stoch[,i], ">=", 0)
    #  add.constraint(linp, sbml$stoch[,i], "<=", abs(sbml$stoch[,biomassf][i]))
    #} 
    #else add.constraint(linp, sbml$stoch[i,], "=", 0)
  }
  set.bounds(linp,lower=lb)
  set.bounds(linp,upper=ub)
  
  #Solve opt problem
  status <- solve(linp)
  print(status)
  #if(status==0) {
    value <- get.objective(linp)
    print(value)
    #get opt fluxes
    flux <- get.variables(linp)
    names(flux) <- sbml$reac
    max_ngam <- max_ngam + flux[[maintenancef]]
  #}
 
  #print(max_ngam)
  
  x <- stoch %*% flux
  names(x) <- rownames(stoch)
  #x[which(x < 0.1)]
  flux[which(flux > 0.1)]
  flux[which(flux < -0.1)]
  
  print(flux[[maintenancef]])
  #get.constraints(linp)
  
  write.lp(linp,"test", NULL) # debug 
  
  return (max_ngam)
}
