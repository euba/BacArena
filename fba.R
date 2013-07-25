library(libSBML)
library(lpSolveAPI)


#
# make stoichiometric matrix
#

read.sbml <- function(sbml_file){
  doc <- readSBML(sbml_file)
  ecoli <- SBMLDocument_getModel(doc)
  Model_getNumSpecies(ecoli)
  Model_getNumReactions(ecoli)  

  spec <- vector(mode = "character", length = 0)
  for(i in seq_len(Model_getNumSpecies(ecoli))){  
    sp  =  Model_getSpecies(ecoli,	i-1);	
    spec[i] <- Species_getId(sp);
  }
  ex <- rep(0, Model_getNumSpecies(ecoli)) 
  ex <- regexpr("_b$",spec) #mark external species/metabolites

  #spec <- as.factor(spec)
  reac <- vector(mode = "character", length = 0)
  for(i in seq_len(Model_getNumReactions(ecoli))){  
    re  =	Model_getReaction(ecoli,	i-1);	
    reac[i] <- Reaction_getId(re);
  }
  stoch <- matrix(data = 0, nrow = length(spec), ncol = length(reac))#, rownames=spec, colnames=reac)
  colnames(stoch) <- reac
  rownames(stoch) <- spec

  # initialize
  lb <- rep(-Inf, Model_getNumReactions(ecoli)) # lower bound
  ub <- rep(Inf, Model_getNumReactions(ecoli))  # upper bound
  
  #
  # fill stoichiometric matrix
  #
  #tmp<-matrix(numeric(0), dim(stoch)[1],0)
  for(i in seq_len(Model_getNumReactions(ecoli))){  
    re  =  Model_getReaction(ecoli,	i-1);	
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
  return(list(stoch=stoch, lb=lb, ub=ub, ex=ex, reac=reac))
}  

fba<-function(substrat, stoch, lb, ub, ex, reac, growth, sub_ex, type){

#Varma and Palsson 1994:
#maximum oxygen utilization rate (15 mmol of 02 per g [dry weight] per h)
#the maximum aerobic glucose utilization rate (10.5 mmol of Glc per g [dryweight] per h), the maximum anaerobic glucose utilization rate (18.5 mmol of Glc per g [dry weight] per h),  
#the non-growth-associated maintenance requirements (7.6 mmol of ATP per g [dry weight] per h), and the growth-associated maintenance requirements (13 mmol of ATP per g of biomass).
  
  
# objective function
c <- rep(0, dim(stoch)[2])
#c[which(colnames(stoch)=="R_Biomass_Ecoli_core_w_GAM")] <- 1 
c[which(colnames(stoch)==get_biomassf(type))] <- 1 

lb <- set_lower_bound(type, substrat, lb)

# linear programming
#linp <- make.lp(0, dim(stoch)[2], verbose = "full")
linp <- make.lp(0, dim(stoch)[2])
lp.control(linp,sense='max')
set.objfn(linp, c)
for(i in 1:nrow(stoch)){
  # only add constraint if it's not an external metabolite!!
  if(ex[i] == -1) add.constraint(linp, stoch[i,], "=", 0) 
}
set.bounds(linp,lower=lb)
set.bounds(linp,upper=ub)

solve(linp)
value <- get.objective(linp)

#get opt fluxes
flux <- get.variables(linp)
names(flux) <- reac

#considering the case that well growing bacs have higher fluxes!!
#print(substrat)
#print("")
#print(flux)
#print("")
#print(lb[which(colnames(stoch)=="R_EX_glc_e_")])
#print("")
pos_uptake <- sapply(names(sapply(substrat, names)),function(x,substrat,growth){
  return(flux[[sub_ex[[x]]]] * growth + substrat[[x]])
},substrat=substrat, growth=growth)
#print(pos_uptake)
#print("")
if(all(pos_uptake>=0) == T) flux <- flux * growth
#if(all(pos_uptake>=0)) print(flux * growth)


rm(linp)
return(flux)

#write.lp(linp,"test", NULL)

}

#tmp <- read.sbml("data/ecoli_core.xml")
#fba("substrat", tmp$stoch, tmp$lb, tmp$ub)
