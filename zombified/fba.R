library(libSBML)
library(lpSolveAPI)
doc <- readSBML("ecoli_core.xml")
ecoli <- SBMLDocument_getModel(doc)
Model_getNumSpecies(ecoli)
Model_getNumReactions(ecoli)

#
# make stoichiometric matrix
#


spec <- vector(mode = "character", length = 0)
for(i in seq_len(Model_getNumSpecies(ecoli))){  
  sp  =	Model_getSpecies(ecoli,	i-1);	
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

# objective function
c <- rep(0, dim(stoch)[2])

#c[which(colnames(stoch)=="R_Biomass_Ecoli_core_w_GAM")] <- 1 
c[which(colnames(stoch)=="R_Biomass_Ecoli_core_N__w_GAM_")] <- 1 
lb[which(colnames(stoch)=="R_ATPM")] <- 7.6
ub[which(colnames(stoch)=="R_ATPM")] <- 7.6

# define growth media
lb[grep("R_EX", colnames(stoch))] <- 0

lb[which(colnames(stoch)=="R_EX_glc_e_")] <- -10
lb[which(colnames(stoch)=="R_EX_h2o_e_")] <- -1000
lb[which(colnames(stoch)=="R_EX_h_e_")] <- -1000
lb[which(colnames(stoch)=="R_EX_o2_e_")] <- -1000
lb[which(colnames(stoch)=="R_EX_pi_e_")] <- -1000

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
get.objective(linp)

write.lp(linp,"test", NULL)

rm(linp)


