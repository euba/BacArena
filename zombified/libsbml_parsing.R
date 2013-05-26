library(libSBML)
doc <- readSBML("ecoli_orth2011_msb201165-s3.xml")
ecoli <- SBMLDocument_getModel(doc)
Model_getNumSpecies(ecoli)
Model_getNumReactions(ecoli)

r <- Model_getReaction(ecoli, 5)
Reaction_getId(r)

m <- Reaction_getListOfReactants(r)

ListOfSpeciesReferences_getElementBySId(m, 1)

ListOfSpeciesReferences_getAllElements(m)

# print all species
for(i in seq_len(Model_getNumSpecies(ecoli))){	
  sp	=	Model_getSpecies(ecoli,	i-1);	
  cat("species:	", Species_getId(sp),	"\n")	
}

# print out all reactions with reactants and stoichiometry
for(i in seq_len(Model_getNumReactions(ecoli))){  
  re	=	Model_getReaction(ecoli,	i-1);	
  cat("reactions:	", Reaction_getId(re),	"\n")	
  for(j in seq_len(Reaction_getNumReactants(re))){
    ed = Reaction_getReactant(re, j-1);
    #print(ed)
    cat("\t", SpeciesReference_getStoichiometry(ed), " ", SimpleSpeciesReference_getSpecies(ed))
  }
  cat(" --> ")
  for(j in seq_len(Reaction_getNumProducts(re))){
    pr = Reaction_getProduct(re, j-1);
    cat("\t", SpeciesReference_getStoichiometry(pr), " ", SimpleSpeciesReference_getSpecies(pr))
  }
  cat("\n")
}

re  =	Model_getReaction(ecoli,	2000);	
cat("reactions:  ", Reaction_getId(re),	"\n")	
for(j in seq_len(Reaction_getNumReactants(re))){
  ed = Reaction_getReactant(re, j-1);
  #print(ed)
  cat("\t", SpeciesReference_getStoichiometry(ed), " ", SimpleSpeciesReference_getSpecies(ed))
}
cat(" --> ")
for(j in seq_len(Reaction_getNumProducts(re))){
  pr = Reaction_getProduct(re, j-1);
  cat("\t", SpeciesReference_getStoichiometry(pr), " ", SimpleSpeciesReference_getSpecies(pr))
}
cat("\n")

#
# make stoichiometric matrix
#
spec <- vector(mode = "character", length = 0)
for(i in seq_len(Model_getNumSpecies(ecoli))){  
  sp	=	Model_getSpecies(ecoli,	i-1);	
  spec[i] <- Species_getId(sp);	
}
#spec <- as.factor(spec)
reac <- vector(mode = "character", length = 0)
for(i in seq_len(Model_getNumReactions(ecoli))){  
  re  =	Model_getReaction(ecoli,	i-1);	
  reac[i] <- Reaction_getId(re);
}
stoch <- matrix(data = 0, nrow = length(spec), ncol = length(reac))#, rownames=spec, colnames=reac)
colnames(stoch) <- reac
rownames(stoch) <- spec


#
# fill stoichiometric matrix
#
tmp<-matrix(numeric(0), dim(stoch)[1],0)
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
  # consider reversible case
  if(Reaction_getReversible(re)){
    tmp <- cbind(tmp, -stoch[,i])
    #colnames(stoch) <- c(colnames(stoch), paste(colnames(stoch[,i]),"rev", sep=""))
  }
}
stoch <- cbind(stoch, tmp)