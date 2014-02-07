# this file is a playground to test how the current classes work in action
# the main goal of this file is to construct a basic framework for BacArena, which can then be merged with diffbac
# it is actually a little bit like diffbac.R, but for the current oop version of BacArena

# load libraries and other R files to have everything in place
library(Rcpp)
library(inline)
library(sybil)
SYBIL_SETTINGS("SOLVER", "clpAPI")
setwd("~/BacArena")
#load class definitions
source(file="class/class_baggage.R")
source(file="class/Arena.R")
source(file="class/Substance.R")
source(file="class/Organism.R")
source(file="class/Bac.R")

#load ecoli core model to play around
load("data/ecore_model.R")
mod <- model

#testing constructor
bac1 = Bac(x=1, y=1, model=mod, growth=1)
bac1 = Bac(x=1, y=1, model=mod, growth=1, objlp="ATPM")
#bac1@model = changeObj(bac1@model, "ATPM")
#org1 <- Organism(x=1, y=2, model=mod)
#org1@model = org1@changeObj(org1@model, "ATPM")

#bac1 = Bac(org1, growth=1) #this does not work, but would be cool if...
#testing constructor
org1 <- Organism(x=1, y=2, model=mod)
#testing constructor
sub1 <- Substance(diffconst=5, n=5, m=5, smax=5)
#testing constructor
arena1 <- Arena(n=1,m=1,iter=1,seed=1,epsilon=1)