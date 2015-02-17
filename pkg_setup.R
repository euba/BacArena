setwd('C:/Users/eugen.bauer/Documents/GitHub/BacArena')
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")

install_github("euba/BacArena", ref="rpkg")
library(BacArena)

############################################################

library(roxygen2)
library(devtools)


roxygenize()
