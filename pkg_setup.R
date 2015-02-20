setwd('C:/Users/eugen.bauer/Documents/GitHub/BacArena')
source(file="R/Arena.R")
source(file="R/Substance.R")
source(file="R/Organism.R")

library(roxygen2)
library(devtools)

install_github("euba/BacArena", ref="rpkg")
library(BacArena)

############################################################


roxygenize()
