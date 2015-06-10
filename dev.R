library(devtools)
library(knitr)

build() #building R package archive

setwd("~/BacArena")
load_all()
check()

build_vignettes()

build_win()


install_github("euba/BacArena", ref="rpkg")
install_local("~/BacArena")

# build vignette
library(rmarkdown)
devtools::use_vignette("BacArena-Introduction")
devtools::use_vignette("my-vignette")
#compile
#Cmd + Shift + K

#
# Rcpp stuff
#
library(Rcpp)
Rcpp.package.skeleton("BacArena") # new dummy package 
compileAttributes() # recreates package skeleton



roxygen2
document()

library(sybil)
library(roxygen2)
ls("package:devtools")

makeNamespace()

library(BacArena)
library(rmarkdown)
