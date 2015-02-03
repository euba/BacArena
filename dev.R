library(devtools)


setwd("~/BacArena")
load_all()
check()

build_vignettes()


install_github("euba/BacArena", ref="rpkg")


# build vignette
library(rmarkdown)
devtools::use_vignette("BacArena-Introduction")
devtools::use_vignette("my-vignette")
#compile
#Cmd + Shift + K

roxygen2

document()

library(sybil)
library(roxygen2)
ls("package:devtools")

makeNamespace()
