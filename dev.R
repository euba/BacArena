library(devtools)


setwd("~/BacArena")
load_all()
check()

build_vignettes()




# build vignette
library(rmarkdown)
devtools::use_vignette("BacArena-Introduction")
devtools::use_vignette("my-vignette")
#compile
#Cmd + Shift + K