#' BacArena: An Agent-Based Modeling Framework for Cellular Communities
#'
#' The BacArena package provides six classes: Arena (subclass Eval), Organism (subclasses Bac, Human) and Substance.
#' Accordingly there are three categories of important functions:
#' Arena, Organism and Substance.
#' 
#' @section Arena functions:
#' The Arena functions ...
#' @section Organism functions:
#' The Organism functions ...
#' @section Substance functions:
#' The Substance functions ...
#'
#' @docType package
#' @name BacArena
NULL
#> NULL


#' @import methods
NULL


.onAttach <- function(libname, pkgname) {
      packageStartupMessage("BacArena paper: https://doi.org/10.1371/journal.pcbi.1005544\n Tutorials: https://bacarena.github.io\n Model import from SBML: https://github.com/euba/BacArena/wiki/Model-import\n Development and help: https://github.com/euba/bacarena")
}
