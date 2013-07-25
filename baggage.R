get_sbml <- function(type){
  if (type=="ecoli") return(ecoli_sbml)
}

get_biomassf <- function(type){
  if (type=="ecoli") return(ecoli_biomassf)
}

get_sub_ex <- function(type){
  if (type=="ecoli") return(ecoli_sub_ex)
}

get_maintenancef <- function(type){
  if (type=="ecoli") return(ecoli_maintenancef)
}

get_lower_bound <- function(type){
  if (type=="ecoli") return(ecoli_lower_bound)
}

get_upper_bound <- function(type){
  if (type=="ecoli") return(ecoli_upper_bound)
}

set_lower_bound <- function(type, substrat, lower_bound){
  if (type=="ecoli") return(ecoli_set_lower_bound(substrat, lower_bound))
}