get_sbml <- function(type){
  if (type=="ecoli") return(ecoli_sbml)
  if (type=="barkeri") return(barkeri_sbml)
}

get_biomassf <- function(type){
  if (type=="ecoli") return(ecoli_biomassf)
  if (type=="barkeri") return(barkeri_biomassf)
}

get_sub_ex <- function(type){
  if (type=="ecoli") return(ecoli_sub_ex)
  if (type=="barkeri") return(barkeri_sub_ex)
}

get_maintenancef <- function(type){
  if (type=="ecoli") return(ecoli_maintenancef)
  if (type=="barkeri") return(barkeri_maintenancef)
}

get_lower_bound <- function(type){
  if (type=="ecoli") return(ecoli_lower_bound)
  if (type=="barkeri") return(barkeri_lower_bound)
}

get_upper_bound <- function(type){
  if (type=="ecoli") return(ecoli_upper_bound)
  if (type=="barkeri") return(barkeri_upper_bound)
}

set_lower_bound <- function(type, substrat){
  if (type=="ecoli") return(ecoli_set_lower_bound(substrat))
  if (type=="barkeri") return(barkeri_set_lower_bound(substrat))
}
