get_sbml <- function(type){
  switch(type,
         ecoli = return(ecoli_sbml),
         barkeri = return(barkeri_sbml))
}

get_biomassf <- function(type){
  switch(type,
         ecoli = return(ecoli_biomassf),
         barkeri = return(barkeri_biomassf))
}

get_sub_ex <- function(type){
  switch(type,
         ecoli = return(ecoli_sub_ex),
         barkeri = return(barkeri_sub_ex))
}

get_maintenancef <- function(type){
  switch(type,
         ecoli = return(ecoli_maintenancef),
         barkeri = return(barkeri_maintenancef))
}

get_lower_bound <- function(type){
  switch(type,
         ecoli = return(ecoli_lower_bound),
         barkeri = return(barkeri_lower_bound))
}

get_upper_bound <- function(type){
  switch(type,
         ecoli = return(ecoli_upper_bound),
         barkeri = return(barkeri_upper_bound))
}

set_lower_bound <- function(type, substrat){
  switch(type,
         ecoli = return(ecoli_set_lower_bound(substrat)),
         barkeri = return(barkeri_set_lower_bound(substrat)))
}
