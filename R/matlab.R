#' @title Read matlab model
#'
#' @description The generic function \code{readMATmod} imports matlab cobra models into sybil model files
#' @export
#' @rdname readMATmod
#'
#' @param file Full path to matlab model file
#' @details Returns sybil model object (time needed: bacterial model ~ 10s, recon2 ~ 60s)
readMATmod <- function(file){
  library(sybil)
  library(R.matlab)
  library(stringr)
  
  #tmp <- readRDS("~/uni/dat/mod/r/recon2v04.RDS")
  #data <- readMat("~/uni/dat/mod/Recon2.v04.mat")
  
  print(system.time(data <- readMat(file)))
  dat.mat <- data[[1]]
  
  mod.var <- dimnames(dat.mat)[[1]]
  
  # 1) model name
  if( "modelID" %in% mod.var ){
    mod.id <- as.character(dat.mat[[which(mod.var=="modelID")]])
  }else{
    mod.id <- names(data)[1]  
  }
  if( "modelName" %in% mod.var ){
    mod.name <- as.character(dat.mat[[which(mod.var=="modelName")]])
  } else{
    mod.name <- mod.id  
  }
  if( "description" %in% mod.var ){
    mod.desc <- as.character(dat.mat[[which(mod.var=="description")]])
  } else{
    mod.desc <- mod.id  
  }
  
  # 2) stoich matrix
  mod.S <- Matrix(dat.mat[[which(mod.var=="S")]], sparse = T)
  
  # 3) rxn
  mod.react_id   <- unlist(dat.mat[[which(mod.var=="rxns")]])
  mod.react_name <- unlist(dat.mat[[which(mod.var=="rxnNames")]])
  if( "rev" %in% mod.var ){
    mod.react_rev  <- as.vector(dat.mat[[which(mod.var=="rev")]]) == TRUE
  }else{
    mod.react_rev  <- as.vector(dat.mat[[which(mod.var=="lb")]]) < 0
  }
  
  # 4) met
  mod.met_id     <- unlist(dat.mat[[which(mod.var=="mets")]])
  mod.met_name   <- unlist(dat.mat[[which(mod.var=="metNames")]])
  
  # 5) genes
  mod.gene_id      <- unlist(dat.mat[[which(mod.var=="genes")]])
  mod.GeneMat    <- dat.mat[[which(mod.var=="rxnGeneMat")]]
  mod.genes <- apply(mod.GeneMat, 1, function(row){ 
    x <- unname(mod.gene_id[which(row != 0)])
    if( length(x)==0 ) "" else x })
  mod.gpr        <- sapply(sapply(dat.mat[[which(mod.var=="grRules")]], unlist), function(entry){
    if ( length(entry) == 0 ) "" else unname(entry)
  })
  # gpr rules needs to be converted (brackets + numbering)
  if( "rules" %in% mod.var ){
    mod.gprRules <- sapply(sapply(dat.mat[[which(mod.var=="rules")]], unlist), function(entry){
      if( length(entry) == 0 ) "" 
      else {
        numbers <- as.numeric(unlist(str_extract_all(entry, "[0-9]+")))
        dict <- as.character(numbers - min(numbers) + 1); names(dict) <- as.character(numbers)
        gsub("\\(([0-9]+)\\)","\\[\\1\\]",str_replace_all(entry, dict)) }
    })
  }else{ # if 'rules' is not present construct own rules from gpr+genes
    mod.gprRules <- sapply(seq_along(mod.gpr), function(i){
      if( mod.gpr[i] == "") ""
      else{
        genes <- unlist(mod.genes[i])
        dict <- paste0("x[",seq_along(genes),"]")
        names(dict) <- genes
        gsub("and","&", gsub("or","|",str_replace_all(mod.gpr[i], dict)))
      }
    })
    
  }
  
  # 6) bounds
  mod.lb <- as.vector(dat.mat[[which(mod.var=="lb")]])
  mod.ub <- as.vector(dat.mat[[which(mod.var=="ub")]])
  
  # 7) compartments
  met_comp <- str_extract_all(mod.met_id, "(?<=\\[)[a-z](?=\\])")
  mod.mod_compart <- unique(unlist(met_comp))
  mod.met_comp    <- match(str_extract_all(mod.met_id, "(?<=\\[)[a-z](?=\\])"), mod.mod_compart)
  
  # 8) subsystems
  sub <- sapply(dat.mat[[which(mod.var=="subSystems")]], unlist)
  sub.unique <- unique(sub)
  mod.subSys <- Matrix(FALSE, nrow = length(sub), ncol = length(sub.unique), sparse=T)
  for(i in 1:length(sub)){
    j <- match(sub[i], sub.unique)
    mod.subSys[i,j] <- TRUE
  }
  colnames(mod.subSys) <- sub.unique
  
  
  # create new model
  model <- modelorg(id = mod.id, name = mod.name)
  model@mod_desc <- mod.desc
  model@S <- mod.S
  model@lowbnd <- mod.lb
  model@uppbnd <- mod.ub
  model@met_id <- mod.met_id
  model@met_name <- mod.met_name
  model@met_num <- length(mod.met_id)
  model@react_id  <- mod.react_id
  model@react_name  <- mod.react_name
  model@react_num <- length(mod.react_id)
  model@react_rev <- mod.react_rev
  model@genes <- mod.genes
  model@gprRules <- mod.gprRules
  model@gpr <- mod.gpr
  model@mod_compart <- mod.mod_compart
  model@met_comp <- mod.met_comp
  model@subSys <- mod.subSys
  
  obj.idx <- which(dat.mat[[which(mod.var == "c")]]!=0)
  if( length(obj.idx) > 0 ){
    model <- changeObjFunc(model, react = obj.idx)  
    print(optimizeProb(model))
  }
  
  return(model)
}
