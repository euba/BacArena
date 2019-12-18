#' @title Read matlab model
#'
#' @description The generic function \code{readMATmod} imports matlab cobra models into sybil model files
#' @export
#' @rdname readMATmod
#'
#' @param file Full path to matlab model file
#' @details Returns sybil model object (time needed: bacterial model ~ 10s, recon2 ~ 60s)
readMATmod <- function(file){

  print(system.time(data <- R.matlab::readMat(file)))
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
  mod.gene_id      <- unname(unlist(dat.mat[[which(mod.var=="genes")]]))
  if( "rxnGeneMat" %in% mod.var ){
    mod.GeneMat    <- dat.mat[[which(mod.var=="rxnGeneMat")]]
    mod.genes <- apply(mod.GeneMat, 1, function(row){ 
      x <- unname(mod.gene_id[which(row != 0)])
      if( length(x)==0 ) "" else x })
    mod.gpr        <- sapply(sapply(dat.mat[[which(mod.var=="grRules")]], unlist), function(entry){
      if ( length(entry) == 0 ) "" else unname(entry)
    })
  }else{ # if gene matrix not present create it
    mod.GeneMat <- Matrix(0,nrow = length(mod.react_id), ncol = length(mod.gene_id))
    mod.genes <- list()
    rules.list <- unlist(dat.mat[[which(mod.var=="rules")]], recursive = F)
    if( length(rules.list) != length(mod.react_id) ) stop("Length of rules not same as length of reactions.")
    for(i in seq_along(mod.react_id)  ){
      rule.tmp <- rules.list[[i]]
      if( length(rule.tmp) == 0 ){
        mod.genes[[i]] <- ""
        next
      } 
      j <- as.numeric(unlist(stringr::str_extract_all(rule.tmp, "(?<=x\\()[0-9]+?(?=\\))")))
      mod.GeneMat[i,j] <- 1
      mod.genes[[i]] <- mod.gene_id[j]
    }
  }
  # gpr rules needs to be converted (brackets + numbering)
  if( "rules" %in% mod.var ){
    mod.gprRules <- sapply(sapply(dat.mat[[which(mod.var=="rules")]], unlist), function(entry){
      if( length(entry) == 0 ) "" 
      else {
        numbers <- as.numeric(unlist(stringr::str_extract_all(entry, "[0-9]+")))
        dict <- as.character(numbers - min(numbers) + 1); names(dict) <- as.character(numbers)
        gsub("\\(([0-9]+)\\)","\\[\\1\\]",stringr::str_replace_all(entry, dict)) }
    })
  }else{ # if 'rules' is not present construct own rules from gpr+genes
    mod.gprRules <- sapply(seq_along(mod.gpr), function(i){
      if( mod.gpr[i] == "") ""
      else{
        genes <- unlist(mod.genes[i])
        dict <- paste0("x[",seq_along(genes),"]")
        names(dict) <- genes
        gsub("and","&", gsub("or","|",stringr::str_replace_all(mod.gpr[i], dict)))
      }
    })
    
  }
  
  # 6) bounds
  mod.lb <- as.vector(dat.mat[[which(mod.var=="lb")]])
  mod.ub <- as.vector(dat.mat[[which(mod.var=="ub")]])
  
  # 7) compartments
  met_comp <- stringr::str_extract_all(mod.met_id, "(?<=\\[)[a-z](?=\\])")
  if( all(sapply(met_comp, length) == 0 ) ){
    met_comp <- stringr::str_extract_all(mod.met_id, "(?<=_)[a-z][0-9]?(?=$)")
  }
  mod.mod_compart <- unique(unlist(met_comp))
  mod.met_comp    <- match(met_comp, mod.mod_compart)
  
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
  model <- sybil::modelorg(id = mod.id, name = mod.name)
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
  model@allGenes <- mod.gene_id
  model@mod_compart <- mod.mod_compart
  model@met_comp <- mod.met_comp
  model@subSys <- mod.subSys
  
  obj.idx <- which(dat.mat[[which(mod.var == "c")]]!=0)
  if( length(obj.idx) > 0 ){
    model <- sybil::changeObjFunc(model, react = obj.idx)  
    print(sybil::optimizeProb(model))
  }
  
  return(model)
}
