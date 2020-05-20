globalVariables(c("colpal1", "colpal2", "colpal3", "Ec_core"))

# status code of linear optimization depending on solver
getLPstat <- function(opt_sol, solver){
  switch(solver,
         glpkAPI =     {solve_ok <- opt_sol$stat==5},
         cplexAPI =    {solve_ok <- opt_sol$stat==1},
         clpAPI =      {solve_ok <- opt_sol$stat==0},
         lpSolveAPI =  {solve_ok <- opt_sol$stat==0},
         sybilGUROBI = {solve_ok <- opt_sol$stat==2},
         stop("Solver not suported!"))
  return(solve_ok)
}

# Diffusion pde solver function
Diff2d <- function (t, y, parms){
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- ReacTran::tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid)$dC
    return (list(dCONC))
  })
}

# advection + diffusion + conservation
AdvecDiffConserv2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    vgrid <- setup.prop.2D(value = 1, y.value=0, grid = gridgeometry.grid2D)
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     v.grid = vgrid, flux.x.down=0, flux.x.up=0, flux.y.up=0, flux.y.down=0)$dC
    return (list(dCONC))
  })
}

AdvecDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    vgrid <- setup.prop.2D(value = 1, y.value=0, grid = gridgeometry.grid2D)
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    #dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, v.grid = vgrid, C.x.down=rep(0.1, gridgeometry.grid2D$x.N))$dC
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     v.grid = vgrid, 
                     C.y.down=rep(1, gridgeometry.grid2D$x.N),
                     flux.x.down=1, flux.x.up=-1, flux.y.up=0, flux.y.down=0)$dC
    return (list(dCONC))
  })
}

ConstBoundDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     C.y.down=rep(boundS, gridgeometry.grid2D$y.N),
                     C.y.up=rep(boundS, gridgeometry.grid2D$y.N),
                     C.x.down=rep(boundS, gridgeometry.grid2D$x.N),
                     C.x.up=rep(boundS, gridgeometry.grid2D$x.N))$dC
    return (list(dCONC))
  })
}

InfluxBoundDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     flux.y.down=rep(-boundS, gridgeometry.grid2D$y.N),
                     flux.y.up=rep(boundS, gridgeometry.grid2D$y.N),
                     flux.x.down=c(0,rep(-boundS, gridgeometry.grid2D$x.N-2),0),  # reduce edge bias
                     flux.x.up=c(0, rep(boundS, gridgeometry.grid2D$x.N-2),0))$dC # "
    return (list(dCONC))
  })
}

ConstBoundAdvecDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     v.grid = diffgeometry.Vgrid,
                     
                     #C.y.down=rep(boundS, gridgeometry.grid2D$y.N),
                     #C.y.up=rep(boundS, gridgeometry.grid2D$y.N))$dC
                     
                     
                     #C.y.down=rep(0, gridgeometry.grid2D$y.N),
                     #C.y.up=rep(1, gridgeometry.grid2D$y.N),
                     flux.y.up=rep(1, gridgeometry.grid2D$y.N),
                     flux.x.down=rep(0, gridgeometry.grid2D$x.N),
                     flux.x.up=rep(0, gridgeometry.grid2D$x.N))$dC
                     #flux.x.down=c(0,rep(-boundS, gridgeometry.grid2D$x.N-2),0),  # reduce edge bias
                     #flux.x.up=c(0, rep(boundS, gridgeometry.grid2D$x.N-2),0))$dC # "
    return (list(dCONC))
  })
}







# function estimates lrw array size paremeter needed to solve stiffy diffusion pde with solver lsodes
estimate_lrw <- function(grid_n, grid_m){
  x=c(10*10, 25*25, 51*51, 61*61, 71*71, 81*81, 91*91, 101*101)
  y=c(3901, 29911, 160000, 230000, 330000, 430000, 580000, 710000)
  lm <- lm(y~x)
  #summary(lm)
  #plot(x,y)
  #abline(coef(lm))
  #abline(coef=c(0, lm$coefficients[2]))
  lrw <- as.numeric(lm$coefficients[2]*grid_n*grid_m + grid_n*grid_m*100)
  #lrw <- ((grid_n*grid_m)*18.5 + 20)*10 -> alternative function
  return(lrw)
}

#' @title Start simulation
#'
#' @description The function \code{openArena} can be used to start a default simulation.
#' @export
#' @rdname openArena
#' @importFrom utils data
#'
#' @return Returns an object of class \code{Eval} which can be used for subsequent analysis steps.
#' @examples
#' \donttest{ 
#' sim <- openArena()
#' evalArena(sim, time=5, phencol = TRUE, 
#'           plot_items=c("Population", "EX_o2(e)", "EX_for(e)", "EX_glc(e)"))
#'}
openArena <- function(){
  data(Ec_core, envir = environment())
  bac = Bac(model=Ec_core, type="E. coli")
  arena <- Arena(n=50, m=50)
  arena <- addOrg(arena, bac, amount=50)
  arena <- addSubs(arena, smax=0.05, unit="mM", difspeed=6.7e-6,
            mediac = c("EX_glc(e)", "EX_o2(e)", "EX_h2o(e)", "EX_pi(e)", "EX_nh4(e)", "EX_h(e)"))
  sim <- simEnv(arena, time=5)  
  plotCurves2(sim, legendpos="left")
  return(sim)
}


#' @title Reset plotting screen
#'
#' @description The function \code{reset_screen} set plotting window to default
#' @export
#' @rdname reset_screen
#'
reset_screen <- function(){
  par(mfrow=c(1,1))
}

#' @title Computer standard deviation upper bound
#'
#' @description Helper function to get upper error bounds in plotting
#' @param y Vector with numbers
#' @export
#' @rdname usd
#'
usd <- function(y){mean(y) + stats::sd(y)}
#' @title Computer standard deviation lower bound
#'
#' @description Helper function to get lower error bounds in plotting
#' @param y Vector with numbers
#' @export
#' @rdname lsd
#'
lsd <- function(y){lb=mean(y)-stats::sd(y); ifelse(lb<0,0,lb)}


#' @title Plot substance curve for several simulations
#'
#' @description The function \code{plotSubCurve} takes a list of simulations and plots the time course of substances with standard deviation.
#' @export
#' @rdname plotSubCurve
#' 
#' @param simlist A list of simulations (eval objects).
#' @param mediac A vector of substances (if not specified most varying substances will be taken.)
#' @param time Vector with two entries defining start and end time.
#' @param scol Vector with colors that should be used.
#' @param ret_data Set true if data should be returned
#' @param num_var Number of varying substances to be shown (if mediac is not specified)
#' @param unit Unit for the substances which should be used for plotting (default: mmol)
#' @param useNames Use substance names instead of ids
#' 
#' @return list of three ggplot object for further formating
#'
plotSubCurve <-function(simlist, mediac=NULL, time=c(NULL,NULL), scol=NULL, unit="mmol", ret_data=FALSE, num_var=10, useNames=FALSE){
  if(is(simlist, "Eval")) simlist <- list(simlist)
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  #if(sum(mediac %in% simlist[[1]]@mediac) != length(mediac)) stop("Substance does not exist in exchange reactions.")
  if(all(!is.null(time)) && (!time[1]<time[2] || !time[2]<length(simlist[[1]]@medlist))) stop("Time interval not valid")
  if(length(mediac)==0) mediac <- names(getVarSubs(simlist[[1]]))[1:num_var] # get the most varying substances (from first sim)
  if(length(mediac) == 0) stop("All substance have a variance of zero.")
  mediac <- intersect(mediac,simlist[[1]]@mediac)
  if(length(mediac)==0) stop("Substance does not exist in exchange reactions.")
  
  all_df <- data.frame()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    if(all(!is.null(time))) time_seq <- seq(time[1],time[2]) else time_seq <- seq_along(object@medlist)
    prelist <- lapply(time_seq, function(j){extractMed(object, j)})
    list <- lapply(prelist, function(x){lapply(x, sum)})
    mat <- matrix(unlist(list), nrow=length(object@media), ncol=length(time_seq))
    rownames(mat) <- object@mediac
    if(length(mediac) > 1){
      mat_nice <- mat[which(rownames(mat) %in% mediac),]
      rownames(mat_nice) <- rownames(mat_nice)
    }else{
      mat_nice <- t(as.matrix(mat[which(rownames(mat) %in% mediac),]))
      rownames(mat_nice) <- mediac
    }
    colnames(mat_nice) <- time_seq
    mat_nice <- reshape2::melt(mat_nice)
    mat_nice$replc <- as.character(i)
    all_df <- rbind(all_df, mat_nice)
  }
  
  all_df$Var2 <- all_df$Var2 * simlist[[1]]@tstep # adjust time to hours 
  if( useNames ) all_df$Var1 <- names(simlist[[1]]@mediac)[match(all_df$Var1, simlist[[1]]@mediac)]
  all_df$Var1 <- gsub("\\(e\\)$","", gsub("\\[e\\]$","", gsub("EX_","", all_df$Var1)))
  
  ylabel = "Amount of substance in"
  switch(unit,
         'mmol'={all_df$value <- all_df$value * 10^{-12}; ylabel=paste(ylabel,"mmol")},
         'umol'={all_df$value <- all_df$value * 10^{-9}; ylabel=paste(ylabel,"umol")},
         'nmol'={all_df$value <- all_df$value * 10^{-6}; ylabel=paste(ylabel,"nmol")},
         'pmol'={all_df$value <- all_df$value * 10^{-3}; ylabel=paste(ylabel,"pmol")},
         'fmol'={all_df$value <- all_df$value * 1; ylabel=paste(ylabel,"fmol")},
         'mM'  ={all_df$value <- all_df$value * 10^{-12}*100/(simlist[[1]]@Lx*simlist[[1]]@Ly); ylabel=paste(ylabel,"mM")},
         stop("Wrong unit for concentration."))
  
  colnames(all_df)[1:2] <- c("sub", "time")
  plot_list <- list()
  
  #unique(all_df$sub) = sort(simlist[[1]]@mediac)
  if(length(simlist)>1){
    q1 <- ggplot2::ggplot(all_df, ggplot2::aes_string(color="sub", y="value", x="time")) + ggplot2::geom_line(size=1) + ggplot2::facet_wrap(~replc)    
    
    q3 <- ggplot2::ggplot(all_df, ggplot2::aes_string(color="sub", y="value", x="time")) +
      ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes_string(fill="sub"), alpha=0.3, size=1) +
      ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel)
    
    q4 <- ggplot2::ggplot(all_df, ggplot2::aes_string(color="sub", y="value", x="time")) +
      ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes_string(fill="sub"), alpha=0.3, size=1) +
      ggplot2::facet_wrap(~sub, scales="free_y") + ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel) #+ 
  plot_list <- list(q1, q3, q4)
  }else{
  q4 <- ggplot2::ggplot(all_df, ggplot2::aes_string(color="sub", x="time", y="value")) +
    ggplot2::geom_line(size=1) + 
    ggplot2::facet_wrap(~sub, scales="free_y") + ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel) #+ 
    plot_list <- list(q4)
  }
  q2 <- ggplot2::ggplot(all_df, ggplot2::aes_string(color="sub", y="value", x="time")) + ggplot2::stat_summary(fun.y = mean, geom="line", size=1) + 
        ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel) + ggplot2::ggtitle("Mean substance curve")
  plot_list[[length(plot_list)+1]] <- q2
  
  if(ret_data) return(all_df) else return(plot_list)
}


#' @title Plot growth curve for several simulations
#'
#' @description The function \code{plotGrowthCurve} takes a list of simulations and plots the time course of species with standard deviation.
#' @export
#' @rdname plotGrowthCurve
#' 
#' @param simlist A list of simulations (eval objects).
#' @param time Vector with two entries defining start and end time
#' @param ret_data Set true if data should be returned
#' @param use_biomass If enabled then biomass is used instead of cell number
#' @param specs List of species for which a growth curve should be shown (default: all)
#'
plotGrowthCurve <-function(simlist, time=c(NULL,NULL), ret_data=FALSE, use_biomass=F, specs=NULL){
  if(is(simlist, "Eval")) simlist <- list(simlist)
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  if(all(!is.null(time)) && (!time[1]<time[2] || !time[2]<length(simlist[[1]]@simlist))) stop("Time interval not valid")
  
  all_df <- data.frame()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    if(all(!is.null(time))) time_seq <- seq(time[1],time[2]) else time_seq <- seq_along(object@simlist)
    if(use_biomass){
      list <- lapply(time_seq, function(i){
        sapply(seq_along(object@specs), function(x){sum(object@simlist[[i]]$biomass[which(object@simlist[[i]]$type==x)])})})
    }else{
      list <- lapply(time_seq, function(i){
        occ <- table(object@simlist[[i]]$type)
        unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)}))}
      )}
    mat_bac  <- do.call(cbind, list)
    rownames(mat_bac) <- names(object@specs)
    colnames(mat_bac) <- time_seq
    mat_bac_m <- reshape2::melt(mat_bac)
    colnames(mat_bac_m) <- c("species", "time", "value")
    mat_bac_m$replc <- as.character(i)
    all_df <- rbind(all_df, mat_bac_m)
  }
  
  all_df$time <- all_df$time * simlist[[1]]@tstep # adjust time to hours
  if( length(specs)>0 ) all_df <- all_df[which(all_df$species %in% specs),]
  
  # test if capacity is reached
  cap <- sapply(seq_along(simlist), function(i){
    sim <- simlist[[i]]
    capacity   <- sim@n*sim@m
    org_number <- lapply(sim@simlist, nrow)
    cap_reached<- which(org_number == capacity)
    if(length(cap_reached)>0) min(cap_reached) else NA
  })
  cap <- cap * simlist[[1]]@tstep
  if(length(cap)!=0){dat_cap <- data.frame(replc=seq_along(simlist), cap=cap)}
  dat_cap <- dat_cap[which(!is.na(dat_cap$cap)),]
  
  plot_list <- list()
  if(length(simlist) > 1){ # first two plots only possible if replicates are available
    q1<-ggplot2::ggplot(all_df, ggplot2::aes_string(color="species", y="value", x="time")) + 
      ggplot2::geom_line(size=1) + ggplot2::facet_wrap(~replc) +
      ggplot2::xlab("Time in h") + ggplot2::ylab(if(use_biomass){"Biomass"} else {"Number of individuals"})
    if(length(cap)!=0){q1 <- q1 + ggplot2::geom_vline(data=dat_cap, ggplot2::aes(xintercept=cap))}
    
    q2<-ggplot2::ggplot(all_df, ggplot2::aes_string(color="species", y="value", x="time")) +
      ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes_string(fill="species"), alpha=0.3) + 
      ggplot2::xlab("Time in h") + ggplot2::ylab(if(use_biomass){"Biomass"} else {"Number of individuals"})
    #if(length(cap)!=0){q2 <- q2 + ggplot2::geom_vline(xintercept=min(cap))}
    
    plot_list <- list(q1, q2)
  }
    
  q3<-ggplot2::ggplot(all_df, ggplot2::aes_string(color="species", y="value", x="time")) +
    ggplot2::stat_summary(fun.y = mean, geom="line", size=1) + 
    ggplot2::xlab("Time in h") + ggplot2::ylab(if(use_biomass){"Biomass"} else {"Number of individuals"})
  q4 <- ggplot2::ggplot(all_df, ggplot2::aes_string("time", "value")) +
    ggplot2::stat_summary(fun.y = mean, geom="bar", width=1, position="fill", ggplot2::aes_string(fill="species"))
  if(length(cap)!=0 & all(!is.na(cap))){q3 <- q3 + ggplot2::geom_vline(xintercept=min(cap))}
  if(length(plot_list)==0){
    if(length(simlist[[1]]@specs)>1) plot_list <- list(q3,q4) else plot_list <- q3
  }else {plot_list[[length(plot_list)+1]] <- q3; plot_list[[length(plot_list)+2]] <- q4}     
  
  if(ret_data) return(all_df) else return(plot_list)
}


#' @title Plot growth curve for several simulations
#'
#' @description The function \code{plotPhenCurve} takes a list of simulations and plots the time course of species with standard deviation.
#' @export
#' @rdname plotPhenCurve
#' 
#' @param simlist A list of simulations (eval objects).
#' @param subs A vector of substance names that are used for phenotype clustering.
#' @param phens If phencurve is given then phens specifies the phenotypes which sould be plotted again.
#' @param time Vector with two entries defining start and end time
#' @param cluster True if phenotypes should be clustered/condensed.
#' @param inAll True if only phenotypes which occur in all replicates should be considered
#' @param col Vector with color that should be used
#' @param with_gc True if growth curve of organisms should be included
#' @param return_dat Should data be returned? (default false)
#'
plotPhenCurve <- function(simlist, subs, phens=NULL, time=c(NULL,NULL), cluster=TRUE, inAll=TRUE, col=colpal3, with_gc=FALSE, return_dat=FALSE){
  if(is(simlist, "Eval")) simlist <- list(simlist)
  if(sum(subs %in% simlist[[1]]@mediac) != length(subs)) stop("Substances invalid.")
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  
  # 1) cluster phenotypes according to relevant substrates
  if(cluster){
    pos <- which(simlist[[1]]@mediac %in% subs)
    simlist_prefac <- lapply(seq_along(simlist), function(j){
      l <- simlist[[j]]
      mediac <- gsub("\\(e\\)","", gsub("EX_","",l@mediac))
      phens <- l@phenotypes
      phenmat <- matrix(0, nrow=length(phens), ncol=length(l@mediac))
      colnames(phenmat) <- mediac
      counter = vector("numeric", length(l@specs))
      names(counter) <- lapply(names(l@specs), function(org_name){org_name})
      new_names = unlist(lapply(names(phens), function(x){
        counter[x] <<- counter[x] + 1
        paste0(x, "_", counter[x],"-sim_",j)
      }))
      rownames(phenmat) <- new_names
      for(i in 1:nrow(phenmat)){
        phenmat[i,] = as.numeric(unlist(strsplit(phens[i],split={})))
      }
      phenmat_bin <- replace(phenmat, phenmat==2, -1)
      for(i in seq_along(l@specs)) { # add inactive phenotypes
        phenmat_bin <- rbind(rep(0,ncol(phenmat_bin)), phenmat_bin)
        rownames(phenmat_bin)[1] <- paste0(names(l@specs)[i],"_",0,"-sim_",j)
      }
      small <- phenmat_bin[,pos]
      if(length(subs)>1) prefac <- apply(small, 1, paste, collapse="_") else prefac <- small
    })
    
    prefac <- unlist(simlist_prefac)
    fac <- factor(prefac, levels=unique(prefac), labels=1:length(unique(prefac)))
    fac2 <- fac
    names(fac2) <- gsub("-sim_.$", "", names(fac)) # remove 'simulation i' tag
    simlist_fac <- split(fac2, unlist(lapply(seq_along(simlist), function(i){rep(i, length(simlist_prefac[[i]]))}))) # split unique pheno ids according for each simulation
  

    # 2) print legend
    group_mem <- unique(cbind(unname(fac),names(prefac)))
    group_mem <- split(group_mem[,2], paste0("P",group_mem[,1]))
    print(group_mem)
    
    group_sub <- unique(cbind(unname(fac),unname(prefac)))
    if(length(subs)>1){
      m_group_sub <- matrix(as.numeric(unlist(strsplit(group_sub[,2], "_"))), byrow=TRUE, ncol=length(subs))
    }else m_group_sub <- matrix(group_sub[,2])
    rownames(m_group_sub) <- paste0("P", group_sub[,1])
    colnames(m_group_sub) <- gsub("\\(lu\\)","", gsub("EX_","",names(simlist[[1]]@mediac)))[pos]
    #print(m_group_sub)
    m_group_sub2 <- m_group_sub
    colnames(m_group_sub2) <- gsub("\\[e\\]","", gsub("EX_","", colnames(m_group_sub2)))
    m_group_sub2[m_group_sub2==-1] <- "-"
    m_group_sub2[m_group_sub2==1] <- "+"
    m_group_sub2[m_group_sub2==0] <- ""
    print(noquote(m_group_sub2))
  }
  
  # 3) get growth of selected phenotypes
  all_df <- data.frame()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    pheno_nr <- table(names(object@phenotypes))
    if(all(!is.null(time))) time_seq <- seq(time[1],time[2]) else time_seq <- seq_along(object@simlist)
    list <- lapply(time_seq, function(t){ # time step
      unlist(lapply(seq_along(object@specs), function(j){ # bac type
        occ <- table(object@simlist[[t]][which(object@simlist[[t]]$type==j),]$phenotype)
        p <- unlist(lapply(seq(0,pheno_nr[[names(object@specs[j])]]), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
        names(p) <- paste0(names(object@specs)[j], "_pheno", seq(0,pheno_nr[[names(object@specs[j])]]))
        p
      }))})
    mat_phen  <- do.call(cbind, list)
    if(cluster){
      mat_groups <- rowsum(mat_phen, group=simlist_fac[[i]]) # phenotype0 singularity
      rownames(mat_groups) <- paste0("P", rownames(mat_groups))
      colnames(mat_groups) <- time_seq
    } else mat_groups <- mat_phen
    #rownames(mat_groups) <- 1:dim(mat_groups)[1]
    mat_groups_m <- reshape2::melt(mat_groups)
    colnames(mat_groups_m) <- c("Cphen", "time", "value")
    mat_groups_m$replc <- as.character(i)
    all_df <- rbind(all_df, mat_groups_m)
  }
  
  all_df$time <- all_df$time * simlist[[1]]@tstep # adjust time to hours
  
  if(inAll){ # only consider phenotypes occurring in all replicates
    t_pr <- table(all_df[which(all_df$value!=0),]$Cphen, all_df[which(all_df$value!=0),]$replc)
    #print(t_pr)  
    phen_inall <- rownames(t_pr)[apply(t_pr, 1, function(row){all(row!=0)})]
    all_df <- all_df[which(all_df$Cphen %in% phen_inall),]
  }
  
  if(with_gc){ # add overall growth curves
    gc_dat <- plotGrowthCurve(simlist, ret_data = TRUE)
    phen_inall <- c(phen_inall, as.character(unique(gc_dat$species)))
    names(gc_dat)[names(gc_dat)=="species"] <- "Cphen"
    all_df <- rbind(all_df, gc_dat)
  }
  
  #avg_occ <- round(rowsum(all_df$value, group=all_df$Cphen) / (length(simlist) * simlist[[1]]@tstep * length(simlist[[1]]@simlist)),3)
  #print("average occurence of phenotypes"); print(avg_occ)
  #sel_phen <- rownames(avg_occ)[which(avg_occ >= min_occ)]
  #all_df <- all_df[which(all_df$Cphen %in% sel_phen),]
  
  p2 <- ggplot2::ggplot(all_df[which(all_df$value!=0),], ggplot2::aes_string(y="time", x="Cphen")) + ggplot2::geom_boxplot(ggplot2::aes_string(fill="Cphen")) + 
    ggplot2::ylab("time [h]") + ggplot2::xlab("Phenotypes") + ggplot2::coord_flip() + ggplot2::scale_fill_manual(values=col)
  
  p1 <- ggplot2::ggplot(all_df[which(all_df$Cphen %in% phen_inall),], ggplot2::aes_string(colour="Cphen", y="value", x="time")) + 
    ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes_string(fill="Cphen"), alpha=0.3, size=1) + 
    ggplot2::scale_fill_manual(values=col) + ggplot2::scale_colour_manual(values=col) +
    ggplot2::xlab("Time in h") + ggplot2::ylab("Number of individuals") + ggplot2::ggtitle("Phenotype growth curve with standard deviation")

  if(return_dat) {if(cluster) return(list(simlist_fac, all_df)) else return(all_df)}
  else(return(list(p1,p2)))
}


#' @title Plot number of phenotypes curve for several simulations
#'
#' @description The function \code{plotPhenNum} takes a list of simulations and plots the time course of the number of phenotypes with standard deviation.
#' @export
#' @rdname plotPhenNum
#' 
#' @param simlist A list of simulations (eval objects).
#' @param size A scaling factor for plot text and line size
#' @param title Title of the plot
#'
plotPhenNum <-function(simlist, title="Phenotype number variation", size=1){
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  pdat = data.frame()
  for(i in 1:length(simlist)){
    pmat = matrix(0,length(simlist[[1]]@specs),length(simlist[[i]]@simlist))
    rownames(pmat) = names(simlist[[1]]@specs)
    for(j in 2:length(simlist[[i]]@simlist)){
      datp = simlist[[i]]@simlist[[j]]
      for(k in 1:length(simlist[[1]]@specs)){
        pmat[k,j] = length(unique(datp[which(datp$type==k),"phenotype"]))
      }
    }
    pdat = rbind(pdat, reshape2::melt(pmat))
  }
  colnames(pdat) = c("spec","time","phens")
  q <- ggplot2::ggplot(pdat, ggplot2::aes(color=pdat$spec, y=pdat$phens, x=pdat$time)) +
    ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes(fill=pdat$spec), alpha=0.3) + 
    ggplot2::xlab("Time") + ggplot2::ylab("Number of phenotypes") + 
    ggplot2::ggtitle(title) + 
    ggplot2::theme_bw(base_size = 30*size) +
    ggplot2::theme(legend.position='none',
          legend.text= ggplot2::element_text(size=14*size),
          legend.key=ggplot2::element_blank(),
          legend.title =ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(size=20*size),
          axis.text.y = ggplot2::element_text(size=20*size),
          axis.title.y = ggplot2::element_text(size=30*size,vjust=0.5),
          #panel.grid.major =ggplot2::element_blank(),
          panel.grid.minor =ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour='black',size=2*size),
          axis.ticks = ggplot2::element_line(size=1*size,color='black'),
          plot.title = ggplot2::element_text(size=30*size)) #15x5           # Position legend in bottom right
  return(q)
}


#' @title Plot number of variation in number of interactions for several simulations
#'
#' @description The function \code{plotInterNum} takes a list of simulations and plots the time course of the number of metabolic interactions with standard deviation.
#' @export
#' @rdname plotInterNum
#' 
#' @param simlist A list of simulations (eval objects).
#' @param size A scaling factor for plot text and line size
#' @param title Title of the plot
#'
plotInterNum <-function(simlist, title="Variation in number of interactions", size=1){
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  imat = matrix(0,length(simlist),length(simlist[[1]]@simlist))
  for(i in 1:length(simlist)){
    for(j in 2:length(simlist[[i]]@simlist)){
      imat[i,j] = length(which(apply(getPhenoMat(simlist[[i]],j-1),2,sum)==3))
    }
  }
  idat = reshape2::melt(imat)
  colnames(idat) = c("rep","time","inter")
  q <- ggplot2::ggplot(idat, ggplot2::aes(y=idat$inter, x=idat$time)) +
    ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", alpha=0.3) + 
    ggplot2::stat_summary(fun.y = mean, geom="line") + 
    ggplot2::xlab("Time") + ggplot2::ylab("Number of phenotypes") + 
    ggplot2::ggtitle(title) + 
    ggplot2::theme_bw(base_size = 30*size) +
    ggplot2::theme(legend.position='none',
          legend.text= ggplot2::element_text(size=14*size),
          legend.key=ggplot2::element_blank(),
          legend.title =ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(size=20*size),
          axis.text.y = ggplot2::element_text(size=20*size),
          axis.title.y = ggplot2::element_text(size=30*size,vjust=0.5),
          #panel.grid.major =ggplot2::element_blank(),
          panel.grid.minor =ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour='black',size=2*size),
          axis.ticks = ggplot2::element_line(size=1*size,color='black'),
          plot.title = ggplot2::element_text(size=30*size)) #15x5           # Position legend in bottom right
  return(q)
}


#' @title Plot abundances of species
#'
#' @description The function \code{plotAbundance} takes a list of simulations and return a boxplot with species abundances 
#' @export
#' @rdname plotAbundance
#' 
#' @param simlist A list of simulations (eval objects).
#' @param time A vector with start and end time to be considered (default: total time)
#' @param col Vector with color that should be used
#' @param return_dat Should plain text mean abundances be returned? (default false)
#' @param use_biomass If enabled then biomass is used instead of cell number
#'
plotAbundance <- function(simlist, time=c(NULL,NULL), col=colpal3, return_dat=F, use_biomass=F){
  if(is(simlist, "Eval")) simlist <- list(simlist)
  all_df <- data.frame()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    if(all(!is.null(time))) time_seq <- seq(time[1],time[2]) else time_seq <- seq_along(object@simlist)
    if(use_biomass){
      list <- lapply(time_seq, function(i){
        sapply(seq_along(object@specs), function(x){sum(object@simlist[[i]]$biomass[which(object@simlist[[i]]$type==x)])})})
    }else{
      list <- lapply(time_seq, function(i){
        occ <- table(object@simlist[[i]]$type)
        unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)}))})}
    mat_bac  <- do.call(cbind, list)
    rownames(mat_bac) <- names(object@specs)
    colnames(mat_bac) <- time_seq
    mat_bac_m <- reshape2::melt(mat_bac)
    colnames(mat_bac_m) <- c("species", "time", "value")
    mat_bac_m$replc <- as.character(i)
    all_df <- rbind(all_df, mat_bac_m)
  }
  if(return_dat){
    abundances<-unlist(lapply(unique(all_df$species), function(s){
      mean(all_df$value[which(all_df$species==s)])}))
    names(abundances) <- unique(all_df$species)
    return(abundances)
  }else{
    q <- ggplot2::ggplot(all_df, ggplot2::aes_string("species", "value")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="species", fill="species"), alpha = 0.2, outlier.size=1) + 
      ggplot2::scale_fill_manual(values=col) + ggplot2::scale_color_manual(values=col) + 
      ggplot2::theme(axis.text.x =ggplot2::element_blank(), legend.title=ggplot2::element_blank(),axis.title.x = ggplot2::element_blank(),axis.title.y = ggplot2::element_blank())
    return(q)
  }
}

#' @title Plot substance variations
#'
#' @description The function \code{plotSubVar} takes a list of simulations and returns a barplot with most varying substances
#' @export
#' @rdname plotSubVar
#' 
#' @param simlist A list of simulations (eval objects).
#' @param metsel A vector with the name of exchange reactions of interest
#'
plotSubVar <- function(simlist, metsel){
  varsub = do.call(rbind,lapply(simlist,function(x,metsel){
    getVarSubs(x)[metsel]
  },metsel=metsel))
  #varsub = log10(varsub)
  vardat = data.frame("mean"=apply(varsub,2,mean),
                      "sd"=apply(varsub,2,stats::sd),
                      "mets"=factor(colnames(varsub),levels=metsel))
  p <- ggplot2::ggplot(vardat, ggplot2::aes(fill=vardat$mets, y=vardat$mean, x=vardat$mets)) +
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggplot2::geom_errorbar(ggplot2::aes(ymax=vardat$mean+vardat$sd,ymin=vardat$mean-vardat$sd), position=ggplot2::position_dodge(width=0.9), width=0.25)
  return(p)
}

#' @title Plot population flux variations
#'
#' @description The function \code{plotFluxVar} takes a list of simulations and metabolites, returning a plot with metabolite fluxes for each species
#' @export
#' @rdname plotFluxVar
#' 
#' @param simlist A list of simulations (eval objects).
#' @param metsel A vector with the name of exchange reactions of interest
#'
plotFluxVar <- function(simlist, metsel){
  if(is(simlist, "Eval")) simlist <- list(simlist)
  
  concdat=data.frame()
  concmean=data.frame()
  for(i in 1:length(metsel)){
    meanmat=matrix(0,ncol=length(simlist),nrow=length(simlist[[1]]@simlist)*length(simlist[[1]]@specs))
    for(j in 1:length(simlist)){
      conc = do.call(cbind,lapply(simlist[[j]]@mfluxlist, function(x,mets){
        unlist(lapply(x,function(x){x[mets]}))
      },mets=metsel[i]))
      conc = ifelse(is.na(conc),0,conc)
      rownames(conc) = unlist(lapply(strsplit(rownames(conc),split="\\."),function(x){x[1]}))
      cdat = reshape2::melt(conc)
      meanmat[,j] = cdat$value
      cdat$rep = rep(j,nrow(cdat))
      cdat$met = rep(metsel[i],nrow(cdat))
      concdat = rbind(concdat,cdat) 
    }
    concmean=rbind(concmean,data.frame(org=cdat$Var1,time=cdat$Var2,mean=apply(meanmat,1,mean),
                                       sd=apply(meanmat,1,stats::sd),met=cdat$met))
  }
  
  msum = vector("numeric",length=length(unique(paste(unique(paste(concmean$org,concmean$met))))))
  names(msum) = unique(paste(concmean$org,concmean$met))
  rownames(concmean) = paste(concmean$org,concmean$met,concmean$time)
  for(i in unique(concmean$time)){
    for(j in names(msum)){
      msum[j] = msum[j] + concmean[paste(j,i),"mean"]
    }
  }
  msum[is.na(msum)] = 0
  rmind = which(rownames(concmean) %in% paste(names(which(msum==0)),unique(concmean$time)))
  if(length(rmind)!=0){concmean=concmean[-rmind,]}
  concmean$time = concmean$time-1
  g <- ggplot2::ggplot(concmean, ggplot2::aes(color=concmean$met, y=concmean$mean, x=concmean$time, shape=concmean$org)) +
    ggplot2::geom_line(ggplot2::aes(x=concmean$time,y=concmean$mean,color=concmean$met),size=0.6) +
    ggplot2::scale_shape_manual(values=1:length(unique(concmean$org))) +
    ggplot2::geom_point(ggplot2::aes(shape=concmean$org),stroke=0.8,size=2.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax=concmean$mean+concmean$sd, ymin=concmean$mean-concmean$sd), width=0.2)
  return(g)
}

#' @title Function to plot usage of substances species wise 
#'
#' @description The generic function \code{plotSubUsage} displays for given substances the quantities of absorption and production for each species
#' @export
#' @rdname plotSubUsage
#'
#' @param simlist An object of class Eval or a list with objects of class Eval.
#' @param subs List of substance names
#' @param cutoff Total values below cutoff will be dismissed
#' @param ret_data Set true if data should be returned
#' @details Returns ggplot objects
plotSubUsage <- function(simlist, subs=vector(), cutoff=1e-2, ret_data=FALSE){
  
  if(is(simlist, "Eval")) simlist <- list(simlist)
  subs = intersect(subs, simlist[[1]]@mediac)
  if(length(subs)==0){ subs <- names(getVarSubs(simlist[[1]], size = 9))
  }else if(sum(subs %in% simlist[[1]]@mediac) != length(subs)) stop("Substance do not exist in arena")
  
  df <- data.frame(spec=as.character(), sub=as.character(), mflux=as.numeric(), time=as.integer())
  
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    for(t in seq_along(object@mfluxlist)){
      for(spec in names(object@specs)){
        available_subs <- intersect(subs,unname(object@specs[[spec]]@medium))
        if(length(available_subs) > 0 && length(names(object@mfluxlist[[t]][[spec]])) > 0 ){
          mflux=object@mfluxlist[[t]][[spec]][which(names(object@mfluxlist[[t]][[spec]]) %in% available_subs )]
          df <- rbind(df, data.frame(spec=spec, sub=names(mflux), mflux=unname(mflux), time=t, replc=i))
        }
      }
    }
  }
  
  if(nrow(df)==0) stop("None of the substances you provided are used by the community.")
  if(!ret_data) df <- df[which(abs(df$mflux) > cutoff),,drop = FALSE] # do not drop if data is used further
  
  
  if(length(unique(df$sub))>1){
    q1 <- ggplot2::ggplot(df, ggplot2::aes_string(x="time", y="mflux")) + ggplot2::geom_line(ggplot2::aes_string(col="spec"), size=1) + 
      ggplot2::facet_wrap(~sub, scales="free_y")+ ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
    
    q2 <- ggplot2::ggplot(df, ggplot2::aes_string("spec", "mflux")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="spec", fill="spec"), alpha=0.2) + 
      ggplot2::facet_wrap(~sub, scales="free_y") + ggplot2::theme(legend.title=ggplot2::element_blank(), axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  }else{ # if less than 2 substances are found than do not use facet_wrap
    q1 <- ggplot2::ggplot(df, ggplot2::aes_string(x="time", y="mflux")) + ggplot2::geom_line(ggplot2::aes_string(col="spec"), size=1) + 
      ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
    
    q2 <- ggplot2::ggplot(df, ggplot2::aes_string("spec", "mflux")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="spec", fill="spec"), alpha=0.2) + 
      ggplot2::theme(legend.title=ggplot2::element_blank(), axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  } 
  if(ret_data) return(df) else return(list(q1, q2))
}


#' @title Function to plot substance usage for every species
#'
#' @description The generic function \code{plotSpecActivity} displays the input/output substances with the highest variance (could also be defiened manually) for all species
#' @export
#' @rdname plotSpecActivity
#'
#' @param simlist An object of class Eval or a list with objects of class Eval.
#' @param subs List of substance names
#' @param var_nr Number of most varying substances to be used (if subs is not specified)
#' @param spec_list List of species names to be considered (default all)
#' @param ret_data Set true if data should be returned
#' @param useNames Use substance names instead of ids
#' @param rm_unused Remove substances which do not change from plot
#' @param cutoff Minimum valu for fluxes to be considered
#' @details Returns ggplot objects
plotSpecActivity <- function(simlist, subs=list(), var_nr=10, spec_list=NULL, ret_data=FALSE, useNames=FALSE, rm_unused=TRUE, cutoff=1e-6){
  
  if(is(simlist, "Eval")) simlist <- list(simlist)
  if(length(subs)==0) {subs_tocheck <- names(getVarSubs(simlist[[1]]))
  }else subs_tocheck <- subs
  if(length(spec_list)>1) #If spec list is provided, then the system takes the index numbers and replaces them with their real names like in line 742
  {for (i in (spec_list)) {spec_list[which(spec_list==i)]<-names(simlist[[1]]@specs)[i]}}
  if(length(spec_list)==0) spec_list <- names(simlist[[1]]@specs)
  df <- data.frame(spec=as.character(), sub=as.character(), mflux=as.numeric(), time=as.integer())
  #browser()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]  
    for(t in seq_along(object@mfluxlist)){
      for(spec in spec_list){
        if(length(intersect(subs_tocheck, unname(object@specs[[spec]]@medium))) > 0 &  length(names(object@mfluxlist[[t]][[spec]])) > 0 ){
          mflux=object@mfluxlist[[t]][[spec]][which(names(object@mfluxlist[[t]][[spec]]) %in% subs_tocheck)]
          df <- rbind(df, data.frame(spec=spec, sub=names(mflux), mflux=unname(mflux), time=t, replc=i))
        }
      }
    }
  }
  
  if(length(subs)==0){ # in case subs is not specified take substances with highest variance
    mflux_var <- unlist(lapply(unique(df$sub), function(sub){
      stats::var(df[which(df$sub==sub),]$mflux)
    }))
    names(mflux_var) <- unique(df$sub)
    mflux_var <- sort(mflux_var, decreasing = TRUE)
    df <- df[which(df$sub %in% names(mflux_var)[1:var_nr]),]
  }
  df$time = df$time-1
  
  if( rm_unused ){
    unused <- subs[sapply(subs_tocheck, function(sub){all(abs(df[df$sub==sub,"mflux"])<=cutoff)})]
    if( length(unused)>0 ) df <- df[!df$sub %in% unused,]
  }

  if( useNames ) df$sub <- names(simlist[[1]]@mediac)[match(df$sub, simlist[[1]]@mediac)]
    
  q1 <- ggplot2::ggplot(df, ggplot2::aes_string(x="time", y="mflux")) + ggplot2::geom_line(ggplot2::aes_string(col="sub"), size=1) + 
        ggplot2::facet_wrap(~spec, scales="free_y") + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  # q1_5 is the same as q2 but contains "+ ggplot2::facet_wrap(~spec, scales="free_y")" to plot the variance for each spec.   
  q1_5 <- ggplot2::ggplot(df, ggplot2::aes_string("sub", "mflux")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="sub", fill="sub"), alpha=0.2) + 
    ggplot2::theme(axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)") + ggplot2::facet_wrap(~spec, scales="free_y")
    
  q2 <- ggplot2::ggplot(df, ggplot2::aes_string("sub", "mflux")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="sub", fill="sub"), alpha=0.2) + 
    ggplot2::theme(axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  #if(length(unique(df$spec)) > 1) q2 <- q2 + ggplot2::facet_wrap(~spec, scales="free_y") Not needed anymore because it is included in q1_5
  
  if(ret_data) return(df) else return(list(q1, q1_5, q2))
}


#' @title Function to plot population level minimum and maximum flux from alternative optimal solutions obtained by FVA
#'
#' @description The generic function \code{plotFVA} plots population level flux results obtained from function fluxVarSim
#' @export
#' @rdname plotFVA
#'
#' @param fvares results of FVA results to plot, obtained from function fluxVarSim
#' @param mediac List with substances.
#' @details Returns ggplot objects
plotFVA = function(fvares, mediac){
  for(i in 1:length(fvares)){
    fvares[[i]] = fvares[[i]][mediac,]
  }
  basedat = data.frame()
  for(i in 1:length(fvares)){
    basedat = rbind(basedat,
                    data.frame(rownames(fvares[[i]]),
                               fvares[[i]],
                               rep(i,nrow(fvares[[i]]))))
  }
  colnames(basedat) = c("met","min","max","time")
  fvag = ggplot2::ggplot(basedat, ggplot2::aes_string(x="time", y="min")) +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin="min",ymax="max",color="met",fill="met"),alpha=0.2) +
    ggplot2::xlab("Time in h") + ggplot2::ylab("Population flux")
  return(fvag)
}


#' @title Function to plot reaction activity for every species
#'
#' @description The generic function \code{plotReaActivity} displays the usage of reactions for all species
#' @export
#' @rdname plotReaActivity
#'
#' @param simlist An object of class Eval or a list with objects of class Eval.
#' @param reactions List of reaction names
#' @param spec_list List of species names to be considered (default all)
#' @param ret_data Set true if data should be returned
#' @details Returns ggplot objects
plotReaActivity <- function(simlist, reactions=list(), spec_list=NULL, ret_data=FALSE){
  
  if(is(simlist, "Eval")) simlist <- list(simlist)
  
  if(length(reactions)==0) stop("You have to define reactions")
  if(length(spec_list)==0) spec_list <- names(simlist[[1]]@specs)
  
  df <- data.frame(spec=as.character(), rea=as.character(), mflux=as.numeric(), time=as.integer())
  
  for(i in seq_along(simlist)){
    object <- simlist[[i]]  
    for(t in seq_along(object@mfluxlist)){
      for(spec in spec_list){
        if(length(intersect(reactions, names(object@mfluxlist[[i]][[spec]]))) > 0 &  length(names(object@mfluxlist[[t]][[spec]])) > 0 ){
          mflux=object@mfluxlist[[t]][[spec]][which(names(object@mfluxlist[[t]][[spec]]) %in% reactions)]
          df <- rbind(df, data.frame(spec=spec, rea=names(mflux), mflux=unname(mflux), time=t, replc=i))
        }
      }
    }
  }
  if(length(unique(df$rea))>=1){
    q1 <- ggplot2::ggplot(df, ggplot2::aes_string(x="time", y="mflux")) + ggplot2::geom_line(ggplot2::aes_string(col="rea"), size=1) + 
      ggplot2::facet_wrap(~spec, scales="free_y") + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
    
    q2 <- ggplot2::ggplot(df, ggplot2::aes_string("rea", "mflux")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="rea", fill="rea"), alpha=0.2) + 
      ggplot2::theme(axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
    if(length(unique(df$spec)) > 2) q2 <- q2 + ggplot2::facet_wrap(~spec, scales="free_y")
    
    df2 <- plyr::ddply(df, c("time","rea"), function(tmp) c("mflux"=sum(tmp$mflux) ))
    q3 <- ggplot2::ggplot(df2, ggplot2::aes_string("time", "rea")) + ggplot2::geom_tile(ggplot2::aes_string(fill = "mflux")) 
  }else{ # if less than 2 substances are found than do not use facet_wrap
    q1 <- ggplot2::ggplot(df, ggplot2::aes_string(x="time", y="mflux")) + ggplot2::geom_line(ggplot2::aes_string(col="spec"), size=1) + 
      ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
    
    q2 <- ggplot2::ggplot(df, ggplot2::aes_string("spec", "mflux")) + ggplot2::geom_boxplot(ggplot2::aes_string(color="spec", fill="spec"), alpha=0.2) + 
      ggplot2::theme(legend.title=ggplot2::element_blank(), axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
    
    df2 <- plyr::ddply(df, c("time","rea"), function(tmp) c("mflux"=sum(tmp$mflux) ))
    q3 <- ggplot2::ggplot(df2, ggplot2::aes_string("time", "rea")) + ggplot2::geom_tile(ggplot2::aes_string(fill = "mflux"))
  } 

  if(ret_data) return(df) else return(list(q1, q2, q3))
}

#' @title Function for investigation of cross feeding patterns of replicate simulations
#'
#' @description The generic function \code{findFeeding3rep} investigates the cross feeding patterns of replicate simulations
#' @export
#' @rdname findFeeding3rep
#' @importFrom igraph V E graph.data.frame layout.circle
#' 
#' @param simlist A list with objects of class Eval.
#' @param time A numeric vector giving the simulation steps which should be plotted. 
#' @param mets Character vector of substance names which should be considered
#' @param plot Should the graph also be plotted?
#' @param mfunction Function by which the replicate simulations should be combined e.g. "mean" or "median"
#' @return Graph (igraph)
findFeeding3rep <- function(simlist, time, mets, plot=TRUE, mfunction="mean"){
  time = time+1
  mfluxmatrep = list()
  for(i in 1:length(simlist)){
    mets = intersect(simlist[[i]]@mediac,as.character(mets))
    mflux = simlist[[i]]@mfluxlist[[time]]
    mfluxmat = do.call(cbind,lapply(mflux,function(x){return(ifelse(is.na(x[mets]),0,x[mets]))}))
    rownames(mfluxmat) = mets
    mfluxmatrep[[i]] = mfluxmat
  }
  mfluxmat = matrix(apply(do.call(cbind,lapply(mfluxmatrep,as.vector)),1,mfunction),
                    nrow=nrow(mfluxmat),ncol=ncol(mfluxmat),
                    dimnames=list(rownames(mfluxmat),colnames(mfluxmat)))
  inter = data.frame()
  for(i in rownames(mfluxmat)){
    x = mfluxmat[i,]
    interact = matrix(0,ncol=2,nrow=1)
    for(j in names(which(x<0))){
      if(length(which(x>0))!=0){interact = rbind(interact,cbind(names(which(x>0)),j))}
    }
    interact = interact[-1,]
    if("character" %in% class(interact)){interact = t(as.matrix(interact))}
    if(nrow(interact)!=0){inter = rbind(inter,data.frame(prod=interact[,1],cons=interact[,2],met=i))}
  }
  if(any(dim(inter)==0)) {
    warning("No crossfeeding found. Try other metaboites or time points.")
    g <- igraph::make_empty_graph()
    return(list(inter,g))
  }
  if (plot) {
  g <- igraph::graph.data.frame(inter[,1:2], directed=TRUE)
  l <- igraph::layout.kamada.kawai(g)
  plot(g,edge.color=grDevices::rainbow(length(unique(inter$met)))[as.numeric(inter$met)],
       edge.width=3,edge.arrow.size=0.8,vertex.color=1:length(igraph::V(g)),layout=l)
  legend("bottomright",legend=unique(inter$met),col=grDevices::rainbow(length(unique(inter$met))), pch=19, cex=0.7)
  return(invisible(list(inter,g)))}
  else return(inter)
}
