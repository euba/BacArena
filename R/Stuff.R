globalVariables(c("Ec_core"))


#color dictionary of 269 maximally distinct colors from all previous colors 
colpal1 <- c("#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B")

# 20 optimally distinct colors
colpal2 <- c("#C48736", "#CE54D1", "#96CED5", "#76D73C", "#403552", "#D4477D", "#5A7E36", "#D19EC4", "#CBC594", "#722A2D", "#D0CD47", "#CF4A31", "#7B6FD0", "#597873", "#6CD3A7", "#484125", "#C17E73", "#688EC1",  "#844081", "#7DD06F")

# K. Kelly (1965): Twenty-two colors of maximum contrast. // Color Eng., 3(6), 1965
colpal3 = c(
  "#FFB300", # Vivid Yellow
  "#803E75", # Strong Purple
  "#FF6800", # Vivid Orange
  "#A6BDD7", # Very Light Blue
  "#C10020", # Vivid Red
  "#CEA262", # Grayish Yellow
  "#817066", # Medium Gray
  "#007D34", # Vivid Green
  "#F6768E", # Strong Purplish Pink
  "#00538A", # Strong Blue
  "#FF7A5C", # Strong Yellowish Pink
  "#53377A", # Strong Violet
  "#FF8E00", # Vivid Orange Yellow
  "#B32851", # Strong Purplish Red
  "#F4C800", # Vivid Greenish Yellow
  "#7F180D", # Strong Reddish Brown
  "#93AA00", # Vivid Yellowish Green
  "#593315", # Deep Yellowish Brown
  "#F13A13", # Vivid Reddish Orange
  "#232C16" # Dark Olive Green
)

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
  #print("test")
  with (as.list(parms), {
    vgrid <- setup.prop.2D(value = 0, y.value=1, grid = gridgeometry.grid2D)
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.x=1, D.y=1, #D.grid = diffgeometry.Dgrid, 
                     v.x = 0, v.y=1, #v.grid = vgrid,
                     #C.y.down=rep(boundS, gridgeometry.grid2D$y.N),
                     #C.y.up=rep(boundS, gridgeometry.grid2D$y.N))$dC
                     
                     #C.y.down=rep(1, gridgeometry.grid2D$x.N),
                     flux.x.down=rep(boundS, gridgeometry.grid2D$y.N),
                     flux.x.up=rep(boundS, gridgeometry.grid2D$y.N))$dC
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
#' evalArena(sim, time=7, phencol = TRUE, 
#'           plot_items=c("Population", "EX_o2(e)", "EX_for(e)",
#'           "EX_glc(e)", "EX_for(e)"))
#'}
openArena <- function(){
  data(Ec_core, envir = environment())
  bac = Bac(model=Ec_core, type="E. coli")
  arena <- Arena(n=50, m=50, stir=F, Lx=0.0125, Ly=0.0125)
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
usd <- function(y){mean(y) + sd(y)}
#' @title Computer standard deviation lower bound
#'
#' @description Helper function to get lower error bounds in plotting
#' @param y Vector with numbers
#' @export
#' @rdname lsd
#'
lsd <- function(y){lb=mean(y)-sd(y); ifelse(lb<0,0,lb)}


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
#'
plotSubCurve <-function(simlist, mediac=NULL, time=c(NULL,NULL), scol=NULL){
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  if(sum(mediac %in% simlist[[1]]@mediac) != length(mediac)) stop("Substance does not exist in exchange reactions.")
  if(all(!is.null(time)) && (!time[1]<time[2] || !time[2]<length(simlist[[1]]@medlist))) stop("Time interval not valid")
  
  if(length(mediac)==0) mediac <- names(getVarSubs(simlist[[1]]))[1:5] # get 5 most varying substances (from first sim)
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
      rownames(mat_nice) <- gsub("\\(e\\)$","", gsub("\\[e\\]$","", gsub("EX_","",rownames(mat_nice))))
    }else{
      mat_nice <- t(as.matrix(mat[which(rownames(mat) %in% mediac),]))
      rownames(mat_nice) <- gsub("\\(e\\)$","", gsub("\\[e\\]$","", gsub("EX_","",mediac)))
    }
    colnames(mat_nice) <- time_seq
    all_df <- rbind(all_df, reshape2::melt(mat_nice))
  }
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + geom_point()
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + stat_smooth() + geom_point()
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + stat_smooth(level=0.999) # confidence level default: 0.95
  
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=Var1), alpha=0.3)
  
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) +
  #  stat_summary(geom="ribbon", fun.ymin="min", fun.ymax="max", aes(fill=Var1), alpha=0.3)
  
  q <- ggplot2::ggplot(all_df, ggplot2::aes(color=Var1, y=value, x=Var2)) +
    ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes(fill=Var1), alpha=0.3) + 
    xlab("time") + ylab("amount of substance [fmol]") + ggtitle("Substance curve with standard deviation") + 
    theme_bw(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour='black',size=2),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5           # Position legend in bottom right
  if(!is.null(scol)) q <- q + scale_fill_manual(values=scol) + scale_color_manual(values=scol)
  print(q)
  
  q <- ggplot2::ggplot(all_df, ggplot2::aes(color=Var1, y=value, x=Var2)) + stat_summary(fun.y = mean, geom="line") + 
    xlab("time") + ylab("amount of substance [fmol]") + ggtitle("Mean substance curve") + 
    theme_bw(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour='black',size=2),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5           # Position legend in bottom right
  if(!is.null(scol)) q <- q + scale_fill_manual(values=scol) + scale_color_manual(values=scol)
  print(q)
}


#' @title Plot growth curve for several simulations
#'
#' @description The function \code{plotGrowthCurve} takes a list of simulations and plots the time course of species with standard deviation.
#' @export
#' @rdname plotGrowthCurve
#' 
#' @param simlist A list of simulations (eval objects).
#' @param time Vector with two entries defining start and end time
#' @param bcol Vector with color that should be used
#'
plotGrowthCurve <-function(simlist, bcol=NULL, time=c(NULL,NULL)){
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  if(all(!is.null(time)) && (!time[1]<time[2] || !time[2]<length(simlist[[1]]@simlist))) stop("Time interval not valid")
  
  all_df <- data.frame()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    if(all(!is.null(time))) time_seq <- seq(time[1],time[2]) else time_seq <- seq_along(object@simlist)
    list <- lapply(time_seq, function(i){
      occ <- table(object@simlist[[i]]$type)
      unlist(lapply(seq_along(object@specs), function(i){ifelse(i %in% names(occ),occ[paste(i)], 0)})) # ugly ;P
    })
    mat_bac  <- do.call(cbind, list)
    rownames(mat_bac) <- names(object@specs)
    colnames(mat_bac) <- time_seq
    all_df <- rbind(all_df, reshape2::melt(mat_bac))
  }
  
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + geom_point()
  q<-ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + #stat_smooth(se = FALSE) + 
    stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", aes(fill=Var1), alpha=0.3) + 
    xlab("Time in h") +
    ylab("Number of individuals") +
    theme_bw(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour='black',size=2),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5           # Position legend in bottom right
  if(!is.null(bcol)) q <- q + scale_fill_manual(values=bcol) + scale_color_manual(values=bcol)
    #xlab("time") + ylab("number organism") + ggtitle("Growth curve with standard deviation")
  print(q)
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
#' @param ret_phengroups True if clustered phenotype groups should be returned. 
#' @param cluster True phenotypes should be clustered/condensed. 
#'
plotPhenCurve <- function(simlist, subs, phens=NULL, time=c(NULL,NULL), ret_phengroups=FALSE, cluster=TRUE){
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
        phenmat_bin <- rbind(phenmat_bin, rep(0,ncol(phenmat_bin)))
        rownames(phenmat_bin)[nrow(phenmat_bin)] <- paste0(names(l@specs)[i],"_",0,"-sim_",j)
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
    colnames(m_group_sub) <- gsub("\\(e\\)","", gsub("EX_","",names(simlist[[1]]@mediac)))[pos]
    print(m_group_sub)
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
      rownames(mat_groups) <- paste0("P", 1:dim(mat_groups)[1])
      colnames(mat_groups) <- time_seq
    } else mat_groups <- mat_phen
    #rownames(mat_groups) <- 1:dim(mat_groups)[1]
    mat_groups_m <- melt(mat_groups)
    colnames(mat_groups_m) <- c("Cphen", "time", "value")
    mat_groups_m$replc <- as.character(i)
    all_df <- rbind(all_df, mat_groups_m)
  }
  #all_df <- all_df[which(all_df$value!=0),]
  
  # 4) plotting
  p <- ggplot(all_df, aes(colour=Var1, y=value, x=Var2)) + #stat_summary(fun.y = mean, geom="line") +
    stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", aes(fill=Var1), alpha=0.3) + 
    scale_fill_manual(values=colpal3) + scale_colour_manual(values=colpal3) +
    xlab("time") + ylab("number organism") + ggtitle("Phenotyp growth curve with standard deviation") + 
    theme_bw(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour='black',size=2, fill=NA),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5   
  print(p)
  
  p <- ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + stat_summary(fun.y = mean, geom="line") +
    scale_colour_manual(values=colpal3) + 
    xlab("time") + ylab("number organism") + ggtitle("Phenotyp growth curve with standard deviation") + 
    theme_bw(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour='black',size=2, fill=NA),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5   
  print(p)
  
  p <- ggplot(all_df, aes(color=replc, group=replc,y=value, x=time)) + geom_line(size=1) +
    scale_colour_manual(values=colpal3) + 
    facet_wrap(~Cphen, nrow=4) +
    #facet_grid(~Cphen) +
    xlab("time") + ylab("number organism") + ggtitle("Phenotype growth curve for each replicate") + 
#    theme_bw(base_size = 30) +
    theme_classic(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=15),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
 #     panel.border = element_rect(colour='black',size=2, fill=NA),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5   
  print(p)

  p <- ggplot(all_df, aes(color=Cphen, group=Cphen,y=value, x=time)) + geom_line(size=1) +
    scale_colour_manual(values=colpal3) + 
    facet_wrap(~replc) +
    xlab("time") + ylab("number organism") + ggtitle("Phenotyp growth curve with standard deviation") + 
    theme_bw(base_size = 30) +
    theme(#legend.position='none',
      legend.text=element_text(size=14),
      legend.key=element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.y = element_text(size=30,vjust=0.5),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour='black',size=2, fill=NA),
      axis.ticks = element_line(size=1,color='black'),
      plot.title = element_text(size=20)) #15x5   
  print(p)
  
    
  
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + geom_point()
  #ggplot(all_df, aes(color=Var1, y=value, x=Var2)) + stat_smooth(level = 0.99) + geom_point()
  if(ret_phengroups & cluster) return(simlist_fac)
}
