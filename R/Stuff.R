globalVariables(c("Ec_core", "colpal1", "colpal2", "colpal3"))


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

#26 distinct colors (Zeileis(2009): Escaping RGBland: Selecting Colors for Statistical Graphics)
colpal4 <- c("#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3", "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4")

# 64 destinct colors
colpal5 <- c("#000000","#00FF00","#0000FF","#FF0000","#01FFFE","#FFA6FE","#FFDB66","#006401","#010067","#95003A","#007DB5","#FF00F6","#FFEEE8","#774D00","#90FB92","#0076FF","#D5FF00","#FF937E","#6A826C","#FF029D","#FE8900","#7A4782","#7E2DD2","#85A900","#FF0056","#A42400","#00AE7E","#683D3B","#BDC6FF","#263400","#BDD393","#00B917","#9E008E","#001544","#C28C9F","#FF74A3","#01D0FF","#004754","#E56FFE","#788231","#0E4CA1","#91D0CB","#BE9970","#968AE8","#BB8800","#43002C","#DEFF74    ","#00FFC6","#FFE502","#620E00","#008F9C","#98FF52","#7544B1","#B500FF","#00FF78","#FF6E41","#005F39","#6B6882","#5FAD4E","#A75740","#A5FFD2","#FFB167","#009BFF","#E85EBE")

# 64 destinc colors
colpal6 <- c("FF0000", "00FF00", "0000FF", "FFFF00", "FF00FF", "00FFFF", "000000", "800000", "008000", "000080", "808000", "800080", "008080", "808080", "C00000", "00C000", "0000C0", "C0C000", "C000C0", "00C0C0", "C0C0C0", "400000", "004000", "000040", "404000", "400040", "004040", "404040", "200000", "002000", "000020", "202000", "200020", "002020", "202020", "600000", "006000", "000060", "606000", "600060", "006060", "606060", "A00000", "00A000", "0000A0", "A0A000", "A000A0", "00A0A0", "A0A0A0", "E00000", "00E000", "0000E0", "E0E000", "E000E0", "00E0E0", "E0E0E0")

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
#' 
#' @return list of three ggplot object for further formating
#'
plotSubCurve <-function(simlist, mediac=NULL, time=c(NULL,NULL), scol=NULL, unit="mmol", ret_data=FALSE, num_var=10){
  if(length(simlist) < 1 | !all(lapply(simlist, class) == "Eval") == TRUE) stop("Simlist is invalid.")
  if(sum(mediac %in% simlist[[1]]@mediac) != length(mediac)) stop("Substance does not exist in exchange reactions.")
  if(all(!is.null(time)) && (!time[1]<time[2] || !time[2]<length(simlist[[1]]@medlist))) stop("Time interval not valid")
  
  if(length(mediac)==0) mediac <- names(getVarSubs(simlist[[1]]))[1:num_var] # get the most varying substances (from first sim)
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
    mat_nice <- reshape2::melt(mat_nice)
    mat_nice$replc <- as.character(i)
    all_df <- rbind(all_df, mat_nice)
  }
  
  
  all_df$Var2 <- all_df$Var2 * simlist[[1]]@tstep # adjust time to hours 
  
  ylabel = "Amount of substance in"
  switch(unit,
         'mmol'={all_df$value <- all_df$value * 10^{-12}; ylabel=paste(ylabel,"mmol")},
         'umol'={all_df$value <- all_df$value * 10^{-9}; ylabel=paste(ylabel,"umol")},
         'nmol'={all_df$value <- all_df$value * 10^{-6}; ylabel=paste(ylabel,"nmol")},
         'pmol'={all_df$value <- all_df$value * 10^{-3}; ylabel=paste(ylabel,"pmol")},
         'fmol'={all_df$value <- all_df$value * 1; ylabel=paste(ylabel,"fmol")},
         'mM'  ={all_df$value <- all_df$value * 10^{-12}/(simlist[[1]]@Lx*simlist[[1]]@Ly); ylabel=paste(ylabel,"mM")},
         stop("Wrong unit for concentration."))
  
  q1 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$Var1, y=all_df$value, x=all_df$Var2)) + ggplot2::geom_line(size=1) + ggplot2::facet_wrap(~all_df$replc)
   
  q2 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$Var1, y=all_df$value, x=all_df$Var2)) + ggplot2::stat_summary(fun.y = mean, geom="line", size=1) + 
        ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel) + ggplot2::ggtitle("Mean substance curve")
   
  q3 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$Var1, y=all_df$value, x=all_df$Var2)) +
        ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes(fill=all_df$Var1), alpha=0.3, size=1) +
        ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel)
   
  q4 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$Var1, y=all_df$value, x=all_df$Var2)) +
        ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes(fill=all_df$Var1), alpha=0.3, size=1) +
        ggplot2::facet_wrap(~all_df$Var1, scales="free_y") + ggplot2::xlab("Time in h") + ggplot2::ylab(ylabel) #+ 
  
  if(ret_data) return(all_df) else return(list(q1, q2, q3, q4))
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
plotGrowthCurve <-function(simlist, bcol=colpal3, time=c(NULL,NULL)){
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
    mat_bac_m <- reshape2::melt(mat_bac)
    colnames(mat_bac_m) <- c("species", "time", "value")
    mat_bac_m$replc <- as.character(i)
    all_df <- rbind(all_df, mat_bac_m)
  }
  
  all_df$time <- all_df$time * simlist[[1]]@tstep # adjust time to hours

  # test if capacity is reached
  cap <- unlist(lapply(simlist, function(sim){
    capacity   <- sim@n*sim@m
    org_number <- lapply(sim@simlist, nrow)
    cap_reached<- which(org_number == capacity)
    if(length(cap_reached)>0) min(cap_reached)
    else NULL
  }))
  cap <- cap * simlist[[1]]@tstep - 1
  if(length(cap)!=0){dat_cap <- data.frame(replc=seq_along(simlist), cap=cap)}
  
  q1<-ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$species, y=all_df$value, x=all_df$time)) + 
    ggplot2::geom_line(size=1) + ggplot2::facet_wrap(~replc) + 
    ggplot2::xlab("Time in h") + ggplot2::ylab("Number of individuals") +
    ggplot2::scale_color_manual(values=bcol)
  if(length(cap)!=0){q1 <- q1 + ggplot2::geom_vline(data=dat_cap, ggplot2::aes(xintercept=cap))}
  
  q2<-ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$species, y=all_df$value, x=all_df$time)) +
    ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes(fill=all_df$species), alpha=0.3) + 
    ggplot2::xlab("Time in h") + ggplot2::ylab("Number of individuals") + ggplot2::scale_color_manual(values=bcol) + ggplot2::scale_fill_manual(values=bcol)
  if(length(cap)!=0){q2 <- q2 + ggplot2::geom_vline(xintercept=min(cap))}
    
  q3<-ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$species, y=all_df$value, x=all_df$time)) +
    ggplot2::stat_summary(fun.y = mean, geom="line", size=1) + 
    ggplot2::xlab("Time in h") + ggplot2::ylab("Number of individuals") + ggplot2::scale_color_manual(values=bcol)
  if(length(cap)!=0){q3 <- q3 + ggplot2::geom_vline(xintercept=min(cap))}
  
  return(list(q1, q2, q3))
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
#' @param col Vector with color that should be used
#'
plotPhenCurve <- function(simlist, subs, phens=NULL, time=c(NULL,NULL), ret_phengroups=FALSE, cluster=TRUE, col=colpal3){
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
  
  all_df <- all_df[which(all_df$value!=0),] # important!
  all_df$time <- all_df$time * simlist[[1]]@tstep # adjust time to hours
  
  # 4) plotting
  p1 <- ggplot2::ggplot(all_df, ggplot2::aes(colour=all_df$Cphen, y=all_df$value, x=all_df$time)) + 
    ggplot2::stat_summary(geom="ribbon", fun.ymin="lsd", fun.ymax="usd", ggplot2::aes(fill=all_df$Cphen), alpha=0.3, size=1) + 
    ggplot2::scale_fill_manual(values=col) + ggplot2::scale_colour_manual(values=col) +
    ggplot2::xlab("time [h]") + ggplot2::ylab("number organism") + ggplot2::ggtitle("Phenotyp growth curve with standard deviation") + 
    ggplot2::theme_bw(base_size = 30) +
    ggplot2::theme(#legend.position='none',
      legend.text= ggplot2::element_text(size=14),
      legend.key=ggplot2::element_blank(),
      legend.title =ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size=20),
      axis.text.y = ggplot2::element_text(size=20),
      axis.title.y = ggplot2::element_text(size=30,vjust=0.5),
      #panel.grid.major =ggplot2::element_blank(),
      panel.grid.minor =ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour='black',size=2, fill=NA),
      axis.ticks = ggplot2::element_line(size=1,color='black'),
      plot.title = ggplot2::element_text(size=20)) #15x5   
  
  p2 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$Cphen, y=all_df$value, x=all_df$time)) + ggplot2::stat_summary(fun.y = mean, geom="line", size=1) +
    ggplot2::stat_summary(fun.y = mean, geom="point", shape=3, size=2) + ggplot2::scale_colour_manual(values=col) + 
    ggplot2::xlab("time [h]") + ggplot2::ylab("number organism") + ggplot2::ggtitle("Phenotyp growth curve with standard deviation") + 
    ggplot2::theme_bw(base_size = 30) +
    ggplot2::theme(#legend.position='none',
      legend.text= ggplot2::element_text(size=14),
      legend.key=ggplot2::element_blank(),
      legend.title =ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size=20),
      axis.text.y = ggplot2::element_text(size=20),
      axis.title.y = ggplot2::element_text(size=30,vjust=0.5),
      #panel.grid.major =ggplot2::element_blank(),
      panel.grid.minor =ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour='black',size=2, fill=NA),
      axis.ticks = ggplot2::element_line(size=1,color='black'),
      plot.title = ggplot2::element_text(size=20)) #15x5   
  
  p3 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$replc, group=all_df$replc,y=all_df$value, x=all_df$time)) + ggplot2::geom_line(size=1) + ggplot2::geom_point(size=2, shape=3) +
    ggplot2::scale_colour_manual(values=col) + 
    ggplot2::facet_wrap(~Cphen, nrow=4) +
    ggplot2::xlab("time [h]") + ggplot2::ylab("number organism") + ggplot2::ggtitle("Comparison of phenotype growth curves") + 
#    theme_bw(base_size = 30) +
    ggplot2::theme_classic(base_size = 30) +
    ggplot2::theme(#legend.position='none',
      legend.text= ggplot2::element_text(size=14),
      legend.key=ggplot2::element_blank(),
      legend.title =ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size=15),
      axis.text.y = ggplot2::element_text(size=20),
      axis.title.y = ggplot2::element_text(size=30,vjust=0.5),
      #panel.grid.major =ggplot2::element_blank(),
      #panel.grid.minor =ggplot2::element_blank(),
 #     panel.border = ggplot2::element_rect(colour='black',size=2, fill=NA),
      axis.ticks = ggplot2::element_line(size=1,color='black'),
      plot.title = ggplot2::element_text(size=20)) #15x5   

  p4 <- ggplot2::ggplot(all_df, ggplot2::aes(color=all_df$Cphen, group=all_df$Cphen,y=all_df$value, x=all_df$time)) + ggplot2::geom_line(size=1) + ggplot2::geom_point(size=2, shape=3) +
    ggplot2::scale_colour_manual(values=col) + 
    ggplot2::facet_wrap(~replc) +
    ggplot2::xlab("time [h]") + ggplot2::ylab("number organism") + ggplot2::ggtitle("Comparison of replicate growth curves") + 
    ggplot2::theme_bw(base_size = 30) +
    ggplot2::theme(#legend.position='none',
      legend.text= ggplot2::element_text(size=14),
      legend.key=ggplot2::element_blank(),
      legend.title =ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size=20),
      axis.text.y = ggplot2::element_text(size=20),
      axis.title.y = ggplot2::element_text(size=30,vjust=0.5),
      #panel.grid.major =ggplot2::element_blank(),
      panel.grid.minor =ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour='black',size=2, fill=NA),
      axis.ticks = ggplot2::element_line(size=1,color='black'),
      plot.title = ggplot2::element_text(size=20)) #15x5   
  
  if(ret_phengroups & cluster) return(simlist_fac)
  else(return(list(p1,p2,p3,p4)))
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
  all_df <- data.frame()
  for(i in seq_along(simlist)){
    object <- simlist[[i]]
    if(all(!is.null(time))) time_seq <- seq(time[1],time[2]) else time_seq <- seq_along(object@simlist)
    if(use_biomass){
      list <- lapply(time_seq, function(i){
        sapply(seq_along(object@specs), function(x){sum(object@simlist[[i]]$growth[which(object@simlist[[i]]$type==x)])})})
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
    abundances<-unlist(lapply(levels(all_df$species), function(s){
      mean(all_df$value[which(all_df$species==s)])}))
    names(abundances) <- levels(all_df$species)
    return(abundances)
  }else{
    q <- ggplot2::ggplot(all_df, ggplot2::aes(factor(all_df$species), all_df$value)) + ggplot2::geom_boxplot(ggplot2::aes(color=factor(all_df$species), fill=factor(all_df$species)), alpha = 0.2, outlier.size=1) + 
      ggplot2::scale_fill_manual(values=col) + ggplot2::scale_color_manual(values=col) + 
      ggplot2::theme(axis.text.x =ggplot2::element_blank(), legend.title=ggplot2::element_blank(),axis.title.x = ggplot2::element_blank(),axis.title.y = ggplot2::element_blank())
    return(q)
  }
}

#' @title Plot substance variations
#'
#' @description The function \code{plotSubVar} takes a list of simulations and return a barplot with most varying substances
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
    ggplot2::scale_shape_manual(values=1:nlevels(concmean$org)) +
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
plotSubUsage <- function(simlist, subs=list(), cutoff=1e-2, ret_data=FALSE){
  
  if(is(simlist, "Eval")) simlist <- list(simlist)
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
  
  if(!ret_data) df <- df[which(abs(df$mflux) > cutoff),,drop = FALSE] # do not drop if date is used further
  
  q1 <- ggplot2::ggplot(df, ggplot2::aes(x=df$time, y=df$mflux)) + ggplot2::geom_line(ggplot2::aes(col=df$spec), size=1) + ggplot2::facet_wrap(~df$sub, scales="free_y")+ ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  
  q2 <- ggplot2::ggplot(df, ggplot2::aes(factor(df$spec), df$mflux)) + ggplot2::geom_boxplot(ggplot2::aes(color=factor(df$spec), fill=factor(df$spec)), alpha=0.2) + 
    ggplot2::facet_wrap(~sub, scales="free_y") + ggplot2::theme(axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
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
#' @details Returns ggplot objects
plotSpecActivity <- function(simlist, subs=list(), var_nr=10, spec_list=NULL, ret_data=FALSE){
  
  if(is(simlist, "Eval")) simlist <- list(simlist)
  if(length(subs)==0) {subs_tocheck <- names(getVarSubs(simlist[[1]]))
  }else subs_tocheck <- subs
  if(length(spec_list)==0) spec_list <- names(simlist[[1]]@specs)
  
  df <- data.frame(spec=as.character(), sub=as.character(), mflux=as.numeric(), time=as.integer())
  
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
    mflux_var <- unlist(lapply(levels(df$sub), function(sub){
      stats::var(df[which(df$sub==sub),]$mflux)
    }))
    names(mflux_var) <- levels(df$sub)
    mflux_var <- sort(mflux_var, decreasing = TRUE)
    df <- df[which(df$sub %in% names(mflux_var)[1:var_nr]),]
  }
  
  q1 <- ggplot2::ggplot(df, ggplot2::aes(x=df$time, y=df$mflux)) + ggplot2::geom_line(ggplot2::aes(col=df$sub), size=1) + ggplot2::facet_wrap(~spec, scales="free_y") + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  
  q2 <- ggplot2::ggplot(df, ggplot2::aes(factor(df$sub), df$mflux)) + ggplot2::geom_boxplot(ggplot2::aes(color=factor(df$sub), fill=factor(df$sub)), alpha=0.2) +  ggplot2::facet_wrap(~df$spec, scales="free_y") +
    ggplot2::theme(axis.text.x =ggplot2::element_blank()) + ggplot2::xlab("") + ggplot2::ylab("mmol/(h*g_dw)")
  
  if(ret_data) return(df) else return(list(q1, q2))
}

