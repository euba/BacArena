########################################################################################################
###################################### GRID ############################################################
########################################################################################################

n <- 20
m <- 20
iter <- 1000
smax <- 70 # substrate start concentration
#seed <- 123 # reproduction of random variables
seed <- sample(1:9999,1)

########################################################################################################
###################################### BACTERIA ########################################################
########################################################################################################

#source(file="ecoli_iAF1260.R")
source(file="ecoli.R")
#source(file="barkeri.R")
source(file="beijerinckii.R")

#bac <- data.frame(x=round(n/2), y=round(m/2),type="ecoli", growth=1) # one cell in the centre
#bac <- rbind(data.frame(x=round(n), y=round(m),type="ecoli", growth=1),data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1))
#bac <- data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1) # one cell in the centre
#bac <- data.frame(x=round(n/2), y=round(m/2),type="Bcoli", growth=1) # one cell in the centre
#bac <- data.frame(x=round(n/2), y=round(m/2),type="beijerinckii", growth=1) # one cell in the centre
#bac <- rbind(data.frame(x=round(n), y=round(m),type="beijerinckii", growth=1),data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1))
#bac <- rbind(data.frame(x=round(n/2), y=round(m/2),type="beijerinckii", growth=1),data.frame(x=round(n/2+1), y=round(m/2+1),type="barkeri", growth=1))
#bac <- data.frame() # empty grid
bac <- rbind(data.frame(x=round(n/2), y=round(m/2),type="beijerinckii", growth=1),data.frame(x=round(n/2+1), y=round(m/2+1),type="Bcoli", growth=1))


########################################################################################################
###################################### SUBSTRATE #######################################################
########################################################################################################

s <- c("acetate","aketoglutarate", "co2", "ethanol", "formiate", "fumarate", "glucose", "h2o", "proton", "lactate","o2", "iphosphate", "pyruvate", "succinate", "h2", "methanol", "methane", "acetone", "butyrate", "butanol")

substrat <- lapply(s, function(x, n, m){
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # random substrate
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m) # downstairs substrate
  #matrix(c(rep(0,(n*m-2*n)/2), rep(smax,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m) # substrate in the middle of our street ohooo
  #matrix(smax,n,m) # homogen substrate distribution
  #matrix(c(smax, rep(0,m*n-1)), nrow=n, ncol=m) # one peak in top left
  matrix(0,n,m) # no substrate distribution
}, n=n, m=m)
names(substrat) <- s
substrat[["iphosphate"]] <- matrix(smax,n,m)
#substrat[["h2o"]] <- matrix(smax,n,m)
#substrat[["proton"]] <- matrix(smax,n,m)
substrat[["glucose"]] <- matrix(smax,n,m)
#substrat[["pyruvate"]] <- matrix(smax,n,m)
#substrat[["h2"]] <- matrix(smax,n,m)
#substrat[["co2"]] <- matrix(smax,n,m)
#substrat[["methanol"]] <- matrix(smax,n,m)
#substrat[["methane"]] <- matrix(0,n,m)
#substrat[["acetate"]] <- matrix(smax,n,m)
#substrat[["butanol"]] <- matrix(0,n,m)
#substrat[["ethanol"]] <- matrix(0,n,m)
#substrat[["aketoglutarate"]] <- matrix(0,n,m)
#substrat[["formiate"]] <- matrix(0,n,m)
#substrat[["o2"]] <- matrix(smax,n,m)
#substrat[["succinate"]] <- matrix(0,n,m)
#substrat[["acetone"]] <- matrix(0,n,m)
#substrat[["butyrate"]] <- matrix(0,n,m)
#substrat[["butanol"]] <- matrix(0,n,m)
#substrat[["lactate"]] <- matrix(0,n,m)
#substrat[["fumarate"]] <- matrix(0,n,m)
#matr<-matrix(0, n, m)
#matr[round(n/2),round(m/2)]=smax
#substrat[["methanol"]] <- matrix(0,n,m) # one cell in the centre
