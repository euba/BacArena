########################################################################################################
###################################### GRID ############################################################
########################################################################################################

n <- 10
m <- 10
iter <- 100
bacs <- 1
smax <- 70


########################################################################################################
###################################### BACTERIA ########################################################
########################################################################################################

#if(!exists("Bcoli_sbml", mode="list"))source(file="ecoli_iAF1260.R")
#if(!exists("ecoli_sbml", mode="list"))source(file="ecoli.R")
#if(!exists("barkeri_sbml", mode="list")) source(file="barkeri.R")
#if(!exists("beijerinckii_sbml", mode="list")) source(file="beijerinckii.R")

source(file="barkeri.R")

#bac <- data.frame(x=round(runif(bacs, min=1, max=n)), y=round(runif(bacs, min=1, max=m)),
#                  type=rep("ecoli", bacs), growth=rep(1, bacs))
#bac <- bac[!duplicated(bac[,1:2]),]
#rownames(bac) <- 1:nrow(bac) #change indices in data.frame

#bac <- data.frame(x=round(n/2), y=round(m/2),type="ecoli", growth=1) # one cell in the centre
#bac <- rbind(data.frame(x=round(n), y=round(m),type="ecoli", growth=1),data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1))
bac <- data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1) # one cell in the centre
#bac <- data.frame(x=round(n/2), y=round(m/2),type="Bcoli", growth=1) # one cell in the centre
#bac <- data.frame(x=round(n/2), y=round(m/2),type="beijerinckii", growth=1) # one cell in the centre
#bac <- rbind(data.frame(x=round(n), y=round(m),type="beijerinckii", growth=1),data.frame(x=round(n/2), y=round(m/2),type="barkeri", growth=1))

# get a color for each bac
bac_color <- as.numeric(as.factor(levels(bac[,3])))
names(bac_color) <- levels(bac[,3])

########################################################################################################
###################################### SUBSTRATE #######################################################
########################################################################################################

s <- c("acetate","aketoglutarate", "co2", "ethanol", "formiate", "fumarate", "glucose", "h2o", "proton", "lactate","o2", "iphosphate", "pyruvate", "succinate", "h2", "methanol", "methane", "acetone", "butyrate", "butanol")

substrat <- lapply(s, function(x, n, m){
  #matrix(runif(n*m,min=0,max=100), nrow=n, ncol=m) # random substrate
  #matrix(c(rep(100, 2*n), rep(0, n*m-2*n)), nrow=n, ncol=m) # downstairs substrate
  #matrix(c(rep(0,(n*m-2*n)/2), rep(10,2*n), rep(0,(n*m-2*n)/2)), nrow=n, ncol=m) # substrate in the middle of our street ohooo
  matrix(smax,n,m) # homogen substrate distribution
}, n=n, m=m)
names(substrat) <- s
#substrat[["iphosphate"]] <- matrix(smax,n,m)
#substrat[["h2o"]] <- matrix(smax,n,m)
#substrat[["proton"]] <- matrix(smax,n,m)
substrat[["pyruvate"]] <- matrix(0,n,m)
substrat[["h2"]] <- matrix(0,n,m)
substrat[["co2"]] <- matrix(0,n,m)
#substrat[["methanol"]] <- matrix(0,n,m)
substrat[["methane"]] <- matrix(0,n,m)
substrat[["acetate"]] <- matrix(0,n,m)
#matr<-matrix(0, n, m)
#matr[round(n/2),round(m/2)]=smax
#substrat[["methanol"]] <- matrix(0,n,m) # one cell in the centre
