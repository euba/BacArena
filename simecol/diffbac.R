# Diffusion and movement

library(simecol)

n <- 100
m <- 100

s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_D_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")

substrat <- lapply(s, function(x, n, m){
  #matrix(runif(n*m), nrow=n, ncol=m)
  #matrix(c(runif(n), rep(0, n*m-n)), nrow=n, ncol=m)
  matrix(c(rep(1, n), rep(0, n*m-n)), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s


diff <- new("rwalkModel",
            main = function(time, init, parms) {
              x<- init$sub             
              # diffusion by mean over neighbourhood
              x <- lapply(x, function(x){
                anb <- eightneighbours(x)
                nb <- neighbours(x)
                mat <- (anb+x) / (nb+1)
                return(mat)
              })
              jpeg(paste("~/BacArena/plot", paste(letters[time-1], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
              image(x$M_co2_b, zlim=c(0,1), col=colorRampPalette(c("white", "black", "red"))(40))
              #par(mfrow=c(3,3))
              #for(i in 1:9){
              #  image(x[[i]], zlim=c(0,1), 
              #    col=colorRampPalette(c("white", "black", rainbow(14)[i]))(40),
              #    main=paste(names(x)[i], paste("step:", time)))
              #}
              #print(x$M_co2_b)
              x2 <- init$bac
              y <- x2
              for (i in 0:(parms$n-1)){
                for(j in 0:(parms$m-1)){
                  if(y[i+1,j+1] != 0){
                    a <- (i + round(runif(1,-1,1))) %% parms$n 
                    b <- (j + round(runif(1,-1,1))) %% parms$m 
                    if(x2[a+1,b+1] == 0){
                      x2[a+1,b+1] <- 1
                      x2[i+1,j+1] <- 0
                      #print(c(i,j,a,b))
                    }
                  }
                }
              }
              #jpeg(paste("~/BacArena/plot", paste(letters[time-1], ".jpeg", sep=""), sep=""), quality = 100, width=600, height=600)
              #image(x2, col=c("white", "black"))
              #print(sum(apply(x2, 1, sum)))
              dev.off()
              return(list(sub=x, bac=x2))
              #return(list(sub=x))
            },
            parms = list(n=n, m=m),
            times = c(from=1, to=27, by=1),
            init = list(sub=substrat, bac=matrix(round(runif(n*m, min=0, max=0.7)), nrow=n, ncol=m)),
            #init = list(sub=substrat, bac=matrix(c(rep(1,n), rep(0,n*m-n)), nrow=n, ncol=m)),
            solver = "iteration"
)

sim(diff)

