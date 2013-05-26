library("simecol")

n <- 10 # grid rows 
m <- 10 # grid cols

s <- c("M_ac_b","M_akg_b", "M_co2_b", "M_etoh_b", "M_for_b", "M_fum_b", "M_glc_D_b", "M_h2o_b", "M_h_b", "M_lac_D_b","M_o2_b", "M_pi_b", "M_pyr_b", "M_succ_b")

substrat <- lapply(s, function(x, n, m){
  matrix(runif(n*m), nrow=n, ncol=m)
}, n=n, m=m)
names(substrat) <- s

#image(substrat$M_co2_b)
#substrat$M_co2_b

#plot(as.data.frame(substrat$M_co2_b))

diff <- new("rwalkModel",
            main = function(time, init, parms) {
              # matrix(round(runif(10*10)-0.4), nrow=10, ncol=10)
              x<- init
                  
              # really nice and really fast but unfortunately not feasible for this problem ..
              #x <- lapply(x, function(x){
              #  apply(x,c(1,2), function(x,m){
              #    print(class(m))
              #  },m=x)})
              #print(x$M_co2_b)
              
              # diffusion by mean over neighbourhood
              x <- lapply(x, function(x){
                anb <- eightneighbours(x)
                nb <- neighbours(x)
                mat <- (anb+x) / (nb+1)
                return(mat)
              })
              image(x$M_co2_b, zlim=c(0,1))
              #print(x$M_co2_b)
              return(x)
            },
            parms = list(substrat=substrat),
            times = c(from=1, to=10, by=1),
            init = substrat,
            solver = "iteration"
)

sim(diff)

plot(sim(diff))


