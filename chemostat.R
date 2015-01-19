source(file="class/Arena.R")
source(file="class/Eval.R")
medium <- c(15, 60)
names(medium) <- c("EX_pi(e)","EX_glc(e)")

bace = Bac(model=ecore, deathrate=0.3, duplirate=1.5, growthlimit=0.05, growtype="exponential")
arena = Arena(n=50, m=50, tstep=0.3, stir=T)
addOrg(arena, bace, amount=1000)
addSubs(arena)
changeSub(arena, smax=medium[1], mediac=names(medium)[1])
changeSub(arena, smax=medium[2], mediac=names(medium)[2])

time=1000
inflow=0.1
outflow=0.063
product <- data.frame()
eval <- Eval(arena)

for(i in 1:time){
  print(i)
  arena <- getArena(simulate(arena, time=1))
  addEval(eval, arena)
  # remove/add concentration to the medium
  for(j in seq_along(names(arena@media))){
    sub <- arena@media[[j]]
    product[sub@name,i] <- sum(sub@diffmat)*outflow
    newconc <- (sum(sub@diffmat) - product[sub@name, i])/(arena@n*arena@m)
    if(sub@name %in% names(medium)){
      newconc <- newconc + (medium[sub@name]*inflow)
    }
    changeSub(arena, smax=newconc, mediac=sub@name)
  }
  # remove bacteria from the medium
  rem_org <- (floor(nrow(arena@orgdat)*outflow))
  if(rem_org!=0){
    arena@orgdat <- arena@orgdat[-(1:rem_org),]
  }
  arena@occmat <- Matrix(dat2mat(arena), sparse=T)
}


#evalArena(eval, plot_items=c('population', "EX_glc(e)"), phencol=T, retdata=F)
#plotCurves(eval, remove=T, retdata = F)

#dev.off()

#rm <- c(which(apply(product, 1, sum)==0), which(rownames(product) %in% c("EX_glc(e)","EX_pi(e)")))
#matplot(t(product[-rm,]), type = c("b"), pch=1) #plot
#matplot(t(product), type = c("l"), col=1:nrow(product)) #plot

times=(1:ncol(product))*arena@tstep
plot(0, 0, xlim=c(0,max(times)), ylim=c(0,max(product)),
     type='n', xlab='time in h', ylab='concentration in mmol',
     main='Product concentrations')
for(i in 1:nrow(product)){
  lines(times, product[i,], col=i)
}
legend('right',legend=rownames(product),col=1:nrow(product),cex=0.6/log10(nrow(product)+1),lwd=1)
plotCurves(eval, remove=T, retdata = F)

