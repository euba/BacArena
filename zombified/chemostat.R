source(file="R/Arena.R")
source(file="R/Eval.R")
medium <- c(15, 60)
names(medium) <- c("EX_pi(e)","EX_glc(e)")

bace = Bac(model=ecore, deathrate=0.3, duplirate=1.5, growthlimit=0.05, growtype="exponential")
arena = Arena(n=50, m=50, tstep=0.3, stir=T)
addOrg(arena, bace, amount=1000)
addSubs(arena)
changeSub(arena, smax=medium[1], mediac=names(medium)[1])
changeSub(arena, smax=medium[2], mediac=names(medium)[2])

time=20
inflow=0.1
outflow=0.063
product <- data.frame()
eval <- Eval(arena)

for(i in 1:time){
  print(i)
  arena <- getArena(simulate(arena, time=1))
  addEval(eval, arena)
  arena_size <- n(arena)*m(arena)
  # remove/add concentration to the medium
  media_subs <- media(arena)
  for(j in seq_along(names(media_subs))){
    sub <- media_subs[[j]]
    conc_matrix <- diffmat(sub)
    product[name(sub),i] <- sum(conc_matrix)*outflow
    newconc <- (sum(conc_matrix) - product[name(sub), i])/(arena_size)
    if(name(sub) %in% names(medium)){
      newconc <- newconc + (medium[name(sub)]*inflow)
    }
    changeSub(arena, smax=newconc, mediac=name(sub))
  }
  # remove bacteria from the medium
  neworgdat <- orgdat(arena)
  rem_org <- (floor(nrow(neworgdat)*outflow))
  if(rem_org!=0){
    neworgdat <- neworgdat[-(1:rem_org),]
  }
  changeOrg(arena, neworgdat)
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

