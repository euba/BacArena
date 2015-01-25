#function for converting Arena object into a position/type data.frame of all organisms involved

setGeneric("pop2dat", function(object){standardGeneric("pop2dat")})
setMethod("pop2dat", "Arena", function(object){
  tmp <- lapply(object@orglist, function(x){return(c(x@x, x@y, x@type))})
  tmp <- t(as.data.frame(tmp))
  popdat <- data.frame("x"=as.numeric(tmp[,1]), "y"=as.numeric(tmp[,2]), "type"=as.factor(tmp[,3]), row.names=1:nrow(tmp))
  return(popdat)
})

#function for converting Arena object into a position/object data.frame of all organisms involved

setGeneric("pop2dat2", function(object){standardGeneric("pop2dat2")})
setMethod("pop2dat2", "Arena", function(object){
  tmp <- lapply(object@orglist, function(x){return(c(x@x, x@y, x))})
  tmp <- t(as.data.frame(tmp))
  popdat <- data.frame("x"=as.numeric(tmp[,1]), "y"=as.numeric(tmp[,2]), "type"=as.factor(tmp[,3]), row.names=1:nrow(tmp))
  return(popdat)
})

#function for converting Arena object into a presence/absence matrix of all organisms involved

setGeneric("pop2mat", function(object){standardGeneric("pop2mat")})
setMethod("pop2mat", "Arena", function(object){
  popdat <- pop2dat(object)
  popmat <- matrix(0, object@n, object@m)
  apply(popdat, 1, function(x, types){
    popmat[as.numeric(x[1]), as.numeric(x[2])] <<- which(types==x[3])
  }, types=levels(popdat$type))
  return(popmat)
})


#function for converting Arena object into a presence/absence matrix of all individuals involved

setGeneric("pop2imat", function(object){standardGeneric("pop2imat")})
setMethod("pop2imat", "Arena", function(object){
  popdat <- pop2dat(object)
  popmat <- matrix(0, object@n, object@m)
  for(i in 1:nrow(popdat)){
    popmat[popdat[i,]$x, popdat[i,]$y] <- i
  }
  return(popmat)
})


#function for random walk (movement) of the whole Arena through the grid space

setGeneric("moveRand", function(object){standardGeneric("moveRand")})
setMethod("moveRand", "Arena", function(object){
  n <- object@n
  m <- object@m
  bmat <- pop2imat(object)
  bmatn <- matrix(NA, nrow=n+2, ncol=m+2) #define environment with boundary conditions
  bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
  bdat <- object@orglist
  for(i in seq_along(bdat)){
    bmatn[2:(n+1), 2:(m+1)] <- bmat #put the values into the environment
    ic = bdat[[i]]@x
    jc = bdat[[i]]@y
    neighbours <- c(bmatn[ic,jc], 
                    bmatn[ic+1,jc], 
                    bmatn[ic+2,jc], 
                    bmatn[ic+2,jc+1],
                    bmatn[ic+2,jc+2], 
                    bmatn[ic+1,jc+2],
                    bmatn[ic,jc+2],
                    bmatn[ic,jc+1])
    neighbours = ifelse(is.na(neighbours),1,neighbours)
    pos <- which(neighbours==0)
    if(length(pos) >= 1){
      if(length(pos)!=1){
        pos = sample(pos, 1)
      }
      switch(pos,
{bmat[ic-1,jc-1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic,jc-1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic+1,jc-1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic+1,jc] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic+1,jc+1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic,jc+1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic-1,jc+1] <- bmat[ic,jc]; bmat[ic,jc] <- 0},
{bmat[ic-1,jc] <- bmat[ic,jc]; bmat[ic,jc] <- 0})
    }else{
      if(length(pos) == 0){
        next
      }
    }
newpos <- which(bmat==i, arr.ind=T)
tmp <- object@occmat[[bdat[[i]]@x,bdat[[i]]@y]]
object@occmat[[bdat[[i]]@x,bdat[[i]]@y]] <- 0
bdat[[i]]@x <- newpos[1]
bdat[[i]]@y <- newpos[2]
object@occmat[[bdat[[i]]@x,bdat[[i]]@y]] <- tmp
  }
eval.parent(substitute(object@orglist <- bdat))
eval.parent(substitute(object@occmat <- object@occmat))
})