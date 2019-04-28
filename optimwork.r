X.sum <- matrix(c(10,-5,-5,20),2)
X.sum
Z.sum <- matrix(sample(1:1000), ncol=2)
Z.sum
rowSums((Z.sum %*% X.sum) * Z.sum)
diag(Z.sum %*% X.sum %*% t(Z.sum))

fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}
optim(c(-1.2,1), fr)

#######################################################################
maxqvarexpand = function(a){
  varf = 0.25^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=20),
                       y=seq(-0.5,(2*pi)+0.5,length.out=20))
  vari = rep(varf,times=400)
  sequem3.1 = mapply(covfD,as.data.frame(t(sequem)))
  sequem3.2 = matrix(sequem3.1,400,9,byrow = TRUE)
  diagcal = rowSums((sequem3.2 %*% svarD) * sequem3.2)
  resdiag = vari - diagcal
  maxres = max(resdiag)
  return(maxres)
}

maxqvarexpand(latpoints)
optim(par = latpoints[1,],maxqvarexpand, method = "Nelder-Mead")
matrix(as.vector(latpoints),9,2)==latpoints
length(latpoints[1,][,1])
mkoy = latpoints[1,]

#############################################################################
maxoverspace = function(x,poin){
  varf = 0.25^2
  b = as.matrix(dist(poin,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(poin)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  varres = varDf(x)
  return(varres)
}

fdejong2 <- function (x) {
  return (cos(x[1])+sin(x[2])-cos(x[3]))
}
optim(c(-0.5,-0.5,-0.5),fdejong2)
dist(c(1,1))
######################################################################

avgvarforone = function(a){
  varf = 0.25^2
  svarD = 1/(varf)
  varDf = function(x){
    j = varf - covfun2d(x,a,varf,2)*svarD*covfun2d(x,a,varf,2)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=20),
                       y=seq(-0.5,(2*pi)+0.5,length.out=20))
  sequem3.5 = mapply(varDf,as.data.frame(t(sequem)))
  average = mean(sequem3.5)
  return(average)
}

optim(c(0,0),avgvarforone,method="L-BFGS-B",lower=c(-0.5,-0.5),upper=c(2*pi,2*pi))

######################################################################

avgvarformany = function(x,numpoints){
  a = matrix(x,nrow = numpoints,ncol = 2)
  varf = 0.6^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
    return(as.vector(i2))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                       y=seq(-0.5,(2*pi)+0.5,length.out=40))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  average = mean(sequem3)
  return(average)
}

as.vector(latpoints)
latpoints
matrix(as.vector(latpoints),nrow=9,ncol=2)
optimpoints.0 = optim(as.vector(latpoints),avgvarformany,numpoints=9,method="L-BFGS-B",lower=rep(-0.5,18),upper=rep(2*pi+0.5,18))
oavgpoints = matrix(optimpoints.0$par,nrow=9,ncol=2)
oavgpoints

avoptvarrun1 = emufuncosxplussiny(0,0.6,2,oavgpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
avoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=oavgpoints))
avoptvarImat = matrix(avoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)})
length(which(avoptvarImat<3))
avoptvarIless = as.vector(which(avoptvarImat<3))
avoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                              y=seq(-0.5,(2*pi)+0.5,length.out=40))[avoptvarIless,]

avoptvarpointIless
avgvarformany.2 = function(x,numpoints,ogpoints,greenspace){
  nwe = matrix(x,nrow = numpoints,ncol = 2)
  nweI = as.vector(mapply(impcosxplussiny,as.data.frame(t(nwe)),
                          MoreArgs = list(z=1,points=ogpoints)))
  nwebadI = sum(((nweI[nweI>3]-3)/10)^2)
  nweIbool = nweI>3
  a = rbind(ogpoints,nwe)
  varf = 0.6^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
    return(as.vector(i2))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = greenspace
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  average = mean(sequem3)
  if(length(which(nweIbool))>0){
    return(average + nwebadI)
  } else{
    return(average)
  }
}
#c(avoptvarpointIless[1:4,1],avoptvarpointIless[1:4,2])
adpoi = avoptvarpointIless[sample(nrow(avoptvarpointIless),size=4),]
adpoi
minpointind = c(which(as.vector(avoptvarI)==sort(as.vector(avoptvarI))[1]),
                which(as.vector(avoptvarI)==sort(as.vector(avoptvarI))[2]),
                which(as.vector(avoptvarI)==sort(as.vector(avoptvarI))[3]),
                which(as.vector(avoptvarI)==sort(as.vector(avoptvarI))[4]))
minpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                         y=seq(-0.5,(2*pi)+0.5,length.out=40))[minpointind,]
minpoint
optimpoints.1 = optim(c(minpoint[,1],minpoint[,2]),avgvarformany.2,numpoints=4,ogpoints=oavgpoints,greenspace=avoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
optimpoints.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optimpoints.1$par,nrow=4,ncol=2))),
                MoreArgs = list(z=1,points=oavgpoints)))
matrix(optimpoints.1$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(optimpoints.1$par,nrow=4,ncol=2),col="black",pch=4)})

oavgpoints.1 = rbind(oavgpoints,matrix(optimpoints.1$par,nrow=4,ncol=2))
oavgpoints.1
avoptvarrun2 = emufuncosxplussiny(0,0.6,2,oavgpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
avoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                   MoreArgs = list(z=1,points=oavgpoints.1))
avoptvarImat.1 = matrix(avoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.1,col="black",pch=19)})

length(which(avoptvarImat.1<3))
avoptvarIless.1 = as.vector(which(avoptvarImat.1<3))
avoptvarpointIless.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                 y=seq(-0.5,(2*pi)+0.5,length.out=40))[avoptvarIless.1,]

avoptvarpointIless.1
#adpoi.1 = avoptvarpointIless.1[sample(nrow(avoptvarpointIless.1),size=2),]
#adpoi.1
minpointind.1 = c(which(as.vector(avoptvarI.1)==sort(as.vector(avoptvarI.1))[1]),
                  which(as.vector(avoptvarI.1)==sort(as.vector(avoptvarI.1))[4]))
minpoint.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                         y=seq(-0.5,(2*pi)+0.5,length.out=40))[minpointind.1,]
minpoint.1
optimpoints.2 = optim(c(minpoint.1[,1],minpoint.1[,2]),avgvarformany.2,numpoints=2,ogpoints=oavgpoints.1,greenspace=avoptvarpointIless.1,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optimpoints.2
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optimpoints.2$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=oavgpoints.1)))
matrix(optimpoints.2$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.1,col="black",pch=19)
                 points(matrix(optimpoints.2$par,nrow=2,ncol=2),col="black",pch=4)})

oavgpoints.2 = rbind(oavgpoints.1,matrix(optimpoints.2$par,nrow=2,ncol=2))
oavgpoints.2
avoptvarrun3 = emufuncosxplussiny(0,0.6,2,oavgpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
avoptvarI.2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=oavgpoints.2))
avoptvarImat.2 = matrix(avoptvarI.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat.2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.2,col="black",pch=19)})

#############################################################################################

quantvarformany = function(x){
  a = matrix(x,nrow = 9,ncol = 2)
  varf = 0.6^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
    return(as.vector(i2))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=20),
                       y=seq(-0.5,(2*pi)+0.5,length.out=20))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  quant = as.vector(quantile(sequem3,0.55))
  return(quant)
}

qoptimpoints.0 = optim(as.vector(latpoints),quantvarformany,method="L-BFGS-B",lower=rep(-0.5,18),upper=rep(2*pi+0.5,18))
qopoints = matrix(qoptimpoints.0$par,nrow=9,ncol=2)
qopoints
qoptimpoints.0

qoptvarrun1 = emufuncosxplussiny(0,0.6,2,qopoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(qoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qopoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
qoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                   MoreArgs = list(z=1,points=qopoints))
qoptvarImat = matrix(qoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(qoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qopoints,col="black",pch=19)})
length(which(qoptvarImat<3))
qoptvarIless = as.vector(which(qoptvarImat<3))
qoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                 y=seq(-0.5,(2*pi)+0.5,length.out=40))[qoptvarIless,]

qoptvarpointIless
quantvarformany.2 = function(x){
  nwe = matrix(x,nrow = 4,ncol = 2)
  nweI = as.vector(mapply(impcosxplussiny,as.data.frame(t(nwe)),
                          MoreArgs = list(z=1,points=qopoints)))
  nweIbool = nweI>3
  if(length(which(nweIbool))>0){
    return(100)
  } else{
    a = rbind(qopoints,nwe)
    varf = 0.6^2
    b = as.matrix(dist(a,diag=TRUE))
    varD = varf*exp(-(b^2)/(2^2))
    covfD = function(x){
      g = mapply(covfun2d,as.data.frame(t(a)),
                 MoreArgs = list(x=x,su=varf,t=2))
      return(as.vector(g))
    }
    covDf = function(x){
      i2 = mapply(covfun2d,as.data.frame(t(a)),
                  MoreArgs = list(y=x,su=varf,t=2))
      return(as.vector(i2))
    }
    svarD = solve(varD)
    varDf = function(x){
      j = varf - covfD(x)%*%svarD%*%covDf(x)
      return(j)
    }
    sequem = qoptvarpointIless
    sequem3 = mapply(varDf,as.data.frame(t(sequem)))
    quant = as.vector(quantile(sequem3,0.95))
    return(quant)
  }
}
qadpoi = qoptvarpointIless[sample(nrow(qoptvarpointIless),size=4),]
qadpoi[,1]
qoptimpoints.1 = optim(c(qadpoi[,1],qadpoi[,2]),quantvarformany.2,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
qoptimpoints.1
matrix(qoptimpoints.1$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(qoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qopoints,col="black",pch=19)
                 points(matrix(qoptimpoints.1$par,nrow=4,ncol=2),col="black",pch=4)})

##############################################################################################

maxvarformany = function(x){
  a = matrix(x,nrow = 9,ncol = 2)
  varf = 0.6^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
    return(as.vector(i2))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=20),
                       y=seq(-0.5,(2*pi)+0.5,length.out=20))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  max = max(sequem3)
  return(max)
}

moptimpoints.0 = optim(c(0,0,0,pi,pi,pi,2*pi,2*pi,2*pi,0,pi,2*pi,0,pi,2*pi,0,pi,2*pi),maxvarformany,method="L-BFGS-B",lower=rep(-0.5,18),upper=rep(2*pi+0.5,18))
mopoints = matrix(moptimpoints.0$par,nrow=9,ncol=2)
mopoints
moptimpoints.0

moptvarrun1 = emufuncosxplussiny(0,0.6,2,mopoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = moptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(moptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(mopoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
moptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                  MoreArgs = list(z=1,points=mopoints))
moptvarImat = matrix(moptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = moptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(moptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(mopoints,col="black",pch=19)})
length(which(moptvarImat<3))
moptvarIless = as.vector(which(moptvarImat<3))
moptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                y=seq(-0.5,(2*pi)+0.5,length.out=40))[moptvarIless,]

moptvarpointIless
maxvarformany.2 = function(x){
  nwe = matrix(x,nrow = 4,ncol = 2)
  nweI = as.vector(mapply(impcosxplussiny,as.data.frame(t(nwe)),
                          MoreArgs = list(z=1,points=mopoints)))
  nwebadI = sum(((nweI[nweI>3]-3)/10)^2)
  nweIbool = nweI>3
  a = rbind(mopoints,nwe)
  varf = 0.6^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
    return(as.vector(i2))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = moptvarpointIless
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  max = max(sequem3)
  if(length(which(nweIbool))>0){
    return(max + nwebadI)
  } else{
    return(max)
  }
}
madpoi = moptvarpointIless[sample(nrow(moptvarpointIless),size=4),]
madpoi
moptimpoints.1 = optim(c(madpoi[,1],madpoi[,2]),maxvarformany.2,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
moptimpoints.1
matrix(moptimpoints.1$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = moptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(moptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(mopoints,col="black",pch=19)
                 points(matrix(moptimpoints.1$par,nrow=4,ncol=2),col="black",pch=4)})




##floor is very high...bad! need a smoother penalty
##first need to scale variance down to implausibility...then choose smooth penalty, like imp^2 - 3 could work
#add nugget through stuff on board
#for variance problem now do it with the 95th quantile....see how different quantiles work
##also try 7 points...see where it puts it for average
#just add bullet points for latex report
#truncation error
#for one not using optim...add implausibility and then reject it if the imp is low