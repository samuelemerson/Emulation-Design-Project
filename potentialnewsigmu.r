greenavgvarformany = function(x,numpoints){
  a = matrix(x,nrow = numpoints,ncol = 2)
  varf = 0.60^2
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
greenoptimpoints.0 = optim(as.vector(latpoints),greenavgvarformany,numpoints=9,method="L-BFGS-B",lower=rep(-0.5,18),upper=rep(2*pi+0.5,18))
greenoavgpoints = matrix(greenoptimpoints.0$par,nrow=9,ncol=2)
greenoavgpoints

greenavoptvarrun1 = emufuncosxplussiny(0,0.60,2,greenoavgpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(greenavoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
greenavoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                        MoreArgs = list(z=1,points=greenoavgpoints))
greenavoptvarImat = matrix(greenavoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(greenavoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints,col="black",pch=19)})
length(which(greenavoptvarImat<3))
greenavoptvarIless = as.vector(which(greenavoptvarImat<3))
greenavoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                 y=seq(-0.5,(2*pi)+0.5,length.out=40))[greenavoptvarIless,]

greenavoptvarpointIless
greenavgvarformany.2 = function(x,numpoints,ogpoints,greenspace){
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
greenadpoi = greenavoptvarpointIless[sample(nrow(greenavoptvarpointIless),size=4),]
greenadpoi
greenminpointind = c(which(as.vector(greenavoptvarI)==sort(as.vector(greenavoptvarI))[1]),
                     which(as.vector(greenavoptvarI)==sort(as.vector(greenavoptvarI))[2]),
                     which(as.vector(greenavoptvarI)==sort(as.vector(greenavoptvarI))[3]),
                     which(as.vector(greenavoptvarI)==sort(as.vector(greenavoptvarI))[4]))
greenminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                       y=seq(-0.5,(2*pi)+0.5,length.out=40))[greenminpointind,]
greenminpoint
greenoptimpoints.1 = optim(c(greenminpoint[,1],greenminpoint[,2]),greenavgvarformany.2,numpoints=4,ogpoints=greenoavgpoints,greenspace=greenavoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
greenoptimpoints.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(greenoptimpoints.1$par,nrow=4,ncol=2))),
                 MoreArgs = list(z=1,points=greenoavgpoints)))
matrix(greenoptimpoints.1$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(greenavoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints,col="black",pch=19)
                 points(matrix(greenoptimpoints.1$par,nrow=4,ncol=2),col="black",pch=4)})

greenoavgpoints.1 = rbind(greenoavgpoints,matrix(greenoptimpoints.1$par,nrow=4,ncol=2))
greenoavgpoints.1
greenavoptvarrun2 = emufuncosxplussiny(0,0.6,2,greenoavgpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(greenavoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
greenavoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=greenoavgpoints.1))
greenavoptvarImat.1 = matrix(greenavoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(greenavoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints.1,col="black",pch=19)})

length(which(greenavoptvarImat.1<3))
greenavoptvarIless.1 = as.vector(which(greenavoptvarImat.1<3))
greenavoptvarpointIless.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                        y=seq(-0.5,(2*pi)+0.5,length.out=40))[greenavoptvarIless.1,]

greenavoptvarpointIless.1

#greenadpoi.1 = greenavoptvarpointIless.1[sample(nrow(greenavoptvarpointIless.1),size=2),]
#greenadpoi.1
greenminpointind.1 = c(which(as.vector(greenavoptvarI.1)==sort(as.vector(greenavoptvarI.1))[1]),
                       which(as.vector(greenavoptvarI.1)==sort(as.vector(greenavoptvarI.1))[3]))
greenminpoint.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                              y=seq(-0.5,(2*pi)+0.5,length.out=40))[greenminpointind.1,]
greenminpoint.1
greenoptimpoints.2 = optim(c(greenminpoint.1[,1],greenminpoint.1[,2]),greenavgvarformany.2,numpoints=2,ogpoints=greenoavgpoints.1,greenspace=greenavoptvarpointIless.1,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
greenoptimpoints.2
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(greenoptimpoints.2$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=greenoavgpoints.1)))
matrix(greenoptimpoints.2$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(greenavoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints.1,col="black",pch=19)
                 points(matrix(greenoptimpoints.2$par,nrow=2,ncol=2),col="black",pch=4)})

greenoavgpoints.2 = rbind(greenoavgpoints.1,matrix(greenoptimpoints.2$par,nrow=2,ncol=2))
greenoavgpoints.2
greenavoptvarrun3 = emufuncosxplussiny(0,0.6,2,greenoavgpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(greenavoptvarrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
greenavoptvarI.2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=greenoavgpoints.2))
greenavoptvarImat.2 = matrix(greenavoptvarI.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = greenavoptvarImat.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(greenavoptvarImat.2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(greenoavgpoints.2,col="black",pch=19)})
