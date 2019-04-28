avoptvarIz0 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                   MoreArgs = list(z=0,points=oavgpoints))
avoptvarImatz0 = matrix(avoptvarIz0,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImatz0, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImatz0)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)})
length(which(avoptvarImatz0<3))
avoptvarIlessz0 = as.vector(which(avoptvarImatz0<3))
avoptvarpointIlessz0 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                 y=seq(-0.5,(2*pi)+0.5,length.out=40))[avoptvarIlessz0,]

avoptvarpointIlessz0
avgvarformany.2z0 = function(x,numpoints,ogpoints,greenspace){
  nwe = matrix(x,nrow = numpoints,ncol = 2)
  nweI = as.vector(mapply(impcosxplussiny,as.data.frame(t(nwe)),
                          MoreArgs = list(z=0,points=ogpoints)))
  nwebadI = sum(((nweI[nweI>3]-3)/10)^2)
  nweIbool = nweI>3
  a = rbind(ogpoints,nwe)
  varf = 0.25^2
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
#adpoiz0 = avoptvarpointIlessz0[sample(nrow(avoptvarpointIlessz0),size=4),]
#adpoiz0
minpointindz0 = c(which(as.vector(avoptvarIz0)==sort(as.vector(avoptvarIz0))[1]),
                  which(as.vector(avoptvarIz0)==sort(as.vector(avoptvarIz0))[2]),
                  which(as.vector(avoptvarIz0)==sort(as.vector(avoptvarIz0))[3]),
                  which(as.vector(avoptvarIz0)==sort(as.vector(avoptvarIz0))[4]))
minpointz0 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                       y=seq(-0.5,(2*pi)+0.5,length.out=40))[minpointindz0,]
minpointz0

optimpoints.1z0 = optim(c(minpointz0[,1],minpointz0[,2]),avgvarformany.2z0,numpoints=4,ogpoints=oavgpoints,greenspace=avoptvarpointIlessz0,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
optimpoints.1z0
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optimpoints.1z0$par,nrow=4,ncol=2))),
                 MoreArgs = list(z=0,points=oavgpoints)))
matrix(optimpoints.1z0$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImatz0, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImatz0)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(optimpoints.1z0$par,nrow=4,ncol=2),col="black",pch=4)})

oavgpoints.1z0 = rbind(oavgpoints,matrix(optimpoints.1z0$par,nrow=4,ncol=2))
oavgpoints.1z0
avoptvarrun2z0 = emufuncosxplussiny(0,0.25,2,oavgpoints.1z0)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun2z0[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun2z0[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.1z0,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
avoptvarI.1z0 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=0,points=oavgpoints.1z0))
avoptvarImat.1z0 = matrix(avoptvarI.1z0,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat.1z0, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat.1z0)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.1z0,col="black",pch=19)})

length(which(avoptvarImat.1z0<3))
avoptvarIless.1z0 = as.vector(which(avoptvarImat.1z0<3))
avoptvarpointIless.1z0 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                   y=seq(-0.5,(2*pi)+0.5,length.out=40))[avoptvarIless.1z0,]

avoptvarpointIless.1z0

minpointind.1z0 = c(which(as.vector(avoptvarI.1z0)==sort(as.vector(avoptvarI.1z0))[1]),
                    which(as.vector(avoptvarI.1z0)==sort(as.vector(avoptvarI.1z0))[2]))
minpoint.1z0 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                       y=seq(-0.5,(2*pi)+0.5,length.out=40))[minpointind.1z0,]
minpoint.1z0
optimpoints.2z0 = optim(c(minpoint.1z0[,1],minpoint.1z0[,2]),avgvarformany.2z0,numpoints=2,ogpoints=oavgpoints.1z0,greenspace=avoptvarpointIless.1z0,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optimpoints.2z0
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optimpoints.2z0$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=0,points=oavgpoints.1z0)))
matrix(optimpoints.2z0$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat.1z0, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat.1z0)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.1z0,col="black",pch=19)
                 points(matrix(optimpoints.2z0$par,nrow=2,ncol=2),col="black",pch=4)})

oavgpoints.2z0 = rbind(oavgpoints.1z0,matrix(optimpoints.2z0$par,nrow=2,ncol=2))
oavgpoints.2z0
avoptvarrun3z0 = emufuncosxplussiny(0,0.25,2,oavgpoints.2z0)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun3z0[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun3z0[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.2z0,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
avoptvarI.2z0 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=0,points=oavgpoints.2z0))
avoptvarImat.2z0 = matrix(avoptvarI.2z0,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat.2z0, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat.2z0)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints.2z0,col="black",pch=19)})
