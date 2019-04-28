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
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

expectedimprovement2D = function(x,B0,sigmau,theta,points){
  a = points
  D = cos(a[,1])+sin(a[,2])
  curmax = max(D)
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  varDw = varD*(1-(10^(-6))) + (10^(-6))*diag(length(points[,1]))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=theta))
    return(as.vector(g)*(1-(10^(-6))))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=theta))
    return(as.vector(i2)*(1-(10^(-6))))
  }
  ju = cbind(D,eD)
  sub = ju[,1]-ju[,2]
  svarD = solve(varDw)
  eDf = function(x){
    h = exf + covfD(x)%*%svarD%*%(sub)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  mu = eDf(x)
  var = varDf(x)
  if(var==0){
    accfun=0
  }
  else{
    z = (mu-curmax)/sqrt(var)
    accfun = (mu-curmax)*pnorm(z) + sqrt(var)*dnorm(z)
  }
  return(as.vector(-accfun))
}

sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                     y=seq(-0.5,(2*pi)+0.5,length.out=40))
sequemacc = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=oavgpoints))
sequemaccmat = matrix(sequemacc,40)
crp.rg <- colorRampPalette(c("purple3","darkblue","blue","cyan","green","yellow","seashell2"))
colors = crp.rg(23)
plot(seq(-0.5,(2*pi)+0.5,length.out=40),seq(-0.5,(2*pi)+0.5,length.out=40),cex=0.2,pch=16,col=colors)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat,
               col = colors,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)})
startpoint = sequem[which(sequemacc==min(sequemacc)),]
startpoint
possnewmax = optim(startpoint,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=oavgpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmax
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat,
               col = colors,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
newmaxpoints = rbind(oavgpoints,possnewax$par)
newmaxpoints
avoptvarrun1acc = emufuncosxplussiny(0,0.6,2,newmaxpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1acc[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1acc[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

sequemacc.1 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpoints))
sequemaccmat.1 = matrix(sequemacc.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.1,
               col = colors,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints,col="black",pch=19)})
startpoint.1 = sequem[which(sequemacc.1==min(sequemacc.1)),]
startpoint.1
possnewmax.1 = optim(startpoint.1,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmax.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints,col="black",pch=19)
                 points(matrix(possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
newmaxpoints.1 = rbind(newmaxpoints,possnewmax.1$par)
newmaxpoints.1
avoptvarrun1acc.1 = emufuncosxplussiny(0,0.6,2,newmaxpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1acc.1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1acc.1[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
sequemacc.2 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpoints.1))
sequemaccmat.2 = matrix(sequemacc.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.2,
               col = colors,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.1,col="black",pch=19)})
startpoint.2 = sequem[which(sequemacc.2==min(sequemacc.2)),]
startpoint.2
possnewmax.2 = optim(startpoint.2,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpoints.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmax.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.2,
               col = topo.colors(26,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.1,col="black",pch=19)
                 points(matrix(possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
newmaxpoints.2 = rbind(newmaxpoints.1,possnewmax.2$par)
newmaxpoints.2
avoptvarrun1acc.2 = emufuncosxplussiny(0,0.6,2,newmaxpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1acc.2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1acc.2[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

#######################################################################################

emufunminuscosxminussiny = function(B0,sigmau,theta,points){
  a = points
  D = -(cos(a[,1])+sin(a[,2]))
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  varDw = varD*(1-(10^(-6))) + (10^(-6))*diag(length(points[,1]))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=theta))
    return(as.vector(g)*(1-(10^(-6))))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=theta))
    return(as.vector(i2)*(1-(10^(-6))))
  }
  ju = cbind(D,eD)
  sub = ju[,1]-ju[,2]
  svarD = solve(varDw)
  eDf = function(x){
    h = exf + covfD(x)%*%svarD%*%(sub)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                       y=seq(-0.5,(2*pi)+0.5,length.out=40))
  sequem2 = mapply(eDf,as.data.frame(t(sequem)))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  sequem2mat = matrix(sequem2,40)
  sequem3mat = matrix(sequem3,40)
  return(list(sequem,sequem2mat,sequem3mat,a,D))
}
minusavoptvarrun1 = emufunminuscosxminussiny(0,0.6,2,oavgpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minusexpectedimprovement2D = function(x,B0,sigmau,theta,points){
  a = points
  D = -cos(a[,1])-sin(a[,2])
  curmax = max(D)
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  varDw = varD*(1-(10^(-6))) + (10^(-6))*diag(length(points[,1]))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=theta))
    return(as.vector(g)*(1-(10^(-6))))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=theta))
    return(as.vector(i2)*(1-(10^(-6))))
  }
  ju = cbind(D,eD)
  sub = ju[,1]-ju[,2]
  svarD = solve(varDw)
  eDf = function(x){
    h = exf + covfD(x)%*%svarD%*%(sub)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  mu = eDf(x)
  var = varDf(x)
  if(var==0){
    accfun=0
  }
  else{
    z = (mu-curmax)/sqrt(var)
    accfun = (mu-curmax)*pnorm(z) + sqrt(var)*dnorm(z)
  }
  return(as.vector(-accfun))
}

sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                     y=seq(-0.5,(2*pi)+0.5,length.out=40))
minussequemacc = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=oavgpoints))
minussequemaccmat = matrix(minussequemacc,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)})
minusstartpoint = oavgpoints[which(minusavoptvarrun1[[5]]==max(minusavoptvarrun1[[5]])),]
minusstartpoint
minuspossnewmax = optim(minusstartpoint,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=oavgpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(minuspossnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints = rbind(oavgpoints,minuspossnewmax$par)
minusnewmaxpoints
minusavoptvarrun1acc = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc.1 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints))
minussequemaccmat.1 = matrix(minussequemacc.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints,col="black",pch=19)})
minusstartpoint.1 = minusnewmaxpoints[which(minusavoptvarrun1acc[[5]]==max(minusavoptvarrun1acc[[5]])),]
minusstartpoint.1
minuspossnewmax.1 = optim(minusstartpoint.1,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints,col="black",pch=19)
                 points(matrix(minuspossnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints.1 = rbind(minusnewmaxpoints,minuspossnewmax.1$par)
minusnewmaxpoints.1
minusavoptvarrun1acc.1 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc.1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc.1[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc.2 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints.1))
minussequemaccmat.2 = matrix(minussequemacc.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat.2,
               col = topo.colors(26,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.1,col="black",pch=19)})
minusstartpoint.2 = minusnewmaxpoints.1[which(minusavoptvarrun1acc.1[[5]]==max(minusavoptvarrun1acc.1[[5]])),]
minusstartpoint.2
minuspossnewmax.2 = optim(minusstartpoint.2,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat.2,
               col = topo.colors(26,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.1,col="black",pch=19)
                 points(matrix(minuspossnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints.2 = rbind(minusnewmaxpoints.1,minuspossnewmax.2$par)
minusnewmaxpoints.2
minusavoptvarrun1acc.2 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc.2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc.2[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc.3 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints.2))
minussequemaccmat.3 = matrix(minussequemacc.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat.3,
               col = topo.colors(26,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.2,col="black",pch=19)})
minusstartpoint.3 = minusnewmaxpoints.2[which(minusavoptvarrun1acc.2[[5]]==max(minusavoptvarrun1acc.2[[5]])),]
minusstartpoint.3
minuspossnewmax.3 = optim(minusstartpoint.3,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat.3,
               col = topo.colors(26,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.2,col="black",pch=19)
                 points(matrix(minuspossnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints.3 = rbind(minusnewmaxpoints.2,minuspossnewmax.3$par)
minusnewmaxpoints.3
minusavoptvarrun1acc.3 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc.3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc.3[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints.3,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

##########################################################################################

minusthroptvarrun1 = emufunminuscosxminussiny(0,0.6,2,othrpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusthroptvarrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusthroptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othrpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc3 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=othrpoints))
minussequemaccmat3 = matrix(minussequemacc3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othrpoints,col="black",pch=19)})
minusstartpoint3 = othrpoints[which(minusthroptvarrun1[[5]]==max(minusthroptvarrun1[[5]])),]
minusstartpoint3
minuspossnewmax3 = optim(minusstartpoint3,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=othrpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othrpoints,col="black",pch=19)
                 points(matrix(minuspossnewmax3$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints3 = rbind(othrpoints,minuspossnewmax3$par)
minusnewmaxpoints3
minusavoptvarrun1acc3 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc3[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc3.1 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3))
minussequemaccmat3.1 = matrix(minussequemacc3.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3,col="black",pch=19)})
minusstartpoint3.1 = minusnewmaxpoints3[which(minusavoptvarrun1acc3[[5]]==max(minusavoptvarrun1acc3[[5]])),]
minusstartpoint3.1
minuspossnewmax3.1 = optim(minusstartpoint3.1,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax3.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3,col="black",pch=19)
                 points(matrix(minuspossnewmax3.1$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints3.1 = rbind(minusnewmaxpoints3,minuspossnewmax3.1$par)
minusnewmaxpoints3.1
minusavoptvarrun1acc3.1 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints3.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc3.1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc3.1[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc3.2 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.1))
minussequemaccmat3.2 = matrix(minussequemacc3.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.1,col="black",pch=19)})
minusstartpoint3.2 = minusnewmaxpoints3.1[which(minusavoptvarrun1acc3.1[[5]]==max(minusavoptvarrun1acc3.1[[5]])),]
minusstartpoint3.2
minuspossnewmax3.2 = optim(minusstartpoint3.2,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax3.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.1,col="black",pch=19)
                 points(matrix(minuspossnewmax3.2$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints3.2 = rbind(minusnewmaxpoints3.1,minuspossnewmax3.2$par)
minusnewmaxpoints3.2
minusavoptvarrun1acc3.2 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints3.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc3.2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc3.2[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc3.3 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.2))
minussequemaccmat3.3 = matrix(minussequemacc3.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.2,col="black",pch=19)})
minusstartpoint3.3 = minusnewmaxpoints3.2[which(minusavoptvarrun1acc3.2[[5]]==max(minusavoptvarrun1acc3.2[[5]])),]
minusstartpoint3.3
minuspossnewmax3.3 = optim(minusstartpoint3.3,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax3.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.2,col="black",pch=19)
                 points(matrix(minuspossnewmax3.3$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints3.3 = rbind(minusnewmaxpoints3.2,minuspossnewmax3.3$par)
minusnewmaxpoints3.3
minusavoptvarrun1acc3.3 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints3.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc3.3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc3.3[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.3,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc3.4 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.3))
minussequemaccmat3.4 = matrix(minussequemacc3.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.4,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.3,col="black",pch=19)})
minusstartpoint3.4 = minusnewmaxpoints3.3[which(minusavoptvarrun1acc3.3[[5]]==max(minusavoptvarrun1acc3.3[[5]])),]
minusstartpoint3.4
minuspossnewmax3.4 = optim(minusstartpoint3.4,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax3.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.4,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.3,col="black",pch=19)
                 points(matrix(minuspossnewmax3.4$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints3.4 = rbind(minusnewmaxpoints3.3,minuspossnewmax3.4$par)
minusnewmaxpoints3.4
minusavoptvarrun1acc3.4 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints3.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc3.4[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc3.4[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.4,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc3.5 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.4))
minussequemaccmat3.5 = matrix(minussequemacc3.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.5,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.4,col="black",pch=19)})
minusstartpoint3.5 = minusnewmaxpoints3.4[which(minusavoptvarrun1acc3.4[[5]]==max(minusavoptvarrun1acc3.4[[5]])),]
minusstartpoint3.5
minuspossnewmax3.5 = optim(minusstartpoint3.5,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints3.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax3.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat3.5,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.4,col="black",pch=19)
                 points(matrix(minuspossnewmax3.5$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints3.5 = rbind(minusnewmaxpoints3.4,minuspossnewmax3.5$par)
minusnewmaxpoints3.5
minusavoptvarrun1acc3.5 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints3.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc3.5[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc3.5[[3]]),
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints3.5,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

#############################################################################################

minustwooptvarrun1 = emufunminuscosxminussiny(0,0.6,2,otwopoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minustwooptvarrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minustwooptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=otwopoints))
minussequemaccmat2 = matrix(minussequemacc2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints,col="black",pch=19)})
minusstartpoint2 = otwopoints[which(minustwooptvarrun1[[5]]==max(minustwooptvarrun1[[5]])),]
minusstartpoint2
minuspossnewmax2 = optim(minusstartpoint2,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=otwopoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints,col="black",pch=19)
                 points(matrix(minuspossnewmax2$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2 = rbind(otwopoints,minuspossnewmax2$par)
minusnewmaxpoints2
minusavoptvarrun1acc2 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2.1 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2))
minussequemaccmat2.1 = matrix(minussequemacc2.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2,col="black",pch=19)})
minusstartpoint2.1 = minusnewmaxpoints2[which(minusavoptvarrun1acc2[[5]]==max(minusavoptvarrun1acc2[[5]])),]
minusstartpoint2.1
minuspossnewmax2.1 = optim(minusstartpoint2.1,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2,col="black",pch=19)
                 points(matrix(minuspossnewmax2.1$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2.1 = rbind(minusnewmaxpoints2,minuspossnewmax2.1$par)
minusnewmaxpoints2.1
minusavoptvarrun1acc2.1 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2.1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2.1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2.2 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.1))
minussequemaccmat2.2 = matrix(minussequemacc2.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.1,col="black",pch=19)})
minusstartpoint2.2 = minusnewmaxpoints2.1[which(minusavoptvarrun1acc2.1[[5]]==max(minusavoptvarrun1acc2.1[[5]])),]
minusstartpoint2.2
minuspossnewmax2.2 = optim(minusstartpoint2.2,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.1,col="black",pch=19)
                 points(matrix(minuspossnewmax2.2$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2.2 = rbind(minusnewmaxpoints2.1,minuspossnewmax2.2$par)
minusnewmaxpoints2.2
minusavoptvarrun1acc2.2 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2.2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2.2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2.3 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.2))
minussequemaccmat2.3 = matrix(minussequemacc2.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.2,col="black",pch=19)})
minusstartpoint2.3 = minusnewmaxpoints2.2[which(minusavoptvarrun1acc2.2[[5]]==max(minusavoptvarrun1acc2.2[[5]])),]
minusstartpoint2.3
minuspossnewmax2.3 = optim(minusstartpoint2.3,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.2,col="black",pch=19)
                 points(matrix(minuspossnewmax2.3$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2.3 = rbind(minusnewmaxpoints2.2,minuspossnewmax2.3$par)
minusnewmaxpoints2.3
minusavoptvarrun1acc2.3 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2.3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2.3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.3,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2.4 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.3))
minussequemaccmat2.4 = matrix(minussequemacc2.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.3,col="black",pch=19)})
minusstartpoint2.4 = minusnewmaxpoints2.3[which(minusavoptvarrun1acc2.3[[5]]==max(minusavoptvarrun1acc2.3[[5]])),]
minusstartpoint2.4
minuspossnewmax2.4 = optim(minusstartpoint2.4,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.3,col="black",pch=19)
                 points(matrix(minuspossnewmax2.4$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2.4 = rbind(minusnewmaxpoints2.3,minuspossnewmax2.4$par)
minusnewmaxpoints2.4
minusavoptvarrun1acc2.4 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2.4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2.4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.4,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2.5 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.4))
minussequemaccmat2.5 = matrix(minussequemacc2.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.4,col="black",pch=19)})
minusstartpoint2.5 = minusnewmaxpoints2.4[which(minusavoptvarrun1acc2.4[[5]]==max(minusavoptvarrun1acc2.4[[5]])),]
minusstartpoint2.5
minuspossnewmax2.5 = optim(minusstartpoint2.5,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.4,col="black",pch=19)
                 points(matrix(minuspossnewmax2.5$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2.5 = rbind(minusnewmaxpoints2.4,minuspossnewmax2.5$par)
minusnewmaxpoints2.5
minusavoptvarrun1acc2.5 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2.5[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2.5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.5,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
minussequemacc2.6 = mapply(minusexpectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.5))
minussequemaccmat2.6 = matrix(minussequemacc2.6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.5,col="black",pch=19)})
minusstartpoint2.6 = minusnewmaxpoints2.5[which(minusavoptvarrun1acc2.5[[5]]==max(minusavoptvarrun1acc2.5[[5]])),]
minusstartpoint2.6
minuspossnewmax2.6 = optim(minusstartpoint2.6,minusexpectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=minusnewmaxpoints2.5,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
minuspossnewmax2.6
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -minussequemaccmat2.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.5,col="black",pch=19)
                 points(matrix(minuspossnewmax2.6$par,nrow=1,ncol=2),col="black",pch=4)})
minusnewmaxpoints2.6 = rbind(minusnewmaxpoints2.5,minuspossnewmax2.6$par)
minusnewmaxpoints2.6
minusavoptvarrun1acc2.6 = emufunminuscosxminussiny(0,0.6,2,minusnewmaxpoints2.6)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = minusavoptvarrun1acc2.6[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(minusavoptvarrun1acc2.6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(minusnewmaxpoints2.6,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
