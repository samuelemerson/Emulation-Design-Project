minmaxfun = function(x){
  return(0.1*(cos(x[1])+sin(x[2]))*(x[1]+x[2]))
}
minmaxfunforlooks = function(x,y){
  return(0.1*(cos(x)+sin(y))*(x+y))
}
x <- seq(0, 3*pi, length= 100)
y <- x
z <- outer(x, y, minmaxfunforlooks)
z[is.na(z)] <- 1

require(lattice)
wireframe(z, drape=T, col.regions=rainbow(100))
op <- par(bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

filled.contour(x,y,z,
               color.palette = colorRampPalette(c('green','yellow','red')))

sequemminmax.2 = expand.grid(x=seq(-0.5,(3*pi)+0.5,length.out=100),
                             y=seq(-0.5,(3*pi)+0.5,length.out=100))
sequemfminmax = mapply(minmaxfun,as.data.frame(t(sequemminmax.2)))
sequemfminmaxmat = matrix(sequemfminmax,100)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=100),
               y = seq(-0.5,(3*pi)+0.5,length.out=100), 
               z = sequemfminmaxmat,
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Real Function ", 0.1*(cos(x[1])+sin(x[2]))*(x[1]+x[2]))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
#################################################################################

emumaxminfun = function(B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(minmaxfun,as.data.frame(t(a))))
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
  sequem = expand.grid(x=seq(-0.5,(3*pi)+0.5,length.out=40),
                       y=seq(-0.5,(3*pi)+0.5,length.out=40))
  sequem2 = mapply(eDf,as.data.frame(t(sequem)))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  sequem2mat = matrix(sequem2,40)
  sequem3mat = matrix(sequem3,40)
  return(list(sequem,sequem2mat,sequem3mat,a,D))
}

minmaxrun1 = emumaxminfun(0,1.8,3,minmaxpoints)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmaxpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmaxpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
expectedimprovement2Dminmax = function(x,B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(minmaxfun,as.data.frame(t(a))))
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

sequemminmax = expand.grid(x=seq(-0.5,(3*pi)+0.5,length.out=40),
                           y=seq(-0.5,(3*pi)+0.5,length.out=40))
sequemaccminmax = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=1.8,theta=3,points=minmaxpoints))
sequemaccmatminmax = matrix(sequemaccminmax,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmaxpoints,col="black",pch=19)})
startpointminmax = sequemminmax[which(sequemaccminmax==min(sequemaccminmax)),]
startpointminmax
possnewmaxminmax = optim(startpointminmax,expectedimprovement2Dminmax,B0=0,sigmau=1.85,theta=3,points=minmaxpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmaxminmax
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmaxpoints,col="black",pch=19)
                 points(matrix(possnewmaxminmax$par,nrow=1,ncol=2),col="black",pch=4)})
newminmaxpoints = rbind(minmaxpoints,possnewmaxminmax$par)
newminmaxpoints
minmaxrun2 = emumaxminfun(0,1.8,3,newminmaxpoints)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmax.2 = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=1.8,theta=3,points=newminmaxpoints))
sequemaccmatminmax.2 = matrix(sequemaccminmax.2,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints,col="black",pch=19)})
startpointminmax.2 = sequemminmax[which(sequemaccminmax.2==min(sequemaccminmax.2)),]
startpointminmax.2
possnewmaxminmax.2 = optim(startpointminmax.2,expectedimprovement2Dminmax,B0=0,sigmau=1.8,theta=3,points=newminmaxpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmaxminmax.2
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints,col="black",pch=19)
                 points(matrix(possnewmaxminmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
newminmaxpoints.2 = rbind(newminmaxpoints,possnewmaxminmax.2$par)
newminmaxpoints.2
minmaxrun3 = emumaxminfun(0,1.8,3,newminmaxpoints.2)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmax.3 = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=1.8,theta=3,points=newminmaxpoints.2))
sequemaccmatminmax.3 = matrix(sequemaccminmax.3,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax.3,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints.2,col="black",pch=19)})
startpointminmax.3 = sequemminmax[which(sequemaccminmax.3==min(sequemaccminmax.3)),]
startpointminmax.3
possnewmaxminmax.3 = optim(startpointminmax.3,expectedimprovement2Dminmax,B0=0,sigmau=1.8,theta=3,points=newminmaxpoints.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmaxminmax.3
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax.3,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints.2,col="black",pch=19)
                 points(matrix(possnewmaxminmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
newminmaxpoints.3 = rbind(newminmaxpoints.2,possnewmaxminmax.3$par)
newminmaxpoints.3
minmaxrun4 = emumaxminfun(0,1.8,3,newminmaxpoints.3)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmax.4 = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=1.8,theta=3,points=newminmaxpoints.3))
sequemaccmatminmax.4 = matrix(sequemaccminmax.4,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmax.4,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpoints.3,col="black",pch=19)})
##################################################################################

expectedimprovement2DMC4 = function(x,B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(minmaxfun,as.data.frame(t(a))))
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
  I = function(e){
    return((max(0,as.vector(mu)+(sqrt(as.vector(var))*e)-curmax))^4)
  }
  return(-sum(sapply(standnorm,I))/10000)
}

filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(minmaxpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmaxpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmaxfourth = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=1.8,theta=3,points=minmaxpoints))
sequemaccmatminmaxfourth = matrix(sequemaccminmaxfourth,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmaxpoints,col="black",pch=19)})
startpointminmaxfourth = sequemminmax[which(sequemaccminmaxfourth==min(sequemaccminmaxfourth)),]
startpointminmaxfourth

possnewmaxminmaxfourth = optim(startpointminmaxfourth,expectedimprovement2DMC4,B0=0,sigmau=1.8,theta=3,points=minmaxpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmaxminmaxfourth
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmaxpoints,col="black",pch=19)
                 points(matrix(possnewmaxminmaxfourth$par,nrow=1,ncol=2),col="black",pch=4)})
newminmaxpointsfourth = rbind(minmaxpoints,possnewmaxminmaxfourth$par)
newminmaxpointsfourth
minmaxrun2fourth = emumaxminfun(0,0.6,1.85,newminmaxpointsfourth)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun2fourth[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpointsfourth,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun2fourth[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpointsfourth,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmaxfourth.2 = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmaxpointsfourth))
sequemaccmatminmaxfourth.2 = matrix(sequemaccminmaxfourth.2,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth.2,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newminmaxpointsfourth,col="black",pch=19)})
startpointminmaxfourth.2 = newminmaxpointsfourth[which(minmaxrun2fourth[[5]]==max(minmaxrun2fourth[[5]])),]
startpointminmaxfourth.2
possnewmaxminmaxfourth.2 = optim(startpointminmaxfourth.2,optimforsmallexpimp,B0=0,sigmau=0.6,theta=1.85,points=newminmaxpointsfourth,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmaxminmaxfourth.2
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth.2,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newminmaxpointsfourth,col="black",pch=19)
                 points(matrix(possnewmaxminmaxfourth.2$par,nrow=1,ncol=2),col="black",pch=4)})
newminmaxpointsfourth.2 = rbind(newminmaxpointsfourth,possnewmaxminmaxfourth.2$par)
newminmaxpointsfourth.2
minmaxrun3fourth = emumaxminfun(0,0.6,1.85,newminmaxpointsfourth.2)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun3fourth[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun3fourth[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpointsfourth.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmaxfourth.3 = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmaxpointsfourth.2))
sequemaccmatminmaxfourth.3 = matrix(sequemaccminmaxfourth.3,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpointsfourth.2,col="black",pch=19)})
startpointminmaxfourth.3 = newminmaxpointsfourth.2[which(minmaxrun3fourth[[5]]==max(minmaxrun3fourth[[5]])),]
startpointminmaxfourth.3
optimforsmallexpimp.4 = function(x,B0,sigmau,theta,points){
  return(1000000000*expectedimprovement2DMC4(x,B0,sigmau,theta,points))
}
possnewmaxminmaxfourth.3 = optim(startpointminmaxfourth.3,optimforsmallexpimp,B0=0,sigmau=0.6,theta=1.85,points=newminmaxpointsfourth.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmaxminmaxfourth.3
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newminmaxpointsfourth.2,col="black",pch=19)
                 points(matrix(possnewmaxminmaxfourth.3$par,nrow=1,ncol=2),col="black",pch=4)})
newminmaxpointsfourth.3 = rbind(newminmaxpointsfourth.2,possnewmaxminmaxfourth.3$par)
newminmaxpointsfourth.3
minmaxrun4fourth = emumaxminfun(0,0.6,1.85,newminmaxpointsfourth.3)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmaxrun4fourth[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmaxrun4fourth[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpointsfourth.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccminmaxfourth.4 = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmaxpointsfourth.3))
sequemaccmatminmaxfourth.4 = matrix(sequemaccminmaxfourth.4,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequemaccmatminmaxfourth.4,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmaxpointsfourth.3,col="black",pch=19)})

##############################################################################################################

#minmax3points = as.matrix(3*pi*maximinLHS(3,2,optimize.on = "grid"))
minmax3points
minmax3run1 = emumaxminfun(0,0.6,1.85,minmax3points)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmax3points,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

sequem3accminmax = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=minmax3points))
sequem3accmatminmax = matrix(sequem3accminmax,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmax3points,col="black",pch=19)})
startpoint3minmax = minmax3points[which(minmax3run1[[5]]==max(minmax3run1[[5]])),]
startpoint3minmax
oldoptimforsmallexpimp = function(x,B0,sigmau,theta,points){
  return(1000000*expectedimprovement2Dminmax(x,B0,sigmau,theta,points))
}
oldoptimforsmallexpimp.2 = function(x,B0,sigmau,theta,points){
  return(10*expectedimprovement2Dminmax(x,B0,sigmau,theta,points))
}
possnewmax3minmax = optim(startpoint3minmax,oldoptimforsmallexpimp.2,B0=0,sigmau=0.6,theta=1.85,points=minmax3points,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmax3minmax
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmax3points,col="black",pch=19)
                 points(matrix(possnewmax3minmax$par,nrow=1,ncol=2),col="black",pch=4)})
newminmax3points = rbind(minmax3points,possnewmax3minmax$par)
newminmax3points
minmax3run2 = emumaxminfun(0,0.6,1.85,newminmax3points)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run2[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmax.2 = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmax3points))
sequem3accmatminmax.2 = matrix(sequem3accminmax.2,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax.2,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points,col="black",pch=19)})
startpoint3minmax.2 = newminmax3points[which(minmax3run2[[5]]==max(minmax3run2[[5]])),]
startpoint3minmax.2
possnewmax3minmax.2 = optim(startpoint3minmax.2,expectedimprovement2Dminmax,B0=0,sigmau=0.6,theta=1.85,points=newminmax3points,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmax3minmax.2
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax.2,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points,col="black",pch=19)
                 points(matrix(possnewmax3minmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
newminmax3points.2 = rbind(newminmax3points,possnewmax3minmax.2$par)
newminmax3points.2
minmax3run3 = emumaxminfun(0,0.6,1.85,newminmax3points.2)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmax.3 = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmax3points.2))
sequem3accmatminmax.3 = matrix(sequem3accminmax.3,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points.2,col="black",pch=19)})
startpoint3minmax.3 = newminmax3points.2[which(minmax3run3[[5]]==max(minmax3run3[[5]])),]
startpoint3minmax.3
possnewmax3minmax.3 = optim(startpoint3minmax.3,expectedimprovement2Dminmax,B0=0,sigmau=0.6,theta=1.85,points=newminmax3points.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmax3minmax.3
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points.2,col="black",pch=19)
                 points(matrix(possnewmax3minmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
newminmax3points.3 = rbind(newminmax3points.2,possnewmax3minmax.3$par)
newminmax3points.3
minmax3run4 = emumaxminfun(0,0.6,1.85,newminmax3points.3)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmax.4 = mapply(expectedimprovement2Dminmax,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmax3points.3))
sequem3accmatminmax.4 = matrix(sequem3accminmax.4,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmax.4,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3points.3,col="black",pch=19)})

###############################################################################################

filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run1[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmax3points,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmax3points,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmaxfourth = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=minmax3points))
sequem3accmatminmaxfourth = matrix(sequem3accminmaxfourth,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmax3points,col="black",pch=19)})
startpoint3minmaxfourth = minmax3points[which(minmax3run1[[5]]==max(minmax3run1[[5]])),]
startpoint3minmaxfourth
optimforsmallexpimp.3 = function(x,B0,sigmau,theta,points){
  return(100000*expectedimprovement2DMC4(x,B0,sigmau,theta,points))
}
possnewmax3minmaxfourth = optim(c(3,6),optimforsmallexpimp.3,B0=0,sigmau=0.6,theta=1.85,points=minmax3points,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmax3minmaxfourth
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(minmax3points,col="black",pch=19)
                 points(matrix(possnewmax3minmaxfourth$par,nrow=1,ncol=2),col="black",pch=4)})
newminmax3pointsfourth = rbind(minmax3points,possnewmax3minmaxfourth$par)
newminmax3pointsfourth
minmax3run2fourth = emumaxminfun(0,0.6,1.85,newminmax3pointsfourth)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run2fourth[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run2fourth[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmaxfourth.2 = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmax3pointsfourth))
sequem3accmatminmaxfourth.2 = matrix(sequem3accminmaxfourth.2,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth.2,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth,col="black",pch=19)})
startpoint3minmaxfourth.2 = newminmax3pointsfourth[which(minmax3run2fourth[[5]]==max(minmax3run2fourth[[5]])),]
startpoint3minmaxfourth.2
optimforsmallexpimp.2 = function(x,B0,sigmau,theta,points){
  return(1000000*expectedimprovement2DMC4(x,B0,sigmau,theta,points))
}
possnewmax3minmaxfourth.2 = optim(c(2,8),optimforsmallexpimp,B0=0,sigmau=0.6,theta=1.85,points=newminmax3pointsfourth,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmax3minmaxfourth.2
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth.2,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newminmax3pointsfourth,col="black",pch=19)
                 points(matrix(possnewmax3minmaxfourth.2$par,nrow=1,ncol=2),col="black",pch=4)})
newminmax3pointsfourth.2 = rbind(newminmax3pointsfourth,possnewmax3minmaxfourth.2$par)
newminmax3pointsfourth.2
minmax3run3fourth = emumaxminfun(0,0.6,1.85,newminmax3pointsfourth.2)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run3fourth[[2]],
               col = rainbow(25,alpha=0.95),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run3fourth[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmaxfourth.3 = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmax3pointsfourth.2))
sequem3accmatminmaxfourth.3 = matrix(sequem3accminmaxfourth.3,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth.2,col="black",pch=19)})
startpoint3minmaxfourth.3 = newminmax3pointsfourth.2[which(minmax3run3fourth[[5]]==max(minmax3run3fourth[[5]])),]
startpoint3minmaxfourth.3
possnewmax3minmaxfourth.3 = optim(startpoint3minmaxfourth.3,optimforsmallexpimp,B0=0,sigmau=0.6,theta=1.85,points=newminmax3pointsfourth.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((3*pi)+0.5,2))
possnewmax3minmaxfourth.3
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newminmax3pointsfourth.2,col="black",pch=19)
                 points(matrix(possnewmax3minmaxfourth.3$par,nrow=1,ncol=2),col="black",pch=4)})
newminmax3pointsfourth.3 = rbind(newminmax3pointsfourth.2,possnewmax3minmaxfourth.3$par)
newminmax3pointsfourth.3
minmax3run4fourth = emumaxminfun(0,0.6,1.85,newminmax3pointsfourth.3)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = minmax3run4fourth[[2]],
               col = rainbow(25,alpha=0.95),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = sqrt(minmax3run4fourth[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequem3accminmaxfourth.4 = mapply(expectedimprovement2DMC4,as.data.frame(t(sequemminmax)),MoreArgs = list(B0=0,sigmau=0.6,theta=1.85,points=newminmax3pointsfourth.3))
sequem3accmatminmaxfourth.4 = matrix(sequem3accminmaxfourth.4,40)
filled.contour(x = seq(-0.5,(3*pi)+0.5,length.out=40),
               y = seq(-0.5,(3*pi)+0.5,length.out=40), 
               z = -sequem3accmatminmaxfourth.4,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newminmax3pointsfourth.3,col="black",pch=19)})

###################################################################################################
