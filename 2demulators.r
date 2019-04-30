covfun2d = function(x,y,su,t){
  b1 = as.matrix(rbind(x,y))
  entr = su*exp(-(dist(b1,diag=FALSE)^2)/(t^2))
  return(entr)
}

latpoints = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
latpoints2 = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
latpoints3 = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))

emufuncosxplussiny = function(B0,sigmau,theta,points){
  a = points
  D = cos(a[,1])+sin(a[,2])
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

cosxplussinyrun1 = emufuncosxplussiny(0,0.6,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
        y = seq(-0.5,(2*pi)+0.5,length.out=40), 
        z = cosxplussinyrun1[[2]])
contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
        y = seq(-0.5,(2*pi)+0.5,length.out=40), 
        z = cosxplussinyrun1[[3]])
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = "x",ylab = "y"))
persp(seq(-0.5,(2*pi)+0.5,length.out=40), 
      seq(-0.5,(2*pi)+0.5,length.out=40), 
      cosxplussinyrun1[[3]], theta = 30, phi = 30, expand = 0.5,
      xlab = "x",ylab = "y",zlab = "z", col = "lightblue")


#######################################################################################################

impcosxplussiny=function(x,z,points){
  a = points
  D = cos(a[,1])+sin(a[,2])
  exf = 0
  varf = 0.6^2
  varep = 0
  vare = 0.0016
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (0.6^2)*exp(-(b^2)/(2^2))
  varDw = varD*(1-(10^(-6))) + (10^(-6))*diag(length(points[,1]))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g)*(1-(10^(-6))))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
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
  Isquared = ((eDf(x)-z)^2)/(varDf(x)+varep+vare)
  I = sqrt(Isquared)
  return(I)
}

sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                     y=seq(-0.5,(2*pi)+0.5,length.out=40))
sequemI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                 MoreArgs = list(z=1,points=latpoints))
sequemImat = matrix(sequemI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)})

persp(seq(-0.5,(2*pi)+0.5,length.out=40), 
      seq(-0.5,(2*pi)+0.5,length.out=40), 
      sequemImat, theta = 30, phi = 30, expand = 0.5, 
      col = "lightblue")

posnewpoints = mapply(impcosxplussiny,as.data.frame(t(latpoints2)),
                      MoreArgs = list(z=1,points=latpoints))
which(posnewpoints<3)
latpointsnew = as.matrix(rbind(latpoints,latpoints2[3,],
                               latpoints2[5,],latpoints2[7,],
                               latpoints2[9,]))


###############################################################################

cosxplussinyrun1new = emufuncosxplussiny(0,0.6,2,latpointsnew)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1new[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = "Emulator Expectation for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x", ylab = "y"))
contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
        y = seq(-0.5,(2*pi)+0.5,length.out=40), 
        z = cosxplussinyrun1new[[2]])
contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
        y = seq(-0.5,(2*pi)+0.5,length.out=40), 
        z = cosxplussinyrun1new[[3]])
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1new[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointsnew,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = "x",ylab = "y"))
persp(seq(-0.5,(2*pi)+0.5,length.out=40), 
      seq(-0.5,(2*pi)+0.5,length.out=40), 
      cosxplussinyrun1new[[2]], theta = 30, phi = 30, expand = 0.5,
      xlab = "x",ylab = "y",zlab = "z", col = "lightblue")

sequemInew = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                 MoreArgs = list(z=1,points=latpointsnew))
sequemImatnew = matrix(sequemInew,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImatnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImatnew)),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointsnew,col="black",pch=19)})

persp(seq(-0.5,(2*pi)+0.5,length.out=40), 
      seq(-0.5,(2*pi)+0.5,length.out=40), 
      sequemImatnew, theta = 30, phi = 30, expand = 0.5, 
      col = "lightblue")

posnewpoints2 = mapply(impcosxplussiny,as.data.frame(t(latpoints3)),
                      MoreArgs = list(z=1,points=latpointsnew))
which(posnewpoints2<3)
latpointsnew2 = as.matrix(rbind(latpointsnew,latpoints3[1,],
                               latpoints3[9,]))

#######################################################################

cosxplussinyrun1new2 = emufuncosxplussiny(0,0.6,2,latpointsnew2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1new2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = "Emulator Expectation for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x", ylab = "y"))
contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
        y = seq(-0.5,(2*pi)+0.5,length.out=40), 
        z = cosxplussinyrun1new2[[2]])
contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
        y = seq(-0.5,(2*pi)+0.5,length.out=40), 
        z = cosxplussinyrun1new2[[3]])
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1new2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointsnew2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = "x",ylab = "y"))
persp(seq(-0.5,(2*pi)+0.5,length.out=40), 
      seq(-0.5,(2*pi)+0.5,length.out=40), 
      cosxplussinyrun1new2[[2]], theta = 30, phi = 30, expand = 0.5,
      xlab = "x",ylab = "y",zlab = "z", col = "lightblue")

sequemInew2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=latpointsnew2))
sequemImatnew2 = matrix(sequemInew2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImatnew2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImatnew2)),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointsnew2,col="black",pch=19)})

persp(seq(-0.5,(2*pi)+0.5,length.out=40), 
      seq(-0.5,(2*pi)+0.5,length.out=40), 
      sequemImatnew2, theta = 30, phi = 30, expand = 0.5, 
      col = "lightblue")

##########################################################################
fdejong <- function (x, y) {
  return (cos(x)+sin(y))
}

x <- seq(0, 2*pi, length= 100)
y <- x
z <- outer(x, y, fdejong)
z[is.na(z)] <- 1

require(lattice)
wireframe(z, drape=T, col.regions=rainbow(100))
op <- par(bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

filled.contour(x,y,z,
               color.palette = colorRampPalette(c('green','yellow','red')))
