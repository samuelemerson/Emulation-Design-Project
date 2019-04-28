robotangles <- function(x)
{
  ##########################################################################
  #
  # ROBOT ARM FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUTS:
  #
  # y = distance from the end of the arm to the origin
  # xx = c(theta1, theta2, theta3, theta4, L1, L2, L3, L4)
  #
  #########################################################################
  
  theta <- c(pi/2,2*pi*x[1],2*pi*x[2],pi/2)
  L     <- c(1,1,1,1)
  
  thetamat <- matrix(rep(theta,times=4), 4, 4, byrow=TRUE)
  thetamatlow <- thetamat
  thetamatlow[upper.tri(thetamatlow)] <- 0
  sumtheta <- rowSums(thetamatlow)
  
  u <- sum(L*cos(sumtheta))
  v <- sum(L*sin(sumtheta))
  
  y <- (u^2 + v^2)^(0.5)
  return(y)
}

robotforlooks = function(x){
  theta <- c(pi/4,0,0,2*pi*x[1])
  L     <- c(1,1,1,x[2])
  
  thetamat <- matrix(rep(theta,times=4), 4, 4, byrow=TRUE)
  thetamatlow <- thetamat
  thetamatlow[upper.tri(thetamatlow)] <- 0
  sumtheta <- rowSums(thetamatlow)
  
  u <- sum(L*cos(sumtheta))
  v <- sum(L*sin(sumtheta))
  
  y <- (u^2 + v^2)^(0.5)
  return(y)
}

robotgrid = expand.grid(x=seq(0,1,length.out=100),
                        y=seq(0,1,length.out=100))
robotvaluesgrid = mapply(robotforlooks,as.data.frame(t(robotgrid)))
zrob = matrix(robotvaluesgrid,100)
require(lattice)
wireframe(zrob, drape=T, col.regions=rainbow(100))
op <- par(bg = "white")
persp(seq(0,1,length.out=100), seq(0,1,length.out=100), zrob, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

filled.contour(x = seq(0,1,length.out=100),
               y = seq(0,1,length.out=100), 
               z = zrob,
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Robot Arm Function ")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
#robot(c(pi/2,1))

robotarmemulator = function(B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(robotforlooks,as.data.frame(t(a))))
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
  sequem = expand.grid(x=seq(0,1,length.out=40),
                       y=seq(0,1,length.out=40))
  sequem2 = mapply(eDf,as.data.frame(t(sequem)))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  sequem2mat = matrix(sequem2,40)
  sequem3mat = matrix(sequem3,40)
  return(list(sequem,sequem2mat,sequem3mat,a,D))
}

avgvarformanyrobotarm = function(x,numpoints){
  a = matrix(x,nrow = numpoints,ncol = 2)
  varf = 0.5^2
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
  sequem = expand.grid(x=seq(0,1,length.out=40),
                       y=seq(0,1,length.out=40))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  average = mean(sequem3)
  return(average)
}

robot7startpoints = (osevpoints)/(2*pi)
robot7startpoints
robot7optimpoints.0 = optim(as.vector(robot7startpoints),avgvarformanyrobotarm,numpoints=7,method="L-BFGS-B",lower=rep(0,18),upper=rep(1,18))
robot7points.0 = matrix(robot7optimpoints.0$par,nrow=7,ncol=2)
robot7points.0

robotarmrun1 = robotarmemulator(3,0.5,0.33,robot7points.0)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.0,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot7points.0,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
EIrobotarm = function(x,B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(robotforlooks,as.data.frame(t(a))))
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

robotsequem = expand.grid(x=seq(0,1,length.out=40),
                          y=seq(0,1,length.out=40))
robotsequemacc = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.0))
robotsequemaccmat = matrix(robotsequemacc,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.0,col="black",pch=19)})
RAstartpoint = robotsequem[which(robotsequemaccmat==min(robotsequemaccmat)),]
RAstartpoint
RApossnewmax = optim(RAstartpoint,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot7points.0,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RApossnewmax
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.0,col="black",pch=19)
                 points(matrix(RApossnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
robot7points.1 = rbind(robot7points.0,RApossnewmax$par)
robot7points.1
robotarmrun2 = robotarmemulator(3,0.5,0.33,robot7points.1)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot7points.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robotsequemacc.1 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.1))
robotsequemaccmat.1 = matrix(robotsequemacc.1,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.1,col="black",pch=19)})
RAstartpoint.1 = robotsequem[which(robotsequemaccmat.1==min(robotsequemaccmat.1)),]
RAstartpoint.1
RApossnewmax.1 = optim(RAstartpoint.1,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot7points.1,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RApossnewmax.1
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.1,col="black",pch=19)
                 points(matrix(RApossnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
robot7points.2 = rbind(robot7points.1,RApossnewmax.1$par)
robot7points.2
robotarmrun3 = robotarmemulator(3,0.5,0.33,robot7points.2)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot7points.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robotsequemacc.2 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.2))
robotsequemaccmat.2 = matrix(robotsequemacc.2,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.2,col="black",pch=19)})
RAstartpoint.2 = robotsequem[which(robotsequemaccmat.2==min(robotsequemaccmat.2)),]
RAstartpoint.2
RApossnewmax.2 = optim(RAstartpoint.2,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot7points.2,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RApossnewmax.2
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.2,col="black",pch=19)
                 points(matrix(RApossnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
robot7points.3 = rbind(robot7points.2,RApossnewmax.2$par)
robot7points.3
robotarmrun4 = robotarmemulator(3,0.5,0.33,robot7points.3)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot7points.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robotsequemacc.3 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.3))
robotsequemaccmat.3 = matrix(robotsequemacc.3,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.3,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.3,col="black",pch=19)})
RAstartpoint.3 = robotsequem[which(robotsequemaccmat.3==min(robotsequemaccmat.3)),]
RAstartpoint.3
RApossnewmax.3 = optim(RAstartpoint.3,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot7points.3,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RApossnewmax.3
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.3,
               col = topo.colors(20),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.3,col="black",pch=19)
                 points(matrix(RApossnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
robot7points.4 = rbind(robot7points.3,RApossnewmax.3$par)
robot7points.4
robotarmrun5 = robotarmemulator(3,0.5,0.33,robot7points.4)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun5[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot7points.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robotsequemacc.4 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.4))
robotsequemaccmat.4 = matrix(robotsequemacc.4,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.4,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.4,col="black",pch=19)})
RAstartpoint.4 = robotsequem[which(robotsequemaccmat.4==min(robotsequemaccmat.4)),]
RAstartpoint.4
RApossnewmax.4 = optim(RAstartpoint.4,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot7points.4,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RApossnewmax.4
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.4,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.4,col="black",pch=19)
                 points(matrix(RApossnewmax.4$par,nrow=1,ncol=2),col="black",pch=4)})
robot7points.5 = rbind(robot7points.4,RApossnewmax.4$par)
robot7points.5
robotarmrun6 = robotarmemulator(3,0.5,0.33,robot7points.5)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun6[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot7points.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robotsequemacc.5 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.5))
robotsequemaccmat.5 = matrix(robotsequemacc.5,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccmat.5,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.5,col="black",pch=19)})

#######################################################################################

EIMC4robot = function(x,B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(robotforlooks,as.data.frame(t(a))))
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

filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarmrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.0,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarmrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.0,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robotsequemaccMC = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot7points.0))
robotsequemaccMCmat = matrix(robotsequemaccMC,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot7points.0,col="black",pch=19)})
MCRAstartpoint = robot7points.0[which(robotarmrun1[[5]]==max(robotarmrun1[[5]])),]
MCRAstartpoint
RAoptimforsmallexpimp = function(x,B0,sigmau,theta,points){
  return(10000*EIMC4robot(x,B0,sigmau,theta,points))
}
MCRApossnewmax = optim(c(1,0.7),EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=robot7points.0,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRApossnewmax
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(robot7points.0,col="black",pch=19)
                 points(matrix(MCRApossnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot7points.1 = rbind(robot7points.0,MCRApossnewmax$par)
MCrobot7points.1
MCrobotarmrun2 = robotarmemulator(3,0.5,0.33,MCrobot7points.1)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarmrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarmrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robotsequemaccMC.1 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.1))
robotsequemaccMCmat.1 = matrix(robotsequemaccMC.1,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.1,col="black",pch=19)})
MCRAstartpoint.1 = robotsequem[which(robotsequemaccMC.1==min(robotsequemaccMC.1)),]
MCRAstartpoint.1
RAoptimforsmallexpimp = function(x,B0,sigmau,theta,points){
  return(10000*EIMC4robot(x,B0,sigmau,theta,points))
}
MCRApossnewmax.1 = optim(MCRAstartpoint.1,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.1,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRApossnewmax.1
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.1,col="black",pch=19)
                 points(matrix(MCRApossnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot7points.2 = rbind(MCrobot7points.1,MCRApossnewmax.1$par)
MCrobot7points.2
MCrobotarmrun3 = robotarmemulator(3,0.5,0.33,MCrobot7points.2)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarmrun3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarmrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robotsequemaccMC.2 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.2))
robotsequemaccMCmat.2 = matrix(robotsequemaccMC.2,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.2,col="black",pch=19)})
MCRAstartpoint.2 = robotsequem[which(robotsequemaccMC.2==min(robotsequemaccMC.2)),]
MCRAstartpoint.2
MCRApossnewmax.2 = optim(MCRAstartpoint.2,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.2,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRApossnewmax.2
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.2,col="black",pch=19)
                 points(matrix(MCRApossnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot7points.3 = rbind(MCrobot7points.2,MCRApossnewmax.2$par)
MCrobot7points.3
MCrobotarmrun4 = robotarmemulator(3,0.5,0.33,MCrobot7points.3)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarmrun4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarmrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robotsequemaccMC.3 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.3))
robotsequemaccMCmat.3 = matrix(robotsequemaccMC.3,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.3,col="black",pch=19)})
MCRAstartpoint.3 = robotsequem[which(robotsequemaccMC.3==min(robotsequemaccMC.3)),]
MCRAstartpoint.3
MCRApossnewmax.3 = optim(MCRAstartpoint.3,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.3,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRApossnewmax.3
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.3,col="black",pch=19)
                 points(matrix(MCRApossnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot7points.4 = rbind(MCrobot7points.3,MCRApossnewmax.3$par)
MCrobot7points.4
MCrobotarmrun5 = robotarmemulator(3,0.5,0.33,MCrobot7points.4)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarmrun5[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarmrun5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robotsequemaccMC.4 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot7points.4))
robotsequemaccMCmat.4 = matrix(robotsequemaccMC.4,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robotsequemaccMCmat.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot7points.4,col="black",pch=19)})

############################################################################################################

robot3startpoints = (othrpoints)/(2*pi)
robot3startpoints
robot3optimpoints.0 = optim(as.vector(robot3startpoints),avgvarformanyrobotarm,numpoints=3,method="L-BFGS-B",lower=rep(0,6),upper=rep(1,6))
robot3points.0 = matrix(robot3optimpoints.0$par,nrow=3,ncol=2)
robot3points.0

robotarm3run1 = robotarmemulator(3,0.5,0.33,robot3points.0)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(robot3points.0,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.0))
robot3sequemaccmat = matrix(robot3sequemacc,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)})
RA3startpoint = robotsequem[which(robot3sequemacc==min(robot3sequemacc)),]
RA3startpoint
RA3possnewmax = optim(RA3startpoint,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.0,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)
                 points(matrix(RA3possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.1 = rbind(robot3points.0,RA3possnewmax$par)
robot3points.1
robotarm3run2 = robotarmemulator(3,0.5,0.33,robot3points.1)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.1 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.1))
robot3sequemaccmat.1 = matrix(robot3sequemacc.1,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.1,col="black",pch=19)})
RA3startpoint.1 = robotsequem[which(robot3sequemacc.1==min(robot3sequemacc.1)),]
RA3startpoint.1
RA3possnewmax.1 = optim(RA3startpoint.1,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.1,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.1
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.1,col="black",pch=19)
                 points(matrix(RA3possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.2 = rbind(robot3points.1,RA3possnewmax.1$par)
robot3points.2
robotarm3run3 = robotarmemulator(3,0.5,0.33,robot3points.2)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.2 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.2))
robot3sequemaccmat.2 = matrix(robot3sequemacc.2,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.2,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.2,col="black",pch=19)})
RA3startpoint.2 = robotsequem[which(robot3sequemacc.2==min(robot3sequemacc.2)),]
RA3startpoint.2
RA3possnewmax.2 = optim(RA3startpoint.2,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.2,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.2
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.2,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.2,col="black",pch=19)
                 points(matrix(RA3possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.3 = rbind(robot3points.2,RA3possnewmax.2$par)
robot3points.3
robotarm3run4 = robotarmemulator(3,0.5,0.33,robot3points.3)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.3 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.3))
robot3sequemaccmat.3 = matrix(robot3sequemacc.3,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.3,col="black",pch=19)})
RA3startpoint.3 = robotsequem[which(robot3sequemacc.3==min(robot3sequemacc.3)),]
RA3startpoint.3
RA3possnewmax.3 = optim(RA3startpoint.3,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.3,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.3
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.3,col="black",pch=19)
                 points(matrix(RA3possnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.4 = rbind(robot3points.3,RA3possnewmax.3$par)
robot3points.4
robotarm3run5 = robotarmemulator(3,0.5,0.33,robot3points.4)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run5[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.4 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.4))
robot3sequemaccmat.4 = matrix(robot3sequemacc.4,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.4,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.4,col="black",pch=19)})
RA3startpoint.4 = robotsequem[which(robot3sequemacc.4==min(robot3sequemacc.4)),]
RA3startpoint.4
RA3possnewmax.4 = optim(RA3startpoint.4,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.4,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.4
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.4,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.4,col="black",pch=19)
                 points(matrix(RA3possnewmax.4$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.5 = rbind(robot3points.4,RA3possnewmax.4$par)
robot3points.5
robotarm3run6 = robotarmemulator(3,0.5,0.33,robot3points.5)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run6[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.5 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.5))
robot3sequemaccmat.5 = matrix(robot3sequemacc.5,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.5,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.5,col="black",pch=19)})
RA3startpoint.5 = robotsequem[which(robot3sequemacc.5==min(robot3sequemacc.5)),]
RA3startpoint.5
RA3possnewmax.5 = optim(RA3startpoint.5,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.5,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.5
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.5,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.5,col="black",pch=19)
                 points(matrix(RA3possnewmax.5$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.6 = rbind(robot3points.5,RA3possnewmax.5$par)
robot3points.6
robotarm3run7 = robotarmemulator(3,0.5,0.33,robot3points.6)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run7[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run7[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.6,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.6 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.6))
robot3sequemaccmat.6 = matrix(robot3sequemacc.6,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.6,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.6,col="black",pch=19)})
RA3startpoint.6 = robotsequem[which(robot3sequemacc.6==min(robot3sequemacc.6)),]
RA3startpoint.6
RA3possnewmax.6 = optim(RA3startpoint.6,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.6,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.6
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.6,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.6,col="black",pch=19)
                 points(matrix(RA3possnewmax.6$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.7 = rbind(robot3points.6,RA3possnewmax.6$par)
robot3points.7
robotarm3run8 = robotarmemulator(3,0.5,0.33,robot3points.7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run8[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.7,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run8[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.7,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.7 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.7))
robot3sequemaccmat.7 = matrix(robot3sequemacc.7,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.7,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.7,col="black",pch=19)})
RA3startpoint.7 = robotsequem[which(robot3sequemacc.7==min(robot3sequemacc.7)),]
RA3startpoint.7
RA3possnewmax.7 = optim(RA3startpoint.7,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.7,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.7
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.7,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.7,col="black",pch=19)
                 points(matrix(RA3possnewmax.7$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.8 = rbind(robot3points.7,RA3possnewmax.7$par)
robot3points.8
robotarm3run9 = robotarmemulator(3,0.5,0.33,robot3points.8)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run9[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.8,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run9[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.8,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.8 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.8))
robot3sequemaccmat.8 = matrix(robot3sequemacc.8,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.8,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.8,col="black",pch=19)})
RA3startpoint.8 = robotsequem[which(robot3sequemacc.8==min(robot3sequemacc.8)),]
RA3startpoint.8
RA3possnewmax.8 = optim(RA3startpoint.8,EIrobotarm,B0=3,sigmau=0.5,theta=0.33,points=robot3points.8,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
RA3possnewmax.8
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.8,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.8,col="black",pch=19)
                 points(matrix(RA3possnewmax.8$par,nrow=1,ncol=2),col="black",pch=4)})
robot3points.9 = rbind(robot3points.8,RA3possnewmax.8$par)
robot3points.9
robotarm3run10 = robotarmemulator(3,0.5,0.33,robot3points.9)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run10[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.9,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run10[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.9,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))

robot3sequemacc.9 = mapply(EIrobotarm,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.9))
robot3sequemaccmat.9 = matrix(robot3sequemacc.9,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.9,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==3,", ",sigma[u]==0.5,", ",theta==0.33)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.9,col="black",pch=19)})
######################################################################################################

filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(robotarm3run1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=robot3points.0))
robot3sequemaccMCmat = matrix(robot3sequemaccMC,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)})
MCRA3startpoint = robotsequem[which(robot3sequemaccMC==min(robot3sequemaccMC)),]
MCRA3startpoint
MCRA3possnewmax = optim(MCRA3startpoint,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=robot3points.0,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(robot3points.0,col="black",pch=19)
                 points(matrix(MCRA3possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.1 = rbind(robot3points.0,MCRA3possnewmax$par)
MCrobot3points.1
MCrobotarm3run2 = robotarmemulator(3,0.5,0.33,MCrobot3points.1)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.1 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.1))
robot3sequemaccMCmat.1 = matrix(robot3sequemaccMC.1,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.1,col="black",pch=19)})
MCRA3startpoint.1 = robotsequem[which(robot3sequemaccMC.1==min(robot3sequemaccMC.1)),]
MCRA3startpoint.1
MCRA3possnewmax.1 = optim(MCRA3startpoint.1,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.1,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax.1
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(MCrobot3points.1,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.2 = rbind(MCrobot3points.1,MCRA3possnewmax.1$par)
MCrobot3points.2
MCrobotarm3run3 = robotarmemulator(3,0.5,0.33,MCrobot3points.2)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.2 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.2))
robot3sequemaccMCmat.2 = matrix(robot3sequemaccMC.2,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.2,col="black",pch=19)})
MCRA3startpoint.2 = robotsequem[which(robot3sequemaccMC.2==min(robot3sequemaccMC.2)),]
MCRA3startpoint.2
MCRA3possnewmax.2 = optim(MCRA3startpoint.2,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.2,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax.2
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.2,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(MCrobot3points.2,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.3 = rbind(MCrobot3points.2,MCRA3possnewmax.2$par)
MCrobot3points.3
MCrobotarm3run4 = robotarmemulator(3,0.5,0.33,MCrobot3points.3)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.3 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.3))
robot3sequemaccMCmat.3 = matrix(robot3sequemaccMC.3,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.3,col="black",pch=19)})
MCRA3startpoint.3 = robotsequem[which(robot3sequemaccMC.3==min(robot3sequemaccMC.3)),]
MCRA3startpoint.3
MCRA3possnewmax.3 = optim(MCRA3startpoint.3,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.3,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax.3
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.3,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(MCrobot3points.3,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.4 = rbind(MCrobot3points.3,MCRA3possnewmax.3$par)
MCrobot3points.4
MCrobotarm3run5 = robotarmemulator(3,0.5,0.33,MCrobot3points.4)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run5[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.4 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.4))
robot3sequemaccMCmat.4 = matrix(robot3sequemaccMC.4,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.4,col="black",pch=19)})
MCRA3startpoint.4 = robotsequem[which(robot3sequemaccMC.4==min(robot3sequemaccMC.4)),]
MCRA3startpoint.4
MCRA3possnewmax.4 = optim(MCRA3startpoint.4,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.4,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax.4
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.4,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(MCrobot3points.4,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.4$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.5 = rbind(MCrobot3points.4,MCRA3possnewmax.4$par)
MCrobot3points.5
MCrobotarm3run6 = robotarmemulator(3,0.5,0.33,MCrobot3points.5)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run6[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.5 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.5))
robot3sequemaccMCmat.5 = matrix(robot3sequemaccMC.5,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.5,col="black",pch=19)})
MCRA3startpoint.5 = robotsequem[which(robot3sequemaccMC.5==min(robot3sequemaccMC.5)),]
MCRA3startpoint.5
MCRA3possnewmax.5 = optim(MCRA3startpoint.5,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.5,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax.5
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.5,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(MCrobot3points.5,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.5$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.6 = rbind(MCrobot3points.5,MCRA3possnewmax.5$par)
MCrobot3points.6
MCrobotarm3run7 = robotarmemulator(3,0.5,0.33,MCrobot3points.6)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run7[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run7[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.6,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.6 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.6))
robot3sequemaccMCmat.6 = matrix(robot3sequemaccMC.6,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.6,col="black",pch=19)})
MCRA3startpoint.6 = robotsequem[which(robot3sequemaccMC.6==min(robot3sequemaccMC.6)),]
MCRA3startpoint.6
MCRA3possnewmax.6 = optim(MCRA3startpoint.6,EIMC4robot,B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.6,method="L-BFGS-B",lower=rep(0,2),upper=rep(1,2))
MCRA3possnewmax.6
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.6,
               col = topo.colors(28,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.6,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.6$par,nrow=1,ncol=2),col="black",pch=4)})
MCrobot3points.7 = rbind(MCrobot3points.6,MCRA3possnewmax.6$par)
MCrobot3points.7
MCrobotarm3run8 = robotarmemulator(3,0.5,0.33,MCrobot3points.7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = MCrobotarm3run8[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.7,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = sqrt(MCrobotarm3run8[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.7,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
robot3sequemaccMC.7 = mapply(EIMC4robot,as.data.frame(t(robotsequem)),MoreArgs = list(B0=3,sigmau=0.5,theta=0.33,points=MCrobot3points.7))
robot3sequemaccMCmat.7 = matrix(robot3sequemaccMC.7,40)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.7,col="black",pch=19)})
