standnorm = rnorm(10000)

expectedimprovement2DMC = function(x,B0,sigmau,theta,points){
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
  I = function(e){
    return(max(0,as.vector(mu)+(sqrt(as.vector(var))*e)-curmax))
  }
  standnormnew = as.vector(mu) + sqrt(as.vector(var))*standnorm
  return(sum(sapply(standnorm,I))/10000)
}

sequemaccMC = mapply(expectedimprovement2DMC,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=oavgpoints))
sequemaccmatMC = matrix(sequemaccMC,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemaccmatMC,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)})


##############################################################################################################

expectedimprovement2DMC4test = function(x,B0,sigmau,theta,points){
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
  I = function(e){
    return((max(0,as.vector(mu)+(sqrt(as.vector(var))*e)-curmax))^4)
  }
  return(-sum(sapply(standnorm,I))/10000)
}

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccMC4test = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=oavgpoints))
sequemaccmatMC4test = matrix(sequemaccMC4test,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(oavgpoints,col="black",pch=19)})
startpointMC4test = oavgpoints[which(avoptvarrun1[[5]]==max(avoptvarrun1[[5]])),]
startpointMC4test
optimforsmallexpimpMC4test = function(x,B0,sigmau,theta,points){
  return(1000000*expectedimprovement2DMC4test(x,B0,sigmau,theta,points))
}
possnewmaxMC4test = optim(startpointMC4test,optimforsmallexpimpMC4test,B0=0,sigmau=0.6,theta=2,points=oavgpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxMC4test
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(possnewmaxMC4test$par,nrow=1,ncol=2),col="black",pch=4)})
newoavgpointsMC4test = rbind(oavgpoints,possnewmaxMC4test$par)
newoavgpointsMC4test
MC4testrun2 = emufuncosxplussiny(0,0.6,1.85,newoavgpointsMC4test)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = MC4testrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newoavgpointsMC4test,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(MC4testrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newoavgpointsMC4test,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccMC4test.2 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newoavgpointsMC4test))
sequemaccmatMC4test.2 = matrix(sequemaccMC4test.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newoavgpointsMC4test,col="black",pch=19)})
startpointMC4test.2 = newoavgpointsMC4test[which(MC4testrun2[[5]]==max(MC4testrun2[[5]])),]
startpointMC4test.2
optimforsmallexpimpMC4test.2 = function(x,B0,sigmau,theta,points){
  return(100000*expectedimprovement2DMC4test(x,B0,sigmau,theta,points))
}
possnewmaxMC4test.2 = optim(c(-0.5,0),optimforsmallexpimpMC4test.2,B0=0,sigmau=0.6,theta=2,points=newoavgpointsMC4test,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxMC4test.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newoavgpointsMC4test,col="black",pch=19)
                 points(matrix(possnewmaxMC4test.2$par,nrow=1,ncol=2),col="black",pch=4)})
newoavgpointsMC4test.2 = rbind(newoavgpointsMC4test,possnewmaxMC4test.2$par)
newoavgpointsMC4test.2
MC4testrun3 = emufuncosxplussiny(0,0.6,1.85,newoavgpointsMC4test.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = MC4testrun3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newoavgpointsMC4test.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(MC4testrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newoavgpointsMC4test.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccMC4test.3 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newoavgpointsMC4test.2))
sequemaccmatMC4test.3 = matrix(sequemaccMC4test.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newoavgpointsMC4test.2,col="black",pch=19)})
startpointMC4test.3 = newoavgpointsMC4test.2[which(MC4testrun3[[5]]==max(MC4testrun3[[5]])),]
startpointMC4test.3
optimforsmallexpimpMC4test.2 = function(x,B0,sigmau,theta,points){
  return(100000*expectedimprovement2DMC4test(x,B0,sigmau,theta,points))
}
possnewmaxMC4test.3 = optim(c(6.5,1),optimforsmallexpimpMC4test.2,B0=0,sigmau=0.6,theta=2,points=newoavgpointsMC4test.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxMC4test.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newoavgpointsMC4test.2,col="black",pch=19)
                 points(matrix(possnewmaxMC4test.3$par,nrow=1,ncol=2),col="black",pch=4)})
newoavgpointsMC4test.3 = rbind(newoavgpointsMC4test.2,possnewmaxMC4test.3$par)
newoavgpointsMC4test.3
MC4testrun4 = emufuncosxplussiny(0,0.6,1.85,newoavgpointsMC4test.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = MC4testrun4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newoavgpointsMC4test.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(MC4testrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newoavgpointsMC4test.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemaccMC4test.4 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newoavgpointsMC4test.3))
sequemaccmatMC4test.4 = matrix(sequemaccMC4test.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatMC4test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newoavgpointsMC4test.3,col="black",pch=19)})

########################################################################################

fivoptvarrun1 = emufuncosxplussiny(0,0.6,2,ofivpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=ofivpoints))
fivsequemaccmatEItest = matrix(fivsequemaccEItest,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)})
fivstartpointEItest = sequem[which(fivsequemaccEItest==min(fivsequemaccEItest)),]
fivstartpointEItest
fivpossnewmaxEItest = optim(fivstartpointEItest,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=ofivpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest = rbind(ofivpoints,fivpossnewmaxEItest$par)
newofivpointsEItest
fivEIoptvarrun2 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun2[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.1 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest))
fivsequemaccmatEItest.1 = matrix(fivsequemaccEItest.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.1,
               col = colors,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest,col="black",pch=19)})
fivstartpointEItest.1 = sequem[which(fivsequemaccEItest.1==min(fivsequemaccEItest.1)),]
fivstartpointEItest.1
fivpossnewmaxEItest.1 = optim(fivstartpointEItest.1,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.1$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.1 = rbind(newofivpointsEItest,fivpossnewmaxEItest.1$par)
newofivpointsEItest.1
fivEIoptvarrun3 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun3[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.2 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.1))
fivsequemaccmatEItest.2 = matrix(fivsequemaccEItest.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.1,col="black",pch=19)})
fivstartpointEItest.2 = sequem[which(fivsequemaccEItest.2==min(fivsequemaccEItest.2)),]
fivstartpointEItest.2
fivpossnewmaxEItest.2 = optim(fivstartpointEItest.2,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.1,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.2$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.2 = rbind(newofivpointsEItest.1,fivpossnewmaxEItest.2$par)
newofivpointsEItest.2
fivEIoptvarrun4 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun4[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.3 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.2))
fivsequemaccmatEItest.3 = matrix(fivsequemaccEItest.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.2,col="black",pch=19)})
fivstartpointEItest.3 = sequem[which(fivsequemaccEItest.3==min(fivsequemaccEItest.3)),]
fivstartpointEItest.3
fivpossnewmaxEItest.3 = optim(fivstartpointEItest.3,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.2,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.3$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.3 = rbind(newofivpointsEItest.2,fivpossnewmaxEItest.3$par)
newofivpointsEItest.3
fivEIoptvarrun5 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun5[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun5[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.4 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.3))
fivsequemaccmatEItest.4 = matrix(fivsequemaccEItest.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.3,col="black",pch=19)})
fivstartpointEItest.4 = sequem[which(fivsequemaccEItest.4==min(fivsequemaccEItest.4)),]
fivstartpointEItest.4
fivpossnewmaxEItest.4 = optim(fivstartpointEItest.4,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.3,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.4$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.4 = rbind(newofivpointsEItest.3,fivpossnewmaxEItest.4$par)
newofivpointsEItest.4
fivEIoptvarrun6 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun6[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun6[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.5 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.4))
fivsequemaccmatEItest.5 = matrix(fivsequemaccEItest.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.4,col="black",pch=19)})
fivstartpointEItest.5 = sequem[which(fivsequemaccEItest.5==min(fivsequemaccEItest.5)),]
fivstartpointEItest.5
fivpossnewmaxEItest.5 = optim(fivstartpointEItest.5,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.4,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.5$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.5 = rbind(newofivpointsEItest.4,fivpossnewmaxEItest.5$par)
newofivpointsEItest.5
fivEIoptvarrun7 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun7[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun7[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.6 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.5))
fivsequemaccmatEItest.6 = matrix(fivsequemaccEItest.6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.5,col="black",pch=19)})
fivstartpointEItest.6 = sequem[which(fivsequemaccEItest.6==min(fivsequemaccEItest.6)),]
fivstartpointEItest.6
fivpossnewmaxEItest.6 = optim(fivstartpointEItest.6,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.5,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.6
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.5,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.6$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.6 = rbind(newofivpointsEItest.5,fivpossnewmaxEItest.6$par)
newofivpointsEItest.6
fivEIoptvarrun8 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.6)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun8[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivEIoptvarrun8[[3]]),
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointsEItest.6,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccEItest.7 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.6))
fivsequemaccmatEItest.7 = matrix(fivsequemaccEItest.7,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.6,col="black",pch=19)})
fivstartpointEItest.7 = sequem[which(fivsequemaccEItest.7==min(fivsequemaccEItest.7)),]
fivstartpointEItest.7
fivpossnewmaxEItest.7 = optim(fivstartpointEItest.7,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newofivpointsEItest.6,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxEItest.7
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.6,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.7$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsEItest.7 = rbind(newofivpointsEItest.6,fivpossnewmaxEItest.7$par)
newofivpointsEItest.7

#############################################################################################

fivoptvarrun1 = emufuncosxplussiny(0,0.6,2,ofivpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=ofivpoints))
fivsequemaccmatMC4test = matrix(fivsequemaccMC4test,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)})
fivstartpointMC4test = sequem[which(fivsequemaccMC4test==min(fivsequemaccMC4test)),]
fivstartpointMC4test
fivpossnewmaxMC4test = optim(fivstartpointMC4test,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=ofivpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test = rbind(ofivpoints,fivpossnewmaxMC4test$par)
newofivpointsMC4test
fivMC4testrun2 = emufuncosxplussiny(0,0.6,2,newofivpointsMC4test)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC4testrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC4testrun2[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test.2 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test))
fivsequemaccmatMC4test.2 = matrix(fivsequemaccMC4test.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test,col="black",pch=19)})
fivstartpointMC4test.2 = sequem[which(fivsequemaccMC4test.2==min(fivsequemaccMC4test.2)),]
fivstartpointMC4test.2
fivpossnewmaxMC4test.2 = optim(fivstartpointMC4test.2,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test.2$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test.2 = rbind(newofivpointsMC4test,fivpossnewmaxMC4test.2$par)
newofivpointsMC4test.2
fivMC4testrun3 = emufuncosxplussiny(0,0.6,2,newofivpointsMC4test.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC4testrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC4testrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test.3 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.2))
fivsequemaccmatMC4test.3 = matrix(fivsequemaccMC4test.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.2,col="black",pch=19)})
fivstartpointMC4test.3 = sequem[which(fivsequemaccMC4test.3==min(fivsequemaccMC4test.3)),]
fivstartpointMC4test.3
fivpossnewmaxMC4test.3 = optim(fivstartpointMC4test.3,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.2,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test.3$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test.3 = rbind(newofivpointsMC4test.2,fivpossnewmaxMC4test.3$par)
newofivpointsMC4test.3
fivMC4testrun4 = emufuncosxplussiny(0,0.6,2,newofivpointsMC4test.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC4testrun4[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC4testrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test.4 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.3))
fivsequemaccmatMC4test.4 = matrix(fivsequemaccMC4test.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.3,col="black",pch=19)})
fivstartpointMC4test.4 = sequem[which(fivsequemaccMC4test.4==min(fivsequemaccMC4test.4)),]
fivstartpointMC4test.4
fivpossnewmaxMC4test.4 = optim(fivstartpointMC4test.4,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.3,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test.4$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test.4 = rbind(newofivpointsMC4test.3,fivpossnewmaxMC4test.4$par)
newofivpointsMC4test.4
fivMC4testrun5 = emufuncosxplussiny(0,0.6,2,newofivpointsMC4test.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC4testrun5[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC4testrun5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test.5 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.4))
fivsequemaccmatMC4test.5 = matrix(fivsequemaccMC4test.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.4,col="black",pch=19)})
fivstartpointMC4test.5 = sequem[which(fivsequemaccMC4test.5==min(fivsequemaccMC4test.5)),]
fivstartpointMC4test.5
fivpossnewmaxMC4test.5 = optim(fivstartpointMC4test.5,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.4,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test.5$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test.5 = rbind(newofivpointsMC4test.4,fivpossnewmaxMC4test.5$par)
newofivpointsMC4test.5
fivMC4testrun6 = emufuncosxplussiny(0,0.6,2,newofivpointsMC4test.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC4testrun6[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC4testrun6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test.6 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.5))
fivsequemaccmatMC4test.6 = matrix(fivsequemaccMC4test.6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.5,col="black",pch=19)})
fivstartpointMC4test.6 = sequem[which(fivsequemaccMC4test.6==min(fivsequemaccMC4test.6)),]
fivstartpointMC4test.6
fivpossnewmaxMC4test.6 = optim(fivstartpointMC4test.6,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.5,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test.6
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.5,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test.6$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test.6 = rbind(newofivpointsMC4test.5,fivpossnewmaxMC4test.6$par)
newofivpointsMC4test.6
fivMC4testrun7 = emufuncosxplussiny(0,0.6,2,newofivpointsMC4test.6)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC4testrun7[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC4testrun7[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC4test.6,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC4test.7 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.6))
fivsequemaccmatMC4test.7 = matrix(fivsequemaccMC4test.7,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.6,col="black",pch=19)})
fivstartpointMC4test.7 = sequem[which(fivsequemaccMC4test.7==min(fivsequemaccMC4test.7)),]
fivstartpointMC4test.7
fivpossnewmaxMC4test.7 = optim(fivstartpointMC4test.7,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC4test.6,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC4test.7
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC4test.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC4test.6,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC4test.7$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC4test.7 = rbind(newofivpointsMC4test.6,fivpossnewmaxMC4test.7$par)
newofivpointsMC4test.7
cos(newofivpointsMC4test.7[12,1])+sin(newofivpointsMC4test.7[12,2])

#################################################################################################


fivlowthetaoptvarrun1 = emufuncosxplussiny(0,0.6,1,ofivpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=ofivpoints))
sequemacclowthetamat = matrix(sequemacclowtheta,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)})
startpointlowtheta = sequem[which(sequemacclowtheta==min(sequemacclowtheta)),]
startpointlowtheta
possnewmaxlowtheta = optim(startpointlowtheta,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=ofivpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta = rbind(ofivpoints,possnewmaxlowtheta$par)
newfivpointslowtheta
fivlowthetaoptvarrun2 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.1 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta))
sequemacclowthetamat.1 = matrix(sequemacclowtheta.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta,col="black",pch=19)})
startpointlowtheta.1 = sequem[which(sequemacclowtheta.1==min(sequemacclowtheta.1)),]
startpointlowtheta.1
possnewmaxlowtheta.1 = optim(startpointlowtheta.1,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.1$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.1 = rbind(newfivpointslowtheta,possnewmaxlowtheta.1$par)
newfivpointslowtheta.1
fivlowthetaoptvarrun3 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.2 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.1))
sequemacclowthetamat.2 = matrix(sequemacclowtheta.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.1,col="black",pch=19)})
startpointlowtheta.2 = sequem[which(sequemacclowtheta.2==min(sequemacclowtheta.2)),]
startpointlowtheta.2
possnewmaxlowtheta.2 = optim(startpointlowtheta.2,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.1,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.2$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.2 = rbind(newfivpointslowtheta.1,possnewmaxlowtheta.2$par)
newfivpointslowtheta.2
fivlowthetaoptvarrun4 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.3 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.2))
sequemacclowthetamat.3 = matrix(sequemacclowtheta.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.2,col="black",pch=19)})
startpointlowtheta.3 = sequem[which(sequemacclowtheta.3==min(sequemacclowtheta.3)),]
startpointlowtheta.3
possnewmaxlowtheta.3 = optim(startpointlowtheta.3,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.2,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.3$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.3 = rbind(newfivpointslowtheta.2,possnewmaxlowtheta.3$par)
newfivpointslowtheta.3
fivlowthetaoptvarrun5 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun5[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.4 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.3))
sequemacclowthetamat.4 = matrix(sequemacclowtheta.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.3,col="black",pch=19)})
startpointlowtheta.4 = sequem[which(sequemacclowtheta.4==min(sequemacclowtheta.4)),]
startpointlowtheta.4
possnewmaxlowtheta.4 = optim(startpointlowtheta.4,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.3,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.4$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.4 = rbind(newfivpointslowtheta.3,possnewmaxlowtheta.4$par)
newfivpointslowtheta.4
fivlowthetaoptvarrun6 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun6[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.5 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.4))
sequemacclowthetamat.5 = matrix(sequemacclowtheta.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.4,col="black",pch=19)})
startpointlowtheta.5 = sequem[which(sequemacclowtheta.5==min(sequemacclowtheta.5)),]
startpointlowtheta.5
possnewmaxlowtheta.5 = optim(startpointlowtheta.5,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.4,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.5$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.5 = rbind(newfivpointslowtheta.4,possnewmaxlowtheta.5$par)
newfivpointslowtheta.5
fivlowthetaoptvarrun7 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun7[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun7[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.6 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.5))
sequemacclowthetamat.6 = matrix(sequemacclowtheta.6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.5,col="black",pch=19)})
startpointlowtheta.6 = sequem[which(sequemacclowtheta.6==min(sequemacclowtheta.6)),]
startpointlowtheta.6
possnewmaxlowtheta.6 = optim(startpointlowtheta.6,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.5,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.6
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.5,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.6$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.6 = rbind(newfivpointslowtheta.5,possnewmaxlowtheta.6$par)
newfivpointslowtheta.6
fivlowthetaoptvarrun6 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun6[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.5 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.4))
sequemacclowthetamat.5 = matrix(sequemacclowtheta.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.4,col="black",pch=19)})
startpointlowtheta.5 = sequem[which(sequemacclowtheta.5==min(sequemacclowtheta.5)),]
startpointlowtheta.5
possnewmaxlowtheta.5 = optim(startpointlowtheta.5,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.4,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.5$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.5 = rbind(newfivpointslowtheta.4,possnewmaxlowtheta.5$par)
newfivpointslowtheta.5
fivlowthetaoptvarrun8 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.6)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun8[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun8[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.6,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.7 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.6))
sequemacclowthetamat.7 = matrix(sequemacclowtheta.7,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.6,col="black",pch=19)})
startpointlowtheta.7 = sequem[which(sequemacclowtheta.7==min(sequemacclowtheta.7)),]
startpointlowtheta.7
possnewmaxlowtheta.7 = optim(startpointlowtheta.7,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.6,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.7
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.6,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.7$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.7 = rbind(newfivpointslowtheta.6,possnewmaxlowtheta.7$par)
newfivpointslowtheta.7
fivlowthetaoptvarrun9 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun9[[2]],
               col = rainbow(20,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.7,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun9[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.7,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.8 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.7))
sequemacclowthetamat.8 = matrix(sequemacclowtheta.8,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.8,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.7,col="black",pch=19)})
startpointlowtheta.8 = sequem[which(sequemacclowtheta.8==min(sequemacclowtheta.8)),]
startpointlowtheta.8
possnewmaxlowtheta.8 = optim(startpointlowtheta.8,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.7,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.8
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.8,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.7,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.8$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.8 = rbind(newfivpointslowtheta.7,possnewmaxlowtheta.8$par)
newfivpointslowtheta.8
fivlowthetaoptvarrun10 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.8)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun10[[2]],
               col = rainbow(20,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.8,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun10[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.8,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.9 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.8))
sequemacclowthetamat.9 = matrix(sequemacclowtheta.9,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.9,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.8,col="black",pch=19)})
startpointlowtheta.9 = sequem[which(sequemacclowtheta.9==min(sequemacclowtheta.9)),]
startpointlowtheta.9
possnewmaxlowtheta.9 = optim(startpointlowtheta.9,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.8,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.9
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.9,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.8,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.9$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.9 = rbind(newfivpointslowtheta.8,possnewmaxlowtheta.9$par)
newfivpointslowtheta.9
fivlowthetaoptvarrun11 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.9)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun11[[2]],
               col = rainbow(20,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.9,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun11[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.9,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.10 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.9))
sequemacclowthetamat.10 = matrix(sequemacclowtheta.10,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.10,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.9,col="black",pch=19)})
startpointlowtheta.10 = sequem[which(sequemacclowtheta.10==min(sequemacclowtheta.10)),]
startpointlowtheta.10
possnewmaxlowtheta.10 = optim(startpointlowtheta.10,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.9,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.10
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.10,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.9,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.10$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.10 = rbind(newfivpointslowtheta.9,possnewmaxlowtheta.10$par)
newfivpointslowtheta.10
fivlowthetaoptvarrun12 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.10)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun12[[2]],
               col = rainbow(20,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.10,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun12[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.10,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.11 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.10))
sequemacclowthetamat.11 = matrix(sequemacclowtheta.11,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.11,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.10,col="black",pch=19)})
startpointlowtheta.11 = sequem[which(sequemacclowtheta.11==min(sequemacclowtheta.11)),]
startpointlowtheta.11
possnewmaxlowtheta.11 = optim(startpointlowtheta.11,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.10,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.11
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.11,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.10,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.11$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.11 = rbind(newfivpointslowtheta.10,possnewmaxlowtheta.11$par)
newfivpointslowtheta.11
fivlowthetaoptvarrun13 = emufuncosxplussiny(0,0.6,1,newfivpointslowtheta.11)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivlowthetaoptvarrun13[[2]],
               col = rainbow(20,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newfivpointslowtheta.11,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivlowthetaoptvarrun13[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newfivpointslowtheta.11,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
sequemacclowtheta.12 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.11))
sequemacclowthetamat.12 = matrix(sequemacclowtheta.12,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.12,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.11,col="black",pch=19)})
startpointlowtheta.12 = sequem[which(sequemacclowtheta.12==min(sequemacclowtheta.12)),]
startpointlowtheta.12
possnewmaxlowtheta.12 = optim(startpointlowtheta.12,expectedimprovement2D,B0=0,sigmau=0.6,theta=1,points=newfivpointslowtheta.11,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
possnewmaxlowtheta.12
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemacclowthetamat.12,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newfivpointslowtheta.11,col="black",pch=19)
                 points(matrix(possnewmaxlowtheta.12$par,nrow=1,ncol=2),col="black",pch=4)})
newfivpointslowtheta.12 = rbind(newfivpointslowtheta.11,possnewmaxlowtheta.12$par)
newfivpointslowtheta.12

###############################################################################################################

fivoptvarlowthetaMC4run1 = emufuncosxplussiny(0,0.6,1,ofivpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run1[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=ofivpoints))
fivsequemaccmatlowthetaMC4test = matrix(fivsequemacclowthetaMC4test,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)})
fivstartpointlowthetaMC4test = sequem[which(fivsequemacclowthetaMC4test==min(fivsequemacclowthetaMC4test)),]
fivstartpointlowthetaMC4test
fivpossnewmaxlowthetaMC4test = optim(fivstartpointlowthetaMC4test,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=ofivpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test = rbind(ofivpoints,fivpossnewmaxlowthetaMC4test$par)
newofivpointslowthetaMC4test
fivoptvarlowthetaMC4run2 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run2[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.1 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test))
fivsequemaccmatlowthetaMC4test.1 = matrix(fivsequemacclowthetaMC4test.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test,col="black",pch=19)})
fivstartpointlowthetaMC4test.1 = sequem[which(fivsequemacclowthetaMC4test.1==min(fivsequemacclowthetaMC4test.1)),]
fivstartpointlowthetaMC4test.1
fivpossnewmaxlowthetaMC4test.1 = optim(fivstartpointlowthetaMC4test.1,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.1
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.1,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.1$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.1 = rbind(newofivpointslowthetaMC4test,fivpossnewmaxlowthetaMC4test.1$par)
newofivpointslowthetaMC4test.1
fivoptvarlowthetaMC4run3 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run3[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.1,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run3[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.1,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.2 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.1))
fivsequemaccmatlowthetaMC4test.2 = matrix(fivsequemacclowthetaMC4test.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.1,col="black",pch=19)})
fivstartpointlowthetaMC4test.2 = sequem[which(fivsequemacclowthetaMC4test.2==min(fivsequemacclowthetaMC4test.2)),]
fivstartpointlowthetaMC4test.2
fivpossnewmaxlowthetaMC4test.2 = optim(fivstartpointlowthetaMC4test.2,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.1,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.2,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.1,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.2$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.2 = rbind(newofivpointslowthetaMC4test.1,fivpossnewmaxlowthetaMC4test.2$par)
newofivpointslowthetaMC4test.2
fivoptvarlowthetaMC4run4 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run4[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run4[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.3 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.2))
fivsequemaccmatlowthetaMC4test.3 = matrix(fivsequemacclowthetaMC4test.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.2,col="black",pch=19)})
fivstartpointlowthetaMC4test.3 = sequem[which(fivsequemacclowthetaMC4test.3==min(fivsequemacclowthetaMC4test.3)),]
fivstartpointlowthetaMC4test.3
fivpossnewmaxlowthetaMC4test.3 = optim(fivstartpointlowthetaMC4test.3,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.2,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.3$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.3 = rbind(newofivpointslowthetaMC4test.2,fivpossnewmaxlowthetaMC4test.3$par)
newofivpointslowthetaMC4test.3
fivoptvarlowthetaMC4run5 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run5[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run5[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.4 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.3))
fivsequemaccmatlowthetaMC4test.4 = matrix(fivsequemacclowthetaMC4test.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.3,col="black",pch=19)})
fivstartpointlowthetaMC4test.4 = sequem[which(fivsequemacclowthetaMC4test.4==min(fivsequemacclowthetaMC4test.4)),]
fivstartpointlowthetaMC4test.4
fivpossnewmaxlowthetaMC4test.4 = optim(fivstartpointlowthetaMC4test.4,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.3,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.4$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.4 = rbind(newofivpointslowthetaMC4test.3,fivpossnewmaxlowthetaMC4test.4$par)
newofivpointslowthetaMC4test.4
fivoptvarlowthetaMC4run6 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run6[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run6[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.5 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.4))
fivsequemaccmatlowthetaMC4test.5 = matrix(fivsequemacclowthetaMC4test.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.4,col="black",pch=19)})
fivstartpointlowthetaMC4test.5 = sequem[which(fivsequemacclowthetaMC4test.5==min(fivsequemacclowthetaMC4test.5)),]
fivstartpointlowthetaMC4test.5
fivpossnewmaxlowthetaMC4test.5 = optim(fivstartpointlowthetaMC4test.5,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.4,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.5$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.5 = rbind(newofivpointslowthetaMC4test.4,fivpossnewmaxlowthetaMC4test.5$par)
newofivpointslowthetaMC4test.5
fivoptvarlowthetaMC4run7 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run7[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run7[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.6 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.5))
fivsequemaccmatlowthetaMC4test.6 = matrix(fivsequemacclowthetaMC4test.6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.5,col="black",pch=19)})
fivstartpointlowthetaMC4test.6 = sequem[which(fivsequemacclowthetaMC4test.6==min(fivsequemacclowthetaMC4test.6)),]
fivstartpointlowthetaMC4test.6
fivpossnewmaxlowthetaMC4test.6 = optim(fivstartpointlowthetaMC4test.6,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.5,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.6
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.5,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.6$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.6 = rbind(newofivpointslowthetaMC4test.5,fivpossnewmaxlowthetaMC4test.6$par)
newofivpointslowthetaMC4test.6
fivoptvarlowthetaMC4run8 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.6)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run8[[2]],
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run8[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.6,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.7 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.6))
fivsequemaccmatlowthetaMC4test.7 = matrix(fivsequemacclowthetaMC4test.7,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.6,col="black",pch=19)})
fivstartpointlowthetaMC4test.7 = sequem[which(fivsequemacclowthetaMC4test.7==min(fivsequemacclowthetaMC4test.7)),]
fivstartpointlowthetaMC4test.7
fivpossnewmaxlowthetaMC4test.7 = optim(fivstartpointlowthetaMC4test.7,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.6,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.7
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.7,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.6,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.7$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.7 = rbind(newofivpointslowthetaMC4test.6,fivpossnewmaxlowthetaMC4test.7$par)
newofivpointslowthetaMC4test.7
fivoptvarlowthetaMC4run9 = emufuncosxplussiny(0,0.6,1,newofivpointslowthetaMC4test.7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarlowthetaMC4run9[[2]],
               col = rainbow(20,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointslowthetaMC4test.7,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarlowthetaMC4run9[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(newofivpointslowthetaMC4test.7,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemacclowthetaMC4test.8 = mapply(expectedimprovement2DMC4test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.7))
fivsequemaccmatlowthetaMC4test.8 = matrix(fivsequemacclowthetaMC4test.8,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.8,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.7,col="black",pch=19)})
fivstartpointlowthetaMC4test.8 = sequem[which(fivsequemacclowthetaMC4test.8==min(fivsequemacclowthetaMC4test.8)),]
fivstartpointlowthetaMC4test.8
fivpossnewmaxlowthetaMC4test.8 = optim(fivstartpointlowthetaMC4test.8,expectedimprovement2DMC4test,B0=0,sigmau=0.6,theta=1,points=newofivpointslowthetaMC4test.7,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxlowthetaMC4test.8
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatlowthetaMC4test.8,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointslowthetaMC4test.7,col="black",pch=19)
                 points(matrix(fivpossnewmaxlowthetaMC4test.8$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointslowthetaMC4test.8 = rbind(newofivpointslowthetaMC4test.7,fivpossnewmaxlowthetaMC4test.8$par)
newofivpointslowthetaMC4test.8
