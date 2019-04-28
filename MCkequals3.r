expectedimprovement2DMC3test = function(x,B0,sigmau,theta,points){
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
    return((max(0,as.vector(mu)+(sqrt(as.vector(var))*e)-curmax))^3)
  }
  return(-sum(sapply(standnorm,I))/10000)
}

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
fivsequemaccMC3test = mapply(expectedimprovement2DMC3test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=ofivpoints))
fivsequemaccmatMC3test = matrix(fivsequemaccMC3test,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)})
fivstartpointMC3test = sequem[which(fivsequemaccMC3test==min(fivsequemaccMC3test)),]
fivstartpointMC3test
fivpossnewmaxMC3test = optim(fivstartpointMC3test,expectedimprovement2DMC3test,B0=0,sigmau=0.6,theta=2,points=ofivpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC3test
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC3test = rbind(ofivpoints,fivpossnewmaxMC3test$par)
newofivpointsMC3test
fivMC3testrun2 = emufuncosxplussiny(0,0.6,2,newofivpointsMC3test)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun2[[2]],
               col = rainbow(30,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC3testrun2[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC3test.2 = mapply(expectedimprovement2DMC3test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test))
fivsequemaccmatMC3test.2 = matrix(fivsequemaccMC3test.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.2,
               col = colors5,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test,col="black",pch=19)})
fivstartpointMC3test.2 = sequem[which(fivsequemaccMC3test.2==min(fivsequemaccMC3test.2)),]
fivstartpointMC3test.2
fivpossnewmaxMC3test.2 = optim(fivstartpointMC3test.2,expectedimprovement2DMC3test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC3test.2
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.2,
               col = colors5,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.2$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC3test.2 = rbind(newofivpointsMC3test,fivpossnewmaxMC3test.2$par)
newofivpointsMC3test.2
fivMC3testrun3 = emufuncosxplussiny(0,0.6,2,newofivpointsMC3test.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC3testrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.2,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC3test.3 = mapply(expectedimprovement2DMC3test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.2))
fivsequemaccmatMC3test.3 = matrix(fivsequemaccMC3test.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.2,col="black",pch=19)})
fivstartpointMC3test.3 = sequem[which(fivsequemaccMC3test.3==min(fivsequemaccMC3test.3)),]
fivstartpointMC3test.3
fivpossnewmaxMC3test.3 = optim(fivstartpointMC3test.3,expectedimprovement2DMC3test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC3test.3
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.3,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.2,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.3$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC3test.3 = rbind(newofivpointsMC3test.2,fivpossnewmaxMC3test.3$par)
newofivpointsMC3test.3
fivMC3testrun4 = emufuncosxplussiny(0,0.6,2,newofivpointsMC3test.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun4[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.3,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC3testrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.3,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC3test.4 = mapply(expectedimprovement2DMC3test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.3))
fivsequemaccmatMC3test.4 = matrix(fivsequemaccMC3test.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.3,col="black",pch=19)})
fivstartpointMC3test.4 = sequem[which(fivsequemaccMC3test.4==min(fivsequemaccMC3test.4)),]
fivstartpointMC3test.4
fivpossnewmaxMC3test.4 = optim(fivstartpointMC3test.4,expectedimprovement2DMC3test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC3test.4
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.4,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.3,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.4$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC3test.4 = rbind(newofivpointsMC3test.3,fivpossnewmaxMC3test.4$par)
newofivpointsMC3test.4
fivMC3testrun5 = emufuncosxplussiny(0,0.6,2,newofivpointsMC3test.4)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun5[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.4,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC3testrun5[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC3test.5 = mapply(expectedimprovement2DMC3test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.4))
fivsequemaccmatMC3test.5 = matrix(fivsequemaccMC3test.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.4,col="black",pch=19)})
fivstartpointMC3test.5 = sequem[which(fivsequemaccMC3test.5==min(fivsequemaccMC3test.5)),]
fivstartpointMC3test.5
fivpossnewmaxMC3test.5 = optim(fivstartpointMC3test.5,expectedimprovement2DMC3test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC3test.5
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.5,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.4,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.5$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC3test.5 = rbind(newofivpointsMC3test.4,fivpossnewmaxMC3test.5$par)
newofivpointsMC3test.5
fivMC3testrun6 = emufuncosxplussiny(0,0.6,2,newofivpointsMC3test.5)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun6[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivMC3testrun6[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.5,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
fivsequemaccMC3test.6 = mapply(expectedimprovement2DMC3test,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.5))
fivsequemaccmatMC3test.6 = matrix(fivsequemaccMC3test.6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==2.5,", ",theta==1.875)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.5,col="black",pch=19)})
fivstartpointMC3test.6 = sequem[which(fivsequemaccMC3test.6==min(fivsequemaccMC3test.6)),]
fivstartpointMC3test.6
fivpossnewmaxMC3test.6 = optim(fivstartpointMC3test.6,expectedimprovement2DMC3test,B0=0,sigmau=0.6,theta=2,points=newofivpointsMC3test.5,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep((2*pi)+0.5,2))
fivpossnewmaxMC3test.6
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.6,
               col = topo.colors(30,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.5,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.6$par,nrow=1,ncol=2),col="black",pch=4)})
newofivpointsMC3test.6 = rbind(newofivpointsMC3test.5,fivpossnewmaxMC3test.6$par)
newofivpointsMC3test.6
fivMC3testrun7 = emufuncosxplussiny(0,0.6,2,newofivpointsMC3test.6)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun7[[2]],
               col = rainbow(25,alpha=0.95),
               levels = seq(-2,2,length.out = 21),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==1.85)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
