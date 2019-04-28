emufunexpsq = function(B0,sigmau,theta){
  a = seq(0.1,0.5,by=0.08)
  D = (-3.5+exp(3.5*a))^2
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length(a))
  b = as.matrix(dist(a))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  covfD = function(x){
    g = covfun(x,a,varf,theta)
    return(g)
  }
  covDf = function(x){
    i = covfun(a,x,varf,theta)
    return(i)
  }
  eDf = function(x){
    h = exf + covfD(x)%*%solve(varD)%*%(D-eD)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%solve(varD)%*%covDf(x)
    return(j)
  }
  sequem = seq(0.05,0.55,length.out=10000)
  sequem2 = c()
  for (i in sequem){
    sequem2 = append(sequem2,eDf(i))
  }
  sequem3 = c()
  for (i in sequem){
    sequem3 = append(sequem3,varDf(i))
  }
  upvar = sequem2 + 3*sqrt(sequem3)
  downvar = sequem2 - 3*sqrt(sequem3)
  return(list(sequem,sequem2,a,D,upvar,downvar))
}

pdf("NewEmulatormodel.pdf")
emufunexpsqrun1 = emufunexpsq(10,2.5,0.1)
plot(emufunexpsqrun1[[1]],emufunexpsqrun1[[2]],type='l',lwd=2,col='blue',
     ylim=c(-1,11.5),xlab="x",
     ylab="f(x)")
lines(emufunexpsqrun1[[1]],emufunexpsqrun1[[5]],type='l',lwd=2,col='red')
lines(emufunexpsqrun1[[1]],emufunexpsqrun1[[6]],type='l',lwd=2,col='red')
points(emufunexpsqrun1[[3]],emufunexpsqrun1[[4]],col='black',cex=1.2,pch=19)
lines(sequemsinexp,realepsq,type='l',lwd=2,col='black')
legend(x=0.075,y=12,legend=c(expression(E[D](f(x))),expression(f(x)),expression(paste("±",3*sqrt("Var"[D](f(x))))),"D"),col=c("blue","black","red","black"),lty=c(1,1,1,NA),pch=c(NA,NA,NA,19),bty="n",lwd=2)
dev.off()
sequemsinexp = seq(0.05,0.55,length.out=10000)
#plot(sequemsinexp,(-3.5+exp(3.5*sequemsinexp))^2,type='l')
#min((-3.5+exp(3.5*sequemsinexp))^2)
realepsq = (-3.5+exp(3.5*sequemsinexp))^2

pdf("emexptit.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main=expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("emvartit.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1[[3]]),
               col = rainbow(13,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main=expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("emimpnotit.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)})
dev.off()
length(which(sequemImat<3))/1600

pdf("emminvarnotit.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(xlab = expression(x[1]),ylab = expression(x[2])))
dev.off()
pdf("emminvarimpnotit.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(avoptvarImat)),
               plot.title = title(xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)})
dev.off()
length(which(avoptvarImat<3))/1600
pdf("emeximpnotit.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

##########################################################################################

pdf("emexptitgoodone.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(main="Emulator Expectation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

pdf("emvartitgoodone.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(main="Emulator Standard Deviation",
                                  xlab = expression(x[1]),ylab = expression(x[2])))
dev.off()

pdf("emEIwave1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("emEIwave2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints,col="black",pch=19)
                 points(matrix(possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("emEIwave3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.2,
               col = topo.colors(26,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.1,col="black",pch=19)
                 points(matrix(possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

#sequemacc.3 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpoints.2))
#sequemaccmat.3 = matrix(sequemacc.3,40)
#filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
y = seq(-0.5,(2*pi)+0.5,length.out=40), 
z = -sequemaccmat.3,
col = topo.colors(26,alpha=1.0),
plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                   xlab = expression(x[1]), ylab = expression(x[2])),
plot.axes = {axis(1,seq(from=0,to=6,by=1))
  axis(2,seq(from=0,to=6,by=1))
  points(newmaxpoints.2,col="black",pch=19)})

#####################################################################################################################

latpointssev = as.matrix(2*pi*maximinLHS(7,2,optimize.on = "grid"))
presbaddesrun1 = emufuncosxplussiny(0,0.6,2,latpointssev)
pdf("emsevbadexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = presbaddesrun1[[2]],
               col = rainbow(24,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointssev,col="black",pch=19)},
               plot.title = title(main = "Emulator Expectation",
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("emsevbadvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(presbaddesrun1[[3]]),
               col = rainbow(17,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointssev,col="black",pch=19)},
               plot.title = title(main = "Emulator Standard Deviation",
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
sevoptvarrun1 = emufuncosxplussiny(0,0.6,2,osevpoints)
pdf("emsevgoodexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarrun1[[2]],
               col = rainbow(35,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Expectation",
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
pdf("emsevgoodvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun1[[3]]),
               col = rainbow(40,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
pdf("emrealfun.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemfmat,
               col = rainbow(25,alpha=0.95),
               plot.title = title(main = expression(paste(f(x[1],x[2])," = ",cos(x[1])+sin(x[2]),sep='')),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

sequemaccsev = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=osevpoints))
sequemaccmatsev = matrix(sequemaccsev,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)})
startpointsev = osevpoints[which(sevoptvarrun1[[5]]==max(sevoptvarrun1[[5]])),]
startpointsev
possnewmaxsev = optim(startpointsev,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=osevpoints,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmaxsev
pdf("emEIwave1sev.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)
                 points(matrix(possnewmaxsev$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpointssev = rbind(osevpoints,possnewmaxsev$par)
newmaxpointssev
sevoptvarrun1acc = emufuncosxplussiny(0,0.6,2,newmaxpointssev)
sequemaccsev.1 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpointssev))
sequemaccmatsev.1 = matrix(sequemaccsev.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev,col="black",pch=19)})
startpointsev.1 = newmaxpointssev[which(sevoptvarrun1acc[[5]]==max(sevoptvarrun1acc[[5]])),]
startpointsev.1
possnewmaxsev.1 = optim(startpointsev.1,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpointssev,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmaxsev.1
pdf("emEIwave2sev.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev.1,
               col = topo.colors(25,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev,col="black",pch=19)
                 points(matrix(possnewmaxsev.1$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpointssev.2 = rbind(newmaxpointssev,possnewmaxsev.1$par)
newmaxpointssev.2
sevoptvarrun1acc.2 = emufuncosxplussiny(0,0.6,2,newmaxpointssev.2)
sequemaccsev.2 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpointssev.2))
sequemaccmatsev.2 = matrix(sequemaccsev.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev.2,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev.2,col="black",pch=19)})
startpointsev.2 = newmaxpointssev.2[which(sevoptvarrun1acc.2[[5]]==max(sevoptvarrun1acc.2[[5]])),]
startpointsev.2
possnewmaxsev.2 = optim(startpointsev.2,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpointssev.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmaxsev.2
pdf("emEIwave3sev.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev.2,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev.2,col="black",pch=19)
                 points(matrix(possnewmaxsev.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpointssev.3 = rbind(newmaxpointssev.2,possnewmaxsev.2$par)
newmaxpointssev.3
sevoptvarrun1acc.3 = emufuncosxplussiny(0,0.6,2,newmaxpointssev.3)
sequemaccsev.3 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpointssev.3))
sequemaccmatsev.3 = matrix(sequemaccsev.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev.3,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev.3,col="black",pch=19)})
startpointsev.3 = newmaxpointssev.3[which(sevoptvarrun1acc.3[[5]]==max(sevoptvarrun1acc.3[[5]])),]
startpointsev.3
possnewmaxsev.3 = optim(startpointsev.3,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpointssev.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmaxsev.3
pdf("emEIwave4sev.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmatsev.3,
               col = topo.colors(27,alpha=1.0),
               plot.title = title(main = "Expected Improvement",
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev.3,col="black",pch=19)
                 points(matrix(possnewmaxsev.3$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpointssev.4 = rbind(newmaxpointssev.3,possnewmaxsev.3$par)
newmaxpointssev.4
sevoptvarrun1acc.4 = emufuncosxplussiny(0,0.6,2,newmaxpointssev.4)
pdf("emsevgoodexpimp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarrun1acc.4[[2]],
               col = rainbow(20,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpointssev.4,col="black",pch=19)},
               plot.title = title(main = "Emulator Expectation",
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

###############################################################################################

realexpthreefive = exp(3.5*sequemsinexp)
pdf("OriginalEmulator.pdf")
plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="x",ylab="f(x)",main=expression("Exponential Model Emulator"))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)
dev.off()

pdf("EmulatorbetaOG.pdf")

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)

dev.off()

pdf("Emulatorbetalow.pdf")
plot(run3.0[[1]],run3.0[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==-3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run3.0[[1]],run3.0[[5]],type='l',lwd=1,col='red')
lines(run3.0[[1]],run3.0[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run3.0[[3]],run3.0[[4]],col='purple',pch=19)

dev.off()

pdf("Emulatorbetahigh.pdf")
plot(run5.0[[1]],run5.0[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==10.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run5.0[[1]],run5.0[[5]],type='l',lwd=1,col='red')
lines(run5.0[[1]],run5.0[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run5.0[[3]],run5.0[[4]],col='purple',pch=19)

dev.off()

pdf("Emulatorvarysigmau.pdf")
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run3.1 = emufunexp(3.5,0.1,0.14)
pdf("Emulatorlowsigma.pdf")
plot(run3.1[[1]],run3.1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==0.1,", ",theta==0.14)))
lines(run3.1[[1]],run3.1[[5]],type='l',lwd=1,col='red')
lines(run3.1[[1]],run3.1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run3.1[[3]],run3.1[[4]],col='purple',pch=19)
dev.off()

run5.1 = emufunexp(3.5,2.9,0.14)
pdf("Emulatorhighsigma.pdf")
plot(run5.1[[1]],run5.1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==2.9,", ",theta==0.14)))
lines(run5.1[[1]],run5.1[[5]],type='l',lwd=1,col='red')
lines(run5.1[[1]],run5.1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run5.1[[3]],run5.1[[4]],col='purple',pch=19)

dev.off()

pdf("Emulatorvarytheta.pdf")
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run3.2 = emufunexp(3.5,1.5,0.01)
pdf("Emulatorlowtheta.pdf")
plot(run3.2[[1]],run3.2[[2]],type='l',lwd=1,col='blue',
     ylim=c(-1.5,8.5),xlab="x",
     ylab="f(x)",main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.01)))
lines(run3.2[[1]],run3.2[[5]],type='l',lwd=1,col='red')
lines(run3.2[[1]],run3.2[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run3.2[[3]],run3.2[[4]],col='purple',pch=19)
dev.off()

run5.2 = emufunexp(3.5,1.5,0.27)
pdf("Emulatorhightheta.pdf")
plot(run5.2[[1]],run5.2[[2]],type='l',lwd=1,col='blue',
     ylim=c(1,6),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.27)))
lines(run5.2[[1]],run5.2[[5]],type='l',lwd=1,col='red')
lines(run5.2[[1]],run5.2[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run5.2[[3]],run5.2[[4]],col='purple',pch=19)

dev.off()

pdf("Emulatorinteractions.pdf")
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)

pdf("Emulatorinter1.pdf")
plot(run2.6[[1]],run2.6[[2]],type='l',lwd=1,col='blue',
     ylim=c(-4,6),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==-3.5,", ",sigma[u]==0.1,", ",theta==0.01)))
lines(run2.6[[1]],run2.6[[5]],type='l',lwd=1,col='red')
lines(run2.6[[1]],run2.6[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run2.6[[3]],run2.6[[4]],col='purple',pch=19)
dev.off()

pdf("Emulatorinter2.pdf")
plot(run3.5[[1]],run3.5[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==2.9,", ",theta==0.27)))
lines(run3.5[[1]],run3.5[[5]],type='l',lwd=1,col='red')
lines(run3.5[[1]],run3.5[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run3.5[[3]],run3.5[[4]],col='purple',pch=19)

dev.off()
#####################################################################################
pdf("2Doriginal.pdf",width = 7.6,height = 7)
cosxplussinyrun1 = emufuncosxplussiny(0,0.6,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("2Dvarybeta0.pdf",width = 7.6,height = 7)
cosxplussinyruntest = emufuncosxplussiny(3.5,0.6,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyruntest[[2]],
               col = rainbow(27,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==3.5,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

pdf("2Doriginalvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1[[3]]),
               col = rainbow(13,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main=expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
cosxplussinyrun4 = emufuncosxplussiny(0,0.4,2,latpoints)
pdf("2Dvarysigmau.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun4[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.4,", ",theta==2)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

cosxplussinyrun2 = emufuncosxplussiny(0,0.6,0.4,latpoints)
pdf("2Dvarythetaexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun2[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==0.4)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("2Dvarythetavar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun2[[3]]),
               col = rainbow(13,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==0.4)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
cosxplussinyrun3 = emufuncosxplussiny(0,0.8,4,latpoints)
pdf("2Dvarythetaexp2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun3[[2]],
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.8,", ",theta==4)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("2Dvarythetavar2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun3[[3]]),
               col = rainbow(35,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.8,", ",theta==4)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()
pdf("Cxoriginal.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = abs(testem),
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,7.5),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(abs(C(x))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

pdf("Cxvarysigma.pdf",width = 7.6,height = 7)
testemvarysig = Cx(0,0.4,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = abs(testemvarysig),
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,7.5),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(abs(C(x))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
testemvarythe = Cx(0,0.6,4,latpoints)
pdf("Cxvarytheta.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = abs(testemvarythe),
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,7.5,max(testemvarythe)),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(abs(C(x))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

###############################################################################################

pdf("EmulatorbetaOGz.pdf")

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator")))
axis(2,at=3.5)
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
abline(h=3.5,col='black',lwd=1)
abline(h=3.5+sqrt(0.025) ,col='black',lty=2)
abline(h=3.5-sqrt(0.025),col='black',lty=2)
points(run1[[3]],run1[[4]],col='purple',pch=19)

dev.off()

sequemdatapoints1D = c()
for (i in run1[[3]]){
  sequemdatapoints1D = append(sequemdatapoints1D,impexp(i,3.5))
}
newdatpointsimp1D = c(impexp(0.333,3.5),impexp(0.367,3.5))

pdf("1DEmulatorimp.pdf")
plot(sequem,sequem2,type='l',lwd=1,col='darkorange',xlab="x",ylab="I(x)",xlim=c(0.068,0.533),
     main=expression(paste(I(x)," for the Exponential Model Emulator")))
axis(2,at=c(3))
abline(h=3,col='black')
points(run1[[3]],sequemdatapoints1D,col='purple',pch=19)
locator()
dev.off()
###################################################################################################

pdf("emrealfunimp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemrealmat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemrealmat)),
               plot.title = title(main=expression(paste(tilde(I)(x)," for ",f(x[1],x[2])," = ",cos(x[1])+sin(x[2]),sep='')),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

pdf("impinitialem.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)})
dev.off()
############################################################################################

pdf("1DEmulatorimpnewpoints.pdf")
plot(sequem,sequem2,type='l',lwd=1,col='darkorange',xlab="x",ylab="I(x)",xlim=c(0.068,0.533),
     main=expression(paste(I(x)," for the Exponential Model Emulator")))
axis(2,at=c(3))
abline(h=3,col='black')
points(run1[[3]],sequemdatapoints1D,col='purple',pch=19)
points(c(0.333,0.367),newdatpointsimp1D,col='purple',pch=4,lwd=2)
locator()
dev.off()

pdf("1DEmulatornewpointsthin.pdf")
newrunex1D = emufunexp2(3.5,1.5,0.14,c(0.1,0.2,0.3,0.333,0.367,0.4,0.5))
plot(newrunex1D[[1]],newrunex1D[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator")))
axis(2,at=3.5)
lines(newrunex1D[[1]],newrunex1D[[5]],type='l',lwd=1,col='red')
lines(newrunex1D[[1]],newrunex1D[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
abline(h=3.5,col='black',lwd=1)
abline(h=3.5+sqrt(0.025) ,col='black',lty=2)
abline(h=3.5-sqrt(0.025),col='black',lty=2)
points(newrunex1D[[3]],newrunex1D[[4]],col='purple',pch=19)
dev.off()

newrunex1Dimp = overallsequem(c(0.1,0.2,0.3,0.333,0.367,0.4,0.5))
sequemdatapoints1D.1 = c()
for (i in newrunex1D[[3]]){
  sequemdatapoints1D.1 = append(sequemdatapoints1D.1,impexp2(i,3.5,c(0.1,0.2,0.3,0.333,0.367,0.4,0.5)))
}
sequemdatapoints1D.1

pdf("1DEmulatorimpnewpointsadded.pdf")
plot(newrunex1Dimp[[1]],newrunex1Dimp[[2]],type='l',lwd=1,col='darkorange',xlab="x",ylab="I(x)",xlim=c(0.068,0.533),
     main=expression(paste(I(x)," for the Exponential Model Emulator")))
axis(2,at=c(3))
abline(h=3,col='black')
points(newrunex1D[[3]],sequemdatapoints1D.1,col='purple',pch=19)
locator()
dev.off()
#############################################################################################

pdf("impinitialemnewpoints.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)
                 points(lpnew[10:13,],col="black",pch=4)})
dev.off()

length(which(sequemImat>=3))/1600

pdf("impwave2emnewpoints.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew,col="black",pch=19)
                 points(lpnew2[14:16,],col="black",pch=4)})
dev.off()

pdf("impwave3emnewpoints.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImnew2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew2,col="black",pch=19)
                 points(lpnew3[17:18,],col="black",pch=4)})
dev.off()

length(which(sImnew2>=3))/1600

##############################################################################################

pdf("impinitialemnewpointsbyeye.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)
                 points(byhp,col="black",pch=4)})
dev.off()
pdf("impwave2emnewpointsbyeye.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImbh, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad,col="black",pch=19)
                 points(byhp2,col="black",pch=4)})
dev.off()

pdf("impwave3emnewpointsbyeye.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImbh2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad2,col="black",pch=19)
                 points(matrix(c(0,6.25,3.5,-0.25),nrow=2,ncol=2),col="black",pch=4)})
dev.off()
############################################################################################

pdf("gridexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridrun[[2]],
               col = rainbow(50,alpha=0.95),
               levels = c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,
                          -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,
                          0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(grid,col="black",pch=19)})
dev.off()
pdf("gridvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(gridrun[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(grid,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

pdf("gridimpnopoints.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(grid,col="black",pch=19)})
dev.off()

pdf("gridimp1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(grid,col="black",pch=19)
                 points(gridbh,col="black",pch=4)})
dev.off()

length(which(gridImat>=3))/1600

pdf("gridimp2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad,col="black",pch=19)
                 points(gridbh2,col="black",pch=4)})
dev.off()

length(which(gridImat2>=3))/1600

pdf("gridimp3.pdf",width = 7.6,height = 7)
gridbh3 = matrix(c(0.5,5.783185,-0.25,-0.25),2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad2,col="black",pch=19)
                 points(gridbh3,col="black",pch=4)})
dev.off()

length(which(gridImat3>=3))/1600

gridad3 = as.matrix(rbind(gridad2,gridbh3))
gridrun4 = emufuncosxplussiny(0,0.6,2,gridad3)
gridI4 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                MoreArgs = list(z=1,points=gridad3))
gridImat4 = matrix(gridI4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad3,col="black",pch=19)})
gridbh4 = matrix(c(1.5,4.95,0.25,0.25),2)
pdf("gridimp4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad3,col="black",pch=19)
                 points(gridbh4,col="black",pch=4)})
dev.off()

length(which(gridImat4>=3))/1600

gridad4 = as.matrix(rbind(gridad3,gridbh4))
gridrun5 = emufuncosxplussiny(0,0.6,2,gridad4)
gridI5 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                MoreArgs = list(z=1,points=gridad4))
gridImat5 = matrix(gridI5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad4,col="black",pch=19)})
gridbh5 = matrix(c(1.5,4.95,2.75,2.75),2)
pdf("gridimp5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad4,col="black",pch=19)
                 points(gridbh5,col="black",pch=4)})
dev.off()

length(which(gridImat5>=3))/1600

gridad5 = as.matrix(rbind(gridad4,gridbh5))
gridrun6 = emufuncosxplussiny(0,0.6,2,gridad5)
gridI6 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                MoreArgs = list(z=1,points=gridad5))
gridImat6 = matrix(gridI6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat6, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad5,col="black",pch=19)})
gridbh6 = matrix(c(0,6.25,6.5,6.5),2)
pdf("gridimp6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat6, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(gridImat)),
               plot.title = title(main=expression(paste(I(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad5,col="black",pch=19)
                 points(gridbh6,col="black",pch=4)})
dev.off()

length(which(gridImat6>=3))/1600

#########################################################################################

pdf("maximinexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsrun[[2]],
               col = rainbow(25,alpha=0.95),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(j,col="black",pch=19)})
dev.off()
pdf("maximinvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(maxminlhsrun[[3]]),
               col = rainbow(13,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(j,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

pdf("maximinimpnopoints.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(j,col="black",pch=19)})
dev.off()

length(which(maxminlhsImat>=3))/1600

pdf("maximinimp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(j,col="black",pch=19)
                 points(jnew[10:13,],col="black",pch=4)})
dev.off()

length(which(maxminlhsImat>=3))/1600

pdf("maximinimp2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat2)),
               plot.title = title(main=expression(paste(I(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew,col="black",pch=19)
                 points(jnew2[14:16,],col="black",pch=4)})
dev.off()

length(which(maxminlhsImat2>=3))/1600

pdf("maximinimp3.pdf",width = 7.6,height = 7)
j4 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j4))){
    j4=r
  }
  else{
    j4=j4
  } 
}
posmaxmin3 = mapply(impcosxplussiny,as.data.frame(t(j4)),
                    MoreArgs = list(z=1,points=jnew2))
which(posmaxmin3<3)
jnew3 = as.matrix(rbind(jnew2,j4[1,],j4[6,]))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat3)),
               plot.title = title(main=expression(paste(I(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew2,col="black",pch=19)
                 points(jnew3[17:18,],col="black",pch=4)})
dev.off()

length(which(maxminlhsImat3>=3))/1600

maxminlhsrun4 = emufuncosxplussiny(0,0.6,2,jnew3)
maxminlhsI4 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=jnew3))
maxminlhsImat4 = matrix(maxminlhsI4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat4)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew3,col="black",pch=19)})
j5 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j5))){
    j5=r
  }
  else{
    j5=j5
  } 
}
posmaxmin4 = mapply(impcosxplussiny,as.data.frame(t(j5)),
                    MoreArgs = list(z=1,points=jnew3))
which(posmaxmin4<3)
jnew4 = as.matrix(rbind(jnew3,j5[3,]))
pdf("maximinimp4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat4)),
               plot.title = title(main=expression(paste(I(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew3,col="black",pch=19)
                 points(jnew4[19,1],jnew4[19,2],col="black",pch=4)})
dev.off()

length(which(maxminlhsImat4>=3))/1600

maxminlhsrun5 = emufuncosxplussiny(0,0.6,2,jnew4)
maxminlhsI5 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=jnew4))
maxminlhsImat5 = matrix(maxminlhsI5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat5)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew4,col="black",pch=19)})
j6 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j6))){
    j6=r
  }
  else{
    j6=j6
  } 
}
posmaxmin5 = mapply(impcosxplussiny,as.data.frame(t(j6)),
                    MoreArgs = list(z=1,points=jnew4))
which(posmaxmin5<3)
jnew5 = as.matrix(rbind(jnew4,j6[4,],j6[6,]))
pdf("maximinimp5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat5)),
               plot.title = title(main=expression(paste(I(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew4,col="black",pch=19)
                 points(jnew5[20:21,],col="black",pch=4)})
dev.off()

length(which(maxminlhsImat5>=3))/1600

maxminlhsrun6 = emufuncosxplussiny(0,0.6,2,jnew5)
maxminlhsI6 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=jnew5))
maxminlhsImat6 = matrix(maxminlhsI6,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat6, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat6)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew5,col="black",pch=19)})
j7 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j7))){
    j7=r
  }
  else{
    j7=j7
  } 
}
posmaxmin6 = mapply(impcosxplussiny,as.data.frame(t(j7)),
                    MoreArgs = list(z=1,points=jnew5))
which(posmaxmin6<3)
jnew6 = as.matrix(rbind(jnew5,j7[3,],j7[6,]))
pdf("maximinimp6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat6, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat6)),
               plot.title = title(main=expression(paste(I(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew5,col="black",pch=19)
                 points(jnew6[22:23,],col="black",pch=4)})
dev.off()

length(which(maxminlhsImat6>=3))/1600

################################################################################################################

pdf("maxoptimvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(moptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(mopoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

######################################################################################################

pdf("avgoptimvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(avoptvarrun1[[3]]),
               col = rainbow(30,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

############################################################################################################

pdf("avgoptimsevexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarrun1[[2]],
               col = rainbow(50,alpha=0.95),
               levels = c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,
                          -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,
                          0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)})
dev.off()
pdf("avgoptimsevvar.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun1[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

pdf("avgoptimsevimpnopoints.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)})
dev.off()

pdf("avgoptimsevimp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat)),
               plot.title = title(main=expression(paste(I(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)
                 points(matrix(optim1points.1$par,nrow=4,ncol=2),col="black",pch=4)})
dev.off()

length(which(sevoptvarImat>=3))/1600

pdf("avgoptimsevimp2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.1)),
               plot.title = title(main=expression(paste(I(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.1,col="black",pch=19)
                 points(matrix(optim1points.2$par,nrow=4,ncol=2),col="black",pch=4)})
dev.off()

length(which(sevoptvarImat.1>=3))/1600

pdf("avgoptimsevimp3.pdf",width = 7.6,height = 7)
length(which(sevoptvarImat.2<3))
sevoptvarIless.2 = as.vector(which(sevoptvarImat.2<3))
sevoptvarpointIless.2 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                    y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevoptvarIless.2,]

sevoptvarpointIless.2
sevminpointind.2 = c(which(as.vector(sevoptvarI.2)==sort(as.vector(sevoptvarI.2))[1]),
                     which(as.vector(sevoptvarI.2)==sort(as.vector(sevoptvarI.2))[3]))
sevminpoint.2 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                            y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevminpointind.2,]
sevminpoint.2
optim1points.3 = optim(c(sevminpoint.2[,1],sevminpoint.2[,2]),avgvarformany.2,numpoints=2,ogpoints=osevpoints.2,greenspace=sevoptvarpointIless.2,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optim1points.3
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim1points.2$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=osevpoints.1)))
matrix(optim1points.2$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.2,col="black",pch=19)
                 points(matrix(optim1points.3$par,nrow=2,ncol=2),col="black",pch=4)})
matrix(optim1points.3$par,nrow=2,ncol=2)

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.2)),
               plot.title = title(main=expression(paste(I(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.2,col="black",pch=19)
                 points(matrix(optim1points.3$par,nrow=2,ncol=2),col="black",pch=4)})
dev.off()

length(which(sevoptvarImat.2>=3))/1600

osevpoints.3 = rbind(osevpoints.2,matrix(optim1points.3$par,nrow=2,ncol=2))
osevpoints.3
sevoptvarrun4 = emufuncosxplussiny(0,0.6,2,osevpoints.3)
sevoptvarI.3 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=osevpoints.3))
sevoptvarImat.3 = matrix(sevoptvarI.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.3)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.3,col="black",pch=19)})
length(which(sevoptvarImat.3<3))
sevoptvarIless.3 = as.vector(which(sevoptvarImat.3<3))
sevoptvarpointIless.3 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                    y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevoptvarIless.3,]

sevoptvarpointIless.3
sevminpointind.3 = c(which(as.vector(sevoptvarI.3)==sort(as.vector(sevoptvarI.3))[5]),
                     which(as.vector(sevoptvarI.3)==sort(as.vector(sevoptvarI.3))[6]))
sevminpoint.3 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                            y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevminpointind.3,]
sevminpoint.3
optim1points.4 = optim(c(sevminpoint.3[,1],sevminpoint.3[,2]),avgvarformany.2,numpoints=2,ogpoints=osevpoints.3,greenspace=sevoptvarpointIless.3,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optim1points.4
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim1points.4$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=osevpoints.3)))
matrix(optim1points.4$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.3)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.3,col="black",pch=19)
                 points(matrix(optim1points.4$par,nrow=2,ncol=2),col="black",pch=4)})
matrix(optim1points.4$par,nrow=2,ncol=2)

length(which(sevoptvarImat.3>=3))/1600


pdf("avgoptimsevimp4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.3)),
               plot.title = title(main=expression(paste(I(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.3,col="black",pch=19)
                 points(matrix(optim1points.4$par,nrow=2,ncol=2),col="black",pch=4)})
dev.off()

length(which(sevoptvarImat.3>=3))/1600

osevpoints.4 = rbind(osevpoints.3,matrix(optim1points.4$par,nrow=2,ncol=2))
osevpoints.4
sevoptvarrun5 = emufuncosxplussiny(0,0.6,2,osevpoints.4)
sevoptvarI.4 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=osevpoints.4))
sevoptvarImat.4 = matrix(sevoptvarI.4,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.4)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.4,col="black",pch=19)})
length(which(sevoptvarImat.4<3))
sevoptvarIless.4 = as.vector(which(sevoptvarImat.4<3))
sevoptvarpointIless.4 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                    y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevoptvarIless.4,]

sevoptvarpointIless.4
sevminpointind.4 = c(which(as.vector(sevoptvarI.4)==sort(as.vector(sevoptvarI.4))[5]),
                     which(as.vector(sevoptvarI.4)==sort(as.vector(sevoptvarI.4))[6]))
sevminpoint.4 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                            y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevminpointind.4,]
sevminpoint.4
c(sevminpoint.4[,1],sevminpoint.4[,2])
optim1points.5 = optim(c(4.789219439,1.489656549,2*pi+0.5,2*pi+0.5),avgvarformany.2,numpoints=2,ogpoints=osevpoints.4,greenspace=sevoptvarpointIless.4,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optim1points.5
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim1points.5$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=osevpoints.4)))
matrix(optim1points.5$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.4)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.4,col="black",pch=19)
                 points(matrix(optim1points.5$par,nrow=2,ncol=2),col="black",pch=4)})


pdf("avgoptimsevimp5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.4, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat5)),
               plot.title = title(main=expression(paste(I(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.4,col="black",pch=19)
                 points(matrix(optim1points.5$par,nrow=2,ncol=2),col="black",pch=4)})
dev.off()

length(which(sevoptvarImat.4>=3))/1600

osevpoints.5 = rbind(osevpoints.4,matrix(optim1points.5$par,nrow=2,ncol=2))
osevpoints.5
sevoptvarrun6 = emufuncosxplussiny(0,0.6,2,osevpoints.5)
sevoptvarI.5 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=osevpoints.5))
sevoptvarImat.5 = matrix(sevoptvarI.5,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.5)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.5,col="black",pch=19)})
length(which(sevoptvarImat.5<3))
sevoptvarIless.5 = as.vector(which(sevoptvarImat.5<3))
sevoptvarpointIless.5 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                    y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevoptvarIless.5,]

sevoptvarpointIless.5
sevminpointind.5 = c(which(as.vector(sevoptvarI.5)==sort(as.vector(sevoptvarI.5))[5]),
                     which(as.vector(sevoptvarI.5)==sort(as.vector(sevoptvarI.5))[6]))
sevminpoint.5 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                            y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevminpointind.5,]
sevminpoint.5
c(sevminpoint.5[,1],sevminpoint.5[,2])
optim1points.6 = optim(c(sevminpoint.5[,1],sevminpoint.5[,2]),avgvarformany.2,numpoints=2,ogpoints=osevpoints.5,greenspace=sevoptvarpointIless.5,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optim1points.6
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim1points.6$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=osevpoints.5)))

matrix(optim1points.6$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.5)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.5,col="black",pch=19)
                 points(matrix(optim1points.6$par,nrow=2,ncol=2),col="black",pch=4)})

pdf("avgoptimsevimp6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.5, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.5)),
               plot.title = title(main=expression(paste(I(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.5,col="black",pch=19)
                 points(matrix(optim1points.6$par,nrow=2,ncol=2),col="black",pch=4)})
dev.off()

length(which(sevoptvarImat.5>=3))/1600

####################################################################################################

pdf("avgoptimsevvar2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun2[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
pdf("avgoptimsevvar3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun3[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
pdf("avgoptimsevvar4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun4[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.3,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
pdf("avgoptimsevvar5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun5[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.4,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()
pdf("avgoptimsevvar6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun6[[3]]),
               col = rainbow(33,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.5,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

################################################################################################

pdf("sineemulator.pdf")
newrun8p = emufunnew8p(0,0.6,0.06,seq(0.1,0.51,length.out = 8))
par(mfrow=c(1,1))
plot(newrun8p[[1]],newrun8p[[2]],type='l',col='blue',xlab="x",ylab="f(x)",ylim=c(-1.5,1.5),main=expression("Sinusoidal Model Emulator"))
lines(newrun8p[[1]],newrun8p[[5]],type='l',col='red')
lines(newrun8p[[1]],newrun8p[[6]],type='l',col='red')
lines(sequemsinexp,realxsin,type='l',col='black')
points(newrun8p[[3]],newrun8p[[4]],col='purple',pch=19)
dev.off()

pdf("EI1D.pdf")
plot(sequem1d,-sequemnew1d,type='l',col='darkgoldenrod',xlab="x",ylab="EI(x)", main = expression("Expected Improvement"))
points(newrun8p[[3]],rep(0,8),col='purple',pch=19)
points(as.vector(possnewax1D[[1]]),-as.vector(possnewax1D[[2]]),col='purple',pch=4,lwd=2)
axis(1,at=c(0.46))
dev.off()

pdf("sineemulator2.pdf")
plot(newrun8p.1[[1]],newrun8p.1[[2]],type='l',col='blue',xlab="x",ylab="f(x)",ylim=c(-1.5,1.5),main=expression("Wave 2 Sinusoidal Model Emulator"))
lines(newrun8p.1[[1]],newrun8p.1[[5]],type='l',col='red')
lines(newrun8p.1[[1]],newrun8p.1[[6]],type='l',col='red')
lines(sequemsinexp,realxsin,type='l',col='black')
points(newrun8p.1[[3]],newrun8p.1[[4]],col='purple',pch=19)
dev.off()

##################################################################################################

crp.rg <- colorRampPalette(c("purple3","darkblue","blue","cyan","green","yellow","darkorange"))
colors = crp.rg(23)
plot(seq(-0.5,(2*pi)+0.5,length.out=40),seq(-0.5,(2*pi)+0.5,length.out=40),cex=0.2,pch=16,col=colors)

pdf("9pointEIwave1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat,
               col = colors,
               plot.title = title(main = expression(paste(EI(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oavgpoints,col="black",pch=19)
                 points(matrix(possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpoints
avoptvarrun1acc[[5]]
dist(matrix(c(0.6713430,0.1776135,0.6729991,1.2450531),2))

colors2 = crp.rg(13)
pdf("9pointEIwave2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.1,
               col = colors2,
               plot.title = title(main = expression(paste(EI(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints,col="black",pch=19)
                 points(matrix(possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpoints.1
avoptvarrun1acc.1[[5]]
dist(matrix(c(0.1776135,-0.2221756,1.2450531,1.4251843),2))

colors3 = crp.rg(26)
pdf("9pointEIwave3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.2,
               col = colors3,
               plot.title = title(main = expression(paste(EI(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.1,col="black",pch=19)
                 points(matrix(possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpoints.2
avoptvarrun1acc.2[[5]]
dist(matrix(c(0.002017447,-0.2221756,1.6255691,1.4251843),2))

sequemacc.3 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpoints.2))
sequemaccmat.3 = matrix(sequemacc.3,40)
colors4 = crp.rg(20)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.3,
               col = colors4,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.2,col="black",pch=19)})
startpoint.3 = sequem[which(sequemacc.3==min(sequemacc.3)),]
startpoint.3
possnewmax.3 = optim(startpoint.3,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpoints.2,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmax.3
pdf("9pointEIwave4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.3,
               col = colors4,
               plot.title = title(main = expression(paste(EI(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.2,col="black",pch=19)
                 points(matrix(possnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpoints.3 = rbind(newmaxpoints.2,possnewmax.3$par)
newmaxpoints.3
avoptvarrun1acc.3 = emufuncosxplussiny(0,0.6,2,newmaxpoints.3)
avoptvarrun1acc.3[[5]]
dist(matrix(c(0.002017447,6.595335567,1.6255691,1.4558197),2))

sequemacc.4 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpoints.3))
sequemaccmat.4 = matrix(sequemacc.4,40)
colorsadd = crp.rg(18)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.4,
               col = colorsadd,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.3,col="black",pch=19)})
startpoint.4 = sequem[which(sequemacc.4==min(sequemacc.4)),]
startpoint.4
possnewmax.4 = optim(startpoint.4,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpoints.3,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmax.4
pdf("9pointEIwave5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.4,
               col = colorsadd,
               plot.title = title(main = expression(paste(EI(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.3,col="black",pch=19)
                 points(matrix(possnewmax.4$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpoints.4 = rbind(newmaxpoints.3,possnewmax.4$par)
newmaxpoints.4
avoptvarrun1acc.4 = emufuncosxplussiny(0,0.6,2,newmaxpoints.4)
avoptvarrun1acc.4[[5]]
dist(matrix(c(0.002017447,6.783185307,1.6255691,0.8442755),2))

sequemacc.5 = mapply(expectedimprovement2D,as.data.frame(t(sequem)),MoreArgs = list(B0=0,sigmau=0.6,theta=2,points=newmaxpoints.4))
sequemaccmat.5 = matrix(sequemacc.5,40)
colors5 = crp.rg(15)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.5,
               col = colors5,
               plot.title = title(main = expression(paste("Expected Improvement for Emulator with ",beta[0]==0,", ",sigma[u]==0.6,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.4,col="black",pch=19)})
startpoint.5 = sequem[which(sequemacc.5==min(sequemacc.5)),]
startpoint.5
possnewmax.5 = optim(startpoint.5,expectedimprovement2D,B0=0,sigmau=0.6,theta=2,points=newmaxpoints.4,method="L-BFGS-B",lower=rep(-0.5,2),upper=rep(2*pi+0.5,2))
possnewmax.5
pdf("9pointEIwave6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -sequemaccmat.5,
               col = colors5,
               plot.title = title(main = expression(paste(EI(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.4,col="black",pch=19)
                 points(matrix(possnewmax.5$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newmaxpoints.5 = rbind(newmaxpoints.4,possnewmax.5$par)
newmaxpoints.5
avoptvarrun1acc.5 = emufuncosxplussiny(0,0.6,2,newmaxpoints.5)
avoptvarrun1acc.5[[5]]

##########################################################################################

pdf("9pointEIexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = avoptvarrun1acc.5[[2]],
               col = rainbow(50,alpha=0.95),
               levels = c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,
                          -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,
                          0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(newmaxpoints.5,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

###################################################################################################

pdf("5pointEIexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarrun1[[2]],
               col = rainbow(30,alpha=0.95),
               levels = seq(-2,2,length.out = 25),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

################################################################################################

colors6 = crp.rg(17)
pdf("5pointEIwave1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest,
               col = colors6,
               plot.title = title(main = expression(paste(EI(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest
fivEIoptvarrun2[[5]]
dist(matrix(c(5.137809,1.145419,1.1455374,1.1453421),2))

pdf("5pointEIwave2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.1,
               col = colors6,
               plot.title = title(main = expression(paste(EI(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.1$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.1
fivEIoptvarrun3[[5]]
dist(matrix(c(5.137809,5.443303,1.1455374,0.5777148),2))

pdf("5pointEIwave3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.2,
               col = colors2,
               plot.title = title(main = expression(paste(EI(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.1,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.2
fivEIoptvarrun4[[5]]
dist(matrix(c(5.137809,4.4776330,1.1455374,0.6381480),2))

pdf("5pointEIwave4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.3,
               col = colors4,
               plot.title = title(main = expression(paste(EI(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.2,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.3$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.3
fivEIoptvarrun5[[5]]
dist(matrix(c(5.8560075,1.1454187,1.4167532,1.1453421),2))

pdf("5pointEIwave5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.4,
               col = colors6,
               plot.title = title(main = expression(paste(EI(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.3,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.4$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.4
fivEIoptvarrun6[[5]]
dist(matrix(c(5.8560075,6.0820602,1.4167532,1.5511127),2))

pdf("5pointEIwave6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.5,
               col = colors5,
               plot.title = title(main = expression(paste(EI(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.4,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.5$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.5
fivEIoptvarrun7[[5]]
dist(matrix(c(6.1676477,6.0820602,1.9423108,1.5511127),2))

pdf("5pointEIwave7.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.6,
               col = colors,
               plot.title = title(main = expression(paste(EI(x)," for Wave 7 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.5,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.6$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.6
fivEIoptvarrun8[[5]]
dist(matrix(c(6.1676477,6.3296800,1.9423108,1.6013315),2))

pdf("5pointEIwave8.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatEItest.7,
               col = colors4,
               plot.title = title(main = expression(paste(EI(x)," for Wave 8 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsEItest.6,col="black",pch=19)
                 points(matrix(fivpossnewmaxEItest.7$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsEItest.7
fivEIoptvarrun9 = emufuncosxplussiny(0,0.6,2,newofivpointsEItest.7)
fivEIoptvarrun9[[5]]
dist(matrix(c(-0.003774649,6.3296800,1.7012676,1.6013315),2))

#######################################################################################################

pdf("5pointemexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivEIoptvarrun9[[2]],
               col = rainbow(30,alpha=0.95),
               levels = seq(-2,2,length.out = 25),
               plot.axes = {axis(1,seq(from=0,to=10,by=2))
                 axis(2,seq(from=0,to=10,by=2))
                 points(newofivpointsEItest.7,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

###########################################################################################################

pdf("MC5pointEIwave1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test,
               col = colors5,
               plot.title = title(main = expression(paste(EI^{3},(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsMC3test
fivMC3testrun2[[5]]
dist(matrix(c(5.137809,5.392545,1.1455374,0.1126828),2))

pdf("MC5pointEIwave2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.2,
               col = colors5,
               plot.title = title(main = expression(paste(EI^{3},(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsMC3test.2
fivMC3testrun3[[5]]
dist(matrix(c(5.137809,0.7185467,1.1455374,0.1761819),2))

pdf("MC5pointEIwave3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.3,
               col = colors3,
               plot.title = title(main = expression(paste(EI^{3},(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.2,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.3$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsMC3test.3
fivMC3testrun4[[5]]
dist(matrix(c(6.1452905,1.1454187,1.6560560,1.1453421),2))

pdf("MC5pointEIwave4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.4,
               col = colors3,
               plot.title = title(main = expression(paste(EI^{3},(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.3,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.4$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsMC3test.4
fivMC3testrun5[[5]]
dist(matrix(c(6.1452905,6.4672573,1.6560560,2.4343509),2))

colors9 = crp.rg(21)
pdf("MC5pointEIwave5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.5,
               col = colors9,
               plot.title = title(main = expression(paste(EI^{3},(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.4,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.5$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsMC3test.5
fivMC3testrun6[[5]]
dist(matrix(c(6.1452905,6.7831853,1.6560560,1.4853566),2))

pdf("MC5pointEIwave6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = -fivsequemaccmatMC3test.6,
               col = colors,
               plot.title = title(main = expression(paste(EI^{3},(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=9,by=1))
                 axis(2,seq(from=0,to=9,by=1))
                 points(newofivpointsMC3test.5,col="black",pch=19)
                 points(matrix(fivpossnewmaxMC3test.6$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
newofivpointsMC3test.6
fivMC3testrun7[[5]]
dist(matrix(c(6.1452905,-0.3132137,1.6560560,1.7409963),2))

##################################################################################################

pdf("MC5pointemexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivMC3testrun7[[2]],
               col = rainbow(25,alpha=0.95),
               levels = seq(-2,2,length.out = 21),
               plot.axes = {axis(1,seq(from=0,to=15,by=2))
                 axis(2,seq(from=0,to=15,by=2))
                 points(newofivpointsMC3test.6,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()               

#################################################################################################################

pdf("truerobotarm.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=100),
               y = seq(0,1,length.out=100), 
               z = zrob,
               col = rainbow(25,alpha=0.95),
               plot.title = title(main = expression(paste("Robot Arm Function ")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
dev.off()

pdf("robotexp.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = robotarm3run1[[2]],
               col = rainbow(25,alpha=0.95),
               levels = seq(2,4,length.out = 21),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
dev.off()

##############################################################################################

pdf("roboteiwave1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat,
               col = colors6,
               plot.title = title(main = expression(paste(EI(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)
                 points(matrix(RA3possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.1
robotarm3run2[[5]]

pdf("roboteiwave2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.1,
               col = colors5,
               plot.title = title(main = expression(paste(EI(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.1,col="black",pch=19)
                 points(matrix(RA3possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.2
robotarm3run3[[5]]

pdf("roboteiwave3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.2,
               col = colors8,
               plot.title = title(main = expression(paste(EI(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.2,col="black",pch=19)
                 points(matrix(RA3possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.3
robotarm3run4[[5]]

pdf("roboteiwave4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.3,
               col = colors5,
               plot.title = title(main = expression(paste(EI(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.3,col="black",pch=19)
                 points(matrix(RA3possnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.4
robotarm3run5[[5]]

pdf("roboteiwave5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.4,
               col = colors9,
               plot.title = title(main = expression(paste(EI(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.4,col="black",pch=19)
                 points(matrix(RA3possnewmax.4$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.5
robotarm3run6[[5]]

pdf("roboteiwave6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.5,
               col = colors6,
               plot.title = title(main = expression(paste(EI(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.5,col="black",pch=19)
                 points(matrix(RA3possnewmax.5$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.6
robotarm3run7[[5]]

colors10 = crp.rg(24)
pdf("roboteiwave7.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.6,
               col = colors10,
               plot.title = title(main = expression(paste(EI(x)," for Wave 7 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.6,col="black",pch=19)
                 points(matrix(RA3possnewmax.6$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.7
robotarm3run8[[5]]

colors11 = crp.rg(19)
pdf("roboteiwave8.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.7,
               col = colors11,
               plot.title = title(main = expression(paste(EI(x)," for Wave 8 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.7,col="black",pch=19)
                 points(matrix(RA3possnewmax.7$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.8
robotarm3run9[[5]]

pdf("roboteiwave9.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccmat.8,
               col = colors10,
               plot.title = title(main = expression(paste(EI(x)," for Wave 9 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.8,col="black",pch=19)
                 points(matrix(RA3possnewmax.8$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
robot3points.9
robotarm3run10[[5]]

##################################################################################################

colors12 = crp.rg(27)
pdf("MCroboteiwave1.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat,
               col = colors12,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Initial Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(robot3points.0,col="black",pch=19)
                 points(matrix(MCRA3possnewmax$par,nrow=1,ncol=2),col="black",pch=4)})

dev.off()

pdf("MCroboteiwave2.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.1,
               col = colors10,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Wave 2 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.1,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.1$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("MCroboteiwave3.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.2,
               col = colors3,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Wave 3 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.2,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.2$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("MCroboteiwave4.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.3,
               col = colors8,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Wave 4 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.3,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.3$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("MCroboteiwave5.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.4,
               col = colors4,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Wave 5 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.4,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.4$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("MCroboteiwave6.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.5,
               col = colors11,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Wave 6 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.5,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.5$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()

pdf("MCroboteiwave7.pdf",width = 7.6,height = 7)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = -robot3sequemaccMCmat.6,
               col = colors8,
               plot.title = title(main = expression(paste(EI^{4},(x)," for Wave 7 Emulator")),
                                  xlab = expression(x[1]), ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.2))
                 axis(2,seq(from=0,to=1,by=0.2))
                 points(MCrobot3points.6,col="black",pch=19)
                 points(matrix(MCRA3possnewmax.6$par,nrow=1,ncol=2),col="black",pch=4)})
dev.off()
