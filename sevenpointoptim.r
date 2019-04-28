sevpoints = as.matrix(2*pi*maximinLHS(7,2,optimize.on = "grid"))
sevpoints
optim1points = optim(as.vector(sevpoints),avgvarformany,numpoints=7,method="L-BFGS-B",lower=rep(-0.5,14),upper=rep(2*pi+0.5,14))
osevpoints = matrix(optim1points$par,nrow=7,ncol=2)
osevpoints
sevoptvarrun1 = emufuncosxplussiny(0,0.6,2,osevpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
sevoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                   MoreArgs = list(z=1,points=osevpoints))
sevoptvarImat = matrix(sevoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)})
length(which(sevoptvarImat<3))
sevoptvarIless = as.vector(which(sevoptvarImat<3))
sevoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                 y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevoptvarIless,]

sevoptvarpointIless

#sevadpoi = sevoptvarpointIless[sample(nrow(sevoptvarpointIless),size=4),]
#sevadpoi
sevminpointind = c(which(as.vector(sevoptvarI)==sort(as.vector(sevoptvarI))[1]),
                   which(as.vector(sevoptvarI)==sort(as.vector(sevoptvarI))[2]),
                   which(as.vector(sevoptvarI)==sort(as.vector(sevoptvarI))[3]),
                   which(as.vector(sevoptvarI)==sort(as.vector(sevoptvarI))[5]))
sevminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                          y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevminpointind,]
sevminpoint
optim1points.1 = optim(c(sevminpoint[,1],sevminpoint[,2]),avgvarformany.2,numpoints=4,ogpoints=osevpoints,greenspace=sevoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
optim1points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim1points.1$par,nrow=4,ncol=2))),
                 MoreArgs = list(z=1,points=osevpoints)))
matrix(optim1points.1$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints,col="black",pch=19)
                 points(matrix(optim1points.1$par,nrow=4,ncol=2),col="black",pch=4)})

osevpoints.1 = rbind(osevpoints,matrix(optim1points.1$par,nrow=4,ncol=2))
osevpoints.1
sevoptvarrun2 = emufuncosxplussiny(0,0.6,2,osevpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
sevoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=osevpoints.1))
sevoptvarImat.1 = matrix(sevoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.1,col="black",pch=19)})
length(which(sevoptvarImat.1<3))
sevoptvarIless.1 = as.vector(which(sevoptvarImat.1<3))
sevoptvarpointIless.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                    y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevoptvarIless.1,]

sevoptvarpointIless.1

#sevadpoi.1 = sevoptvarpointIless.1[sample(nrow(sevoptvarpointIless.1),size=2),]
#sevadpoi.1
sevminpointind.1 = c(which(as.vector(sevoptvarI.1)==sort(as.vector(sevoptvarI.1))[1]),
                     which(as.vector(sevoptvarI.1)==sort(as.vector(sevoptvarI.1))[2]),
                     which(as.vector(sevoptvarI.1)==sort(as.vector(sevoptvarI.1))[3]),
                     which(as.vector(sevoptvarI.1)==sort(as.vector(sevoptvarI.1))[4]))
sevminpoint.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                            y=seq(-0.5,(2*pi)+0.5,length.out=40))[sevminpointind.1,]
sevminpoint.1
c(sevminpoint.1[,1],sevminpoint.1[,2])
c(-0.5,-0.5,2*pi+0.5,2*pi+0.5,-0.5,2*pi+0.5,-0.5,2*pi +0.5)
optim1points.2 = optim(c(sevminpoint.1[,1],sevminpoint.1[,2]),avgvarformany.2,numpoints=4,ogpoints=osevpoints.1,greenspace=sevoptvarpointIless.1,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
optim1points.2
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim1points.2$par,nrow=4,ncol=2))),
                 MoreArgs = list(z=1,points=osevpoints.1)))
matrix(optim1points.2$par,nrow=4,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sevoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.1,col="black",pch=19)
                 points(matrix(optim1points.2$par,nrow=4,ncol=2),col="black",pch=4)})

osevpoints.2 = rbind(osevpoints.1,matrix(optim1points.2$par,nrow=4,ncol=2))
osevpoints.2
sevoptvarrun3 = emufuncosxplussiny(0,0.6,2,osevpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sevoptvarrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sevoptvarrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osevpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
sevoptvarI.2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=osevpoints.2))
sevoptvarImat.2 = matrix(sevoptvarI.2,40)
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
                 points(osevpoints.2,col="black",pch=19)})
#####################################################################################################

twopoints = as.matrix(2*pi*maximinLHS(2,2,optimize.on = "grid"))
twopoints
optim2points = optim(as.vector(twopoints),avgvarformany,numpoints=2,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
otwopoints = matrix(optim2points$par,nrow=2,ncol=2)
otwopoints
twooptvarrun1 = emufuncosxplussiny(0,0.6,2,otwopoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = twooptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(twooptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
twooptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=otwopoints))
twooptvarImat = matrix(twooptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = twooptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(twooptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints,col="black",pch=19)})
#####################################################################################

thrpoints = as.matrix(2*pi*maximinLHS(3,2,optimize.on = "grid"))
thrpoints
optim3points = optim(as.vector(thrpoints),avgvarformany,numpoints=3,method="L-BFGS-B",lower=rep(-0.5,6),upper=rep(2*pi+0.5,6))
othrpoints = matrix(optim3points$par,nrow=3,ncol=2)
othrpoints
throptvarrun1 = emufuncosxplussiny(0,0.6,2,othrpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = throptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(throptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othrpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
throptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=othrpoints))
throptvarImat = matrix(throptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = throptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(throptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othrpoints,col="black",pch=19)})

#####################################################################################################


foupoints = as.matrix(2*pi*maximinLHS(4,2,optimize.on = "grid"))
foupoints
optim4points = optim(as.vector(foupoints),avgvarformany,numpoints=4,method="L-BFGS-B",lower=rep(-0.5,8),upper=rep(2*pi+0.5,8))
ofoupoints = matrix(optim4points$par,nrow=4,ncol=2)
ofoupoints
fouoptvarrun1 = emufuncosxplussiny(0,0.6,2,ofoupoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fouoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fouoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofoupoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
fouoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=ofoupoints))
fouoptvarImat = matrix(fouoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fouoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(fouoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofoupoints,col="black",pch=19)})

################################################################################

fivpoints = as.matrix(2*pi*maximinLHS(5,2,optimize.on = "grid"))
fivpoints
optim5points = optim(as.vector(fivpoints),avgvarformany,numpoints=5,method="L-BFGS-B",lower=rep(-0.5,10),upper=rep(2*pi+0.5,10))
ofivpoints = matrix(optim5points$par,nrow=5,ncol=2)
ofivpoints
fivoptvarrun1 = emufuncosxplussiny(0,0.6,2,ofivpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofivpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
fivoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=ofivpoints))
fivoptvarImat = matrix(fivoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(fivoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofivpoints,col="black",pch=19)})

##############################################################################

sixpoints = as.matrix(2*pi*maximinLHS(6,2,optimize.on = "grid"))
sixpoints
optim6points = optim(as.vector(sixpoints),avgvarformany,numpoints=6,method="L-BFGS-B",lower=rep(-0.5,12),upper=rep(2*pi+0.5,12))
osixpoints = matrix(optim6points$par,nrow=6,ncol=2)
osixpoints
sixoptvarrun1 = emufuncosxplussiny(0,0.6,2,osixpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sixoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(sixoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osixpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
sixoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=osixpoints))
sixoptvarImat = matrix(sixoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sixoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sixoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(osixpoints,col="black",pch=19)})

#############################################################################

eigpoints = as.matrix(2*pi*maximinLHS(8,2,optimize.on = "grid"))
eigpoints
optim8points = optim(as.vector(eigpoints),avgvarformany,numpoints=8,method="L-BFGS-B",lower=rep(-0.5,16),upper=rep(2*pi+0.5,16))
oeigpoints = matrix(optim8points$par,nrow=8,ncol=2)
oeigpoints
eigoptvarrun1 = emufuncosxplussiny(0,0.6,2,oeigpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eigoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(eigoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeigpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
eigoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=oeigpoints))
eigoptvarImat = matrix(eigoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eigoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(eigoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeigpoints,col="black",pch=19)})

##############################################################################################

tenpoints = as.matrix(2*pi*maximinLHS(10,2,optimize.on = "grid"))
tenpoints
optim10points = optim(as.vector(tenpoints),avgvarformany,numpoints=10,method="L-BFGS-B",lower=rep(-0.5,20),upper=rep(2*pi+0.5,20))
otenpoints = matrix(optim10points$par,nrow=10,ncol=2)
otenpoints
tenoptvarrun1 = emufuncosxplussiny(0,0.6,2,otenpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(tenoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
tenoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=otenpoints))
tenoptvarImat = matrix(tenoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(tenoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints,col="black",pch=19)})
