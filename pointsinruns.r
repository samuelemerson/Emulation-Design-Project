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
length(which(tenoptvarImat<3))
tenoptvarIless = as.vector(which(tenoptvarImat<3))
tenoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                  y=seq(-0.5,(2*pi)+0.5,length.out=40))[tenoptvarIless,]

tenoptvarpointIless
tenminpointind = c(which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[1]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[2]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[3]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[4]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[5]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[6]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[7]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[8]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[9]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[10]))
tenminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                          y=seq(-0.5,(2*pi)+0.5,length.out=40))[tenminpointind,]
tenminpoint
optim10points.1 = optim(c(tenminpoint[,1],tenminpoint[,2]),avgvarformany.2,numpoints=10,ogpoints=otenpoints,greenspace=tenoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,20),upper=rep(2*pi+0.5,20))
optim10points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim10points.1$par,nrow=10,ncol=2))),
                 MoreArgs = list(z=1,points=otenpoints)))
matrix(optim10points.1$par,nrow=10,ncol=2)
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
                 points(otenpoints,col="black",pch=19)
                 points(matrix(optim10points.1$par,nrow=10,ncol=2),col="black",pch=4)})

otenpoints.1 = rbind(otenpoints,matrix(optim10points.1$par,nrow=10,ncol=2))
otenpoints.1
tenoptvarrun2 = emufuncosxplussiny(0,0.6,2,otenpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(tenoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
tenoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=otenpoints.1))
tenoptvarImat.1 = matrix(tenoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(tenoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.1,col="black",pch=19)})
length(which(tenoptvarImat.1<3))

##############################################################################

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

length(which(fivoptvarImat<3))
fivoptvarIless = as.vector(which(fivoptvarImat<3))
fivoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                  y=seq(-0.5,(2*pi)+0.5,length.out=40))[fivoptvarIless,]

fivoptvarpointIless

fivminpointind = c(which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[1]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[2]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[3]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[4]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[5]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[6]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[7]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[8]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[9]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[10]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[11]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[12]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[13]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[14]),
                   which(as.vector(fivoptvarI)==sort(as.vector(fivoptvarI))[15]))
fivminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                          y=seq(-0.5,(2*pi)+0.5,length.out=40))[fivminpointind,]
fivminpoint
optim5points.1 = optim(c(fivminpoint[,1],fivminpoint[,2]),avgvarformany.2,numpoints=15,ogpoints=ofivpoints,greenspace=fivoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,30),upper=rep(2*pi+0.5,30))
optim5points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim5points.1$par,nrow=15,ncol=2))),
                 MoreArgs = list(z=1,points=ofivpoints)))
matrix(optim5points.1$par,nrow=15,ncol=2)
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
                 points(ofivpoints,col="black",pch=19)
                 points(matrix(optim5points.1$par,nrow=15,ncol=2),col="black",pch=4)})

ofivpoints.1 = rbind(ofivpoints,matrix(optim5points.1$par,nrow=15,ncol=2))
ofivpoints.1
fivoptvarrun2 = emufuncosxplussiny(0,0.6,2,ofivpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fivoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofivpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
fivoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=ofivpoints.1))
fivoptvarImat.1 = matrix(fivoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fivoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(fivoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofivpoints.1,col="black",pch=19)})
length(which(fivoptvarImat.1<3))

############################################################################

fifpoints = as.matrix(2*pi*maximinLHS(15,2,optimize.on = "grid"))
fifpoints
optim15points = optim(as.vector(fifpoints),avgvarformany,numpoints=15,method="L-BFGS-B",lower=rep(-0.5,30),upper=rep(2*pi+0.5,30))
ofifpoints = matrix(optim15points$par,nrow=15,ncol=2)
ofifpoints
fifoptvarrun1 = emufuncosxplussiny(0,0.6,2,ofifpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fifoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fifoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofifpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
fifoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=ofifpoints))
fifoptvarImat = matrix(fifoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fifoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(fifoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofifpoints,col="black",pch=19)})

length(which(fifoptvarImat<3))
fifoptvarIless = as.vector(which(fifoptvarImat<3))
fifoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                  y=seq(-0.5,(2*pi)+0.5,length.out=40))[fifoptvarIless,]

fifoptvarpointIless

fifminpointind = c(which(as.vector(fifoptvarI)==sort(as.vector(fifoptvarI))[1]),
                   which(as.vector(fifoptvarI)==sort(as.vector(fifoptvarI))[2]),
                   which(as.vector(fifoptvarI)==sort(as.vector(fifoptvarI))[3]),
                   which(as.vector(fifoptvarI)==sort(as.vector(fifoptvarI))[4]),
                   which(as.vector(fifoptvarI)==sort(as.vector(fifoptvarI))[5]))
fifminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                          y=seq(-0.5,(2*pi)+0.5,length.out=40))[fifminpointind,]
fifminpoint
optim15points.1 = optim(c(fifminpoint[,1],fifminpoint[,2]),avgvarformany.2,numpoints=5,ogpoints=ofifpoints,greenspace=fifoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,10),upper=rep(2*pi+0.5,10))
optim15points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim15points.1$par,nrow=5,ncol=2))),
                 MoreArgs = list(z=1,points=ofifpoints)))
matrix(optim15points.1$par,nrow=5,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fifoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(fifoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofifpoints,col="black",pch=19)
                 points(matrix(optim15points.1$par,nrow=5,ncol=2),col="black",pch=4)})

ofifpoints.1 = rbind(ofifpoints,matrix(optim15points.1$par,nrow=5,ncol=2))
ofifpoints.1
fifoptvarrun2 = emufuncosxplussiny(0,0.6,2,ofifpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fifoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(fifoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofifpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
fifoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=ofifpoints.1))
fifoptvarImat.1 = matrix(fifoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = fifoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(fifoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(ofifpoints.1,col="black",pch=19)})
length(which(fifoptvarImat.1<3))

###################################################################################################

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

length(which(twooptvarImat<3))
twooptvarIless = as.vector(which(twooptvarImat<3))
twooptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                  y=seq(-0.5,(2*pi)+0.5,length.out=40))[twooptvarIless,]

twooptvarpointIless

twominpointind = c(which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[1]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[2]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[3]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[4]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[5]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[6]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[7]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[8]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[9]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[10]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[11]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[12]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[13]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[14]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[15]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[16]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[17]),
                   which(as.vector(twooptvarI)==sort(as.vector(twooptvarI))[18]))
twominpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                          y=seq(-0.5,(2*pi)+0.5,length.out=40))[twominpointind,]
twominpoint
optim2points.1 = optim(c(twominpoint[,1],twominpoint[,2]),avgvarformany.2,numpoints=18,ogpoints=otwopoints,greenspace=twooptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,36),upper=rep(2*pi+0.5,36))
optim2points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim2points.1$par,nrow=18,ncol=2))),
                 MoreArgs = list(z=1,points=otwopoints)))
matrix(optim2points.1$par,nrow=18,ncol=2)
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
                 points(otwopoints,col="black",pch=19)
                 points(matrix(optim2points.1$par,nrow=18,ncol=2),col="black",pch=4)})

otwopoints.1 = rbind(otwopoints,matrix(optim2points.1$par,nrow=18,ncol=2))
otwopoints.1
twooptvarrun2 = emufuncosxplussiny(0,0.6,2,otwopoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = twooptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(twooptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
twooptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=otwopoints.1))
twooptvarImat.1 = matrix(twooptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = twooptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(twooptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otwopoints.1,col="black",pch=19)})
length(which(twooptvarImat.1<3))

##################################################################################################


eighteenpoints = as.matrix(2*pi*maximinLHS(18,2,optimize.on = "grid"))
eighteenpoints
optim18points = optim(as.vector(eighteenpoints),avgvarformany,numpoints=18,method="L-BFGS-B",lower=rep(-0.5,36),upper=rep(2*pi+0.5,36))
oeighteenpoints = matrix(optim18points$par,nrow=18,ncol=2)
oeighteenpoints
eighteenoptvarrun1 = emufuncosxplussiny(0,0.6,2,oeighteenpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eighteenoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(eighteenoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeighteenpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
eighteenoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=oeighteenpoints))
eighteenoptvarImat = matrix(eighteenoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eighteenoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(eighteenoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeighteenpoints,col="black",pch=19)})

length(which(eighteenoptvarImat<3))
eighteenoptvarIless = as.vector(which(eighteenoptvarImat<3))
eighteenoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                       y=seq(-0.5,(2*pi)+0.5,length.out=40))[eighteenoptvarIless,]

eighteenoptvarpointIless

eighteenminpointind = c(which(as.vector(eighteenoptvarI)==sort(as.vector(eighteenoptvarI))[1]),
                        which(as.vector(eighteenoptvarI)==sort(as.vector(eighteenoptvarI))[2]))
eighteenminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                               y=seq(-0.5,(2*pi)+0.5,length.out=40))[eighteenminpointind,]
eighteenminpoint
optim18points.1 = optim(c(eighteenminpoint[,1],eighteenminpoint[,2]),avgvarformany.2,numpoints=2,ogpoints=oeighteenpoints,greenspace=eighteenoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,4),upper=rep(2*pi+0.5,4))
optim18points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim18points.1$par,nrow=2,ncol=2))),
                 MoreArgs = list(z=1,points=oeighteenpoints)))
matrix(optim18points.1$par,nrow=2,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eighteenoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(eighteenoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeighteenpoints,col="black",pch=19)
                 points(matrix(optim18points.1$par,nrow=2,ncol=2),col="black",pch=4)})

oeighteenpoints.1 = rbind(oeighteenpoints,matrix(optim18points.1$par,nrow=2,ncol=2))
oeighteenpoints.1
eighteenoptvarrun2 = emufuncosxplussiny(0,0.6,2,oeighteenpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eighteenoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(eighteenoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeighteenpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
eighteenoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                           MoreArgs = list(z=1,points=oeighteenpoints.1))
eighteenoptvarImat.1 = matrix(eighteenoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = eighteenoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(eighteenoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(oeighteenpoints.1,col="black",pch=19)})
length(which(eighteenoptvarImat.1<3))

##########################################################################################################

#10,10,10

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
length(which(tenoptvarImat<3))
tenoptvarIless = as.vector(which(tenoptvarImat<3))
tenoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                  y=seq(-0.5,(2*pi)+0.5,length.out=40))[tenoptvarIless,]

tenoptvarpointIless

tenminpointind = c(which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[1]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[2]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[3]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[4]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[5]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[6]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[7]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[8]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[9]),
                   which(as.vector(tenoptvarI)==sort(as.vector(tenoptvarI))[10]))
tenminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                          y=seq(-0.5,(2*pi)+0.5,length.out=40))[tenminpointind,]
tenminpoint
optim10points.1 = optim(c(tenminpoint[,1],tenminpoint[,2]),avgvarformany.2,numpoints=10,ogpoints=otenpoints,greenspace=tenoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,20),upper=rep(2*pi+0.5,20))
optim10points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim10points.1$par,nrow=10,ncol=2))),
                 MoreArgs = list(z=1,points=otenpoints)))
matrix(optim10points.1$par,nrow=10,ncol=2)
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
                 points(otenpoints,col="black",pch=19)
                 points(matrix(optim10points.1$par,nrow=10,ncol=2),col="black",pch=4)})

otenpoints.1 = rbind(otenpoints,matrix(optim10points.1$par,nrow=10,ncol=2))
otenpoints.1
tenoptvarrun2 = emufuncosxplussiny(0,0.6,2,otenpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(tenoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
tenoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=otenpoints.1))
tenoptvarImat.1 = matrix(tenoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(tenoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.1,col="black",pch=19)})
length(which(tenoptvarImat.1<3))
tenoptvarIless.1 = as.vector(which(tenoptvarImat.1<3))
tenoptvarpointIless.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                    y=seq(-0.5,(2*pi)+0.5,length.out=40))[tenoptvarIless.1,]

tenoptvarpointIless.1

tenminpointind.1 = c(which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[1]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[2]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[3]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[4]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[5]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[6]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[7]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[8]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[9]),
                     which(as.vector(tenoptvarI.1)==sort(as.vector(tenoptvarI.1))[10]))
tenminpoint.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                            y=seq(-0.5,(2*pi)+0.5,length.out=40))[tenminpointind.1,]
tenminpoint.1
optim10points.2 = optim(c(tenminpoint.1[,1],tenminpoint.1[,2]),avgvarformany.2,numpoints=10,ogpoints=otenpoints.1,greenspace=tenoptvarpointIless.1,method="L-BFGS-B",lower=rep(-0.5,20),upper=rep(2*pi+0.5,20))
optim10points.2
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim10points.2$par,nrow=10,ncol=2))),
                 MoreArgs = list(z=1,points=otenpoints.1)))
matrix(optim10points.2$par,nrow=10,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(tenoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.1,col="black",pch=19)
                 points(matrix(optim10points.2$par,nrow=10,ncol=2),col="black",pch=4)})

otenpoints.2 = rbind(otenpoints.1,matrix(optim10points.2$par,nrow=10,ncol=2))
otenpoints.2
tenoptvarrun3 = emufuncosxplussiny(0,0.6,2,otenpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(tenoptvarrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
tenoptvarI.2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=1,points=otenpoints.2))
tenoptvarImat.2 = matrix(tenoptvarI.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = tenoptvarImat.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(tenoptvarImat.2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(otenpoints.2,col="black",pch=19)})
length(which(tenoptvarImat.2<3))

#############################################################################################################

#13,12,5

thirteenpoints = as.matrix(2*pi*maximinLHS(13,2,optimize.on = "grid"))
thirteenpoints
optim13points = optim(as.vector(thirteenpoints),avgvarformany,numpoints=13,method="L-BFGS-B",lower=rep(-0.5,26),upper=rep(2*pi+0.5,26))
othirteenpoints = matrix(optim13points$par,nrow=13,ncol=2)
othirteenpoints
thirteenoptvarrun1 = emufuncosxplussiny(0,0.6,2,othirteenpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(thirteenoptvarrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
thirteenoptvarI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=othirteenpoints))
thirteenoptvarImat = matrix(thirteenoptvarI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(thirteenoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints,col="black",pch=19)})
length(which(thirteenoptvarImat<3))
thirteenoptvarIless = as.vector(which(thirteenoptvarImat<3))
thirteenoptvarpointIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                       y=seq(-0.5,(2*pi)+0.5,length.out=40))[thirteenoptvarIless,]

thirteenoptvarpointIless

thirteenminpointind = c(which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[1]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[2]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[3]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[4]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[5]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[6]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[7]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[8]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[9]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[10]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[11]),
                        which(as.vector(thirteenoptvarI)==sort(as.vector(thirteenoptvarI))[12]))
thirteenminpoint = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                               y=seq(-0.5,(2*pi)+0.5,length.out=40))[thirteenminpointind,]
thirteenminpoint
optim13points.1 = optim(c(thirteenminpoint[,1],thirteenminpoint[,2]),avgvarformany.2,numpoints=12,ogpoints=othirteenpoints,greenspace=thirteenoptvarpointIless,method="L-BFGS-B",lower=rep(-0.5,24),upper=rep(2*pi+0.5,24))
optim13points.1
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim13points.1$par,nrow=12,ncol=2))),
                 MoreArgs = list(z=1,points=othirteenpoints)))
matrix(optim13points.1$par,nrow=12,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(thirteenoptvarImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints,col="black",pch=19)
                 points(matrix(optim13points.1$par,nrow=12,ncol=2),col="black",pch=4)})

othirteenpoints.1 = rbind(othirteenpoints,matrix(optim13points.1$par,nrow=12,ncol=2))
othirteenpoints.1
thirteenoptvarrun2 = emufuncosxplussiny(0,0.6,2,othirteenpoints.1)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(thirteenoptvarrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints.1,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
thirteenoptvarI.1 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                           MoreArgs = list(z=1,points=othirteenpoints.1))
thirteenoptvarImat.1 = matrix(thirteenoptvarI.1,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(thirteenoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints.1,col="black",pch=19)})
length(which(thirteenoptvarImat.1<3))
thirteenoptvarIless.1 = as.vector(which(thirteenoptvarImat.1<3))
thirteenoptvarpointIless.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                         y=seq(-0.5,(2*pi)+0.5,length.out=40))[thirteenoptvarIless.1,]

thirteenoptvarpointIless.1

thirteenminpointind.1 = c(which(as.vector(thirteenoptvarI.1)==sort(as.vector(thirteenoptvarI.1))[1]),
                          which(as.vector(thirteenoptvarI.1)==sort(as.vector(thirteenoptvarI.1))[2]),
                          which(as.vector(thirteenoptvarI.1)==sort(as.vector(thirteenoptvarI.1))[3]),
                          which(as.vector(thirteenoptvarI.1)==sort(as.vector(thirteenoptvarI.1))[4]),
                          which(as.vector(thirteenoptvarI.1)==sort(as.vector(thirteenoptvarI.1))[5]))
thirteenminpoint.1 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                                 y=seq(-0.5,(2*pi)+0.5,length.out=40))[thirteenminpointind.1,]
thirteenminpoint.1
optim13points.2 = optim(c(thirteenminpoint.1[,1],thirteenminpoint.1[,2]),avgvarformany.2,numpoints=5,ogpoints=othirteenpoints.1,greenspace=thirteenoptvarpointIless.1,method="L-BFGS-B",lower=rep(-0.5,10),upper=rep(2*pi+0.5,10))
optim13points.2
as.vector(mapply(impcosxplussiny,as.data.frame(t(matrix(optim13points.2$par,nrow=5,ncol=2))),
                 MoreArgs = list(z=1,points=othirteenpoints.1)))
matrix(optim13points.2$par,nrow=5,ncol=2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarImat.1, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(thirteenoptvarImat.1)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints.1,col="black",pch=19)
                 points(matrix(optim13points.2$par,nrow=5,ncol=2),col="black",pch=4)})

othirteenpoints.2 = rbind(othirteenpoints.1,matrix(optim13points.2$par,nrow=5,ncol=2))
othirteenpoints.2
thirteenoptvarrun3 = emufuncosxplussiny(0,0.6,2,othirteenpoints.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(thirteenoptvarrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
thirteenoptvarI.2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                           MoreArgs = list(z=1,points=othirteenpoints.2))
thirteenoptvarImat.2 = matrix(thirteenoptvarI.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = thirteenoptvarImat.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(thirteenoptvarImat.2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(othirteenpoints.2,col="black",pch=19)})
length(which(thirteenoptvarImat.2<3))
