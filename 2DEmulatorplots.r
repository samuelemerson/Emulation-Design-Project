install.packages('lhs')
library(lhs)

########################################################################
pdf("2DEmulatorplots.pdf",width = 7.6,height = 7)

cosxplussinyrun1 = emufuncosxplussiny(0,0.25,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1[[3]]),
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==2.0)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))

cosxplussinyrun2 = emufuncosxplussiny(0,0.25,0.1,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==0.1)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun2[[3]]),
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==0.1)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))

cosxplussinyrun3 = emufuncosxplussiny(0,0.25,5,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==5)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun3[[3]]),
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.25,", ",theta==5)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))

cosxplussinyrun4 = emufuncosxplussiny(0,0.01,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun4[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation for ",beta[0]==0,", ",sigma[u]==0.01,", ",theta==2)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun4[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.01,", ",theta==2)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))

cosxplussinyrun5 = emufuncosxplussiny(0,0.65,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun5[[3]]),
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.65,", ",theta==2)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))

cosxplussinyrun6 = emufuncosxplussiny(0,0.65,5,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun6[[3]]),
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation for ",beta[0]==0,", ",sigma[u]==0.65,", ",theta==5)),
                                  xlab = expression(x[1]), ylab = expression(x[2])))

dev.off()
