lp2 = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
lp3 = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
lp4 = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))

#####################################################################

pdf("2DIbyhand.pdf",width = 7.6,height = 7)

f = function(x){
  r = cos(x[1])+sin(x[2])
  return(r)
}
realfI = function(x,z){
  r = (((cos(x[1])+sin(x[2]))-z)^2)/0.0016
  return(r)
}
sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                     y=seq(-0.5,(2*pi)+0.5,length.out=40))
sequemf = mapply(f,as.data.frame(t(sequem)))
sequemfmat = matrix(sequemf,40)
sequemreal = mapply(realfI,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1))
sequemrealmat = matrix(sequemreal,40)
length(which(sequemrealmat<3))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemfmat,
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Real Function ", cos(x[1])+sin(x[2]))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemrealmat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemrealmat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

length(which(sequemrealmat>=3))/1600

######################################################################

cosxplussinyrun1 = emufuncosxplussiny(0,0.6,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cosxplussinyrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cosxplussinyrun1[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)})
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)
                 points(lp2,col="black",pch=4)})
posnewp = mapply(impcosxplussiny,as.data.frame(t(lp2)),
                      MoreArgs = list(z=1,points=latpoints))
which(posnewp<3)
lpnew = as.matrix(rbind(latpoints,lp2[1,],lp2[3,],lp2[6,],lp2[8,]))
#####################################################################

cxpsyr1new = emufuncosxplussiny(0,0.6,2,lpnew)
sInew = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=lpnew))
sImnew = matrix(sInew,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cxpsyr1new[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cxpsyr1new[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImnew)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew,col="black",pch=19)})
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImnew)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew,col="black",pch=19)
                 points(lp3,col="black",pch=4)})
posnewp2 = mapply(impcosxplussiny,as.data.frame(t(lp3)),
                       MoreArgs = list(z=1,points=lpnew))
which(posnewp2<3)
lpnew2 = as.matrix(rbind(lpnew,lp3[5,],lp3[6,],lp3[9,]))

#####################################################################

cxpsyr1new2 = emufuncosxplussiny(0,0.6,2,lpnew2)
sInew2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=lpnew2))
sImnew2 = matrix(sInew2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cxpsyr1new2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cxpsyr1new2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImnew2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImnew2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(lpnew2,col="black",pch=19)})
posnewp3 = mapply(impcosxplussiny,as.data.frame(t(lp4)),
                  MoreArgs = list(z=1,points=lpnew2))
which(posnewp3<3)
lpnew3 = as.matrix(rbind(lpnew2,lp4[6,],lp4[8,]))

#########################################################################

byhp = matrix(c(0.45,0.2,6.25,2,0.2,3),nrow=3,ncol=2)
byhpad = as.matrix(rbind(latpoints,byhp))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)})
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)
                 points(byhp,col="black",pch=4)})

##################################################################################

cxpsyr1bh = emufuncosxplussiny(0,0.6,2,byhpad)
sIbh = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=byhpad))
sImbh = matrix(sIbh,40)
byhp2 = matrix(c(1.5,5,6.25,0,0,4.5),nrow=3,ncol=2)
byhpad2 = as.matrix(rbind(byhpad,byhp2))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cxpsyr1bh[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cxpsyr1bh[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImbh, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad,col="black",pch=19)})
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImbh, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad,col="black",pch=19)
                 points(byhp2,col="black",pch=4)})
######################################################################

cxpsyr1bh2 = emufuncosxplussiny(0,0.6,2,byhpad2)
sIbh2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
              MoreArgs = list(z=1,points=byhpad2))
sImbh2 = matrix(sIbh2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = cxpsyr1bh2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(cxpsyr1bh2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sImbh2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(byhpad2,col="black",pch=19)})

dev.off()
