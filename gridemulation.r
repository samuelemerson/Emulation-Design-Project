pdf("2Dinitialpointchoices.pdf",width = 7.6,height = 7)

sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                     y=seq(-0.5,(2*pi)+0.5,length.out=40))
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
               plot.title = title(expression(paste("Emulator Standard Deviation")),
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

#######################################################################
grid = as.matrix(expand.grid(x=seq(0.5,(2*pi)-0.5,length.out=3),
                   y=seq(0.5,(2*pi)-0.5,length.out=3)))
gridrun = emufuncosxplussiny(0,0.6,2,grid)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridrun[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(gridrun[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(grid,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
gridI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                 MoreArgs = list(z=1,points=grid))
gridImat = matrix(gridI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(grid,col="black",pch=19)})
gridbh = matrix(c(0.5,0.5,5.783185,5.783185,2,4.5,2,4.5),4)
gridad = as.matrix(rbind(grid,gridbh))
gridrun2 = emufuncosxplussiny(0,0.6,2,gridad)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(gridrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
gridI2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=gridad))
gridImat2 = matrix(gridI2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad,col="black",pch=19)})
gridbh2 = matrix(c(1.5,4.95,1.5,1.5),2)
gridad2 = as.matrix(rbind(gridad,gridbh2))
gridrun3 = emufuncosxplussiny(0,0.6,2,gridad2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(gridrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
gridI3 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                MoreArgs = list(z=1,points=gridad2))
gridImat3 = matrix(gridI3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = gridImat3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sImbh)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(gridad2,col="black",pch=19)})
########################################################################
j = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j))){
    j=r
  }
  else{
    j=j
  } 
}
maxminlhsrun = emufuncosxplussiny(0,0.6,2,j)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsrun[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(maxminlhsrun[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(j,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
maxminlhsI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=j))
maxminlhsImat = matrix(maxminlhsI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(j,col="black",pch=19)})
j2 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j2))){
    j2=r
  }
  else{
    j2=j2
  } 
}
posmaxmin = mapply(impcosxplussiny,as.data.frame(t(j2)),
                   MoreArgs = list(z=1,points=j))
which(posmaxmin<3)
jnew = as.matrix(rbind(j,j2[1,],j2[2,],j2[4,],j2[6,]))
maxminlhsrun2 = emufuncosxplussiny(0,0.6,2,jnew)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(maxminlhsrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
maxminlhsI2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                MoreArgs = list(z=1,points=jnew))
maxminlhsImat2 = matrix(maxminlhsI2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew,col="black",pch=19)})
j3 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(j3))){
    j3=r
  }
  else{
    j3=j3
  } 
}
posmaxmin2 = mapply(impcosxplussiny,as.data.frame(t(j3)),
                   MoreArgs = list(z=1,points=jnew))
which(posmaxmin2<3)
jnew2 = as.matrix(rbind(jnew,j3[3,],j3[6,],j3[8,]))
maxminlhsrun3 = emufuncosxplussiny(0,0.6,2,jnew2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(maxminlhsrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
maxminlhsI3 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=jnew2))
maxminlhsImat3 = matrix(maxminlhsI3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = maxminlhsImat3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(maxminlhsImat3)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(jnew2,col="black",pch=19)})

length(which(maxminlhsImat3>=3))/1600

##########################################################################

k = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(k))){
    k=r
  }
  else{
    k=k
  } 
}
kmaxminlhsrun = emufuncosxplussiny(0,0.6,2,k)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = kmaxminlhsrun[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(kmaxminlhsrun[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
kmaxminlhsI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=k))
kmaxminlhsImat = matrix(kmaxminlhsI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = kmaxminlhsImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(kmaxminlhsImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k,col="black",pch=19)})
k2 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(k2))){
    k2=r
  }
  else{
    k2=k2
  } 
}
kposmaxmin = mapply(impcosxplussiny,as.data.frame(t(k2)),
                   MoreArgs = list(z=1,points=k))
which(kposmaxmin<3)
knew = as.matrix(rbind(k,k2[8,]))
kmaxminlhsrun2 = emufuncosxplussiny(0,0.25,2,knew)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = kmaxminlhsrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(kmaxminlhsrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(knew,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
kmaxminlhsI2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=knew))
kmaxminlhsImat2 = matrix(kmaxminlhsI2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = kmaxminlhsImat2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(kmaxminlhsImat2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(knew,col="black",pch=19)})
k3 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(k3))){
    k3=r
  }
  else{
    k3=k3
  } 
}
kposmaxmin2 = mapply(impcosxplussiny,as.data.frame(t(k3)),
                    MoreArgs = list(z=1,points=knew))
which(kposmaxmin2<3)
knew2 = as.matrix(rbind(knew,k3[9,]))
kmaxminlhsrun3 = emufuncosxplussiny(0,0.25,2,knew2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = kmaxminlhsrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(kmaxminlhsrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(knew2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
kmaxminlhsI3 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=knew2))
kmaxminlhsImat3 = matrix(kmaxminlhsI3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = kmaxminlhsImat3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(kmaxminlhsImat3)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(knew2,col="black",pch=19)})

dev.off()

#################################################################################

pdf("2Dinitialpointchoicesz=0.pdf",width = 7.6,height = 7)

sequemreal.z0 = mapply(realfI,as.data.frame(t(sequem)),
                    MoreArgs = list(z=0))
sequemreal.z0mat = matrix(sequemreal.z0,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemfmat,
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Real Function ", sin(x[1])+cos(x[2]))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemreal.z0mat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemreal.z0mat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

k.z0 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(k.z0))){
    k.z0=r
  }
  else{
    k.z0=k.z0
  } 
}
k.z0maxminlhsrun = emufuncosxplussiny(0,0.6,2,k.z0)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = k.z0maxminlhsrun[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(k.z0maxminlhsrun[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k.z0,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
k.z0maxminlhsI = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=0,points=k.z0))
k.z0maxminlhsImat = matrix(k.z0maxminlhsI,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = k.z0maxminlhsImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(k.z0maxminlhsImat)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k.z0,col="black",pch=19)})
k.z02 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(k.z02))){
    k.z02=r
  }
  else{
    k.z02=k.z02
  } 
}
k.z0posmaxmin = mapply(impcosxplussiny,as.data.frame(t(k.z02)),
                    MoreArgs = list(z=0,points=k.z0))
which(k.z0posmaxmin<3)
k.z0new = as.matrix(rbind(k.z0,k.z02[2,],k.z02[5,],k.z02[9,]))
k.z0maxminlhsrun2 = emufuncosxplussiny(0,0.25,2,k.z0new)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = k.z0maxminlhsrun2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(k.z0maxminlhsrun2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k.z0new,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
k.z0maxminlhsI2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=0,points=k.z0new))
k.z0maxminlhsImat2 = matrix(k.z0maxminlhsI2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = k.z0maxminlhsImat2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(k.z0maxminlhsImat2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k.z0new,col="black",pch=19)})
k.z03 = matrix(rep(0,18),nrow=9,ncol=2)
for (i in seq(100000)){
  r = as.matrix(2*pi*maximinLHS(9,2,optimize.on = "grid"))
  if(min(dist(r))>min(dist(k.z03))){
    k.z03=r
  }
  else{
    k.z03=k.z03
  } 
}
k.z0posmaxmin2 = mapply(impcosxplussiny,as.data.frame(t(k.z03)),
                     MoreArgs = list(z=0,points=k.z0new))
which(k.z0posmaxmin2<3)
k.z0new2 = as.matrix(rbind(k.z0new,k.z03[4,],k.z03[5,],k.z03[6,]))
k.z0maxminlhsrun3 = emufuncosxplussiny(0,0.25,2,k.z0new2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = k.z0maxminlhsrun3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(k.z0maxminlhsrun3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k.z0new2,col="black",pch=19)},
               plot.title = title(main = expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
k.z0maxminlhsI3 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                      MoreArgs = list(z=0,points=k.z0new2))
k.z0maxminlhsImat3 = matrix(k.z0maxminlhsI3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = k.z0maxminlhsImat3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(k.z0maxminlhsImat3)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(k.z0new2,col="black",pch=19)})

dev.off()



#pdf plots 7*8 and then scale up by 50%...or down by 50% to fix white line problem
#generate latin hypercube with 100000000000 or bigger do waves etc.
#generate 1000 hypercubes, clculate emulator variance of these
#choose average one and one with min max variance...show wave one plots
#consider wave 2 choices 
#optim