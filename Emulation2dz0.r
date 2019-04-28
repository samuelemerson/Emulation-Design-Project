filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemImatnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemImatnew2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))
###############################################################################################

cosxplussinyrun1 = emufuncosxplussiny(0,0.25,2,latpoints)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = cosxplussinyrun1[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = "Emulator Expectation for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x", ylab = "y"))
contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
        y = seq(-0.5,(2*pi)+0.5,length.out=100), 
        z = cosxplussinyrun1[[2]])
contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
        y = seq(-0.5,(2*pi)+0.5,length.out=100), 
        z = cosxplussinyrun1[[3]])
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = cosxplussinyrun1[[3]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)},
               plot.title = title(main = "Emulator Variance for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x",ylab = "y"))
persp(seq(-0.5,(2*pi)+0.5,length.out=100), 
      seq(-0.5,(2*pi)+0.5,length.out=100), 
      cosxplussinyrun1[[2]], theta = 30, phi = 30, expand = 0.5,
      xlab = "x",ylab = "y",zlab = "z", col = "lightblue")

sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=100),
                     y=seq(-0.5,(2*pi)+0.5,length.out=100))
sequemI0 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                 MoreArgs = list(z=0,points=latpoints))
sequemI0mat = matrix(sequemI0,100)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemI0mat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

persp(seq(-0.5,(2*pi)+0.5,length.out=100), 
      seq(-0.5,(2*pi)+0.5,length.out=100), 
      sequemI0mat, theta = 30, phi = 30, expand = 0.5, 
      col = "lightblue")

posnewpoints0 = mapply(impcosxplussiny,as.data.frame(t(latpoints2)),
                      MoreArgs = list(z=0,points=latpoints))
which(posnewpoints0<3)
latpointsnew0 = as.matrix(rbind(latpoints,latpoints2[5,],
                               latpoints2[7,],latpoints2[8,],
                               latpoints2[9,]))

#########################################################################

cosxplussinyrun1new0 = emufuncosxplussiny(0,0.25,2,latpointsnew0)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = cosxplussinyrun1new0[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = "Emulator Expectation for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x", ylab = "y"))
contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
        y = seq(-0.5,(2*pi)+0.5,length.out=100), 
        z = cosxplussinyrun1new0[[2]])
contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
        y = seq(-0.5,(2*pi)+0.5,length.out=100), 
        z = cosxplussinyrun1new0[[3]])
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = cosxplussinyrun1new0[[3]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointsnew0,col="black",pch=19)},
               plot.title = title(main = "Emulator Variance for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x",ylab = "y"))
persp(seq(-0.5,(2*pi)+0.5,length.out=100), 
      seq(-0.5,(2*pi)+0.5,length.out=100), 
      cosxplussinyrun1new[[2]], theta = 30, phi = 30, expand = 0.5,
      xlab = "x",ylab = "y",zlab = "z", col = "lightblue")

sequemI0new = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=0,points=latpointsnew0))
sequemI0matnew = matrix(sequemI0new,100)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemI0matnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

persp(seq(-0.5,(2*pi)+0.5,length.out=100), 
      seq(-0.5,(2*pi)+0.5,length.out=100), 
      sequemI0matnew, theta = 30, phi = 30, expand = 0.5, 
      col = "lightblue")

posnewpoints20 = mapply(impcosxplussiny,as.data.frame(t(latpoints3)),
                       MoreArgs = list(z=0,points=latpointsnew0))
which(posnewpoints20<3)
latpointsnew20 = as.matrix(rbind(latpointsnew0,latpoints3[3,],
                                latpoints3[6,],latpoints3[8,],
                                latpoints3[9,]))

####################################################################

cosxplussinyrun1new20 = emufuncosxplussiny(0,0.25,2,latpointsnew20)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = cosxplussinyrun1new20[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = "Emulator Expectation for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x", ylab = "y"))
contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
        y = seq(-0.5,(2*pi)+0.5,length.out=100), 
        z = cosxplussinyrun1new20[[2]])
contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
        y = seq(-0.5,(2*pi)+0.5,length.out=100), 
        z = cosxplussinyrun1new20[[3]])
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = cosxplussinyrun1new20[[3]],
               col = rainbow(25,alpha=0.75),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpointsnew20,col="black",pch=19)},
               plot.title = title(main = "Emulator Variance for B0=0,Sigmau=0.25,Theta=2.0",
                                  xlab = "x",ylab = "y"))
persp(seq(-0.5,(2*pi)+0.5,length.out=100), 
      seq(-0.5,(2*pi)+0.5,length.out=100), 
      cosxplussinyrun1new20[[2]], theta = 30, phi = 30, expand = 0.5,
      xlab = "x",ylab = "y",zlab = "z", col = "lightblue")

sequemI0new2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=0,points=latpointsnew20))
sequemI0matnew2 = matrix(sequemI0new2,100)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemI0matnew2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

persp(seq(-0.5,(2*pi)+0.5,length.out=100), 
      seq(-0.5,(2*pi)+0.5,length.out=100), 
      sequemI0matnew2, theta = 30, phi = 30, expand = 0.5, 
      col = "lightblue")

##################################################################################

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemI0mat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemI0matnew, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=100),
               y = seq(-0.5,(2*pi)+0.5,length.out=100), 
               z = sequemI0matnew2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))
