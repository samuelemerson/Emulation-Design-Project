realfI = function(x,z){
  r = (((cos(x[1])+sin(x[2]))-z)^2)/0.0016
  return(r)
}

sequemreal = mapply(realfI,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1))
sequemrealmat = matrix(sequemreal,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemrealmat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemrealmat)),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"))

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sequemImat, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(sequemImat)),
               plot.title = title(main = "Implausability measure I",
                                  xlab = "x",ylab = "y"),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(latpoints,col="black",pch=19)})

posnewpoints = mapply(impcosxplussiny,as.data.frame(t(latpoints2)),
                      MoreArgs = list(z=1,points=latpoints))
which(posnewpoints<3)
latpointsnew = as.matrix(rbind(latpoints,latpoints2[3,]))
