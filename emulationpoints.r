layout(matrix(c(1,3,2,4), 2, 2, byrow = TRUE))
par(mar=c(4,4,2,2))

run1 = emufunexp(3.5,1.5,0.14)
plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14",cex.main=0.65)
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)
abline(h=3.5,col='black')
abline(h=(3.5+2*sqrt(0.025)),col='black',lty=2)
abline(h=(3.5-2*sqrt(0.025)),col='black',lty=2)
#axis(2,at=3.5,cex.axis=0.5)

sequem = seq(0.05,0.55,length.out=10000)
sequem2 = c()
for (i in sequem){
  sequem2 = append(sequem2,impexp(i,3.5))
}
plot(sequem,sequem2,type='l',col='green',xlab='x',ylab='I')
abline(h=3,col='darkgoldenrod')
points(c(0.317,0.396),c(impexp(0.317,3.5),impexp(0.396,3.5)),col='purple',pch=19)
#axis(1,at=0.317,cex.axis=0.5,tick=FALSE,line=-5)
#axis(1,at=0.396,cex.axis=0.5,tick=FALSE,line=-5)


#layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
#par(mar=c(4,4,2,2))

newrun = emufunexp2(3.5,1.5,0.14,c(0.1,0.2,0.3,0.317,0.396,0.4,0.5))
plot(newrun[[1]],newrun[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14",cex.main=0.65)
lines(newrun[[1]],newrun[[5]],type='l',lwd=1,col='red')
lines(newrun[[1]],newrun[[6]],type='l',lwd=1,col='red')
points(newrun[[3]],newrun[[4]],col='purple',pch=19)
abline(h=3.5,col='black')
abline(h=(3.5+2*sqrt(0.025)),col='black',lty=2)
abline(h=(3.5-2*sqrt(0.025)),col='black',lty=2)

newpointsrun1 = overallsequem(c(0.1,0.2,0.3,0.317,0.396,0.4,0.5))  
plot(newpointsrun1[[1]],newpointsrun1[[2]],type='l',col='green',xlab='x',ylab='I',ylim=c(0,14))
abline(h=3,col='darkgoldenrod')
########################################################################

layout(matrix(c(1,3,2,4), 2, 2, byrow = TRUE))
par(mar=c(4,4,2,2))

run1 = emufunexp(3.5,1.5,0.14)
plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14",cex.main=0.65)
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)
abline(h=3.5,col='black')
abline(h=(3.5+2*sqrt(0.025)),col='black',lty=2)
abline(h=(3.5-2*sqrt(0.025)),col='black',lty=2)

plot(sequem,sequem2,type='l',col='green',xlab='x',ylab='I')
abline(h=3,col='darkgoldenrod')
points(c(0.342,0.369),c(impexp(0.342,3.5),impexp(0.369,3.5)),col='purple',pch=19)
#axis(1,at=0.317,cex.axis=0.5,tick=FALSE,line=-5)
#axis(1,at=0.396,cex.axis=0.5,tick=FALSE,line=-5)

newrun = emufunexp2(3.5,1.5,0.14,c(0.1,0.2,0.3,0.342,0.369,0.4,0.5))
plot(newrun[[1]],newrun[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14",cex.main=0.65)
lines(newrun[[1]],newrun[[5]],type='l',lwd=1,col='red')
lines(newrun[[1]],newrun[[6]],type='l',lwd=1,col='red')
points(newrun[[3]],newrun[[4]],col='purple',pch=19)
abline(h=3.5,col='black')
abline(h=(3.5+2*sqrt(0.025)),col='black',lty=2)
abline(h=(3.5-2*sqrt(0.025)),col='black',lty=2)

newpointsrun1 = overallsequem(c(0.1,0.2,0.3,0.342,0.369,0.4,0.5))  
plot(newpointsrun1[[1]],newpointsrun1[[2]],type='l',col='green',xlab='x',ylab='I',ylim=c(0,14))
abline(h=3,col='darkgoldenrod')
#points(c(0.342,0.369),c(impexp2(0.342,3.5,c(0.1,0.2,0.3,0.342,0.369,0.4,0.5)),impexp2(0.369,3.5,c(0.1,0.2,0.3,0.342,0.369,0.4,0.5))),col='purple',pch=19)

##########################################################################

layout(matrix(c(1,3,2,4), 2, 2, byrow = TRUE))
par(mar=c(4,4,2,2))

run1 = emufunexp(3.5,1.5,0.14)
plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14",cex.main=0.65)
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)
abline(h=3.5,col='black')
abline(h=(3.5+2*sqrt(0.025)),col='black',lty=2)
abline(h=(3.5-2*sqrt(0.025)),col='black',lty=2)

plot(sequem,sequem2,type='l',col='green',xlab='x',ylab='I')
locator()
abline(h=3,col='darkgoldenrod')
points(0.3664,impexp(0.3664,3.5),col='purple',pch=19)
#axis(1,at=0.317,cex.axis=0.5,tick=FALSE,line=-5)
#axis(1,at=0.396,cex.axis=0.5,tick=FALSE,line=-5)

newrun = emufunexp2(3.5,1.5,0.14,c(0.1,0.2,0.3,0.3664,0.4,0.5))
plot(newrun[[1]],newrun[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14",cex.main=0.65)
lines(newrun[[1]],newrun[[5]],type='l',lwd=1,col='red')
lines(newrun[[1]],newrun[[6]],type='l',lwd=1,col='red')
points(newrun[[3]],newrun[[4]],col='purple',pch=19)
abline(h=3.5,col='black')
abline(h=(3.5+2*sqrt(0.025)),col='black',lty=2)
abline(h=(3.5-2*sqrt(0.025)),col='black',lty=2)

newpointsrun1 = overallsequem(c(0.1,0.2,0.3,0.3664,0.4,0.5))  
plot(newpointsrun1[[1]],newpointsrun1[[2]],type='l',col='green',xlab='x',ylab='I',ylim=c(0,14))
abline(h=3,col='darkgoldenrod')

