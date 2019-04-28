pdf("Emulatorplotssidebyside.pdf")

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

run1 = emufunexp(3.5,1.5,0.14)
plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run3.0 = emufunexp(-3.5,1.5,0.14)
plot(run3.0[[1]],run3.0[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==-3.5,", ",sigma[u]==1.5,", ",theta==0.14)),
     cex.main=0.75)
lines(run3.0[[1]],run3.0[[5]],type='l',lwd=1,col='red')
lines(run3.0[[1]],run3.0[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run3.0[[3]],run3.0[[4]],col='purple',pch=19)

run5.0 = emufunexp(10.5,1.5,0.14)
plot(run5.0[[1]],run5.0[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="x",
     ylab="f(x)",
     main=expression(paste("Exponential Model Emulator with ",beta[0]==10.5,", ",sigma[u]==1.5,", ",theta==0.14)),
     cex.main=0.75)
lines(run5.0[[1]],run5.0[[5]],type='l',lwd=1,col='red')
lines(run5.0[[1]],run5.0[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run5.0[[3]],run5.0[[4]],col='purple',pch=19)

######################################################################

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)),
     cex.main=0.75)
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
lines(sequemsinexp,realexpthreefive,type='l',lwd=1,col='black')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.1 = emufunexp(3.5,0.8,0.14)
plot(run2.1[[1]],run2.1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==0.8,", ",theta==0.14)),
     cex.main=0.75)
lines(run2.1[[1]],run2.1[[5]],type='l',lwd=1,col='red')
lines(run2.1[[1]],run2.1[[6]],type='l',lwd=1,col='red')
points(run2.1[[3]],run2.1[[4]],col='purple',pch=19)

run3.1 = emufunexp(3.5,0.1,0.14)
plot(run3.1[[1]],run3.1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==0.1,", ",theta==0.14)),
     cex.main=0.75)
lines(run3.1[[1]],run3.1[[5]],type='l',lwd=1,col='red')
lines(run3.1[[1]],run3.1[[6]],type='l',lwd=1,col='red')
points(run3.1[[3]],run3.1[[4]],col='purple',pch=19)

run5.1 = emufunexp(3.5,2.9,0.14)
plot(run5.1[[1]],run5.1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==2.9,", ",theta==0.14)),
     cex.main=0.75)
lines(run5.1[[1]],run5.1[[5]],type='l',lwd=1,col='red')
lines(run5.1[[1]],run5.1[[6]],type='l',lwd=1,col='red')
points(run5.1[[3]],run5.1[[4]],col='purple',pch=19)

########################################################################

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)),
     cex.main=0.75)
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.2 = emufunexp(3.5,1.5,0.075)
plot(run2.2[[1]],run2.2[[2]],type='l',lwd=1,col='blue',
     ylim=c(-0.6,7.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.075)),
     cex.main=0.75)
lines(run2.2[[1]],run2.2[[5]],type='l',lwd=1,col='red')
lines(run2.2[[1]],run2.2[[6]],type='l',lwd=1,col='red')
points(run2.2[[3]],run2.2[[4]],col='purple',pch=19)

run3.2 = emufunexp(3.5,1.5,0.01)
plot(run3.2[[1]],run3.2[[2]],type='l',lwd=1,col='blue',
     ylim=c(-1.5,8.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.01)),
     cex.main=0.75)
lines(run3.2[[1]],run3.2[[5]],type='l',lwd=1,col='red')
lines(run3.2[[1]],run3.2[[6]],type='l',lwd=1,col='red')
points(run3.2[[3]],run3.2[[4]],col='purple',pch=19)

run5.2 = emufunexp(3.5,1.5,0.27)
plot(run5.2[[1]],run5.2[[2]],type='l',lwd=1,col='blue',
     ylim=c(1,6),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.27)),
     cex.main=0.75)
lines(run5.2[[1]],run5.2[[5]],type='l',lwd=1,col='red')
lines(run5.2[[1]],run5.2[[6]],type='l',lwd=1,col='red')
points(run5.2[[3]],run5.2[[4]],col='purple',pch=19)

######################################################################

layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.3 = emufunexp(-3.5,0.1,0.14)
plot(run2.3[[1]],run2.3[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==-3.5,", ",sigma[u]==0.1,", ",theta==0.14)))
lines(run2.3[[1]],run2.3[[5]],type='l',lwd=1,col='red')
lines(run2.3[[1]],run2.3[[6]],type='l',lwd=1,col='red')
points(run2.3[[3]],run2.3[[4]],col='purple',pch=19)

#######################################################################

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run3.4 = emufunexp(10.5,1.5,0.27)
plot(run3.4[[1]],run3.4[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==10.5,", ",sigma[u]==1.5,", ",theta==0.27)),
     cex.main=0.75)
lines(run3.4[[1]],run3.4[[5]],type='l',lwd=1,col='red')
lines(run3.4[[1]],run3.4[[6]],type='l',lwd=1,col='red')
points(run3.4[[3]],run3.4[[4]],col='purple',pch=19)

run5.4 = emufunexp(10.5,1.5,0.01)
plot(run5.4[[1]],run5.4[[2]],type='l',lwd=1,col='blue',
     ylim=c(0,15),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==10.5,", ",sigma[u]==1.5,", ",theta==0.01)),
     cex.main=0.75)
lines(run5.4[[1]],run5.4[[5]],type='l',lwd=1,col='red')
lines(run5.4[[1]],run5.4[[6]],type='l',lwd=1,col='red')
points(run5.4[[3]],run5.4[[4]],col='purple',pch=19)

#####################################################################

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.5 = emufunexp(3.5,0.1,0.01)
plot(run2.5[[1]],run2.5[[2]],type='l',lwd=1,col='blue',
     ylim=c(1,6),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==0.1,", ",theta==0.01)),
     cex.main=0.75)
lines(run2.5[[1]],run2.5[[5]],type='l',lwd=1,col='red')
lines(run2.5[[1]],run2.5[[6]],type='l',lwd=1,col='red')
points(run2.5[[3]],run2.5[[4]],col='purple',pch=19)

run3.5 = emufunexp(3.5,2.9,0.27)
plot(run3.5[[1]],run3.5[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==2.9,", ",theta==0.27)),
     cex.main=0.75)
lines(run3.5[[1]],run3.5[[5]],type='l',lwd=1,col='red')
lines(run3.5[[1]],run3.5[[6]],type='l',lwd=1,col='red')
points(run3.5[[3]],run3.5[[4]],col='purple',pch=19)

######################################################################

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,3,2))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',
     ylim=c(0.5,6.5),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==3.5,", ",sigma[u]==1.5,", ",theta==0.14)))
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.6 = emufunexp(-3.5,0.1,0.01)
plot(run2.6[[1]],run2.6[[2]],type='l',lwd=1,col='blue',
     ylim=c(-4,6),xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==-3.5,", ",sigma[u]==0.1,", ",theta==0.01)),
     cex.main=0.75)
lines(run2.6[[1]],run2.6[[5]],type='l',lwd=1,col='red')
lines(run2.6[[1]],run2.6[[6]],type='l',lwd=1,col='red')
points(run2.6[[3]],run2.6[[4]],col='purple',pch=19)

run3.6 = emufunexp(10.5,2.9,0.27)
plot(run3.6[[1]],run3.6[[2]],type='l',lwd=1,col='blue',
     xlab="Rate parameter value x",
     ylab="Concentration of f(x)",
     main=expression(paste("Exponential Model Emulation with ",beta[0]==10.5,", ",sigma[u]==2.9,", ",theta==0.27)),
     cex.main=0.75)
lines(run3.6[[1]],run3.6[[5]],type='l',lwd=1,col='red')
lines(run3.6[[1]],run3.6[[6]],type='l',lwd=1,col='red')
points(run3.6[[3]],run3.6[[4]],col='purple',pch=19)

dev.off()
