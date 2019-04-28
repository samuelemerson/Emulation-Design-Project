expfun = function(x){
  return(exp(x*3.5))
}

covfun = function(x,y,su,t){
  entr = su*exp(-(abs(x-y)^2)/(t^2))
  return(entr)
}

emufunexp = function(B0,sigmau,theta){
  a = seq(0.1,0.5,by=0.1)
  D = expfun(a)
  exf = B0
  varf = sigmau^2
  eD = c(exf,exf,exf,exf,exf)
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

###################################################################################

#par(mfrow=c(2,3))
pdf("Emulatorplots.pdf")

run1 = emufunexp(3.5,1.5,0.14)
plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

#run2.0 = emufunexp(0.0,1.5,0.14)
#plot(run2.0[[1]],run2.0[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Bo=0.0")
#points(run2.0[[3]],run2.0[[4]],col='purple',pch=19)
#lines(run2.0[[1]],run2.0[[5]],type='l',col='red')
#lines(run2.0[[1]],run2.0[[6]],type='l',col='red')

run3.0 = emufunexp(-3.5,1.5,0.14)
plot(run3.0[[1]],run3.0[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulation with Beta0=-3.5,Sigmau=1.5,Theta=0.14")
lines(run3.0[[1]],run3.0[[5]],type='l',lwd=1,col='red')
lines(run3.0[[1]],run3.0[[6]],type='l',lwd=1,col='red')
points(run3.0[[3]],run3.0[[4]],col='purple',pch=19)

#run4.0 = emufunexp(7.0,1.5,0.14)
#plot(run4.0[[1]],run4.0[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=7.0")
#points(run4.0[[3]],run4.0[[4]],col='purple',pch=19)
#lines(run4.0[[1]],run4.0[[5]],type='l',col='red')
#lines(run4.0[[1]],run4.0[[6]],type='l',col='red')

run5.0 = emufunexp(10.5,1.5,0.14)
plot(run5.0[[1]],run5.0[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulation when Beta0=10.5,Sigmau=1.5,Theta=0.14")
lines(run5.0[[1]],run5.0[[5]],type='l',lwd=1,col='red')
lines(run5.0[[1]],run5.0[[6]],type='l',lwd=1,col='red')
points(run5.0[[3]],run5.0[[4]],col='purple',pch=19)

####################################################################################

#par(mfrow=c(2,3))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.1 = emufunexp(3.5,0.8,0.14)
plot(run2.1[[1]],run2.1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=0.8,Theta=0.14")
lines(run2.1[[1]],run2.1[[5]],type='l',lwd=1,col='red')
lines(run2.1[[1]],run2.1[[6]],type='l',lwd=1,col='red')
points(run2.1[[3]],run2.1[[4]],col='purple',pch=19)

run3.1 = emufunexp(3.5,0.1,0.14)
plot(run3.1[[1]],run3.1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=0.1,Theta=0.14")
lines(run3.1[[1]],run3.1[[5]],type='l',lwd=1,col='red')
lines(run3.1[[1]],run3.1[[6]],type='l',lwd=1,col='red')
points(run3.1[[3]],run3.1[[4]],col='purple',pch=19)

#run4.1 = emufunexp(3.5,2.2,0.14)
#plot(run4.1[[1]],run4.1[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Sigmau=2.2")
#points(run4.1[[3]],run4.1[[4]],col='purple',pch=19)
#lines(run4.1[[1]],run4.1[[5]],type='l',col='red')
#lines(run4.1[[1]],run4.1[[6]],type='l',col='red')

run5.1 = emufunexp(3.5,2.9,0.14)
plot(run5.1[[1]],run5.1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=2.9,Theta=0.14")
lines(run5.1[[1]],run5.1[[5]],type='l',lwd=1,col='red')
lines(run5.1[[1]],run5.1[[6]],type='l',lwd=1,col='red')
points(run5.1[[3]],run5.1[[4]],col='purple',pch=19)

###################################################################################

#par(mfrow=c(2,3))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.2 = emufunexp(3.5,1.5,0.075)
plot(run2.2[[1]],run2.2[[2]],type='l',lwd=1,col='blue',ylim=c(-0.6,7.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.075")
lines(run2.2[[1]],run2.2[[5]],type='l',lwd=1,col='red')
lines(run2.2[[1]],run2.2[[6]],type='l',lwd=1,col='red')
points(run2.2[[3]],run2.2[[4]],col='purple',pch=19)

run3.2 = emufunexp(3.5,1.5,0.01)
plot(run3.2[[1]],run3.2[[2]],type='l',lwd=1,col='blue',ylim=c(-1.5,8.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.01")
lines(run3.2[[1]],run3.2[[5]],type='l',lwd=1,col='red')
lines(run3.2[[1]],run3.2[[6]],type='l',lwd=1,col='red')
points(run3.2[[3]],run3.2[[4]],col='purple',pch=19)

#run4.2 = emufunexp(3.5,1.5,0.205)
#plot(run4.2[[1]],run4.2[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Theta=0.205")
#points(run4.2[[3]],run4.2[[4]],col='purple',pch=19)
#lines(run4.2[[1]],run4.2[[5]],type='l',col='red')
#lines(run4.2[[1]],run4.2[[6]],type='l',col='red')

run5.2 = emufunexp(3.5,1.5,0.27)
plot(run5.2[[1]],run5.2[[2]],type='l',lwd=1,col='blue',ylim=c(1,6),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.27")
lines(run5.2[[1]],run5.2[[5]],type='l',lwd=1,col='red')
lines(run5.2[[1]],run5.2[[6]],type='l',lwd=1,col='red')
points(run5.2[[3]],run5.2[[4]],col='purple',pch=19)

###########################################################

#par(mfrow=c(2,3))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.3 = emufunexp(-3.5,0.1,0.14)
plot(run2.3[[1]],run2.3[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=-3.5,Sigmau=0.1,Theta=0.14")
lines(run2.3[[1]],run2.3[[5]],type='l',lwd=1,col='red')
lines(run2.3[[1]],run2.3[[6]],type='l',lwd=1,col='red')
points(run2.3[[3]],run2.3[[4]],col='purple',pch=19)

#run3.3 = emufunexp(10.5,2.9,0.14)
#plot(run3.3[[1]],run3.3[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=10.5,Sigmau=2.9")
#points(run3.3[[3]],run3.3[[4]],col='purple',pch=19)
#lines(run3.3[[1]],run3.3[[5]],type='l',col='red')
#lines(run3.3[[1]],run3.3[[6]],type='l',col='red')

#run4.3 = emufunexp(-3.5,2.9,0.14)
#plot(run4.3[[1]],run4.3[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=-3.5,Sigmau=2.9")
#points(run4.3[[3]],run4.3[[4]],col='purple',pch=19)
#lines(run4.3[[1]],run4.3[[5]],type='l',col='red')
#lines(run4.3[[1]],run4.3[[6]],type='l',col='red')

#run5.3 = emufunexp(10.5,0.1,0.14)
#plot(run5.3[[1]],run5.3[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=10.5,Sigmau=0.1")
#points(run5.3[[3]],run5.3[[4]],col='purple',pch=19)
#lines(run5.3[[1]],run5.3[[5]],type='l',col='red')
#lines(run5.3[[1]],run5.3[[6]],type='l',col='red')

#############################################################

#par(mfrow=c(2,3))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

#run2.4 = emufunexp(-3.5,1.5,0.01)
#plot(run2.4[[1]],run2.4[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=-3.5,Theta=0.01")
#points(run2.4[[3]],run2.4[[4]],col='purple',pch=19)
#lines(run2.4[[1]],run2.4[[5]],type='l',col='red')
#lines(run2.4[[1]],run2.4[[6]],type='l',col='red')

run3.4 = emufunexp(10.5,1.5,0.27)
plot(run3.4[[1]],run3.4[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=10.5,Sigmau=1.5,Theta=0.27")
lines(run3.4[[1]],run3.4[[5]],type='l',lwd=1,col='red')
lines(run3.4[[1]],run3.4[[6]],type='l',lwd=1,col='red')
points(run3.4[[3]],run3.4[[4]],col='purple',pch=19)

#run4.4 = emufunexp(-3.5,1.5,0.27)
#plot(run4.4[[1]],run4.4[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=-3.5,Theta=0.27")
#points(run4.4[[3]],run4.4[[4]],col='purple',pch=19)
#lines(run4.4[[1]],run4.4[[5]],type='l',col='red')
#lines(run4.4[[1]],run4.4[[6]],type='l',col='red')

run5.4 = emufunexp(10.5,1.5,0.01)
plot(run5.4[[1]],run5.4[[2]],type='l',lwd=1,col='blue',ylim=c(0,15),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=10.5,Sigmau=1.5,Theta=0.01")
lines(run5.4[[1]],run5.4[[5]],type='l',lwd=1,col='red')
lines(run5.4[[1]],run5.4[[6]],type='l',lwd=1,col='red')
points(run5.4[[3]],run5.4[[4]],col='purple',pch=19)

###########################################################

#par(mfrow=c(2,3))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.5 = emufunexp(3.5,0.1,0.01)
plot(run2.5[[1]],run2.5[[2]],type='l',lwd=1,col='blue',ylim=c(1,6),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=0.1,Theta=0.01")
lines(run2.5[[1]],run2.5[[5]],type='l',lwd=1,col='red')
lines(run2.5[[1]],run2.5[[6]],type='l',lwd=1,col='red')
points(run2.5[[3]],run2.5[[4]],col='purple',pch=19)

run3.5 = emufunexp(3.5,2.9,0.27)
plot(run3.5[[1]],run3.5[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=2.9,Theta=0.27")
lines(run3.5[[1]],run3.5[[5]],type='l',lwd=1,col='red')
lines(run3.5[[1]],run3.5[[6]],type='l',lwd=1,col='red')
points(run3.5[[3]],run3.5[[4]],col='purple',pch=19)

#run4.5 = emufunexp(3.5,0.1,0.27)
#plot(run4.5[[1]],run4.5[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Sigmau=0.1,Theta=0.27")
#points(run4.5[[3]],run4.5[[4]],col='purple',pch=19)
#lines(run4.5[[1]],run4.5[[5]],type='l',col='red')
#lines(run4.5[[1]],run4.5[[6]],type='l',col='red')

#run5.5 = emufunexp(3.5,2.9,0.01)
#plot(run5.5[[1]],run5.5[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Sigmau=2.9,Theta=0.01")
#points(run5.5[[3]],run5.5[[4]],col='purple',pch=19)
#lines(run5.5[[1]],run5.5[[5]],type='l',col='red')
#lines(run5.5[[1]],run5.5[[6]],type='l',col='red')

###########################################################

#par(mfrow=c(2,3))

plot(run1[[1]],run1[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(run1[[1]],run1[[5]],type='l',lwd=1,col='red')
lines(run1[[1]],run1[[6]],type='l',lwd=1,col='red')
points(run1[[3]],run1[[4]],col='purple',pch=19)

run2.6 = emufunexp(-3.5,0.1,0.01)
plot(run2.6[[1]],run2.6[[2]],type='l',lwd=1,col='blue',ylim=c(-4,6),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=-3.5,Sigmau=0.1,Theta=0.01")
lines(run2.6[[1]],run2.6[[5]],type='l',lwd=1,col='red')
lines(run2.6[[1]],run2.6[[6]],type='l',lwd=1,col='red')
points(run2.6[[3]],run2.6[[4]],col='purple',pch=19)

run3.6 = emufunexp(10.5,2.9,0.27)
plot(run3.6[[1]],run3.6[[2]],type='l',lwd=1,col='blue',xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=10.5,Sigmau=2.9,Theta=0.27")
lines(run3.6[[1]],run3.6[[5]],type='l',lwd=1,col='red')
lines(run3.6[[1]],run3.6[[6]],type='l',lwd=1,col='red')
points(run3.6[[3]],run3.6[[4]],col='purple',pch=19)

#run4.6 = emufunexp(-3.5,2.9,0.01)
#plot(run4.6[[1]],run4.6[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=-3.5,Sigmau=2.9,Theta=0.01")
#points(run4.6[[3]],run4.6[[4]],col='purple',pch=19)
#lines(run4.6[[1]],run4.6[[5]],type='l',col='red')
#lines(run4.6[[1]],run4.6[[6]],type='l',col='red')

#run5.6 = emufunexp(10.5,0.1,0.27)
#plot(run5.6[[1]],run5.6[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="B0=10.5,Sigmau=0.1,Theta=0.27")
#points(run5.6[[3]],run5.6[[4]],col='purple',pch=19)
#lines(run5.6[[1]],run5.6[[5]],type='l',col='red')
#lines(run5.6[[1]],run5.6[[6]],type='l',col='red')

dev.off()
####################################################################

newfun = function(x){
  return(3*x*sin((5*pi*(x-0.1))/0.4))
}

emufunnew = function(B0,sigmau,theta){
  a = seq(0.1,0.5,length.out = 10)
  D = newfun(a)
  exf = B0
  varf = sigmau^2
  eD = c(exf,exf,exf,exf,exf)
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

newrun1 = emufunnew(0,0.6,0.06)
par(mfrow=c(1,1))
plot(newrun1[[1]],newrun1[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)")
points(newrun1[[3]],newrun1[[4]],col='purple',pch=19)
lines(newrun1[[1]],newrun1[[5]],type='l',col='red')
lines(newrun1[[1]],newrun1[[6]],type='l',col='red')


#############################################################

emufuncos = function(B0,sigmau,theta){
  a = seq(0,2*pi,by = (pi/4))
  D = cos(a)
  exf = B0
  varf = sigmau^2
  eD = c(exf,exf,exf,exf,exf,exf,exf,exf,exf)
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
  sequem = seq(-0.5,(2*pi)+0.5,length.out=10000)
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

cosrun1 = emufuncos(0,0.25,0.65)
plot(cosrun1[[1]],cosrun1[[2]],type='l',col='blue',xlim = c(-0.5,(2*pi)+0.5),ylim = c(-1.3,1.4),xlab="Rate parameter value x",ylab="Concentartion of f(x)")
points(cosrun1[[3]],cosrun1[[4]],col='purple',pch=19)
lines(cosrun1[[1]],cosrun1[[5]],type='l',col='red')
lines(cosrun1[[1]],cosrun1[[6]],type='l',col='red')
#e = 0.05
############################################################

impexp=function(x,z){
  a = seq(0.1,0.5,by=0.1)
  D = expfun(a)
  exf = 3.5
  varf = 1.5^2
  varep = 0
  vare = 0.025
  eD = c(exf,exf,exf,exf,exf)
  b = as.matrix(dist(a))
  varD = (1.5^2)*exp(-(b^2)/(0.14^2))
  covfD = function(x){
    g = covfun(x,a,varf,0.14)
    return(g)
  }
  covDf = function(x){
    i = covfun(a,x,varf,0.14)
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
  Isquared = ((eDf(x)-z)^2)/(varDf(x)+varep+vare)
  I = sqrt(Isquared)
  return(I)
}

sequem = seq(0.05,0.55,length.out=10000)
sequem2 = c()
for (i in sequem){
  sequem2 = append(sequem2,impexp(i,3.5))
}

plot(sequem,sequem2,type='l',col='darkorange',xlab='x',ylab='I')
abline(h=3,col='black')
locator()

#####################################################################

emufunexp2 = function(B0,sigmau,theta,points){
  D = expfun(points)
  exf = B0
  varf = sigmau^2
  eD = rep(exf,times=length(D))
  b = as.matrix(dist(points))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  covfD = function(x){
    g = covfun(x,points,varf,theta)
    return(g)
  }
  covDf = function(x){
    i = covfun(points,x,varf,theta)
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
  return(list(sequem,sequem2,points,D,upvar,downvar))
}

newrun = emufunexp2(3.5,1.5,0.14,c(0.1,0.2,0.3,0.321,0.395,0.4,0.5))
plot(newrun[[1]],newrun[[2]],type='l',lwd=1,col='blue',ylim=c(0.5,6.5),xlab="Rate parameter value x",ylab="Concentration of f(x)",main="Exponential Model Emulator with Beta0=3.5,Sigmau=1.5,Theta=0.14")
lines(newrun[[1]],newrun[[5]],type='l',lwd=1,col='red')
lines(newrun[[1]],newrun[[6]],type='l',lwd=1,col='red')
points(newrun[[3]],newrun[[4]],col='purple',pch=19)

impexp2=function(x,z,points){
  D = expfun(points)
  exf = 3.5
  varf = 1.5^2
  varep = 0
  vare = 0.025
  eD = rep(exf,times=length(D))
  b = as.matrix(dist(points))
  varD = (1.5^2)*exp(-(b^2)/(0.14^2))
  covfD = function(x){
    g = covfun(x,points,varf,0.14)
    return(g)
  }
  covDf = function(x){
    i = covfun(points,x,varf,0.14)
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
  Isquared = ((eDf(x)-z)^2)/(varDf(x)+varep+vare)
  I = sqrt(Isquared)
  return(I)
}

overallsequem = function(points){
  sequem = seq(0.05,0.55,length.out=10000)
  newsequem = c()
  for (i in sequem){
    newsequem = append(newsequem,impexp2(i,3.5,points))
  }
  return(list(sequem,newsequem))
}

newpointsrun1 = overallsequem(c(0.1,0.2,0.3,0.317,0.396,0.4,0.5))  
plot(newpointsrun1[[1]],newpointsrun1[[2]],type='l',col='green',xlab='x',ylab='I')
abline(h=2,col='yellow')
locator()

#fill.contour
#add points
#change code to 2d
#test 2d function sin+cos for eg
#look up expand grid function....probably for sequem (might need as.matrix)