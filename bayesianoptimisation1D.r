newfun = function(x){
  return(3*x*sin((5*pi*(x-0.1))/0.4))
}

emufunnew8p = function(B0,sigmau,theta,points){
  a = points
  D = newfun(a)
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length(a))
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

newfun2 = function(x){
  return(-newfun(x))
}

optim(0.4,newfun2,method="L-BFGS-B",lower=0.05,upper=0.55)
realxsin = 3*sequemsinexp*sin((5*pi*(sequemsinexp - 0.1))/0.4)


newrun8p = emufunnew8p(0,0.6,0.06,seq(0.1,0.51,length.out = 8))
par(mfrow=c(1,1))
plot(newrun8p[[1]],newrun8p[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",ylim=c(-1.5,1.5))
points(newrun8p[[3]],newrun8p[[4]],col='purple',pch=19)
lines(newrun8p[[1]],newrun8p[[5]],type='l',col='red')
lines(newrun8p[[1]],newrun8p[[6]],type='l',col='red')
lines(sequemsinexp,realxsin,type='l',col='black')

expectedimprovement = function(x,B0,sigmau,theta,points){
  a = points
  D = newfun(a)
  curmax = max(D)
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length(points))
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
  mu = eDf(x)
  var = varDf(x)
  til = 0.0
  if(var==0){
    accfun=0
  }
  else{
    z = (mu-curmax-til)/sqrt(var)
    accfun = (mu-curmax-til)*pnorm(z) + sqrt(var)*dnorm(z)
  }
  return(as.vector(-accfun))
}
#expectedimprovement(0.45,0,0.6,0.06,seq(0.1,0.5,length.out = 8))

sequem1d = seq(0.05,0.55,length.out=10000)
sequemnew1d = c()
for (i in sequem1d){
  sequemnew1d = append(sequemnew1d,expectedimprovement(i,B0=0,sigmau=0.6,theta=0.06,points=seq(0.1,0.51,length.out = 8)))
}
#sequemnew1d
plot(sequem1d,-sequemnew1d,type='l',col='darkgoldenrod',xlab="x",ylab="EI(x)", main = expression("Expected Improvement"))
points(newrun8p[[3]],rep(0,8),col='purple',pch=19)
points(as.vector(possnewax1D[[1]]),-as.vector(possnewax1D[[2]]),col='purple',pch=4,lwd=2)
axis(1,at=c(0.46))

startpoint1D = seq(0.1,0.51,length.out = 8)[which(newrun8p[[4]]==max(newrun8p[[4]]))]
startpoint1D
possnewax1D = optim(startpoint1D,expectedimprovement,B0=0,sigmau=0.6,theta=0.06,points=seq(0.1,0.51,length.out = 8),method="L-BFGS-B",lower=0.05,upper=0.55)
possnewax1D

newmaxpoints1D = c(seq(0.1,0.51,length.out = 8),as.vector(possnewax1D[[1]]))
newmaxpoints1D
newrun8p.1 = emufunnew8p(0,0.6,0.06,newmaxpoints1D)
plot(newrun8p.1[[1]],newrun8p.1[[2]],type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)")
points(newrun8p.1[[3]],newrun8p.1[[4]],col='purple',pch=19)
lines(newrun8p.1[[1]],newrun8p.1[[5]],type='l',col='red')
lines(newrun8p.1[[1]],newrun8p.1[[6]],type='l',col='red')
lines(sequemsinexp,realxsin,type='l',col='black')

dist(c(newfun(0.4579478),1.382109))
dist(c(newfun(0.4514286),1.382109))
#sequemnew1d.1 = c()
#for (i in sequem1d){
  #sequemnew1d.1 = append(sequemnew1d.1,expectedimprovement(i,B0=0,sigmau=0.6,theta=0.06,points=newmaxpoints))
#}
#plot(sequem1d,-sequemnew1d.1,type='l',col='orange',xlab="x",ylab="EI(x)")
#startpoint.1 = newmaxpoints[which(newrun8p.1[[4]]==max(newrun8p.1[[4]]))]
#startpoint.1
#possnewax.1 = optim(startpoint.1,expectedimprovement,B0=0,sigmau=0.6,theta=0.06,points=newmaxpoints,method="L-BFGS-B",lower=0.05,upper=0.55)
#possnewax
