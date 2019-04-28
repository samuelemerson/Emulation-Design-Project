expfun = function(x){
  return(exp(x*3.5))
}
D = c(expfun(0.1),expfun(0.2),expfun(0.3),expfun(0.4),expfun(0.5))

exf = 3.5

varf = 1.5^2

eD = c(exf,exf,exf,exf,exf)

covfun = function(x,y,su,t){
  entr = su*exp(-(abs(x-y)^2)/(t^2))
  return(entr)
}
a = seq(0.1,0.5,by=0.1)
b = covfun(0.1,a,varf,0.14)
c = covfun(0.2,a,varf,0.14)
d = covfun(0.3,a,varf,0.14)
e = covfun(0.4,a,varf,0.14)
f = covfun(0.5,a,varf,0.14)

varD = matrix(c(b,c,d,e,f),nrow=5,ncol=5,byrow=TRUE)

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

par(mfrow=c(2,3))
plot(sequem,sequem2,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=3.5")
points(a,D,col='purple',pch=19)
lines(sequem,upvar,type='l',col='red')
lines(sequem,downvar,type='l',col='red')


exfnew = 0.0

eDnew = c(exfnew,exfnew,exfnew,exfnew,exfnew)

eDfnew = function(x){
  h = exfnew + covfD(x)%*%solve(varD)%*%(D-eDnew)
  return(as.vector(h))
}

sequem2new = c()
for (i in sequem){
  sequem2new = append(sequem2new,eDfnew(i))
}

upvarnew = sequem2new + 3*sqrt(sequem3)
downvarnew = sequem2new - 3*sqrt(sequem3)

plot(sequem,sequem2new,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=0.0")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnew,type='l',col='red')
lines(sequem,downvarnew,type='l',col='red')

exfnew2 = -3.5

eDnew2 = c(exfnew2,exfnew2,exfnew2,exfnew2,exfnew2)

eDfnew2 = function(x){
  h = exfnew2 + covfD(x)%*%solve(varD)%*%(D-eDnew2)
  return(as.vector(h))
}

sequem2new2 = c()
for (i in sequem){
  sequem2new2 = append(sequem2new2,eDfnew2(i))
}

upvarnew2 = sequem2new2 + 3*sqrt(sequem3)
downvarnew2 = sequem2new2 - 3*sqrt(sequem3)

plot(sequem,sequem2new2,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=-3.5")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnew2,type='l',col='red')
lines(sequem,downvarnew2,type='l',col='red')

exfnew3 = 7.0

eDnew3 = c(exfnew3,exfnew3,exfnew3,exfnew3,exfnew3)

eDfnew3 = function(x){
  h = exfnew3 + covfD(x)%*%solve(varD)%*%(D-eDnew3)
  return(as.vector(h))
}

sequem2new3 = c()
for (i in sequem){
  sequem2new3 = append(sequem2new3,eDfnew3(i))
}

upvarnew3 = sequem2new3 + 3*sqrt(sequem3)
downvarnew3 = sequem2new3 - 3*sqrt(sequem3)

plot(sequem,sequem2new3,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=7.0")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnew3,type='l',col='red')
lines(sequem,downvarnew3,type='l',col='red')

exfnew4 = 10.5

eDnew4 = c(exfnew4,exfnew4,exfnew4,exfnew4,exfnew4)

eDfnew4 = function(x){
  h = exfnew4 + covfD(x)%*%solve(varD)%*%(D-eDnew4)
  return(as.vector(h))
}

sequem2new4 = c()
for (i in sequem){
  sequem2new4 = append(sequem2new4,eDfnew4(i))
}

upvarnew4 = sequem2new4 + 3*sqrt(sequem3)
downvarnew4 = sequem2new4 - 3*sqrt(sequem3)

plot(sequem,sequem2new4,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=0.0")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnew4,type='l',col='red')
lines(sequem,downvarnew4,type='l',col='red')

#########################################################

par(mfrow=c(2,3))
plot(sequem,sequem2,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="sigmau=1.5")
points(a,D,col='purple',pch=19)
lines(sequem,upvar,type='l',col='red')
lines(sequem,downvar,type='l',col='red')

varfnew = 0.8^2

bnew = covfun(0.1,a,varfnew,0.14)
cnew = covfun(0.2,a,varfnew,0.14)
dnew = covfun(0.3,a,varfnew,0.14)
enew = covfun(0.4,a,varfnew,0.14)
fnew = covfun(0.5,a,varfnew,0.14)

varDnew = matrix(c(bnew,cnew,dnew,enew,fnew),nrow=5,ncol=5,byrow=TRUE)

covfDnew = function(x){
  g = covfun(x,a,varfnew,0.14)
  return(g)
}
covDfnew = function(x){
  i = covfun(a,x,varfnew,0.14)
  return(i)
}

eDfnewnew = function(x){
  h = exf + covfDnew(x)%*%solve(varDnew)%*%(D-eD)
  return(as.vector(h))
}

varDfnew = function(x){
  j = varfnew - covfDnew(x)%*%solve(varDnew)%*%covDfnew(x)
  return(j)
}


sequem2newnew = c()
for (i in sequem){
  sequem2newnew = append(sequem2newnew,eDfnewnew(i))
}

sequem3new = c()
for (i in sequem){
  sequem3new = append(sequem3new,varDfnew(i))
}

upvarnewnew = sequem2newnew + 3*sqrt(sequem3new)
downvarnewnew = sequem2newnew - 3*sqrt(sequem3new)

plot(sequem,sequem2newnew,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="sigmau=0.8")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnewnew,type='l',col='red')
lines(sequem,downvarnewnew,type='l',col='red')

varfnew2 = 0.1^2

bnew2 = covfun(0.1,a,varfnew2,0.14)
cnew2 = covfun(0.2,a,varfnew2,0.14)
dnew2 = covfun(0.3,a,varfnew2,0.14)
enew2 = covfun(0.4,a,varfnew2,0.14)
fnew2 = covfun(0.5,a,varfnew2,0.14)

varDnew2 = matrix(c(bnew2,cnew2,dnew2,enew2,fnew2),nrow=5,ncol=5,byrow=TRUE)

covfDnew2 = function(x){
  g = covfun(x,a,varfnew2,0.14)
  return(g)
}
covDfnew2 = function(x){
  i = covfun(a,x,varfnew2,0.14)
  return(i)
}

eDfnewnew2 = function(x){
  h = exf + covfDnew2(x)%*%solve(varDnew2)%*%(D-eD)
  return(as.vector(h))
}

varDfnew2 = function(x){
  j = varfnew2 - covfDnew2(x)%*%solve(varDnew2)%*%covDfnew2(x)
  return(j)
}


sequem2newnew2 = c()
for (i in sequem){
  sequem2newnew2 = append(sequem2newnew2,eDfnewnew2(i))
}

sequem3new2 = c()
for (i in sequem){
  sequem3new2 = append(sequem3new2,varDfnew2(i))
}

upvarnewnew2 = sequem2newnew2 + 3*sqrt(sequem3new2)
downvarnewnew2 = sequem2newnew2 - 3*sqrt(sequem3new2)

plot(sequem,sequem2newnew2,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="sigmau=0.1")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnewnew2,type='l',col='red')
lines(sequem,downvarnewnew2,type='l',col='red')

varfnew3 = 2.2^2

bnew3 = covfun(0.1,a,varfnew3,0.14)
cnew3 = covfun(0.2,a,varfnew3,0.14)
dnew3 = covfun(0.3,a,varfnew3,0.14)
enew3 = covfun(0.4,a,varfnew3,0.14)
fnew3 = covfun(0.5,a,varfnew3,0.14)

varDnew3 = matrix(c(bnew3,cnew3,dnew3,enew3,fnew3),nrow=5,ncol=5,byrow=TRUE)

covfDnew3 = function(x){
  g = covfun(x,a,varfnew3,0.14)
  return(g)
}
covDfnew3 = function(x){
  i = covfun(a,x,varfnew3,0.14)
  return(i)
}

eDfnewnew3 = function(x){
  h = exf + covfDnew3(x)%*%solve(varDnew3)%*%(D-eD)
  return(as.vector(h))
}

varDfnew3 = function(x){
  j = varfnew3 - covfDnew3(x)%*%solve(varDnew3)%*%covDfnew3(x)
  return(j)
}


sequem2newnew3 = c()
for (i in sequem){
  sequem2newnew3 = append(sequem2newnew3,eDfnewnew3(i))
}

sequem3new3 = c()
for (i in sequem){
  sequem3new3 = append(sequem3new3,varDfnew3(i))
}

upvarnewnew3 = sequem2newnew3 + 3*sqrt(sequem3new3)
downvarnewnew3 = sequem2newnew3 - 3*sqrt(sequem3new3)

plot(sequem,sequem2newnew3,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="sigmau=2.2")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnewnew3,type='l',col='red')
lines(sequem,downvarnewnew3,type='l',col='red')

varfnew4 = 2.9^2

bnew4 = covfun(0.1,a,varfnew4,0.14)
cnew4 = covfun(0.2,a,varfnew4,0.14)
dnew4 = covfun(0.3,a,varfnew4,0.14)
enew4 = covfun(0.4,a,varfnew4,0.14)
fnew4 = covfun(0.5,a,varfnew4,0.14)

varDnew4 = matrix(c(bnew4,cnew4,dnew4,enew4,fnew4),nrow=5,ncol=5,byrow=TRUE)

covfDnew4 = function(x){
  g = covfun(x,a,varfnew4,0.14)
  return(g)
}
covDfnew4 = function(x){
  i = covfun(a,x,varfnew4,0.14)
  return(i)
}

eDfnewnew4 = function(x){
  h = exf + covfDnew4(x)%*%solve(varDnew4)%*%(D-eD)
  return(as.vector(h))
}

varDfnew4 = function(x){
  j = varfnew4 - covfDnew4(x)%*%solve(varDnew4)%*%covDfnew4(x)
  return(j)
}


sequem2newnew4 = c()
for (i in sequem){
  sequem2newnew4 = append(sequem2newnew4,eDfnewnew4(i))
}

sequem3new4 = c()
for (i in sequem){
  sequem3new4 = append(sequem3new4,varDfnew4(i))
}

upvarnewnew4 = sequem2newnew4 + 3*sqrt(sequem3new4)
downvarnewnew4 = sequem2newnew4 - 3*sqrt(sequem3new4)

plot(sequem,sequem2newnew4,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="sigmau=2.9")
points(a,D,col='purple',pch=19)
lines(sequem,upvarnewnew4,type='l',col='red')
lines(sequem,downvarnewnew4,type='l',col='red')


#########################################################

par(mfrow=c(1,2))
plot(sequem,sequem2,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="correlation length=0.14")
points(a,D,col='purple',pch=19)
lines(sequem,upvar,type='l',col='red')
lines(sequem,downvar,type='l',col='red')

bnewnew = covfun(0.1,a,varf,0.07)
cnewnew = covfun(0.2,a,varf,0.07)
dnewnew = covfun(0.3,a,varf,0.07)
enewnew = covfun(0.4,a,varf,0.07)
fnewnew = covfun(0.5,a,varf,0.07)

varDnewnew = matrix(c(bnewnew,cnewnew,dnewnew,enewnew,fnewnew),nrow=5,ncol=5,byrow=TRUE)

covfDnewnew = function(x){
  g = covfun(x,a,varf,0.07)
  return(g)
}
covDfnewnew = function(x){
  i = covfun(a,x,varf,0.07)
  return(i)
}

eDfnewnewnew = function(x){
  h = exf + covfDnewnew(x)%*%solve(varDnewnew)%*%(D-eD)
  return(as.vector(h))
}

varDfnewnew = function(x){
  j = varf - covfDnewnew(x)%*%solve(varDnewnew)%*%covDfnewnew(x)
  return(j)
}


sequem2newnew2 = c()
for (i in sequem){
  sequem2 = append(sequem2,eDf(i))
}

sequem3 = c()
for (i in sequem){
  sequem3 = append(sequem3,varDf(i))
}

upvar = sequem2 + 3*sqrt(sequem3)
downvar = sequem2 - 3*sqrt(sequem3)

par(mfrow=c(2,3))
plot(sequem,sequem2,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)",main="Bo=3.5")
points(a,D,col='purple',pch=19)
lines(sequem,upvar,type='l',col='red')
lines(sequem,downvar,type='l',col='red')


##########################################################

newfun = function(x){
  return(3*x*sin((5*pi*(x-0.1))/0.4))
}
a2 = seq(0.1,0.5,length.out = 10)
D2 = newfun(even)
exf2 = 0
varf2 = 0.6^2
eD2 = c(exf2,exf2,exf2,exf2,exf2)
b2 = covfun(a2[1],a2,varf2,0.06)
c2 = covfun(a2[2],a2,varf2,0.06)
d2 = covfun(a2[3],a2,varf2,0.06)
e2 = covfun(a2[4],a2,varf2,0.06)
f2 = covfun(a2[5],a2,varf2,0.06)
g2 = covfun(a2[6],a2,varf2,0.06)
h2 = covfun(a2[7],a2,varf2,0.06)
i2 = covfun(a2[8],a2,varf2,0.06)
j2 = covfun(a2[9],a2,varf2,0.06)
k2 = covfun(a2[10],a2,varf2,0.06)
varD2 = matrix(c(b2,c2,d2,e2,f2,g2,h2,i2,j2,k2),nrow=10,ncol=10,byrow=TRUE)
covfD2 = function(x){
  g = covfun(x,a2,varf2,0.06)
  return(g)
}
covDf2 = function(x){
  i = covfun(a2,x,varf2,0.06)
  return(i)
}
eDf2 = function(x){
  h = exf2 + covfD2(x)%*%solve(varD2)%*%(D2-eD2)
  return(as.vector(h))
}
varDf2 = function(x){
  j = varf2 - covfD2(x)%*%solve(varD2)%*%covDf2(x)
  return(j)
}
sequem4 = c()
for (i in sequem){
  sequem4 = append(sequem4,eDf2(i))
}
sequem5 = c()
for (i in sequem){
  sequem5 = append(sequem5,varDf2(i))
}
upvar2 = sequem4 + 3*sqrt(sequem5)
downvar2 = sequem4 - 3*sqrt(sequem5)
plot(sequem,sequem4,type='l',col='blue',xlab="Rate parameter value x",ylab="Concentartion of f(x)")
points(a2,D2,col='purple',pch=19)
lines(sequem,upvar2,type='l',col='red')
lines(sequem,downvar2,type='l',col='red')

