#pdf("Nansproducedfixed.pdf",width = 7.6,height = 7)

vartestexpand = function(a){
  varf = 0.25^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=2))
    return(as.vector(i2))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=20),
                       y=seq(-0.5,(2*pi)+0.5,length.out=20))
  sequem3 = mapply(varDf,as.data.frame(t(sequem)))
  sequem3mat = matrix(sequem3,20)
  return(sequem3mat)
}

hyp = function(x){
  return(2*pi*maximinLHS(9,2,optimize.on = "grid"))
}

hypl = lapply(seq(100),hyp)
hypvarexpand = lapply(hypl,vartestexpand)
maxvarexpand = lapply(hypvarexpand,max)
optsampexpand = min(unlist(maxvarexpand))
optpointsexpand = hypl[[which(maxvarexpand==optsampexpand)]]

length(optpointsexpand[,1])

varrunexpand = emufuncosxplussiny(0,0.25,2,optpointsexpand)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = varrunexpand[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(varrunexpand[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(optpointsexpand,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
varIexpand = mapply(impcosxplussiny,as.data.frame(t(sequem)),
               MoreArgs = list(z=1,points=optpointsexpand))
varImatexpand = matrix(varIexpand,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = varImatexpand, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(varImatexpand)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(optpointsexpand,col="black",pch=19)})

length(which(varIexpand<3))
possIless = as.vector(which(varIexpand<3))
pointpossIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                             y=seq(-0.5,(2*pi)+0.5,length.out=40))[possIless,]

pointpossIless
randIless = list()
for(i in seq(100)){
  randIless[[i]] = matrix(rbind(optpointsexpand,as.matrix(pointpossIless[sample(nrow(pointpossIless),
                                         size=4),])),13,2)
}
#randIless

randvarexpand = lapply(randIless,vartestexpand)
randmaxvarexpand = lapply(randvarexpand,max)
randoptsampexpand = min(unlist(randmaxvarexpand))
randoptpointsexpand = randIless[[which(randmaxvarexpand==randoptsampexpand)]]
#randexpand.2 = as.matrix(rbind(optpointsexpand,as.matrix(randoptpointsexpand)))
#randexpand.2
testrand = matrix(c(4.2092002,2.3556705,2.4087874,4.2934502,
                    3.1514164,3.5063800,1.2706833,5.9837297,
                    0.3593197,0.7364953,4.1318984,1.5406247,
                    5.4470498,5.3897040,5.5997081,0.4011520,
                    1.7559731,3.4023346,5.4759469,2.8614701,
                    0.4337417,3.0482185,6.0,2.9,
                    5.6,2.5),13,2,byrow = TRUE)
randoptpointsexpand
#randexpand.21

randvarrunexpand = emufuncosxplussiny(0,0.25,2,randoptpointsexpand)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = randvarrunexpand[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(randvarrunexpand[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(randoptpointsexpand,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
#which(randvarrunexpand[[3]]<0)
#randvarIexpand = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    #MoreArgs = list(z=1,points=randexpand.2))
#randvarImatexpand = matrix(randvarIexpand,40)
#filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               #y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               #z = randvarImatexpand, 
               #color.palette = colorRampPalette(c('green','yellow','red')),
               #levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(randvarImatexpand)),
               #plot.title = title(main = expression(paste("Implausability measure I")),
                                  #xlab = expression(x[1]), 
                                  #ylab = expression(x[2])),
               #plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 #axis(2,seq(from=0,to=6,by=1))
                 #points(randexpand.2,col="black",pch=19)})
#dev.off()

########################################################################################################

meanvarexpand = lapply(hypvarexpand,mean)
moptsampexpand = min(unlist(meanvarexpand))
moptpointsexpand = hypl[[which(meanvarexpand==moptsampexpand)]]

mvarrunexpand = emufuncosxplussiny(0,0.25,2,moptpointsexpand)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = mvarrunexpand[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(mvarrunexpand[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(moptpointsexpand,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
mvarIexpand = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=moptpointsexpand))
mvarImatexpand = matrix(mvarIexpand,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = mvarImatexpand, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(mvarImatexpand)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(moptpointsexpand,col="black",pch=19)})

##########################################################################

quickvarexpand = function(a){
  varf = 0.25^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=20),
                       y=seq(-0.5,(2*pi)+0.5,length.out=20))
  vari = rep(varf,times=400)
  sequem3.1 = mapply(covfD,as.data.frame(t(sequem)))
  sequem3.2 = matrix(sequem3.1,400,9,byrow = TRUE)
  diagcal = rowSums((sequem3.2 %*% svarD) * sequem3.2)
  resdiag = vari - diagcal
  resdiagmat = matrix(resdiag,20)
  return(resdiagmat)
}

hyplte = lapply(seq(100),hyp)
qhypvarexpand = lapply(hyplte,quickvarexpand)
qmaxvarexpand = lapply(qhypvarexpand,max)
qoptsampexpand = min(unlist(qmaxvarexpand))
qoptpointsexpand = hyplte[[which(qmaxvarexpand==qoptsampexpand)]]
qoptpointsexpand

qvarrunexpand = emufuncosxplussiny(0,0.25,2,qoptpointsexpand)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qvarrunexpand[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(qvarrunexpand[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qoptpointsexpand,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
qvarIexpand = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                    MoreArgs = list(z=1,points=qoptpointsexpand))
qvarImatexpand = matrix(qvarIexpand,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qvarImatexpand, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(qvarImatexpand)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qoptpointsexpand,col="black",pch=19)})

length(which(qvarIexpand<3))
qpossIless = as.vector(which(qvarIexpand<3))
qpointpossIless = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                             y=seq(-0.5,(2*pi)+0.5,length.out=40))[qpossIless,]

qpointpossIless

qrandIless = list()
for(i in seq(100)){
  qrandIless[[i]] = matrix(rbind(qoptpointsexpand,as.matrix(qpointpossIless[sample(nrow(qpointpossIless),
                                                                                size=4),])),13,2)
}
qrandIless

quickvarexpand.2 = function(a,sequem){
  varf = 0.25^2
  b = as.matrix(dist(a,diag=TRUE))
  varD = varf*exp(-(b^2)/(2^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=2))
    return(as.vector(g))
  }
  svarD = solve(varD)
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  vari = rep(varf,times=length(sequem))
  sequem3.1 = mapply(covfD,as.data.frame(t(sequem)))
  sequem3.2 = matrix(sequem3.1,length(sequem),length(a[,1]),byrow = TRUE)
  diagcal = rowSums((sequem3.2 %*% svarD) * sequem3.2)
  resdiag = vari - diagcal
  return(resdiag)
}
qhypvarexpand.2 = lapply(qrandIless,quickvarexpand.2,sequem = qpointpossIless)
qmaxvarexpand.2 = lapply(qhypvarexpand.2,max)
qoptsampexpand.2 = min(unlist(qmaxvarexpand.2))
qoptpointsexpand.2 = qrandIless[[which(qmaxvarexpand.2==qoptsampexpand.2)]]
qoptpointsexpand.2

qvarrunexpand.2 = emufuncosxplussiny(0,0.25,2,qoptpointsexpand.2)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qvarrunexpand.2[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(qvarrunexpand.2[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qoptpointsexpand.2,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
qvarIexpand.2 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                     MoreArgs = list(z=1,points=qoptpointsexpand.2))
qvarImatexpand.2 = matrix(qvarIexpand.2,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qvarImatexpand.2, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(qvarImatexpand.2)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qoptpointsexpand.2,col="black",pch=19)})

length(which(qvarIexpand.2<3))
which()
qpossIless.3 = as.vector(which(qvarIexpand.2<3))
qpointpossIless.3 = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                              y=seq(-0.5,(2*pi)+0.5,length.out=40))[qpossIless.3,]

qpointpossIless.3
qrandIless.3 = list()
for(i in seq(100)){
  qrandIless.3[[i]] = matrix(rbind(qoptpointsexpand.2,as.matrix(qpointpossIless.3[sample(nrow(qpointpossIless.3),
                                                                                   size=2),])),15,2)
}

qhypvarexpand.3 = lapply(qrandIless.3,quickvarexpand.2,sequem = qpointpossIless.3)
qmaxvarexpand.3 = lapply(qhypvarexpand.3,max)
qoptsampexpand.3 = min(unlist(qmaxvarexpand.3))
qoptpointsexpand.3 = qrandIless.3[[which(qmaxvarexpand.3==qoptsampexpand.3)]]
qoptpointsexpand.3

qvarrunexpand.3 = emufuncosxplussiny(0,0.25,2,qoptpointsexpand.3)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qvarrunexpand.3[[2]],
               col = rainbow(25,alpha=0.75),
               plot.title = title(main = expression(paste("Emulator Expectation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = sqrt(qvarrunexpand.3[[3]]),
               col = rainbow(25,alpha=0.95),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qoptpointsexpand.3,col="black",pch=19)},
               plot.title = title(expression(paste("Emulator Standard Deviation")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
qvarIexpand.3 = mapply(impcosxplussiny,as.data.frame(t(sequem)),
                       MoreArgs = list(z=1,points=qoptpointsexpand.3))
qvarImatexpand.3 = matrix(qvarIexpand.3,40)
filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = qvarImatexpand.3, 
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,9,11,14,16,18,20,max(qvarImatexpand.3)),
               plot.title = title(main = expression(paste("Implausability measure I")),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])),
               plot.axes = {axis(1,seq(from=0,to=6,by=1))
                 axis(2,seq(from=0,to=6,by=1))
                 points(qoptpointsexpand.3,col="black",pch=19)})
