Cx = function(B0,sigmau,theta,points){
  a = points
  D = cos(a[,1])+sin(a[,2])
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=theta))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=theta))
    return(as.vector(i2))
  }
  ju = cbind(D,eD)
  sub = ju[,1]-ju[,2]
  svarD = solve(varD)
  eDf = function(x){
    h = exf + covfD(x)%*%svarD%*%(sub)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  diff = function(x){
    pol = (eDf(x)-f(x))/(sqrt(varDf(x)))
    return(pol)
  }
  sequem = expand.grid(x=seq(-0.5,(2*pi)+0.5,length.out=40),
                       y=seq(-0.5,(2*pi)+0.5,length.out=40))
  sequem4 = mapply(diff,as.data.frame(t(sequem)))
  sequem4mat = matrix(sequem4,40)
  return(sequem4mat)
}

testem = Cx(0,0.6,2,latpoints)

filled.contour(x = seq(-0.5,(2*pi)+0.5,length.out=40),
               y = seq(-0.5,(2*pi)+0.5,length.out=40), 
               z = abs(testem),
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,7.5),
               plot.title = title(main = expression(abs(C(x[T]))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

###################################################################################
Cx2 = function(B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(minmaxfun,as.data.frame(t(a))))
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=theta))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=theta))
    return(as.vector(i2))
  }
  ju = cbind(D,eD)
  sub = ju[,1]-ju[,2]
  svarD = solve(varD)
  eDf = function(x){
    h = exf + covfD(x)%*%svarD%*%(sub)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  diff = function(x){
    pol = (eDf(x)-minmaxfun(x))/(sqrt(varDf(x)))
    return(pol)
  }
  sequem = expand.grid(x=seq(4.5,15.5,length.out=40),
                       y=seq(4.5,15.5,length.out=40))
  sequem4 = mapply(diff,as.data.frame(t(sequem)))
  sequem4mat = matrix(sequem4,40)
  return(sequem4mat)
}

minmaxpoints = matrix(c(1,4.5,8,1,4.5,8,1,4.5,8,1,1,1,4.5,4.5,4.5,8,8,8),nrow=9,ncol=2)
minmaxpoints
testem1 = Cx2(0,1.8,3,minmaxpoints)

filled.contour(x = seq(0,3*pi,length.out=40),
               y = seq(0,3*pi,length.out=40), 
               z = abs(testem1),
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,7.5),
               plot.axes = {axis(1,seq(from=0,to=9,by=2))
                 axis(2,seq(from=0,to=9,by=2))
                 points(minmaxpoints,col="black",pch=19)},
               plot.title = title(main = expression(abs(C(x[T]))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))

#########################################################################################

Cx3 = function(B0,sigmau,theta,points){
  a = points
  D = as.vector(mapply(robotforlooks,as.data.frame(t(a))))
  exf = B0
  varf = sigmau^2
  eD = rep(exf,length.out = length(D))
  b = as.matrix(dist(a,diag=TRUE))
  varD = (sigmau^2)*exp(-(b^2)/(theta^2))
  covfD = function(x){
    g = mapply(covfun2d,as.data.frame(t(a)),
               MoreArgs = list(x=x,su=varf,t=theta))
    return(as.vector(g))
  }
  covDf = function(x){
    i2 = mapply(covfun2d,as.data.frame(t(a)),
                MoreArgs = list(y=x,su=varf,t=theta))
    return(as.vector(i2))
  }
  ju = cbind(D,eD)
  sub = ju[,1]-ju[,2]
  svarD = solve(varD)
  eDf = function(x){
    h = exf + covfD(x)%*%svarD%*%(sub)
    return(as.vector(h))
  }
  varDf = function(x){
    j = varf - covfD(x)%*%svarD%*%covDf(x)
    return(j)
  }
  diff = function(x){
    pol = (eDf(x)-robotforlooks(x))/(sqrt(varDf(x)))
    return(pol)
  }
  sequem = expand.grid(x=seq(0,1,length.out=40),
                       y=seq(0,1,length.out=40))
  sequem4 = mapply(diff,as.data.frame(t(sequem)))
  sequem4mat = matrix(sequem4,40)
  return(sequem4mat)
}
robotogpoints = matrix(c(0.1,0.5,0.9,0.1,0.5,0.9,0.1,0.5,0.9,0.1,0.1,0.1,0.5,0.5,0.5,0.9,0.9,0.9),nrow=9,ncol=2)
robotogpoints
testem2 = Cx3(3,0.5,0.33,robotogpoints)
filled.contour(x = seq(0,1,length.out=40),
               y = seq(0,1,length.out=40), 
               z = abs(testem2),
               color.palette = colorRampPalette(c('green','yellow','red')),
               levels = c(0,0.5,1,1.5,2,2.5,2.75,3,3.25,4,5,6,7,7.5),
               plot.axes = {axis(1,seq(from=0,to=1,by=0.1))
                 axis(2,seq(from=0,to=1,by=0.1))
                 points(robotogpoints,col="black",pch=19)},
               plot.title = title(main = expression(abs(C(x[T]))),
                                  xlab = expression(x[1]), 
                                  ylab = expression(x[2])))
