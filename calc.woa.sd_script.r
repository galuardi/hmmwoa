# calculate sd on a grid
load()

library(raster)
library(fields)

fmat = list()

for(i in 1:57){
  # for(j in 1:57){
    r = flip(raster(t(dat[,,i,1])),2)
    plot(r, col = tim.colors(100))
    f1 = focal(r, w=matrix(1/9,nrow=3,ncol=3), fun=sd)
    plot(f1, add=T)
    fmat[[i]] = f1
  # }
}
cellStats(stack(fmat), 'summary')

pdf(height = 8, width = 10, file = 'woasd.pdf')

for(i in 1:57){
  par(mfrow=c(1,2))
  par(mar = c(2,2,2,6))
  r = flip(raster(t(dat[,,i,1])),2)
  plot(r, col = tim.colors(100))
  title(paste0(depth[i], ' m'))
  plot(fmat[[i]])
  title('focal sd')
}

dev.off()

# compare to imager
gausskern <-
  function(siz, sigma, muadv = 0){
    x = 1:round(siz);
    mu = c(mean(x), mean(x)) + muadv;
    fx = (matrix(exp((-0.5*(x-mu[1])/sigma)^2))/(sqrt(2*pi)*sigma));
    options(digits=5)
    fx = exp(-.5*((x-mu[1])/sigma)^2)/sqrt((2*pi)*sigma)
    fy = exp(-.5*((x-mu[2])/sigma)^2)/sqrt((2*pi)*sigma)
    fx[!is.finite(fx)] = 0
    #fy = (matrix(exp((-0.5*((x-mu[2])/sigma))^2))/(sqrt(2*pi)*sigma));
    fy[!is.finite(fy)] = 0
    kern = (fx%*%t(fy))
    kern = kern/(sum(sum(kern,na.rm=T),na.rm=T))
    kern[is.nan(kern)]=0
    kern
  }


# Just an example... 
d1 = 10
gk = gausskern(d1, 3) 
gk = (array(gk, dim = c(10,10,1,1)))


## Compare focal w/ convolve

r = flip(raster(t(dat[,,1,1])),2)
f1 = focal(r, w = gk[1:9,1:9,1,1], fun = mean, na.rm=T) #, fun = sum, na.rm=T
f1 = focal(r, w=matrix(1,nrow=9,ncol=9), fun = mean, na.rm=T)
f1 = t(as.matrix(flip(f1,2)))
ssti = as.cimg((fliplr(dat[,,1,1])))
ssti[is.na(ssti)] = 1e-15
sstc = fliplr(as.matrix(convolve(ssti, gk)))
sstc[sstc==1e-15] = NA

par(mfrow=c(2,2))
image.plot(dat[,,1,1])
contour(dat[,,1,1], add=T, col='white')
title('woa Jan, 0 depth')
image.plot(f1, col = tim.colors(100))
contour(f1, add=T)
title('FOCAL')
image.plot(sstc)
contour(sstc, add=T, col ='white')
title('CONV')
image.plot(sstc-f1, zlim = c(-2,2))
contour(f1, add=T)
contour(sstc, add=T, col ='white')
title('COMPARE')

library(spatilfil)

K <- convKernel(sigma = 3, k = 'gaussian')
Mfil <- applyFilter(x = M, kernel = convKernel(sigma = 1.4, k = 'gaussian'))

