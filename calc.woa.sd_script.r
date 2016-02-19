# calculate sd on a grid

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

