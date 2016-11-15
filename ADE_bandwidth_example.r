library(adehabitatHR)
library(raster)

Sigma <- matrix(c(10,3,3,2),2,2)
pts = mvrnorm(n = 1000, rep(0, 2), Sigma)
pts = as.data.frame(pts)
names(pts) = c('x','y')
coordinates(pts) = ~x+y
ud = kernelUD(pts)
ud@h

# now make a raster with values
rud = raster(ud)

rudf = as.data.frame(as(rud, 'SpatialPointsDataFrame'))
rudf[,1] = rudf[,1]/max(rudf[,1])

n = 1000 # samples to draw
samp = rudf[sample(1:nrow(rudf), n, prob = rudf[,1], replace = T),2:3]

coordinates(samp) = ~x+y
sampud = kernelUD(samp)

# compare h values

rbind(ud@h$h, sampud@h$h)

# if less samples are drawn, h increases...

n = 10 # samples to draw
samp = rudf[sample(1:nrow(rudf), n, prob = rudf[,1], replace = T),2:3]

coordinates(samp) = ~x+y
sampud = kernelUD(samp)

# compare h values

rbind(ud@h$h, sampud@h$h)
