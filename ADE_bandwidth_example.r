
#----------------
# REAL LIKELIHOOD EXAMPLE FOR BANDWIDTH CALCS
# from day 2 of blue259 likelihood
# used here to compare 2 different bandwidth approaches
#----------------

# load example likelihood raster, r
load('example_likelihood_raster.RData')

# ADE bandwidth approach
rudf = as.data.frame(as(r, 'SpatialPointsDataFrame'))
rudf[,1] = rudf[,1]/max(rudf[,1])

iter <- seq(10,10010,by=100)
hvec <- rep(NA, length(iter))
for (n in iter){
  #n = 10 # samples to draw
  samp = rudf[sample(1:nrow(rudf), n, prob = rudf[,1], replace = T),2:3]
  coordinates(samp) = ~x+y
  sampud = kernelUD(samp)
  
  # compare h values
  hvec[which(iter == n)] <- sampud@h$h
  #rbind(ud@h$h, sampud@h$h)
  
}

# ESRI bandwidth approach
library(aspace)

r.pts <- rasterToPoints(r, spatial=TRUE)

rdf <- as.data.frame(r.pts)
rdf[,1] <- rdf[,1] / max(rdf[,1], na.rm=T)

source('../hmmwoa/R/esri_h.r')
h <- esri_h(rdf)

pdf('bandwidth_example.pdf',height=12,width=8)
par(mfrow=c(2,1))
plot(r)
world(add=TRUE)
plot(hvec~iter,type='l',xlab='# iterations',ylab='bandwidth, h')
points(h~1,pch=16,col='blue')
dev.off()

#----------------
# the question remains, how do we get back to the movement parameters from these kernel metrics?
# here's an example using movement kernel for resident behavior, K2
D2 <- c(3,1)
K2 <- gausskern(D2[1], D2[2], muadv = 0) # normally we then convert this into an image
k2r <- raster(K2, xmn=-1, xmx=1, ymn=-1, ymx=1)

# ADE bandwidth approach
rudf = as.data.frame(as(k2r, 'SpatialPointsDataFrame'))
rudf[,1] = rudf[,1]/max(rudf[,1])

iter <- seq(10,10010,by=100)
hvec <- rep(NA, length(iter))
for (n in iter){
  #n = 10 # samples to draw
  samp = rudf[sample(1:nrow(rudf), n, prob = rudf[,1], replace = T),2:3]
  coordinates(samp) = ~x+y
  sampud = kernelUD(samp)
  
  # compare h values
  hvec[which(iter == n)] <- sampud@h$h
  #rbind(ud@h$h, sampud@h$h)
  
}

plot(hvec~iter,type='l',xlab='# iterations',ylab='bandwidth, h')

# h ~ 0.1-0.2 for resident kernel 
# how do we go from here to a new set of movement parameters?



#------------
## END



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


