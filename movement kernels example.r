#----------------
# the question remains, how do we get back to the movement parameters from these kernel metrics?
# here's an example using the movement kernels, K1 = migratory K2 = resident
D1 <- c(9,10); D2 <- c(3,1) # let's call these units degrees, for now
K1 <- gausskern(D1[1], D1[2], muadv = 0) # normally we then convert this into an image
K2 <- gausskern(D2[1], D2[2], muadv = 0) # normally we then convert this into an image
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
k1r <- raster(K1, xmn=(0-D1[1]/2), xmx=(0+D1[1]/2), ymn=(0-D1[1]/2), ymx=(0+D1[1]/2),crs=crs)
k2r <- raster(K2, xmn=(0-D2[1]/2), xmx=(0+D2[1]/2), ymn=(0-D2[1]/2), ymx=(0+D2[1]/2),crs=crs)

# ADE bandwidth approach
rdf1 = as.data.frame(as(k1r, 'SpatialPointsDataFrame'))
rdf1[,1] = rdf1[,1]/max(rdf1[,1])

rdf2 = as.data.frame(as(k2r, 'SpatialPointsDataFrame'))
rdf2[,1] = rdf2[,1]/max(rdf2[,1])

iter <- seq(10,10010,by=100)
hvec1 <- rep(NA, length(iter))
hvec2 <- rep(NA, length(iter))

for (n in iter){
  #n = 10 # samples to draw
  samp = rdf1[sample(1:nrow(rdf1), n, prob = rdf1[,1], replace = T),2:3]
  coordinates(samp) = ~x+y
  sampud = kernelUD(samp)
  
  # compare h values
  hvec1[which(iter == n)] <- sampud@h$h
  #rbind(ud@h$h, sampud@h$h)
  
  samp = rdf2[sample(1:nrow(rdf2), n, prob = rdf2[,1], replace = T),2:3]
  coordinates(samp) = ~x+y
  sampud = kernelUD(samp)
  
  hvec2[which(iter == n)] <- sampud@h$h
  
  
}
rdf1[,1] <- rdf1[,1] / max(rdf1[,1])
rdf2[,1] <- rdf2[,1] / max(rdf2[,1])
h1 <- esri_h(rdf1)
h2 <- esri_h(rdf2)

plot(hvec1~iter,type='l',xlab='# iterations',ylab='bandwidth, h', ylim=c(0,1.6))
lines(hvec2~iter, col='blue')
points(h1~1,pch=16)
points(h2~1, pch=16,col='blue')

# area returns area of each cell in km2, sum it, sqrt and / 111 to get back to diffusion radius of that kernel area
D1.new <- D1; D2.new <- D2
D1.new[1] <- sqrt(raster::cellStats(raster::area(k1r), 'sum')) / 111 
D2.new[1] <- sqrt(raster::cellStats(raster::area(k2r), 'sum')) / 111 

# now what about bandwidth?


# h ~ 0.1-0.2 for resident kernel 
# how do we go from here to a new set of movement parameters?


# here's an example using the movement kernels, K1 = migratory K2 = resident
D1 <- c(9,10); D2 <- c(9,1) # let's call these units degrees, for now
K1 <- gausskern(D1[1], D1[2], muadv = 0) # normally we then convert this into an image
K2 <- gausskern(D2[1], D2[2], muadv = 0) # normally we then convert this into an image
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
k1r <- raster(K1, xmn=(0-D1[1]/2), xmx=(0+D1[1]/2), ymn=(0-D1[1]/2), ymx=(0+D1[1]/2),crs=crs)
k2r <- raster(K2, xmn=(0-D2[1]/2), xmx=(0+D2[1]/2), ymn=(0-D2[1]/2), ymx=(0+D2[1]/2),crs=crs)
k3r <- raster(K3, xmn=(0-D1[1]/2), xmx=(0+D1[1]/2), ymn=(0-D1[1]/2), ymx=(0+D1[1]/2),crs=crs)

par(mfrow=c(1,2))
plot(k1r); title('D1 <- c(9,10)')
plot(k2r); title('D1 <- c(9,1)')


rd = rnorm(1000,10,2)
plot(rd)
plot(density(rdf[,1]))
lines(density(rdf[,1], bw = sd(rd)), col = 2)
rug(rdf[,1])
plot(density(rdf[,1], bw = sd(rd)), col = 2)
lines(density(rdf[,1]))
