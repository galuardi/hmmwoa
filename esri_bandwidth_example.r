
# example bandwidth calculation from a likelihood raster
# this one is from ESRI: http://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/how-kernel-density-works.htm

library(aspace)

set.seed(123)
r <- raster(matrix(abs(rnorm(9)), nrow=3, ncol=3), xmn=-1, xmx=1, ymn=-1, ymx=1)
r.pts <- rasterToPoints(r, spatial=TRUE)

sdd <- aspace::calc_sdd(centre.xy=c(0,0), calccentre=FALSE,
         weighted=TRUE, weights=as.data.frame(r.pts)[,1], points=as.data.frame(r.pts)[,c(2,3)], verbose=TRUE)$SDD.radius

mc <- mean_centre(weighted=TRUE, weights=as.data.frame(r.pts)[,1], points=as.data.frame(r.pts)[,c(2,3)])
Dm <- median(distances(centre.xy = mc[,2:3], destmat = as.data.frame(r.pts)[,c(2,3)]))

term2 <- sqrt(1 / log(2)) * Dm
n <- cellStats(r, 'sum')

if(sdd < term2){
  bnd <- .9 * sdd * n^(-0.2)
} else{
  bnd <- .9 * term2 * n^(-0.2)
}

