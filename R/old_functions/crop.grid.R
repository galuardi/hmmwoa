crop.grid <-
function(xmin, xmax, ymin, ymax, ingrid){
  bxlim = c(xmin, xmax)
  bylim = c(ymin, ymax)
  blon = which(ingrid$lon>=bxlim[1]&ingrid$lon<=bxlim[2])
  blat = which(ingrid$lat>=bylim[1]&ingrid$lat<=bylim[2])
  bath.out=list()
  bath.out$lon=ingrid$lon[blon]
  bath.out$lat=ingrid$lat[blat]
  bath.out$data=t(ingrid$data[blat, blon])
  bath.out
}
