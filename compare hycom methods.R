
# put in two L's: L.ohc and L.prof
# these are from the mako257 example and are already currently loaded
str(L.ohc); str(L.prof)

# get SPOT positions for as many times as possible
# SPOT also already loaded
str(spot)

df <- data.frame(matrix(ncol=7, nrow=T))
names(df) <- list('dateOnly','lon','lat','extr','mc.lon','mc.lat','rmse')
spotDates <- as.Date(spot$date)

for (t in 1:5){
  spot.idx <- which(spotDates %in% dateVec[t])
  if(length(spot.idx > 0)){
    spot.i <- spot[min(spot.idx),]
    r.ex <- raster::extract(L.ohc[[t]], cbind(spot.i$lon,spot.i$lat))
    r.pts <- rasterToPoints(L.ohc[[t]], spatial=TRUE)
    
    # generate mean weighted centre at each time T
    mc <- mean_centre(weighted=TRUE, weights=as.data.frame(r.pts)[,1], points=as.data.frame(r.pts)[,c(2,3)])
    
    # compare centre of L to day's spot position and do RMSE
    df[t,1] <- paste(year(spot.i$date),'-',month(spot.i$date),'-',day(spot.i$date),sep='')
    df[t,2:3] <- c(spot.i$lon, spot.i$lat)
    df[t,4] <- r.ex
    df[t,5:6] <- mc[,2:3]
  }
   
}
