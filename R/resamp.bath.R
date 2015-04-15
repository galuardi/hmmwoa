resamp.bath <-
function(sstL, bathy){
	require(raster)
   xmin = min(sstL$lon)
   xmax = max(sstL$lon)
   ymin = min(sstL$lat)
   ymax = max(sstL$lat)	
   
   b2 = crop.grid(xmin, xmax, ymin, ymax, bathy)
   
	s = raster(nrow=nrow(sstL$sstL),ncol=ncol(sstL$sstL)) # dimensions of the SST 
	 bath2=raster(ncol = ncol(b2[[3]]), nrow = nrow(b2[[3]]))
	 
	 #extent(bath2)=extent(c(min(b2$lon), max(b2$lon), min(b2$lat), max(b2$lat)))
	 
	 bath2 = setValues(bath2, as.vector(t(b2[[3]])))
	 
	 s <- resample(bath2, s, method='bilinear')
	 
	 extent(s)=extent(c(min(sstL$lon)-360, max(sstL$lon)-360, min(sstL$lat), max(sstL$lat)))
	 
	 #plot(s, col = bath.colors(100), zlim=c(-250,0))
	 
	bath2 = list(lon = sstL$lon-360, lat = sstL$lat,data = matrix(s@data@values, nrow=nrow(sstL$sstL),ncol=ncol(sstL$sstL), byrow = T)) 
	bath2
 }
