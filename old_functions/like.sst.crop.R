like.sst.crop <-
function(zgrid=zgrid, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, lonvec = NULL, latvec = NULL, datax = vec, sigma = 2){
	vec = as.numeric(datax)
	data.lon = as.numeric(vec[8])
	data.lat = as.numeric(vec[9])
	data.sst = as.numeric(vec[11])
	lat = latvec
	lon = lonvec

	xmin = data.lon-3 #if(is.null(xmin)==T) 
	 xmax = data.lon+3 #if(is.null(xmax)==T)
	ymin = data.lat-3 #if(is.null(ymin)==T) 
	 ymax = data.lat+3 #if(is.null(xmax)==T)
	  bxlim = c(xmin, xmax)
	  bylim = c(ymin, ymax)
	  blon = which(lon>=bxlim[1]&lon<=bxlim[2])
	  blat = which(lat>=bylim[1]&lat<=bylim[2])
	  bout=zgrid[blon, blat]
	lout = zgrid*0
	lout[blon, blat] = like.sst(bout, datax = data.sst, sigma = sigma)
	lout
}
