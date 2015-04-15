pixtomap <-
function (M,px,py){
	lgt = length(px);
	latlong = cbind((numeric(lgt)+1), py, (numeric(lgt)+1), px) %*% M;
	lat  = latlong[,2];
	long = latlong[,1];
	cbind(long,lat)
}
