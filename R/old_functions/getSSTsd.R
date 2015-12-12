getSSTsd <-
function(psat, surf = -11){
	tzidx = psat$Z>=(surf)
	tzidx[is.na(tzidx)]=F
	temp = psat$T
	temp[tzidx==F] = NA
	tsd = apply(temp,1,sd,na.rm=T)
	tsd[tsd==0] = NA
	tsd
}
