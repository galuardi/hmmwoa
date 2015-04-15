normalise <-
function(ingrid){
	normConst = sum(ingrid,na.rm=T)
	ingrid/normConst
}
