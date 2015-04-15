make.btrack.alldates <-
function(track, psat){
	require(date)
	alldates  = seq(psat$day0, psat$dayT)
	alldates.dmy  = date.mdy(alldates)
	alldates.dmy = cbind(alldates.dmy[[2]],alldates.dmy[[1]],alldates.dmy[[3]])

	bdates = mdy.date(track[,2],track[,3],track[,1])

	bidx = match(bdates, alldates)

	output = as.data.frame(matrix(NA,nrow = length(alldates), ncol = ncol(track)))
	names(output) = names(track)
	
	output[,1:3] = alldates.dmy
	output[bidx,4:ncol(track)] = track[,4:ncol(track)]
	output
}
