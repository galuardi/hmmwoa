log.like.sst <-
function(ingrid, datax = as.numeric(vec[11]), sigma = as.numeric(vec[12])){
if(is.na(sigma)) sigma = 5
tempL = ingrid*0  #matrix(0,dim(ingrid)[1],dim(ingrid)[2])
tdim = dim(ingrid)
tempL = dlnorm( log(ingrid), log(datax),log(sigma), log =F)
tempL
}
