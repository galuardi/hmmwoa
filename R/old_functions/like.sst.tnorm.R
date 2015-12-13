like.sst.tnorm <-
function(ingrid, datax = as.numeric(vec[11]), sigma = as.numeric(vec[12])){
if(is.na(sigma)) sigma = 5
tempL = ingrid*0  #matrix(0,dim(ingrid)[1],dim(ingrid)[2])
tdim = dim(ingrid)
ingrid2 = ingrid
ingrid2[is.nan(ingrid2)] = 0
tempL = matrix(dtnorm(ingrid2, datax, sigma, datax-2, 40), tdim[1], tdim[2])
tempL[tempL==0] = NaN
tempL
}
