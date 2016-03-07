load('~/Documents/WHOI/RData/WhiteSharks/2013/121325/integrate example.RData')

lik.pdt = array(NA, dim=c(dim(dat.i)[1], dim(dat.i)[2], length(depIdx)))

# works
b=2
lik.pdt[,,b] = likint2(dat.i[,,depIdx[b]], sd.i[,,depIdx[b]], df[b,1], df[b,2])
image.plot(lik.pdt[,,b])

# doesn't work
b=3
lik.pdt[,,b] = likint2(dat.i[,,depIdx[b]], sd.i[,,depIdx[b]], df[b,1], df[b,2])
image.plot(lik.pdt[,,b])

# add TINY difference anywhere, works again
lik.pdt[,,b] = likint2(dat.i[,,depIdx[b]], sd.i[,,depIdx[b]], df[b,1]+.001, df[b,2])
image.plot(lik.pdt[,,b])

