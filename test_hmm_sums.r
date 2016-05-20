# fix the last day's weird probability sum...

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f2 = hmm.filter_test(g,L,K1,K2,P)
# res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
# fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s2 = hmm.smoother_test(f2, K1, K2, P, plot = F)
apply(s2, c(2), sum)

# PLOT IT IF YOU WANT TO SEE LIMITS (CI)
sres = apply(s2[1,,,], 2:3, sum, na.rm=T)
fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
distr = s2
meanlat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)


# Things to look at... 
# L
# L.locs
# 2:(T-1); (T-1:2) in filter/smooth.. need to leave the first/last day alone while updating all others. but, we want to use first/last days info going forward and backward...


