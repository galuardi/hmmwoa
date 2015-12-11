# should be kernel 1 mean, kernel 1 sd, kernel 2 mean, kernel 2 sd, ?, ? Bmode??
#-------------------------------------------#
#run run_lydia_example_bg first!
#-------------------------------------------#
library(imager)

# sstf = open.ncdf('C:/Users/ben/Google Drive/SCDNR/ANALYSIS/SSTNEW (1)/37078-SST/oisst.nc')
# sst = get.var.ncdf(sstf, 'sst')
# sstr = stack('C:/Users/ben/Google Drive/SCDNR/ANALYSIS/SSTNEW (1)/37078-SST/oisst.nc')
# source('C:/Users/ben/Google Drive/Camrin-WOA/sphmmInR/sim/sstdb.r')
source('C:/Users/ben/Google Drive/Camrin-WOA/sphmmInR/spathmm/sphmmfuns.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\misc_funs.r')
source('C:\\Users\\ben\\Documents\\GitHub\\hmmwoa\\R\\normalise.r')

gausskern <-
	function(siz, sigma, muadv = 0){
	x = 1:round(siz);
	mu = c(mean(x), mean(x)) + muadv;
	fx = (matrix(exp((-0.5*(x-mu[1])/sigma)^2))/(sqrt(2*pi)*sigma));
	options(digits=5)
	fx = exp(-.5*((x-mu[1])/sigma)^2)/sqrt((2*pi)*sigma)
	fy = exp(-.5*((x-mu[2])/sigma)^2)/sqrt((2*pi)*sigma)
	fx[!is.finite(fx)] = 0
	#fy = (matrix(exp((-0.5*((x-mu[2])/sigma))^2))/(sqrt(2*pi)*sigma));
	fy[!is.finite(fy)] = 0
	kern = (fx%*%t(fy))
	kern = kern/(sum(sum(kern,na.rm=T),na.rm=T))
	kern[is.nan(kern)]=0
	kern
}


# Just an example... 
d1 = 10
gk = gausskern(d1, 3) 
gk = (array(gk, dim = c(10,10,1,1)))

ssti = as.cimg(t(as.matrix(sstr[[1]])))
plot(convolve(ssti, gk) )

ssta = as.cimg(array(values(sstr), dim = c(20,21,1,14)))

str(idply(ssta, 'c', fun = function(x) convolve(x, gk)))

#### Older Function based on Pedersen 2008, using gaussian kernel convolution
hmmfilter1 <- function(sstL, fmat, ks=29, sigma=1){
allpost = sstL[[3]]*NaN
lon = L$lon
lat = L$lat
# First position
xidx = which.min((fmat[1,8]-(lon))^2);yidx = which.min((fmat[1,9]-lat)^2);
allpost[xidx,yidx,1] = 1
# Last Position
len = nrow(fmat)
xidx = which.min((fmat[len,8]-(lon))^2);yidx = which.min((fmat[len,9]-lat)^2);
allpost[xidx,yidx,len] = 1

post = allpost[,,1]
post[is.nan(post)] = 0

#=======================================================#
   # P <- lowpass(post,radius =50)  # rimage package

#=======================================================#
# EBImage package
#=======================================================#

datediff = c(1,diff(as.numeric(mdy.date(fmat[,2],fmat[,3],fmat[,1]))))

#n = ks+(datediff[1]-1)*14
f = gausskern(ks[1], sigma[1])

 # f = makeBrush(n, shape='disc', step=FALSE) 
 f = f/sum(f)

post = Image(post)
P <- filter2(post, f)
#=======================================================#
   
post = P@.Data
post = normalise(post)
normaliser = numeric(nrow(fmat)-1)+1
#-------------------------------------#
# Filter loop
#-------------------------------------#
time1 = Sys.time()
for(i in 2:(nrow(fmat)-1)){#:
print(i)
 
if(fmat$behav[i-1]==1){
		sig = sigma[1]; kern = ks[1];
	 }else{
		sig = sigma[2]; kern = ks[2]
	 }
 
   post[is.nan(post)] = 0
   # n = ks+(datediff[i]-1)*14
   # f = makeBrush(ks, shape='disc', step=FALSE) 
   # n = ks+(datediff[1]-1)*14
   f = gausskern(kern, sig)
   f = f/sum(f)
   P <- filter2(post, f)
   P = P@.Data
   #allP[,,i] = P
   post = P*sstL[[3]][,,i]
   normaliser[i-1] = sum(post,na.rm=T)
   post = normalise(post)
   allpost[,,i] = post
}

print(Sys.time()-time1)
 allpost
}

#-------------------------------------------------------------------------------#
## edits to function start here
#-------------------------------------------------------------------------------#
## Fixed parameter values
# par0=c(8.908,10.27,1.152,0.0472,0.707,0.866) # what units are these?
par0=c(8.908,10.27,3,1,0.707,0.866) # what units are these?
D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6] # not sure what these parameters are.. look like the diagonal of a 2x2 transition matrix.  

# Probably need to express kernel movement in terms of pixels per time step. The sparse matrix work likely renders this unnecessary, but going back to gausskern, it is. For example, if we have .25 degree and daily time step, what would the speed of the fish be when moving fast? 4 pixels/day?

K1 = as.cimg(gausskern(D1[1], D1[2], muadv = 0))
K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

# make all NA's very tiny for the convolution
# the previous steps may have taken care of this...
# L[L==0] = 1e-15
# L[is.na(L)] = 1e-15

# add a 'skip' index for missing days in the L.. 

#### Function based on Pedersen 2011, originally using sparse matrix multiplication.. editing for gaussian kernel convolution
# watch out for matrix dimensionality. even though dimensions lined up, data did not.. 
hmm.filter2 <- function(g,L,K1,K2,P){
  require(imager) # convolution function
  require(magic) # has a rotate function.. and isn't matlab
  ## Filter data to estimate locations and behaviour
  
  T <- dim(L)[1] # dimension of time 
  row <- dim(g$lon)[1] # nrows
  col <- dim(g$lon)[2] # ncols
  m <- 2 # Number of behavioural states
  
  pred <- array(0,dim=c(m,T,col,row)) # empty array for prediction step. ordering in col before row emulates lon before lat
  phi  <- array(0,dim=c(m,T,col,row)) # posterior (final) step array
  # Start in resident state at the known initial location
  phi[2,1,,]  <- L[1,,] # first position is known
  pred[2,1,,] <- L[1,,] # first position is known
  psi <- rep(0,T-1) # sum of the probability of both states at each step
  # Start filter iterations
  for(t in 2:T){
	# replace this part with older workflow using a gaussian kernel.. 
  # p1 <- as.vector(phi[1,t-1,,])
    # p2 <- as.vector(phi[2,t-1,,])
    # q1 <- as.vector(p1%*%K1)
    # q2 <- as.vector(p2%*%K2)
    
	p1 = as.cimg(t(phi[1,t-1,,]))
	p2 = as.cimg(t(phi[2,t-1,,]))
	q1 = convolve(p1, K1)
	q2 = convolve(p2, K2)
	
	# q1 = arot(t(as.matrix(q1)),3)
	# q2 = arot(t(as.matrix(q2)),3)
	q1 = t(as.matrix(q1))
	q2 = t(as.matrix(q2))
	
# 	par(mfrow=c(1,2))
# 	image(q1)
# 	image(q2)
	
    # pred[1,t,,] <- matrix(P[1,1]*q1+P[2,1]*q2,row,col)
    # pred[2,t,,] <- matrix(P[1,2]*q1+P[2,2]*q2,row,col)
	
	# multiply by transition probability 
	pred[1,t,,] <- P[1,1]*q1+P[2,1]*q2
    pred[2,t,,] <- P[1,2]*q1+P[2,2]*q2
   
  sumL = sum(L[t,,])  
  if(sumL > 0){
    post1 <- pred[1,t,,]*L[t,,]
    post2 <- pred[2,t,,]*L[t,,]
  }else{
    post1 <- pred[1,t,,]
    post2 <- pred[2,t,,]
  }

    psi[t-1] <- sum(as.vector(post1), na.rm=T) + sum(as.vector(post2), na.rm=T)
	
	# remove NaNs... 
	# normalise (divide here by sum, not max)
# 	post1 <- normalise(post1)
# 	post2 <- normalise(post2)
# 	post1[is.nan(post1)] = 0
# 	post2[is.nan(post2)] = 0
	
    phi[1,t,,] <- post1/(psi[t-1]+1e-15)
    phi[2,t,,] <- post2/(psi[t-1]+1e-15)
  }
  list(phi=phi,pred=pred,psi=psi)
}

f = hmm.filter2(g,L,K1,K2,P)
res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, res/max(res), zlim = c(.05,1)

# Next up.... 
hmm.smoother <- function(f,K1,K2,P){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  smooth <- array(0,dim=dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  for(t in T:2){
    RAT <- smooth[,t,,]/(f$pred[,t,,]+1e-15)
#     Rp1 <- as.vector(K1 %*% as.vector(RAT[1,,]))
#     Rp2 <- as.vector(K2 %*% as.vector(RAT[2,,]))
    
    
    p1 = as.cimg(t(RAT[1,,]))
    Rp1 <- convolve(p1, K1)
    p2 = as.cimg(t(RAT[2,,]))
    Rp2 <- convolve(p2, K2)
    
    Rp1 = t(as.matrix(Rp1))
    Rp2 = t(as.matrix(Rp2))
    
    par(mfrow=c(1,2))
    image.plot(Rp1)
    image.plot(Rp2)
    
    post1 <- matrix(P[1,1]*Rp1 + P[1,2]*Rp2,row,col)
    post2 <- matrix(P[2,1]*Rp1 + P[2,2]*Rp2,row,col)
    post1 <- post1 * f$phi[1,t-1,,]
    post2 <- post2 * f$phi[2,t-1,,]
    fac <- sum(as.vector(post1)) + sum(as.vector(post2))
    smooth[1,t-1,,] <- post1/fac
    smooth[2,t-1,,] <- post2/fac
  }
  smooth
}

s = hmm.smoother(f, K1, K2, P)

sres = apply(s[1,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

sres = apply(s[2,,,],2:3,sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

# calculate track

calc.track(s, g)  # dimensions flipped...

# switch the dimensions in the calc.track.r... gives a weird output.. ON FIN LAND!
# this is either 1) right and we have to deal with the L and K elements or 2) the dimensions need adjusting..

distr = s 
meanlat <- apply(apply(distr,c(2,4),sum)*repmat(t(as.matrix(g$lat[,1])),T,1),1,sum)
meanlon <- apply(apply(distr,c(2,3),sum)*repmat(t(as.matrix(g$lon[1,])),T,1),1,sum)

plot(meanlon, meanlat)
plot(countriesLow, add = T)

sres = apply(s,c(3,4), sum, na.rm=T)
image.plot(lon, lat, sres/max(sres), zlim = c(.01,1))
lines(meanlon, meanlat, pch=19, col=2)
plot(countriesLow, add = T)

sr = raster(sres/max(sres),xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
spot = read.csv('C:/Users/benjamin.galuardi/Google Drive/Camrin-WOA/hmmwoa_files/121325-SPOT.csv')

plot(sr)
plot(countriesLow, add = T)
lines(meanlon, meanlat, pch=19, col=2)
lines(spot$Longitude, spot$Latitude, typ='o', pch=19)


