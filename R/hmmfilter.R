hmmfilter <-
function(sstL, fmat, ks=29, sigma=1){
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
