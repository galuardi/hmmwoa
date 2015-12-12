hmmsmooth <-
function(allpost, fmat, ks, sigma){
time1 = Sys.time()	
	allpost[is.nan(allpost)] = 0
	allpost2 = allpost*0
	post = allpost[,,dim(allpost)[3]]
	#post[is.nan(post)] = 0
	# f = gausskern(ks, sigma)
	# S <- filter2(post,f)
	# S <- S@.Data
	# post = S
	# post = normalise(post)
	allpost2[,,dim(allpost)[3]] = post
	allpost2[,,1] = allpost[,,1]

	#time1 = Sys.time()

for(i in rev(2:dim(allpost)[3])){#

if(fmat$behav[i]==1){
		sig = sigma[1]; kern = ks[1];
	 }else{
		sig = sigma[2]; kern = ks[2]
	 }
	 
   print(i)
   post[is.nan(post)] = 0
   #n = ks+(datediff[1]-1)*14
   # f = makeBrush(ks, shape='disc', step=FALSE) 
   f = gausskern(kern, sig)
   f = f/sum(f)
   rat = filter2(allpost[,,i-1],f)
   ratio = (allpost2[,,i]/(rat@.Data+1e-20))
   S <- filter2(ratio,f)
   #S <- filter2(allpost[,,i-1],f)
   S <- S@.Data
   #allS[,,i] = S
   post = S*allpost[,,i-1]
   allpost2[,,i-1] = normalise(post)
}
print(Sys.time()-time1)
allpost2
}
