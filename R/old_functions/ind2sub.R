ind2sub <-
function(siz, ndx){
	#nout = max(nargout,1);
	nout = 2
	siz = as.double(siz);
	len = end(siz)[1]

	if (length(siz)<=nout){
	  siz = c(siz, numeric(nout-length(siz)));
	}else{
	  siz = c(siz[1:(nout-1)], prod(siz[nout:len]));
	}
	n = length(siz);
	k = c(1, cumprod(siz[1:(len-1)]));
	varargout=numeric(n)
	for(i in rev(1:n)){
	  vi = rem(ndx-1, k[i]) + 1;         
	  vj = (ndx - vi)/k[i] + 1; 
	  varargout[i] = vj; 
	  ndx = vi;     
	}
	varargout
}
