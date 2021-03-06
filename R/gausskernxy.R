gausskernxy <-
function(siz, xsigma, ysigma, muadv = 0){
x = 1:round(siz);
mu = c(mean(x), mean(x)) + muadv;
fx = (matrix(exp((-0.5*(x-mu[1])/xsigma)^2))/(sqrt(2*pi)*xsigma));

options(digits=5)
fx = exp(-.5*((x-mu[1])/xsigma)^2)/sqrt((2*pi)*xsigma)
fy = exp(-.5*((x-mu[2])/ysigma)^2)/sqrt((2*pi)*ysigma)

fx[!is.finite(fx)] = 0
#fy = (matrix(exp((-0.5*((x-mu[2])/sigma))^2))/(sqrt(2*pi)*sigma));
fy[!is.finite(fy)] = 0
kern = (fx%*%t(fy))
kern = kern/(sum(sum(kern,na.rm=T),na.rm=T))
kern[is.nan(kern)]=0
kern
}
