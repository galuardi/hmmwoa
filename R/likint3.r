likint3 <- function(w, wsd, minT, maxT){
  midT = (maxT+minT)/2
  Tsd = (maxT-minT)/4
  widx = w>=minT&w<=maxT&!is.na(w)
  widxv = as.vector(widx)
  wdf = data.frame(w = as.vector(w[widx]), wsd = as.vector(wsd[widx]))
  wdf$wsd[is.na(wdf$wsd)] = 0
  # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
  wint = apply(wdf, 1, function(x)integrate(dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd*2)$value) 
  w = w*0  
  w[widx] = wint
  w
} 