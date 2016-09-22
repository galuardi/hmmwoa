# use dnorm to geenrate a likelihood
# use focal to get a sd field
# likint3 <- function(w, wsd, minT, maxT){
liksrss <- function(obs, srss, srsd){
  # midT = (maxT + minT) / 2
  # Tsd = (maxT - minT) / 4
  #d = obs-srss
  # widx = w >= minT & w <= maxT & !is.na(w)
  sdf = data.frame(sr = as.vector(srss), srsd = as.vector(srsd))
  sdf$srsd[is.na(sdf$srsd)] = 0
  # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
  # wint = apply(wdf, 1, function(x) integrate(dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd * 2)$value) 
  res = dnorm(obs, sdf$sr, sdf$srsd)
  srssout = srss
  values(srssout) = res
  # w = w * 0
  # w[widx] = wint
  # w
  srssout
} 