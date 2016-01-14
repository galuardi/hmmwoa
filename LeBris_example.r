# build an example with some real data
# objective here is to properly characterize error in the PDT
# matching routine between depth-temp data from the tag
# and WOA data.

library(locfit)
library(fields)
library(plyr)

# this is real data
minT <- c(21.4,21.6,18.2,17.4,15.6,14.2,11.8,8.4)
maxT <- c(23,22.4,18.6,17.8,16.6,15.4,12.6,8.8)
deps <- c(0,64,320,456,576,640,752,952)
midT <- (maxT + minT) / 2

stdDepth <- c(0.0,2.5,7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5,
              62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5,112.5,137.5,162.5,187.5,212.5,
              237.5,262.5,287.5,312.5,337.5,362.5,387.5,412.5,437.5,462.5,487.5,525.0,575.0,
              625.0,675.0,725.0,775.0,825.0,875.0,925.0,975.0, 1025.0, 1075.0, 1125.0, 1175.0,
              1225.0, 1275.0, 1325.0, 1375.0, 1425.0, 1475.0)

# do the regression
fit <- locfit(midT ~ deps)
n = length(deps)
# find the standard depth levels that correspond to the depths the tag data is measured at
depIdx <- findInterval(deps, stdDepth)
woaDep <- stdDepth[depIdx] 
# make predictions based on the regression model earlier for the temperature at standard WOA depth levels
pred = predict(fit, newdata=woaDep, se=T, get.data=T)

# see the regression
plot(fit, get.data=TRUE)
# and add SE based on model predictions at new depth levels
lines(woaDep, pred$fit+pred$se.fit*sqrt(n), col=2, lty=2)
lines(woaDep, pred$fit-pred$se.fit*sqrt(n), col=2, lty=2)
#points(woaDep, pred$fit+pred$se.fit*sqrt(n), col=2)
#points(woaDep, pred$fit-pred$se.fit*sqrt(n), col=2)

# fake temperature surface
woa1 = matrix(1:100/3, 10,10)
woa2 = matrix(1:100/3,10,10)
woa = rbind(woa1[1:5,],woa2[1:5,])
image.plot(woa,xlab='fake longitude',ylab='fake latitude', main='ocean temperature surface')

# put everything in one spot
out <- cbind(woaDep,fitT=pred$fit,sd_lwr=pred$fit-pred$se.fit*sqrt(n),sd_upr=pred$fit+pred$se.fit*sqrt(n))

# ** now how do we implement LeBris's likelihood routine via integration from sd_lwr to sd_upr?

#-------------------------------------------------------------------------# 
# like this! 
#-------------------------------------------------------------------------# 
# do the regression at min and max
fit.low <- locfit(minT ~ deps)
fit.high <- locfit(maxT ~ deps)
n = length(deps)
# find the standard depth levels that correspond to the depths the tag data is measured at
# depIdx <- findInterval(deps, stdDepth)

# use the which.min
depIdx = apply(as.data.frame(deps), 1, FUN=function(x) which.min((x-stdDepth)^2))
woaDep <- stdDepth[depIdx] 

# make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
pred = predict(fit, newdata=stdDepth, se=T, get.data=T)
pred.low = predict(fit.low, newdata=stdDepth, se=T, get.data=T)
pred.high = predict(fit.high, newdata=stdDepth, se=T, get.data=T)
plot(fit, get.data=TRUE, lwd=2)

# and add SE based on model predictions at new depth levels
lines(stdDepth, pred$fit+pred$se.fit*sqrt(n), col=2, lty=2)
lines(stdDepth, pred$fit-pred$se.fit*sqrt(n), col=2, lty=2)

lines(fit.low, col = 4, lwd=2)           
lines(fit.high, col = 4, lwd=2)

lines(stdDepth, pred.high$fit+pred.high$se.fit*sqrt(n), col=4, lty=2)
lines(stdDepth, pred.low$fit-pred.low$se.fit*sqrt(n), col=4, lty=2)

n =length(stdDepth)

# data frame for next step
df = data.frame(low=pred.low$fit[depIdx]-pred.low$se.fit[depIdx]*sqrt(n), high=pred.high$fit[depIdx]+pred.high$se.fit[depIdx]*sqrt(n), row.names = deps)

#---------------------------------------------------------------------------------------#
# WOA sd.. 
#---------------------------------------------------------------------------------------#
# fixed sd of the grid
woasd = .7

# slow derivation of sd by surrounding cells
for(i in 1:nrow(woa)){
  for(j in 1:ncol(woa)){
        
  }
}

#---------------------------------------------------------------------------------------#
# integration would be something like this for each WOA layer..
# note to self: read up on integration.. 
#---------------------------------------------------------------------------------------#

# I think a sequence might hae to be set between the min and max... I tried just adding dnorm(min) and dnorm(max) and it is not right
# for example:

# image.plot(dnorm(df[1,1], woa, .7)+dnorm(df[1,2], woa, .7))

# sequence
dt = seq(df[1,1], df[1,2], length=10)

# likelihood array for one depth
lik0 = aaply(dt, 1, .fun = function(x) dnorm(x, woa, .7))
int0 = apply(lik0, 2:3, sum)

# plot it
par(mfrow=c(1,2))
image.plot(woa)
image.plot(int0, col = terrain.colors(100))

    





