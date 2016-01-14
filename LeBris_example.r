# build an example with some real data
# objective here is to properly characterize error in the PDT
# matching routine between depth-temp data from the tag
# and WOA data.

library(locfit)

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


