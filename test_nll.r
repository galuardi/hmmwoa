# par0 <- c(log(10), log(10), log(0.5), log(0.5), log(0.95/0.05), log(0.95/0.05))

par0 = c(40, 10, 10, 5, .707, .866)

fit <- nlm(get.nll.fun, par0, g, L, dt)
D1 <- exp(fit$estimate[1:2])
D2 <- exp(fit$estimate[3:4])
p <- 1/(1+exp(-fit$estimate[5:6])) 

## Try with larger Diffusion
t <- Sys.time()
par0 = c(100, 300, 25, 10, .707, .866) # transitions from Pedersen 2011
fit <- nlm(f = get.nll.fun, p = par0, g.mle, L.mle)
Sys.time() - t


## Try with larger Diffusion and logit transformed transition prob (.95)
t <- Sys.time()
par0 = c(100, 300, 25, 10, .95, .95) # transitions from Pedersen 2011
fit <- nlm(f = get.nll.fun, p = par0, g.mle, L.mle)
Sys.time() - t

# 
# > fit
# $minimum
# [1] 3857.8
# 
# $estimate
# [1]   40.00000 1745.00624   10.00000    1.22522    0.89162    0.71539
# 
# $gradient
# [1]  0.0000e+00 -8.9099e-07  0.0000e+00  8.2211e-04 -2.5525e-03  3.4684e-03
# 
# $code
# [1] 1
# 
# $iterations
# [1] 46
# 
# > par0
# [1] 40.000 10.000 10.000  5.000  0.707  0.866

par0 = c(40, 10, 10, 5, .707, .866)
p1 = p2 = exp(log(0.95/0.05))/(1+exp(log(0.95/0.05))) # inverse logit for trans prob terms
logpar = c(par0[1:4], p1, p2) 
t1 = Sys.time()
logfit <- nlm(f = get.nll.fun, p = logpar, g.mle, L.mle)
runtime = Sys.time()-t1
runtime


D1 <- exp(fit$estimate[1:2])
D2 <- exp(fit$estimate[3:4])
p <- 1/(1+exp(-fit$estimate[5:6]))

# use optim and logit transformed transition parameters
par0 = c(100/2.6, 300/2.6, 25/2.6, 10/2.6, .95, .95)
# par0 = c(40, 10, 10, 5, .707, .866)
p1 = p2 = exp(log(0.95/0.05))/(1+exp(log(0.95/0.05))) # inverse logit for trans prob terms
logpar = c(par0[1:4], p1, p2) 
ofit <- optim(par = c(logpar), get.nll.fun, hessian = F)
1/(1+exp(-ofit$par[5:6]))

par1 = c(abs(ofit$par[1:4])*2.6, 1/(1+exp(-ofit$par[5:6])))

        
