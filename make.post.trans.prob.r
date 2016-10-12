# random vector 
sv = c(1,1,1,2,1,2,2,2,2,2,1,1,1,2,1,2,1,2,1,1,2,1,1,1,1,2,2,1,2,1,2,1,2,2,2,2,2,1,2,1,2,2,2,2,2)

# difference it to find transitions
svd = rbind(sv, c(0,diff(sv)))

# plot it 
plot(svd[1,], typ = 'l', col = 2, lwd= 3)

# get a new transition probability matrix
r1 = rev(table(svd[2, svd[1,]==1])/sum(svd[1,]==1))
r2 = rev(table(svd[2, svd[1,]==2])/sum(svd[1,]==2))
P1 = matrix(rbind(r1, r2), 2,2)

