source("C:/Users/aper/Documents/PH.R")
###############################################################################################################################################################
#Test

obj <- PH(alpha=c(0.2,0.3,0.5),S=matrix(c(-1,0,0,0.2,-1,0,0,0.3,-2),nrow=3))
obj

y <- 1:100
weights <- dexp(1:100,1/10)



obj <- fit(obj,y,weights,100,every=10,scale=1)



obj2 <- ph(alpha=c(0.2,0.3,0.5),S=matrix(c(-1,0,0,0.2,-1,0,0,0.3,-2),nrow=3))
obj2 <- fit(obj2,y,weight=weights,stepsEM=100)


plot(y,weights/sum(weights))
curve(dens(obj,x),add=T,col='red')
curve(dens(obj2,x),add=T,col="green")

plot(obj)



obj3 <- PH(structure = "erlang généralisée",dimension=4)
y <- 1:500
weights <- dgamma(1:500,3,1/25)

obj3 <- fit(obj3,y,weights,1000,every=100)

plot(obj3)

obj4 <- pmle(obj3)
plot(obj4)








#benchmark              
library(microbenchmark)
bench <- microbenchmark(fit(obj,y,weights,stepsEM=10),fit(obj2,y,weight=weights,stepsEM=10,every=10))
bench
plot(bench)


