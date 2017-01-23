
library(ConsensusClusterParallelPlus)
library(ggplot2)

m1 <- data.frame(X=rnorm(100), Y=rnorm(100))
m2 <- data.frame(X=rnorm(n=100, 10, 2), Y=rnorm(n=100,10,2))
m3 <- data.frame(X=rnorm(n=100, -2, 2), Y=rnorm(n=100,5,2))
m <- rbind(m1,m2,m3)

mplot <- cbind(data.frame(Class=c(rep(1,100),rep(2,100),rep(3,100))), m)
qplot(data=mplot, x=X, y=Y, col=Class)


res0 <- ConsensusClusterPlus(d=t(m), maxK=5, distance="euclidean", innerLinkage="ward.D2")
mplot0 <- cbind(data.frame(Class=res0[[3]]$consensusClass), m)
qplot(data=mplot0, x=X, y=Y, col=Class)


resTest <- ConsensusClusterParallelPlus(d=t(m), maxK=5, cores=4, distance="euclidean", innerLinkage="ward.D2")
mplotTest <- cbind(data.frame(Class=resTest[[5]]$consensusClass), m)
qplot(data=mplotTest, x=X, y=Y, col=Class)
