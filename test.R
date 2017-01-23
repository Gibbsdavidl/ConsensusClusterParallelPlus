
library(ConsensusClusterParallelPlus)
library(ggplot2)

m1 <- data.frame(X=rnorm(2000), Y=rnorm(100))
m2 <- data.frame(X=rnorm(n=2000, 10, 2), Y=rnorm(n=100,10,2))
m3 <- data.frame(X=rnorm(n=2000, -2, 2), Y=rnorm(n=100,5,2))
m <- rbind(m1,m2,m3)

print(system.time(res0 <- ConsensusClusterPlus(d=t(m), maxK=5, distance="euclidean", innerLinkage="ward.D2", reps=5000)))

print(system.time(resTest <- ConsensusClusterParallelPlus(d=t(m), maxK=5, cores=4, distance="euclidean", innerLinkage="ward.D2", reps=5000)))


#mplot <- cbind(data.frame(Class=c(rep(1,100),rep(2,100),rep(3,100))), m)
#qplot(data=mplot, x=X, y=Y, col=Class)

#system.time(res0 <- ConsensusClusterPlus(d=t(m), maxK=5, distance="euclidean", innerLinkage="ward.D2", reps=5000))
#mplot0 <- cbind(data.frame(Class=res0[[3]]$consensusClass), m)
#qplot(data=mplot0, x=X, y=Y, col=Class)

#system.time(resTest <- ConsensusClusterParallelPlus(d=t(m), maxK=5, cores=4, distance="euclidean", innerLinkage="ward.D2", reps=5000))
#mplotTest <- cbind(data.frame(Class=resTest[[3]]$consensusClass), m)
#qplot(data=mplotTest, x=X, y=Y, col=Class)
