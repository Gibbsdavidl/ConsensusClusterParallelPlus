### R code from vignette source 'ConsensusClusterPlus.Rnw'

###################################################
### code chunk number 1: ConsensusClusterPlus.Rnw:37-41
###################################################
library(ALL)
data(ALL)
d=exprs(ALL)
d[1:5,1:5]


###################################################
### code chunk number 2: ConsensusClusterPlus.Rnw:48-50
###################################################
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:5000],]


###################################################
### code chunk number 3: ConsensusClusterPlus.Rnw:55-56
###################################################
d = sweep(d,1, apply(d,1,median,na.rm=T))


###################################################
### code chunk number 4: ConsensusClusterPlus.Rnw:70-74
###################################################
library(ConsensusClusterPlus)
title=tempdir()
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")


###################################################
### code chunk number 5: ConsensusClusterPlus.Rnw:77-79
###################################################
cat(sprintf("\\graphicspath{{%s}}", paste(gsub("[\\]","/",title),"/",sep="")))
cat("\n")


###################################################
### code chunk number 6: ConsensusClusterPlus.Rnw:84-96
###################################################
#consensusMatrix - the consensus matrix.  
#For .example, the top five rows and columns of results for k=2:
results[[2]][["consensusMatrix"]][1:5,1:5]

#consensusTree - hclust object 
results[[2]][["consensusTree"]]

#consensusClass - the sample classifications
results[[2]][["consensusClass"]][1:5]

#ml - consensus matrix result
#clrs - colors for cluster  


###################################################
### code chunk number 7: ConsensusClusterPlus.Rnw:104-105
###################################################
icl = calcICL(results,title=title,plot="png")


###################################################
### code chunk number 8: ConsensusClusterPlus.Rnw:109-110
###################################################
icl[["clusterConsensus"]]


###################################################
### code chunk number 9: ConsensusClusterPlus.Rnw:113-114
###################################################
icl[["itemConsensus"]][1:5,]


###################################################
### code chunk number 10: ConsensusClusterPlus.Rnw:125-126
###################################################
cat("\\includegraphics[width=60mm]{consensus001.png}",sep="")


###################################################
### code chunk number 11: ConsensusClusterPlus.Rnw:135-137
###################################################
cat("\\includegraphics[width=60mm]{consensus002.png}",sep="")
cat("\\includegraphics[width=60mm]{consensus003.png}",sep="")


###################################################
### code chunk number 12: ConsensusClusterPlus.Rnw:140-142
###################################################
cat("\\includegraphics[width=60mm]{consensus004.png}",sep="")
cat("\\includegraphics[width=60mm]{consensus005.png}",sep="")


###################################################
### code chunk number 13: ConsensusClusterPlus.Rnw:148-149
###################################################
cat("\\includegraphics[width=60mm]{consensus007.png}",sep="")


###################################################
### code chunk number 14: ConsensusClusterPlus.Rnw:156-157
###################################################
cat("\\includegraphics[width=60mm]{consensus008.png}",sep="")


###################################################
### code chunk number 15: ConsensusClusterPlus.Rnw:167-168
###################################################
cat("\\includegraphics[width=60mm]{consensus009.png}",sep="")


###################################################
### code chunk number 16: ConsensusClusterPlus.Rnw:178-179
###################################################
cat("\\includegraphics[width=60mm]{icl003.png}",sep="")


###################################################
### code chunk number 17: ConsensusClusterPlus.Rnw:188-189
###################################################
cat("\\includegraphics[width=60mm]{icl001.png}",sep="")


###################################################
### code chunk number 18: ConsensusClusterPlus.Rnw:203-206
###################################################
#example of providing a custom distance matrix as input:
#dt = as.dist(1-cor(d,method="pearson"))
#ConsensusClusterPlus(dt,maxK=4,reps=100,pItem=0.8,pFeature=1,title="example2",distance="pearson",clusterAlg="hc")


###################################################
### code chunk number 19: ConsensusClusterPlus.Rnw:209-212
###################################################
#example of providing a custom distance function:
#myDistFunc = function(x){ dist(x,method="manhattan")}
#ConsensusClusterPlus(d,maxK=4,reps=100,pItem=0.8,pFeature=1,title="example3",distance="myDistFunc",clusterAlg="pam")


###################################################
### code chunk number 20: ConsensusClusterPlus.Rnw:216-223
###################################################
#library(cluster)
#dianaHook = function(this_dist,k){
  #tmp = diana(this_dist,diss=TRUE)
  #assignment = cutree(tmp,k)
  #return(assignment)  
#}
#ConsensusClusterPlus(d,clusterAlg="dianaHook",distance="pearson",...)


