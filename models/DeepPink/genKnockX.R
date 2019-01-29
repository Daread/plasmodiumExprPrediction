
setwd("/media/yanglu/TOSHIBA/data/featuresAndResponseDataframeCSVs/2018_6_15_fiveFoldCSVs")

source('/media/yanglu/TOSHIBA/deepknockoff/src/RANK/gknockoff_all.R')
source('/media/yanglu/TOSHIBA/deepknockoff/src/RANK/isee_all.R')


library(knockoff)
require(graphics)
library(igraph)
library(MASS) 
library(Matrix)
library(pracma)

regfactor = "log"
npermu = 5
sis.use = 0
bia.cor = 0


fdr = 0.20

dataTypeList=c("LRnoMotifs","LRwithMotifs","LSnoMotifs","LSwithMotifs","LTnoMotifs","LTwithMotifs")


for (dataType in dataTypeList) {

xURL = paste(dataType,"/","X.csv",sep="")
yURL = paste(dataType,"/","Y.csv",sep="")

xKnockURL = paste(dataType,"/","X_knockoff.csv",sep="")
print(xURL)

X <- read.csv(xURL, header = FALSE, sep = ",")
print(nrow(X))
print(ncol(X))

Y <- read.csv(yURL, header = FALSE, sep = ",")
print(nrow(Y))
print(ncol(Y))
y=as.numeric(unlist(Y))


#result = knockoff.filter(X, y, fdr=fdr, knockoffs= create.second_order, statistic=stat.glmnet_lambdasmax)

#knockoffs = function(X) create.fixed(X, method='equi')
#result = knockoff.filter(X, y, knockoffs=knockoffs)

#Xnew=cbind(result$X,result$Xk)

obj = isee(X, regfactor, npermu, sis.use, bia.cor) 
Omega = obj$Omega.isee

Xnew = gknockoffX(X, Omega)

print(nrow(Xnew))
print(ncol(Xnew))
print('--------------------------------------------------------')
write.table(Xnew, file=xKnockURL,row.names=FALSE, col.names=FALSE, sep=",")

}







