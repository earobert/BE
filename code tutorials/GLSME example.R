# GLSME example (I think the bootstrap tutorial is better)

library(GLSME)
library(mvSLOUCH)
library(ape)
library(ouch)
n<-5 ## number of species
apetree<-rtree(n)
phyltree<-ape2ouch(apetree) ##mvslouch requires ouch format
### Correct the names of the internal node labels.
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-
  as.character(1:(phyltree@nnodes-phyltree@nterm))
### Define Brownian motion parameters to be able to simulate data under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=2,ncol=1),Sxx=rbind(c(1,0),c(0.2,1)))
### Now simulate the data and remove the values corresponding to the internal nodes.
xydata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
xydata<-xydata[(nrow(xydata)-n+1):nrow(xydata),]
x<-xydata[,1]
y<-xydata[,2]
yerror<-diag((rnorm(n,mean=0,sd=0.1))^2) #create error matrix
y<-rmvnorm(1,mean=y,sigma=yerror)[1,]
xerror<-diag((rnorm(n,mean=0,sd=0.1))^2) #create error matrix
x<-rmvnorm(1,mean=x,sigma=xerror)[1,]
GLSME(y=y, CenterPredictor=TRUE, D=cbind(rep(1, n), x), Vt=vcv(apetree),
      Ve=yerror, Vd=list("F",vcv(apetree)), Vu=list("F", xerror))
      
plot(x,y)


vignette("using-lsmeans", package="lsmeans")
