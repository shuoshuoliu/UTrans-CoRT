rm(list=ls())
set.seed(0, kind = "L'Ecuyer-CMRG")
library(glmtrans)
library(pROC)

source("help2.R")

np=500 #p in the paper
n=100 #target size
h=30 #maximum difference divided by p
nK=10

family="gaussian"

s.vector=c(rep(0.4,3),rep(0.5,3),-rep(0.6,4)) #true signal variables
s=length(s.vector)
wk <- c(s.vector, rep(0, np-s)) #true beta
b0=c(0,wk) #true beta plus intercept

itr=30
ans=matrix(NA,itr,2)
final=matrix(NA,nK,4)
ans2=ans
final2=final

for (k in 1:4){
  cv=2*k+3
  print(cv)
  for (i in 1:itr){
    # data generation: Ka is # of transferable source data
    D.training <- simu(family=family, s=s,s.vector=s.vector,type = "all",Ka=8,K =10,n.target =n,h=h,p =np)
    D.test <- simu(family=family, s=s,s.vector=s.vector,type = "target",Ka=8,K =10,n.target =n,h=h,p =np)
    
    # proposal
    # mode="data" is for CoRT, mode="test" is for UTrans
    result=UTrans(D.training,family = family,valid.nfolds=cv,mode="data",type="lasso")
    beta.lasso=result$b1
    
    ######error#######
    r3=sqrt(sum((beta.lasso - b0)^2))# prop lasso
    ans[i,]=c(1,r3)
    if (i>1){print(colMeans(ans[1:i,]))}
    
  }
  final[k,]=c(colMeans(ans),apply(ans,2,sd))
  print(final[k,])
}


