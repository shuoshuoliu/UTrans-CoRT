rm(list=ls())
set.seed(0, kind = "L'Ecuyer-CMRG")
library(glmtrans)
library(pROC)

source("help2.R")

np=500 #p in the paper
n=100 #target size
h=10 #maximum difference divided by p
nK=10

family="binomial"

s.vector=c(0.5,-0.4,0.7,-0.3,0.8)
s=length(s.vector)
wk <- c(s.vector, rep(0, np-s)) #true beta
b0=c(0,wk) #true beta plus intercept

itr=50
ans=matrix(NA,itr,4)
final=matrix(NA,nK,8)


for (k in 5:nK){
  print(k)
  for (i in 1:itr){
    # data generation
    D.training <- simu(family=family, s=s,s.vector=s.vector,type = "all",K =k,n.target =n,h=h,p =np)
    D.test <- simu(family=family, s=s,s.vector=s.vector, type = "target",n.target =n,h=h,p =np)
    
    #transfer GLM-Trans
    fit.gaussian <- glmtrans(target = D.training$target, source = D.training$source,detection.info =FALSE,
                             transfer.source.id = "all",family=family)

    
    # lasso with target only
    library(glmnet)
    fit.lasso <- cv.glmnet(x = D.training$target$x, y = D.training$target$y,
                           family = family)

    
    # proposal
    result=AUTrans(D.training,family = family)
    beta.lasso=result$b1
    beta.scad=result$b2
    
    ######error#######
    r1=sqrt(sum((fit.gaussian$beta - b0)^2)) #glmtrans
    r2=sqrt(sum((coef(fit.lasso) - b0)^2)) # lasso target only
    r3=sqrt(sum((beta.lasso - b0)^2))# prop lasso
    r4=sqrt(sum((beta.scad - b0)^2)) # prop scad
    ans[i,]=c(r1,r2,r3,r4)
    if (i>1){print(colMeans(ans[1:i,]))}
  }
  final[k,]=c(colMeans(ans),apply(ans,2,sd))
  print(final[k,])
}

