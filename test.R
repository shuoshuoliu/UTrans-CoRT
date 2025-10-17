library(glmnet); library(MASS)

##################################################################################
#################     Implementation        ###################
##################################################################################
# Data Structure. 
#y: vector of response; 
#x1: matrix of covariates with nuisance parameter; 
#x2: matrix of covariates with testing parameter
# p1: dimension of nuisance parameter
# estimation based on model in the null hypothesis 

detect_test=function(x1,x2,y,n0,p1,family){
  X=cbind(x2,x1)
  dim(X)
  nn=nrow(x1)
  n1=nn-n0
  xx1=tail(x1,n=n0)
  yy1=tail(y,n=n0)
  if (family=="gaussian"){
    #beta_K=coef(cv.ncvreg(x1,y,penalty="SCAD",family="gaussian"))
    beta_K=coef(cv.ncvreg(X,y,penalty="SCAD",family="gaussian"))
    beta_K=c(beta_K[1],tail(beta_K,p1))
    u0=beta_K[1]+head(x1,n1)%*%beta_K[-1]
    # beta_K= coef(cv.ncpen(y.vec=y,x.mat=x1,intercept=FALSE,family="gaussian", penalty="scad"))$beta
    # u0=head(x1,n1)%*%beta_K
    
    un1=0; rn1=0;
    for (II in 1:n1){
      for (JJ in 1:n1){
        if (II != JJ){
          un1=un1+(y[II]-u0[II])*(y[JJ]-u0[JJ])*(t(x2[II,])%*%x2[JJ,])     # x2:   matrix of covariates with testing parameter
          rn1=rn1+(y[II]-u0[II])^2*((y[JJ]-u0[JJ])^2)*((t(x2[II,])%*%x2[JJ,])^2)
        }
      }
    }
    un1=un1/n1; rn1=rn1/(n1*(n1-1));
    test_stat=abs(un1)/sqrt(2*rn1);   # test_stat is the calculated test statistics
    recj=0
    if (test_stat>qnorm(0.995)){recj=1}  # when recj=0: fail to reject H0 (transferable); when recj=1: reject H0
  } else if (family=="binomial"){
    # beta_K=coef(cv.ncvreg(x1,y,penalty="SCAD",family="binomial"))
    # pii=exp(beta_K[1]+head(x1,n1)%*%beta_K[-1])/(1+exp(beta_K[1]+head(x1,n1)%*%beta_K[-1]))
    beta_K= coef(cv.ncpen(y.vec=y,x.mat=x1,intercept=FALSE,family="binomial", penalty="scad"))$beta
    pii=exp(head(x1,n1)%*%beta_K)/(1+exp(head(x1,n1)%*%beta_K))
    
    un1=0; rn1=0;
    for (II in 1:n1){
      for (JJ in 1:n1){
        if(II!=JJ){
          un1=un1+(y[II]-pii[II])*(y[JJ]-pii[JJ])*(t(x2[II,])%*%x2[JJ,])  # x2: matrix of covariates with testing parameter
          rn1=rn1+(y[II]-pii[II])^2*((y[JJ]-pii[JJ])^2)*((t(x2[II,])%*%x2[JJ,])^2)
        }
      }
    }
    un1=un1/n1; rn1=rn1/(n1*(n1-1)); 
    test_stat=abs(un1)/sqrt(2*rn1);   # test_stat is the calculated test statistics
    recj=0
    if (test_stat>qnorm(0.995)){recj=1}
  } else if (family=="poisson"){
    # beta_K=coef(cv.ncvreg(x1,y,penalty="SCAD",family="poisson"))
    # u0=exp(beta_K[1]+head(x1,n1)%*%beta_K[-1])
    beta_K= coef(cv.ncpen(y.vec=y,x.mat=x1,intercept=FALSE,family="poisson", penalty="scad"))$beta
    u0=exp(head(x1,n1)%*%beta_K)
    
    un1=0; rn1=0;
    for (II in 1:n1){
      for (JJ in 1:n1){
        if (II != JJ){
          un1=un1+(y[II]-u0[II])*(y[JJ]-u0[JJ])*(t(x2[II,])%*%x2[JJ,])     # x2:   matrix of covariates with testing parameter
          rn1=rn1+(y[II]-u0[II])^2*((y[JJ]-u0[JJ])^2)*((t(x2[II,])%*%x2[JJ,])^2)
        }
      }
    }
    un1=un1/n1; rn1=rn1/(n1*(n1-1));
    test_stat=abs(un1)/sqrt(2*rn1);   # test_stat is the calculated test statistics
    recj=0
    if (test_stat>qnorm(0.995)){recj=1}  # when recj=0: fail to reject H0; when recj=1: reject H0
  }
  return(recj)
}

# find_mode <- function(x) {
#   u <- unique(x)
#   tab <- tabulate(match(x, u))
#   u[tab == max(tab)]
# }
# 
# detect_test=function(x1,x2,y,n0,p1,family){
#   ans=replicate(4,detect_test0(x1,x2,y,n0,p1,family))
#   recj=find_mode(ans)
#   return(recj)
# }