library(ncvreg)
library(glmnet)
library(foreach)
library(parallel)
cl <- makeCluster(4)
source("test.R")

cal_auc=function(x,y){
  #x:predicted, y:true
  n=length(x)
  
  if (sum(x)==n|sum(x)==0){
    x[1]=abs(x[1]-1)
  }
  
  res=auc(x,y, levels = c(0, 1), direction = "<")
  return(res)
}

get.source=function(D.source){
  #D.source: 1st col=source name; 2nd=y; others=x
  Data.s=NULL
  D.s=split(D.source, f = D.source[,1])
  
  source <- sapply(1:length(D.s), function(k){
    dt=D.s[[k]]
    x <- sapply(dt[,3:ncol(dt)], as.numeric)
    y <- as.numeric(dt[,2])
    list(x = x, y = y)
  }, simplify = FALSE)
  
  return(source)
}

allComb=function(x, k){
  dx = data.frame()
  dy = c()
  
  for (i in 1:k){
    out1=x[[i]]$x
    out2=x[[i]]$y
    dx=rbind(dx,out1)
    dy=c(dy,out2)
  }
  return(list(x=dx,y=dy))
}

propComb=function(source,target,family){ # stack all
  #n: sample size of source data
  transfer.source.id=1:length(source)
  
  X <- as.matrix(foreach(j = transfer.source.id, .combine = "rbind") %do% {
    source[[j]]$x})
  X=rbind(X,target$x)
    Y <- foreach(j = transfer.source.id, .combine = "c") %do% {
      source[[j]]$y}
    Y=c(Y,target$y)

  ans=list(X=X,Y=Y)
  return(ans)
}

propComb2=function(x,K,nv){
  #nv: vector containing sample sizes of source data 
  n.sum=c(0,cumsum(nv))
  dx = NULL
  
  for (i in 1:K){
    aa=matrix(0,nrow(x),ncol(x))
    aa[(n.sum[i]+1):n.sum[i+1],]=x[(n.sum[i]+1):n.sum[i+1],]
    dx=cbind(dx,aa)
  }
  
  final=cbind(dx,x)
  return(final)
}

simu <- function(family = c("gaussian", "binomial", "poisson"), s.vector,type = c("all", "source", "target"), 
                 cov.type = 1, h, K, n.target = 200, n.source = rep(200, K), 
                 s = length(s.vector), p = 500, Ka=K,sparsity="s",sigma=0.5) {
  #sigma: strength in the AR covariance
  family <- match.arg(family)
  target <- NULL
  source <- NULL
  
  Sigma <- outer(1:p, 1:p, function(x,y){
    sigma^(abs(x-y))
  })
  
  type <- match.arg(type)
  sig.strength <- 0.5
  
  if (family == "gaussian" || family == "binomial") {
    if(type == "all" || type == "target") {
      wk <- c(s.vector, rep(0, p-s))
      R <- chol(Sigma)
      target <- list(x = NULL, y = NULL)
      
      target$x <- matrix(rnorm(n.target*p), nrow = n.target) %*% R
      
      if (family == "gaussian") {
        target$y <- as.numeric(target$x %*% wk + rnorm(n.target))
      } else if (family == "binomial") {
        pr <- 1/(1+exp(-target$x %*% wk))
        target$y <- sapply(1:n.target, function(i){sample(0:1, size = 1, prob = c(1-pr[i], pr[i]))})
      }
    }
    
    if(type == "all" || type == "source") {
      if (cov.type == 1) {
        eps <- rnorm(p, sd = 0.3)
        Sigma <- Sigma + eps %*% t(eps)
        R <- chol(Sigma)
      }
      
      source <- sapply(1:K, function(k){
        if (k <= Ka){
          if (sparsity=="p"){
            wk <- c(s.vector, rep(0, p-s)) + h/p*sample(c(-1,1), size = p, replace = TRUE)
          }else if (sparsity=="s"){
            wk <- c(s.vector+ h/p*sample(c(-1,1), size = s, replace = TRUE), rep(0, p-s)) 
          }
        } else {
          sig.index <- c(s+1:s, sample((2*s+1):p, s))
          wk <- rep(0, p)
          wk[sig.index] <- 0.5#c(s.vector,s.vector)
          wk <- wk + 2*h/p*sample(c(-1,1), size = p, replace = TRUE)
        }
        
        if (cov.type == 1) { # normal
          x <- matrix(rnorm(n.source[k]*p), nrow = n.source[k]) %*% R
        } else if (cov.type == 2) { # normal with CS
          temp=matrix(0.5,p,p)
          diag(temp)=1
          x=mvrnorm(n.source[k],rep(0,p),temp)%*% chol(temp)
        } else if (cov.type == 3) { # t dist
          x <- matrix(rt(n.source[k]*p, df = 4), nrow = n.source[k])%*% R
        } else if (cov.type == 4) { # beta dist
          myMix <- UnivarMixingDistribution(Norm(mean=2*sqrt(5)/5, sd=1/sqrt(5)), 
                                            Norm(mean=-2*sqrt(5)/5, sd=1/sqrt(5)),
                                            mixCoeff=c(0.5,0.5))
          rmyMix <- r(myMix)
          x=matrix(rmyMix(n.source[k]*p), nrow = n.source[k])%*% R
        }
        
        if (family == "gaussian") {
          y <- as.numeric(0.5*I(k > Ka) + x %*% wk + rnorm(n.source[k]))
        } else if (family == "binomial") {
          pr <- 1/(1+exp(-0.5*I(k > Ka)-x %*% wk))
          y <- sapply(1:n.source[k], function(i){
            sample(0:1, size = 1, prob = c(1-pr[i], pr[i]))
          })
        }
        list(x = x, y = y)
      }, simplify = FALSE)
    }
    
  } else if (family=="poisson"){ # model == "poisson
    if(type == "all" || type == "target") {
      R <- chol(Sigma)
      wk <- c(s.vector, rep(0, p-s))
      # 
      target$x <- matrix(rnorm(n.target*p), nrow = n.target) %*% R
      target$x[target$x > 0.5] <- 0.5
      target$x[target$x < -0.5] <- -0.5
      lambda <- as.numeric(exp(target$x %*% wk))
      target$y <- rpois(n.target, lambda)
    }
    
    if(type == "all" || type == "source") {
      if (cov.type == 1) {
        eps <- rnorm(p, sd = 0.3)
        Sigma <- Sigma + eps %*% t(eps)
        
        R <- chol(Sigma)
      }
      
      source <- sapply(1:K, function(k){
        if (k <= Ka){
          wk <- c(s.vector, rep(0, p-s)) + h/p*sample(c(-1,1), size = p, replace = TRUE)
        } else {
          sig.index <- c(s+1:s, sample((2*s+1):p, s))
          wk <- rep(0, p)
          wk[sig.index] <- 0.5#c(s.vector,s.vector)
          wk <- wk + 2*h/p*sample(c(-1,1), size = p, replace = TRUE)
        }
        if (cov.type == 1) {
          x <- matrix(rnorm(n.source[k]*p), nrow = n.source[k]) %*% R
        } else if (cov.type == 2) {
          x <- matrix(rt(n.source[k]*p, df = 4), nrow = n.source[k])
        }
        x[x > 0.5] <- 0.5
        x[x < -0.5] <- -0.5
        lambda <- as.numeric(exp(0.5*I(k > Ka) + x %*% wk))
        y <- rpois(n.source[k], lambda)
        
        list(x = x, y = y)
      }, simplify = FALSE)
    }
  } 
  if (type == "all") {
    return(list(target = target, source = source))
  } else if (type == "target") {
    return(list(target = target))
  } else {
    return(list(source = source))
  }
}


library(caret)
source_d <- function(target, source, type="lasso",valid.nfolds = 9, mode=c("test","data"),
                     family = c("gaussian", "binomial", "poisson"),...) {
  if (mode=="test"){ #detect by testing
    
    ans=rep(NA,length(source))
    for (k in 1:length(source)){
      x1=rbind(source[[k]]$x,target$x) #nuisance
      x2=rbind(source[[k]]$x,matrix(0,nrow(target$x),ncol(target$x)))
      y=c(source[[k]]$y,target$y)
      ans[k]=detect_test(x1,x2,y=y,nrow(target$x),ncol(target$x),family=family)
    }
    transfer.source.id=which(ans==0)
    cat("Detection by test:", transfer.source.id)
    cat("\n")
    obj=list(transfer.source.id = transfer.source.id)
  } else if (mode=="data"){
    n=length(target$y)
    
    if(nrow(target$x)<50){
      valid.nfolds=round(nrow(target$x)/5)
      valid.nfolds=ifelse(valid.nfolds%%2==0,valid.nfolds-1,valid.nfolds)
    }
    
    family <- match.arg(family)
    folds <- createFolds(target$y, valid.nfolds)
    penalty="lasso"
    
    loss.cv <- t(sapply(1:valid.nfolds, function(i){
      source.loss <- sapply(1:length(source), function(k){
        wa=coef(cv.ncvreg(as.matrix(rbind(target$x[-folds[[i]], , drop = F],source[[k]]$x)),
                          c(target$y[-folds[[i]]],source[[k]]$y),penalty=penalty,family=family,cluster=cl))
        
        loss(wa, as.matrix(target$x[folds[[i]], , drop = F]), target$y[folds[[i]]], family=family)
      })
      
      wa.target=coef(cv.ncvreg(as.matrix(target$x[-folds[[i]], , drop = F]),target$y[-folds[[i]]],penalty=penalty,family=family,cluster=cl))
      target.loss <- loss(wa.target, as.matrix(target$x[folds[[i]], , drop = F]), target$y[folds[[i]]], family=family)
      c(source.loss, target.loss)
    }))
    
    source.loss <- colMeans(loss.cv)[1:(ncol(loss.cv)-1)]
    target.valid.loss <- colMeans(loss.cv)[ncol(loss.cv)]
    
    n_col=ncol(loss.cv)-1
    transfer.source.id=rep(NA,n_col)
    for (j in 1:n_col){
      transfer.source.id[j]=sum(loss.cv[,j]<=loss.cv[,(n_col+1)])
    }
    transfer.source.id <- which(transfer.source.id >=(valid.nfolds+1)/2)
    cat("Detection by data:", transfer.source.id)
    cat("\n")
    obj <- list(transfer.source.id = transfer.source.id, source.loss = source.loss, target.valid.loss = target.valid.loss)
  }
  return(obj)
}

loss <- function(wa, x.valid, y.valid, family, tau = NULL) {
  if (family == "gaussian") {
    mean((y.valid - x.valid %*% wa[-1] - wa[1])^2)
  } else if (family == "binomial"){
    xb <- x.valid %*% wa[-1] + wa[1]
    as.numeric(- t(y.valid) %*% xb + sum(log(1+exp(xb))))/length(y.valid)
  } else if (family == "poisson") {
    xb <- x.valid %*% wa[-1] + wa[1]
    as.numeric(- t(y.valid) %*% xb + sum(exp(xb)) + sum(lgamma(y.valid+1)))/length(y.valid)
  } else { # family == "huber
    xb <- x.valid %*% wa[-1] + wa[1]
    res <- y.valid - x.valid %*% wa[-1] - wa[1]
    mean((res^2)/2*I(abs(res)<=tau) + (tau*abs(res)-(tau^2)/2)*I(abs(res)>tau))
  }
}

UTrans=function(Data,transfer.id=NULL,family = c("gaussian", "binomial", "poisson"),mode=c("test","data"),
                type=c("lasso","scad"),valid.nfolds){
  target=Data$target
  source=Data$source
  np=ncol(target$x)
  
  if (is.null(transfer.id)){
    res=source_d(target=target, source=source,family=family,mode=mode,valid.nfolds=valid.nfolds)
    transfer.source.id=res$transfer.source.id
  }else{
    transfer.source.id=transfer.id
  }
  
  kk=length(transfer.source.id)
  
  if (kk>0){ # when some sources are transferable
    X <- as.matrix(foreach(j = transfer.source.id, .combine = "rbind") %do% {
      source[[j]]$x})
    X=rbind(X,target$x)
    
    Y <- foreach(j = transfer.source.id, .combine = "c") %do% {
      source[[j]]$y}
    Y=c(Y,target$y)
    
    nk=rep(NA,kk)
    for (i in 1:kk){
      ind=transfer.source.id[i]
      nk[i]=length(source[[ind]]$y)
    }
    XX2=propComb2(X,kk,nk)
  } else if (kk==0){
    X=target$x
    XX2=target$x
    Y=target$y
  }
  
  if (type=="lasso"){
    lasso.fit=cv.ncvreg(XX2,Y,penalty="lasso",family=family)
    res=list(transfer.source.id=transfer.source.id,dx=X,dy=Y,x=XX2,b1=c(coef(lasso.fit)[1],tail(coef(lasso.fit),np)))
  } else if (type=="scad"){
    scad.fit=cv.ncvreg(XX2,Y,penalty="SCAD",family=family)
    res=list(transfer.source.id=transfer.source.id,dx=X,dy=Y,x=XX2,b2=c(coef(scad.fit)[1],tail(coef(scad.fit),np)))
  } else {
    lasso.fit=cv.ncvreg(XX2,Y,penalty="lasso",family=family)
    scad.fit=cv.ncvreg(XX2,Y,penalty="SCAD",family=family)
    #b1:lasso   b2:scad
    res=list(b1=c(coef(lasso.fit)[1],tail(coef(lasso.fit),np)),b2=c(coef(scad.fit)[1],tail(coef(scad.fit),np)),
             transfer.source.id=transfer.source.id,dx=X,dy=Y,x=XX2)
  }
  return(res)
}













AUTrans=function(D.training,family = c("gaussian", "binomial", "poisson"),type="both"){
  target=D.training$target
  source=D.training$source
  np=ncol(target$x)
  
  transfer.source.id=1:length(source)
  kk=length(transfer.source.id)
  
  X <- as.matrix(foreach(j = transfer.source.id, .combine = "rbind") %do% {
    source[[j]]$x})
  X=rbind(X,target$x)
  Y <- foreach(j = transfer.source.id, .combine = "c") %do% {
    source[[j]]$y}
  Y=c(Y,target$y)
  
  nk=rep(NA,kk)
  for (i in 1:kk){
    ind=transfer.source.id[i]
    nk[i]=length(source[[ind]]$y)
  }
  XX2=propComb2(X,kk,nk)
  
  if (type=="lasso"){
    if (family=="poisson"){ # ncvreg has some error in poisson regression
      lasso.fit <- cv.glmnet(x=XX2, y =Y, family = family, alpha = 1,parallel=TRUE)
    } else {
      lasso.fit=cv.ncvreg(XX2,Y,penalty="lasso",family=family,cluster=cl)
    }
    res=list(b1=c(coef(lasso.fit)[1],tail(coef(lasso.fit),np)),transfer.source.id=transfer.source.id,dx=X,dy=Y,x=XX2)
  } else if (type=="scad"){
    scad.fit=cv.ncvreg(XX2,Y,penalty="SCAD",family=family,cluster=cl)
    res=list(b2=c(coef(scad.fit)[1],tail(coef(scad.fit),np)),transfer.source.id=transfer.source.id,dx=X,dy=Y,x=XX2)
  } else{
    lasso.fit=cv.ncvreg(XX2,Y,penalty="lasso",family=family,cluster=cl)
    scad.fit=cv.ncvreg(XX2,Y,penalty="SCAD",family=family,cluster=cl)
    res=list(b1=c(coef(lasso.fit)[1],tail(coef(lasso.fit),np)),b2=c(coef(scad.fit)[1],tail(coef(scad.fit),np)),
             transfer.source.id=transfer.source.id,dx=X,dy=Y,x=XX2)
  }
  return(res)
}
