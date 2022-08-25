library(glmnet)
library(foreach)
library(dplyr)
library(magrittr)

# # lasso with regularisation parameter selected by cv, rolling window. Lag = 4
# 
# # cross validation (b/w period T1 n T2)
#
# # Assuming the data size is n*240, the moving window size for estimation is 72 and for evaluation is 61.

y <- dat[, targetVar] %>%
  set_colnames("y")
lambdaChoises <- 10^(seq(0.5,-3,len=100)) # lambda choices, selection on CV

# up till above

X <- lag.xts(dat, h:(h+3)) # lag=4

predErr <-
  foreach(t = 1:72, .combine = "cbind", .inorder = F,.packages=c("magrittr", "glmnet","zoo")) %dopar% {
    fit <- glmnet(X[(4+h):(T1+t-1),],y[(4+h):(T1+t-1)], # (p+h):T1 instead of 1:T1 bc first p obs's are missing #?
                  lambda=lambdaChoises, family="gaussian", alpha=1, standardize=F, intercept=T,
                  thresh=1e-15, maxit=1e07) # choose smaller thresh if nr of nonzero coef exceeds winSize
    pred <- predict.glmnet(fit, coredata(X[(T1+t),]))
    as.numeric((pred - as.numeric(y[T1+t]))^2)
  } # endforeach

cvScore <- apply(predErr,1,mean) # MSE for each lambda in cross validation period
optLam <- lambdaChoises[which.min(cvScore)]


LASSOlambda[horizon,targetVar] <- optLam

# evaluatoin
eval1 <- foreach(t = 1:61, .packages=c("magrittr", "glmnet","zoo")) %dopar% { 
  fit <- glmnet(X[(T2+t-60):(T2+t-1),], y[(T2+t-60):(T2+t-1)], lambda = optLam,
                family = "gaussian", alpha = 1, standardize = F, intercept=T,
                thresh=1e-15, maxit = 1e07)
  pred <- predict.glmnet(fit, coredata(X[(T2+t),]))
  err <- as.numeric((pred - y[T2+t])^2)
  coefs <- as.numeric(fit$beta)
  list(err, coefs)
}
predErr1 <- unlist(sapply(eval1, function(foo) foo[1]))
coefTracker <- matrix(unlist(sapply(eval1, function(foo) foo[2])),
                      nrow=61, ncol=ncol(X), byrow=T) # 61x960 matrix

coefTracker[coefTracker == 0] <- 0
coefTracker[coefTracker != 0] <- 1 # 1 if param is selected (non-zero)

msfeLasso <- mean(predErr1)

LASSOnonzero[horizon,targetVar] <- sum(coefTracker)/61 

# number of non-zero coef
LASSOV <- round(LASSOnonzero[horizon,targetVar])

## find the all significant regressors
ave_coefTracker<-matrix(0,1,960)
sum_coefTracker<-matrix(0,1,960)
if(horizon == 1){
  Lasso_coef_names[[var]] = matrix(NA, nrow=5, ncol = 960, 
                                   dimnames = list(c(paste("h=",hChoises,sep="")), 
                                   seq(1,960,by=1)))
  LASSORF_importance[[var]] = matrix(NA, nrow=5, ncol = 960, 
                                     dimnames = list(c(paste("h=",hChoises,sep="")), 
                                                     seq(1,960,by=1)))
}


for(j in 1:960){

  sum_coefTracker[j]<-sum(coefTracker[,j])
  ave_coefTracker[j]<-sum(coefTracker[,j])/61
  
}

colnames(sum_coefTracker) = colnames(X)
colnames(ave_coefTracker) = colnames(X)

ave_t = as.matrix(ave_coefTracker) 


# bubble sort

list = c()
for(n in 1:LASSOV){
  Lasso_coef_names[[var]][horizon,n] <- names(ave_t[,which.max(ave_t)])
  list = append(list, names(ave_t[,which.max(ave_t)]))
  ave_t = ave_t[,-which.max(ave_t)]
  ave_t = t(as.matrix(ave_t))
}

Val<-merge.xts(y,X)
Val<-as_tibble(Val)

Val1<-select(Val, y,list)

blocksizechoices <-(seq(1,18,by=1))

#Validation
MSEs<-c()
predErr<-c()
for(i in 1:length(blocksizechoices)){
  for(t in 1 :72){
    fit <-rangerts::rangerts(y~., data=Val1[(4+h):(T1+t-1),], num.trees=50, mtry=ceiling(LASSOV/3), replace=T, seed=711, bootstrap.ts = "moving", block.size=blocksizechoices[i] )
    pred <- predict(fit, Val1[(T1+t),-1])
    predErr[t]<-as.numeric((pred$predictions - as.numeric(y[T1+t]))^2)  
  }
  
  MSEs[i]<-mean(predErr)
} 
optmblocksize<-blocksizechoices[which.min(MSEs)]


# evaluation
predErr2<-c()
for(t in 1:61){
  fit <- rangerts::rangerts(y~., data=Val1[(T2+t-60):(T2+t-1),], num.trees=50, mtry=ceiling(LASSOV/3), replace=T, seed=711, bootstrap.ts = "moving", block.size=optmblocksize )
  pred <- predict(fit, Val1[(T2+t),-1])
  predErr2[t]<- as.numeric((pred$predictions - y[T2+t])^2)
}

msfeLassoRF <- mean(predErr2)

MSFEs[[horizon]]["LASSORF", targetVar] <- msfeLassoRF

c = matrix(NA, nrow = 61, ncol = 960)
for(t in 1:61){
  fit <- rangerts::rangerts(y~., data=Val[(T2+t-60):(T2+t-1),], num.trees=500, mtry=320, replace=T, seed=711, importance = "impurity" , bootstrap.ts = "moving", block.size=optmblocksize )
  c[t,] <- fit$variable.importance
}

LASSORF_importance[[var]][horizon,] = colMeans(c)
colnames(LASSORF_importance[[var]]) = colnames(X)


rm(coefTracker,X, y, cvScore,lambdaChoises,msfeLassoRF, optLam,predErr1,predErr2,optmblocksize,eval1,eval2,Val1)
