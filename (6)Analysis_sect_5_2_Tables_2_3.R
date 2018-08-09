################################
## File 6
##    Analysis for section 5.2,
##    table 2 and 
##    table 3 (table 1 is not 
##    a table of estimation 
##    results).
##
## Dependencies:
##    - dem.panel.RData
##    - ResultsTable3.RData
################################

library(reshape)
library(plyr)
library(survey)
library(mgcv)
library(gbm)
library(BayesTree)
library(MASS)
library(foreach)
library(doMC)
library(AUC)
library(xtable)


# Clear working space
rm(list=ls())


###################
#### Load data ####
###################
setwd("~/Dropbox/TreeBasedMethods/ReplicationMaterials/")
load("dem.panel.RData")



#### Setting up covariates/formulas


## Recoding year to be binary indicator
dem.panel$y2000<-dem.panel$y2002 <- dem.panel$y2004 <- dem.panel$y2006<-NA
dem.panel$y2000[dem.panel$year==2000] <- 1
dem.panel$y2000[dem.panel$year!=2000] <- 0
dem.panel$y2002[dem.panel$year==2002] <- 1
dem.panel$y2002[dem.panel$year!=2002] <- 0
dem.panel$y2004[dem.panel$year==2004] <- 1
dem.panel$y2004[dem.panel$year!=2004] <- 0
dem.panel$y2006[dem.panel$year==2006] <- 1
dem.panel$y2006[dem.panel$year!=2006] <- 0



## Group the variables here into relevant categories ##
base.vars <- paste(c("camp.length","deminc",
                     "base.poll*deminc","as.factor(year)",
                     "base.und*deminc", "office"
), collapse = "+")
treat.hist <- paste(c("d.gone.neg.l1","d.gone.neg.l2",
                      "d.neg.frac.l3","week*camp.length"
),collapse = "+")
n.cov.hist <- paste(c("s(dem.polls.l1, k = 3)", "neg.rep.l1", "neg.rep.l2", "r.neg.frac.l3",
                      "num.rep.l1", "num.dem.l1",
                      "undother.l1"), collapse = "+")
i.cov.hist <- paste(c("neg.rep.l1", "neg.rep.l2", "r.neg.frac.l2",
                      "s(dem.polls.l1)"), collapse = "+")
xvars.numer.trees<-c("d.gone.neg.l1", "d.gone.neg.l2", "d.neg.frac.l3"
                     , "week", "camp.length", "base.poll", "y2000" ,"y2002"
                     , "y2004", "y2006", "base.und", "office",  "deminc")
xvars.denom.trees<-c("d.gone.neg.l1", "d.gone.neg.l2", "d.neg.frac.l3", "week",
                     "camp.length", "base.poll", "y2000", "y2002", "y2004", "y2006",
                     "base.und", "office", "dem.polls.l1", "dem.polls.l2", "neg.rep.l1",
                     "neg.rep.l2", "r.neg.frac.l2", "r.neg.frac.l3", "num.rep.l1",
                     "num.rep.l2", "num.dem.l1", "num.dem.l2", "undother.l1",
                     "undother.l2", "deminc", "dem.contrib.l1", "dem.contrib.l2",
                     "rep.contrib.l1", "rep.contrib.l2")
i.denom <- as.formula(paste("d.gone.neg ~ ",
                            paste(treat.hist, base.vars,
                                  i.cov.hist, sep = "+")))
n.denom <- as.formula(paste("d.gone.neg ~ ",
                            paste(treat.hist, base.vars,
                                  n.cov.hist, sep = "+")))
numer <- as.formula(paste("d.gone.neg ~ ",
                          paste(treat.hist, base.vars,
                                sep = "+")))
denom.trees<-as.formula(paste("d.gone.neg ~ ",
                              paste(xvars.denom.trees, collapse = "+")))
numer.trees<-as.formula(paste("d.gone.neg ~ ",
                              paste(xvars.numer.trees, collapse = "+")))


## non.comp are the extremely uncompetitive races that fail common
## support from baseline. Thus, this is an analysis is for the
## population of reasonably competive races.
non.comp <- unique(subset(dem.panel, base.poll > 70 | base.poll < 30
                          | is.na(base.poll) ,
                          "race"))[,1]
dem.panel <- subset(dem.panel, !(race %in% non.comp))

## Positive probability subset remove the weeks where there is either
## no chance of staying positive (too many ads to not have negative
## ones), or no chance of going negative (race is too safe, or there
## are no ads)
dem.panel$pos.prob.subset <- (dem.panel$week !=0 &
                                dem.panel$week > dem.panel$first.week +2 &
                                dem.panel$num.dem != 0 )

## Identifying cases with missingness on the explanatory variables.
dem.panel$missing <- aaply(is.na(dem.panel[,c(xvars.denom.trees, xvars.numer.trees)]), 1, any)

### Set up places for the weights and predicted values are to be placed
dem.panel$gam.sw <- dem.panel$gbm.sw <- dem.panel$bart.sw <- NA
dem.panel$gam.pscore <- dem.panel$gbm.pscore <- dem.panel$bart.pscore <- NA


#### Cross validation study

## Extract a test set
set.seed(1234)

## subset down to positive probability subset
subData <- dem.panel[dem.panel$pos.prob.subset,]
rows <- dim(subData)[1]
test <- sample(1:rows, floor(1/10*rows), replace=FALSE)
length(test)
train <-c(1:rows)[!c(1:rows)%in%test]
length(train)
dem.panel.test <- subData[test,]
dim(dem.panel.test)
dem.panel.train <- subData[train,]
dim(dem.panel.train)

### Set up grid for doing 10-fold cross validation only within the training set
min.obs=1
shrink=c(0.005)
trees=seq(1, 2000, by=1)
depth=c(3, 5)
param.grid = expand.grid("min.obs"=min.obs, "shrink"=shrink, "trees"=trees, "depth"=depth)
nSweep<-dim(param.grid)[1]


#### oneFold function

## This function conducts a 10-fold cross validation using one row (x) from param.grid
oneFold <- function(.folds=10, thisData, x){
  set.seed(123456)
  param.grid.row <- param.grid[x,]
  folds.index <- sample(as.numeric(cut(1:nrow(thisData),.folds)))
  outPreds <- rep(NA, nrow(thisData))
  names(outPreds) <- rownames(thisData)
  for(this.fold in 1:.folds){
    inSample <- folds.index!=this.fold
    outSample <- folds.index==this.fold
    gbm.dem<-gbm(denom.trees
                 , data=thisData[inSample,]
                 , interaction.depth=param.grid.row$depth, n.tree=param.grid.row$trees
                 , n.minobsinnode=param.grid.row$min.obs, shrinkage = param.grid.row$shrink
                 , distribution="bernoulli")
    predGbm <- predict.gbm(gbm.dem, newdata=thisData[outSample,]
                           , interaction.depth=param.grid.row$depth, n.tree=param.grid.row$trees
                           , n.minobsinnode=param.grid.row$min.obs, shrinkage = param.grid.row$shrink
                           ,type="response")
    outPreds[rownames(thisData[outSample,])] <- predGbm
  }
  print(x)
  outPreds
}

registerDoMC(detectCores())
cvList <- llply(1:4000, oneFold, .folds=10, thisData=dem.panel.train, .parallel=TRUE)

rmse.cvStudy <- laply(cvList, function(z) mean(c(z-dem.panel.train[, "d.gone.neg"])^2))
bestParam.cvStudy <- param.grid[which(rmse.cvStudy==min(rmse.cvStudy)),] ## this is the parameter grid we will use
bestParam.cvStudy

## Now we will do a small tournament predicting the true test set
yvar <- c("d.gone.neg")
missing.train <- aaply(is.na(dem.panel.train[,xvars.denom.trees]), 1, any)
missing.test <- aaply(is.na(dem.panel.test[,xvars.denom.trees]), 1, any)


## Calculate out-of-sample predictive power for gbm, bart, and gam
dem.panel.test$bart.pred <- dem.panel.test$gbm.pred <- dem.panel.test$gam.pred <- NA

gbm.dem<-gbm(denom.trees, data=dem.panel.train
             , interaction.depth=bestParam.cvStudy$depth, n.tree=bestParam.cvStudy$trees
             , n.minobsinnode=bestParam.cvStudy$min.obs, shrinkage = bestParam.cvStudy$shrink,
             distribution="bernoulli")
dem.panel.test$gbm.pred <- predict.gbm(gbm.dem, newdata=dem.panel.test, n.tree=bestParam.cvStudy$trees
                                       , n.minobsinnode=bestParam.cvStudy$min.obs, shrinkage = bestParam.cvStudy$shrink, type="response")
dem.panel.test$gbm.pred[missing.test] <- NA

## Bart using default parameters
bart.dem <- bart(x.train=dem.panel.train[!missing.train, xvars.denom.trees],
                 y.train=as.factor(dem.panel.train[!missing.train, yvar]),
                 x.test=dem.panel.test[!missing.test, xvars.denom.trees], ndpost=20000, nskip=10000)
dem.panel.test[!missing.test, "bart.pred"] <- pnorm(colMeans(bart.dem$yhat.test))


## Gam
dem.yesInc <- gam(i.denom, data = dem.panel.train[dem.panel.train$deminc==1,] , family ="binomial")
dem.nonInc <- gam(n.denom, data = dem.panel.train[dem.panel.train$deminc==0,] , family ="binomial")
predGam.yesInc <- predict.gam(dem.yesInc, newdata=dem.panel.test[dem.panel.test$deminc==1,], type="response")
predGam.nonInc <- predict.gam(dem.nonInc, newdata=dem.panel.test[dem.panel.test$deminc==0,], type="response")
dem.panel.test$gam.pred[dem.panel.test$deminc==1] <-  predGam.yesInc
dem.panel.test$gam.pred[dem.panel.test$deminc==0] <-  predGam.nonInc
dem.panel.test$gam.pred[missing.test] <- NA

## Calculate the errors for the test set
errors.gbm <- c(dem.panel.test[dem.panel.test$pos.prob.subset,"d.gone.neg"]-dem.panel.test$gbm.pred[dem.panel.test$pos.prob.subset])
errors.gam <- c(dem.panel.test[dem.panel.test$pos.prob.subset,"d.gone.neg"]-dem.panel.test$gam.pred[dem.panel.test$pos.prob.subset])
errors.bart <- c(dem.panel.test[dem.panel.test$pos.prob.subset,"d.gone.neg"]-dem.panel.test$bart.pred[dem.panel.test$pos.prob.subset])


outMat<-matrix(NA, nrow=3, ncol=5)

outMat[1,1] <- 1
outMat[2,1] <- cor(dem.panel.test$gbm.pred, dem.panel.test$gam.pred, use="pairwise.complete")
outMat[3,1] <- cor(dem.panel.test$bart.pred, dem.panel.test$gam.pred, use="pairwise.complete")
outMat[1,2] <- sqrt(mean(errors.gam[dem.panel.test$pos.prob.subset]^2,na.rm=TRUE))
outMat[2,2] <- sqrt(mean(errors.gbm[dem.panel.test$pos.prob.subset]^2, na.rm=TRUE))
outMat[3,2] <- sqrt(mean(errors.bart[dem.panel.test$pos.prob.subset]^2, na.rm=TRUE))

sqrt(mean(errors.bart^2, na.rm=TRUE))


outMat[1,3] <- auc(roc(dem.panel.test$gam.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))
outMat[2,3] <-auc(roc(dem.panel.test$gbm.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))
outMat[3,3] <-auc(roc(dem.panel.test$bart.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))


outMat[1,4] <- auc(sensitivity(dem.panel.test$gam.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))
outMat[2,4] <- auc(sensitivity(dem.panel.test$gbm.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))
outMat[3,4] <- auc(sensitivity(dem.panel.test$bart.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))


outMat[1,5] <- auc(specificity(dem.panel.test$gam.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))
outMat[2,5] <- auc(specificity(dem.panel.test$gbm.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))
outMat[3,5] <- auc(specificity(dem.panel.test$bart.pred[dem.panel.test$pos.prob.subset], labels=as.factor(dem.panel.test$d.gone.neg[dem.panel.test$pos.prob.subset])))



## Now calculate the fit statistics for a 10-fold CV within the training set alone
withinTraining <- function(.folds=1, thisData, params){
  set.seed(54321)
  folds.index <- sample(as.numeric(cut(1:nrow(thisData),.folds)))
  outPredsGAM <- outPredsGBM <- outPredsBART <- rep(NA, nrow(thisData))
  names(outPredsGAM) <- names(outPredsGBM) <- names(outPredsBART) <- rownames(thisData)
  for(this.fold in 1:.folds) {
    inSample <- folds.index!=this.fold
    outSample <- folds.index==this.fold
    gbm.dem<-gbm(denom.trees
                 , data=thisData[inSample,]
                 , interaction.depth=params$depth, n.tree=params$trees
                 , n.minobsinnode=params$min.obs, shrinkage = params$shrink
                 , distribution="bernoulli")
    predGbm <- predict.gbm(gbm.dem, newdata=thisData[outSample,]
                           , interaction.depth=params$depth, n.tree=params$trees
                           , n.minobsinnode=params$min.obs, shrinkage = params$shrink
                           ,type="response")
    outPredsGBM[rownames(thisData[outSample,])] <- predGbm
    
    bart.dem <- bart(x.train=thisData[inSample, xvars.denom.trees],
                     y.train=as.factor(thisData[inSample, yvar]),
                     x.test=thisData[outSample, xvars.denom.trees], ndpost=20000, nskip=10000
                     ,verbose=FALSE)
    outPredsBART[rownames(thisData[outSample,])] <- pnorm(colMeans(bart.dem$yhat.test))
    
    dem.yesInc <- gam(i.denom, data = thisData[thisData$deminc==1 & inSample,] , family ="binomial")
    dem.nonInc <- gam(n.denom, data = thisData[thisData$deminc==0 & inSample,] , family ="binomial")
    predGam.yesInc <- predict.gam(dem.yesInc, newdata=thisData[thisData$deminc==1 & outSample,]
                                  , type="response")
    predGam.nonInc <- predict.gam(dem.nonInc, newdata=thisData[thisData$deminc==0 & outSample,]
                                  , type="response")
    outPredsGAM[thisData$deminc==1 & outSample] <-  predGam.yesInc
    outPredsGAM[thisData$deminc==0 & outSample] <-  predGam.nonInc
  }
  list(gam=outPredsGAM, gbm=outPredsGBM, bart=outPredsBART)
}

##Warning: the following can take some time to run. 
cvWithinTraining<-withinTraining(.folds=10
                                 , thisData=dem.panel.train[dem.panel.train$pos.prob.subset & !missing.train,]
                                 , params=bestParam.cvStudy)

dem.panel.train$gam.pred <- dem.panel.train$bart.pred <- dem.panel.train$gbm.pred <- NA
dem.panel.train$gam.pred[dem.panel.train$pos.prob.subset & !missing.train] <- cvWithinTraining$gam
dem.panel.train$gbm.pred[dem.panel.train$pos.prob.subset & !missing.train] <- cvWithinTraining$gbm
dem.panel.train$bart.pred[dem.panel.train$pos.prob.subset & !missing.train] <- cvWithinTraining$bart

table(is.na(dem.panel.train$bart.pred))

## Calculate the errors for the test set
errors.gbm.train <- c(dem.panel.train[dem.panel.train$pos.prob.subset,"d.gone.neg"]-dem.panel.train$gbm.pred[dem.panel.train$pos.prob.subset])
errors.gam.train <- c(dem.panel.train[dem.panel.train$pos.prob.subset,"d.gone.neg"]-dem.panel.train$gam.pred[dem.panel.train$pos.prob.subset])
errors.bart.train <- c(dem.panel.train[dem.panel.train$pos.prob.subset,"d.gone.neg"]-dem.panel.train$bart.pred[dem.panel.train$pos.prob.subset])


outMat2<-matrix(NA, nrow=3, ncol=5)

outMat2[1,1] <- 1
outMat2[2,1] <- cor(dem.panel.train$gbm.pred, dem.panel.train$gam.pred, use="pairwise.complete")
outMat2[3,1] <- cor(dem.panel.train$bart.pred, dem.panel.train$gam.pred, use="pairwise.complete")
outMat2[1,2] <- sqrt(mean(errors.gam.train[dem.panel.train$pos.prob.subset]^2 , na.rm=TRUE))
outMat2[2,2] <- sqrt(mean(errors.gbm.train[dem.panel.train$pos.prob.subset]^2 , na.rm=TRUE))
outMat2[3,2] <- sqrt(mean(errors.bart.train[dem.panel.train$pos.prob.subset]^2 , na.rm=TRUE))

outMat2[1,3] <- auc(roc(dem.panel.train$gam.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))
outMat2[2,3] <-auc(roc(dem.panel.train$gbm.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))
outMat2[3,3] <-auc(roc(dem.panel.train$bart.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))


outMat2[1,4] <- auc(sensitivity(dem.panel.train$gam.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))
outMat2[2,4] <- auc(sensitivity(dem.panel.train$gbm.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))
outMat2[3,4] <- auc(sensitivity(dem.panel.train$bart.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))


outMat2[1,5] <- auc(specificity(dem.panel.train$gam.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))
outMat2[2,5] <- auc(specificity(dem.panel.train$gbm.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))
outMat2[3,5] <- auc(specificity(dem.panel.train$bart.pred[dem.panel.train$pos.prob.subset], labels=as.factor(dem.panel.train$d.gone.neg[dem.panel.train$pos.prob.subset])))


##Table 2
xtable(outMat2, digits=3)
xtable(outMat, digits=3)
dim(dem.panel.train)
dim(dem.panel.test)



###########################
#### pan.prod function ####
###########################
## This function was originally located in panel-utils.R
## I have moved it here
pan.prod <- function(x,ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cumprod(na.omit(x))
    return(xt)
  }))
}


#######################
#### iptw function ####
#######################
## This function was originally located in panel-utils.R
## I have moved it here, and trimmed it down to focus only on the binomial model

## ones is an argument to manually set certain observations as having
## weight 1. This is needed due to common support or positivity
## problems. Basically, in certain ranges of the covariates, there are
## either random or structural zeros. For negative advertising, a
## strutural zero occurs when there are no ads. Random zeros generally
## occur in the extremes of the polls and number of ads
## covariates. Extremely uncompetitive races never see negativity and
## extreme competitive race (in the number of ads) almost always go
## negative.

iptw <- function(denominator, numerator, data, id, time, subset,
                 ones = NULL) {
  
  if (missing(ones)) {
    ones <- rep(FALSE, nrow(data))
  }
  
  ## Subset via the subset.data.frame function
  if (missing(subset)) {
    r <- TRUE
  } else {
    e <- substitute(subset)
    r <- eval(e, data, parent.frame())
    if (!is.logical(r)) {
      stop("'subset' must evaluate to logical")
    }
    r <- r & !is.na(r)
    if (sum(r) == 0) {
      stop("no observations in the subset")
    }
  }
  
  odata <- data
  data <- data[r,]
  ones <- ones[r]
  idmat <- data[,c(id, time)]
  d.fit <- rep(NA, nrow(data))
  n.fit <- rep(NA, nrow(data))
  names(d.fit) <- rownames(data)
  names(n.fit) <- rownames(data)
  trvar <- all.vars(denominator)[1]
  tr <- data[, trvar]
  
  set.seed(99)
  dem <- gam(denominator, data = data[!ones,], family = "binomial")
  nom <- gam(numerator, data = data[!ones,], family = "binomial")
  
  d.nms <- rownames(model.matrix(dem))
  n.nms <- rownames(model.matrix(nom))
  
  d.fit[d.nms] <- dem$fitted.values
  n.fit[n.nms] <- nom$fitted.values
  
  d.pr.tr <- ifelse(tr == 1, d.fit, 1-d.fit)
  n.pr.tr <- ifelse(tr == 1, n.fit, 1-n.fit)
  
  sw <- n.pr.tr/d.pr.tr
  names(sw) <- rownames(data)
  sw[ones] <- 1
  sw <- pan.prod(sw, idmat[,1])
  
  out <- list(sw = sw, d.pscore = d.fit, n.pscore = n.fit)
  return(out)
}


###############################
#### Calculate GAM weights ####
###############################
dem.wmod.ni <- iptw(denominator = n.denom, numerator = numer,
                    data = dem.panel, id = "race", time = "week",
                    subset = deminc == 0,
                    ones  = !dem.panel$pos.prob.subset)
dem.wmod.inc <- iptw(denominator = i.denom, numerator = numer,
                     data = dem.panel, id = "race", time = "week",
                     subset = deminc == 1,
                     ones  = !dem.panel$pos.prob.subset)
### Put weights into the main panel
dem.panel$gam.sw[dem.panel$deminc == 1] <- dem.wmod.inc$sw
dem.panel$gam.sw[dem.panel$deminc == 0] <- dem.wmod.ni$sw
dem.panel$gam.pscore[dem.panel$deminc == 1] <- dem.wmod.inc$d.pscore
dem.panel$gam.pscore[dem.panel$deminc == 0] <- dem.wmod.ni$d.pscore
hist(dem.panel$gam.sw)



###############################
#### Calculate GBM weights ####
###############################

iptw.gbm <- function(denominator, numerator, data, id, time, ones = NULL) {
  
  if (missing(ones)) {
    ones <- rep(FALSE, nrow(data))
  }
  
  idmat <- data[,c(id, time)]
  d.fit <- rep(NA, nrow(data))
  n.fit <- rep(NA, nrow(data))
  names(d.fit) <- rownames(data)
  names(n.fit) <- rownames(data)
  trvar <- all.vars(denominator)[1]
  tr <- data[, trvar]
  
  set.seed(99)
  dem<-gbm(denominator, distribution="bernoulli", data=data[!ones,]
           , interaction.depth=3, n.tree=828, n.minobsinnode=1, shrinkage=0.005)
  nom<-gbm(numerator,  distribution="bernoulli"
           , data=data[!ones,]
           , interaction.depth=3, n.tree=828, shrinkage=0.005, n.minobsinnode=1)
  d.nms <- rownames(data[!ones,])
  n.nms <- rownames(data[!ones,])
  
  d.fit[d.nms]<-plogis(dem$fit)
  n.fit[n.nms] <- plogis(nom$fit)
  
  d.pr.tr <- ifelse(tr == 1, d.fit, 1-d.fit)
  n.pr.tr <- ifelse(tr == 1, n.fit, 1-n.fit)
  
  sw <- n.pr.tr/d.pr.tr
  names(sw) <- rownames(data)
  sw[ones] <- 1
  sw <- pan.prod(sw, idmat[,1])
  
  out <- list(sw = sw, d.pscore = d.fit, n.pscore = n.fit)
  return(out)
}

gbm.fit <- iptw.gbm(denominator = denom.trees, numerator = numer.trees,
                    data = dem.panel, id = "race", time = "week",
                    ones  = !dem.panel$pos.prob.subset)


dem.panel$gbm.sw<-gbm.fit$sw
dem.panel$gbm.pscore<-gbm.fit$d.pscore


###############################
#### Calculate BART weights ####
###############################
iptw.bart <- function(y, xvars.denominator, xvars.numerator, data, id, time, ones = NULL, missing) {
  
  if (missing(ones)) {
    ones <- rep(FALSE, nrow(data))
  }
  
  idmat <- data[,c(id, time)]
  d.fit <- rep(NA, nrow(data))
  n.fit <- rep(NA, nrow(data))
  names(d.fit) <- rownames(data)
  names(n.fit) <- rownames(data)
  trvar <- y
  tr <- data[, trvar]
  set.seed(99)
  dem <- bart(x.train=data[!ones & !missing,xvars.denominator],
              y.train=data[!ones & !missing, y], ndpost=10000, nskip=10000,verbose=TRUE)
  nom <- bart(x.train=data[!ones & !missing, xvars.numerator],
              y.train=data[!ones & !missing, y], ndpost=10000, nskip=10000,verbose=TRUE)
  
  
  d.nms <- rownames(data[!ones & !missing,])
  n.nms <- rownames(data[!ones & !missing,])
  
  d.fit[d.nms]<-pnorm(colMeans(dem$yhat.train))
  n.fit[n.nms] <- pnorm(colMeans(nom$yhat.train))
  
  d.pr.tr <- ifelse(tr == 1, d.fit, 1-d.fit)
  n.pr.tr <- ifelse(tr == 1, n.fit, 1-n.fit)
  
  sw <- n.pr.tr/d.pr.tr
  names(sw) <- rownames(data)
  sw[ones] <- 1
  sw <- pan.prod(sw, idmat[,1])
  
  out <- list(sw = sw, d.pscore = d.fit, n.pscore = n.fit)
  return(out)
}

bart.fit <- iptw.bart(y="d.gone.neg", xvars.denominator=xvars.denom.trees,
                      xvars.numerator=xvars.numer.trees, data=dem.panel, id = "race", time = "week",
                      ones  = !dem.panel$pos.prob.subset, missing=dem.panel$missing)

dem.panel$bart.sw <- bart.fit$sw
dem.panel$bart.pscore <- bart.fit$d.pscore

hist(dem.panel$bart.sw)
plot(dem.panel$bart.pscore[!dem.panel$missing], dem.panel$gam.pscore[!dem.panel$missing])
cor(dem.panel$bart.pscore, dem.panel$gam.pscore, use="pairwise.complete")

####### Re-shaping for the analysis

dem.final <- dem.panel[dem.panel$week == -1, ]
miss <- apply(is.na(dem.final[, c("gam.sw", "gbm.sw", "bart.sw")]), 1, any)
dim(dem.final)
dem.final <- dem.final[!miss,]


cor(dem.final[,c("gam.sw",  "bart.sw", "gbm.sw")], use="pairwise.complete")

## Compare weights

par(mfrow=c(2,3), mgp=c(1,.1,0), mar=c(2.5,2.5,3,1), tcl=0, din=c(10,4))
plot(dem.final[,c("gam.sw",  "gbm.sw")], pch=19, xlab="GBM weights", ylab="GAM weights", main="")
title(line=.5, main="GBM vs. GAM (r=0.309)")
abline(lm(gbm.sw~gam.sw, data=dem.final), lty=2)
plot(dem.final[,c("gam.sw",  "bart.sw")], pch=19, xlab="BART weights", ylab="GAM weights", main="")
title(line=.5, main="BART vs. GAM (r=0.589)")
abline(lm(bart.sw~gam.sw, data=dem.final), lty=2)
plot(dem.final[,c("gbm.sw",  "bart.sw")], pch=19, xlab="BART weights", ylab="GBM weights", main="")
title(line=.5, main="BART vs. GBM (r=0.738)")
abline(lm(gbm.sw~bart.sw, data=dem.final), lty=2)
plot(density(dem.final[,c("gam.sw")]), xlab="", main="")
rug(dem.final[,c("gam.sw")])
title(line=.5, main="GAM weights")
plot(density(dem.final[,c("gbm.sw")], na.rm=TRUE), xlab="", main="")
rug(dem.final[,c("gbm.sw")])
title(line=.5, main="GBM weights")
plot(density(dem.final[,c("bart.sw")], na.rm=TRUE), xlab="", main="")
rug(dem.final[,c("bart.sw")])
title(line=.5, main="BART weights")




#######################
### Analysis Models ###
#######################

## First, we perform the regressions (iptw, naive, and adjusted),
## switching the incumbent indicator to get the effect for each
## group. The SEs for the IPTW are asympototically correct, but I
## bootstrap the SEs below just to be sure.

# adj.vars <- paste(c("r.neg.frac", "poll.mean", "I(poll.mean^2)",
#                     "I(d.total/10)", "I(r.total/10)"
# ),
# collapse = "+")
causal.des <- svydesign(ids = ~ race, weights = ~ gam.sw, data =
                          dem.final)
causal.des.gbm <- svydesign(ids = ~ race, weights = ~ gbm.sw, data =
                              dem.final)
causal.des.bart <- svydesign(ids = ~ race, weights = ~ bart.sw, data =
                               dem.final)

## Calculating effects for non-incumbents
## Formula for IPTW and Naive: only baseline covariates
share.form <- paste("demprcnt ~ ",
                    "deminc * d.neg.rec + deminc * I(d.neg.dur - d.neg.rec)+as.factor(year)+",
                    base.vars)
share.form <- as.formula(share.form)

## incumbent vote share ##
inc.form <- paste("demprcnt ~ ",
                  "I(deminc == 0) * d.neg.rec   + I(deminc == 0) * I(d.neg.dur - d.neg.rec) +",
                  gsub("deminc", "I(deminc == 0)", base.vars))
inc.form <- as.formula(inc.form)

dim(dem.final)
naive.share <- lm(share.form, data = dem.final)
summary(naive.share)
gam.share <-svyglm(share.form, causal.des)
summary(gam.share)
gbm.share <-svyglm(share.form, causal.des.gbm)
summary(gbm.share)
bart.share <-svyglm(share.form, causal.des.bart)
summary(bart.share)

naive.share.inc <- lm(inc.form, data = dem.final)
summary(naive.share.inc)
gam.share.inc <- svyglm(inc.form, causal.des)
summary(gam.share.inc)
gbm.share.inc <- svyglm(inc.form, causal.des.gbm)
summary(gbm.share.inc)
bart.share.inc <- svyglm(inc.form, causal.des.bart)
summary(bart.share.inc)



#### BOOTSTRAPPED SEs ######
## bootstrap the race
sims <-500

## set up for bootstrap samples.  We are sampling by race.
race.list <- unique(dem.final$race)
race.lengths <- c()
race.rows <- list()
for (i in 1:length(race.list)) {
  race.rows[[i]] <- rownames(dem.panel[which(dem.panel$race==race.list[i]),])
  race.lengths <- c(race.lengths, nrow(dem.panel[which(dem.panel$race ==race.list[i]),]))
}
names(race.rows) <- race.list
names(race.lengths) <- race.list

registerDoMC(detectCores())

out<-foreach(i=1:sims) %dopar% {
  cat(i, " ")
  if (i %% 10 == 0) cat("\n")
  set.seed(i)
  races <- sample(race.list, replace = TRUE,
                  size = length(race.list))
  races <- as.character(races)
  star <- unlist(race.rows[races])
  bootid <- rep(1:length(races), times = race.lengths[races])
  dem.star <- dem.panel[star,]
  dem.star$bootid <- bootid
  dem.wmod.ni <- iptw(denominator = n.denom, numerator = numer,
                      data = dem.star, id = "bootid", time = "week",
                      subset = deminc == 0,
                      ones  = !dem.star$pos.prob.subset)
  dem.wmod.inc <- iptw(denominator = i.denom, numerator = numer,
                       data = dem.star, id = "bootid", time = "week",
                       subset = deminc == 1,
                       ones  = !dem.star$pos.prob.subset)
  ### Put weights into the main panel
  dem.star$gam.sw[dem.star$deminc == 1] <- dem.wmod.inc$sw
  dem.star$gam.sw[dem.star$deminc == 0] <- dem.wmod.ni$sw
  gbm.fit <- iptw.gbm(denominator = denom.trees, numerator = numer.trees,
                      data = dem.star, id = "bootid", time = "week",
                      ones  = !dem.star$pos.prob.subset)
  dem.star$gbm.sw<-gbm.fit$sw
  bart.fit <- iptw.bart(y="d.gone.neg", xvars.denominator=xvars.denom.trees,
                        xvars.numerator=xvars.numer.trees, data=dem.star, id = "bootid", time = "week",
                        ones  = !dem.star$pos.prob.subset, missing=dem.star$missing)
  dem.star$bart.sw <- bart.fit$sw
  ## Reorganize the data
  dem.star <- dem.star[dem.star$week == -1,]
  .miss <- apply(is.na(dem.star[, c("gam.sw", "gbm.sw", "bart.sw")]), 1, any)
  dem.star <- dem.star[!.miss,]
  ## Set up the designs
  causal.des <- svydesign(ids = ~ bootid, weights = ~ gam.sw, data = dem.star)
  causal.des.gbm <- svydesign(ids = ~ bootid, weights = ~ gbm.sw, data =dem.star)
  causal.des.bart <- svydesign(ids = ~ bootid, weights = ~ bart.sw, data =dem.star)
  ## Fit the models
  naive.share.star <- lm(share.form, data = dem.star)
  gam.share.star <-svyglm(share.form, causal.des)
  gbm.share.star <-svyglm(share.form, causal.des.gbm)
  bart.share.star <-svyglm(share.form, causal.des.bart)
  naive.inc.star <- lm(inc.form, data = dem.star)
  gam.inc.star <- svyglm(inc.form, causal.des)
  gbm.inc.star <- svyglm(inc.form, causal.des.gbm)
  bart.inc.star <- svyglm(inc.form, causal.des.bart)
  ## Output the results
  print(i)
  list(naive.share=naive.share.star$coef, naive.inc=naive.inc.star$coef, gam.share=gam.share.star$coef,
       gam.inc=gam.inc.star$coef, gbm.share=gbm.share.star$coef, gbm.inc=gbm.inc.star$coef, bart.share=bart.share.star$coef, bart.inc=bart.inc.star$coef)
}
####


#####################################################
## FOR REPLICATION PURPOSES, YOU MUST LOAD OBJECT
load("ResultsTable3.RData")
## This is because the bart function does not take a random seed.
#####################################################


## ALL SUBSEQUENT LINES SHOULD RUN AFTER THE LOAD.
str(out[[1]])
naive.share.boots<-laply(out, function(x){unlist(x[["naive.share"]])})
naive.inc.boots<-laply(out, function(x){unlist(x[["naive.inc"]])})
gam.share.boots<-laply(out, function(x){unlist(x[["gam.share"]])})
gam.inc.boots<-laply(out, function(x){unlist(x[["gam.inc"]])})
gbm.share.boots<-laply(out, function(x){unlist(x[["gbm.share"]])})
gbm.inc.boots<-laply(out, function(x){unlist(x[["gbm.inc"]])})
bart.share.boots<-laply(out, function(x){unlist(x[["bart.share"]])})
bart.inc.boots<-laply(out, function(x){unlist(x[["bart.inc"]])})

colnames(naive.share.boots)
colnames(naive.inc.boots)

regOut <- function(x,name){
  results <- rbind(round(mean(x),3), paste(round(quantile(x, .025),3), round(quantile(x, .975),3),sep=", "), paste(round(quantile(x, .05),3), round(quantile(x, .95),3),sep=", "))
  rownames(results) <- c(name,"95% Bootstrapped CI","90% Bootstrapped CI")
  return(results)
}

# Table 3
tab_3 <- cbind(rbind(regOut(naive.share.boots[,"d.neg.rec"],"No weights"),
                     regOut(gam.share.boots[,"d.neg.rec"],"GAM weights"),
                     regOut(gbm.share.boots[,"d.neg.rec"],"GBM weights"),
                     regOut(bart.share.boots[,"d.neg.rec"],"BART weights")),
               rbind(regOut(naive.inc.boots[,"d.neg.rec"],"No weights"),
                     regOut(gam.inc.boots[,"d.neg.rec"],"GAM weights"),
                     regOut(gbm.inc.boots[,"d.neg.rec"],"GBM weights"),
                     regOut(bart.inc.boots[,"d.neg.rec"],"BART weights")))
colnames(tab_3) <- c("Democratic\\Non-incumbents","Democratic Incumbents")
print(tab_3)

