#!/usr/bin/env Rscript
#########################################
## File 4:
##    Simulations for section 5.1 
##
## WARNING: This script was 
## designed to run on a cluster
## with over 400 cores available
## using open MPI and the 
## the LSF scheduler. 
## The scheduler batch file is
## given in "contest.job".
##
## Dependencies:
##  - File (2)Synth_data_for_5_1.R
##  - File SynthData.RData
##  - File CVFunction.R
#######################################

library(plyr)
library(MASS)
library(lattice)
library(latticeExtra)
#library(doMC)
library(doMPI)

##Cross validation function
source("CVFunction.R")
source("(2)Synth_data_for_5_1.R")

## Load simulated data
load("SynthData.RData")

## Start cluster and worker loop
cl  <- startMPIcluster()
registerDoMPI(cl)
#registerDoMC(20)

library(randomForest)
library(gbm)
library(rpart)
library(BayesTree)

### Model estimation

#Parameters to CV
krls.par <- list(derivative=FALSE
                 ,vcov=FALSE)
nn.par <- list(size=seq(1,15,by=2)
               ,decay=c(0.1)
               ,trace=FALSE
               ,linout=TRUE
               ,maxit=9.5e3)
gam.par <- list()
rpart.par <- list(minbucket=c(3,5,10)
                  ,maxdepth=c(3,5)
                  ,cp=c(0.1,0.01,0.001)
                  ,maxsurrogate=0)
rf.par <- list(ntree=c(50,100,500,1000,2500,5000,10000)
               ,mtry=c(5,10,15)
               ,nodesize=c(3,5))
gbm.par <- list(n.trees=c(seq(1,2500,by=3))
                ,interaction.depth=c(3,5)
                ,shrinkage=c(.01)
                ,bag.fraction=1
                ,n.minobsnode=1)
bart.par <- list(verbose=FALSE) #Use default hyperpriors


##Additive Model
#KRLS

#print(dim(add.data.tr))
#print(dim(add.data.test))
#debug(cross.val)

cat("Working on KRLS additive DGP.\n")
add.dgp.krls <- alply(add.data.tr,3
                      ,cross.val
                      ,.data.test=add.data.test
                      ,.formula=Y~.
                      ,.model="krls"
                      ,.params=krls.par
                      ,.parallel=TRUE
                      ,.paropts=list(.packages=c("plyr")
                                     ,.export=c("add.data.tr","add.data.test")
                      )
)

add.dgp.krls.rmse <- laply(add.dgp.krls,function(x)x$rmse)
add.dgp.krls.model  <- add.dgp.krls[[which.min(add.dgp.krls.rmse)]]$mod.res

#Neural Net
cat("Working on NN additive DGP.\n")
add.dgp.nn <- alply(add.data.tr,3
                    ,cross.val
                    ,.data.test=add.data.test
                    ,.formula=Y~.
                    ,.model="nnet.formula"
                    ,.params=nn.par
                    ,.paropts=list(.packages=c("plyr")
                                   ,.export=c("add.data.tr","add.data.test")
                    )
                    ,.parallel=TRUE
)

add.dgp.nn.rmse <- laply(add.dgp.nn,function(x)x$rmse)
add.dgp.nn.model  <- add.dgp.nn[[which.min(add.dgp.nn.rmse)]]$mod.res


#Generalized additive model
cat("Working on GAM additive DGP.\n")
add.dgp.gam <- alply(add.data.tr,3
                     ,cross.val
                     ,.data.test=add.data.test
                     ,.formula=as.formula(paste("Y~s(X1,k=20)+s(X2,k=20)+s(X3,k=20)+s(X10,k=20)+X4+X5+X6+X7+X8+X9+",paste(colnames(add.data.tr)[-c(1:11)],sep="",collapse="+"),sep=""))
                     ,.model="gam"
                     ,.params=gam.par
                     ,.parallel=TRUE
)

add.dgp.gam.rmse <- laply(add.dgp.gam,function(x)x$rmse)
add.dgp.gam.model  <- add.dgp.gam[[which.min(add.dgp.gam.rmse)]]$mod.res

#Single tree
cat("Working on CART additive DGP.\n")
add.dgp.cart <- alply(add.data.tr,3
                      ,cross.val
                      ,.data.test=add.data.test
                      ,.formula=Y~.
                      ,.model="rpart"
                      ,.params=rpart.par
                      ,.parallel=TRUE
)
add.dgp.cart.rmse <- laply(add.dgp.cart,function(x)x$rmse)
add.dgp.cart.model  <- add.dgp.cart[[which.min(add.dgp.cart.rmse)]]$mod.res

#CV & estmation Random Forest
cat("Working on RF additive DGP.\n")
add.dgp.rf<- alply(add.data.tr,3
                   ,cross.val
                   ,.data.test=add.data.test
                   ,.formula=Y~.
                   ,.model="randomForest.formula"
                   ,.params=rf.par
                   ,.parallel=TRUE
)
add.dgp.rf.rmse <- laply(add.dgp.rf,function(x)x$rmse)
add.dgp.rf.model  <- add.dgp.rf[[which.min(add.dgp.rf.rmse)]]$mod.res

#CV & estimation GBM
cat("Working on GBM additive DGP.\n")
add.dgp.gbm<- alply(add.data.tr,3
                    ,cross.val
                    ,.data.test=add.data.test
                    ,.formula=Y~.
                    ,.model="gbm"
                    ,.params=gbm.par
                    ,.parallel=TRUE
)
add.dgp.gbm.rmse <- laply(add.dgp.gbm,function(x)x$rmse)
add.dgp.gbm.model  <- add.dgp.gbm[[which.min(add.dgp.gbm.rmse)]]$mod.res



#CV & estimation BART
cat("Working on BART additive DGP.\n")
add.dgp.bart <-  alply(add.data.tr,3
                       ,cross.val
                       ,.data.test=add.data.test
                       ,.formula=Y~.
                       ,.model="bart"
                       ,.params=bart.par
                       ,.parallel=TRUE
)
add.dgp.bart.rmse <- laply(add.dgp.bart,function(x)x$rmse)
add.dgp.bart.model  <- add.dgp.bart[[which.min(add.dgp.bart.rmse)]]$mod.res



#RMSE across folds

rmse.dataAdd <- do.call(rbind,list(
  data.frame(rmse=add.dgp.cart.rmse,model="Single Tree (c.v.)")
  ,data.frame(rmse=add.dgp.rf.rmse,model="Random Forest (c.v.)")
  ,data.frame(rmse=add.dgp.gbm.rmse,model="GBM (c.v.)")
  ,data.frame(rmse=add.dgp.bart.rmse,model="BART")
  ,data.frame(rmse=add.dgp.krls.rmse,model="KRLS (c.v.)")
  ,data.frame(rmse=add.dgp.nn.rmse,model="Neural Net (c.v.)")
  ,data.frame(rmse=add.dgp.gam.rmse,model="GAM")
))



##Additive and Multiplicative Model
#KRLS


cat("Working on KRLS additive multi DGP.\n")
addmult.dgp.krls <- alply(addmulti.data.tr,3
                          ,cross.val
                          ,.data.test=addmulti.data.test
                          ,.formula=Y~.
                          ,.model="krls"
                          ,.params=krls.par
                          ,.parallel=TRUE
)
addmult.dgp.krls.rmse <- laply(addmult.dgp.krls,function(x)x$rmse)
addmult.dgp.krls.model  <- addmult.dgp.krls[[which.min(addmult.dgp.krls.rmse)]]$mod.res


#Neural Net

cat("Working on NN additive multi DGP.\n")
addmult.dgp.nn <- alply(addmulti.data.tr,3
                        ,cross.val
                        ,.data.test=addmulti.data.test
                        ,.formula=Y~.
                        ,.model="nnet.formula"
                        ,.params=nn.par
                        ,.parallel=TRUE
)
addmult.dgp.nn.rmse <- laply(addmult.dgp.nn,function(x)x$rmse)
addmult.dgp.nn.model  <- addmult.dgp.nn[[which.min(addmult.dgp.nn.rmse)]]$mod.res


#Generalized additive model

cat("Working on gam addtive multi DGP.\n")
addmult.dgp.gam <- alply(addmulti.data.tr,3
                         ,cross.val
                         ,.data.test=addmulti.data.test
                         ,.formula=as.formula(paste("Y~s(X1,k=20)+s(X2,k=20)+s(X10,k=20)+X3+X4+X5+X6+X7+X8+X9+",paste(colnames(add.data.tr)[-c(1:11)],sep="",collapse="+"),sep=""))
                         ,.model="gam"
                         ,.params=gam.par
                         ,.parallel=TRUE
)
addmult.dgp.gam.rmse <- laply(addmult.dgp.gam,function(x)x$rmse)
addmult.dgp.gam.model  <- addmult.dgp.gam[[which.min(addmult.dgp.gam.rmse)]]$mod.res

#Single tree
cat("Working on cart addtive multi DGP.\n")
addmult.dgp.cart <- alply(addmulti.data.tr,3
                          ,cross.val
                          ,.data.test=addmulti.data.test
                          ,.formula=Y~.
                          ,.model="rpart"
                          ,.params=rpart.par
                          ,.parallel=TRUE
)
addmult.dgp.cart.rmse <- laply(addmult.dgp.cart,function(x)x$rmse)
addmult.dgp.cart.model  <- addmult.dgp.cart[[which.min(addmult.dgp.gam.rmse)]]$mod.res

#CV & estmation Random Forest
cat("Working on RF addtive multi DGP.\n")
addmult.dgp.rf<- alply(addmulti.data.tr,3
                       ,cross.val
                       ,.data.test=addmulti.data.test
                       ,.formula=Y~.
                       ,.model="randomForest.formula"
                       ,.params=rf.par
                       ,.parallel=TRUE
)
addmult.dgp.rf.rmse <- laply(addmult.dgp.rf,function(x)x$rmse)
addmult.dgp.rf.model  <- addmult.dgp.rf[[which.min(addmult.dgp.rf.rmse)]]$mod.res

#CV & estimation GBM
cat("Working on GBM addtive multi DGP.\n")
addmult.dgp.gbm<- alply(addmulti.data.tr,3
                        ,cross.val
                        ,.data.test=addmulti.data.test
                        ,.formula=Y~.
                        ,.model="gbm"
                        ,.params=gbm.par
                        ,.parallel=TRUE
)
addmult.dgp.gbm.rmse <- laply(addmult.dgp.gbm,function(x)x$rmse)
addmult.dgp.gbm.model  <- addmult.dgp.gbm[[which.min(addmult.dgp.gbm.rmse)]]$mod.res


#CV & estimation BART
cat("Working on BART addtive multi DGP.\n")
addmult.dgp.bart <-  alply(addmulti.data.tr,3
                           ,cross.val
                           ,.data.test=addmulti.data.test
                           ,.formula=Y~.
                           ,.model="bart"
                           ,.params=bart.par
                           ,.parallel=TRUE
)
addmult.dgp.bart.rmse <- laply(addmult.dgp.bart,function(x)x$rmse)
addmult.dgp.bart.model  <- addmult.dgp.bart[[which.min(addmult.dgp.bart.rmse)]]$mod.res


#RMSE across folds

rmse.dataAddMult <- do.call(rbind,list(
  data.frame(rmse=addmult.dgp.cart.rmse,model="Single Tree (c.v.)")
  ,data.frame(rmse=addmult.dgp.rf.rmse,model="Random Forest (c.v.)")
  ,data.frame(rmse=addmult.dgp.gbm.rmse,model="GBM (c.v.)")
  ,data.frame(rmse=addmult.dgp.bart.rmse,model="BART")
  ,data.frame(rmse=addmult.dgp.krls.rmse,model="KRLS (c.v.)")
  ,data.frame(rmse=addmult.dgp.nn.rmse,model="Neural Net (c.v.)")
  ,data.frame(rmse=addmult.dgp.gam.rmse,model="GAM")
))


##Complicated DGP
#KRLS
cat("Working on KRLS comp DGP.\n")
comp.dgp.krls <- alply(comp.data.tr,3
                       ,cross.val
                       ,.data.test=comp.data.test
                       ,.formula=Y~.
                       ,.model="krls"
                       ,.params=krls.par
                       ,.parallel=TRUE
)
comp.dgp.krls.rmse <- laply(comp.dgp.krls,function(x)x$rmse)
comp.dgp.krls.model  <- comp.dgp.krls[[which.min(comp.dgp.krls.rmse)]]$mod.res

#Neural Net
cat("Working on NN comp DGP.\n")
comp.dgp.nn <- alply(comp.data.tr,3
                     ,cross.val
                     ,.data.test=comp.data.test
                     ,.formula=Y~.
                     ,.model="nnet.formula"
                     ,.params=nn.par
                     ,.parallel=TRUE
)
comp.dgp.nn.rmse <- laply(comp.dgp.nn,function(x)x$rmse)
comp.dgp.nn.model  <- comp.dgp.nn[[which.min(comp.dgp.nn.rmse)]]$mod.res



#Generalized additive model
cat("Working on GAM comp DGP.\n")
comp.dgp.gam <- alply(comp.data.tr,3
                      ,cross.val
                      ,.data.test=comp.data.test
                      ,.formula=as.formula(paste("Y~s(X1,k=20)+s(X2,k=20)+s(X10,k=20)+X3+X4+X5+X6+X7+X8+X9+",paste(colnames(add.data.tr)[-c(1:11)],sep="",collapse="+"),sep=""))
                      ,.model="gam"
                      ,.params=gam.par
                      ,.parallel=TRUE
)
comp.dgp.gam.rmse <- laply(comp.dgp.gam,function(x)x$rmse)
comp.dgp.gam.model  <- comp.dgp.gam[[which.min(comp.dgp.gam.rmse)]]$mod.res

#Single tree
cat("Working on cart comp DGP.\n")
comp.dgp.cart <- alply(comp.data.tr,3
                       ,cross.val
                       ,.data.test=comp.data.test
                       ,.formula=Y~.
                       ,.model="rpart"
                       ,.params=rpart.par
                       ,.parallel=TRUE
)
comp.dgp.cart.rmse <- laply(comp.dgp.cart,function(x)x$rmse)
comp.dgp.cart.model  <- comp.dgp.cart[[which.min(comp.dgp.cart.rmse)]]$mod.res


#CV & estmation Random Forest
cat("Working on RF comp DGP.\n")
comp.dgp.rf<- alply(comp.data.tr,3
                    ,cross.val
                    ,.data.test=comp.data.test
                    ,.formula=Y~.
                    ,.model="randomForest.formula"
                    ,.params=rf.par
                    ,.parallel=TRUE
)
comp.dgp.rf.rmse <- laply(comp.dgp.rf,function(x)x$rmse)
comp.dgp.rf.model  <- comp.dgp.rf[[which.min(comp.dgp.rf.rmse)]]$mod.res


#CV & estimation GBM
cat("Working on GBM comp DGP.\n")
comp.dgp.gbm<- alply(comp.data.tr,3
                     ,cross.val
                     ,.data.test=comp.data.test
                     ,.formula=Y~.
                     ,.model="gbm"
                     ,.params=gbm.par
                     ,.parallel=TRUE
)
comp.dgp.gbm.rmse <- laply(comp.dgp.gbm,function(x)x$rmse)
comp.dgp.gbm.model  <- comp.dgp.gbm[[which.min(comp.dgp.gbm.rmse)]]$mod.res



#CV & estimation BART
cat("Working on bart comp DGP.\n")
comp.dgp.bart <-  alply(comp.data.tr,3
                        ,cross.val
                        ,.data.test=comp.data.test
                        ,.formula=Y~.
                        ,.model="bart"
                        ,.params=bart.par
                        ,.parallel=TRUE
)
comp.dgp.bart.rmse <- laply(comp.dgp.bart,function(x)x$rmse)
comp.dgp.bart.model  <- comp.dgp.bart[[which.min(comp.dgp.bart.rmse)]]$mod.res




#RMSE across folds

rmse.dataComp <- do.call(rbind,list(
  data.frame(rmse=comp.dgp.cart.rmse,model="Single Tree (c.v.)")
  ,data.frame(rmse=comp.dgp.rf.rmse,model="Random Forest (c.v.)")
  ,data.frame(rmse=comp.dgp.gbm.rmse,model="GBM (c.v.)")
  ,data.frame(rmse=comp.dgp.bart.rmse,model="BART")
  ,data.frame(rmse=comp.dgp.krls.rmse,model="KRLS (c.v.)")
  ,data.frame(rmse=comp.dgp.nn.rmse,model="Neural Net (c.v.)")
  ,data.frame(rmse=comp.dgp.gam.rmse,model="GAM")
))


##Collect all rmse data.frames and plot using lattice
#rmse.dataMult$DGP <- "Multiplicative DGP"
rmse.dataAdd$DGP  <- "Additive DGP"
rmse.dataAddMult$DGP  <- "Additive and Multiplicative DGP"
rmse.dataComp$DGP  <- "Complicated DGP"
rmse.data <- do.call(rbind,list(
  rmse.dataAdd
  ,rmse.dataAddMult
  ,rmse.dataComp))

rmse.data$rel.rmse  <- ave(rmse.data$rmse,rmse.data$DGP,FUN=function(x)x/min(x)) 
rmse.data$model <- revalue(rmse.data$model, c("GAM"="GAM (informed)" ))
rmse.data$model  <- factor(rmse.data$model,levels=rev(c(
  "Neural Net (c.v.)"
  ,"Single Tree (c.v.)"
  ,"KRLS (c.v.)"
  ,"GAM (informed)"
  ,"BART"
  ,"Random Forest (c.v.)"
  ,"GBM (c.v.)"
)))

## Write results in csv
write.csv(rmse.data,file="rmsedata(finegrid).csv",row.names=FALSE)

## Save results of best crossvalidated models
save(add.dgp.cart.model,add.dgp.rf.model,add.dgp.gbm.model,add.dgp.bart.model
     ,addmult.dgp.cart.model,addmult.dgp.rf.model,addmult.dgp.gbm.model,addmult.dgp.bart.model
     ,comp.dgp.cart.model,comp.dgp.rf.model,comp.dgp.gbm.model,comp.dgp.bart.model
     ,file="BestModels(finegrid).RData")


## Prepare data for Figure 5
x10.sim <- seq(quantile(comp.data.test[,11],0)
               ,quantile(comp.data.test[,11],1)
               ,length.out=10)
x1.min <- quantile(comp.data.test[,2],0.25)
x1.max <- quantile(comp.data.test[,2],0.75)
x2.min <- quantile(comp.data.test[,3],0.25)
x2.max <- quantile(comp.data.test[,3],0.75)
sim.data.lo1.lo2 <- data.frame(Y = NA
                               ,X1=x1.min
                               ,X2=x2.min
                               ,X3=mean(comp.data.test[,4])
                               ,X4=mean(comp.data.test[,5])
                               ,X5=mean(comp.data.test[,6])
                               ,X6=0
                               ,X7=0
                               ,X8=1
                               ,X9=mean(comp.data.test[,10])
                               ,X10=x10.sim)
rem.vars <-  as.data.frame(t(replicate(10,colMeans(comp.data.test[,-c(1:11)]))))
sim.data.lo1.lo2 <- cbind(sim.data.lo1.lo2,rem.vars)
sim.data.lo1.hi2 <- mutate(sim.data.lo1.lo2,X2=x2.max)
sim.data.hi1.lo2 <- mutate(sim.data.lo1.lo2,X1=x1.max)
sim.data.hi1.hi2 <- mutate(sim.data.lo1.hi2,X1=x1.max)

all.sim.data <- rbind(sim.data.lo1.lo2,sim.data.lo1.hi2,sim.data.hi1.lo2,sim.data.hi1.hi2)


## Define and Execute prediction function
pred.df <- data.frame(model=c("rpart","randomForest.formula","gbm","bart")
                      ,inputs=c("comp.dgp.cart.model","comp.dgp.rf.model","comp.dgp.gbm.model","comp.dgp.bart.model")
                      ,stringsAsFactors = FALSE)

pred.fun <- function(model,inputs){
  return(cross.val(comp.data.test,all.sim.data,Y~.,as.character(model), getAnywhere(inputs)$objs[[1]],1,TRUE))
}

all.preds <- mdply(pred.df
                   ,.fun=pred.fun
                   ,.inform=TRUE)
all.preds$X1 <- rep(c("X1: Low","X1: High"),each=20) 
all.preds$X2 <- rep(c("X2: Low","X2: High"),each=10,length.out=40) 
all.preds$X10 <- rep(x10.sim,4)
all.preds <- rbind(data.frame(model="real"
                              ,inputs="none"
                              ,predictions=with(all.sim.data,ifelse(X10<2.5,c(X1 - X1^2 - X2^2  - 15*X1*X2*X10 + poly(X10,3,raw=TRUE)%*%c(10,-5,0.9)),1750+350*X10)
                              )
                              ,X1=rep(c("X1: Low","X1: High"),each=20) 
                              ,X2= rep(c("X2: Low","X2: High"),each=10,length.out=40)
                              ,X10=rep(x10.sim,4))
                   ,all.preds)

## Save effect predictions
save(all.preds,file="FXPreds.RData")

