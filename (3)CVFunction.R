###Cross validation function
cross.val <- function( .data.train
                       ,.data.test
                       ,.formula
                       ,.model
                       ,.params
                       ,.folds=5
                       ,.predict=FALSE
                       #,.parallel=FALSE
                       #,.n.cores=22
){
  set.seed(831213)
  switch(.model
         ,"krls"=library(KRLS)
         ,"gam"=library(mgcv)
         ,"nnet.formula"=library(nnet)
         ,"rpart"=library(rpart)
         ,"randomForest.formula"=library(randomForest)
         ,"gbm"=library(gbm)
         ,"bart"=library(BayesTree)
  )
  # .modelName <- deparse(substitute(.model)) #Get the name of the model, as string
  .modelName <- .model
  .model <- getAnywhere(.model)$objs[[1]] #catch the actual function object of the model; see body of apply().
  .data.train  <- as.data.frame(.data.train)
  .data.test  <- as.data.frame(.data.test)
  .params$fold <- 1:.folds
  param.grid <- do.call(expand.grid,.params) #get all possible paramater combinations
  
  #The following splits the data into .folds groups,
  #initializes train=0 flags, cycles through groups
  #and produces a list with data.frames with train=1 for members of 
  #the corresponding group, and train=0 for all other obs.
  .data.test$train <- 0
  if(.predict){
    .data.train$train <- 1
    data.folds <- list(.data.train)
  }else{
    fold.index  <- sample(as.numeric(cut(1:nrow(.data.train),.folds)))
    .data.train$train  <- NA
    data.folds  <- llply(unique(fold.index)
                         ,function(x,data=.data.train){
                           temp.data  <- data
                           temp.data[,"train"]  <- as.numeric(fold.index!=x)
                           return(temp.data)
                         })
  }
  
  #This will generate predictions based on an estimated model object 
  pred.function  <- function(model.obj,prediction=FALSE){
    
    model.preds <- switch(model.obj$modelName #different predict methods for different models
                          #form adequate test set predictions for each model
                          ,bart = model.obj$results$yhat.test.mean
                          ,gbm =  predict(model.obj$results,newdata=model.obj$data.test,n.trees=model.obj$inputs$n.trees)
                          ,krls = predict(model.obj$results,newdata=as.matrix(model.obj$data.test[,!grepl("Y",colnames(model.obj$data.test))]))$fit
                          ,predict(model.obj$results,newdata=model.obj$data.test) #default value.
    )
    rmse <- sqrt(mean((model.preds-model.obj$data.test[,"Y"])^2)) #calculate rmse
    if(prediction)
    {
      return(data.frame(preds=model.preds))
    }else{
      return(data.frame(rmse=rmse)) 
    }
  }
  
  #This will form the argument list and execute the model function:
  model.wrapper <- function(inputs,.data.folds=data.folds,mod.name=.modelName,mod.fun=.model,form=.formula,pred=FALSE){ 
    set.seed(831213)
    inputs <- llply(inputs,function(x)x) #create model list
    train.data  <-  subset(.data.folds[[inputs$fold]],train==1,select=-c(train)) #get correct train data (kth fold)
    test.data <- subset(.data.folds[[inputs$fold]],train==0,select=-c(train))
    inputs$formula <- form
    outcome.name <- as.character(form)[2]
    outcome.index  <- which(colnames(test.data)==outcome.name)
    
    
    inputs$data <- train.data #Create data argument using train.data
    
    if(mod.name=="bart"){ #call to bart is different, so form special arguments     
      inputs$x.train <- train.data[,-outcome.index]
      inputs$y.train  <- train.data[,outcome.index]
      inputs$x.test <- test.data[,-outcome.index]
    }
    if(mod.name=="rpart"){ #call to rpart takes parameters as list named `control'  
      inputs$control  <- inputs[-grep("formula|data",names(inputs))]
    }
    if(mod.name=="krls"){ #and call to krls is different also
      #cat("\nI'm here!\n",file="modname.txt",append=TRUE)
      inputs$X  <- as.matrix(train.data[,-outcome.index])
      inputs$y  <- as.matrix(train.data[,outcome.index])
    }
    #    cat("\n\t",.modelName,"\t\n",file="modname.txt",append=TRUE)
    #Catch unused argument error, and correct:
    #cat("\n-------------\n\n",file="data0.txt",append=TRUE)
    #cat("\n",names(inputs),"\n\n",file="data0.txt",append=TRUE)
    model.res <- tryCatch(do.call(mod.fun,inputs) #estimate model passing arguments as list
                          ,error=function(x){
                            err.msg  <- conditionMessage(x)
                            if(grepl("argument|invalid",err.msg,ignore.case=TRUE)){ #Won't work if errors are given in locale other than English
                              inputs <- inputs[names(inputs)%in%names(formals(mod.fun))]
                              #cat("\n-------------\n\n",file="data1.txt",append=TRUE)
                              #cat("\n",names(inputs),"\n\n",file="data1.txt",append=TRUE)
                              return(do.call(mod.fun,inputs))
                            }else{stop(x)}
                          })
    model.obj.res <- list(results=model.res,modelName=mod.name,data.test=test.data,inputs=inputs)
    
    pred.res <- pred.function(model.obj.res,pred)
    if(pred){
      return(unlist(pred.res$preds)) 
    }else{
      return(c(rmse=pred.res$rmse))
    }
  }
  
  all.data <- list(data.folds[[1]])
  all.data[[1]]$train = 1
  all.data[[1]] <- rbind(all.data[[1]],.data.test)
  
  #Conduct parameter sweep, and return list
  #with parameters and rmse for each parameter combination and fold:
  if(!.predict){
    if(.modelName!="krls"&&.modelName!="gam"&&.modelName!="bart")
    {
      cv.results.rmse <- adply(param.grid
                               ,1
                               ,model.wrapper
      )
      
      un.param.index <- grep("^fold$|^rmse$",names(cv.results.rmse)) #get index of things that aren't parameters
      
      #Find best model inputs given train data
      #Get median rmse per param combination across data folds:
      res.temp <- ddply(cv.results.rmse
                        ,.variables=names(cv.results.rmse)[-un.param.index]
                        ,.fun=function(x)c(ave.rmse=median(x$rmse)))
      best.inputs <- as.list(res.temp[which.min(res.temp$ave.rmse),-un.param.index])
    }else{
      best.inputs = as.list(param.grid[1,])
    }
    best.inputs$fold <- 1
    #Fit best model to entire train data
    oos.preds <- model.wrapper(best.inputs,all.data)
  }else{
    #Fit best model to entire train data
    oos.preds <- model.wrapper(.params,all.data,pred=TRUE)
  }
  
  
  #Obtain out-of-sample RMSE for best model given train data 
  #oos.preds  <- pred.function(res.temp.obj)
  if(.predict){
    return(data.frame(predictions=oos.preds))
  }else{
    return(list(rmse=oos.preds,mod.res=best.inputs))
  }
}
