###########################
## File 2:
##    Create data used in 
##    section 5.1
##
## Dependencies: none
##  
###########################


library(plyr)
library(MASS)


set.seed(831213)

## Function to create datasets:
data.gen  <- function(DGP,n=500,w.rand.error=TRUE,mean.vec){
  X1 <- rgamma(n,8,2)
  X2 <- rgamma(n,10,1)
  X3.5 <-  mvrnorm(n,c(2,3,6),diag(c(1.5,0.5,3.3)))
  X6.8 <- t(rmultinom(n,size = 1,prob = c(1/3,1/3,1/3)))
  X9.10 <- mvrnorm(n,mu = c(-0.3,2),Sigma = matrix(c(1.5,0.685,0.685,5.5),2))  
  X11.40 <- mvrnorm(n,mu = mean.vec,Sigma=diag(30))
  
  ### Create outcomes
  Y  <- switch(DGP
               ,"additive" = X1  +  X2 + X3.5[,1] + X9.10[,2]
               ,"multiplicative" =  X1*X2*X3.5[,1]*X9.10[,2] 
               ,"add+multip" = 1.5 + (3.2)*X1 + 0.1*X2 - (2.8)*X9.10[,2] + 1.2*X1*X2 -
                 1.2*X1*X9.10[,2] - 0.2*X2*X9.10[,2] + 2.4*X1*X2*X9.10[,2]
               ,"complicated" =  ifelse(X9.10[,2]<2.5,c(X1 - X1^2 - X2^2  - 15*X1*X2*X9.10[,2] + poly(X9.10[,2],3,raw=TRUE)%*%c(10,-5,0.9)),1750+350*X9.10[,2])
  )
  if(w.rand.error){
    Y  <- Y + sd(Y)*rnorm(n)
  }
  
  ### Create dataset
  temp.data  <-  data.frame(Y=Y
                            ,X1=X1
                            ,X2=X2
                            ,X3=X3.5[,1]
                            ,X4=X3.5[,2]
                            ,X5=X3.5[,3]
                            ,X6=X6.8[,1]
                            ,X7=X6.8[,2]
                            ,X8=X6.8[,3]
                            ,X9=X9.10[,1]
                            ,X10=X9.10[,2])
  temp.data <- cbind(temp.data,as.data.frame(X11.40))
  return(as.matrix(temp.data))
}

#Common mean vector for irrelevant variables
irr_mean_vec <- sample(2:10, 30, replace=TRUE)

## 100 train datasets
#Additive data
add.data.tr  <- raply(100,.expr = data.gen("additive",mean.vec=irr_mean_vec))
add.data.tr  <- aperm(add.data.tr,c(2,3,1))

#Additive+Multiplicative data
addmulti.data.tr  <- raply(100,.expr = data.gen("add+multip",mean.vec=irr_mean_vec))
addmulti.data.tr  <- aperm(addmulti.data.tr,c(2,3,1))

#Complicated data
comp.data.tr  <- raply(100,.expr = data.gen("complicated",mean.vec=irr_mean_vec))
comp.data.tr  <- aperm(comp.data.tr,c(2,3,1))

# Common test sets
add.data.test  <- data.gen("additive",n=3e3,mean.vec=irr_mean_vec)
addmulti.data.test  <- data.gen("add+multip",n=3e3,mean.vec=irr_mean_vec)
comp.data.test  <- data.gen("complicated",n=3e3,mean.vec=irr_mean_vec)

## Save data
save(add.data.test,addmulti.data.test,comp.data.test
     ,add.data.tr,addmulti.data.tr,comp.data.tr
     ,file="SynthData.RData")
