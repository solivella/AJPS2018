#######################################
## File 5:
##    Figure 4 and SI-1
##
## Dependencies: 
##  - File (4)Analysis_section_5_1.R
##      (see file notes in Readme.txt)
##  - File rmsedata(finegrid).csv
##  - File SynthData.RData
##  - File FXPreds.RData
#######################################


library(lattice)
library(latticeExtra)

rmse.data <- read.csv("rmsedata(finegrid).csv")
load("SynthData.RData")
load("FXPreds.RData")

## Figure 4
bwplot(model~rel.rmse|DGP
       ,data=rmse.data
       ,DGP=rmse.data$DGP
       ,index.cond=list(c(2,1,3))
       ,rmse = rmse.data$rmse
       ,layout=c(3,1)
       ,xlab="Relative RMSE"
       ,pch="|"
       #,xlim=c(-0.1,2.5)
       ,scales=list(x=list(relation="free"
                           ,limits=list(c(1,1.5),c(1,1.4),c(1,1.42))
       ))
       ,par.settings=standard.theme("pdf", color=FALSE)
       ,panel = function(x, y, rmse, DGP,subscripts, ...) { 
         panel.bwplot(x, y, do.out=FALSE, fill=c(rep("gray80",3),rep(NA,2),"gray80",NA),...) 
         meds <- tapply(x, y, median) 
         ylocs <- seq_along(meds) 
         panel.segments(meds, ylocs - 1/4, 
                        meds, ylocs + 1/4, 
                        lwd = 3, col = "black")
         switch (as.character(unique(DGP[subscripts])),
                 "Additive DGP" = panel.abline(v=sqrt(mean((add.data.test[,1]-mean(add.data.test[,1]))^2))/min(rmse[subscripts]),lty=2),
                 "Additive and Multiplicative DGP" = panel.abline(v=sqrt(mean((addmulti.data.test[,1]-mean(addmulti.data.test[,1]))^2))/min(rmse[subscripts]),lty=2),
                 "Complicated DGP" = panel.abline(v=sqrt(mean((comp.data.test[,1]-mean(comp.data.test[,1]))^2))/min(rmse[subscripts]),lty=2)
         )
         #panel.abline(v=1,lty=3)
       }
)


## Figure SI-1
myPlot <- xyplot(predictions~X10|as.factor(X1)+as.factor(X2)
                 ,groups=model
                 ,data=all.preds
                 #,scales=list(y="free")
                 #,ylim=list(c())
                 ,key=list(title="Model"
                           ,space="right"
                           ,lines=list(lwd=c(3.5,rep(2,4))
                                       ,col=c("black"
                                              ,"#377eb8"
                                              ,"#4daf4a"
                                              ,"#984ea3"
                                              ,"#ff7f00"
                                       ))
                           ,points=list(cex=c(1.5,rep(1.1,4)),pch=c(19,6,0,5,1))
                           ,text=list(c("Real DGP","Single Tree (c.v)","Random Forest (c.v.)","GBM (c.v)","BART")))
                 ,par.settings=standard.theme("pdf",color=FALSE)
                 ,ylab="Outcome"
                 ,xlab="Predictor X10"
                 ,main=""
                 ,panel=function(x,y,groups,subscripts,...){
                   panel.superpose(x,y,groups=groups,subscripts=subscripts
                                   ,panel.groups=function(x,y,col.line,type,lwd,...,lty,cex,pch){
                                     panel.xyplot(x,y,type='b'
                                                  ,lty=1
                                                  ,lwd=c(3.5,rep(2,4))
                                                  ,cex=c(1.5,rep(1.1,4))
                                                  ,groups=groups
                                                  ,pch=c(19,6,0,5,1)
                                                  ,col.line=c("black"
                                                              ,"#377eb8"
                                                              ,"#4daf4a"
                                                              ,"#984ea3"
                                                              ,"#ff7f00"
                                                  ),...)
                                     
                                   }
                   )
                 }
                 
)
useOuterStrips(myPlot)