#####################################
## File 7
##    Figure 5
##
## Dependencies:
##    - (7)Analysis_sect_5_3.R
##    - MRP_BART_Expansion_Preds.csv
##
######################################

library(plyr)
library(lattice)
library(latticeExtra)

##Load data
census_sub_ext <- read.csv("MRP_BART_Expansion_Preds.csv")

## Figure 7
par(mfrow=c(1,2))
colors.edu2 <- rev(gray.colors(5,0.1,0.8,alpha=0.4))

## Post-stratified
##State X Ethinicity x income x age x sex x edu
withsexedu <- ddply(census_sub_ext
                    ,c("stt","eth","inc","age","sex","edu")
                    ,function(x){
                      dento <- sum(x[,"x"])
                      denvo <- sum(x[,"x"]*x[,"turn2008"])
                      numto <- sum(x[,"turn2008"]*x[,"x"])
                      numvo  <- sum(x[,"vote2008"]*x[,"turn2008"]*x[,"x"])
                      return(data.frame(predsto=numto/dento,predsvo=numvo/denvo,size=dento))
                    }) 
withsexedu$eth <- factor(withsexedu$eth,labels=c("White","Black","Latino","Other"))
withsexedu$sex <- factor(withsexedu$sex,labels=c("Male","Female"))

withsexedu <- withsexedu[sample(1:nrow(withsexedu),nrow(withsexedu)),]

myPlot <- xyplot(predsto~predsvo|sex+eth
                 ,data=withsexedu
                 ,ylim=c(0,1)
                 ,xlim=c(-0.05,1)
                 ,size=withsexedu$size
                 ,par.settings=standard.theme("pdf",color=FALSE)
                 ,ylab="Turnout"
                 ,xlab="McCain Vote '08"
                 ,main="State x Ethnicity x Income x Age x Sex x Edu"
                 ,panel=function(x,y,subscripts
                                 ,size,...){
                   panel.xyplot(x,y
                                ,col=colors.edu2[withsexedu$edu[subscripts]]
                                ,cex=sqrt(size[subscripts])/50
                                ,pch=19
                   )
                   panel.abline(h=0.5,col="gray70",lty=2)
                   panel.abline(v=0.5,col="gray70",lty=2)
                   
                 }
)
useOuterStrips(myPlot)
