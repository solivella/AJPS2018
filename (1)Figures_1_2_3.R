##############################
## File 1:
## Figures 1, 2 and 3
##############################

library(lattice)
library(grid)
library(latticeExtra)
library(rpart)
library(caret)

set.seed(234)
mytheme <- standard.theme("pdf",color=FALSE)
X1 <- sample(1:10, 600, replace=TRUE)
X2 <- sample(1:10, 600, replace=TRUE)

## Figure 1
freq_table <- xtabs(~X2+X1)
levelplot(freq_table
          ,col.regions = gray(15:0/15)
          ,colorkey=list(at=0:13))
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text("Obs.\nFrequency", 0.2, 1.1, hjust=0.5, vjust=1)
trellis.unfocus()


y <- (sqrt(X1*X2)*(5.5-X1)^2/((5.5-X1)*(5.5-X2)))/100
y.norm<- y + rnorm(600, 0, .35)
y.pred1 <- predict(lm(y.norm~X1+X2))
y.pred2 <- predict(lm(y.norm~poly(X1,degree=2,raw=TRUE)*poly(X2,degree=2,raw=TRUE)))#+X1*X2+X1^2+X2^2))
sim_data <- data.frame(y.norm=y.norm,X1=X1,X2=X2)
this <- data.frame(rbind(aggregate(y, list(X1, X2), mean),
                         aggregate(y.norm, list(X1, X2), mean),
                         aggregate(y.pred1, list(X1, X2), mean),
                         aggregate(y.pred2, list(X1, X2), mean)))
this$type <- rep(c("True data generating process","Naive nonparametric model",
                   "Simple linear model", "Complex linear model"),each=99)
colnames(this)[1:3] <- c("X1","X2","y")

## Figure 2
wireframe(y~X1*X2|type, this, xlim=c(1,10), ylim=c(1,10), zlim=c(-.2,.7), drape=TRUE, zlab="y",
          scales = list(arrows=FALSE, col="black",font=10, cex=.5), auto.key=FALSE,
          par.settings=mytheme, colorkey=FALSE , layout=c(2,2), index.cond=list(c(3,1,4,2)),
          xlab=expression(x[1]), ylab=expression(x[2]))

## Figure 3
set.seed(831213)
t_ctrl <- trainControl("cv",15)
param_grid <- expand.grid(.cp=seq(0.001,0.01,length.out=10))
t_rpart <- train(y.norm~X1+X2
                 ,data=sim_data
                 ,method="rpart"
                 ,trControl=t_ctrl
                 #,control=list(minsplit=20)
                 ,tuneGrid=param_grid)

#left panel
regions <- xyplot(1:10~1:10
       ,ylim=c(1,10)
       ,xlim=c(1,10)
       ,xlab=expression(X[1])
       ,ylab=expression(X[2])
       ,par.settings = list(axis.line = list(col = 0))
       ,scales=list(col=1,tck=c(1,0),x=list(at=1:10,cex=0.7),y=list(at=1:10,cex=0.7))
       ,panel = function(x,y,subscripts,...){
         panel.abline(v=1,h=1,lwd=1.2)
         panel.segments(x0=c(8.5,8.5,8.5,8.5,8.5,2.5,1,2.5,2.5,3.5,3.5,5.5,5.5), y0=c(1,5.5,6.5,7.5,4.5,1,1.5,6.5,4.5,4.5,5.5,4.5,5.5),
                        x1=c(8.5,10,10,10,10,2.5,2.5,8.5,8.5,3.5,8.5,5.5,5.5), y1=c(10,5.5,6.5,7.5,4.5,10,1.5,6.5,4.5,6.5,5.5,5.5,6.5),
                        lwd=1.5,col="black")
                        
       })

#central panel
node_locs <- plot(t_rpart$finalModel,compress=TRUE,uniform=TRUE,margin=0.05)
tree_lines <-  rpart:::rpart.branch(node_locs$x,node_locs$y
                                    ,as.numeric(row.names(t_rpart$finalModel$frame))
                                    ,1)
panel.text.mine <- function(x,y,labels,...){
  labels <- ifelse(is.na(labels),"",labels)
  panel.text(x,y,labels,...)
}
tree <- xyplot(y~x
       ,data=node_locs
       ,xlab=NULL
       ,ylab=NULL
       ,scales=list(draw=FALSE)
       ,par.settings = list(axis.line = list(col = 0))
       ,tree_lines = tree_lines
       ,tree_model = t_rpart$finalModel
       ,panel=function(x,y,tree_lines,tree_model,...){
         lab_orig <- labels(tree_model)
         nodes <- as.numeric(row.names(tree_model$frame))
         left_childs <-  match(2L * nodes, nodes)
         panel.lines(tree_lines$x,tree_lines$y,lwd=1.5,col="black")
         rpart:::text.rpart(t_rpart$finalModel,digits=2
                            ,FUN=panel.text.mine
                            ,cex=0.6
                            ,fancy=FALSE
                            )
         
         #panel.text()
       })


#right panel
pred_data <- expand.grid(X1=1:10,X2=1:10)
pred_data$z_vals <- predict(t_rpart,pred_data)
surf_pred <- wireframe(z_vals~X1*X2, pred_data,xlim=c(1,10), ylim=c(1,10),zlim=c(-.2,.7), drape=TRUE, zlab=expression(hat(y)),
          scales = list(arrows=FALSE, col="black",font=10, cex=.5,x=list(at=1:10),y=list(at=1:10)), auto.key=FALSE,
          par.settings=mytheme, colorkey=FALSE,
          xlab=expression(X[1]), ylab=expression(X[2]))
figure3 <- c(regions,tree,surf_pred,layout=c(3,1))
update(figure3,xlab=list(label=c(expression(X[1]),"","")),ylab=expression(X[2])
       ,scales=list(at=list(c(1:10),NULL,NULL))
       )
