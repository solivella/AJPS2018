############################
## File 7:
##    Code for analysis in
##    Section 5.3: replicating
##    and extendding GG 2013
##    using BART.
##
## Dependencies:
##    - File mapReduce_1.2.6.tar.gz
##    - MRP-20120715.RData
##    - census-pums-pop-2000-04-08.dat
##    - cps2000-04-08-DKs.dat
##    - votechoice2000-04-08.dat
##
## Note:
##    The first part of this file
##    along with the *.dat files
##    in the dependencies
##    were provided as part of the
##    replication materials for
##    Ghitza and Gelman 2013.
###################################

library(BayesTree)
library(plyr)
library(arm)

##FOR WINDOWS USERS:
## Prior to running lines 35--38
## download Rtools from https://cran.r-project.org/bin/windows/Rtools/
## in order to install mapReduce from source.
##ALL OTHER OS USERS:
##continue with code execution.
if(!require(mapReduce,quietly=TRUE))
{
  install.packages("mapReduce_1.2.6.tar.gz",repos=NULL,type="source")
  library(mapReduce)
}


set.seed(831213)

####Read data
load("MRP-20120715.RData")
census <- read.table("census-pums-pop-2000-04-08.dat"
                     ,sep="\t"
                     ,header=TRUE)
dat.stt <- read.table("state-stats.dat"
                      ,sep="\t"
                      ,header=TRUE)
data.to <- read.table("cps2000-04-08-DKs.dat"
                      ,sep="\t"
                      ,header=TRUE)
data.vo <- read.table("votechoice2000-04-08.dat"
                      ,sep="\t"
                      ,header=TRUE)


####Create auxiliary functions defined by Ghitza and Gelman
WMean <- function(a, w=rep(1,length(a)), subset=rep(TRUE,length(a))) {
  keep <- !is.na(a) & !is.na(w) & !is.na(subset) & subset
  return(sum((w*a)[keep])/sum(w[keep]))
}

FindDelta <- function(delta, a, w, x0)
  abs(x0-sum(invlogit(logit(a) + delta)*w))

CorrectWeighted <- function(a, w, x0) {
  delta <- optimize(FindDelta, interval=c(-5,5), a, w, x0)$minimum
  corrected <- invlogit(logit(a) + delta)
  return(list(delta=delta, corrected=corrected))
}

logit <- function (a) log(a/(1-a))

#### Replication
## Prepare data
#CPS
data.to.full <- M.cps[[2]]@frame
unobs.ind.to <- which(rowSums(data.to.full[,1])==0)
y.to <- data.to.full[-unobs.ind.to,1]
y.to <- round(y.to)
data.to.agg <- data.to.full[-unobs.ind.to,c(6:9)]
data.to.1 <- cbind(1,do.call(cbind,alply(data.to.agg,2,function(x)rep(unlist(x),y.to[,1]))))
data.to.0 <- cbind(0,do.call(cbind,alply(data.to.agg,2,function(x)rep(unlist(x),y.to[,2]))))
data.to.sub <- rbind(data.to.1,data.to.0)
colnames(data.to.sub) <- c("turnout", "stt","eth","inc","age")

#Vote intention data
data.vo.full <- M.vot[[2]]@frame
unobs.ind.vo <- which(rowSums(data.vo.full[,1])==0)
y.vo <- data.vo.full[-unobs.ind.vo,1]
y.vo <- round(y.vo)
data.vo.agg <- data.vo.full[-unobs.ind.vo,c(6:9)]
data.vo.1 <- cbind(1,do.call(cbind,alply(data.vo.agg,2,function(x)rep(unlist(x),y.vo[,1]))))
data.vo.0 <- cbind(0,do.call(cbind,alply(data.vo.agg,2,function(x)rep(unlist(x),y.vo[,2]))))
data.vo.sub <- rbind(data.vo.1,data.vo.0)
colnames(data.vo.sub) <- c("vote", "stt","eth","inc","age")


data.to.sub <- numcolwise(factor)(as.data.frame(data.to.sub))
data.vo.sub <- numcolwise(factor)(as.data.frame(data.vo.sub))

#Census data
census_rep <- subset(census
                 ,(inc!=-1&stt!=-1&eth!=-1&age!=-1)
                 ,select=c(stt
                           ,eth
                           ,inc
                           ,age
                           ,wtd2008
                 ))
census.sub_rep <- aggregate(census_rep[,5]
                        ,census_rep[,1:4]
                        ,FUN=sum)
census.sub_rep <- subset(census.sub_rep,stt!=2)

# state-level data
dat.stt$z.inc2000 <- rescale(dat.stt$inc2000)
dat.stt$z.inc2004 <- rescale(dat.stt$inc2004)
dat.stt$z.inc2007 <- rescale(dat.stt$inc2007)
dat.stt$z.rep1996 <- rescale(dat.stt$rep1996)
dat.stt$z.rep2000 <- rescale(dat.stt$rep2000)
dat.stt$z.rep2004 <- rescale(dat.stt$rep2004)
dat.stt$z.rep2008 <- rescale(dat.stt$rep2008)
dat.stt$z.trn1996 <- rescale(dat.stt$vote1996/dat.stt$pop1996)
dat.stt$z.trn2000 <- rescale(dat.stt$vote2000/dat.stt$pop2000)
dat.stt$z.trn2004 <- rescale(dat.stt$vote2004/dat.stt$pop2004)
dat.stt$z.trn2008 <- rescale(dat.stt$vote2008/dat.stt$pop2007)
dat.stt$stt <- 1:nrow(dat.stt)

##Estimate BARTs
#Turnout
mrp.bart.data.to <- makeind(data.to.sub[,-1])
mrp.bart.sim.to <- expand.grid(stt=1:51
                               ,eth=1:4
                               ,inc=1:5
                               ,age=1:4)
mrp.bart.sim.to <- as.data.frame(numcolwise(factor)(mrp.bart.sim.to))
mrp.bart.sim.full.to <- makeind(mrp.bart.sim.to)


system.time(mrp.bart.to <- bart(x.train=mrp.bart.data.to
                                ,y.train=data.to.sub[,1]
                                ,x.test=mrp.bart.sim.full.to
                                ,binaryOffset=-0.5))

save(mrp.bart.to,file="MRP_BART_TO_Rep.RData")

mrp.bart.meds.to <- colMeans(mrp.bart.to$yhat.test)
mrp.bart.sim.to$preds <-  mrp.bart.meds.to

#Vote Intention
mrp.bart.data.vo <- makeind(data.vo.sub[,-1])
mrp.bart.sim.vo <- expand.grid(stt=c(1:51)[-2]
                               ,eth=1:4
                               ,inc=1:5
                               ,age=1:4)
mrp.bart.sim.vo <- as.data.frame(numcolwise(factor)(mrp.bart.sim.vo))
mrp.bart.sim.full.vo <- makeind(mrp.bart.sim.vo)



system.time(mrp.bart.vo <- bart(x.train=mrp.bart.data.vo
                                ,y.train=data.vo.sub[,1]
                                ,x.test=mrp.bart.sim.full.vo))
save(mrp.bart.vo,file="MRP_BART_VO_Rep.RData")

mrp.bart.meds.vo <- colMeans(mrp.bart.vo$yhat.test)
mrp.bart.sim.vo$preds <-  mrp.bart.meds.vo


##Full predicted dataset
mrp.bart.sim  <- subset(mrp.bart.sim.to,stt!=2)

mrp.bart.sim$predsto <- pnorm(mrp.bart.sim$preds)
mrp.bart.sim$predsvo <- pnorm(mrp.bart.meds.vo)

census.sub_rep$predsto  <- mrp.bart.sim$predsto
census.sub_rep$predsvo <- mrp.bart.sim$predsvo

## Correct based on actual turnout and vote estimates at state level
census.sub_rep$vote2008 <- census.sub_rep$turn2008 <- NA
for (i in 1:51) {
  cat(paste(i, "of", 51, "\n"))
  ok <- census.sub_rep$stt==i
  census.sub_rep$turn2008[ok] <- CorrectWeighted(a=census.sub_rep$predsto[ok]
                                             , w=census.sub_rep$x[ok]
                                             , x=dat.stt[as.character(i)
                                                         , "vote2008"])$corrected
  census.sub_rep$vote2008[ok] <- CorrectWeighted(a=census.sub_rep$predsvo[ok]
                                             , w=census.sub_rep$x[ok]/sum(census.sub_rep$x[ok])
                                             , x=dat.stt[as.character(i)
                                                         , "rep2008"])$corrected
}


##Comparison between models BART and MRP models
D$grp.plot <- apply(D[, c("stt", "eth", "inc", "age")], 1, paste, collapse="_")
D.plot <- as.data.frame(mapReduce(data=D, map=grp.plot,
                                  stt=unique(stt),
                                  eth=unique(eth),
                                  inc=unique(inc),
                                  age=unique(age),
                                  vote2008=sum(vote2008*turn2008*pop2008) / sum(turn2008*pop2008),
                                  turn2008=sum(turn2008*pop2008) / sum(pop2008),
                                  pop2008=sum(pop2008)))
D.plot <- D.plot[with(D.plot,order(stt,eth,inc,age)),]
D.plot <- subset(D.plot,stt!=2)

#Correlation of predicted values:
round(cor(D.plot$turn2008,sttethincagepos$predsto,use="complete"),3)

round(cor(D.plot$vote2008,sttethincagepos$predsvo,use="complete"),3)


#### Extension
##Prepare data

#Turnout
data.to.sub <- subset(data.to,file=="2008-cps"&is.na(vote)==FALSE
                      ,select=c(vote
                                ,stt
                                ,eth
                                ,inc
                                ,age
                                ,sex
                                ,edu
                                ,mar
                                ,kid
                      )
)

#Vote intention
data.vo.sub <- subset(data.vo,file=="2008-pew"&is.na(rvote)==FALSE
                      ,select=c(rvote
                                ,stt
                                ,eth
                                ,inc
                                ,age
                                ,sex
                                ,edu
                                ,mar
                                ,kid
                      )
)
data.to.sub <- data.frame(vote=data.to.sub[,1],numcolwise(factor,exclude=FALSE)(data.to.sub[,-1]))
data.vo.sub <- data.frame(rvote=data.vo.sub[,1],numcolwise(factor,exclude=FALSE)(data.vo.sub[,-1]))

#Census data
census_ext <- subset(census
                 ,(inc!=-1&stt!=-1&eth!=-1&age!=-1&sex!=-1&edu!=-1&mar!=-1&kid!=-1)
                 ,select=c(stt
                           ,eth
                           ,inc
                           ,age
                           ,sex
                           ,edu
                           ,mar
                           ,kid
                           ,wtd2008
                 ))

census.sub_ext <- aggregate(census_ext[,"wtd2008"]
                        ,census_ext[,1:8]
                        ,FUN=sum)
census.sub_ext <- subset(census.sub_ext,stt!=2)


##Estimate BARTs of extension
#Turnout
mrp.bart.data.to <- makeind(data.to.sub[,-1])
mrp.bart.sim.to <- expand.grid(stt=1:51
                               ,eth=c(1:4)
                               ,inc=c(1:5,NA)
                               ,age=c(1:4)
                               ,sex=1:2
                               ,edu=c(1:5)
                               ,mar=c(1:2)
                               ,kid=c(1:2,NA))
mrp.bart.sim.to <- as.data.frame(numcolwise(factor
                                            ,exclude=FALSE)(mrp.bart.sim.to))
mrp.bart.sim.full.to <- makeind(mrp.bart.sim.to)

system.time(mrp.bart.to <- bart(x.train=mrp.bart.data.to
                                ,y.train=data.to.sub[,1]
                                ,x.test=mrp.bart.sim.full.to
                                ,power=1.5))
save(mrp.bart.to,file="MRP_BART_TO_Exp.RData")

mrp.bart.meds.to <- colMeans(mrp.bart.to$yhat.test)
mrp.bart.sim.to$preds <-  mrp.bart.meds.to

#Vote Intention
mrp.bart.data.vo <- makeind(data.vo.sub[,-1])
mrp.bart.sim.vo <- expand.grid(stt=c(1:51)[-2]
                               ,eth=c(1:4,NA)
                               ,inc=c(1:5,NA)
                               ,age=c(1:4,NA)
                               ,sex=1:2
                               ,edu=c(1:5,NA)
                               ,mar=c(1:2,NA)
                               ,kid=c(1:2,NA))
mrp.bart.sim.vo <- as.data.frame(numcolwise(factor
                                            ,exclude=FALSE)(mrp.bart.sim.vo))
mrp.bart.sim.full.vo <- makeind(mrp.bart.sim.vo)

system.time(mrp.bart.vo <- bart(x.train=mrp.bart.data.vo
                                ,y.train=data.vo.sub[,1]
                                ,x.test=mrp.bart.sim.full.vo
                                ,power=1.5))#From first system time.
save(mrp.bart.vo,file="MRP_BART_VO_Exp.RData")

mrp.bart.meds.vo <- colMeans(mrp.bart.vo$yhat.test)
mrp.bart.sim.vo$preds <-  mrp.bart.meds.vo

##Full predicted dataset
mrp.bart.sim  <- subset(mrp.bart.sim.to,stt!=2,drop=TRUE)
names(mrp.bart.sim)[9] <- "predsto"

mrp.bart.sim <- merge(mrp.bart.sim,mrp.bart.sim.vo)
names(mrp.bart.sim)[10] <- "predsvo"

census.sub_ext  <- merge(census.sub_ext,mrp.bart.sim)
census.sub_ext$predsto <- pnorm(census.sub_ext$predsto-0.5)
census.sub_ext$predsvo <- pnorm(census.sub_ext$predsvo)


## Correct based on actual state-level turnout and vote estimates
census.sub_ext$vote2008 <- census.sub_ext$turn2008 <- NA
for (i in 1:51) {
  cat(paste(i, "of", 51, "\n"))
  ok <- census.sub_ext$stt==i
  census.sub_ext$turn2008[ok] <- CorrectWeighted(a=census.sub_ext$predsto[ok], w=census.sub_ext$x[ok], x=dat.stt[as.character(i), "vote2008"])$corrected
  census.sub_ext$vote2008[ok] <- CorrectWeighted(a=census.sub_ext$predsvo[ok], w=census.sub_ext$x[ok]/sum(census.sub_ext$x[ok]), x=dat.stt[as.character(i), "rep2008"])$corrected
}

#Write results of prediction exercise
write.csv(census.sub_ext, "MRP_BART_Expansion_Preds.csv",row.names=FALSE)


