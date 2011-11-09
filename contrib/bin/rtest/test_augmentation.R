library(qtl)

set.seed(1000)              

version = mqm_version()
cat("R/qtl=",version$RQTL)
cat("R-MQM=",version$RMQM)
cat("MQM=",version$MQM)

testaugmentation <- function(cross, ...){
  crossML <- mqmaugment(cross, ...)
  crossIMP <- mqmaugment(cross,strategy="impute", ...)
  crossDROP <- mqmaugment(cross,strategy="drop", ...)

  res1 <- mqmscan(crossML)
  res2 <- mqmscan(crossIMP)
  res3 <- mqmscan(crossDROP)

  plot(res1,res2,res3,lty=c(1,2,3),col=c("black","blue","red"))
  legend("topleft",c("MinProb","Imputation","Drop"),lty=c(1,2,3),col=c("black","blue","red"))
  list(res1,res2,res3)
}

stabilitytest <- function(cross, ...){
  crossML1 <- mqmaugment(cross, minprob=0.5,maxaugind=16,...)
  crossML2 <- mqmaugment(cross, minprob=0.05,maxaugind=16,...)
  crossML3 <- mqmaugment(cross, minprob=0.005,maxaugind=16,...)
  res1 <- mqmscan(crossML1)
  res2 <- mqmscan(crossML2)
  res3 <- mqmscan(crossML3)

  plot(res1,res2,res3,lty=c(1,2,3),col=c("black","blue","red"))
  legend("topleft",c("minprob 0.5","minprob 0.05","minprob 0.005"),lty=c(1,2,3),col=c("black","blue","red"))
  list(res1,res2,res3)
}

data(multitrait)
multimissing <- simulatemissingdata(multitrait,25)
#r <- stabilitytest(multimissing,verbose=TRUE)
r <- testaugmentation(multimissing)
if(!round(r[[1]][3,3],3)==0.764) stop("Multitrait ML dataaugmentation error")
if(!round(r[[1]][3,3],3)==round(r[[2]][3,3],3)) stop("Multitrait ML compared versus IMP error")
if(!round(r[[3]][3,3],3)==0.844) stop("Multitrait DROP dataaugmentation error")

data(hyper)
r <- testaugmentation(hyper,maxaugind=32,minprob=0.75)
if(!round(r[[1]][3,3],3)==0.568) stop("Hyper ML dataaugmentation error")
if(!round(r[[1]][3,3],3)==round(r[[2]][3,3],3)) stop("Hyper ML compared versus IMP error")
if(!round(r[[3]][3,3],3)==0.252) stop("Hyper DROP dataaugmentation error")


data(listeria)
r <- testaugmentation(listeria)
if(!round(r[[1]][3,3],3)==0.352) stop("Listeria ML dataaugmentation error")
if(!round(r[[1]][3,3],3)==round(r[[2]][3,3],3)) stop("Listeria ML compared versus IMP error")
if(!round(r[[3]][3,3],3)==0.034) stop("Listeria DROP dataaugmentation error")


cat("testaugmentation.R, tests succesfully run!")
