# Test augmentation with MQM
#
# Note: the full version of this test has moved to ./contrib/bin/rtest, 
#       as it takes a long time to run. The full version can be run with:
#
#   cd contrib/bin
#   rm CMakeCache.txt ; cmake -DTEST_R=TRUE
#   make testR

library(qtl)

set.seed(1000)              

version = mqm_version()
cat("R/qtl=",version$RQTL)
cat("R-MQM=",version$RMQM)
cat("MQM=",version$MQM)

testaugmentation <- function(cross, ...){
  crossML <- mqmaugment(cross, ...)

  res1 <- mqmscan(crossML,logtransform=TRUE)
  list(res1)
}

data(listeria)
r <- testaugmentation(listeria)
if(!round(r[[1]][3,3],3)==0.307) stop("Listeria ML dataaugmentation error")

cat("testaugmentation.R, tests succesfully run!")
