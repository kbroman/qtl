# library(foreach)
# library(doParallel)
# 
# cl <- makeCluster(detectCores()-1)
# registerDoParallel(cl)
# 
# 
# detectCores()
# 
# a <- 0
# 
# ls <- foreach(icount(100)) %dopar% {
# 
# 	b <- a + 1
# 
# }
# StopCluster(cl)