#####################################################################
#
# scanMQMall.R
# Alias: scanall, cimall, mqmall
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Feb, 2009
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/qtl package
# Contains: scanall, cimall, mqmall
#
######################################################################

cimall <- function(...){
	scanall(...,Funktie=cim)
}

mqmall <- function(...){
	scanall(...,Funktie=mqm)
}

scanall <- function(cross= NULL,Funktie=scanone,multiC=TRUE,n.clusters=2,b_size=10,FF=0,...,plot=FALSE,verbose=FALSE){

	
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		start <- proc.time()
		n.pheno <- nphe(cross)
		if(verbose) {
                  ourline()
                  cat("Starting R/QTL multitrait analysis\n")
                  cat("Number of phenotypes:",n.pheno,"\n")
                  cat("Batchsize:",b_size," & n.clusters:",n.clusters,"\n")
                  ourline()	
                }
		
		result <- NULL	 	#BATCH result variable
		res <- NULL			#GLOBAL result variable
		
		all_data <- fill.geno(cross)
		
		bootstraps <- 1:n.pheno
		batches <- length(bootstraps) %/% b_size
		last.batch.num <- length(bootstraps) %% b_size
		if(last.batch.num > 0){
			batches = batches+1
		}
		
		#INIT TIME VARS
		SUM <- 0
		AVG <- 0
		LEFT <- 0

		#TEST FOR SNOW CAPABILITIES
		if(multiC && n.clusters>1 && suppressWarnings(require(snow,quietly=TRUE))) {
                  if(verbose) cat("INFO: Library snow found using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
			for(x in 1:(batches)){
				start <- proc.time()
                                if(verbose) {
                                  ourline()
                                  cat("INFO: Starting with batch",x,"/",batches,"\n")				
                                  ourline()
                                }
				if(x==batches && last.batch.num > 0){
					boots <- bootstraps[((b_size*(x-1))+1):((b_size*(x-1))+last.batch.num)]
				}else{
					boots <- bootstraps[((b_size*(x-1))+1):(b_size*(x-1)+b_size)]
				}	
				cl <- makeCluster(n.clusters)
				clusterEvalQ(cl, require(qtl, quietly=TRUE))
				result <- parLapply(cl,boots, fun=snowCoreALL,all_data=all_data,Funktie=Funktie,verbose=verbose,...)
				stopCluster(cl)
				if(plot){
					temp <- result
					class(temp) <- c(class(temp),"MQMmulti")
					plotMQMnice(temp)
				}
				res <- c(res,result)
				end <- proc.time()
				SUM <- SUM + (end-start)[3]
				AVG <- SUM/x
				LEFT <- AVG*(batches-x)
                                if(verbose) {
                                  cat("INFO: Done with batch",x,"/",batches,"\n")	
                                  cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
                                  cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
                                  cat("INFO: Average time per batch:",round((AVG), digits=3)," per trait:",round((AVG %/% b_size), digits=3),"seconds\n")
                                  cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
                                  ourline()
                                }
			}
		}else{
			if(verbose) cat("INFO: Library snow not found, so going into singlemode.\n")
			for(x in 1:(batches)){
				start <- proc.time()
				if(verbose) {
                                  ourline()
                                  cat("INFO: Starting with batch",x,"/",batches,"\n")				
                                  ourline()
                                }
				if(x==batches && last.batch.num > 0){
					boots <- bootstraps[((b_size*(x-1))+1):((b_size*(x-1))+last.batch.num)]
				}else{
					boots <- bootstraps[((b_size*(x-1))+1):(b_size*(x-1)+b_size)]
				}	
				result <- lapply(boots, FUN=snowCoreALL,all_data=all_data,Funktie=Funktie,verbose=verbose,...)
				if(plot){
					temp <- result
					class(temp) <- c(class(temp),"MQMmulti")
					plotMQMnice(temp)
				}
				res <- c(res,result)				
				end <- proc.time()
				SUM <- SUM + (end-start)[3]
				AVG <- SUM/x
				LEFT <- AVG*(batches-x)
                                if(verbose) {
                                  cat("INFO: Done with batch",x,"/",batches,"\n")	
                                  cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
                                  cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
                                  cat("INFO: Average time per batch:",round((AVG), digits=3)," per trait:",round((AVG %/% b_size), digits=3),"seconds\n")
                                  cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
                                  ourline()
                                }
			}
		}
		if(FF){
			if(verbose) cat(rownames(res[[1]]),"\n",res[[1]][,1],"\n",res[[1]][,2],"\n",file="out.frank")
			for(i in 1:length(res)){
				if(verbose) cat("INFO: Saving trait",i,"in frankformat\n")
				qtl <- res[[i]]
				if(verbose) cat(colnames(qtl)[3],qtl[,3],"\n",file="out.frank",append = TRUE)
			}
		}
		#Return the results
		class(res) <- c(class(res),"MQMmulti")
		#All done now plot the results
		end <- proc.time()
		SUM <- SUM + (end-start)[3]
		AVG <- SUM/n.pheno	
                if(verbose) {
                  cat("------------------------------------------------------------------\n")
                  cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")		
                  cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
                  cat("------------------------------------------------------------------\n")
                }
		res
	}else{
		stop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}


# end of scanMQMall.R
