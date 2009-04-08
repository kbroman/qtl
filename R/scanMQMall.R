#####################################################################
#
# scanMQMall.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: scanMQM
#
######################################################################

######################################################################
#
# scanMQMall: Contains scanMQMall routine and the plot.MQMall routine
#
######################################################################

scanMQMall <- function(cross= NULL,cofactors = NULL,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.clusters=2,b_size=10,FF=0,plot=TRUE,verbose=TRUE,...){

	
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		start <- proc.time()
		n.pheno <- nphe(cross)
		cat("------------------------------------------------------------------\n")
		cat("Starting MQM multitrait analysis\n")
		cat("Number of phenotypes:",n.pheno,"\n")
		cat("Batchsize:",b_size," & n.clusters:",n.clusters,"\n")
		cat("------------------------------------------------------------------\n")		
		
		result <- NULL	 	#BATCH result variable
		res <- NULL			#GLOBAL result variable
		
		all_data <- fill.geno(cross)
		#Some tests from scanMQM repeated here so they are not hidden when using snow
		if((step.min+step.size) > step.max){
				ourstop("Surrent Step settings (step.min/step.max) would crash the algorithm")
		}
		if(step.min>0){
				ourstop("Step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
				ourstop("Step.size needs to be larger than 1")
		}
		
		
		bootstraps <- 1:n.pheno
		batches <- length(bootstraps) %/% b_size
		last.batch.num <- length(bootstraps) %% b_size
		if(last.batch.num > 0){
			batches = batches+1
		}
		SUM <- 0
		AVG <- 0
		LEFT <- 0
		#TEST FOR SNOW CAPABILITIES
		if("snow" %in% installed.packages()[1:dim(installed.packages())[1]]){
			cat("INFO: Library snow found using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
			library(snow)
			for(x in 1:(batches)){
				start <- proc.time()
				ourline()
				cat("INFO: Starting with batch",x,"/",batches,"\n")				
				ourline()
				if(x==batches && last.batch.num > 0){
					boots <- bootstraps[((b_size*(x-1))+1):((b_size*(x-1))+last.batch.num)]
				}else{
					boots <- bootstraps[((b_size*(x-1))+1):(b_size*(x-1)+b_size)]
				}	
				cl <- makeCluster(n.clusters)
				clusterEvalQ(cl, library(MQMpackage))
				result <- parLapply(cl,boots, snowCoreALL,all_data=all_data,cofactors=cofactors,...)
				stopCluster(cl)
				if(plot){
					temp <- result
					class(temp) <- c(class(temp),"MQMmulti")
					plot.MQMnice(temp)
				}
				res <- c(res,result)
				end <- proc.time()
				SUM <- SUM + (end-start)[3]
				AVG <- SUM/x
				LEFT <- AVG*(batches-x)
				cat("INFO: Done with batch",x,"/",batches,"\n")	
				cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
				cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
				cat("INFO: Average time per batch:",round((AVG), digits=3)," per trait:",round((AVG %/% b_size), digits=3),"seconds\n")
				cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
				ourline()
			}
		}else{
			cat("INFO: Library snow not found, so going into singlemode.\n")
			for(x in 1:(batches)){
				start <- proc.time()
				ourline()
				cat("INFO: Starting with batch",x,"/",batches,"\n")				
				ourline()
				if(x==batches && last.batch.num > 0){
					boots <- bootstraps[((b_size*(x-1))+1):((b_size*(x-1))+last.batch.num)]
				}else{
					boots <- bootstraps[((b_size*(x-1))+1):(b_size*(x-1)+b_size)]
				}	
				result <- lapply(cl,boots, snowCoreALL,all_data=all_data,cofactors=cofactors,...)
				if(plot){
					temp <- result
					class(temp) <- c(class(temp),"MQMmulti")
					plot.MQMnice(temp)
				}
				res <- c(res,result)				
				end <- proc.time()
				SUM <- SUM + (end-start)[3]
				AVG <- SUM/x
				LEFT <- AVG*(batches-x)
				cat("INFO: Done with batch",x,"/",batches,"\n")	
				cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
				cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
				cat("INFO: Average time per batch:",round((AVG), digits=3)," per trait:",round((AVG %/% b_size), digits=3),"seconds\n")
				cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
				ourline()
			}
		}
		if(FF){
			cat(rownames(res[[1]]),"\n",res[[1]][,1],"\n",res[[1]][,2],"\n",file="out.frank")
			for(i in 1:length(res)){
				cat("INFO: Saving trait",i,"in frankformat\n")
				qtl <- res[[i]]
				cat(colnames(qtl)[3],qtl[,3],"\n",file="out.frank",append = T)
			}
		}
		#Return the results
		class(res) <- c(class(res),"MQMmulti")
		#All done now plot the results
		end <- proc.time()
		SUM <- (end-start)[3]
		AVG <- SUM/n.pheno	
		cat("------------------------------------------------------------------\n")
		cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")		
		cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
		cat("------------------------------------------------------------------\n")	
		res
	}else{
		stop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}

snowCoreALL <- function(x,all_data,cofactors,...){
	b <- proc.time()
	num_traits <- nphe(all_data)
	cat("------------------------------------------------------------------\n")
	cat("INFO: Starting analysis of trait (",x,"/",num_traits,")\n")
	cat("------------------------------------------------------------------\n")
	result <- scanMQM(all_data,cofactors=cofactors,pheno.col=x,plot=F,verbose=F,...)
	e <- proc.time()
	cat("------------------------------------------------------------------\n")
	cat("INFO: Done with the analysis of trait (",x,"/",num_traits,")\n")	
	cat("INFO: Calculation of trait",x,"took:",round((e-b)[3], digits=3)," seconds\n")
	cat("------------------------------------------------------------------\n")
	result
}


# end of scanMQMall.R
