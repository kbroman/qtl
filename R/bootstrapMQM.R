#####################################################################
#
# bootstrapMQM.R
#
# copyright (c) 2009, Danny Arends
# last modified Jun, 2009
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
# Contains: bootstrap - Main function for bootstrap analysis
#           MQMpermObject - Helperfunction to create permObjects (R/QTL format)
#           bootstrapmqm, bootstrapcim
#           FDRpermutation
#
######################################################################

FDRpermutation <- function(cross=NULL, Funktie=scanall, thresholds=c(1,2,3,4,5,7,10,15,20), n.perm = 10, verbose=TRUE, ...){
	if(verbose){cat("Calculation of FDR estimate of threshold in multitrait analysis.\n")}
	results <- NULL
	above_in_real_res <- NULL
	res <- Funktie(cross,...)
	for(threshold in thresholds){
		above_in_real <- 0
		for(x in 1:nphe(cross)){
			above_in_real = above_in_real + sum(res[[x]][,3] > threshold)
		}
		above_in_real_res <- c(above_in_real_res,above_in_real)
	}
	perm <- cross
	if(verbose){cat("QTL's above threshold:",above_in_real,"\n")}
	above_in_perm_res <- rep(0,length(thresholds))
	for(x in 1:n.perm){
		if(verbose){cat("Starting permutation",x,"\n")}
		perm_res <- NULL
		
		neworder <- sample(nind(cross))
		for(chr in 1:nchr(cross)){
			perm$geno[[1]]$data <- perm$geno[[1]]$data[neworder,]
		}
		res <- Funktie(perm,...)
		for(threshold in thresholds){
			above_in_perm <- 0
			for(y in 1:nphe(cross)){
				above_in_perm = above_in_perm + sum(res[[y]][,3] > threshold)
			}
			perm_res <- c(perm_res,above_in_perm)
			#if(verbose){cat("Permutation",x,"QTL's above threshold:",above_in_perm,"\n")}
		}
		above_in_perm_res <- above_in_perm_res+perm_res
	}
	above_in_perm_res <- above_in_perm_res/n.perm
	results <- cbind(above_in_real_res,above_in_perm_res,above_in_perm_res/above_in_real_res)
	rownames(results) <- thresholds
	results
}



bootstrapmqm <- function(...){
	bootstrap(...,Funktie=mqm)
}

bootstrapcim <- function(...){
	bootstrap(...,Funktie=cim)
}

######################################################################
#
# bootstrap: Shuffles phenotype or does parametric bootstrapping of scanMQM
#
######################################################################

bootstrap <- function(cross= NULL,Funktie=scanone,pheno.col=1,multiC=TRUE,n.run=10,b_size=10,file="MQM_output.txt",n.clusters=2,method=0,plot=FALSE,verbose=FALSE,...)
{
	
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		#Echo back the cross type
		if(verbose) {
                  cat("------------------------------------------------------------------\n")
                  cat("Starting bootstrap analysis\n")
                  cat("Number of bootstrapping runs:",n.run,"\n")
                  cat("Batchsize:",b_size," & n.clusters:",n.clusters,"\n")
                  cat("------------------------------------------------------------------\n")		
                  cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
                }
		b <- proc.time()		
		if(!method){
			if(verbose) cat("INFO: Shuffleling traits between individuals.\n")
		}else{
			if(verbose) cat("INFO: Parametric bootstrapping\nINFO: Calculating new traits for each individual.\n")
		}

		#Set the Phenotype under intrest as the first
		cross$pheno[[1]] <- cross$pheno[[pheno.col]]

		if(n.clusters > b_size){
				ourstop("Please have more items in a batch then clusters assigned per batch")
		}

		#Scan the original
		cross <- fill.geno(cross)
		res0 <- lapply(1, FUN=snowCoreALL,all_data=cross,Funktie=Funktie,verbose=verbose,...)
		
		#Setup bootstraps by generating a list of random numbers to set as seed for each bootstrap
		bootstraps <- runif(n.run)
		batches <- length(bootstraps) %/% b_size
		last.batch.num <- length(bootstraps) %% b_size
		results <- NULL
		if(last.batch.num > 0){
			batches = batches+1
		}
		SUM <- 0
		AVG <- 0
		LEFT <- 0
		#TEST FOR SNOW CAPABILITIES
#		if(("snow" %in% installed.packages()[1:dim(installed.packages())[1]]) && multiC){
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
				res <- parLapply(cl,boots, fun=snowCoreBOOT,all_data=cross,Funktie=Funktie,method=method,verbose=verbose,...)
				stopCluster(cl)
				results <- c(results,res)
				if(plot){
					temp <- c(res0,results)
					class(temp) <- c(class(temp),"MQMmulti")
					plotMQMboot(temp)
				}
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
				res <- lapply(boots, FUN=snowCoreBOOT,all_data=cross,Funktie=Funktie,method=method,verbose=verbose,...)
				results <- c(results,res)	
				if(plot){
					temp <- c(res0,results)
					class(temp) <- c(class(temp),"MQMmulti")
					plotMQMboot(temp)
				}
				end <- proc.time()
				SUM <- SUM + (end-start)[3]
				AVG <- SUM/x
				LEFT <- AVG*(batches-x)
                                if(verbose) {
                                  cat("INFO: Done with batch",x,"/",batches,"\n")	
                                  cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
                                  cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
                                  cat("INFO: Average time per batch:",round((AVG), digits=3),",per run:",round((AVG %/% b_size), digits=3),"seconds\n")
                                  cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")				
                                  ourline()
                                }
			}
		}
		res <- c(res0,results)
		#Set the class of the result to MQMmulti (so we can use our plotting routines)
		class(res) <- c(class(res),"MQMmulti")
		e <- proc.time()
		SUM <- (e-b)[3]
		AVG <- SUM/(n.run+1)	
                if(verbose) {
                  cat("INFO: Done with MQM bootstrap analysis\n")
                  cat("------------------------------------------------------------------\n")
                  cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")		
                  cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
                  cat("------------------------------------------------------------------\n")
                }
		res
	}else{
		ourstop("Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}

MQMpermObject <- function(MQMbootresult = NULL){
	if(class(MQMbootresult)[2] == "MQMmulti"){
		result <- NULL
		names <- NULL
		for(i in 2:length(MQMbootresult)) {
			result <- rbind(result,max(MQMbootresult[[i]][,3]))
			names <- c(names,i-1)
		}
		result <- as.matrix(result)
		rownames(result) <- names
		result <- cbind(result,result,result)
		class(result) <- c("scanoneperm",class(result))
		result
	}else{
		ourstop("PLease supply a valid resultobject (MQMmulti).")
	}
}

#result <- bootstrap(cross)
# tiff(object, file="namemeplease.tiff" res=300, unit="in", width=6, height=6)

# end of bootstrapMQM.R
