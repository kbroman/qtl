#####################################################################
#
# bootstrapMQM.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: bootstrapMQM - Main function for bootstrap analysis
# Contains: MQMpermObject - Helperfunction to create permObjects (R/QTL format)
# Contains: snowCoreBOOT - Helperfunction to execute bootstrapping on a R-Term
#
######################################################################

######################################################################
#
# bootstrapMQM: Shuffles phenotype or does parametric bootstrapping of scanMQM
#
######################################################################

#setwd("D:/")
#library(qtl)
#library(MQMpackage)
#cross <- read.cross("csv","","Test.csv")

bootstrapMQM <- function(cross= NULL,cofactors = NULL,pheno.col=1,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.run=10,b_size=10,file="MQM_output.txt",n.clusters=2,parametric=0,plot=TRUE,...)
{
	
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		#Echo back the cross type
		cat("------------------------------------------------------------------\n")
		cat("Starting MQM bootstrap analysis\n")
		cat("Number of bootstrapping runs:",n.run,"\n")
		cat("Batchsize:",b_size," & n.clusters:",n.clusters,"\n")
		cat("------------------------------------------------------------------\n")		
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
		b <- proc.time()		
		if(!parametric){
			cat("INFO: Shuffleling traits between individuals.\n")
		}else{
			cat("INFO: Parametric bootstrapping\nINFO: Calculating new traits for each individual.\n")
		}

		#Set the Phenotype under intrest as the first
		cross$pheno[[1]] <- cross$pheno[[pheno.col]]

		#Some tests from scanMQM repeated here so they are not hidden when using snow
		if((step.min+step.size) > step.max){
				ourstop("Current Step setting would crash the algorithm")
		}
		if(step.min>0){
				ourstop("Step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
				ourstop("Step.size needs to be larger than 1")
		}
		if(n.clusters > b_size){
				ourstop("Please have more items in a batch then clusters assigned per batch")
		}

		#Scan the original
		cross <- fill.geno(cross)
		res0 <- lapply(1, snowCoreALL,all_data=cross,cofactors=cofactors,...)
		
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
		if(("snow" %in% installed.packages()[1:dim(installed.packages())[1]])){
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
				res <- parLapply(cl,boots, snowCoreBOOT,all_data=cross,cofactors=cofactors,parametric=parametric,...)
				stopCluster(cl)
				results <- c(results,res)
				if(plot){
					temp <- c(res0,results)
					class(temp) <- c(class(temp),"MQMmulti")
					plot.MQMboot(temp)
				}
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
			for(x in batch:(batches)){
				start <- proc.time()
				ourline()
				cat("INFO: Starting with batch",x,"/",batches,"\n")				
				ourline()			
				if(x==batches && last.batch.num > 0){
					boots <- bootstraps[((b_size*(x-1))+1):((b_size*(x-1))+last.batch.num)]
				}else{
					boots <- bootstraps[((b_size*(x-1))+1):(b_size*(x-1)+b_size)]
				}	
				res <- lapply(boots, snowCoreBOOT,all_data=cross,cofactors=cofactors,parametric=parametric,...)
				if(plot){
					temp <- c(res0,results)
					class(temp) <- c(class(temp),"MQMmulti")
					plot.MQMboot(temp)
				}
				results <- c(results,res)				
				end <- proc.time()
				SUM <- SUM + (end-start)[3]
				AVG <- SUM/x
				LEFT <- AVG*(batches-x)
				cat("INFO: Done with batch",x,"/",batches,"\n")	
				cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
				cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
				cat("INFO: Average time per batch:",round((AVG), digits=3),",per run:",round((AVG %/% b_size), digits=3),"seconds\n")
				cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")				
				ourline()
			}
		}
		res <- c(res0,results)
		#Set the class of the result to MQMmulti (so we can use our plotting routines)
		class(res) <- c(class(res),"MQMmulti")
		e <- proc.time()
		SUM <- (e-b)[3]
		AVG <- SUM/(n.run+1)	
		cat("INFO: Done with MQM bootstrap analysis\n")
		cat("------------------------------------------------------------------\n")
		cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")		
		cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
		cat("------------------------------------------------------------------\n")	
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

snowCoreBOOT <- function(x,all_data,cofactors,parametric,...){
	b <- proc.time()
	if(!parametric){
		neworder <- sample(nind(all_data))			
		all_data$pheno[[1]] <- all_data$pheno[[1]][neworder]
	}else{	
		pheno <- all_data$pheno[[1]]
		variance <- var(pheno,na.rm = TRUE)
		for(j in 1:nind(all_data)) {
			all_data$pheno[[1]][j] <- runif(1)*(variance^0.5)
		}
	}
	result <- scanMQM(all_data,cofactors=cofactors,pheno.col=1,plot=F,verbose=F,...)
	e <- proc.time()
	cat("------------------------------------------------------------------\n")
	cat("INFO: Done with bootstrap\n")
	cat("INFO: Calculation took:",round((e-b)[3], digits=3)," seconds\n")
	cat("------------------------------------------------------------------\n")
	result
}

#result <- bootstrapMQM(cross)

# end of bootstrapMQM.R
