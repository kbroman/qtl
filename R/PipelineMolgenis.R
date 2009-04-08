#####################################################################
#
# PipelineMolgenis.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: PipelineMolgenis
#
######################################################################

######################################################################
#
# PipelineMolgenis:
#
######################################################################

PipelineMolgenis <- function(DBmarkerID,DBtraitID,name="MQMResults",DBpath,each=0,n.clusters=2,...){
	cat("------------------------------------------------------------------\n")
	cat("Starting Molgenis <-> MQM <-> Molgenis automated pipeline\n")
	cat("INFO: Molgenisserver:",DBpath,"\n")
	cat("INFO: Genotype info-tableID:",DBmarkerID," (DBmarkerID)\n")
	cat("INFO: Phenotype values-tableID:",DBtraitID," (DBtraitID)\n")
	cat("INFO: Results will be stored in a table named:",name,"\n")	
	cat("------------------------------------------------------------------\n")
	cat("INFO: Starting data retrieval.\n")	
	cat("INFO: Please be patient while all data is retrieved.\n")	
	cat("------------------------------------------------------------------\n")
	start <- proc.time()
	all_data <- CrossFromMolgenis(DBmarkerID,DBtraitID,trait=0,DBpath,verbose=F)
	all_data <- fill.geno(all_data)
	num_traits <- nphe(all_data)
	summary(all_data)
	end <- proc.time()
	cat("------------------------------------------------------------------\n")
	cat("INFO: Data retrieval finished in",round((end-start)[3], digits = 3),"seconds\n")
	cat("------------------------------------------------------------------\n")
	PRE <- (end-start)[3]
	SUM <- 0
	AVG <- 0
	LEFT <- 0
	#TEST FOR SNOW CAPABILITIES
	if(("snow" %in% installed.packages()[1:dim(installed.packages())[1]])){
		start <- proc.time()
		cat("INFO: Library snow found using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
		outcome <- NULL
		library(snow)
		cl <- makeCluster(n.clusters)
		clusterEvalQ(cl, library(MQMpackage))
		outcome <- parLapply(cl,1:num_traits, snowCore,each=each,all_data=all_data,name=name,DBpath=DBpath)
		stopCluster(cl)
		end <- proc.time()
		SUM <- SUM + (end-start)[3]
		AVG <- SUM/num_traits
		cat("------------------------------------------------------------------\n")		
		cat("INFO: Elapsed time:",(SUM+PRE%/%3600),":",(SUM+PRE%%3600)%/%60,":",round(SUM+PRE%%60, digits=0),"(Hour:Min:Sec) (",round(PRE, digits=3),",",round(SUM, digits=3),")\n")		
		cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
		cat("------------------------------------------------------------------\n")	
	}else{
		for(x in 1:num_traits){
			start <- proc.time()
			cat("\n\n------------------------------------------------------------------\n")
			cat("INFO: Starting analysis of trait (",x,"/",num_traits,")",names(all_data$pheno)[x],"\n")
			cat("------------------------------------------------------------------\n")
			cat("INFO: Scanning for QTL's\n")
			cat("------------------------------------------------------------------\n")
			if(each>1){
				cof <- MQMCofactorsEach(all_data,each)
				result <- scanMQM(all_data,cof,pheno.col=x,plot=T,verbose=F,...)
			}else{
				result <- scanMQM(all_data,pheno.col=x,plot=T,verbose=F,...)
			}
			cat("------------------------------------------------------------------\n")			
			cat("INFO: Finished scanning for QTL's\n")	
			cat("INFO: Uploading calculated QTL's to Molgenis\n")
			cat("------------------------------------------------------------------\n")
				ResultsToMolgenis(result, name,(x-1),DBpath, verbose=F)
			end <- proc.time()		
			SUM <- SUM + (end-start)[3]
			AVG <- SUM/x
			LEFT <- AVG*(num_traits-x)
			cat("------------------------------------------------------------------\n")
			cat("INFO: Finished uploading of QTL's\n")
			cat("------------------------------------------------------------------\n")
			cat("INFO: Calculation of trait",x,"took:",round((end-start)[3], digits=3),"seconds\n")
			cat("INFO: Elapsed time:",(SUM+PRE%/%3600),":",(SUM+PRE%%3600)%/%60,":",round(SUM+PRE%%60, digits=0),"(Hour:Min:Sec) (",round(PRE, digits=3),",",round(SUM, digits=3),")\n")
			cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
			cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
			cat("------------------------------------------------------------------\n")	
		}
	}
}

snowCore <- function(x,each,all_data,name,DBpath,...){
	num_traits <- nphe(all_data)
	b <- NULL
	e <- NULL
	b <- proc.time()
	r_string <- NULL
	r_string <- paste("------------------------------------------------------------------\n")
	r_string <- paste(r_string,"INFO: Starting analysis of trait (",x,"/",num_traits,")\n")
	r_string <- paste(r_string,"------------------------------------------------------------------\n")
	if(each>1){
		cof <- MQMCofactorsEach(all_data,each)
		result <- scanMQM(all_data,cof,pheno.col=x,plot=F,verbose=F,...)
	}else{
		result <- scanMQM(all_data,pheno.col=x,plot=F,verbose=F,...)
	}
	try(ResultsToMolgenis(result, name,(x-1),DBpath, verbose=F),TRUE)
	e <- proc.time()
	r_string <- paste("------------------------------------------------------------------\n")
	r_string <- paste(r_string,"INFO: Done with the analysis of trait (",x,"/",num_traits,")\n")	
	r_string <- paste(r_string,"INFO: Calculation of trait",x,"took:",round((e-b)[3], digits=3)," seconds\n")
	r_string <- paste(r_string,"------------------------------------------------------------------\n")
	r_string
}

#x=1:num_traits
#paste("script_",x,".R",sep="")


