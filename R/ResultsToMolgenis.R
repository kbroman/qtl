#####################################################################
#
# ResultsToMolgenis.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: CrossToMolgenis
#
######################################################################

######################################################################
#
# ResultsToMolgenis: Puts results from scanone scanMQM and scanALL into a molgenis database
#
######################################################################

ResultsToMolgenis <- function(intervalQTLmap=NULL,name="MQMresultsTest",Trait_num=0,DBpath=NULL,verbose=TRUE){
	library("RCurl")
	if(!("RCurl" %in% names( getLoadedDLLs()))){
		ourstop("Please install the package RCurl from bioconductor to use the molgenis interface\n")
	}
	if(is.null(DBpath)){
		ourstop("Please provide a valid DBpath\n")
	}else{
		#Set the path to molgenis
		molgenispath <- paste(DBpath,"/api/R/",sep="")
		
		#Source the interface
		if(!exists(".MOLGENIS")){
			source(molgenispath)
		}
	}

	#get data from server
	if(is.null(intervalQTLmap)){
		ourstop("Please supply a QTL interval map\n")
	}
	if(any(class(intervalQTLmap) == "scanone")){
		#cat("INFO: Valid object from scanone, containing 1 phenotype\n")
		num_pheno <- 1
	}
	if(any(class(intervalQTLmap) == "MQMmulti")){
		#cat("INFO: Valid object from MultiQTL scan, containing ",length(intervalQTLmap)," phenotypes\n")
		num_pheno <- length(intervalQTLmap)
	}
	investi <- find.investigation(.verbose=verbose)
	ourcat("INFO: Found",dim(investi)[1],"investigations in the current database\n",a=verbose)
	if("MQMQTL" %in% investi$name){
		ourcat("INFO: Found MQMQTL investigation in the current database\n",a=verbose)
		num <- find.investigation(name="MQMQTL",.verbose=verbose)
	}else{
		ourcat("INFO: Created a new instance of an MQMQTL investigation in the current database\n",a=verbose)
		num <- add.investigation(name="MQMQTL",.verbose=verbose)
	}
	
	#Get all the markers
	markers <- find.marker(.verbose=verbose)
	ourcat("INFO: Found",dim(markers)[1],"markers in the current database\n",a=verbose)
	cnt <- 0
	if(num_pheno == 1){	
		intervalQTLmap = list(intervalQTLmap)
	}
	for(j in 1:num_pheno){
		for(i in 1:dim(intervalQTLmap[[j]])[1]) {
			if(!rownames(intervalQTLmap[[j]])[i] %in% markers$name){
				add.marker(name=rownames(intervalQTLmap[[j]])[i],chr=intervalQTLmap[[j]][i,"chr"],cm=intervalQTLmap[[j]][i,"pos (Cm)"],investigation_id=num$id)
				cnt=cnt+1
			}
		}
	}
	ourcat("INFO: Added",cnt,"markers to the current database\n",a=verbose)	
	#Markers are inside molgenis, now we need to get the QTL's in
	
	aaa <- find.data(name=name,.verbose=verbose)
	colnam <- colnames(intervalQTLmap[[j]])[3]
	trait_name <- find.trait(name=substr(colnam,5,nchar(colnames(intervalQTLmap[[1]])[3])),.verbose=verbose)
	if(!dim(aaa)[1]){
		#No find, so we'll create one
		ourcat("INFO: Not matrix named",name,"found in the current database\n",a=verbose)
		ourcat("INFO: Creating:",name,"in the current database\n",a=verbose)
		aaa <- add.data(name = name,investigation_id=num$id,rowtype="Marker",coltype=trait_name$type,totalrows=1,totalcols=1,valuetype="Decimal",.verbose=verbose)
	}else{
		ourcat("INFO: Matrix named",name,"found in the current database\n",a=verbose)
	}

	for(j in 1:num_pheno){
		colnam <- colnames(intervalQTLmap[[j]])[3]
		colnam <- substr(colnam,5,nchar(colnames(intervalQTLmap[[j]])[3]))
		if(colnam ==""){
			colnam = paste("unknown",j,sep="")
		}
		names <- NULL
		values <- NULL
		rowindex <- NULL
		colindex <- Trait_num
		for(i in 1:dim(intervalQTLmap[[j]])[1]) {
			names <- c(names,rownames(intervalQTLmap[[j]])[i])
			values <- c(values,intervalQTLmap[[j]][i,3])
			rowindex <- c(rowindex,i-1)
		}
		ourcat("INFO: Trying to upload a trait to column:",colindex,"\n",a=verbose)  
		add.decimaldataelement(data_id=aaa$id, col_name=colnam, row_name=names, rowindex=rowindex, colindex=colindex, value=values,.verbose=verbose)
		ourcat("INFO: Uploaded",length(values)," QTL estimates\n",a=verbose)
	}
}

# end of ResultsToMolgenis.R
