#####################################################################
#
# CrossFromMolgenis.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: CrossFromMolgenis
#
######################################################################

######################################################################
#
# CrossFromMolgenis:
#
######################################################################

ResultsFromMolgenis <- function(DBqtlName=NULL,DBpath=NULL,verbose=TRUE){
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
	if(is.null(DBqtlName)){
		ourstop("Please supply either a qtl DB name or a qtl DB identifier\n")
	}else{
		dataM <- find.data(name=DBqtlName,.verbose=verbose)
		returnObj <- NULL
		tempObj <- NULL
		data <- find.datamatrix(id=dataM$id,.verbose=verbose)
		marker_info <- find.marker(.verbose=verbose)
		matchV <- match(rownames(data),marker_info$name)
		marker_info <- marker_info[matchV,]
		
		if(dim(data)[2] == 1){
			#We need to give back a single QTL scanone object
			returnObj <- cbind(as.numeric(marker_info$chr),as.numeric(marker_info$cm),as.numeric(data[,1]))
			colnames(returnObj) <- c("chr","cm","QTL")
			returnObj <- returnObj[order(returnObj[,"chr"],returnObj[,"cm"]),]
			returnObj <- as.data.frame(returnObj)
			class(returnObj) <- c("scanone",class(returnObj))
		}else{	
			for(i in 1:dim(data)[2]) {
				tempObj <- cbind(as.numeric(marker_info$chr),as.numeric(marker_info$cm),as.numeric(data[,i]))
				colnames(tempObj) <- c("chr","cm","QTL")
				tempObj <- tempObj[order(tempObj[,"chr"],tempObj[,"cm"]),]
				class(tempObj) <- c(class(tempObj),"scanone")
				returnObj[[i]] <- tempObj
			}
			class(returnObj) <- c(class(returnObj),"MQMmulti")
		}
		returnObj
	}
}

# end of ResultsFromMolgenis.R
