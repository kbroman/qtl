#####################################################################
#
# MQMLicence.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: MQMLicence
#
######################################################################

MQMLicence <- function(){
    text <- paste("This software is distributed under the terms of \n")
	text <- paste("Developed by: R.C. Janssen\nImplementation into R: D. Arends\n")
	MQMwindow(text)
}

MQMwindow <- function(text){
    outFile <- tempfile()
    outConn <- file(outFile, open = "w")
    writeLines(paste(text, sep = ""), outConn)
    close(outConn)
    file.show(outFile, delete.file = TRUE)
}

ourcat <- function(...,a=TRUE){
	if(a){
		cat(...)
	}else{
		return
	}
}

ourstop <- function(...){
	library(tools)
	parents <- sys.calls()
	ourline()
	cat("ERROR: Detected in function: '",as.character(as.call(parents[[1]])),"'\n",sep="")	
	ourline()
	db <- Rd_db("MQMpackage")
	db <- lapply(db, function(txt) Rd_parse(text = txt))
	data <- lapply(db, "[[", "data")
	p_Rd <- paste(as.character(parents[[1]]),".Rd",sep="")
	#if(!is.null(data[[p_Rd]]$tags) && "examples" %in% data[[p_Rd]]$tags){
	#	ourline()
	#	cat("ERROR: Some examples on how-to use this function from the helpfiles\n")
	#	cat("ERROR: To show help type: \"?",as.character(parents[[1]]),"\"\n",sep="")
	#	ourline()
	#	aa <- which(data[[p_Rd]]$tags=="examples")
	#	cat(data[[p_Rd]]$vals[[aa]])
	#	ourline()
	#}
	stop(...)
}

ourline <- function(){
	cat("------------------------------------------------------------------\n")		
}




## Extract the metadata.



