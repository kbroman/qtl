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

CrossFromMolgenis <- function(DBmarkerID=0,DBtraitID=0,trait=0,DBpath=NULL,verbose=T){
	library("RCurl")
	if(!("RCurl" %in% names( getLoadedDLLs()))){
		ourstop("Please install the package RCurl from bioconductor to use the molgenis interface\n")
	}

	if(DBmarkerID==0 || DBtraitID==0){
		ourstop("Please provide valid ID's for DBmarker and DBtrait\n")
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
	#marker_data <- find.datamatrix(id=DBmarkerID) #Markerdata
	#trait_data <- find.datamatrix(id=DBtraitID) #Traitdata
	#prepare data for the cross object (we should match which matches which)
	marker_row <- find.data(id=DBmarkerID,.verbose=verbose)["rowtype"]
	marker_col <- find.data(id=DBmarkerID,.verbose=verbose)["coltype"]
	trait_row <- find.data(id=DBtraitID,.verbose=verbose)["rowtype"]
	trait_col <- find.data(id=DBtraitID,.verbose=verbose)["coltype"]
	#Checks
	if(trait_row != "Individual" && trait_col != "Individual"){
		ourstop("No Individuals found in DBtraitID")
	}
	if(marker_row != "Marker" && marker_col != "Marker"){
		ourstop("No markers found in DBmarkerID")
	}	
	
	m_data_url <- paste(DBpath,"/find.datamatrix?id=",DBmarkerID,"&download=all",sep="")
	#cat(m_data_url,"\n")
	marker_data <- read.table(m_data_url,sep="\t")
	if(trait==0){
		t_data_url <- paste(DBpath,"/find.datamatrix?id=",DBtraitID,"&download=all",sep="")
	}else{
		if(trait_col != "Individual"){
			#traits are in the cols
			t_data_url <- paste(DBpath,"/find.datamatrix?id=",DBtraitID,"&download=some&coff=",trait,"&clim=1&roff=0&rlim=1000000",sep="")
		}else{
			#traits are in the rows
			t_data_url <- paste(DBpath,"/find.datamatrix?id=",DBtraitID,"&download=some&coff=0&clim=1000000&roff=",trait,"&rlim=1",sep="")
		}
	}
	#cat(t_data_url,"\n")	
	trait_data <- read.table(t_data_url,sep="\t")
	marker_info <- find.marker(.verbose=verbose)
	marker_info_reduced <- marker_info[,c(1,2,10)]
	
	if(marker_row != "Marker"){
		#cat("INFO: Flipping markerset\n")
		marker_data <- t(marker_data)
	}

	if(trait_col != "Individual"){
		#cat("INFO: Flipping traitset\n")
		trait_data <- t(trait_data)
	}
	ourcat("INFO: Number of individuals in Marker set:",dim(marker_data)[2],"\n",a=verbose)
	ourcat("INFO: Number of individuals in Phenotype set:",dim(trait_data)[2],"\n",a=verbose)
	
	#We assume that if we have IND in markers = IND in trait that individuals match
	if(dim(marker_data)[2] > dim(trait_data)[2]){
		ourcat("INFO: Scaling down the markerset\n",a=verbose)
		matchV <- na.omit(match(colnames(trait_data),colnames(marker_data)))
		marker_data <- marker_data[,matchV]
		matchV <- na.omit(match(colnames(marker_data),colnames(trait_data)))
		trait_data <- trait_data[,matchV]
	}else{
		ourcat("INFO: Scaling down the traitset\n",a=verbose)	
		matchV <- na.omit(match(colnames(marker_data),colnames(trait_data)))
		trait_data <- trait_data[,matchV]
		matchV <- na.omit(match(colnames(trait_data),colnames(marker_data)))
		marker_data <- marker_data[,matchV]
	}
	ourcat("INFO: Number of individuals in Marker set:",dim(marker_data)[2],"\n",a=verbose)
	ourcat("INFO: Number of individuals in Phenotype set:",dim(trait_data)[2],"\n",a=verbose)
	
	#Parse data towards the R/QTL format we need to convert all AA/AB/BB etc to 1,2,3
	for(i in 1:dim(marker_data)[1]) {
		for(j in 1:dim(marker_data)[2]) {
			repl <- FALSE
		    if(is.na(marker_data[i,j])){
				marker_data[i,j] <- NA
				repl <- TRUE
			}
			#RIL
			if(!repl && as.character(marker_data[i,j]) == 'A'){
				marker_data[i,j] <- 1
				repl <- TRUE
			}
			if(!repl && as.character(marker_data[i,j]) == 'B'){
				marker_data[i,j] <- 2
				repl <- TRUE
			}
			#BC
			if(!repl && as.character(marker_data[i,j]) == 'AA'){
				marker_data[i,j] <- 1
				repl <- TRUE
			}
			if(!repl && as.character(marker_data[i,j]) == 'AB'){
				marker_data[i,j] <- 2
				repl <- TRUE
			}
			#F2 intercross
			if(!repl && as.character(marker_data[i,j]) == 'AA'){
				marker_data[i,j] <- 1
				repl <- TRUE
			}
			if(!repl && as.character(marker_data[i,j]) == 'AB'){
				marker_data[i,j] <- 2
				repl <- TRUE
			}
			
			if(!repl && as.character(marker_data[i,j]) == 'BB'){
				marker_data[i,j] <- 3
				repl <- TRUE
			}
			if(!repl){
				marker_data[i,j] <- NA
			}
		}
	}
	chr <- NULL
	loc <- NULL
	names <- NULL
	for(i in 1:dim(marker_data)[1]){
		chr <- c(chr,marker_info[which(rownames(marker_data)[i]== marker_info$name),"chr"])
		loc <- c(loc,marker_info[which(rownames(marker_data)[i]== marker_info$name),"cm"])
		names <- c(names,marker_info[which(rownames(marker_data)[i]== marker_info$name),"name"])
	}
	#FIX for NA in chromosome
	remFromChr <- NULL
	for(i in 1:length(chr)) {
		if(is.na(chr[i])){
			ourcat("INFO: Gonna remove marker #",i,"Which is prob:",names(marker_data[,1])[i],"\n",a=verbose)
			remFromChr <- c(remFromChr,i)

		}
	}
	if(!is.null(remFromChr)){
		chr <- chr[-remFromChr]
		loc <- loc[-remFromChr]
		names <- names[-remFromChr]
		marker_data <- marker_data[-remFromChr,]
	}
	
	#All in expected format, So we can begin filling our cross object
	cross <- NULL
	for(i in sort(unique(chr))){
		# for each chromosome do
		matrix <- NULL
		map <- NULL
		namez <- NULL
		length(cross$geno) <- i
		for(j in which(chr==i)){
				#For all markers on the chromosome do
				#cat("INFO: marker:",j,names[j],"on chr",i,"at",loc[j],"\n")
				matrix <- rbind(matrix,marker_data[j,])
				map <- rbind(map,loc[j])
				namez <- rbind(namez,names[j])
		}
		mapi <- cbind(namez,map)	
		resort <- na.omit(match(sort(as.double(mapi[,2])),as.double(mapi[,2])))
		map <- map[resort]
		namez <- namez[resort]
		matrix <- matrix[resort,]

		#We got everything so lets start adding it to the cross object

		cross$geno[[i]]$data <- t(matrix)
		colnames(cross$geno[[i]]$data)<- namez
		cross$geno[[i]]$map <- as.numeric(t(map))
		names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
		#Type of the chromosome should be retrieved from the database
		class(cross$geno[[i]])[1] <- "A"
	}
	names(cross$geno) <- unique(chr)
	cross$pheno <- as.data.frame(t(trait_data[,]))
	#Make it into a crossobject get the kind of cross from the database (RISELF/RISIBL/F2/BC)
	class(cross)[1] <- "riself"
	class(cross)[2] <- "cross"
	cross
}

# end of CrossFromMolgenis.R
