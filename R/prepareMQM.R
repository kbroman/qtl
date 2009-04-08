#####################################################################
#
# prepareMQM.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: prepareMQM, printMQMin
#
######################################################################

MQMLoadanalysis <- function(file="MQM_output.txt"){
	data <- read.table(file)
	class(data) <- c("scanone",class(data))
	plot(data)
	data
}

prepareMQM <- function(cross, name,cofactors=NULL,dominance='n',RemLorML=0){
#print files needed by the MQM algorithm (Genotype & Phenotype, also calls 
#the construction function of the run-input file)

 f1mar <- NULL
 f2mar <- t(pull.geno(cross))
 n.ind <- nind(cross)
 for(i in 1:dim(f2mar)[1]) {
   f1mar <- rbind(f1mar,12)
   for(j in 1:dim(f2mar)[2]) {
    if(is.na(f2mar[i,j])){
	     f2mar[i,j] <- '-'
	}
	if(as.character(f2mar[i,j]) == '1'){
        f2mar[i,j] <- 'A'
    }
    if(as.character(f2mar[i,j]) == '2'){
        f2mar[i,j] <- 'H'
    }
    if(as.character(f2mar[i,j]) == '3'){
        f2mar[i,j] <- 'B'
    }
    if(as.character(f2mar[i,j]) == '4'){
        f2mar[i,j] <- 'C'
    }
    if(as.character(f2mar[i,j]) == '5'){
        f2mar[i,j] <- 'D'
    }	
	}
 }
 rownames(f1mar) = rownames(f2mar)
 pheno <- cross$pheno
 #check for missing phenotypes
 dropped <- NULL
 for(i in 1:dim(pheno)[1]) {
		if(is.na(pheno[i,1])){
			cat("INFO: Dropped individual ",i," with missing phenotype.\n")
			dropped <- c(dropped,i) 
			n.ind = n.ind-1
		}
	}
 #throw em out
    if(!is.null(dropped)){
			f2mar<- f2mar[,-dropped]  
			pheno <- pheno[-dropped,]
		}
 filename <- paste(name,"_F2.MAR.TXT", sep="")
 write.table(f2mar, file = filename, col.names=FALSE, quote=FALSE)

 filename <- paste(name,"_F2.QUA.TXT",sep="")
 write.table(t(pheno), file = filename, row.names=FALSE,col.names=FALSE, quote=FALSE)
 MQM_in <- printMQMin(cross,name,cofactors,dominance,RemLorML) 
}

printMQMin <- function(cross,name,cofactors=NULL,dominance='n',RemLorML=0){
#print mqm_in.txt file needed by the MQM algorithm
	pheno <- cross$pheno	
	n.individuals <- nind(cross)
	for(i in 1:dim(pheno)[1]) {
		if(is.na(pheno[i,1])){
			n.individuals = n.individuals-1
		}
	}	
	n.fam <- 1
	n.mark <- sum(nmar(cross))
	info <- NULL
	info <- rbind(info,c("Nind=",n.individuals))
	info <- rbind(info,c("Nfam=",n.fam))
	info <- rbind(info,c("Nmark=",n.mark))
	
	write.table(info, file = "mqm_in.txt", append = FALSE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	# Printing the markers per chromosome
	n.chr <- nchr(cross)
	cur.marker <- 0
    map <- pull.map(cross)
	result <- NULL
	for(i in 1:n.chr) {
		n.markers <- length(map[[i]])
		for(j in 1:n.markers){
		    if(cur.marker %in% cofactors){
				result <- rbind(result,c(i,row.names(as.matrix(map[[i]]))[j],as.matrix(map[[i]])[j],"*"))
			}else{
			    result <- rbind(result,c(i,row.names(as.matrix(map[[i]]))[j],as.matrix(map[[i]])[j]," "))
			}
			cur.marker <- cur.marker+1
		}
	}
	write.table(result, file = "mqm_in.txt", append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	settings <- NULL;
	settings <- rbind(settings,c("Cross=","f2"))
	settings <- rbind(settings,c("FileF1=","-"))
	settings <- rbind(settings,c("FileF2=",paste(name,"_F2.MAR.TXT",sep="")))
	settings <- rbind(settings,c("FileY=",paste(name,"_F2.QUA.TXT",sep="")))
	settings <- rbind(settings,c("Dominance=",dominance))
	settings <- rbind(settings,c("RemLorML=",RemLorML))
	settings <- rbind(settings,c("defset=","y"))
	settings <- rbind(settings,c("real_simu=","0"))
	settings <- rbind(settings,c("perm_simu=","1 0"))
	write.table(as.matrix(settings), file = "mqm_in.txt", append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)	
}

RunTest <- function(testset = "T1", exe="V_1.exe"){
#Runs an executable versus a, and compares the output to the output from runPrelimTest 
 
 shell("del mqm_out.txt", intern=TRUE) 
 #read old output
 v0.output <- read.table(paste("V_0/",testset,"_OUT",sep=""),comment.char = ":")
 #setup the In-File for testing
 outcome <- shell(paste("copy V_0\\",testset,"_mqm_in.txt mqm_in.txt",sep=""), intern=TRUE)
 outcome <- shell(exe, intern=TRUE)
 vNew_out <- readMQMout()
 error <- 0
 for(i in 1:dim(v0.output)[1]) {
   for(j in 1:dim(v0.output)[2]) {
    if(abs(v0.output[i,j] - vNew_out[i,j]) > 0.01){
	#print(paste("WARNING No match between old and new version at: ",i,j,"",sep=""))
	error <- 1;
	}
  }
 }
 if(error == 0){
   print(paste("No error detected between V_0.exe and ",exe,sep=""))
 }else{
   print(paste("Errors detected between V_0.exe and ",exe,sep=""))
 }
}

RunPrelimTest <- function(testset = "T1",qtl = NULL,cofactors=NULL,dominance='n',RemLorML=0){
#Runs MQM (original version on a testset)
#Places files in an "V_0" directory (which should be setup in advance by the user

   library(qtl)
   setwd("D:/MQM/compiled")
   library(qtl)
   shell("del mqm_out.txt", intern=TRUE)    # make sure the output from the previous runs are gone
   data(map10)
   cross <- sim.cross(map10,qtl,n=100)
   prepareMQM(cross,testset,cofactors,dominance,RemLorML)
   outcome <- shell("V_0.exe", intern=TRUE)
   write.table(readMQMout(), file = paste("V_0/",testset,"_OUT",sep=""), row.names=FALSE,col.names=FALSE, quote=FALSE)
   shell(paste("copy mqm_in.txt V_0\\",testset,"_mqm_in.txt",sep=""), intern=TRUE)   #Stores the runFILE & output to a V_0 directory
   shell("del mqm_out.txt", intern=TRUE)   											#make sure
   shell("del mqm_in.txt", intern=TRUE)    											# also replace the input file (for the new run
}


GenerateTestSets <- function(){
#Run once version to setup some testfiles
#Executes "V_0.exe" to get output for each testfile (so we can compare newer versions)

   RunPrelimTest(testset="T1")
   RunPrelimTest(testset="T1c",cofactors=c(10,50,100))   
   RunPrelimTest(testset="T2",qtl=c(1,30,1,0))
   RunPrelimTest(testset="T2c",qtl=c(1,30,1,0),cofactors=c(10,50,100))
   RunPrelimTest(testset="T3",qtl=c(19,30,0,1))
   RunPrelimTest(testset="T4",qtl=c(19,500,0,1))
   RunPrelimTest(testset="T5",qtl=rbind(c(19,5,0,1),c(19,45,0,1)))
   RunPrelimTest(testset="T6",qtl=rbind(c(19,5,0,1),c(19,45,0,-1)))

}

loadOAT <- function(){
	setwd("D:/test/Oat")
	n.ind=50
	trait <- read.table("trait.qua")
	markers <- read.table("rflp.mar",colClasses=c("character","character"),row.names=1)
	info <- read.table("mqm_R.txt",colClasses=c("integer","character","double"),fill=T)
	n.mark=dim(markers)[1]
	all.mark = dim(info)[1]
	mar <- NULL
	for(m in 1:n.mark){
		colmn <- NULL
		for(i in 1:n.ind){
				num <- as.integer(substr(markers[m,1],i,i))
				f <- 0
				if(num == 3){
					num <- 3
					f <- 1
				}
				if(!f && num == 1){
					num <- 1
					f <- 1
				}
				if(!f && num == 2){
					num <- 2
					f <- 1
				}
				if(!f && num == 0){
					num <- NA
				}

				colmn <- c(colmn,num)
		}
		mar <- cbind(mar,colmn)
	}
	colnames(mar) <- rownames(markers)
	mar <- t(mar)

	rownames(info) <- info[,2]
	selected_markers <- substr(tolower(rownames(markers)),2,10)
	all_markers <- substr(tolower(info[,2]),2,10)
	conv <- which(selected_markers %in% all_markers)
	new_markers <- NULL
	new_info <- NULL
	for (i in conv){
		num <- NULL
		if(!is.na(which(selected_markers[i]==all_markers)==0&&1)){
			num <- which(selected_markers[i]==all_markers)
			new_markers <- rbind(new_markers,mar[i,])
			new_info <- rbind(new_info,info[num,])
		}
	}
	rownames(new_markers) <- rownames(new_info)
	info <- new_info
	mar <- new_markers
	chr <- info[,1]
	cross <- NULL
	cnt <- 1
	for(i in sort(unique(chr))){
		matrix <- NULL
		map <- NULL
		names <- NULL
		cat("Chromosome:",i,"\n")
		for(j in which(chr==i)){
			cat("Marker:",j,"->",info[j,2],"\n")	
			#For all markers on the chromosome do
			matrix <- rbind(matrix,mar[j,])
			map <- rbind(map,info[j,3])
			names <- c(names,info[j,2])
		}
		map <- as.numeric(map)
		rownames(matrix) <- names
		names(map) <- names
		order <- NULL
		num <- length(map)
		while(num>0){
			order <- c(which(as.double(sort(map)[num])==as.double(map)),order)
			if(length(which(as.double(sort(map)[num])==as.double(map))) != 0){
				num = num-length(which(as.double(sort(map)[num])==as.double(map)))
			}else{
				num = num-1
			}
		}
		matrix <- matrix[order,]
		map <- map[order]
		cat(map,"\n")
		#We got everything so lets start adding it to the cross object
		length(cross$geno) <- length(cross$geno)+1
		if(!is.null(dim(matrix)[1])){
			cross$geno[[cnt]]$data <- as.matrix(t(matrix))
		}else{
			cat("INFO SINGLE MARKER ON A CHROMOSOME\n")
			cat(cnt,"->",names,"\n")
			cross$geno[[cnt]]$data <- as.matrix(matrix)
			colnames(cross$geno[[cnt]]$data) <- names
		}
		cat("DIMS:",nrow(cross$geno[[cnt]]$data),ncol(cross$geno[[cnt]]$data),"\n")
		cross$geno[[cnt]]$map <- map
		#Type of the chromosome should be retrieved from the database
		class(cross$geno[[cnt]])[1] <- "A"
		cnt = cnt +1
	}
	names(cross$geno) <- unique(chr)
	cross$pheno <- as.data.frame(t(trait[1:n.ind]))
	#Make it into a crossobject get the kind of cross from the database (RISELF/RISIBL/F2/BC)
	class(cross)[1] <- "f2"
	class(cross)[2] <- "cross"
	cross
}

readMQMout <- function(cross = NULL, file = "mqm_out.txt", plot = TRUE,chr = 1){
#reads the output from the MQM algorithm
   data <-read.table(file, quote=":")
   data <-data[-dim(data)[1],]
   if(plot){
       plot(rownames(data),data[,1],xlab="Markers",ylab="QTL",main="MQM", type='n')
       lines(rownames(data),data[,3],col="red",cex=0.5,pch=20)
       lines(rownames(data),data[,4],col="blue",cex=0.5,pch=20)  
   }
   rownames(data) <- paste(data[,1],data[,2])
   data <- data[,-4]
   colnames(data) <- c("chr","pos","lod")
   data
   #should be pushed to the cross object
}


MQMfind.marker <- function(cross=NULL,scanMQM=NULL,perm=NULL,alpha=0.05){
	
	chr <- summary(scanMQM,alpha=alpha,perms=perm,pvalues=F)$'chr'
	pos <- summary(scanMQM,alpha=alpha,perms=perm,pvalues=F)$'pos (Cm)'
	cat("INFO: Found",length(chr),"markers with alpha <",alpha,".\n")
	ret <- NULL
	for(i in 1:length(chr)){
		ret <- rbind(ret,cbind(find.marker(cross,chr=chr[i],pos=pos[i]),as.integer(chr[i]),as.double(pos[i])))
	}
	colnames(ret) <- c("marker","chr","pos (Cm)")
	ret
}

loadMOUSE <- function(pheno=c(1:10),gr1=NULL,gr2=NULL){
	
	setwd("D:/test/Illumina")
	library(MQMpackage)
	if(gr1==gr2){
		ourstop("Cannot compare groups with themselves.\n")
	}
	selected <- NULL
	#We select from a ordered file, the first cols (1:19) are stemcells, the progenitor, RED and White
	if(gr1 == "S" || gr2 == "S"){
		selected <- c(selected,1:19)
	}
	if(gr1 == "P" || gr2 == "P"){
		selected <- c(selected,20:33)
	}
	if(gr1 == "R" || gr2 == "R"){
		selected <- c(selected,34:48)
	}
	if(gr1 == "W" || gr2 == "W"){
		selected <- c(selected,49:62)
	}
	notselected <- NULL
	for(i in 1:62){
		if(!(i %in% selected)){
			notselected <- c(notselected,i)
		}
	}
	genotypes <- read.csv("BXD.geno",sep="\t")
	#phenotypes <- read.csv("Experiment.pheno",sep=",")
	#save(phenotypes,file="pheno.Rdata")
	load("pheno.Rdata")
	phenotypes <- phenotypes[,-notselected]
	conversion <- read.csv("Conversion2.csv",sep=";")
	n.ind <- dim(phenotypes)[2]
	trait <- t(phenotypes[1,])
	traits <- t(phenotypes[pheno,])
	genomatrix <- NULL
	for (i in 1:n.ind){
		id <- which(rownames(trait)[i]==conversion[,2])
		bxd <- as.character(conversion[id,1])
		bxd <- gsub(" ","",bxd)
		#cat(id," ",bxd,"\n")
		genomatrix <- rbind(genomatrix,as.character(genotypes[,bxd]))
		rownames(genomatrix)[i] <- as.character(bxd)
	}
	for(i in 1:dim(genomatrix)[1]){
		for(j in 1:dim(genomatrix)[2]){
			f <- 0
			if(!f && genomatrix[i,j] == "B"){
				genomatrix[i,j] <- as.integer(1)
				f <- 1
			}
			if(!f && genomatrix[i,j] == "H"){
				genomatrix[i,j] <- NA
				f <- 1
			}
			if(!f && genomatrix[i,j] == "D"){
				genomatrix[i,j] <- as.integer(2)
				f <- 1
			}
			if(!f && genomatrix[i,j] == "U"){
				genomatrix[i,j] <- NA
				f <- 1
			}		
		}
	}
	colnames(genomatrix) <- genotypes[,2]
	chr <- as.character(genotypes[,1])
	names(chr) <- genotypes[,2]
	Cm <- as.double(genotypes[,3])
	names(Cm) <- genotypes[,2]

	cross <- NULL
	for(i in unique(chr)){
		sum <- 0
		asum <- 0
		matrix <- NULL
		map <- NULL
		n <- NULL
		cm <- 0
		for(j in which(chr==i)){
			if(Cm[j] > (cm+2)){
				matrix <- rbind(matrix,as.integer(genomatrix[,j]))
				map <- rbind(map,Cm[j])
				n <- c(n,names(Cm[j]))
				sum <- sum+1
				cm <- Cm[j]
			}else{
				asum <- asum +1
			}
		}
		rownames(matrix) <- n
		rownames(map) <- n
		
		#Everything is okay now
		cross$geno[[i]]$data <- as.matrix(t(matrix))
		cross$geno[[i]]$map <- as.numeric(t(map))
		names(cross$geno[[i]]$map) <- n
		if(i != "X"){
			class(cross$geno[[i]])[1] <- "A"
		}else{
			class(cross$geno[[i]])[1] <- "X"
		}
	}
	cross$pheno <- as.data.frame(traits)
	class(cross)[1] <- "riself"
	class(cross)[2] <- "cross"
	cross
}


# end of prepareMQM.R
