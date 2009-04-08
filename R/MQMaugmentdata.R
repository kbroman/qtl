#####################################################################
#
# MQMaugmentdata.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: MQMaugment
#
######################################################################

######################################################################
#
# MQMaugmentdata:
#
######################################################################
#setwd("D:/")
#library(qtl)
#dyn.load("scanMQM.dll")
#qtl <- c(3,15,3,7)							# QTL at chromosome 3
#data(map10)									# Mouse genome
#cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross
#data(listeria)

MQMaugment <- function(cross= NULL,pheno.col=1,maxaug=1000,maxiaug=10,neglect=10,verbose=TRUE){
	start <- proc.time()
	library(qtl)
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.")
		return 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		if(class(cross)[1] == "f2"){
			ctype = 1
		}
		if(class(cross)[1] == "bc"){
			ctype = 2
		}
		if(class(cross)[1] == "riself"){
			ctype = 3
		#	stop("Somethings still wrong in the algorithm, please analyse RIL as BC.")
		}
		ourcat("INFO: Received a valid cross file type:",class(cross)[1],".\n",a=verbose)
		n.ind <- nind(cross)
		n.chr <- nchr(cross)
		n.aug <- maxaug
		ourcat("INFO: Number of individuals:",n.ind,".\n",a=verbose)
		ourcat("INFO: Number of chr:",n.chr,".\n",a=verbose)
		phenonaam <- colnames(cross$pheno)[pheno.col]
		geno <- NULL
		chr <- NULL
		dist <- NULL
		out.qtl <- NULL	
		for(i in 1:n.chr) {
			geno <- cbind(geno,cross$geno[[i]]$data)
			chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
			dist <- c(dist,cross$geno[[i]]$map)
		}
		pheno <- cross$pheno
		n.mark <- ncol(geno)
		ourcat("INFO: Number of markers:",n.mark,".\n",a=verbose)
		#Check for na genotypes and replace them with a 9
		for(i in 1:n.ind) {
			for(j in 1:n.mark) {
				if(is.na(geno[i,j])){
					geno[i,j] <- 9;
				}
			}
		}
		#check for missing phenotypes
		dropped <- NULL
		for(i in 1:dim(pheno)[1]) {
			if(is.na(pheno[i,pheno.col])){
				ourcat("INFO: Dropped individual ",i ," with missing phenotype.\n",a=verbose)
				dropped <- c(dropped,i) 
				n.ind = n.ind-1
			}
		}
		#throw em out
		if(!is.null(dropped)){
			geno <- geno[-dropped,]  
			pheno <- pheno[-dropped,]
		}
		if(pheno.col != 1){
			
			ourcat("INFO: Selected phenotype ",pheno.col," -> ",phenonaam,".\n",a=verbose)
			ourcat("INFO: # of phenotypes in object ",nphe(cross),".\n",a=verbose)
			if(nphe(cross) < pheno.col || pheno.col < 1){
				ourstop("No such phenotype at column index:",pheno.col,"in cross object.\n")
			}			
		}
		result <- .C("R_augdata",
				as.integer(geno),
				as.double(dist),
				as.double(pheno[,pheno.col]),
				augGeno=as.integer(rep(0,n.mark*maxaug)),
				augPheno=as.double(rep(0,maxaug)),
				augIND=as.integer(rep(0,maxiaug*n.ind)),
				as.integer(n.ind),
				as.integer(n.aug),
				as.integer(n.mark),
				as.integer(1),    # 1 phenotype
				as.integer(maxaug),
				as.integer(maxiaug),
				as.double(neglect),
				as.integer(chr),
				as.integer(ctype)
				)
		n.ind = result[[7]]
		n.aug = result[[8]]	
		markONchr <- 0
		markdone <- 0
		for(c in 1:n.chr){
			#print(paste("Cromosome",c,"\n",sep=""))
			matri <- NULL
			matri2 <- NULL
			markONchr <- dim(cross$geno[[c]]$data)[2]
			#print(paste("# markers",markONchr,"\n",sep=""))
			for(j in markdone:(markdone+markONchr-1)){
				#print(paste("Start",markdone,":End",(markdone+markONchr-1),"\n",sep=""))
				ind2 <- NULL
				pheno <- NULL
				ind2 = result[[4]][(1+(j*maxaug)):(n.aug+(j*maxaug))]
				matri <- rbind(matri,ind2)
			}
			pheno <- as.matrix(result[[5]][1:n.aug])
			matri <- t(matri)
			if(markdone==0){
				colnames(matri) <- colnames(geno)[markdone:(markdone+markONchr)]
			}else{
				#print(paste("Markdone",markdone,"End",(markdone+markONchr-1)))
				colnames(matri) <- colnames(geno)[(markdone+1):(markdone+markONchr)]			
			}
			cross$geno[[c]]$data <- matri
			colnames(pheno) = phenonaam
			cross$pheno <- as.data.frame(pheno)
			#Store extra information (needed by the MQM algorithm) which individual was which original etc..
			cross$extra$Nind <- n.ind
			cross$extra$Naug <- n.aug
			cross$extra$augIND <- result[[6]][1:n.aug]
			markdone <- (markdone+markONchr)  
		}
		#RETURN THE RESULTS
		end <- proc.time()
		cat("INFO: DATA-Augmentation took: ",round((end-start)[3], digits=3)," seconds\n")		
		cross
	}else{
		ourstop("Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

MQMlogPheno <- function(cross= NULL,pheno.col=NULL,verbose=TRUE){
	#Helperfunction to logtransform a specific phenotype specified by the Phenot parameter
	library(qtl)
	if(is.null(cross)){
		ourstop("ERROR: No cross file. Please supply a valid cross object.")
		return 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		if(!is.null(pheno.col)){
			pheno <- NULL
			logpheno <- NULL
			pheno <- cross$pheno[[pheno.col]]
			logpheno <- log(pheno)
			cross$pheno[[pheno.col]] <- logpheno
			ourcat("INFO: Phenotype:",names(cross$pheno)[pheno.col],".\n",a=verbose)
			ourcat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"Variation:",var(pheno,na.rm = TRUE),".\n",a=verbose)
			ourcat("INFO: After LOG transformation Mean:",mean(logpheno,na.rm = TRUE),"Variation:",var(logpheno,na.rm = TRUE),".\n",a=verbose)
		}else{
			n.pheno <- nphe(cross)
			for(i in 1:n.pheno) {
				pheno <- cross$pheno[[i]]
				logpheno <- log(pheno)
				cross$pheno[[i]] <- logpheno
				ourcat("INFO: Phenotype:",names(cross$pheno)[i],".\n",a=verbose)
				ourcat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"Variation:",var(pheno,na.rm = TRUE),".\n",a=verbose)
				ourcat("INFO: After LOG transformation Mean:",mean(logpheno,na.rm = TRUE),"Variation:",var(logpheno,na.rm = TRUE),".\n",a=verbose)				
			}
		}
		cross
	}else{
		ourstop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

#cross_good <- MQMaugment(cross)
#listeria_good <- MQMaugment(listeria,maxaug=1000,maxiaug=10,neglect=10)

# end of MQMaugment.R
