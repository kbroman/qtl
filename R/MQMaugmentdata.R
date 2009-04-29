#####################################################################
#
# MQMaugmentdata.R
#
# copyright (c) 2009, Danny Arends and Karl W. Broman
# last modified Apr, 2009
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
# Contains: MQMaugment
#
######################################################################

######################################################################
#
# MQMaugmentdata:
#
######################################################################
MQMaugment <- function(cross, pheno.col=1, maxaug=1000, maxiaug=10, neglect=10, verbose=FALSE){
	start <- proc.time()

        crosstype <- class(cross)[1]
        if(crosstype != "f2" && crosstype != "bc" && crosstype != "riself")
          ourstop("Currently only F2 / BC / RIL by selfing crosses can be analyzed by MQM.")

        if(class(cross)[1] == "f2"){
          ctype = 1
        }
        if(class(cross)[1] == "bc"){
          ctype = 2
        }
        if(class(cross)[1] == "riself"){
          ctype = 3
        }

        ourcat("INFO: Received a valid cross file type:", crosstype,".\n",a=verbose)

        n.ind <- nind(cross)
        n.chr <- nchr(cross)
        n.aug <- maxaug
        ourcat("INFO: Number of individuals:",n.ind,".\n",a=verbose)
        ourcat("INFO: Number of chr:",n.chr,".\n",a=verbose)

        
        if(length(pheno.col) > 1) {
           warning("Only one phenotype may be considered; using the first one.")
           pheno.col <- pheno.col[1]
         }
        if(is.character(pheno.col)) {
          num <- find.pheno(cross, pheno.col)
          if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
          pheno.col <- num
        }
        phenonaam <- colnames(cross$pheno)[pheno.col]

        if(pheno.col != 1){
			
          ourcat("INFO: Selected phenotype ",pheno.col," -> ",phenonaam,".\n",a=verbose)
          ourcat("INFO: # of phenotypes in object ",nphe(cross),".\n",a=verbose)
            if(nphe(cross) < pheno.col || pheno.col < 1){
              ourstop("No such phenotype at column index:",pheno.col,"in cross object.\n")
            }			
        }


		out.qtl <- NULL	

        geno <- pull.geno(cross)
        chr <- rep(1:nchr(cross), nmar(cross))
        dist <- unlist(pull.map(cross))

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
		result <- .C("R_augdata",
				as.integer(geno),
				as.double(dist),
				as.double(pheno[,pheno.col]),
				augGeno=as.integer(rep(0,n.mark*maxaug)),
				augPheno=as.double(rep(0,maxaug)),
				augIND=as.integer(rep(0,maxiaug*n.ind)),
				nind=as.integer(n.ind),
				naug=as.integer(n.aug),
				as.integer(n.mark),
				as.integer(1),    # 1 phenotype
				as.integer(maxaug),
				as.integer(maxiaug),
				as.double(neglect),
				as.integer(chr),
				as.integer(ctype),
                                as.integer(verbose),
				PACKAGE="qtl")
		n.ind = result$nind
		n.aug = result$naug
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
				ind2 <- result[[4]][(1+(j*maxaug)):(n.aug+(j*maxaug))]
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
			colnames(pheno) <- phenonaam
			cross$pheno <- as.data.frame(pheno)
			#Store extra information (needed by the MQM algorithm) which individual was which original etc..
			cross$extra$Nind <- n.ind
			cross$extra$Naug <- n.aug
			cross$extra$augIND <- result[[6]][1:n.aug]
			markdone <- (markdone+markONchr)  
		}
		#RETURN THE RESULTS
		end <- proc.time()
		ourcat("INFO: DATA-Augmentation took: ",round((end-start)[3], digits=3)," seconds\n", a=verbose)		
		cross

}

# end of MQMaugment.R
