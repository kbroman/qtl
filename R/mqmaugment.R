#####################################################################
#
# mqmaugment.R
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
# Contains: mqmaugment
#
######################################################################

######################################################################
#
# mqmaugment:
#
######################################################################
mqmaugment <- function(cross, pheno.col=1, maxind=1000, maxaugind=1, neglect=10, verbose=FALSE){
	start <- proc.time()
  maxaug = maxind   # maxaug is the maximum of individuals to augment to
  maxiaug = maxaugind

        # check for supported crosses and set ctype
        crosstype <- class(cross)[1]
        if(crosstype != "f2" && crosstype != "bc" && crosstype != "riself")
          ourstop("Currently only F2 / BC / RIL by selfing crosses can be analyzed by MQM.")

        if(crosstype == "f2"){
          ctype = 1
        }
        if(crosstype == "bc"){
          ctype = 2
        }
        if(crosstype == "riself"){
          ctype = 3
        }

        if(verbose) cat("INFO: Received a valid cross file type:", crosstype,".\n")

        # check whether the X chromosome should be dropped
        # (backcross with all one sex should be fine)
        chrtype <- sapply(cross$geno, class)
        if(any(chrtype == "X") && (crosstype=="f2" || 
                 length(getgenonames(crosstype, "X", "full",
                                     getsex(cross), attributes(cross))) != 2)) { # drop X chr
          warning("MQM not yet available for the X chromosome; omitting chr ",
                  paste(names(cross$geno)[chrtype == "X"], collapse=" "))
          cross <- subset(cross, chr=(chrtype != "X"))
        }

        n.ind <- nind(cross)
        n.chr <- nchr(cross)
        n.aug <- maxaug
        if(verbose) {
          cat("INFO: Number of individuals:",n.ind,".\n")
          cat("INFO: Number of chr:",n.chr,".\n")
        }
       
        # select the phenotype 
        if(length(pheno.col) > 1) {
           warning("Only one phenotype in pheno.col may be considered; using the first one.")
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
          if(verbose) {
            cat("INFO: Selected phenotype ",pheno.col," -> ",phenonaam,".\n")
            cat("INFO: # of phenotypes in object ",nphe(cross),".\n")
            }
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
		if(verbose) cat("INFO: Number of markers:",n.mark,".\n")

		# Check for NA genotypes and replace them with a 9
        geno[is.na(geno)] <- 9
        if(ctype==3) { 
          # RIL
          if(any(geno==3)) { # have 3's, so replace 2's with missing values
            if(any(geno==2)) {
              n2 <- sum(geno==2)
              geno[geno==2] <- 9
              warning("Removed ", n2, " het genotypes")
            }
          }
          else {
            if(any(geno==2)) # no 3's; replace 2's with 3's.
              geno[geno==2] <- 3
          }
        } # end if(RIL)

		# check for missing phenotypes and drop
		dropped <- NULL
		for(i in 1:dim(pheno)[1]) {
			if(is.na(pheno[i,pheno.col])){
				if(verbose) cat("INFO: Dropped individual ",i ," with missing phenotype.\n")
				dropped <- c(dropped,i) 
				n.ind = n.ind-1
			}
		}
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
		if(verbose) cat("INFO: DATA-Augmentation took: ",round((end-start)[3], digits=3)," seconds\n")
		cross

}

# end of mqmaugment.R
