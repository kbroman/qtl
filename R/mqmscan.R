#####################################################################
#
# mqmscan.R
#
# copyright (c) 2009, Danny Arends
# last modified Jun, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: mqmscan
#
######################################################################

######################################################################
#
# mqmscan: main scanning function to the MQMpackage
#
######################################################################
	
mqmscan <- function(cross,cofactors,pheno.col=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,file="MQM_output.txt",doLOG=0,est.map=0,dominance=0,plot=FALSE,forceRIL=0,verbose=FALSE){
    start <- proc.time()

	n.run=0
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
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

                        # check genotypes
                        g <- as.numeric(pull.geno(cross))
                        g <- sort(unique(g[!is.na(g)]))
                        if(max(g)==2) { # convert genotypes from 1/2 to 1/3
                          for(i in seq(along=cross$geno)) 
                            cross$geno[[i]]$data[!is.na(cross$geno[[i]]$data) & cross$geno[[i]]$data==2] <- 3
                        }

		}
		n.ind <- nind(cross)
		n.chr <- nchr(cross)
		if(verbose) {
                  cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
                  cat("INFO: Number of individuals: ",n.ind,"\n")
                  cat("INFO: Number of chromosomes: ",n.chr,"\n")
                }
		geno <- NULL
		chr <- NULL
		dist <- NULL
		out.qtl <- NULL	
		for(i in 1:n.chr) {
			geno <- cbind(geno,cross$geno[[i]]$data)
			chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
			dist <- c(dist,cross$geno[[i]]$map)
		}
		if(alfa <=0 || alfa >= 1){
			ourstop("Alfa must be between 0 and 1.\n")
		}
		#CHECK if the phenotype exists
		if (length(pheno.col) > 1){
			ourstop("For multiple phenotype analysis use the function: 'scanall'.\n")	
		}
		if(pheno.col != 1){
                  if(verbose) {
                    cat("INFO: Selected phenotype ",pheno.col,".\n")
                    cat("INFO: Number of phenotypes in object ",nphe(cross),".\n")
                  }
                  if(nphe(cross) < pheno.col || pheno.col < 1){
                    ourstop("No such phenotype in cross object.\n")
                  }			
		}
		pheno <- cross$pheno[[pheno.col]]
		if(var(pheno,na.rm = TRUE)> 1000){
			if(doLOG == 0){
				if(verbose) cat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"variation:",var(pheno,na.rm = TRUE),".\n")
				warning("INFO: Perhaps we should LOG-transform this phenotype, please set parameter: doLOG=1 to correct this error")
			}
		}
		if(doLOG != 0){
				#transform the cross file
				cross <- transformPheno(cross,pheno.col,transf=log)
				pheno <- cross$pheno[[pheno.col]]
		}
		n.mark <- ncol(geno)
		if(verbose) cat("INFO: Number of markers:",n.mark,"\n")
		Fril.replaced <- 0

		for(i in 1:n.ind) {
			for(j in 1:n.mark) {
				if(is.na(geno[i,j])){
					ourstop("Missing genotype information, please estimate unknown data, before running mqmscan.\n")
					geno[i,j] <- 9
				}else{
					if(forceRIL && ctype != 2 && geno[i,j]==2){
					#We have a 2 (AB) change it to a 3
					geno[i,j] <- 3
					Fril.replaced <- Fril.replaced+1
					}
				}
				
			}
		}
		if(verbose && forceRIL==1 && Fril.replaced > 0){
			 cat("INFO: Changed ",Fril.replaced," AB markers into BB markers.\n")
		}
		#check for missing phenotypes
		dropped <- NULL
		for(i in 1:length(pheno)) {
			if(is.na(pheno[i])){
			  if(verbose) cat("INFO: Dropped individual ",i," with missing phenotype.\n")
			  dropped <- c(dropped,i) 
			  n.ind = n.ind-1
			}
		}
		#throw em out
		if(!is.null(dropped)){
			geno <- geno[-dropped,]  
			pheno <- pheno[-dropped]
		}
		
		#CHECK for previously augmented dataset
		if(!is.null(cross$extra)){
		#	ourcat("INFO: previously augmented dataset.\n",a=verbose)			
		#	ourcat("INFO: Individuals before augmentation",cross$extra$Nind,".\n",a=verbose)
			extra1 <- cross$extra$Nind
		#	ourcat("INFO: Individuals after augmentation",cross$extra$augIND,".\n",a=verbose)
			extra2 <- cross$extra$augIND
		}else{
			#No augmentation so just set extra1 to be Nind (Naug internally of mqmscan)
			extra1 <- n.ind
			extra2 <- 0:n.ind
		}
		
		#CHECK if we have cofactors, so we can do backward elimination
		backward <- 0;
		if(missing(cofactors)){
			if(verbose) cat("INFO: No cofactors, setting cofactors to 0\n")
			cofactors = rep(0,n.mark)
		}else{
			if(length(cofactors) != n.mark){
				if(verbose) cat("ERROR: # Cofactors != # Markers\n")		
			}else{
				if(verbose) cat("INFO:",length(cofactors),"Cofactors received to be analyzed\n")
				if((sum(cofactors) > n.ind-10 && dominance==0)){
					ourstop("INFO: Cofactors don't look okay for use without dominance\n")
				}
				if((sum(cofactors)*2 > n.ind-10 && dominance==1)){
					ourstop("INFO: Cofactors don't look okay for use with dominance\n")
				}				
				if(sum(cofactors) > 0){
					if(verbose) cat("INFO: Doing backward elimination of selected cofactors.\n")
					backward <- 1;
					n.run <- 0;
				}else{
					backward <- 0;
					ourstop("Are u trying to give an empty cofactor list ???")
				}
			}
		}

		if((step.min+step.size) > step.max){
			ourstop("Current Step setting would crash the algorithm")
		}
		if(step.min>0){
			ourstop("step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
			ourstop("Step.size needs to be larger than 1")
		}
		max.cm.on.map <- max(unlist(pull.map(cross)))
		if(step.max < max.cm.on.map){
				ourstop("Markers outside of the mapping at ",max.cm.on.map," Cm, please set parameter step.max larger than this value.")		
		}
		qtlAchromo <- length(seq(step.min,step.max,step.size))
		if(verbose) cat("INFO: Number of locations per chromosome: ",qtlAchromo, "\n")
		end.1 <- proc.time()
		result <- .C("R_mqmscan",
				as.integer(n.ind),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
                as.integer(geno),
				as.integer(chr),
				DIST=as.double(dist),
				as.double(pheno),
				COF=as.integer(cofactors),
				as.integer(backward),
				as.integer(REMLorML),
				as.double(alfa),
				as.integer(em.iter),
				as.double(windowsize),
				as.double(step.size),
				as.double(step.min),
				as.double(step.max),
				as.integer(n.run),
				as.integer(extra1),
				as.integer(extra2),
				QTL=as.double(rep(0,2*n.chr*qtlAchromo)),
				as.integer(est.map),
				as.integer(ctype),
				as.integer(dominance),
				as.integer(verbose),
			    PACKAGE="qtl")
		end.2 <- proc.time()				
		# initialize output object
		qtl <- NULL
		info <- NULL
		names <- NULL
		for(i in 1:(n.chr*qtlAchromo)) {
			#Store the result in the qtl object
			qtl <- rbind(qtl,c(ceiling(i/qtlAchromo),rep(seq(step.min,step.max,step.size),n.chr)[i],result$QTL[i]))
			info <- rbind(info,result$QTL[(n.chr*qtlAchromo)+i])
			#make names in the form: C L
			names <- c(names,paste("C",ceiling(i/qtlAchromo),"L",rep(seq(step.min,step.max,step.size),n.chr)[i],sep=""))
		}
		if(plot){
			if(est.map && backward){
				op <- par(mfrow = c(3,1))
			}else{
				if(est.map || backward){
					op <- par(mfrow = c(2,1))
				}else{
					op <- par(mfrow = c(1,1))
				}
			}
			if(est.map){
				new.map <- pull.map(cross)
				aa <- nmar(cross)
				sum <- 1
				for(i in 1:length(aa)) {
					for(j in 1:aa[[i]]) {
						new.map[[i]][j] <- result$DIST[sum]
						sum <- sum+1
					}
				}
				if(verbose) cat("INFO: Viewing the user supplied map versus genetic map used during analysis.\n")
				plot.map(pull.map(cross), new.map,main="Supplied map versus re-estimated map")
			}
			if(backward){
				if(!est.map){
					new.map <- pull.map(cross)
				}
				aa <- nmar(cross)			
				sum <- 1
				model.present <- 0
				qc <- NULL
				qp <- NULL
				qn <- NULL
				for(i in 1:length(aa)) {
					for(j in 1:aa[[i]]) {
						#cat("INFO ",sum," ResultCOF:",result$COF[sum],"\n")
						if(result$COF[sum] != 48){
							if(verbose) cat("MODEL: Marker",sum,"from model found, CHR=",i,",POSITION=",as.double(unlist(new.map)[sum])," Cm\n")
							qc <- c(qc, as.character(names(cross$geno)[i]))
							qp <- c(qp, as.double(unlist(new.map)[sum]))
							qn <- c(qn, substr(names(unlist(new.map))[sum],3,nchar(names(unlist(new.map))[sum])))
							model.present <- 1
						}
						sum <- sum+1
					}
				}
				why <- sim.geno(cross)
				if(!is.null(qc) && model.present){
					QTLmodel <- makeqtl(why, qc, qp, qn, what="draws")
					plot(QTLmodel)
				}else{
					op <- par(mfrow = c(1,1))
				}
			}
		}
		rownames(qtl) <- names
		qtl <- cbind(qtl,1/(min(info))*(info-min(info)))
		qtl <- cbind(qtl,1/(min(info))*(info-min(info))*qtl[,3])
		colnames(qtl) = c("chr","pos (Cm)",paste("QTL",colnames(cross$pheno)[pheno.col]),"Info","QTL*INFO")
		#Convert to data/frame and scan.one object so we can use the standard plotting routines
		qtl <- as.data.frame(qtl)
		class(qtl) <- c("scanone",class(qtl)) 
		if(verbose) cat("INFO: Saving output to file: ",file, "\n")
		#write.table(qtl,file)
		#Reset plotting and return the results
		if(plot){
			info.c <- qtl
			#Check for errors in the information content IF err we can't do a second plot
			e <- 0
			for(i in 1:ncol(qtl)){
				if(is.na(info.c[i,5])){
					e<- 1
				}
				if(is.infinite(info.c[i,5])){
					e<- 1
				}
				if(is.null(info.c[i,5])){
					e<- 1
				}
			}
			#No error do plot 2
			if(!e){
				plotMQMone(qtl)
			}else{
				plot(qtl,lwd=1)
				grid(max(qtl$chr),5)
				labels <- paste("QTL",colnames(cross$pheno)[pheno.col])
				legend("topright", labels,col=c("black"),lty=c(1))
			}
		}
		#Reset the plotting window to contain 1 plot (fot the next upcomming pots
		if(plot){
		  op <- par(mfrow = c(1,1))
		}
		
		for( x in 1:nchr(cross)){
			#Remove pseudomarkers from the dataset and scale to the chromosome
			to.remove <- NULL
			chr.length <- max(cross$geno[[x]]$map)
			markers.on.chr <- which(qtl[,1]==x)
			to.remove <- markers.on.chr[which(qtl[markers.on.chr,2] > chr.length+step.size)]
			to.remove <- c(to.remove,markers.on.chr[which(qtl[markers.on.chr,2] < 0)])
			qtl <- qtl[-to.remove,]
		}
		
		end.3 <- proc.time()
		if(verbose) cat("INFO: Calculation time (R->C,C,C-R): (",round((end.1-start)[3], digits=3), ",",round((end.2-end.1)[3], digits=3),",",round((end.3-end.2)[3], digits=3),") (in seconds)\n")
		qtl
	}else{
		ourstop("Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

mqm <- function(...){
	mqmscan(...)
}

# end of mqmscan.R
