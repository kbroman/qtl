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
# mqmscan: main scanning function to the mqmpackage
#
######################################################################
	
mqmscan <- function(cross,cofactors,pheno.col=1,model=c("Additive","Dominance"),method=c("REML","ML"),
                    cof.significance=0.02,em.iter=1000,window.size=25.0,step.size=5.0,
                    step.min=-20.0,step.max=max(unlist(pull.map(cross))),logtransform = FALSE, estimate.map = FALSE,plot=FALSE,verbose=FALSE){
  
  start <- proc.time()
  method <- match.arg(method)
  model <- match.arg(model)
  step.max <- as.integer(step.max+step.size)
  #Because iirc we cannot pass booleans from R to C
  REMLorML <- 0             #We code :0 -> Restricted Maximum Likelyhood
  if(method=="ML"){
    REMLorML <- 1           #and 1 -> Maximum Likelyhood
  }
  dominance <- 0            #We code :0 -> Additive model (no_dominance)
  if(model=="Dominance"){
    dominance <- 1          #and 1 -> Dominance model
  }
  if(estimate.map){
    estimate.map <- 1
  }else{
    estimate.map <- 0
  }
	n.run <- 0
	if(is.null(cross)){
		stop("No cross file. Please supply a valid cross object.") 
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
		if(cof.significance <=0 || cof.significance >= 1){
			stop("cof.significance must be between 0 and 1.\n")
		}
    if(any(is.na(geno))){
			stop("Missing genotype information, please estimate unknown data, before running mqmscan.\n")
		}
		#CHECK if the phenotype exists
		if (length(pheno.col) > 1){
      warning("For a multiple phenotype analysis use the function: 'mqmall' for improved performance.\n")
      ##DANNY: HERE just call MQMall and pass the parameters
      cross$pheno <- cross$pheno[,pheno.col]   #Scale down the triats
      if(missing(cofactors)) cofactors <- rep(0,sum(nmar(cross)))
      result <- mqmall( cross,cofactors=cofactors,method=method,model=model,
                        cof.significance=cof.significance,step.min=step.min,step.max=step.max,step.size=step.size,window.size=window.size,
                        logtransform=logtransform, estimate.map = estimate.map,plot=plot, verbose=verbose)
			return(result)
		}
		if(pheno.col != 1){
      if(verbose) {
        cat("INFO: Selected phenotype ",pheno.col,".\n")
        cat("INFO: Number of phenotypes in object ",nphe(cross),".\n")
      }
      if(nphe(cross) < pheno.col || pheno.col < 1){
        stop("No such phenotype in cross object.\n")
      }			
		}

		if(any(rownames(installed.packages())=="nortest")){
			library(nortest)
			if(pearson.test(cross$pheno[[pheno.col]])$p.value < 0.05){
				warning("Trait might not be normal (pearsons normallity test)\n")
			}
		}
		pheno <- cross$pheno[[pheno.col]]
		if(var(pheno,na.rm = TRUE)> 1000){
			if(!logtransform){
				if(verbose) cat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"variation:",var(pheno,na.rm = TRUE),".\n")
				warning("INFO: Perhaps we should LOG-transform this phenotype, please set parameter: logtransform=1 to correct this error")
			}
		}
		if(logtransform){
				#transform the cross file
				cross <- transformPheno(cross,pheno.col,transf=log)
				pheno <- cross$pheno[[pheno.col]]
		}
		n.mark <- ncol(geno)
		if(verbose) cat("INFO: Number of markers:",n.mark,"\n")

		#check for missing phenotypes
		dropped <- NULL
		for(i in 1:length(pheno)) {
			if(is.na(pheno[i])){
			  if(verbose) cat("INFO: Dropped individual ",i," with missing phenotype.\n")
			  dropped <- c(dropped,i) 
			  n.ind = n.ind-1
			}
		}
		#Throw out missing phenotypes from phenotype vector and genotype matrix
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
				if(verbose) cat("INFO:",sum(cofactors!=0),"Cofactors received to be analyzed\n")
				if((sum(cofactors) > n.ind-10 && dominance==0)){
					stop("INFO: Cofactors don't look okay for use without dominance\n")
				}
				if((sum(cofactors)*2 > n.ind-10 && dominance==1)){
					stop("INFO: Cofactors don't look okay for use with dominance\n")
				}				
				if(sum(cofactors) > 0){
					if(verbose) cat("INFO: Doing backward elimination of selected cofactors.\n")
					backward <- 1;
					n.run <- 0;
				}else{
					backward <- 0;
					cofactors = rep(0,n.mark)
				}
			}
		}

		if((step.min+step.size) > step.max){
			stop("Current Step setting would crash the algorithm")
		}
		if(step.min>0){
			stop("step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
			stop("Step.size needs to be larger than 1")
		}
		max.cm.on.map <- max(unlist(pull.map(cross)))
		if(step.max < max.cm.on.map){
				stop("Markers outside of the mapping at ",max.cm.on.map," Cm, please set parameter step.max larger than this value.")		
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
				as.double(cof.significance),
				as.integer(em.iter),
				as.double(window.size),
				as.double(step.size),
				as.double(step.min),
				as.double(step.max),
				as.integer(n.run),
				as.integer(extra1),
				as.integer(extra2),
				QTL=as.double(rep(0,2*n.chr*qtlAchromo)),
				as.integer(estimate.map),
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
			#make names in the form: cX.locXX
			names <- c(names,paste("c",ceiling(i/qtlAchromo),".loc",rep(seq(step.min,step.max,step.size),n.chr)[i],sep=""))
		}
		if(plot){
			if(estimate.map && backward){
				op <- par(mfrow = c(3,1))
			}else{
				if(estimate.map || backward){
					op <- par(mfrow = c(2,1))
				}else{
					op <- par(mfrow = c(1,1))
				}
			}
			if(estimate.map){
				new.map <- pull.map(cross)
				chrmarkers <- nmar(cross)
				sum <- 1
				for(i in 1:length(chrmarkers)) {
					for(j in 1:chrmarkers[[i]]) {
						new.map[[i]][j] <- result$DIST[sum]
						sum <- sum+1
					}
				}
				if(verbose) cat("INFO: Viewing the user supplied map versus genetic map used during analysis.\n")
				plot.map(pull.map(cross), new.map,main="Supplied map versus re-estimated map")
			}
    }
    if(backward){
      if(!estimate.map){
        new.map <- pull.map(cross)
      }
      chrmarkers <- nmar(cross)			
      sum <- 1
      model.present <- 0
      qc <- NULL
      qp <- NULL
      qn <- NULL
      for(i in 1:length(chrmarkers)) {
        for(j in 1:chrmarkers[[i]]) {
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
      if(!is.null(qc) && model.present){
		why <- sim.geno(cross,n.draws=1)
        QTLmodel <- makeqtl(why, qc, qp, qn, what="draws")
		if(plot){
			plot(QTLmodel)
        }
      }
    }
	rownames(qtl) <- names
	qtl <- cbind(qtl,1/(min(info))*(info-min(info)))
	qtl <- cbind(qtl,1/(min(info))*(info-min(info))*qtl[,3])
	colnames(qtl) = c("chr","pos (Cm)",paste("LOD",colnames(cross$pheno)[pheno.col]),"info","LOD*info")
	#Convert to data/frame and scan.one object so we can use the standard plotting routines
	qtl <- as.data.frame(qtl)
	if(backward && !is.null(qc) && model.present){
	  attr(qtl,"mqmmodel") <- QTLmodel
	}
	class(qtl) <- c("scanone",class(qtl)) 
	for( x in 1:nchr(cross)){
		#Remove pseudomarkers from the dataset and scale to the chromosome
		to.remove <- NULL
		chr.length <- max(cross$geno[[x]]$map)
		markers.on.chr <- which(qtl[,1]==x)
		to.remove <- markers.on.chr[which(qtl[markers.on.chr,2] > chr.length+step.size)]
		to.remove <- c(to.remove,markers.on.chr[which(qtl[markers.on.chr,2] < 0)])
		qtl <- qtl[-to.remove,]
	}		
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
			mqmplotone(qtl,main=paste(colnames(cross$pheno)[pheno.col],"at significance=",cof.significance))
		}else{
			plot(qtl,main=paste(colnames(cross$pheno)[pheno.col],"at significance=",cof.significance),lwd=1)
			grid(max(qtl$chr),5)
			labels <- paste("QTL",colnames(cross$pheno)[pheno.col])
			legend("topright", labels,col=c("black"),lty=c(1))
		}
	}
	#Reset the plotting window to contain 1 plot (fot the next upcomming pots
	if(plot){
	  op <- par(mfrow = c(1,1))
	}
	end.3 <- proc.time()   
	if(verbose) cat("INFO: Calculation time (R->C,C,C-R): (",round((end.1-start)[3], digits=3), ",",round((end.2-end.1)[3], digits=3),",",round((end.3-end.2)[3], digits=3),") (in seconds)\n")
	qtl
	}else{
		stop("Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

mqm <- function(...){
	mqmscan(...)
}

# end of mqmscan.R
