#####################################################################
#
# mqmutil.R
#
# copyright (c) 2009, Danny Arends
# last modified Nov, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: ourstop, ourline
#
######################################################################

ourstop <- function(...){
	stop(...)
}

ourline <- function(){
	cat("------------------------------------------------------------------\n")		
}


mqmextractmarkers <- function(mqmresult){
  cleanedresult <- NULL
  for(x in 1:nrow(mqmresult)){
    if(!grepl(".loc",rownames(mqmresult)[x])){
      cleanedresult <- rbind(cleanedresult,mqmresult[x,])
    }
  }
  cleanedresult
}

mqmplot_directedqtl <- function(cross, mqmresults, draw = TRUE){
  onlymarkers <- mqmextractmarkers(mqmresults)
  eff <- effectscan(sim.geno(cross),draw=F)
  eff <- eff[-which(eff[,1]=="X"),]
  onlymarkers[,3] <- onlymarkers[,3]*(eff[,3]/abs(eff[,3]))
  if(draw) plot(ylim=c(-max(mqmresults[,3]),max(mqmresults[,3])),onlymarkers)
  onlymarkers
}

estimatemarkerlod <- function(interresults){
	for(x in 1:nrow(interresults)){
		if(is.na(interresults[x,3])){
			y <- x
			while(is.na(interresults[y,3])){
				y <- y + 1
			}
			nY <- interresults[y,3]
			nX <- interresults[y,2]
			distp = interresults[x,2] - pX
			distn = nX - interresults[x,2]
			disttot = distn+distp
			interresults[x,3] <- (((nY-pY)/disttot) * distp) + pY
			interresults[x,4] <- 1
			interresults[x,5] <- interresults[x,3]
		}
		pY <- interresults[x,3]		
		pX <- interresults[x,2]		
	}
	interresults
}


addmarkerstointervalmap <- function(cross,intervalresult,verbose=FALSE){
	map <- pull.map(cross)
	newres <- NULL
	intervalmaploc <- 1
	n <- NULL
	for(chr in 1:length(map)){
		for(mar in 1:length(map[[chr]])){
			if(verbose) cat(chr,"Placing marker: ",names(map[[chr]])[mar]," at ",map[[chr]][mar],"\t",intervalresult[intervalmaploc,2],"\n")
			if((class(map[[chr]])=="A")){
			while(intervalresult[intervalmaploc,2] < map[[chr]][mar]  || intervalresult[intervalmaploc,1] < chr ){
				newres <- rbind(newres,intervalresult[intervalmaploc,])
				n <- c(n,rownames(intervalresult)[intervalmaploc])
				intervalmaploc <- intervalmaploc+1
			}
			if(intervalresult[intervalmaploc,2] == map[[chr]][mar]){
				newres <- rbind(newres,c(chr,map[[chr]][mar],intervalresult[intervalmaploc,3],intervalresult[intervalmaploc,4],intervalresult[intervalmaploc,5]))
				while(intervalresult[intervalmaploc,2] == map[[chr]][mar]){
					intervalmaploc <- intervalmaploc+1
				}
			}else{
				newres <- rbind(newres,c(chr,map[[chr]][mar],NA,NA,NA))
			}
			
			n <- c(n,names(map[[chr]])[mar])
			colnames(newres) <- colnames(intervalresult)
			}
		}
	}
	if(intervalmaploc <= nrow(intervalresult)){
		newres <- rbind(newres,intervalresult[intervalmaploc,])
		n <- c(n,rownames(intervalresult)[intervalmaploc])
		intervalmaploc <- intervalmaploc+1		
	}
	rownames(newres) <- n
	newres <- estimatemarkerlod(newres)
	newres
}

mqmtestnormal <- function(cross, pheno.col=1){
	returnval <- FALSE
	if(pheno.col <0 || pheno.col > nphe(cross)){
		stop("No such phenotype (pheno.col = ",pheno.col,")")
	}
	if(!is.numeric(cross$pheno[[pheno.col]])){
		stop("Please supply a numeric trait (pheno.col = ",pheno.col," is not numeric)")
	}
	if(any(rownames(installed.packages())=="nortest")){
		require(nortest)
		if(pearson.test(cross$pheno[[pheno.col]])$p.value < 0.05){
			cat("Trait distribution not normal\n")
			returnval<- FALSE
		}else{
			cat("Trait distribution normal\n")
			returnval<- TRUE
		}
		returnval
	}else{
		cat("Please install package: nortest to enable testing of normality\n")
	}
}

mqmgetmodel <- function(scanresult){
	if(!is.null(scanresult)){
		model <- attr(scanresult,"mqmmodel")
		model
	}else{
		stop("Please supply a scan result made by using mqm qith cofactors")
	}
}

# end of mqmutil.R

