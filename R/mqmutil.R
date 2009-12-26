#####################################################################
#
# mqmutil.R
#
# Copyright (c) 2009, Danny Arends
#
# Modified by Pjotr Prins and Karl Broman
#
# 
# first written Februari 2009
# last modified December 2009
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
# Contains: ourstop
#           ourline
#           mqmextractmarkers
#           estimatemarkerlod
#           addmarkerstointervalmap
#           mqmtestnormal
#           mqmgetmodel
#           
#
#####################################################################





ourstop <- function(...){
	stop(...)
}

ourline <- function(){
	cat("------------------------------------------------------------------\n")		
}


# Return the real markers in the set (remove fake ones)
mqmextractmarkers <- function(mqmresult){
  result <- NULL
  for(x in 1:nrow(mqmresult)){
    # for every marker...
    marker = mqmresult[x,]
    found = grep('.loc',rownames(marker))
    if (length(found)==0) {
      result <- rbind(result,marker)
    }
  }
  class(result) <- class(mqmresult)
  result
}


estimatemarkerlod <- function(interresults){
  #For an okay return, with all markers filled every REAL marker has to be surrounded by interval markers
  #It does skip markers untill we reach the next pseudomarker. When one of the assumptions fails
  #we return, so there could be markers without a LOD score
  if(all(is.na(interresults[,3]))) return (interresults) 
	pY <- interresults[1,3]		
	pX <- interresults[1,2]
  if(is.na(pY) || is.na(pX)) return (interresults) #The first marker needs to be a interval marker
  for(x in 2:nrow(interresults)){
		if(is.na(interresults[x,3])){
			y <- x
			while(y <= nrow(interresults) && is.na(interresults[y,3]) && interresults[y,1] == interresults[x,1]){
				y <- y + 1
			}
			nY <- interresults[y,3]     
			nX <- interresults[y,2]
      if(is.na(nY) || is.na(nX)) return (interresults) #The next marker also needs to be a interval marker
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
	class(newres) <- c("scanone",class(newres))
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

