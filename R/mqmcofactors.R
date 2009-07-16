#####################################################################
#
# MQMCofactors.R
#
# copyright (c) 2009, Danny Arends
# last modified Jun, 2009
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
# Contains: MQMCofactors, MQMCofactorsEach
#
######################################################################

######################################################################
#
# MQMCofactors: Prepares a cofactor list to use with mqmscan
# MQMCofactorsEach: Prepares a cofactor list to use with mqmscan
#
######################################################################

MQMCofactors <- function(cross,cofactors,sexfactors,verbose=FALSE){
	if(missing(cross))
		ourstop("No cross file. Please supply a valid cross object.")

	if(missing(cofactors))
		ourstop("Cofactors to set. Please supply a list of markers to serve as cofactors.")
		
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL
	individuals <- nind(cross)
        if(verbose) {
          cat("INFO: Found",individuals,"individuals in the cross object.\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-10),"leaves 10 Degrees of Freedom, (No Dominance).\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-10)/2,"leaves 10 Degrees of Freedom (Dominance).\n")
        }
	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	
	n.mark <- ncol(geno)

	if(max(cofactors) > n.mark){
		ourstop("Trying to set a non-existent marker as a cofactor.")	  
	}

	if(min(cofactors) <= 0){
		ourstop("Trying to set a non-existent marker as a cofactor.")	  
	}
	
	if(!missing(sexfactors)){
		if(max(sexfactors) > n.mark){
			ourstop("Trying to set a non-existent marker as a sexfactor.")
		}
		if(min(sexfactors) <= 0){
			ourstop("Trying to set a non-existent marker as a sexfactor.")
		}
	}
        else sexfactors <- NULL

    cofactorlist <- rep(0,n.mark)
	for(i in 1:length(cofactors)) {
	  cofactorlist[cofactors[i]]=1
	}
	if(!is.null(sexfactors)){
      for(i in 1:length(sexfactors)) {
	    cofactorlist[sexfactors[i]]=2
	  }
	}
    if(sum(cofactorlist) > (individuals-10)){
		ourstop("Trying to set: ",sum(cofactorlist)," markers as cofactor. This leaves less than 10 Degrees of Freedom.\n")
		return
	}
    cofactorlist
}

MQMCofactorsEach <- function(cross,each = 3,verbose=FALSE){
	if(missing(cross))
          ourstop("No cross file. Please supply a valid cross object.")

	individuals <- nind(cross)
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL

	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	n.mark <- ncol(geno)
	
        if(verbose) {
          cat("INFO: Found",individuals,"individuals in the cross object.\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-10)," (each =",ceiling(sum(n.mark)/(individuals-10)),") leaves 10 Degrees of Freedom (no Dominance).\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-10)/2," (each =",ceiling(sum(n.mark)/(individuals-10))*2,") leaves 10 Degrees of Freedom (Dominance).\n")
        }

	if(each > n.mark){
      ourstop("Not enough markers to place cofactors at.")
	  return 
	}	
	
    cofactorlist <- rep(0,n.mark)
	for(i in 1:n.mark) {
		if(i%%as.integer(each)==0){
			cofactorlist[i] = 1
		}
	}
    if(sum(cofactorlist) > (individuals-10)){
		ourstop("Trying to set: ",ceiling(sum(n.mark)/each)," markers as cofactor. This leaves less than 10 Degrees of Freedom.\n")
		return
	}
	cofactorlist
}

#a <- MQMCofactors(cross,c(10,20,30,40,50,60,70,80),c(186,187))

# end of cofactorsMQM.R
