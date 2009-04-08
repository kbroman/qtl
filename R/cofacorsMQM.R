#####################################################################
#
# MQMCofactors.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: MQMCofactors
#
######################################################################

######################################################################
#
# MQMCofactors: Prepares a cofactor list to use with scanMQM
#
######################################################################

#setwd("D:/")
#library(qtl)
#dyn.load("scanMQM.dll")
#cross <- read.cross("csv","","Test.csv")

MQMCofactors <- function(cross= NULL,cofactors = NULL,sexfactors=NULL,verbose=TRUE){
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.")
	}
	if(is.null(cofactors)){
		ourstop("Cofactors to set. Please supply a list of markers to serve as cofactors.")
	}
		
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL
	individuals <- nind(cross)
	ourcat("INFO: Found",individuals,"individuals in the cross object.\n",a=verbose)
	ourcat("INFO: Mamimum amount of cofactors",(individuals-5),"leaves 5 Degrees of Freedom.\n",a=verbose)
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
	
	if(!is.null(sexfactors)){
		if(max(sexfactors) > n.mark){
			ourstop("Trying to set a non-existent marker as a sexfactor.")
		}
		if(min(sexfactors) <= 0){
			ourstop("Trying to set a non-existent marker as a sexfactor.")
		}
	}	

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

MQMCofactorsEach <- function(cross = NULL,each = 3,verbose=TRUE){
	if(is.null(cross)){
      ourstop("No cross file. Please supply a valid cross object.")
	  return 
	}

	individuals <- nind(cross)
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL

	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	n.mark <- ncol(geno)
	
	ourcat("INFO: Found",individuals,"individuals in the cross object.\n",a=verbose)
	ourcat("INFO: Mamimum amount of cofactors",(individuals-10)," (each =",ceiling(sum(n.mark)/(individuals-10)),") leaves 10 Degrees of Freedom.\n",a=verbose)


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
