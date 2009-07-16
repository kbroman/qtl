#####################################################################
#
# mqmsnow.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Apr, 2009
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
# Contains: snowCoreALL, snowCoreBOOT, snowCore
#
######################################################################


snowCoreALL <- function(x,all.data,Funktie,verbose=FALSE,...){
	b <- proc.time()
	result <- NULL
	num.traits <- nphe(all.data)
        if(verbose) {
          cat("------------------------------------------------------------------\n")
          cat("INFO: Starting analysis of trait (",x,"/",num.traits,")\n")
          cat("------------------------------------------------------------------\n")
        }
	if("cofactors" %in% names(formals(Funktie))){
		if(exists("cofactors")){
			result <- Funktie(cross=all.data,cofactors=cofactors,pheno.col=x,verbose=FALSE,...)
		}
	}else{
		result <- Funktie(cross=all.data,pheno.col=x,...)
	}
	colnames(result)[3] <- paste("lod",names(all.data$pheno)[x])
	e <- proc.time()
        if(verbose) {
          cat("------------------------------------------------------------------\n")
          cat("INFO: Done with the analysis of trait (",x,"/",num.traits,")\n")	
          cat("INFO: Calculation of trait",x,"took:",round((e-b)[3], digits=3)," seconds\n")
          cat("------------------------------------------------------------------\n")
        }
	result
}


snowCoreBOOT <- function(x,all.data,Funktie,bootmethod,verbose=FALSE,...){
	b <- proc.time()
	result <- NULL
	if(!bootmethod){
		#random permutation
		neworder <- sample(nind(all.data))			
		all.data$pheno[[1]] <- all.data$pheno[[1]][neworder]
	}else{
		#parametric permutation
		pheno <- all.data$pheno[[1]]
		variance <- var(pheno,na.rm = TRUE)
		for(j in 1:nind(all.data)) {
			all.data$pheno[[1]][j] <- rnorm(1)*(variance^0.5)
		}
	}
	if("cofactors" %in% names(formals(Funktie))){
		if(exists("cofactors")){
			result <- Funktie(cross=all.data,cofactors=cofactors,pheno.col=1,verbose=FALSE,...)
		}else{
			result <- Funktie(cross=all.data,pheno.col=1,verbose=FALSE,...)
		}
	}else{
		if("plot" %in% names(formals(Funktie))){
			result <- Funktie(cross=all.data,pheno.col=1,...)
		}else{
			result <- Funktie(cross=all.data,pheno.col=1,...)
		}
	}
	e <- proc.time()
        if(verbose) {
          cat("------------------------------------------------------------------\n")
          cat("INFO: Done with bootstrap\n")
          cat("INFO: Calculation took:",round((e-b)[3], digits=3)," seconds\n")
          cat("------------------------------------------------------------------\n")
        }
	result
}

# end of mqmsnow.R
