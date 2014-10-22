#####################################################################
#
# mqmsnow.R
#
# Copyright (c) 2009-2010, Danny Arends
#
# Modified by Karl Broman and Pjotr Prins
#
#
# first written Februari 2009
# last modified April 2010
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
# Contains: snowCoreALL
#           snowCoreBOOT
#
#
#####################################################################






snowCoreALL <- function(x,all.data,scanfunction,cofactors,verbose=FALSE,...){
    b <- proc.time()
    result <- NULL
    num.traits <- nphe(all.data)
    if(verbose) {
        cat("------------------------------------------------------------------\n")
        cat("INFO: Starting analysis of trait (",x,"/",num.traits,")\n")
        cat("------------------------------------------------------------------\n")
    }
    if(!is.null(cofactors)) {
        result <- scanfunction(cross=all.data,cofactors=cofactors,pheno.col=x,verbose=verbose,...)
    }
    else{
        result <- scanfunction(cross=all.data,pheno.col=x,...)
    }
    colnames(result)[3] <- paste("LOD",names(all.data$pheno)[x])
    e <- proc.time()
    if(verbose) {
        cat("------------------------------------------------------------------\n")
        cat("INFO: Done with the analysis of trait (",x,"/",num.traits,")\n")
        cat("INFO: Calculation of trait",x,"took:",round((e-b)[3], digits=3)," seconds\n")
        cat("------------------------------------------------------------------\n")
    }
    result
}


snowCoreBOOT <- function(x,all.data,scanfunction,bootmethod,cofactors,verbose=FALSE,...){
    b <- proc.time()
    result <- NULL
    if(!bootmethod){
        #random permutation
        neworder <- sample(nind(all.data))
        all.data$pheno[[1]] <- all.data$pheno[[1]][neworder]
    }else{
        #parametric permutation
        all.data$pheno[[1]] <- rnorm(nind(all.data))
    }
    if("cofactors" %in% names(formals(scanfunction))){
        if(!is.null(cofactors)){
            result <- scanfunction(cross=all.data,cofactors=cofactors,pheno.col=1,verbose=FALSE,...)
        }else{
            result <- scanfunction(cross=all.data,pheno.col=1,verbose=FALSE,...)
        }
    }else{
        if("plot" %in% names(formals(scanfunction))){
            result <- scanfunction(cross=all.data,pheno.col=1,...)
        }else{
            result <- scanfunction(cross=all.data,pheno.col=1,...)
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
