#####################################################################
#
# transformPheno.R
#
# Copyright (c) 2009, Danny Arends
#
# Modified by Karl Broman
#
#
# first written Februari 2009
# last modified April 2009
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
# Contains: transformPheno
#
#####################################################################


transformPheno <- function(cross, pheno.col=1, transf=log, ...)
{
    #Helperfunction to transform a specific phenotype specified by the pheno.col parameter
    # by default, a log transformation is used, though one may use any function

    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(any(is.na(num))) {
            if(sum(is.na(num))>1)
                stop("Couldn't identify phenotypes ", paste(paste("\"", pheno.col[is.na(num)], "\"", sep=""),
                                                            collapse=" "))
            else
                stop("Couldn't identify phenotype \"", pheno.col[is.na(num)], "\"")
        }
        pheno.col <- num
    }

    if(any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")

    for(i in pheno.col)
        cross$pheno[,i] <- transf(cross$pheno[,i], ...)

    cross
}

# end of transformPheno.R
