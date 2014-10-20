######################################################################
#
# add.cim.covar.R
#
# copyright (c) 2007-8, Karl W Broman
# last modified Aug, 2008
# first written Mar, 2007
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
# Contains: add.cim.covar
#
######################################################################

######################################################################
# Indicate the locations of the selected marker covariates
# in a plot of CIM results (as obtained by plot.scanone)
#
# The chr and gap arguments must be identical to those used in the
# call to plot.scanone.
######################################################################
add.cim.covar <-
    function(cimresult, chr, gap=25, ...)
{
    cimcovar <- attr(cimresult, "marker.covar.pos")

    if(!missing(chr)) cimresult <- subset(cimresult, chr=chr)
    if(nrow(cimcovar) == 0) return(invisible(NULL))
    chr <- as.character(unique(cimresult[,1]))

    dots <- list(...)
    ndots <- names(dots)

    u <- par("usr")
    if(length(chr)==1) {
        if(!("col" %in% ndots) && !("pch" %in% ndots))
            points(cimcovar[,2], u[3], xpd=TRUE, col="red", pch=16, ...)
        else if(!("col" %in% ndots))
            points(cimcovar[,2], u[3], xpd=TRUE, col="red", ...)
        else if(!("pch" %in% ndots))
            points(cimcovar[,2], u[3], xpd=TRUE, pch=16, ...)
        else
            points(cimcovar[,2], u[3], xpd=TRUE, ...)
    }
    else {
        begend <- matrix(unlist(tapply(cimresult[,2],cimresult[,1],range)),ncol=2,byrow=TRUE)
        rownames(begend) <- chr
        begend <- begend[as.character(chr),,drop=FALSE]
        len <- begend[,2]-begend[,1]
        start <- c(0,cumsum(len+gap))-c(begend[,1],0)
        start <- start[-length(start)]
        names(start) <- chr

        for(i in 1:nrow(cimcovar)) {
            x <- start[cimcovar[i,1]] + cimcovar[i,2]
            if(!("col" %in% ndots) && !("pch" %in% ndots))
                points(x, u[3], xpd=TRUE, col="red", pch=16, ...)
            else if(!("col" %in% ndots))
                points(x, u[3], xpd=TRUE, col="red", ...)
            else if(!("pch" %in% ndots))
                points(x, u[3], xpd=TRUE, pch=16, ...)
            else
                points(x, u[3], xpd=TRUE, ...)
        }
    }

    invisible(cimcovar)
}

# end of add.cim.covar.R
