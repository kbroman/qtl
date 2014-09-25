######################################################################
#
# countXO.R
#
# copyright (c) 2008-2013, Karl W Broman
# last modified Sep, 2013
# first written Feb, 2008
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
# Contains: countXO
#
######################################################################

######################################################################
#
# countXO: Count number of obligate crossovers for each individual
#          on individual chromosomes or overall
#
# if bychr=TRUE, return matrix with no. obligate crossovers for each
#                individual on each chromosome
#        =FALSE, return vector with total no. crossovers across the
#                selected chromosomes
######################################################################

countXO <-
    function(cross, chr, bychr=FALSE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    # pull out relevant chromosome
    if(!missing(chr)) cross <- subset(cross,chr=chr)
    chr.name <- names(cross$geno)

    type <- class(cross)[1]
    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    nxo <- matrix(0, ncol=n.chr, nrow=n.ind)
    id <- getid(cross)
    if(is.null(id)) id <- 1:n.ind
    dimnames(nxo) <- list(id, chr.name)

    for(i in 1:n.chr) {

        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") xchr <- TRUE
        else xchr <- FALSE

        # which type of cross is this?
        if(type == "f2") {
            if(!xchr) # autosomal
                func <- "R_countXO_f2"
            else func <- "R_countXO_bc"        # X chromsome
        }
        else if(type == "bc" || type=="riself" || type=="risib" || type=="dh" || type=="haploid") func <- "R_countXO_bc"
        else if(type == "4way") func <- "R_countXO_4way"
        else if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="bgmagic16") {
            func <- "R_countXO_ril48"
            if(xchr)
                warning("countXO not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else
            stop("ripple not available for cross ", type)

        # data to be input
        genodat <- cross$geno[[i]]$data
        genodat[is.na(genodat)] <- 0
        n.mar <- ncol(genodat)

        if(n.mar > 1) {
            z <- .C(func,
                    as.integer(n.ind),
                    as.integer(n.mar),
                    as.integer(genodat),
                    oblxo=as.integer(rep(0,n.ind)),
                    PACKAGE="qtl")

            nxo[,i] <- z$oblxo
        }
    }

    if(!bychr) nxo <- apply(nxo, 1, sum)

    nxo
}

# end of countXO.R
