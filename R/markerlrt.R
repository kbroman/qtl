######################################################################
#
# markerlrt.R
#
# copyright (c) 2010, Karl W Broman
# last modified Jul, 2010
# first written Jul, 2010
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
# Contains: markerlrt
#
######################################################################

######################################################################
#
# markerlrt: General likelihood ratio test to assess linkage between
#            all pairs of markers
#
######################################################################

markerlrt <-
    function(cross)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    Geno <- pull.geno(cross)
    Geno[is.na(Geno)] <- 0
    maxg <- max(as.numeric(Geno))
    n.ind <- nrow(Geno)
    n.mar <- ncol(Geno)
    mar.names <- colnames(Geno)

    z <- .C("R_markerlrt",
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            as.integer(maxg),          # maximum observed genotype
            lod = as.double(rep(0,n.mar*n.mar)),
            PACKAGE="qtl")

    cross$rf <- matrix(z$lod,ncol=n.mar)
    dimnames(cross$rf) <- list(mar.names,mar.names)
    attr(cross$rf, "onlylod") <- TRUE

    cross
}

# end of markerlrt.R
