#####################################################################
#
# inferFounderHap.R
#
# copyright (c) 2011, Karl W Broman
# last modified Dec, 2011
# first written Dec, 2011
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
# Contains: inferFounderHap, restoreMWrilGeno
#
# This is for reconstructing the founder haplotypes in inbred lines
# by a crude method using groups of adjacent SNPs
#
######################################################################

inferFounderHap <-
    function(cross, chr, max.n.markers=15)
{
    if(!missing(chr))
        cross <- subset(cross, chr=chr)
    if(nchr(cross) > 1) {
        thechr <- names(cross$geno)[1]
        cross <- subset(cross, chr=thechr)
        warning("inferFounderHap is only for one chromosome; considering ", thechr)
    }

    # pull out genotypes for RIL and founders
    offspringGen <- restoreMWrilGeno(cross)
    founderGen <- cross$founderGeno

    # drop markers with any missing data in the founders
    nomissing <- apply(founderGen, 2, function(a) !any(is.na(a)))
    names(nomissing) <- colnames(offspringGen)
    if(!any(nomissing))
        stop("No markers with complete founder genotypes")
    offspringGen <- offspringGen[,nomissing,drop=FALSE]
    founderGen <- founderGen[,nomissing,drop=FALSE]

    longbits <- .Machine$sizeof.long*8
    if(max.n.markers > longbits-1) {
        max.n.markers <- longbits-1
        warning("We can't use max.n.markers > ", longbits-1,
                ", so we're taking max.n.markers = ", longbits-1)
    }
    n.mar <- ncol(offspringGen)
    if(max.n.markers > n.mar) max.n.markers <- n.mar
    max.offset <- ceiling((max.n.markers-1)/2)

    n.ind <- nrow(offspringGen)
    n.founders <- nrow(founderGen)
    if(n.mar != ncol(founderGen))
        stop("ncol(offspringGen) != ncol(founderGen)")
    if(any(!is.na(offspringGen) & offspringGen != 0 & offspringGen != 1))
        stop("offspringGen should be NA, 0 or 1")
    if(any(!is.na(founderGen) & founderGen != 0 & founderGen != 1))
        stop("founderGen should be NA, 0 or 1")

    offspringGen[is.na(offspringGen)] <- -1

    z <- .C("R_inferFounderHap",
            as.integer(n.mar),
            as.integer(n.founders),
            as.integer(n.ind),
            as.integer(founderGen),
            as.integer(offspringGen),
            as.integer(max.offset),
            hap=as.integer(rep(0,n.mar * n.ind)),
            PACKAGE="qtl")
    z$hap[z$hap <= 0] <- NA

    fullhap <- matrix(ncol=length(nomissing), nrow=n.ind)
    fullhap[,nomissing] <- matrix(z$hap, ncol=n.mar, nrow=n.ind)
    colnames(fullhap) <- names(nomissing)

    fullhap
}


restoreMWrilGeno <-
    function(cross)
{
    g <- pull.geno(cross)
    f <- cross$founderGeno
    uf <- unique(as.numeric(f[!is.na(f)]))
    f[is.na(f)] <- missingval <- min(uf)-1
    g[is.na(g)] <- 0

    n.mar <- ncol(g)
    n.ind <- nrow(g)
    n.str <- nrow(f)
    if(n.mar != ncol(f))
        stop("no. genotypes inconsistent between offspring and founders.")

    crosses <- cross$cross
    if(ncol(crosses) != n.str || nrow(crosses) != n.ind)
        stop("Incompatiability in cross$cross dimension.")

    gen <-
        .C("R_restoreMWrilGeno",
           as.integer(n.ind),
           as.integer(n.mar),
           as.integer(n.str),
           as.integer(f),
           gen=as.integer(g),
           as.integer(crosses),
           as.integer(missingval),
           PACKAGE="qtl")$gen
    gen[gen==missingval] <- NA
    gen <- matrix(gen, nrow=n.ind, ncol=n.mar)
    colnames(gen) <- colnames(g)
    gen
}


# end of inferFounderHap.R
