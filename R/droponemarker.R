######################################################################
#
# droponemarker.R
#
# copyright (c) 2010, Karl W Broman
# last modified Oct, 2010
# first written Oct, 2001
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
# Contains: droponemarker
#
######################################################################

######################################################################
#
# droponemarker: Drop one marker at a time from a genetic map, to
#                evaluate the change in log likelihood and in map
#                length, in order to identify problematic markers
#
######################################################################

droponemarker <-
    function(cross, chr, error.prob=0.0001,
             map.function=c("haldane", "kosambi", "c-f", "morgan"),
             m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE, verbose=TRUE)
{
    if(!("cross" %in% class(cross)))
        stop("Input must have class \"cross\".")

    if(!missing(chr)) cross <- subset(cross, chr=chr)
    if(any(nmar(cross) < 3)) {
        if(all(nmar(cross) < 3))
            stop("No chromosomes with at least three markers\n")
        todrop <- names(cross$geno)[nmar(cross) < 3]
        tokeep <- names(cross$geno)[nmar(cross) > 2]
        warning("Dropping chr with <3 markers: ", paste(todrop, collapse=", "))
        cross <- subset(cross, chr=tokeep)
    }

    map.function <- match.arg(map.function)
    if(verbose) cat(" -Re-estimating map\n")
    origmap <- est.map(cross, error.prob=error.prob, map.function=map.function,
                       m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp)

    cross <- replace.map(cross, origmap)
    origmaptab <- pull.map(cross, as.table=TRUE)
    origmaptab <- cbind(origmaptab, LOD=rep(NA, nrow(origmaptab)))
    if(is.matrix(origmap[[1]])) {
        origmaptab <- cbind(origmaptab, Ldiff.female=rep(NA, nrow(origmaptab)),
                            Ldiff.male=rep(NA, nrow(origmaptab)))
        sexsp <- TRUE
    }
    else {
        origmaptab <- cbind(origmaptab, Ldiff=rep(NA, nrow(origmaptab)))
        sexsp <- FALSE
    }

    for(i in names(cross$geno)) {
        if(sexsp) {
            Lf <- diff(range(origmap[[i]][1,]))
            Lm <- diff(range(origmap[[i]][2,]))
        }
        else L <- diff(range(origmap[[i]]))

        if(verbose) cat(" -Chromosome", i, "\n")
        mnames <- markernames(cross, chr=i)
        temp <- subset(cross, chr=i)
        for(j in seq(along=mnames)){
            if(verbose > 1) cat(" ---Marker", j, "of", length(mnames), "\n")
            markerll <- markerloglik(cross, mnames[j], error.prob)

            newmap <- est.map(drop.markers(temp, mnames[j]), error.prob=error.prob,
                              map.function=map.function, m=m, p=p, maxit=maxit, tol=tol,
                              sex.sp=sex.sp)

            if(sexsp) {
                origmaptab[mnames[j],4] <- -(attr(origmap[[i]], "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10)
                origmaptab[mnames[j],5] <- Lf - diff(range(newmap[[1]][1,]))
                origmaptab[mnames[j],6] <- Lm - diff(range(newmap[[1]][2,]))
            }
            else {
                origmaptab[mnames[j],3] <- -(attr(origmap[[i]], "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10)
                origmaptab[mnames[j],4] <- L - diff(range(newmap[[1]]))
            }
        }
    }

    class(origmaptab) <- c("scanone", "data.frame")
    origmaptab$chr <- factor(origmaptab$chr, levels=unique(origmaptab$chr))
    origmaptab
}


# end of droponemarker.R
