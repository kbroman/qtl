#####################################################################
#
# scanoneboot.R
#
# copyright (c) 2007-2011, Karl W Broman
# last modified Mar, 2011
# first written Apr, 2007
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
# Contains: scanoneboot, summary.scanoneboot, print.scanoneboot
#
######################################################################

######################################################################
# scanoneboot: function to get bootstrap-based confidence interval for
#              QTL location
######################################################################
scanoneboot <-
    function(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"), upper=FALSE,
             ties.random=FALSE, start=NULL, maxit=4000, tol=1e-4,
             n.boot=1000, verbose=FALSE)
{
    if(!missing(chr)) cross <- subset(cross, chr)
    if(nchr(cross) != 1) { # scan just one chromosome
        warning("Considering just the first chromosome (", names(cross$geno)[1], ").")
        cross <- subset(cross, names(cross$geno)[1])
    }

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    if(length(pheno.col) > 1)
        stop("pheno.col should indicate a single phenotype")
    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }
    if(pheno.col < 1 || pheno.col > nphe(cross))
        stop("pheno.col should be between 1 and ", nphe(cross))

    # do scan with actual data
    out <- scanone(cross, pheno.col=pheno.col, model=model, method=method,
                   addcovar=addcovar, intcovar=intcovar, weights=weights,
                   use=use, upper=upper, ties.random=ties.random, start=start,
                   maxit=maxit, tol=tol)

    maxlod <- max(out[,3],na.rm=TRUE)
    w <- which(!is.na(out[,3]) & out[,3] == maxlod)

    results <- rep(NA, n.boot)
    n.ind <- nind(cross)
    n.prnt <- floor(n.boot/20)

    for(i in 1:n.boot) {
        temp <- subset(cross, ind=sample(n.ind, replace=TRUE))
        out <- scanone(temp, pheno.col=pheno.col, model=model, method=method,
                       addcovar=addcovar, intcovar=intcovar, weights=weights,
                       use=use, upper=upper, ties.random=ties.random, start=start,
                       maxit=maxit, tol=tol)
        mx <- max(out[,3],na.rm=TRUE)
        w <- out[!is.na(out[,3]) & out[,3]==mx,2]
        if(length(w) > 1) w <- sample(w,1)
        results[i] <- w

        if(verbose && ((i-1) %% n.prnt) == 0)
            cat("replicate", i, "\n")
    }

    attr(results, "results") <- out
    class(results) <- "scanoneboot"
    results
}


# summary function for scanoneboot output
summary.scanoneboot <-
    function(object, prob=0.95, expandtomarkers=FALSE, ...)
{
    lo <- (1-prob)/2
    results <- attr(object, "results")
    o <- max(results)
    qu <- quantile(object, c(lo, 1-lo))

    wh1 <- which(results[,2] <= qu[1])
    wh1 <- wh1[length(wh1)]
    wh2 <- which(results[,2] >= qu[2])
    wh2 <- wh2[1]

    if(expandtomarkers) {
        markerpos <- (1:nrow(results))[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", rownames(results))]
        if(any(markerpos <= wh1))
            wh1 <- max(markerpos[markerpos <= wh1])
        if(any(markerpos >= wh2))
            wh2 <- min(markerpos[markerpos >= wh2])
    }

    rbind(results[wh1,], o, results[wh2,])
}


# print function for scanoneboot output
print.scanoneboot <-
    function(x, ...)
{
    print(as.numeric(x))
}

# end of scanoneboot.R
