######################################################################
#
# map_construction.R
#
# copyright (c) 2008-2013, Karl W Broman
# last modified Dec, 2013
# first written Oct, 2008
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
# Contains: formLinkageGroups, orderMarkers, orderMarkers.sub
#
######################################################################

######################################################################
# formLinkageGroups
#
# Use the estimated recombination fractions between pairs of markers
# (and LOD scores for a test of rf = 1/2) to partition the markers
# into a set of linkage groups.
#
# Two markers are placed in the same linkage group if rf <= max.rf
# and LOD >= min.lod.  The transitive property (if A is linked to B
# and B is linked to C then A is linked to C) is used to close the
# groups.
#
# If reorgMarkers=FALSE, the output is a data frame with two columns:
# the initial chromosome assignments and the linkage group assigments
# determined from the pairwise recombination fractions.
#
# If reorgMarkers=TRUE, the output is an experimental cross object,
# with the data reorganized according to the inferred linkage groups.
#
# The linkage groups are sorted by the number of markers they contain
# (from largest to smallest).
#
# If verbose=TRUE, tracing information is printed.
######################################################################
formLinkageGroups <-
    function(cross, max.rf=0.25, min.lod=3, reorgMarkers=FALSE,
             verbose=FALSE)
{
    if(!("rf" %in% names(cross))) {
        warning("Running est.rf.")
        cross <- est.rf(cross)
    }
    n.mar <- nmar(cross)
    tot.mar <- totmar(cross)
    rf <- cross$rf
    diagrf <- diag(rf)
    if(ncol(rf) != tot.mar)
        stop("dimension of recombination fractions inconsistent with no. markers in cross.")

    onlylod <- attr(cross$rf, "onlylod")
    if(!is.null(onlylod) && onlylod) { # results of markerlrt()
        if(!missing(max.rf))
            warning("max.rf ignored, as markerlrt() was used.")
        max.rf <- Inf
    }

    marnam <- colnames(rf)
    chrstart <- rep(names(cross$geno), n.mar)

    lod <- rf
    lod[lower.tri(rf)] <- t(rf)[lower.tri(rf)]
    rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
    diag(rf) <- 1
    diag(lod) <- 0

    ingrp <- 1:tot.mar
    for(i in 1:tot.mar) {
        if(verbose) {
            if(tot.mar > 100)
                if(i==round(i,-2)) cat(i,"of", tot.mar, "\n")
                else
                    if(i==round(i,-1)) cat(i,"of", tot.mar, "\n")
        }
        wh <- (rf[,i]<=max.rf & lod[,i] > min.lod)

        if(any(wh) && length(unique(ingrp[c(i, which(wh))]))>1) {
            oldgrp <- ingrp[wh]
            ingrp[wh] <- ingrp[i]
            u <- unique(oldgrp[oldgrp != ingrp[i]])
            ingrp[!is.na(match(ingrp, u))] <- ingrp[i]
        }
    }
    tab <- sort(table(ingrp), decreasing=TRUE)
    u <- names(tab)

    revgrp <- ingrp
    for(i in seq(along=u))
        revgrp[ingrp==u[i]] <- i

    if(reorgMarkers) {
        cross <- clean(cross)
        chrtype <- rep(sapply(cross$geno, class), n.mar)
        crosstype <- class(cross)[1]
        g <- pull.geno(cross)

        cross$geno <- vector("list", max(revgrp))
        names(cross$geno) <- 1:max(revgrp)

        for(i in 1:max(revgrp)) {
            cross$geno[[i]]$data <- g[,revgrp==i,drop=FALSE]

            cross$geno[[i]]$map <- seq(0, by=10, length=tab[i])
            if(crosstype=="4way") {
                cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
                                             cross$geno[[i]]$map)
                colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
            }
            else
                names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)

            thechrtype <- unique(chrtype[revgrp==i])
            if(length(thechrtype) > 1)
                warning("Problem with linkage group ", i, ": A or X?\n", paste(thechrtype, collapse=" "))
            else
                class(cross$geno[[i]]) <- thechrtype

        }

        mname <- markernames(cross)
        m <- match(mname, marnam)
        rf <- rf[m,m]
        lod <- lod[m,m]
        rf[upper.tri(rf)] <- lod[upper.tri(lod)]
        diag(rf) <- diagrf[m]
        cross$rf <- rf

        return(cross)
    }
    else {
        result <- data.frame(origchr=factor(chrstart, levels=names(cross$geno)),
                             LG=factor(revgrp, levels=1:max(revgrp)), stringsAsFactors=TRUE)
        rownames(result) <- marnam
        return(result)
    }
}

######################################################################
# orderMarkers
#
# For each of the selected chromosomes, construct a new genetic map
# from scratch.  Marker order is determined in by an expedient and not
# necessarily good algorithm, with orders compared by the number of
# obligate crossovers.
######################################################################
orderMarkers <-
    function(cross, chr, window=7, use.ripple=TRUE, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             maxit=4000, tol=1e-4, sex.sp=TRUE, verbose=FALSE)
{
    map.function <- match.arg(map.function)
    if(!missing(chr)) chr <- matchchr(chr, names(cross$geno))
    else chr <- names(cross$geno)

    n.mar <- nmar(cross)

    if(verbose > 1) verbose.sub <- TRUE else verbose.sub <- FALSE

    for(i in chr) {
        if(verbose && length(chr) > 1) cat(" - Chr", i,"\n")
        if(n.mar[i] > 2) {
            neworder <- orderMarkers.sub(cross, i, window=window, use.ripple=use.ripple,
                                         verbose=verbose.sub)
            cross <- switch.order(cross, i, neworder, error.prob=error.prob,
                                  map.function=map.function, maxit=maxit, tol=tol,
                                  sex.sp=sex.sp)
        }
    }

    cross
}

######################################################################
# orderMarkers.sub
# For the markers on a chromosome, use a greedy algorithm to order
# the markers de novo, possibly running ripple() to refine the order.
######################################################################
orderMarkers.sub <-
    function(cross, chr, window=7, use.ripple=TRUE, verbose=FALSE)
{
    if(missing(chr)) chr <- names(cross$geno)[1]
    if(length(chr) > 1) {
        if(length(grep("^-", chr)) > 0)
            stop("Need to give a single chromosome name.")
        warning("Need to give a single chromosome name; using just the first")
        chr <- chr[1]
    }

    if(length(matchchr(chr, names(cross$geno)))>1)
        stop("Chr ", chr, " not found.")

    cross <- subset(cross, chr=chr)
    names(cross$geno)[1] <- "1"

    n.mar <- totmar(cross)
    if(n.mar < 3) return(1:n.mar)

    if(use.ripple && n.mar <= window) { # just use ripple
        rip <- summary(ripple(cross, chr=1, window=window, verbose=FALSE))

        nxo <- rip[1:2,ncol(rip)]
        if(nxo[1] <= nxo[2]) return(1:n.mar)
        else return(rip[2,1:(ncol(rip)-1)])
    }

    nt <- ntyped(cross, "mar")

    # start with the most typed markers and move down
    themar <- order(nt, decreasing=TRUE)
    marnam <- markernames(cross)
    if(n.mar > 3) {
        for(i in 4:n.mar)
            cross <- movemarker(cross, marnam[themar[i]], 2)
    }

    # create matrix of orders to test
    makeorders <-
        function(n)
        {
            orders <- matrix(ncol=n, nrow=n)
            for(k in 1:n) { orders[k,n-k+1] <- n; orders[k,-(n-k+1)] <- 1:(n-1) }
            orders
        }

    # simple switch of marker order on chr 1
    simpleswitch <-
        function(cross, neworder)
        {
            cross$geno[[1]]$data <- cross$geno[[1]]$data[,neworder]
            if(is.matrix(cross$geno[[1]]$map))
                cross$geno[[1]]$map <- cross$geno[[1]]$map[,neworder]
            else
                cross$geno[[1]]$map <- cross$geno[[1]]$map[neworder]
            cross
        }

    # work on marker 3
    if(verbose) cat(" --- Adding marker 3 of", n.mar, "\n")
    orders <- makeorders(3)
    nxo <- rep(NA, nrow(orders))
    nxo[1] <- sum(countXO(cross, 1))
    temp <- cross
    for(kk in 2:nrow(orders)) {
        temp$geno[[1]]$data <- temp$geno[[1]]$data[,orders[kk,]]
        nxo[kk] <- sum(countXO(cross, 1))
    }
    wh <- which(nxo==min(nxo))
    if(length(wh) > 1) wh <- sample(wh, 1)
    if(wh > 1)
        cross <- simpleswitch(cross, orders[wh,])

    # rest of the markers
    if(n.mar > 3) {
        for(k in 4:n.mar) {
            if(verbose) cat(" --- Adding marker", k, "of", n.mar, "\n")
            cross <- movemarker(cross, marnam[themar[k]], 1)
            orders <- makeorders(k)
            nxo <- rep(NA, nrow(orders))
            nxo[1] <- sum(countXO(cross, 1))
            temp <- cross
            for(kk in 2:nrow(orders)) {
                temp$geno[[1]]$data <- cross$geno[[1]]$data[,orders[kk,]]
                nxo[kk] <- sum(countXO(temp, 1))
            }

            wh <- which(nxo==min(nxo))
            if(length(wh) > 1) wh <- sample(wh, 1)
            if(wh > 1)
                cross <- simpleswitch(cross, orders[wh,])
        }
    }

    if(use.ripple) {
        dif <- -8
        while(dif < 0) {
            rip <- summary(ripple(cross, chr=1, window=window, verbose=FALSE))
            dif <- diff(rip[1:2, ncol(rip)])
            if(dif < 0)
                cross <- simpleswitch(cross, rip[2,1:(ncol(rip)-1)])
        }
    }

    match(colnames(cross$geno[[1]]$data), marnam)
}


# end of map_construction.R
