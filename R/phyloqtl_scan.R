######################################################################
# phyloqtl_scan.R
#
# copyright (c) 2009-2012, Karl W Broman
# last modified Oct, 2012
# first written May, 2009
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
# Single-QTL scan to map QTL to a phylogenetic tree
#
# Part of the R/qtl package
# Contains: scanPhyloQTL, max.scanPhyloQTL, plot.scanPhyloQTL,
#           inferredpartitions
#
######################################################################

scanPhyloQTL <-
    function(crosses, partitions, chr, pheno.col=1,
             model=c("normal","binary"), method=c("em","imp","hk"),
             addcovar, maxit=4000, tol=0.0001, useAllCrosses=TRUE,
             verbose=FALSE)
{
    if(missing(chr)) chr <- names(crosses[[1]]$geno)
    model <- match.arg(model)
    method <- match.arg(method)

    thecrosses <- names(crosses)
    taxa <- sort(unique(unlist(strsplit(thecrosses, ""))))

    thecrosses <- names(crosses) <- checkPhyloCrosses(thecrosses, taxa)

    if(missing(partitions)) {
        # generate all partitions
        partitions <- genAllPartitions(length(taxa), taxa)
    }
    checkPhyloPartition(partitions, taxa)

    crossmat <- qtlByPartition(thecrosses, partitions)

    dimnames(crossmat) <- list(thecrosses, partitions)

    if(!missing(addcovar)) {
        if(!is.list(addcovar) || length(addcovar) != length(crosses))
            stop("addcovar must be a list with the same length as crosses (", length(crosses), ")")
        n.ind <- sapply(crosses, nind)
        if(is.matrix(addcovar[[1]])) {
            nind.addcovar <- sapply(addcovar, nrow)
            n.addcovar <- sapply(addcovar, ncol)
        }
        else {
            nind.addcovar <- sapply(addcovar, length)
            n.addcovar <- 1
        }
        if(any(nind.addcovar != n.ind)) {
            err <- paste("crosses: ", paste(n.ind, collapse=" "), "\n",
                         "addcovar:", paste(n.addcovar, collapse=" "), "\n")

            stop("Mismatch between no. individuals in addcovar and crosses.\n", err)
        }
        if(length(unique(n.addcovar)) > 1)
            stop("Mismatch in no. add've covariates: ", paste(n.addcovar, collapse=" "))
    }

    # check that the marker maps are all exactly the same
    n.chr <- sapply(crosses, nchr)
    if(!all(n.chr == n.chr[1]))
        stop("Different numbers of chromosomes")
    chrnam1 <- chrnames(crosses[[1]])
    for(j in 2:length(crosses)) {
        chrnam2 <- chrnames(crosses[[j]])
        if(!all(chrnam1 == chrnam2))
            stop("Different chromosome names")
    }
    n.mar1 <- nmar(crosses[[1]])
    for(j in 2:length(crosses)) {
        n.mar2 <- nmar(crosses[[j]])
        if(!all(n.mar1 == n.mar2))
            stop("Different numbers of markers")
    }
    mn1 <- markernames(crosses[[1]])
    for(j in 2:length(crosses)) {
        mn2 <- markernames(crosses[[j]])
        if(!all(mn1 == mn2))
            stop("Different marker names")
    }
    mp1 <- unlist(pull.map(crosses[[1]]))
    for(j in 2:length(crosses)) {
        mp2 <- unlist(pull.map(crosses[[j]]))
        if(!all(mp1 == mp2))
            stop("Different marker positions")
    }

    out <- vector("list", length(partitions))

    for(i in seq(along=partitions)) {
        if(verbose) cat("Partition", i, "of", length(partitions), "\n")
        cm <- crossmat[,i]

        # if all crosses need to be flipped, don't flip any of them
        if(!any(cm > 0)) cm <- -cm

        if(!useAllCrosses) {
            whcm <- which(cm != 0)
            x <- crosses[whcm]
            cm <- cm[whcm]
        }
        else {
            o <- order(abs(cm), decreasing=TRUE)
            cm <- cm[o]
            whcm <- seq(along=cm)
            x <- crosses[o]
            if(!missing(addcovar))
                addcovar <- addcovar[o]
        }

        # flip crosses if necessary
        if(any(cm < 0))
            for(j in which(cm < 0))
                x[[j]] <- flipcross(x[[j]])

        # combine the crosses
        xx <- x[[1]]
        n.phe <- nphe(x[[1]])
        if(length(x) > 1) {
            for(j in 2:length(x)) {
                xx <- c(xx, x[[j]])
                xx$pheno <- xx$pheno[,1:n.phe,drop=FALSE]
            }
        }

        # create cross indicators (as additive covariates)
        if(length(x)==1) alladdcovar<-NULL
        else {
            ni <- sapply(x, nind)
            alladdcovar <- matrix(0, ncol=length(x)-1, nrow=sum(ni))
            thestart <- cumsum(c(1,ni))
            end <- cumsum(ni)
            for(j in 2:length(x))
                alladdcovar[thestart[j]:end[j],j-1] <- 1
        }
        if(!missing(addcovar)) {
            theaddcovar <- as.matrix(addcovar[[whcm[1]]])
            if(length(whcm) > 1)
                for(j in 2:length(whcm))
                    theaddcovar <- rbind(theaddcovar, as.matrix(addcovar[[whcm[j]]]))

            # quick check to be sure that it's not a column with all one value
            if(ncol(theaddcovar)>1 || length(unique(theaddcovar))>1)
                alladdcovar <- cbind(alladdcovar, theaddcovar)
        }

        # ind with no QTL effect
        ind.noqtl <- rep(cm == 0, sapply(x, nind))

        # do the scan
        out[[i]] <- scanone(xx, chr=chr, pheno.col=pheno.col, addcovar=alladdcovar,
                            model=model, method=method, maxit=maxit, tol=tol,
                            ind.noqtl=ind.noqtl)
    }

    # just one partition
    if(length(out) == 1) return(out[[1]])

    # multiple partitions
    result <- out[[1]]
    for(j in 2:length(out))
        result <- c(result, out[[j]])
    colnames(result)[-(1:2)] <- partitions
    class(result) <- c("scanPhyloQTL", "scanone", "data.frame")
    result
}


max.scanPhyloQTL <-
    function(object, chr, format=c("postprob", "lod"), ...)
{
    format <- match.arg(format)

    if(!missing(chr))
        object <- subset.scanone(object, chr=chr)

    mx <- summary.scanPhyloQTL(object, format=format)

    if(nrow(mx) > 1) {
        nc <- ncol(mx)
        mx <- mx[!is.na(mx[,nc]) & mx[,nc]==max(mx[,nc], na.rm=TRUE),,drop=FALSE]
    }

    class(mx) <- c("summary.scanPhyloQTL", "summary.scanone", "data.frame")

    mx
}

summary.scanPhyloQTL <-
    function(object, format=c("postprob", "lod"), threshold, ...)
{
    format <- match.arg(format)

    themax <- apply(object[,-(1:2)], 2, tapply, object[,1], max, na.rm=TRUE)
    if(length(unique(object[,1]))==1) {
        themax <- as.data.frame(matrix(themax, nrow=1), stringsAsFactors=TRUE)
        names(themax) <- colnames(object)[-(1:2)]
    }
    wh <- apply(themax, 1, function(a) { a <- which(a==max(a)); if(length(a) > 1) a <- sample(a, 1); a })
    whpos <- rep(NA, length(wh))
    names(whpos) <- unique(object[,1])
    for(i in seq(along=whpos)) {
        temp <- object[object[,1]==names(whpos)[i],,drop=FALSE]
        temp2 <- which(temp[,wh[i]+2]==max(temp[,wh[i]+2],na.rm=TRUE))
        if(length(temp2) > 1) temp2 <- sample(temp2, 1)
        whpos[i] <- temp[temp2,2]
        names(whpos)[i] <- rownames(temp)[temp2]
    }
    if(format=="lod") {
        out <- data.frame(chr=unique(object[,1]), pos=whpos, themax,
                          loddif=apply(themax, 1, function(a) -diff(sort(a, decreasing=TRUE)[1:2])),
                          inferred=colnames(object)[wh+2],
                          maxlod=apply(themax, 1, max),
                          stringsAsFactors=TRUE)
    }
    else {
        out <- data.frame(chr=unique(object[,1]), pos=whpos, themax,
                          inferred=colnames(object)[wh+2],
                          maxlod=apply(themax, 1, max),
                          stringsAsFactors=TRUE)
        temp <- out[,-c(1:2, ncol(out)-0:1)]
        out[,-c(1:2, ncol(out)-0:1)] <- t(apply(temp, 1, function(a) 10^a/sum(10^a)))
    }
    colnames(out)[(1:ncol(themax))+2] <- colnames(themax)

    rownames(out) <- names(whpos)

    if(!missing(threshold))
        out <- out[out$maxlod >= threshold,,drop=FALSE]

    class(out) <- c("summary.scanPhyloQTL", "summary.scanone", "data.frame")
    out
}


# plot results of scanPhyloQTL
plot.scanPhyloQTL <-
    function(x, chr, incl.markers=TRUE, col, xlim, ylim, lwd=2,
             gap=25, mtick=c("line", "triangle"),
             show.marker.names=FALSE, alternate.chrid=FALSE,
             legend=TRUE, ...)
{
    mtick <- match.arg(mtick)

    if(!missing(chr)) x <- subset(x, chr=chr)

    if(missing(col)) {
        col <- c("black","blue","red","green","orange","brown","gray","cyan","magenta")
        if(ncol(x)-2 > length(col))
            stop("Please give a list of colors")
        col <- col[1:(ncol(x)-2)]
    }

    if(missing(ylim))
        ylim <- c(0, max(apply(x[,-(1:2)], 2, max)))

    dots <- list(...)

    if(missing(xlim)) {
        if("ylab" %in% names(dots))
            plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                         ylim=ylim, lwd=lwd, gap=gap, mtick=mtick,
                         show.marker.names=show.marker.names, alternate.chrid=alternate.chrid, ...)
        else
            plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                         ylim=ylim, lwd=lwd, gap=gap, mtick=mtick,
                         show.marker.names=show.marker.names, alternate.chrid=alternate.chrid,
                         ylab="LOD score", ...)
    }
    else {
        if("ylab" %in% names(dots))
            plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                         ylim=ylim, xlim=xlim, lwd=lwd, gap=gap, mtick=mtick,
                         show.marker.names=show.marker.names, alternate.chrid=alternate.chrid, ...)
        else
            plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                         ylim=ylim, xlim=xlim, lwd=lwd, gap=gap, mtick=mtick,
                         show.marker.names=show.marker.names, alternate.chrid=alternate.chrid,
                         ylab="LOD score", ...)
    }

    if(ncol(x) > 3)
        for(i in 2:(ncol(x)-2))
            plot.scanone(x, col=col[i], lodcolumn=i, add=TRUE, ...)

    if(is.character(legend) || legend) {
        if(is.character(legend))
            legend(legend, legend=colnames(x)[-(1:2)], col=col, lwd=lwd)
        else
            legend("topright", legend=colnames(x)[-(1:2)], col=col, lwd=lwd)
    }

    invisible()
}

inferredpartitions <-
    function(output, chr, lodthreshold, probthreshold=0.9)
{
    if(!("scanPhyloQTL" %in% class(output)))
        stop("Argument 'output' must be of class \"scanPhyloQTL\", as output by scanPhyloQTL.")

    if(missing(chr)) {
        chr <- output[1,1]
        warning("Missing chromosome; using ", chr)
    }
    else if(!any(output[,1]==chr))
        stop("Chromosome \"", chr, "\" not found.")

    if(missing(lodthreshold)) {
        warning("No lodthreshold given; using lodthreshold=0.")
        lodthreshold=0
    }

    output <- output[output[,1]==chr,]
    output[,1] <- as.factor(as.character(output[,1]))
    output <- summary(output, format="postprob")

    if(output$maxlod < lodthreshold) return("null")

    prob <- sort(output[,3:(ncol(output)-2)], decreasing=TRUE)
    cs <- cumsum(as.numeric(prob))
    wh <- min(which(cs >= probthreshold))
    names(prob)[1:wh]
}

# end of phyloqtl_scan.R
