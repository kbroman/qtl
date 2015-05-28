######################################################################
#
# add_threshold.R
#
# copyright (c) 2006-2014, Karl W Broman
# last modified Apr, 2014
# first written Dec, 2006
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
# Contains: add.threshold, xaxisloc.scanone
#
######################################################################

######################################################################
# add.threshold: function to add threshold lines to a plot
# created by plot.scanone()
#
# out: scanone output used to create the plot
# chr: chromosomes that were plotted
#
# perms: scanone permutation results
# alpha: significance level (a single number)
#
# gap:   the gap between chromosomes, as specified in the call to plot.scanone
#
# ...:   extra arguments passed to abline or segments to control line types/widths/colors
#        (e.g., specify lty=2 to give dashed lines)
######################################################################

add.threshold <-
    function(out, chr, perms, alpha=0.05, lodcolumn=1, gap=25, ...)
{
    if(missing(out))
        stop("You must provide scanone output, so we can get chromosome lengths.")
    if(!any(class(out) == "scanone"))
        stop("out should have class \"scanone\".")
    if(!missing(chr))
        out <- subset(out, chr=chr)

    if(missing(perms))
        stop("You must specify permutation results, to get the thresholds")
    if(length(alpha) > 1 || alpha<0 || alpha>1)
        stop("alpha should have length 1 and be between 0 and 1.")

    thr <- summary(perms, alpha=alpha)


    if(!is.list(thr)) {
        if(any(lodcolumn < 1 | lodcolumn > length(thr)))
            stop("lodcolumn should be between 1 and ", length(thr))
        abline(h=thr[lodcolumn], ...)
    }
    else {
        if(any(lodcolumn < 1 | lodcolumn > length(thr$A)))
            stop("lodcolumn should be between 1 and ", length(thr$A))
        a <- thr$A[lodcolumn]
        x <- thr$X[lodcolumn]
        noX <- FALSE
        xchr <- attr(perms, "xchr")

        L <- tapply(out[,2], out[,1], function(a) diff(range(a)))
        L <- L[!is.na(L)]
        xchr <- xchr[match(names(L), names(xchr))]

        if(all(xchr))
            abline(h=x, ...)
        else if(all(!xchr))
            abline(h=a, ...)
        else {

            start <- c(0,cumsum(L+gap))
            end <- start+ c(L,0)
            wh <- which(!xchr)
            if(length(wh)==1 || all(diff(wh)==1))
                segments(start[min(wh)], a, end[max(wh)], a, ...)
            else
                segments(start[wh], a, start[wh+1], a, ...)
            wh <- which(xchr)
            if(length(wh)==1 || all(diff(wh)==1))
                segments(start[min(wh)], x, end[max(wh)], x, ...)
            else
                segments(start[wh], x, start[wh+1], x, ...)
        }
    }
    invisible()
}

# xaxisloc.scanone
# find x-axis locations for a plot of scanone output
xaxisloc.scanone <-
    function(out, thechr, thepos, chr, gap=25)
{
    if(missing(out))
        stop("You must provide scanone output, so we can get chromosome lengths.")
    if(!any(class(out) == "scanone"))
        stop("out should have class \"scanone\".")
    if(!missing(chr))
        out <- subset(out, chr=chr)
    chr <- unique(out[,1])

    if(length(thechr) != 1) {
        if(length(thepos)==1)
            thepos <- rep(thepos, length(thechr))
        else if(length(thechr) != length(thepos))
            stop("If thechr and thepos have length>1 they must both have the same length")
    }
    else {
        if(length(thepos)!=1)
            thechr <- rep(thechr, length(thepos))
    }

    if(length(thechr) > 1) {
        res <- rep(NA, length(thechr))
        for(i in seq(along=thechr))
            res[i] <- xaxisloc.scanone(out, thechr[i], thepos[i], chr, gap)
        return(res)
    }

    L <- tapply(out[,2], out[,1], function(a) diff(range(a, na.rm=TRUE)))
    Lmin <- tapply(out[,2], out[,1], min, na.rm=TRUE)
    start <- c(0,cumsum(L+gap))

    start[which(as.character(thechr)==chr)]+(thepos - Lmin[as.character(thechr)==chr])
}

# end of add_threshold.R
