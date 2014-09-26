######################################################################
# plotModel.R
#
# copyright (c) 2008-9, Karl W Broman
# last modified January, 2009
# first written Apr, 2008
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
# Contains: plotModel
#
######################################################################

######################################################################
# plotModel
#
# plot a QTL model
######################################################################
plotModel <-
    function(qtl, formula, circrad.rel=0.25, circrad.abs,
             cex.name=1, chronly=FALSE, order, ...)
{
    if(missing(qtl))
        stop("Must provide qtl object or a vector of QTL names")

    if(is.character(qtl)) { # qtl names
        if(missing(formula)) formula <- NULL
    }
    else {
        if(missing(formula)) {
            if("formula" %in% names(attributes(qtl)))
                formula <- attr(qtl, "formula")
            else {
                if("formula" %in% names(qtl))
                    formula <- qtl$formula
                else formula <- NULL
            }
        }
        if(length(qtl)>0) {

            if(chronly) qtl <- qtl$chr
            else {
                if("name" %in% names(qtl))
                    qtl <- qtl$name
                else {
                    if("pos" %in% names(qtl))
                        qtl <- paste(qtl$chr, roundqtlpos(qtl$pos, 1), sep="@")
                    else qtl <- qtl$chr
                }
            }
        }
    }

    nqtl <- length(qtl)

    if(!is.null(formula)) {
        if(is.character(formula))
            formula <- as.formula(formula)

        theterms <- attr(terms(formula), "factors")[-1,]
        rn <- rownames(theterms)
        g <- grep("^[Qq][0-9]+$", rn)
        qtlnum <- as.numeric(substr(rn[g], 2, nchar(rn[g])))
        if(any(qtlnum < 1 | qtlnum > nqtl))
            stop("QTL in formula must be numbered between 1 and ", nqtl)

        cn <- colnames(theterms)
        intxn <- cn[grep("^[Qq][0-9]+:[Qq][0-9]+$", cn)]
    }
    else intxn <- NULL

    plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n",
         xaxs="i", yaxs="i", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),
         ...)

    if(!is.null(qtl) && length(qtl) > 0) {
        radloc <- rev(seq(pi/2, 2.5*pi, len=nqtl+1)[-1])
        xloc <- cos(radloc)
        yloc <- sin(radloc)

        rad <- seq(0, 2*pi, length=100)
        if(missing(circrad.abs)) {
            if(length(xloc) == 1)
                circrad.abs <- circrad.rel*2
            else
                circrad.abs <- circrad.rel * sqrt(diff(xloc[1:2])^2 + diff(yloc[1:2])^2)
            if(circrad.abs > 0.45) circrad.abs <- 0.45
        }
        # use the smallest of the two relative lengths
        u <- par("usr")
        pin <- par("pin")
        pin <- pin/c(diff(u[1:2]), diff(u[3:4]))
        circrad.abs<- circrad.abs/(pin/min(pin))

        if(!missing(order)) {
            if(length(order) != length(qtl))
                stop("order should have length ", nqtl)
            if(!all(sort(order) == seq(along=qtl)))
                stop("order should be a permutaton of 1,2,...,", nqtl)
            xloc <- xloc[order(order)]
            yloc <- yloc[order(order)]
        }

        if(length(intxn) > 0) {
            for(i in intxn) {
                thisint <- strsplit(i, ":")[[1]]
                theseqtl <- as.numeric(substr(thisint, 2, nchar(thisint)))
                lines(xloc[theseqtl], yloc[theseqtl], lwd=2)
            }
        }

        for(i in seq(along=xloc))
            polygon(xloc[i]+cos(rad)*circrad.abs[1], yloc[i]+sin(rad)*circrad.abs[2], col="white", border="black", lwd=2)

        text(xloc, yloc, qtl, cex=cex.name)
    }

    invisible()
}

roundqtlpos <-
    function (x, digits = 1)
{
    if (digits < 1)
        stop("This is intended for the case digits >= 1.")
    y <- as.character(round(x, digits))
    z <- strsplit(y, "\\.")
    sapply(z, function(a, digits) {
        if (length(a) == 1)
            b <- paste(a[1], ".", paste(rep("0", digits), collapse = ""),
                       sep = "")
        else {
            if (nchar(a[2]) == digits)
                b <- paste(a, collapse = ".")
            else b <- paste(a[1], ".", a[2], paste(rep("0", digits -
                                                       nchar(a[2])), collapse = "."), sep = "")
        }
    }, digits)
}

# end of plotModel.R
