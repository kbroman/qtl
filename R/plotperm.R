#####################################################################
#
# plotperm.R
#
# copyright (c) 2007-2014, Karl W Broman
# last modified Oct, 2014
# first written Dec, 2007
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
# Contains: plot.scanoneperm, plot.scantwoperm, plot.scanoneboot
#
######################################################################

############################################################
# plot.scanoneboot
#
# plot histogram of the results from scanoneboot
############################################################
plot.scanoneboot <-
    function(x, ...)
{
    results <- attr(x, "results")
    markerpos <- results[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", rownames(results)),2]

    hideplot.scanoneboot <-
        function(x, breaks, xlim, main, xlab, ...)
        {
            if(missing(breaks)) breaks <- ceiling(2*sqrt(length(x)))
            if(missing(xlim)) xlim <- range(results[,2])
            if(missing(main)) main <- ""
            if(missing(xlab)) xlab <- "QTL position (cM)"

            hist(x, xlim=xlim, xlab=xlab, main=main, breaks=breaks, ...)
        }

    hideplot.scanoneboot(x, ...)
    rug(markerpos)
}

############################################################
# plot.scanoneperm
#
# plot histogram of the permutation results from scanone
############################################################
plot.scanoneperm <-
    function(x, lodcolumn=1, ...)
{
    # subroutine for hiding arguments in ...
    hideplot.scanoneperm <-
        function(A, X, breaks, xlab, xlim, main, ...)
        {
            if(missing(xlab)) xlab <- "maximum LOD score"

            if(missing(X)) {
                if(missing(breaks)) breaks <- ceiling(2*sqrt(length(A)))
                if(missing(xlim)) xlim <- c(0, max(A))
                if(missing(main)) main <- ""

                hist(A, breaks=breaks, xlab=xlab, xlim=xlim, main=main, ...)
            }
            else {
                mfrow <- par("mfrow")
                on.exit(par(mfrow=mfrow))
                par(mfrow=c(2,1))

                if(missing(breaks)) {
                    breaks.missing <- TRUE
                    breaks <- seq(0, max(c(A,X)), length=2*sqrt(length(A)))
                }
                else breaks.missing <- FALSE

                if(missing(xlim)) xlim <- c(0, max(c(A,X)))
                if(missing(main)) {
                    main <- "Autosomes"
                    main.missing <- TRUE
                }
                else {
                    main.missing <- FALSE
                    if(length(main)>1) { main2 <- main[2]; main <- main[1] }
                    else main2 <- main
                }

                hist(A, xlim=xlim, breaks=breaks, xlab=xlab, main=main, ...)
                rug(A)

                if(breaks.missing)
                    breaks <- seq(0, max(c(A,X)), length=2*sqrt(length(X)))
                if(main.missing) main <- "X chromosome"
                else main <- main2
                hist(X, xlim=xlim, breaks=breaks, xlab=xlab, main=main, ...)
                rug(X)
            }
        }

    # now to the actual code
    if(is.list(x)) { # separate X chr results
        if(lodcolumn < 1 || lodcolumn > ncol(x$A))
            stop("lodcolumn should be between 1 and ", ncol(x$A))
        A <- as.numeric(x$A[,lodcolumn])
        X <- as.numeric(x$X[,lodcolumn])
        hideplot.scanoneperm(A, X, ...)
    }
    else {
        if(lodcolumn < 1 || lodcolumn > ncol(x))
            stop("lodcolumn should be between 1 and ", ncol(x$A))
        A <- as.numeric(x[,lodcolumn])
        hideplot.scanoneperm(A, ...)
    }
}


############################################################
# plot.scantwoperm
#
# plot histogram of the permutation results from scantwo
############################################################
plot.scantwoperm <-
    function(x, lodcolumn=1, include_rug=TRUE, ...)
{
    hideplot.scantwoperm <-
        function(x, include_rug=TRUE, xlim, breaks, xlab, main, ...)
        {
            if(missing(xlim)) xlim <- c(0, max(unlist(x)))
            if(missing(xlab)) xlab <- "maximum LOD score"
            if(missing(main)) main.missing <- TRUE
            else { main.missing <- FALSE; main.input <- main }


            mfcol <- par("mfcol")
            on.exit(par(mfcol=mfcol))

            if("AA" %in% names(x)) {

                if(missing(breaks))
                    breaks <- seq(0, max(unlist(x)), len=ceiling(4*sqrt(length(x[[1]][[1]]))+1))

                par(mfrow=c(3,6))

                for(i in seq(along=x)) {
                    for(j in seq(along=x[[i]])) {
                        if(main.missing) main <- paste(names(x)[i], names(x[[i]])[j])
                        else {
                            if(length(main.input) >= i) main <- main.input[i]
                            else if(length(main.input) == 1) main <- main.input
                            else main <- ""
                        }

                        hist(x[[i]][[j]], xlim=xlim, breaks=breaks, xlab=xlab, main=main,...)
                        if(include_rug) rug(x[[i]][[j]])
                    }
                }
                invisible(return(NULL))
            }

            if(missing(breaks))
                breaks <- seq(0, max(unlist(x)), len=ceiling(4*sqrt(length(x[[1]]))+1))


            par(mfcol=c(3,2))

            for(i in seq(along=x)) {
                if(main.missing) main <- names(x)[i]
                else {
                    if(length(main.input) >= i) main <- main.input[i]
                    else if(length(main.input) == 1) main <- main.input
                    else main <- ""
                }

                hist(x[[i]], xlim=xlim, breaks=breaks, xlab=xlab, main=main,...)
                if(include_rug) rug(x[[i]])
            }
        }

    if(length(lodcolumn) > 1)
        stop("Select just one lod column")
    x <- x[,lodcolumn]
    hideplot.scantwoperm(x, include_rug=include_rug, ...)
}


# end of plotperm.R
