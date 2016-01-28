######################################################################
#
# plot.scantwo.R
#
# copyright (c) 2001-2016, Karl W Broman, Hao Wu and Brian Yandell
# last modified Jan, 2016
# first written Nov, 2001
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
# Hao Wu (The Jackson Lab) wrote the initial code
#
# Part of the R/qtl package
# Contains: plot.scantwo
#
######################################################################

plot.scantwo <-
    function(x, chr, incl.markers = FALSE, zlim, lodcolumn=1,
             lower = c("full", "add", "cond-int", "cond-add", "int"),
             upper = c("int", "cond-add", "cond-int", "add", "full"),
             nodiag = TRUE,
             contours = FALSE, main, zscale = TRUE, point.at.max=FALSE,
             col.scheme = c("viridis", "redblue","cm","gray","heat","terrain","topo"),
             gamma=0.6, allow.neg=FALSE, alternate.chrid=FALSE, ...)
{
    if(!any(class(x) == "scantwo"))
        stop("Input should have class \"scantwo\".")

    col.scheme <- match.arg(col.scheme)

    if(length(dim(x$lod)) > 2) { # results from multiple phenotypes
        if(length(lodcolumn) > 1) {
            warning("Argument lodcolumn should be of length 1.")
            lodcolumn <- lodcolumn[1]
        }

        if(lodcolumn < 0 || lodcolumn > dim(x$lod)[3])
            stop("Argument lodcolumn misspecified.")
        x$lod <- x$lod[,,lodcolumn]
    }

    if(!missing(chr)) x <- subset(x, chr=chr)
    if(nrow(x$lod)==0) {
        warning("Empty scantwo object.")

        return(invisible(NULL))
    }


    chr <- as.character(unique(x$map[,1]))

    addpair <- attr(x, "addpair")
    if(!is.null(addpair) && addpair) {
        lower <- "full"
        upper <- "add"
        if(missing(zlim)) {
            if(allow.neg)
                zlim <- rep(max(abs(x$lod), na.rm=TRUE), 2)
            else
                zlim <- rep(max(x$lod, na.rm=TRUE), 2)
        }
        addpair <- TRUE
    }
    else addpair <- FALSE

    if(length(lower)==1 && lower == "fv1") lower <- "cond-int"
    if(length(lower)==1 && lower == "av1") lower <- "cond-add"
    if(length(upper)==1 && upper == "fv1") upper <- "cond-int"
    if(length(upper)==1 && upper == "av1") upper <- "cond-add"

    lower <- match.arg(lower)
    upper <- match.arg(upper)
    if(!any(class(x) == "scantwo"))
        stop("Input variable is not an object of class scantwo!")
    lod <- x$lod
    map <- x$map

    # backward compatibility for previous version of R/qtl
    if(!("scanoneX" %in% names(x))) {
        warning("It would be best to re-run scantwo() with the R/qtl version 0.98 or later.\n")
        scanoneX <- NULL
    }
    else scanoneX <- x$scanoneX

    # if incl.markers is FALSE, drop positions
    #     for which third column of map is 0
    if(!incl.markers && any(map[, 3] == 0)) {
        o <- (map[, 3] == 1)
        lod <- lod[o, o]
        map <- map[o, ]
        if(!is.null(scanoneX)) scanoneX <- scanoneX[o]
    }

    if(all(diag(lod) < 1e-14) && (lower == "cond-int" || lower=="cond-add") )
        stop("Need to run scantwo with run.scanone=TRUE.")

    oldlod <- lod

    lo <- lower.tri(lod)
    up <- upper.tri(lod)

    # grab the interaction LOD scores
    if(upper=="int")
        lod[up] <- t(oldlod)[up] - oldlod[up]
    if(lower=="int")
        lod[lo] <- oldlod[lo] - t(oldlod)[lo]

    if(lower=="add")
        lod[lo] <- t(oldlod)[lo]

    if(upper=="full")
        lod[up] <- t(oldlod)[up]

    # get conditional LOD scores
    if(lower=="cond-int" || lower=="cond-add") {
        if(lower=="cond-add")
            lod[lo] <- t(oldlod)[lo]

        thechr <- map$chr
        uchr <- unique(thechr)
        thechr <- factor(as.character(thechr), levels=as.character(uchr))
        uchr <- factor(as.character(uchr), levels=levels(thechr))
        xchr <- tapply(map$xchr, thechr, function(a) a[1])

        maxo <- tapply(diag(lod), thechr, max, na.rm=TRUE)
        if(any(xchr) && !is.null(scanoneX)) {
            maxox <- tapply(scanoneX, thechr, max, na.rm=TRUE)
            maxo[xchr] <- maxox[xchr]
        }

        n.chr <- length(chr)
        for(i in 1:n.chr) {
            pi <- which(thechr==uchr[i])
            for(j in i:n.chr) {
                pj <- which(thechr==uchr[j])
                temp <- lod[pj,pi] - max(maxo[c(i,j)])
                temp[!is.na(temp) & temp<0] <- 0
                if(i==j) lod[pj,pi][lower.tri(temp)] <- temp[lower.tri(temp)]
                else lod[pj,pi] <- temp
            }
        }
    }

    if(upper=="cond-int" || upper=="cond-add") {
        if(upper=="cond-int")
            lod[up] <- t(oldlod)[up]

        thechr <- map$chr
        uchr <- unique(thechr)
        thechr <- factor(as.character(thechr), levels=as.character(uchr))
        uchr <- factor(as.character(uchr), levels=levels(thechr))
        xchr <- tapply(map$xchr, thechr, function(a) a[1])

        maxo <- tapply(diag(lod), thechr, max, na.rm=TRUE)
        if(any(xchr) && !is.null(scanoneX)) {
            maxox <- tapply(scanoneX, thechr, max, na.rm=TRUE)
            maxo[xchr] <- maxox[xchr]
        }

        n.chr <- length(chr)
        for(i in 1:n.chr) {
            pi <- which(thechr==uchr[i])
            for(j in i:n.chr) {
                pj <- which(thechr==uchr[j])
                temp <- lod[pi,pj] - max(maxo[c(i,j)])
                temp[!is.na(temp) & temp<0] <- 0
                if(i==j) lod[pi,pj][upper.tri(temp)] <- temp[upper.tri(temp)]
                else lod[pi,pj] <- temp
            }
        }
    }


    if(nodiag) diag(lod) <- 0

    # deal with bad LOD score values
    if(any(is.na(lod))) {
        u <- is.na(lod)
        n <- sum(u)
        warning(n, " LOD scores NA, set to 0")
        lod[u] <- 0
    }
    if(!allow.neg && any(!is.na(lod) & lod < -1e-6)) {
        u <- !is.na(lod) & lod < 0
        n <- sum(u)
        warning(n, " LOD scores <0, set to 0")
        lod[u] <- 0
    }
    if(any(!is.na(lod) & lod == Inf)) {
        u <- !is.na(lod) & lod == Inf
        n <- sum(u)
        warning(n, " LOD scores =Inf, set to maximum observed value")
        lod[u] <- max(lod[!is.na(lod) & lod < Inf])
    }

    if(missing(zlim)) { # no given zlim
        # calculate the zlim for interactive and full LODs
        if(allow.neg) {
            zlim.int <- max(abs(lod[row(lod) < col(lod)]))
            zlim.jnt <- max(abs(lod[row(lod) >= col(lod)]))
        }
        else {
            zlim.int <- max(lod[row(lod) < col(lod)])
            zlim.jnt <- max(lod[row(lod) >= col(lod)])
        }
    }
    else {
        zlim.jnt <- zlim[1]
        if(length(zlim) < 2)
            zlim.int <- zlim[1]
        else
            zlim.int <- zlim[2]
    }

    # rescale the data in upper triangle based on zlims.jnt
    lod[row(lod) < col(lod)] <- lod[row(lod) < col(lod)] * zlim.jnt/zlim.int

    # make sure LOD values are below (0,zlim.jnt) or update zlim.jnt
    #  if(max(lod) > zlim.jnt) {
    #    warning("LOD values out of range; updating zlim.")
    #    temp <- max(lod)
    #    zlim.int <- zlim.int * temp/zlim.jnt
    #    zlim.jnt <- temp
    #  }

    # save old par parameters, to restore them on exit
    if(zscale) {
        old.mar <- par("mar")
        old.mfrow <- par("mfrow")
        old.las <- par("las")
        on.exit(par(las = old.las, mar = old.mar, mfrow = old.mfrow))
    }
    else {
        old.las <- par("las")
        on.exit(par(las = old.las))
    }
    par(las = 1)
    dots <- list(...)
    if(zscale) {
        if("layout" %in% names(dots))
            layout(dots[["layout"]][[1]],dots[["layout"]][[2]])
        else
            layout(cbind(1, 2), c(6, 1))

        if("mar1" %in% names(dots))
            par(mar=dots[["mar1"]])
        else
            par(mar = c(5, 4, 4, 2) + 0.1)
    }
    if( gamma < 0 && col.scheme == "redblue")
        stop( "gamma must be non-negative" )
    cols <- switch(col.scheme,
                   gray = if( gamma <= 0) rev(gray(seq(0,1,len=256)))
                   else rev(gray(log(seq(1,exp(gamma),len=256))/gamma)),
                   heat = heat.colors(256),
                   terrain = terrain.colors(256),
                   topo = topo.colors(256),
                   cm = cm.colors(256),
                   redblue = rev(rainbow(256, start = 0, end = 2/3)),
                   viridis = viridis_qtl(256)
                   )
    if(col.scheme=="redblue") {
        # convert colors using gamma=0.6 (which will no longer be available in R)
        rgbval <- (col2rgb(cols)/255)^0.6
        cols <- rgb(rgbval[1,], rgbval[2,], rgbval[3,])
    }

    if(allow.neg) {
        lo <- -zlim.jnt
        lo.int <- -zlim.int
    }
    else lo.int <- lo <- 0

    if("xlab" %in% names(dots)) {
        xlab <- dots$xlab
        if("ylab" %in% names(dots))
            ylab <- dots$ylab
        else ylab <- xlab
    }
    else {
        if("ylab" %in% names(dots))
            xlab <- ylab <- dots$ylab
        else {
            if(length(chr) > 1)
                xlab <- ylab <- "Chromosome"
            else
                xlab <- ylab <- "Location (cM)"
        }
    }

    if(length(chr) > 1)
        image(1:ncol(lod), 1:nrow(lod), lod, ylab = ylab,
              xlab = xlab, zlim = c(lo, zlim.jnt), col = cols,
              xaxt = "n", yaxt = "n")
    else
        image(map[,2], map[,2], lod, ylab = ylab,
              xlab = xlab, zlim = c(lo, zlim.jnt), col = cols)

    # plot point at maximum, if requested
    if(point.at.max) {
        temp <- lod
        temp[upper.tri(temp)] <- -50
        temp[diag(temp)] <- -50
        wh <- which(temp == max(temp), arr.ind=TRUE)
        if(length(chr) > 1)
            points(wh,rev(wh),pch=4,lwd=2)
        else {
            points(map[wh[,1],2],map[wh[,2],2],pch=4,lwd=2,col="blue")
            points(map[wh[,2],2],map[wh[,1],2],pch=4,lwd=2,col="blue")
        }
    }

    # add contours if requested
    if(any(contours > 0)) {
        if(is.logical(contours))
            contours = 1.5
        tmp = lod
        tmp[row(lod) < col(lod)] <- NA
        if(length(chr) > 1)
            thepos <- 1:ncol(lod)
        else thepos <- map[,2]
        contour(thepos, thepos, tmp, add = TRUE,drawlabels=FALSE,
                levels = max(tmp,na.rm=TRUE) - contours, col = "blue",
                lwd = 2)
        tmp = lod
        tmp[row(lod) > col(lod)] <- NA
        contour(thepos, thepos, tmp, add = TRUE,drawlabels=FALSE,
                levels = max(tmp,na.rm=TRUE) - contours * zlim.jnt/zlim.int,
                col = "blue", lwd = 2)
    }

    if(length(chr) > 1) {
        # calculate how many markers in each chromesome
        n.mar <- NULL
        for(i in 1:length(chr)) n.mar[i] <- sum(map[, 1] == chr[i])

        # plot lines at the chromosome boundaries
        if(length(chr) > 1)
            wh <- c(0.5, cumsum(n.mar) + 0.5)
        abline(v = wh, xpd = FALSE)
        abline(h = wh, xpd = FALSE)

        # add chromesome numbers
        a <- par("usr")
        placement <- (wh[-1] + wh[-length(wh)])/2
        if(!alternate.chrid || length(chr)<2) {
            for(i in 1:length(chr)) {
                axis(side=1, at=placement[i], labels=chr[i])
                axis(side=2, at=placement[i], labels=chr[i])
            }
        }
        else {
            odd <- seq(1, length(chr), by=2)
            even <- seq(2, length(chr), by=2)
            for(i in odd) {
                axis(side=1, at=placement[i], labels="")
                axis(side=2, at=placement[i], labels="")
                axis(side=1, at=placement[i], labels=chr[i], line=-0.4, tick=FALSE)
                axis(side=2, at=placement[i], labels=chr[i], line=-0.4, tick=FALSE)
            }

            for(i in even) {
                axis(side=1, at=placement[i], labels="")
                axis(side=2, at=placement[i], labels="")
                axis(side=1, at=placement[i], labels=chr[i], line=+0.4, tick=FALSE)
                axis(side=2, at=placement[i], labels=chr[i], line=+0.4, tick=FALSE)
            }
        }
    }
    else {
        u <- par("usr")
        abline(v=u[1:2], h=u[3:4])
    }

    # add title
    if(!missing(main))
        title(main = main)

    if(zscale) {
        # plot the colormap
        dots <- list(...)
        if("mar2" %in% names(dots))
            par(mar=dots[["mar2"]])
        else
            par(mar = c(5, 2, 4, 2) + 0.1)

        colorstep <- (zlim.jnt-lo)/255
        image(x = 1:1, y = seq(lo, zlim.jnt, colorstep), z = matrix(1:256, 1, 256),
              zlim = c(1, 256), ylab = "", xlab = "",
              xaxt = "n", yaxt = "n", col = cols)

        # make sure there's a box around it
        u <- par("usr")
        abline(v = u[1:2], xpd = FALSE)
        abline(h = u[3:4], xpd = FALSE)
        if(any(contours) > 0) {
            for(i in seq(length(contours))) {
                segments(mean(u[1:2]),
                         max(lod[row(lod) > col(lod)]) - contours[i],
                         u[2], max(lod[row(lod) > col(lod)]) - contours[i],
                         xpd = FALSE, col = "blue", lwd = 2)
                segments(u[1], max(lod[row(lod) < col(lod)]) - contours[i] *
                         zlim.jnt/zlim.int, mean(u[1:2]),
                         max(lod[row(lod) < col(lod)]) - contours[i] * zlim.jnt / zlim.int,
                         xpd = FALSE, col = "blue", lwd = 2)
            }
        }

        # figure out how big the axis labels should be
        fin <- par("fin")[1] # figure width in inches
        pin <- par("pin")[1] # plot width in inches
        mai <- par("mai")[2] # margin width in inches
        # note: pin + 2*mai = fin
        xlen.mar <- mai/pin * diff(u[1:2])

        # axis for full LODs
        yloc <- pretty(c(lo, zlim.jnt), 4)
        yloc <- yloc[yloc >= u[3] & yloc <= u[4]]
        #    segments(u[2], yloc, u[2] + xlen.mar/4, yloc, xpd = TRUE)
        #    text(u[2] + xlen.mar/3, yloc, as.character(yloc), xpd = TRUE, adj = 0)
        axis(side=4, at=yloc, labels=yloc)

        # axis for int've LODs
        yloc <- pretty(c(lo.int, zlim.int), 4)
        yloc.rev <- yloc * zlim.jnt/zlim.int
        yloc <- yloc[yloc.rev >= u[3] & yloc.rev <= u[4]]
        yloc.rev <- yloc.rev[yloc.rev >= u[3] & yloc.rev <= u[4]]
        #    segments(u[1], yloc.rev, u[1] - xlen.mar/4, yloc.rev, xpd = TRUE)
        #    text(u[1] - xlen.mar/3, yloc.rev, as.character(yloc), xpd = TRUE, adj = 1)
        if(!addpair)
            axis(side=2, at=yloc.rev, labels=yloc)
    }

    invisible()
}


# end of plot.scantwo.R
