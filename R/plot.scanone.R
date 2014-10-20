#####################################################################
#
# plot.scanone.R
#
# copyright (c) 2001-2014, Karl W Broman
# last modified Apr, 2014
# first written Feb, 2001
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
# Contains: plot.scanone,
#
######################################################################

######################################################################
#
# plot.scanone: plot output from scanone
#
######################################################################

plot.scanone <-
    function(x,x2,x3,chr,lodcolumn=1,incl.markers=TRUE,xlim, ylim,
             lty=1,col=c("black","blue","red"),lwd=2,add=FALSE,gap=25,
             mtick=c("line", "triangle"), show.marker.names=FALSE,
             alternate.chrid=FALSE, bandcol=NULL, type="l", cex=1,
             pch=1, bg="transparent", bgrect=NULL, ...)
{
    if(!any(class(x) == "scanone") ||
       (!missing(x2) && !any(class(x2) == "scanone")) ||
       (!missing(x3) && !any(class(x3) == "scanone")))
        stop("Input should have class \"scanone\".")

    if(!is.factor(x$chr)) x$chr <- factor(x$chr, levels=unique(x$chr))

    dots <- list(...)

    # handle special arguments to be used in lines()
    if(length(type)==1) type <- rep(type, 3)
    if(length(cex)==1) cex <- rep(cex, 3)
    if(length(pch)==1) pch <- rep(pch, 3)
    if(length(bg)==1) bg <- rep(bg, 3)

    mtick <- match.arg(mtick)

    if(length(dim(x))!=2)
        stop("Argument x must be a matrix or data.frame.")
    if(!missing(x2) && length(dim(x2))!=2)
        stop("Argument x2 must be a matrix or data.frame.")
    if(!missing(x3) && length(dim(x3))!=2)
        stop("Argument x3 must be a matrix or data.frame.")

    if(length(lodcolumn)==1)
        lodcolumn <- rep(lodcolumn,3)[1:3]
    else if(length(lodcolumn)==2) {
        if(missing(x2)) x2 <- x
        lodcolumn <- lodcolumn[c(1,2,3)]
    }
    else {
        if(missing(x2)) x2 <- x
        if(missing(x3)) x3 <- x
    }
    lodcolumn <- lodcolumn+2

    second <- third <- TRUE
    if(missing(x2) && missing(x3))
        second <- third <- FALSE
    if(missing(x3))
        third <- FALSE
    if(missing(x2))
        second <- FALSE

    # rename things and turn into data frames
    if(lodcolumn[1] > ncol(x) ||
       (second && lodcolumn[2] > ncol(x2)) ||
       (third && lodcolumn[3] > ncol(x3)))
        stop("Argument lodcolumn misspecified.")

    out <- x[,c(1:2,lodcolumn[1])]
    if(second) out2 <- x2[,c(1:2,lodcolumn[2])]
    if(third) out3 <- x3[,c(1:2,lodcolumn[3])]
    if(length(lty)==1) lty <- rep(lty,3)
    if(length(lwd)==1) lwd <- rep(lwd,3)
    if(length(col)==1) col <- rep(col,3)

    # pull out desired chromosomes
    if(missing(chr) || length(chr)==0)
        chr <- unique(as.character(out[,1]))
    else chr <- matchchr(chr, unique(out[,1]))

    out <- out[!is.na(match(out[,1],chr)),]
    if(second) out2 <- out2[!is.na(match(out2[,1],chr)),]
    if(third) out3 <- out3[!is.na(match(out3[,1],chr)),]

    onechr <- FALSE
    if(length(chr) == 1) {
        gap <- 0
        onechr <- TRUE
    }

    # beginning and end of chromosomes
    temp <- out
    begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
    rownames(begend) <- unique(out[,1])
    begend <- begend[as.character(chr),,drop=FALSE]
    len <- begend[,2]-begend[,1]

    # locations to plot start of each chromosome
    if(!onechr)
        start <- c(0,cumsum(len+gap))-c(begend[,1],0)
    else start <- 0

    maxx <- sum(len+gap)-gap
    if(all(is.na(out[,3]))) maxy <- 1
    else maxy <- max(out[,3],na.rm=TRUE)
    if(second) maxy <- max(c(maxy,out2[,3]),na.rm=TRUE)
    if(third) maxy <- max(c(maxy,out3[,3]),na.rm=TRUE)

    # graphics parameters
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd=FALSE,las=1)
    on.exit(par(xpd=old.xpd,las=old.las))

    # make frame of plot
    if(missing(ylim)) ylim <- c(0,maxy)
    if(missing(xlim)) {
        if(onechr) xlim <- c(0,max(out[,2]))
        else xlim <- c(-gap/2,maxx+gap/2)
    }

    if(!add) {
        if(onechr) {

            if("ylab" %in% names(dots)) {
                if("xlab" %in% names(dots)) {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",...)
                }
                else {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",
                         xlab="Map position (cM)",...)
                }
            }
            else {
                if("xlab" %in% names(dots)) {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",
                         ylab=dimnames(out)[[2]][3],
                         ...)
                }
                else {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",
                         xlab="Map position (cM)",ylab=dimnames(out)[[2]][3],
                         ...)
                }
            }
        }
        else {
            if("ylab" %in% names(dots)) {
                if("xlab" %in% names(dots)) {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",xaxt="n",
                         xaxs="i", ...)
                }
                else {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",xaxt="n",
                         xlab="Chromosome", xaxs="i", ...)
                }
            }
            else {
                if("xlab" %in% names(dots)) {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",xaxt="n",
                         ylab=dimnames(out)[[2]][3], xaxs="i",
                         ...)
                }
                else {
                    plot(0,0,ylim=ylim,xlim=xlim,type="n",xaxt="n",
                         xlab="Chromosome",ylab=dimnames(out)[[2]][3], xaxs="i",
                         ...)
                }
            }
        }
    }

    if(!add && !is.null(bgrect)) {
        u <- par("usr")
        rect(u[1], u[3], u[2], u[4], col=bgrect)
    }
    if(!add && !onechr && !is.null(bandcol)) {
        u <- par("usr")
        for(i in seq(2, by=2, length(chr))) {
            rect(min(out[out[,1]==chr[i],2]) + start[i]-gap/2, u[3],
                 max(out[out[,1]==chr[i],2]) + start[i]+gap/2, u[4], border=bandcol, col=bandcol)
        }
        abline(h=u[3:4])
    }

    # initialize xtick and xtickmark
    xtick <- NULL
    xticklabel <- NULL
    for(i in 1:length(chr)) {
        # plot first out
        x <- out[out[,1]==chr[i],2]+start[i]
        y <- out[out[,1]==chr[i],3]
        if(length(x)==1) {
            g <- max(gap/10,2)
            x <- c(x-g,x,x+g)
            y <- rep(y,3)
        }
        lines(x,y,lwd=lwd[1],lty=lty[1],col=col[1], type=type[1], cex=cex[1], pch=pch[1], bg=bg[1])
        # plot chromosome number
        if(!add && !onechr) {
            tloc <- mean(c(min(x),max(x)))
            xtick <- c(xtick, tloc)
            xticklabel <- c(xticklabel, as.character(chr[i]))
        }

        # plot second out
        if(second) {
            x <- out2[out2[,1]==chr[i],2]+start[i]
            y <- out2[out2[,1]==chr[i],3]
            if(length(x)==1) {
                g <- max(gap/10,2)
                x <- c(x-g,x,x+g)
                y <- rep(y,3)
            }
            lines(x,y,lty=lty[2],col=col[2],lwd=lwd[2], type=type[2], cex=cex[2], pch=pch[2], bg=bg[2])
        }

        if(third) {
            x <- out3[out3[,1]==chr[i],2]+start[i]
            y <- out3[out3[,1]==chr[i],3]
            if(length(x)==1) {
                g <- max(gap/10,2)
                x <- c(x-g,x,x+g)
                y <- rep(y,3)
            }
            lines(x,y,lty=lty[3],col=col[3],lwd=lwd[3], type=type[3], cex=cex[3], pch=pch[3], bg=bg[3])
        }

        # plot lines or triangles at marker positions
        if(!add) {
            nam <- dimnames(out)[[1]][out[,1]==chr[i]]
            wh.genoprob <- grep("^c.+\\.loc-*[0-9]+", nam)
            if(length(wh.genoprob)==0) wh.genoprob <- seq(along=nam)
            else wh.genoprob <- (seq(along=nam))[-wh.genoprob]
            pos <- out[out[,1]==chr[i],2][wh.genoprob]+start[i]

            if(incl.markers) {
                if(mtick=="line")
                    rug(pos, 0.02, quiet=TRUE)
                else {
                    a <- par("usr")
                    points(pos, rep(a[3]+diff(a[3:4])*0.01, length(pos)), pch=17, cex=1.5)
                }
            }

            if(show.marker.names) {
                a <- par("usr")
                text(pos, rep(a[3]+diff(a[3:4])*0.03, length(pos)), nam[wh.genoprob],
                     srt=90, adj=c(0,0.5))
            }
        }
    }

    # draw the axis
    if(!add && !onechr) {
        if(!alternate.chrid || length(xtick) < 2) {
            for(i in seq(along=xtick))
                axis(side=1, at=xtick[i], labels=xticklabel[i])
        }
        else {
            odd <- seq(1, length(xtick), by=2)
            even <- seq(2, length(xtick), by=2)
            for(i in odd) {
                axis(side=1, at=xtick[i], labels="")
                axis(side=1, at=xtick[i], labels=xticklabel[i], line=-0.4, tick=FALSE)
            }
            for(i in even) {
                axis(side=1, at=xtick[i], labels="")
                axis(side=1, at=xtick[i], labels=xticklabel[i], line=+0.4, tick=FALSE)
            }
        }
    }

    if(!add && !is.null(bgrect)) {
        u <- par("usr")
        rect(u[1], u[3], u[2], u[4])
    }

    invisible()
}

# end of plot.scanone.R
