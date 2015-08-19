######################################################################
#
# plot.R
#
# copyright (c) 2000-2015, Karl W Broman
#       [modifications of plot.cross from Brian Yandell]
# last modified Aug, 2015
# first written Mar, 2000
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
# Contains: plotMissing, plotMap, plot.cross, plotGeno, plotInfo,
#           plotPXG, plotPheno
#
######################################################################

plotMissing <- plot.missing <-
    function(x, chr, reorder=FALSE, main="Missing genotypes",
             alternate.chrid=FALSE, ...)
{
    cross <- x
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if(!missing(chr)) cross <- subset(cross,chr=chr)

    # get full genotype data into one matrix
    Geno <- pull.geno(cross)

    # reorder the individuals according to their phenotype
    o <- 1:nrow(Geno)
    if(reorder) {
        # if reorder is a number, use the corresponding phenotype
        if(is.numeric(reorder)) {
            if(reorder < 1 || reorder > nphe(cross))
                stop("reorder should be TRUE, FALSE, or an integer between 1 and", nphe(cross))

            o <- order(cross$pheno[,reorder])
        }

        # otherwise, order according to the sum of the numeric phenotypes
        else {
            wh <- sapply(cross$pheno, is.numeric)
            o <- order(apply(cross$pheno[,wh,drop=FALSE],1,sum))
        }
    }

    # make matrix with  0 where genotype data is missing
    #                   1 where data is not missing
    #                 0.5 where data is partially missing
    type <- class(cross)[1]
    g <- t(Geno[o,])
    g[is.na(g)] <- 0
    if(type == "bc" || type=="risib" || type=="riself" || type=="bc")
        g[g > 0] <- 1
    else if(type=="f2") {
        g[g > 0 & g < 4] <- 1
        g[g > 3] <- 0.5
    }
    else if(type=="4way") {
        g[g > 0 & g < 5] <- 1
        g[g > 4] <- 0.5
    }
    else {
        g[g > 0] <- 1
    }

    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd=TRUE,las=1)
    on.exit(par(xpd=old.xpd,las=old.las))

    colors <- c("#000000", "gray80", "#FFFFFF")

    # plot grid with black pixels where there is missing data
    image(1:nrow(g),1:ncol(g),g,ylab="Individuals",xlab="Markers",col=colors,zlim=c(0,1))

    # plot lines at the chromosome boundaries
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)
    a <- c(0.5,cumsum(n.mar)+0.5)

    # the following makes the lines go slightly above the plotting region
    b <- par("usr")
    segments(a,b[3],a,b[4]+diff(b[3:4])*0.02)

    # this line adds a line above the image
    #     (the image function seems to leave it out)
    abline(h=0.5+c(0,ncol(g)),xpd=FALSE)

    # add chromosome numbers
    a <- par("usr")
    wh <- cumsum(c(0.5,n.mar))
    x <- 1:n.chr
    for(i in 1:n.chr)
        x[i] <- mean(wh[i+c(0,1)])

    thechr <- names(cross$geno)
    if(!alternate.chrid || length(thechr) < 2) {
        for(i in seq(along=x))
            axis(side=3, at=x[i], thechr[i], tick=FALSE, line=-0.5)
    }
    else {
        odd <- seq(1, length(x), by=2)
        even <- seq(2, length(x), by=2)
        for(i in odd) {
            axis(side=3, at=x[i], labels=thechr[i], line=-0.75, tick=FALSE)
        }
        for(i in even) {
            axis(side=3, at=x[i], labels=thechr[i], line=+0, tick=FALSE)
        }
    }

    title(main=main)
    invisible()
}

geno.image <-
    function(x, chr, reorder=FALSE, main="Genotype data",
             alternate.chrid=FALSE, ...)
{
    cross <- x
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if(!missing(chr)) cross <- subset(cross,chr=chr)

    type <- class(cross)[1]

    # revise X chromosome data
    if(type=="bc" || type=="f2") {
        chrtype <- sapply(cross$geno, class)
        if(any(chrtype=="X")) {
            for(i in which(chrtype=="X"))
                cross$geno[[i]]$data <- reviseXdata(type, "simple", getsex(cross),
                                                    geno=cross$geno[[i]]$data,
                                                    cross.attr=attributes(cross))
        }
    }

    # get full genotype data into one matrix
    Geno <- pull.geno(cross)

    # colors to use
    maxgeno <- max(Geno, na.rm=TRUE)
    if(type != "4way") {
        thecolors <- c("white", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
        thebreaks <- seq(-0.5, 5.5, by=1)
    }
    else {
        if(maxgeno <= 5) {
            thecolors <- c("white", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
            thebreaks <- seq(-0.5, 5.5, by=1)
        }
        else {
            thecolors <- c("white",  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                           "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
            thebreaks <- seq(-0.5, 10.5, by=1)
        }
    }
    thecolors <- thecolors[1:(maxgeno+1)]
    thebreaks <- thebreaks[1:(maxgeno+2)]

    # reorder the individuals according to their phenotype
    o <- 1:nrow(Geno)
    if(reorder) {
        # if reorder is a number, use the corresponding phenotype
        if(is.numeric(reorder)) {
            if(reorder < 1 || reorder > nphe(cross))
                stop("reorder should be TRUE, FALSE, or an integer between 1 and", nphe(cross))

            o <- order(cross$pheno[,reorder])
        }

        # otherwise, order according to the sum of the numeric phenotypes
        else {
            wh <- sapply(cross$pheno, is.numeric)
            o <- order(apply(cross$pheno[,wh,drop=FALSE],1,sum))
        }
    }

    g <- t(Geno[o,])
    g[is.na(g)] <- 0

    # make matrix with  0 where genotype data is missing
    #                   1 where data is not missing
    #                 0.5 where data is partially missing
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd=TRUE,las=1)
    on.exit(par(xpd=old.xpd,las=old.las))

    # plot grid with black pixels where there is missing data
    plot_image_sub <-
        function(g, ylab="Individuals",xlab="Markers", col=thecolors, ...)
        {
            if(length(thebreaks) != length(col)+1)
                stop("Must have one more break than color\n",
                     "length(breaks) = ", length(thebreaks),
                     "\nlength(col) = ", length(col))
            image(1:nrow(g),1:ncol(g), g, col=col, xlab=xlab, ylab=ylab,
                  breaks=thebreaks, ...)
        }
    plot_image_sub(g, ...)

    # plot lines at the chromosome boundaries
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)
    a <- c(0.5,cumsum(n.mar)+0.5)

    # the following makes the lines go slightly above the plotting region
    b <- par("usr")
    segments(a,b[3],a,b[4]+diff(b[3:4])*0.02)

    # this line adds a line above the image
    #     (the image function seems to leave it out)
    abline(h=0.5+c(0,ncol(g)),xpd=FALSE)

    # add chromosome numbers
    a <- par("usr")
    wh <- cumsum(c(0.5,n.mar))
    x <- 1:n.chr
    for(i in 1:n.chr)
        x[i] <- mean(wh[i+c(0,1)])
    thechr <- names(cross$geno)
    if(!alternate.chrid || length(thechr) < 2) {
        for(i in seq(along=x))
            axis(side=3, at=x[i], thechr[i], tick=FALSE, line=-0.5)
    }
    else {
        odd <- seq(1, length(x), by=2)
        even <- seq(2, length(x), by=2)
        for(i in odd)
            axis(side=3, at=x[i], labels=thechr[i], line=-0.75, tick=FALSE)
        for(i in even)
            axis(side=3, at=x[i], labels=thechr[i], line=+0, tick=FALSE)
    }

    title(main=main)
    invisible()
}

plotMap <- plot.map <-
    function(x, map2, chr, horizontal=FALSE, shift=TRUE,
             show.marker.names=FALSE, alternate.chrid=FALSE, ...)
{
    dots <- list(...)
    if("main" %in% names(dots)) {
        themain <- dots$main
        usemaindefault <- FALSE
    }
    else usemaindefault <- TRUE

    if("xlim" %in% names(dots)) {
        xlim <- dots$xlim
        usexlimdefault <- FALSE
    }
    else usexlimdefault <- TRUE

    if("ylim" %in% names(dots)) {
        ylim <- dots$ylim
        useylimdefault <- FALSE
    }
    else useylimdefault <- TRUE

    if("xlab" %in% names(dots))
        xlab <- dots$xlab
    else {
        if(horizontal)
            xlab <- "Location (cM)"
        else
            xlab <- "Chromosome"
    }

    if("ylab" %in% names(dots))
        ylab <- dots$ylab
    else {
        if(horizontal)
            ylab <- "Chromosome"
        else
            ylab <- "Location (cM)"
    }

    map <- x
    # figure out if the input is a cross (containing a map)
    #    or is the map itself
    if(any(class(map) == "cross"))
        map <- pull.map(map)
    if(!missing(map2) && any(class(map2) == "cross"))
        map2 <- pull.map(map2)

    if(!any(class(map)=="map")  || (!missing(map2) && !any(class(map2) == "map")))
        warning("Input should have class \"cross\" or \"map\".")

    if(!missing(map2) && is.matrix(map[[1]]) != is.matrix(map2[[1]]))
        stop("Maps must be both sex-specific or neither sex-specific.")

    if(!missing(chr)) {
        map <- map[matchchr(chr, names(map))]
        if(!missing(map2)) map2 <- map2[matchchr(chr, names(map2))]
    }

    sex.sp <- FALSE

    if(is.matrix(map[[1]])) { # sex-specific map
        one.map <- FALSE
        sex.sp <- TRUE
        if(!missing(map2)) {
            if(is.logical(map2)) {
                horizontal <- map2
                map2 <- lapply(map,function(a) a[2,])
                map <- lapply(map,function(a) a[1,])
            }
            else {
                Map1 <- lapply(map,function(a) a[1,,drop=TRUE])
                Map2 <- lapply(map,function(a) a[2,,drop=TRUE])
                Map3 <- lapply(map2,function(a) a[1,,drop=TRUE])
                Map4 <- lapply(map2,function(a) a[2,,drop=TRUE])
                old.mfrow <- par("mfrow")
                on.exit(par(mfrow=old.mfrow))

                par(mfrow=c(2,1))
                class(Map1) <- class(Map2) <- class(Map3) <- class(Map4) <- "map"
                plotMap(Map1,Map3,horizontal=horizontal,shift=shift,
                        show.marker.names=show.marker.names,alternate.chrid=alternate.chrid)
                plotMap(Map2,Map4,horizontal=horizontal,shift=shift,
                        show.marker.names=show.marker.names,alternate.chrid=alternate.chrid)
                return(invisible(NULL))
            }
        }
        else {
            map2 <- lapply(map,function(a) a[2,])
            map <- lapply(map,function(a) a[1,])
        }
    }
    else { # single map
        # determine whether a second map was given
        if(!missing(map2))
            one.map <- FALSE
        else one.map <- TRUE
    }

    if(one.map) {

        n.chr <- length(map)
        if(!show.marker.names) {  # locations of chromosomes
            chrpos <- 1:n.chr
            thelim <- range(chrpos)+c(-0.5, 0.5)
        }
        else {
            chrpos <- seq(1, n.chr*2, by=2)
            thelim <- range(chrpos)+c(-0.35, 2.35)
        }

        if(shift) map <- lapply(map, function(a) a-a[1])
        maxlen <- max(unlist(lapply(map,max)))

        if(horizontal) {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd=TRUE, las=1)
            on.exit(par(xpd=old.xpd,las=old.las))

            if(usexlimdefault) xlim <- c(0,maxlen)
            if(useylimdefault) ylim <- rev(thelim)
            plot(0,0,type="n",xlim=xlim, ylim=ylim,yaxs="i",
                 xlab=xlab, ylab=ylab, yaxt="n")
            a <- par("usr")

            for(i in 1:n.chr) {
                segments(min(map[[i]]), chrpos[i], max(map[[i]]), chrpos[i])
                segments(map[[i]], chrpos[i]-0.25, map[[i]], chrpos[i]+0.25)

                if(show.marker.names)
                    text(map[[i]], chrpos[i]+0.35, names(map[[i]]), srt=90, adj=c(1,0.5))
            }

            # add chromosome labels
            if(!alternate.chrid || length(chrpos) < 2) {
                for(i in seq(along=chrpos))
                    axis(side=2, at=chrpos[i], labels=names(map)[i])
            }
            else {
                odd <- seq(1, length(chrpos), by=2)
                even <- seq(2, length(chrpos), by=2)
                for(i in odd) {
                    axis(side=2, at=chrpos[i], labels="")
                    axis(side=2, at=chrpos[i], labels=names(map)[i], line=-0.4, tick=FALSE)
                }
                for(i in even) {
                    axis(side=2, at=chrpos[i], labels="")
                    axis(side=2, at=chrpos[i], labels=names(map)[i], line=+0.4, tick=FALSE)
                }
            }
        }
        else {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd=TRUE,las=1)
            on.exit(par(xpd=old.xpd,las=old.las))

            if(usexlimdefault) xlim <- thelim
            if(useylimdefault) ylim <- c(maxlen, 0)
            plot(0,0,type="n",ylim=ylim,xlim=xlim,xaxs="i",
                 xlab=xlab, ylab=ylab, xaxt="n")

            a <- par("usr")

            for(i in 1:n.chr) {
                segments(chrpos[i], min(map[[i]]), chrpos[i], max(map[[i]]))
                segments(chrpos[i]-0.25, map[[i]], chrpos[i]+0.25, map[[i]])

                if(show.marker.names)
                    text(chrpos[i]+0.35, map[[i]], names(map[[i]]), adj=c(0,0.5))

            }
            # add chromosome labels
            if(!alternate.chrid || length(chrpos) < 2) {
                for(i in seq(along=chrpos))
                    axis(side=1, at=chrpos[i], labels=names(map)[i])
            }
            else {
                odd <- seq(1, length(chrpos), by=2)
                even <- seq(2, length(chrpos), by=2)
                for(i in odd) {
                    axis(side=1, at=chrpos[i], labels="")
                    axis(side=1, at=chrpos[i], labels=names(map)[i], line=-0.4, tick=FALSE)
                }
                for(i in even) {
                    axis(side=1, at=chrpos[i], labels="")
                    axis(side=1, at=chrpos[i], labels=names(map)[i], line=+0.4, tick=FALSE)
                }
            }
        }
        if(usemaindefault)
            title(main="Genetic map")
        else if(themain != "")
            title(main=themain)

    }
    else {
        map1 <- map

        # check that maps conform
        if(is.matrix(map2[[1]]))
            stop("Second map appears to be a sex-specific map.")
        if(length(map1) != length(map2))
            stop("Maps have different numbers of chromosomes.")
        if(any(names(map1) != names(map2))) {
            cat("Map1: ", names(map1), "\n")
            cat("Map2: ", names(map2), "\n")
            stop("Maps have different chromosome names.")
        }

        if(shift) {
            map1 <- lapply(map1,function(a) a-a[1])
            map2 <- lapply(map2,function(a) a-a[1])
        }

        n.mar1 <- sapply(map1, length)
        n.mar2 <- sapply(map2, length)
        markernames1 <- lapply(map1, names)
        markernames2 <- lapply(map2, names)
        if(any(n.mar1 != n.mar2)) {
            if(show.marker.names) {
                warning("Can't show marker names because of different numbers of markers.")
                show.marker.names <- FALSE
            }
        }
        else if(any(unlist(markernames1) != unlist(markernames2))) {
            if(show.marker.names) {
                warning("Can't show marker names because markers in different orders.")
                show.marker.names <- FALSE
            }
        }

        n.chr <- length(map1)
        maxloc <- max(c(unlist(lapply(map1,max)),unlist(lapply(map2,max))))

        if(!show.marker.names) {  # locations of chromosomes
            chrpos <- 1:n.chr
            thelim <- range(chrpos)+c(-0.5, 0.5)
        }
        else {
            chrpos <- seq(1, n.chr*2, by=2)
            thelim <- range(chrpos)+c(-0.4, 2.4)
        }

        if(!horizontal) {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd=TRUE,las=1)
            on.exit(par(xpd=old.xpd,las=old.las))

            if(usexlimdefault) xlim <- thelim
            if(useylimdefault) ylim <- c(maxloc, 0)
            plot(0,0,type="n",ylim=ylim,xlim=xlim, xaxs="i",
                 xlab=xlab, ylab=ylab, xaxt="n")

            a <- par("usr")

            for(i in 1:n.chr) {
                if(max(map2[[i]]) < max(map1[[i]]))
                    map2[[i]] <- map2[[i]] + (max(map1[[i]])-max(map2[[i]]))/2

                else
                    map1[[i]] <- map1[[i]] + (max(map2[[i]])-max(map1[[i]]))/2

                segments(chrpos[i]-0.3, min(map1[[i]]), chrpos[i]-0.3, max(map1[[i]]))
                segments(chrpos[i]+0.3, min(map2[[i]]), chrpos[i]+0.3, max(map2[[i]]))

                # lines between markers
                wh <- match(markernames1[[i]], markernames2[[i]])
                for(j in which(!is.na(wh)))
                    segments(chrpos[i]-0.3, map1[[i]][j], chrpos[i]+0.3, map2[[i]][wh[j]])
                if(any(is.na(wh)))
                    segments(chrpos[i]-0.4, map1[[i]][is.na(wh)], chrpos[i]-0.2, map1[[i]][is.na(wh)])
                wh <- match(markernames2[[i]], markernames1[[i]])
                if(any(is.na(wh)))
                    segments(chrpos[i]+0.4, map2[[i]][is.na(wh)], chrpos[i]+0.2, map2[[i]][is.na(wh)])

                if(show.marker.names)
                    text(chrpos[i]+0.35, map2[[i]], names(map2[[i]]), adj=c(0,0.5))

            }
            # add chromosome labels
            if(!alternate.chrid || length(chrpos) < 2) {
                for(i in seq(along=chrpos))
                    axis(side=1, at=chrpos[i], labels=names(map1)[i])
            }
            else {
                odd <- seq(1, length(chrpos), by=2)
                even <- seq(2, length(chrpos), by=2)
                for(i in odd) {
                    axis(side=1, at=chrpos[i], labels="")
                    axis(side=1, at=chrpos[i], labels=names(map1)[i], line=-0.4, tick=FALSE)
                }
                for(i in even) {
                    axis(side=1, at=chrpos[i], labels="")
                    axis(side=1, at=chrpos[i], labels=names(map1)[i], line=+0.4, tick=FALSE)
                }
            }
        }
        else {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd=TRUE,las=1)
            on.exit(par(xpd=old.xpd,las=old.las))

            if(usexlimdefault) xlim <- c(0,maxloc)
            if(useylimdefault) ylim <- rev(thelim)
            plot(0,0,type="n",xlim=xlim,ylim=ylim,
                 xlab=xlab, ylab=ylab, yaxt="n", yaxs="i")

            a <- par("usr")

            for(i in 1:n.chr) {

                if(max(map2[[i]]) < max(map1[[i]]))
                    map2[[i]] <- map2[[i]] + (max(map1[[i]])-max(map2[[i]]))/2
                else
                    map1[[i]] <- map1[[i]] + (max(map2[[i]])-max(map1[[i]]))/2

                segments(min(map1[[i]]), chrpos[i]-0.3, max(map1[[i]]), chrpos[[i]]-0.3)
                segments(min(map2[[i]]), chrpos[i]+0.3, max(map2[[i]]), chrpos[[i]]+0.3)

                # lines between markers
                wh <- match(markernames1[[i]], markernames2[[i]])
                for(j in which(!is.na(wh)))
                    segments(map1[[i]][j], chrpos[i]-0.3, map2[[i]][wh[j]], chrpos[i]+0.3)
                if(any(is.na(wh)))
                    segments(map1[[i]][is.na(wh)], chrpos[i]-0.4, map1[[i]][is.na(wh)], chrpos[i]-0.2)
                wh <- match(markernames2[[i]], markernames1[[i]])
                if(any(is.na(wh)))
                    segments(map2[[i]][is.na(wh)], chrpos[i]+0.4, map2[[i]][is.na(wh)], chrpos[i]+0.2)

                if(show.marker.names)
                    text(map2[[i]], chrpos[i]+0.35, names(map2[[i]]), srt=90, adj=c(1,0.5))

            }
            # add chromosome labels
            if(!alternate.chrid || length(chrpos) < 2) {
                for(i in seq(along=chrpos))
                    axis(side=2, at=chrpos[i], labels=names(map1)[i])
            }
            else {
                odd <- seq(1, length(chrpos), by=2)
                even <- seq(2, length(chrpos), by=2)
                for(i in odd) {
                    axis(side=2, at=chrpos[i], labels="")
                    axis(side=2, at=chrpos[i], labels=names(map1)[i], line=-0.4, tick=FALSE)
                }
                for(i in even) {
                    axis(side=2, at=chrpos[i], labels="")
                    axis(side=2, at=chrpos[i], labels=names(map1)[i], line=+0.4, tick=FALSE)
                }
            }

        }
        if(usemaindefault) {
            if(!sex.sp) title(main="Comparison of genetic maps")
            else title(main="Genetic map")
        }
        else if(themain != "")
            title(main=themain)
    }
    invisible()
}


plot.cross <-
    function (x, auto.layout = TRUE, pheno.col,
              alternate.chrid=TRUE, ...)
{
    # look to see whether this should really be shipped to plotMap
    if("map" %in% class(auto.layout) &&
       ("map" %in% class(x) || "cross" %in% class(x))) {
        plotMap(x, auto.layout, alternate.chrid=alternate.chrid, ...)
        return(invisible())
    }

    # make sure this is a cross
    if(!("cross" %in% class(x)))
        stop("Input should have class \"cross\".")

    old.yaxt <- par("yaxt")
    old.mfrow <- par("mfrow")
    on.exit(par(yaxt = old.yaxt, mfrow = old.mfrow))

    n.phe <- nphe(x)
    if(missing(pheno.col)) pheno.col <- 1:n.phe
    if(is.character(pheno.col)) {
        temp <- match(pheno.col, names(x$pheno))
        if(any(is.na(temp)))
            warning("Some phenotypes not found:",
                    paste(pheno.col[is.na(temp)], collapse=" "))
        pheno.col <- temp[!is.na(temp)]
    }
    else pheno.col <- (1:nphe(x))[pheno.col]

    n.plot = length(pheno.col) + 2

    # automatically choose row/column structure for the plots
    if(auto.layout) {
        nr <- ceiling(sqrt(n.plot))
        nc <- ceiling((n.plot)/nr)
        par(mfrow = c(nr, nc))
    }

    plotMissing(x,alternate.chrid=alternate.chrid)
    plotMap(x,alternate.chrid=alternate.chrid)

    for(i in pheno.col) plotPheno(x, pheno.col=i)

    invisible()
}


##################################################r####################
#
# plotGeno: Plot genotypes for a specified chromosome, with likely
#           genotyping errors indicated.
#
######################################################################

plotGeno <- plot.geno <-
    function(x, chr, ind, include.xo=TRUE, horizontal=TRUE,
             cutoff=4, min.sep=2, cex=1.2, ...)
{
    cross <- x
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(missing(chr)) chr <- names(cross$geno)[1]
    cross <- subset(cross,chr=chr)
    if(nchr(cross) > 1)
        cross <- subset(cross,chr=names(cross$geno)[1])

    if(!missing(ind)) {
        if(is.null(getid(cross))) cross$pheno$id <- 1:nind(cross)
        if(!is.logical(ind)) ind <- unique(ind)
        cross <- subset(cross, ind=ind)
    }
    id <- getid(cross)
    if(is.null(id)) id <- 1:nind(cross)
    use.id <- TRUE

    type <- class(cross)[1]

    old.las <- par("las")
    on.exit(par(las=old.las))
    par(las=1)

    if(!("errorlod" %in% names(cross$geno[[1]]))) {
        warning("First running calc.errorlod.")
        cross <- calc.errorlod(cross,error.prob=0.01)
    }

    # indicators for apparent errors
    errors <- matrix(0,ncol=ncol(cross$geno[[1]]$data),
                     nrow=nrow(cross$geno[[1]]$data))
    dimnames(errors) <- dimnames(cross$geno[[1]]$data)

    top <- top.errorlod(cross,names(cross$geno)[1],cutoff,FALSE)
    if(length(top) > 0)
        for(i in 1:nrow(top))
            errors[match(top[i,2],id),as.character(top[i,3])] <- 1

    # map, data, errors
    map <- cross$geno[[1]]$map
    if(is.matrix(map)) map <- map[1,] # if sex-specific map
    L <- diff(range(map))
    min.d <- L*min.sep/100
    d <- diff(map)
    d[d < min.d] <- min.d
    map <- cumsum(c(0,d))
    cross$geno[[1]]$map <- map

    n.ind <- nrow(errors)

    color <- c("white","gray60","black","green","orange","red")

    # revise X chr data for backcross/intercross
    data <- cross$geno[[1]]$data
    chrtype <- class(cross$geno[[1]])
    if(chrtype=="X" && (type=="f2" || type=="bc"))
        data <- reviseXdata(type, sexpgm=getsex(cross), geno=data, cross.attr=attributes(cross), force=TRUE)

    if(include.xo) {
        if(type != "4way") { # find crossover locations
            xoloc <- locateXO(cross)
            xoloc <- data.frame(ind=rep(1:length(xoloc),sapply(xoloc,length)),
                                loc=unlist(xoloc), stringsAsFactors=TRUE)
        }
        else { # 4-way cross
            mcross <- dcross <- cross
            class(mcross)[1] <- class(dcross)[1] <- "bc"
            mcross$geno[[1]]$data[!is.na(data) & data==1 | data==3 | data==5] <- 1
            mcross$geno[[1]]$data[!is.na(data) & data==2 | data==4 | data==6] <- 2
            mcross$geno[[1]]$data[!is.na(data) & data==7 | data==8 | data==9 | data==10] <- NA

            dcross$geno[[1]]$data[!is.na(data) & data==1 | data==2 | data==7] <- 1
            dcross$geno[[1]]$data[!is.na(data) & data==3 | data==4 | data==8] <- 2
            dcross$geno[[1]]$data[!is.na(data) & data==5 | data==6 | data==9 | data==10] <- NA

            mxoloc <- locateXO(mcross)
            mxoloc <- data.frame(ind=rep(1:length(mxoloc),sapply(mxoloc,length)),
                                 loc=unlist(mxoloc), stringsAsFactors=TRUE)

            dxoloc <- locateXO(dcross)
            dxoloc <- data.frame(ind=rep(1:length(dxoloc),sapply(dxoloc,length)),
                                 loc=unlist(dxoloc), stringsAsFactors=TRUE)
        }
    }

    # check for 'main' in the ...
    args <- list(...)
    if("main" %in% names(args))
        themain <- args$main
    else
        themain <- paste("Chromosome",names(cross$geno)[1])

    # check for 'xlim' and 'ylim'
    if("xlim" %in% names(args)) thexlim <- args$xlim
    else thexlim <- NULL
    if("ylim" %in% names(args)) theylim <- args$ylim
    else theylim <- NULL

    if(type=="4way") {
        jit <- 0.15
        mdata <- data
        ddata <- data

        # mom's allele
        mdata[!is.na(data) & (data==1 | data==3 | data==5)] <- 1
        mdata[!is.na(data) & (data==2 | data==4 | data==6)] <- 2
        mdata[!is.na(data) & (data==7 | data==8)] <- NA

        # dad's allele
        ddata[!is.na(data) & (data==1 | data==2 | data==7)] <- 1
        ddata[!is.na(data) & (data==3 | data==4 | data==8)] <- 2
        ddata[!is.na(data) & (data==5 | data==6)] <- NA

        if(horizontal) {
            if(is.null(thexlim)) thexlim <- c(0, max(map))
            if(is.null(theylim)) theylim <- c(n.ind+1, 0)
            plot(0,0,type="n",xlab="Location (cM)",ylab="Individual",
                 main=themain,
                 ylim=theylim,xlim=thexlim, yaxt="n", yaxs="i")
            segments(0, 1:n.ind-jit, max(map), 1:n.ind-jit)
            segments(0, 1:n.ind+jit, max(map), 1:n.ind+jit)

            if(use.id) axis(side=2, at=1:n.ind, labels=id)
            else axis(side=2, at=1:n.ind)

            # A alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=1] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind-jit,pch=21,col="black", bg=color[1],cex=cex)

            # B alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=2] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind-jit,pch=21,col="black", bg=color[3],cex=cex)

            # 9/10 genotypes
            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=9] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind-jit,pch=21,col="black", bg=color[4],cex=cex)

            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=10] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind-jit,pch=21,col="black", bg=color[5],cex=cex)

            # C alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=1] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind+jit,pch=21,col="black", bg=color[1],cex=cex)

            # D alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=2] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind+jit,pch=21,col="black", bg=color[3],cex=cex)

            # 9/10 genotypes
            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=9] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind+jit,pch=21,col="black", bg=color[4],cex=cex)

            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=10] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind+jit,pch=21,col="black", bg=color[5],cex=cex)

            # plot map
            u <- par("usr")
            segments(map,u[3],map,u[3]-1/2)
            segments(map,u[4],map,u[4]+1/2)

            if(any(errors != 0)) {
                ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
                points(x,ind-jit,pch=0,col=color[6],cex=cex+0.4,lwd=2)
                points(x,ind+jit,pch=0,col=color[6],cex=cex+0.4,lwd=2)
            }

            if(include.xo) {
                points(mxoloc$loc,mxoloc$ind-jit,pch=4,col="blue",lwd=2)
                points(dxoloc$loc,dxoloc$ind+jit,pch=4,col="blue",lwd=2)
            }
        }
        else {
            if(is.null(theylim)) theylim <- c(max(map), 0)
            if(is.null(thexlim)) thexlim <- c(0, n.ind+1)

            plot(0,0,type="n",ylab="Location (cM)",xlab="Individual",
                 main=themain,
                 xlim=thexlim,ylim=theylim, xaxt="n", xaxs="i")

            segments(1:n.ind-jit, 0, 1:n.ind-jit, max(map))
            segments(1:n.ind+jit, 0, 1:n.ind+jit, max(map))


            if(use.id) axis(side=1, at=1:n.ind, labels=id)
            else axis(side=1, at=1:n.ind)

            # A alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=1] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind-jit,y,pch=21,col="black",bg=color[1],cex=cex)

            # B alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=2] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind-jit,y,pch=21,col="black",bg=color[3],cex=cex)

            # 9/10 genotypes
            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=9] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind-jit,y,pch=21,col="black", bg=color[4],cex=cex)

            tind <- rep(1:n.ind,length(map));tind[is.na(mdata)] <- NA
            ind <- tind; ind[!is.na(mdata) & mdata!=10] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind-jit,y,pch=21,col="black", bg=color[5],cex=cex)

            # C alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=1] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind+jit,y,pch=21,col="black", bg=color[1],cex=cex)

            # D alleles
            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=2] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind+jit,y,pch=21,col="black", bg=color[3],cex=cex)

            # 9/10 genotypes
            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=9] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind+jit,y,pch=21,col="black", bg=color[4],cex=cex)

            tind <- rep(1:n.ind,length(map));tind[is.na(ddata)] <- NA
            ind <- tind; ind[!is.na(ddata) & ddata!=10] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind+jit,y,pch=21,col="black", bg=color[5],cex=cex)

            # plot map
            u <- par("usr")
            segments(u[1],map,(u[1]+1)/2,map)
            segments(u[2],map,(n.ind+u[2])/2,map)

            if(any(errors != 0)) {
                ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
                points(ind-jit,y,pch=0,col=color[6],cex=cex+0.4,lwd=2)
                points(ind+jit,y,pch=0,col=color[6],cex=cex+0.4,lwd=2)
            }

            if(include.xo) {
                points(mxoloc$ind-jit,mxoloc$loc,pch=4,col="blue",lwd=2)
                points(dxoloc$ind+jit,dxoloc$loc,pch=4,col="blue",lwd=2)
            }

        }


    }
    else {

        if(horizontal) {
            if(is.null(thexlim)) thexlim <- c(0, max(map))
            if(is.null(theylim)) theylim <- c(n.ind+0.5,0.5)

            plot(0,0,type="n",xlab="Location (cM)",ylab="Individual",
                 main=themain,
                 ylim=theylim,xlim=thexlim, yaxt="n")
            segments(0, 1:n.ind, max(map), 1:n.ind)
            if(use.id) axis(side=2, at=1:n.ind, labels=id)
            else axis(side=2)

            # AA genotypes
            tind <- rep(1:n.ind,length(map));tind[is.na(data)] <- NA
            ind <- tind; ind[!is.na(data) & data!=1] <- NA
            x <- rep(map,rep(n.ind,length(map)))
            points(x,ind,pch=21,col="black", bg=color[1],cex=cex)

            # AB genotypes
            ind <- tind; ind[!is.na(data) & data!=2] <- NA
            if(type=="f2" || (type=="bc" && chrtype=="X"))
                points(x,ind,pch=21,col="black", bg=color[2],cex=cex)
            else points(x,ind,pch=21,col="black", bg=color[3],cex=cex)

            if(type=="f2" || (type=="bc" && chrtype=="X")) {
                # BB genotypes
                ind <- tind; ind[!is.na(data) & data!=3] <- NA
                points(x,ind,pch=21,col="black", bg=color[3],cex=cex)
            }

            if(type=="f2") {
                # not BB (D in mapmaker/qtl) genotypes
                ind <- tind; ind[!is.na(data) & data!=4] <- NA
                points(x,ind,pch=21,col="black", bg=color[4],cex=cex)

                # not AA (C in mapmaker/qtl) genotypes
                ind <- tind; ind[!is.na(data) & data!=5] <- NA
                points(x,ind,pch=21,col="black", bg=color[5],cex=cex)
            }

            # plot map
            u <- par("usr")
            segments(map,u[3],map,u[3]-1/2)
            segments(map,u[4],map,u[4]+1/2)

            if(any(errors != 0)) {
                ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
                points(x,ind,pch=0,col=color[6],cex=cex+0.4,lwd=2)
            }

            if(include.xo) points(xoloc$loc,xoloc$ind,pch=4,col="blue",lwd=2)
        }
        else {
            if(is.null(theylim)) theylim <- c(max(map), 0)
            if(is.null(thexlim)) thexlim <- c(0.5,n.ind+0.5)
            plot(0,0,type="n",ylab="Location (cM)",xlab="Individual",
                 main=themain,
                 xlim=thexlim,ylim=theylim, xaxt="n")
            segments(1:n.ind,0,1:n.ind,max(map))
            if(use.id) axis(side=1, at=1:n.ind, labels=id)
            else axis(side=1)

            # AA genotypes
            tind <- rep(1:n.ind,length(map));tind[is.na(data)] <- NA
            ind <- tind; ind[!is.na(data) & data!=1] <- NA
            y <- rep(map,rep(n.ind,length(map)))
            points(ind,y,pch=21,col="black", bg="white",cex=cex)

            # AB genotypes
            ind <- tind; ind[!is.na(data) & data!=2] <- NA
            if(type=="f2" || (type=="bc" && chrtype=="X"))
                points(ind,y,pch=21,col="black", bg=color[2],cex=cex)
            else points(ind,y,pch=21,col="black", bg=color[3],cex=cex)

            if(type=="f2" || (type=="bc" && chrtype=="X")) {
                # BB genotypes
                ind <- tind; ind[!is.na(data) & data!=3] <- NA
                points(ind,y,pch=21,col="black", bg=color[3],cex=cex)
            }

            if(type=="f2") {
                # not BB genotypes
                ind <- tind; ind[!is.na(data) & data!=4] <- NA
                points(ind,y,pch=21,col="black", bg=color[4],cex=cex)

                # not AA genotypes
                ind <- tind; ind[!is.na(data) & data!=5] <- NA
                points(ind,y,pch=21,col="black", bg=color[5],cex=cex)
            }

            # plot map
            u <- par("usr")
            segments(u[1],map,(u[1]+1)/2,map)
            segments(u[2],map,(n.ind+u[2])/2,map)

            if(any(errors != 0)) {
                ind <- rep(1:n.ind,length(map));ind[errors!=1]<-NA
                points(ind,y,pch=0,col=color[6],cex=cex+0.4,lwd=2)
            }

            if(include.xo) points(xoloc$ind,xoloc$loc,pch=4,col="blue",lwd=2)

        }
    }
    invisible()
}

######################################################################
#
# plotInfo: Plot the proportion of missing information in the
#            genotype data.
#
######################################################################
plotInfo <- plot.info <-
    function(x,chr,method=c("entropy","variance","both"), step=1,
             off.end=0, error.prob=0.001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             alternate.chrid=FALSE, fourwaycross=c("all", "AB", "CD"),
             include.genofreq=FALSE, ...)
{
    cross <- x
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    method <- match(match.arg(method),c("entropy","variance","both"))-1
    map.function <- match.arg(map.function)

    if(!missing(chr)) cross <- subset(cross,chr=chr)

    n.chr <- nchr(cross)
    results <- NULL

    cross <- calc.genoprob(cross, step=step, error.prob=error.prob,
                           off.end=off.end, map.function=map.function)

    gap <- 25

    # 4-way cross: can consider just A/B or just C/D
    fourwaycross <- match.arg(fourwaycross)

    n.ind <- nind(cross)
    if(include.genofreq) theprob <- NULL
    for(i in 1:n.chr) {

        n.gen <- dim(cross$geno[[i]]$prob)[3]
        n.pos <- ncol(cross$geno[[i]]$prob)

        # 4-way cross: can consider just A/B or just C/D
        if(n.gen==4 && (fourwaycross!="all")) {
            att <- attributes(cross$geno[[i]]$prob)
            if(fourwaycross == "AB") {
                cross$geno[[i]]$prob[,,1] <- cross$geno[[i]]$prob[,,1] + cross$geno[[i]]$prob[,,3]
                cross$geno[[i]]$prob[,,2] <- cross$geno[[i]]$prob[,,2] + cross$geno[[i]]$prob[,,4]
            }
            else {
                cross$geno[[i]]$prob[,,1] <- cross$geno[[i]]$prob[,,1] + cross$geno[[i]]$prob[,,2]
                cross$geno[[i]]$prob[,,2] <- cross$geno[[i]]$prob[,,3] + cross$geno[[i]]$prob[,,4]
            }
            cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,,1:2]
            n.gen <- 2

            for(j in seq(along=att))
                if(names(att)[j] != "dim" && names(att)[j] != "dimnames")
                    attr(cross$geno[[i]]$prob, names(att)[j]) <- att[[j]]
        }

        # calculate information (between 0 and 1)
        info <- .C("R_info",
                   as.integer(n.ind),
                   as.integer(n.pos),
                   as.integer(n.gen),
                   as.double(cross$geno[[i]]$prob),
                   info1=as.double(rep(0,n.pos)),
                   info2=as.double(rep(0,n.pos)),
                   as.integer(method),
                   PACKAGE="qtl")

        if(method != 1) { # rescale entropy version
            if(n.gen==3) maxent <- 1.5*log(2)
            else maxent <- log(n.gen)
            info$info1 <- -info$info1/maxent
        }
        if(method != 0) { # rescale variance version
            maxvar <- c(0.25,0.5,1.25)[n.gen-1]
            info$info2 <- info$info2/maxvar
        }

        # reconstruct map
        if("map" %in% names(attributes(cross$geno[[i]]$prob)))
            map <- attr(cross$geno[[i]]$prob,"map")
        else {
            stp <- attr(cross$geno[[i]]$prob, "step")
            oe <- attr(cross$geno[[i]]$prob, "off.end")

            if("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
                stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
            else stpw <- "fixed"
            map <- create.map(cross$geno[[i]]$map,stp,oe,stpw)
        }

        if(is.matrix(map)) map <- map[1,]

        z <- data.frame(chr=rep(names(cross$geno)[i],length(map)),
                        pos=as.numeric(map),
                        "Missing information"=info$info1,
                        "Missing information"=info$info2, stringsAsFactors=TRUE)

        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
        rownames(z) <- w
        results <- rbind(results, z)

        if(include.genofreq) {
            p <- cross$geno[[i]]$prob
            if(class(cross$geno[[i]])=="X")
                p <- reviseXdata(class(cross)[1], expandX="full", sexpgm=getsex(cross), prob=p,
                                 cross.attr=attributes(cross))
            p <- apply(p, 2:3, mean, na.rm=TRUE)
            if(is.null(theprob))
                theprob <- p
            else {
                m1 <- match(colnames(p), colnames(theprob))
                m2 <- match(colnames(theprob), colnames(p))
                if(!any(is.na(m1)) && !any(is.na(m2))) {
                    theprob <- rbind(theprob, p[,colnames(theprob)])
                }
                else { # some differences in column names
                    if(all(!is.na(m1))) { # need to add some new columns to p
                        noldcol <- ncol(p)
                        nnewcol <- sum(is.na(m2))
                        for(i in 1:nnewcol)
                            p <- cbind(p, 0)
                        colnames(p)[-(1:noldcol)] <- colnames(theprob)[is.na(m2)]
                        theprob <- rbind(theprob, p[,colnames(theprob)])
                    }
                    else if(all(!is.na(m2))) { # need to add some new columns to theprob
                        noldcol <- ncol(theprob)
                        nnewcol <- sum(is.na(m1))
                        for(i in 1:nnewcol)
                            theprob <- cbind(theprob, 0)
                        colnames(theprob)[-(1:noldcol)] <- colnames(p)[is.na(m1)]
                        theprob <- rbind(theprob, p[,colnames(theprob)])
                    }
                    else { # need to add some new columns to each
                        noldcol <- ncol(theprob)
                        nnewcol <- sum(is.na(m1))
                        oldcolnam <- colnames(theprob)
                        for(i in 1:nnewcol)
                            theprob <- cbind(theprob, 0)
                        colnames(theprob)[-(1:noldcol)] <- colnames(p)[is.na(m1)]

                        noldcol <- ncol(p)
                        nnewcol <- sum(is.na(m2))
                        for(i in 1:nnewcol)
                            p <- cbind(p, 0)
                        colnames(p)[-(1:noldcol)] <- oldcolnam[is.na(m2)]

                        theprob <- rbind(theprob, p[,colnames(theprob)])
                    }
                }
            }
        } # end of if(include.genofreq)

    }
    colnames(results)[3:4] <- c("misinfo.entropy","misinfo.variance")

    if(method==0) results <- results[,-4]
    if(method==1) results <- results[,-3]

    if(include.genofreq)
        results <- cbind(results, theprob)

    class(results) <- c("scanone","data.frame")

    # check whether main was included as an argument
    args <- list(...)
    if("main" %in% names(args)) hasmain <- TRUE
    else hasmain <- FALSE

    # check whether gap was included as an argument
    if(!("gap" %in% names(args))) {
        if(method==0) {
            if(hasmain)
                plot.scanone(results,ylim=c(0,1),gap=gap,
                             alternate.chrid=alternate.chrid,...)
            else
                plot.scanone(results,ylim=c(0,1),gap=gap,
                             main="Missing information",
                             alternate.chrid=alternate.chrid,...)
        }
        else if(method==1) {
            if(hasmain)
                plot.scanone(results,ylim=c(0,1),gap=gap,
                             alternate.chrid=alternate.chrid,...)
            else
                plot.scanone(results,ylim=c(0,1),gap=gap,
                             main="Missing information",
                             alternate.chrid=alternate.chrid,...)
        }
        else if(method==2) {
            if(hasmain)
                plot.scanone(results,results,lodcolumn=1:2,ylim=c(0,1),gap=gap,
                             alternate.chrid=alternate.chrid,...)
            else
                plot.scanone(results,results,lodcolumn=1:2,ylim=c(0,1),gap=gap,
                             main="Missing information",
                             alternate.chrid=alternate.chrid,...)
        }
    }
    else { # gap was included in ...
        if(method==0) {
            if(hasmain)
                plot.scanone(results,ylim=c(0,1),
                             alternate.chrid=alternate.chrid,...)
            else
                plot.scanone(results,ylim=c(0,1),
                             main="Missing information",
                             alternate.chrid=alternate.chrid,...)
        }
        else if(method==1) {
            if(hasmain)
                plot.scanone(results,ylim=c(0,1),
                             alternate.chrid=alternate.chrid,...)
            else
                plot.scanone(results,ylim=c(0,1),
                             main="Missing information",
                             alternate.chrid=alternate.chrid,...)
        }
        else if(method==2) {
            if(hasmain)
                plot.scanone(results,results,lodcolumn=1:2,ylim=c(0,1),
                             alternate.chrid=alternate.chrid,...)
            else
                plot.scanone(results,results,lodcolumn=1:2,ylim=c(0,1),
                             main="Missing information",
                             alternate.chrid=alternate.chrid,...)
        }
    }

    invisible(results)
}


# plot phenotypes against one or more markers
plotPXG <- plot.pxg <-
    function(x, marker, pheno.col = 1, jitter = 1, infer = TRUE,
             pch, ylab, main, col, ...)
{
    cross <- x
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    type <- class(cross)[1]

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("plotPXG can take just one phenotype; only the first will be used")
    }

    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }

    if(pheno.col < 1 | pheno.col > nphe(cross))
        stop("pheno.col values should be between 1 and the no. phenotypes")

    if(!is.numeric(cross$pheno[,pheno.col]))
        stop("phenotype \"", colnames(cross$pheno)[pheno.col], "\" is not numeric.")

    if(missing(pch)) pch <- par("pch")
    if(missing(ylab)) ylab <-  colnames(cross$pheno)[pheno.col]

    oldlas <- par("las")
    on.exit(par(las = oldlas))
    par(las = 1)

    # find chromosomes containing the markers
    o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data),
                marker)
    if(length(marker)==1) o <- matrix(o,nrow=1)
    if(!all(apply(o,1,any))) {
        oo <- apply(o,1,any)
        stop("Marker ", marker[!oo], " not found")
    }
    n.mark <- length(marker)
    o <- apply(o, 1, which)
    chr <- names(cross$geno)[o]
    uchr <- unique(chr)

    cross <- subset(cross, chr=uchr)
    map <- pull.map(cross)
    pos <- NULL
    for(i in seq(length(chr))) pos[i] <- map[[chr[i]]][marker[i]]
    chrtype <- sapply(cross$geno, class)
    names(chrtype) <- names(cross$geno)
    chrtype <- chrtype[chr]

    # if X chromosome and backcross or intercross, get sex/direction data
    if(any(chrtype == "X") && (type == "bc" || type == "f2"))
        sexpgm <- getsex(cross)
    else sexpgm <- NULL

    # number of possible genotypes
    gen.names <- list()
    for(i in seq(length(chr)))
        gen.names[[i]] <- getgenonames(type, chrtype[i], "full", sexpgm, attributes(cross))
    n.gen <- sapply(gen.names, length)

    jitter <- jitter/10
    if(any(n.gen == 2)) jitter <- jitter * 0.75

    # function to determine whether genotype is fully known
    tempf <-
        function(x, type)
        {
            tmp <- is.na(x)
            if(type=="f2") tmp[!is.na(x) & x>3] <- TRUE
            if(type=="4way") tmp[!is.na(x) & x>4] <- TRUE
            tmp
        }

    # if infer=TRUE, fill in genotype data by a single imputation
    if(infer) {
        which.missing <- tempf(cross$geno[[chr[1]]]$data[, marker[1]],type)
        if(n.mark > 1)
            for(i in 2:n.mark)
                which.missing <- which.missing | tempf(cross$geno[[chr[i]]]$data[,marker[i]],type)
        which.missing <- as.numeric(which.missing)

        cross <- fill.geno(cross, method = "imp")
    }
    else which.missing <- rep(1,nind(cross))

    # data to plot
    x <- cross$geno[[chr[1]]]$data[, marker[1]]
    if(n.mark > 1)
        for(i in 2:n.mark)
            x <- cbind(x, cross$geno[[chr[i]]]$data[, marker[i]])
    else x <- as.matrix(x)
    y <- cross$pheno[, pheno.col]

    if(!infer) { # replace partially informative genotypes with NAs
        if(type == "f2") x[x > 3] <- NA
        if(type == "4way") x[x > 4] <- NA
        if(sum(!is.na(x)) == 0)
            stop("Can't use infer=FALSE as there are no fully informative genotypes")
    }

    # in case of X chromosome, recode some genotypes
    if(any(chrtype == "X") && (type == "bc" || type == "f2")) {
        ix = seq(n.mark)[chrtype == "X"]
        for(i in ix)
            x[, i] <- as.numeric(reviseXdata(type, "full", sexpgm,
                                             geno = as.matrix(x[, i]),
                                             cross.attr=attributes(cross)))
    }

    # save all of the data, returned invisibly
    data <- as.data.frame(x, stringsAsFactors=TRUE)
    names(data) <- marker
    for(i in marker) data[[i]] <- ordered(data[[i]])
    data$pheno <- y
    data$inferred <- which.missing

    # re-code the multi-marker genotypes
    if(n.mark > 1) {
        for(i in 2:n.mark)
            x[, 1] <- n.gen[i] * (x[, 1] - 1) + x[, i]
    }
    x <- x[, 1]

    observed <- sort(unique(x))

    # amount of jitter
    u <- runif(nind(cross), -jitter, jitter)
    r <- (1 - 2 * jitter)/2

    # genotype names
    if(n.mark == 1)
        gnames <- gen.names[[1]]
    else {
        gnames <- array(gen.names[[n.mark]], c(prod(n.gen), n.mark))
        for(i in (n.mark - 1):1) {
            tmpi <- rep(gen.names[[i]], rep(prod(n.gen[(i + 1):n.mark]),
                                            n.gen[i]))
            if(i > 1)
                tmpi <- rep(tmpi, prod(n.gen[1:(i - 1)]))
            gnames[, i] <- tmpi
        }
        gnames <- apply(gnames, 1, function(x) paste(x, collapse = "\n"))
    }

    # create plot
    plot(x + u, y, xlab = "Genotype", ylab = ylab, type = "n",
         main = "", xlim = c(1 - r + jitter, length(gnames) + r + jitter), xaxt = "n")

    # marker names at top
    if(missing(main))
        mtext(paste(marker, collapse = "\n"),
              line=0.5, cex = ifelse(n.mark==1, 1.2, 0.8))
    else
        title(main=main)

    abline(v = 1:prod(n.gen), col = "gray", lty = 3)

    if(length(pch) == 1)
        pch = rep(pch, length(x))
    if(infer) {
        points((x + u)[which.missing == 1], y[which.missing ==
                       1], col = "red", pch = pch[which.missing == 1])
        points((x + u)[which.missing == 0], y[which.missing ==
                       0], pch = pch[which.missing == 0])
    }
    else points(x + u, y, pch = pch)
    sux = sort(unique(x))

    # add confidence intervals
    me <- se <- array(NA, length(gnames))
    me[sux] <- tapply(y, x, mean, na.rm = TRUE)
    se[sux] <- tapply(y, x, function(a) sd(a, na.rm = TRUE)/sqrt(sum(!is.na(a))))
    thecolors <- c("black", "blue", "red", "purple", "green", "orange")
    if(missing(col)) {
        col <- thecolors[1:n.gen[n.mark]]
        if(n.gen[n.mark] == 3)
            col <- c("blue", "purple", "red")
        else if(n.gen[n.mark] == 2)
            col <- c("blue", "red")
    }

    ng <- length(gnames)
    segments(seq(ng) + jitter * 2, me, seq(ng) +
             jitter * 4, me, lwd = 2, col = col)
    segments(seq(ng) + jitter * 3, me - se, seq(ng) +
             jitter * 3, me + se, lwd = 2, col = col)
    segments(seq(ng) + jitter * 2.5, me - se, seq(ng) +
             jitter * 3.5, me - se, lwd = 2, col = col)
    segments(seq(ng) + jitter * 2.5, me + se, seq(ng) +
             jitter * 3.5, me + se, lwd = 2, col = col)

    # add genotypes below
    u <- par("usr")
    cxaxis <- par("cex.axis")

    segments(seq(along=gnames), u[3], seq(along=gnames), u[3] - diff(u[3:4]) *
             0.015, xpd = TRUE)
    axis(side=1, at=seq(along=gnames), labels=gnames,
         cex=ifelse(n.mark==1, cxaxis, cxaxis*0.8),
         tick=FALSE, line = (length(marker)-1)/2)

    invisible(data)
}

plotPheno <- plot.pheno <-
    function(x, pheno.col=1, ...)
{
    if(!any(class(x) == "cross"))
        stop("Input should have class \"cross\".")

    if(LikePheVector(pheno.col, nind(x), nphe(x))) {
        x$pheno <- cbind(pheno.col, x$pheno)
        pheno.col <- 1
    }

    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("Ignoring all but the first element in pheno.col.")
    }
    if(is.character(pheno.col)) {
        num <- find.pheno(x, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }

    if(pheno.col < 1 | pheno.col > nphe(x))
        warning("pheno.col should be between 1 and ", nphe(x))

    phe <- x$pheno[,pheno.col]
    u <- length(unique(phe))
    if(u==2 || (u < 10  && nind(x) > 50))
        phe <- as.factor(phe)

    plot_pheno_sub <-
        function(phe, xlab=paste("phe", pheno.col),
                 main=colnames(x$pheno)[pheno.col], col="white",
                 breaks=ceiling(2*sqrt(nind(x))),
                 las=1,
                 ...)
        {
            if(is.factor(phe)) {
                barplot(table(phe), xlab=xlab, main=main, col=col, las=las, ...)
            }
            else {
                phe <- as.numeric(phe)[1:nind(x)]
                hist(phe, breaks = breaks,
                     xlab = xlab, main = main, las=las, ...)
            }
        }

    plot_pheno_sub(phe, ...)
}

# end of plot.R
