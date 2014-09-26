######################################################################
#
# errorlod.R
#
# copyright (c) 2001-2013, Karl W Broman
# last modified Sep, 2013
# first written Apr, 2001
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
# Contains: calc.errorlod, plotErrorlod, top.errorlod
#
######################################################################

######################################################################
#
# calc.errorlod: Calculate LOD scores indicating likely genotyping
#                errors.
#
######################################################################

calc.errorlod <-
    function(cross, error.prob=0.01,
             map.function=c("haldane","kosambi","c-f","morgan"),
             version=c("new","old"))
{
    version <- match.arg(version)

    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    origcross <- cross

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 0.5) {
        error.prob <- 0.5
        warning("error.prob shouldn't be > 0.5.")
    }

    # map function
    map.function <- match.arg(map.function)

    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    type <- class(cross)[1]

    # calculate genotype probabilities one chromosome at a time
    for(i in 1:n.chr) {

        chr.type <- class(cross$geno[[i]])
        if(type=="bc" || type=="risib" || type=="riself" || type=="dh" || type=="haploid")
            cfunc <- "calc_errorlod_bc"
        else if(type=="f2" || type=="bcsft") {
            if(chr.type!="X") cfunc <- "calc_errorlod_f2"
            else cfunc <- "calc_errorlod_bc"
        }
        else if(type=="4way") cfunc <- "calc_errorlod_4way"
        else if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="bgmagic16")
            cfunc <- paste("calc_errorlod_", type, sep="")
        else
            stop("calc.errorlod not available for cross type ", type, ".")

        # skip chromosomes with only 1 marker
        if(n.mar[i] < 2) next

        if(version=="old") {
            if((!("prob" %in% names(cross$geno[[i]])) ||
                abs(attr(cross$geno[[i]]$prob,"error.prob")
                    - error.prob) > 1e-9)) {
                # need to run calc.genoprob
                cross <- calc.genoprob(cross,error.prob=error.prob,
                                       map.function=map.function)
            }

            Pr <- cross$geno[[i]]$prob
            u <- grep("^loc-*[0-9]+",colnames(Pr))
            if(length(u) > 0) Pr <- Pr[,-u,]
        }
        else { # new version
            cross <- calc.genoprob.special(cross,error.prob=error.prob,
                                           map.function=map.function)
            Pr <- cross$geno[[i]]$prob
        }

        nm <- dim(Pr)[2]
        dat <- cross$geno[[i]]$data
        dat[is.na(dat)] <- 0

        z <- .C(cfunc,
                as.integer(n.ind),
                as.integer(nm),
                as.integer(dat),
                as.double(error.prob),
                as.double(Pr),
                errlod=as.double(rep(0,n.ind*nm)),
                PACKAGE="qtl")

        errlod <- array(z$errlod, dim=dim(Pr)[1:2])
        if(version=="new") errlod[dat==0] <- -99

        dimnames(errlod) <- list(NULL,colnames(cross$geno[[i]]$data))
        origcross$geno[[i]]$errorlod <- errlod

        # attribute set to the error.prob value used, for later
        #     reference.
        attr(origcross$geno[[i]]$errorlod,"error.prob") <- error.prob
        attr(origcross$geno[[i]]$errorlod,"map.function") <- map.function
    }

    origcross
}





######################################################################
#
# plotErrorlod
#
######################################################################

plotErrorlod <- plot.errorlod <-
    function(x, chr, ind, breaks=c(-Inf,2,3,4.5,Inf),
             col=c("white","gray85","hotpink","purple3"),
             alternate.chrid=FALSE, ...)
{
    if(!any(class(x) == "cross"))
        stop("Input should have class \"cross\".")

    if(length(breaks) != length(col)+1)
        stop("Length of breaks should be length(col)+1.")
    if(length(breaks) != length(col)+1)
        col <- col[1:(length(breaks)+1)]

    cross <- x
    if(!missing(chr)) cross <- subset(cross,chr=chr)

    use.id <- FALSE
    if(!missing(ind)) {
        if(is.null(getid(cross))) cross$pheno$id <- 1:nind(cross)
        cross <- subset(cross,ind=ind)
        use.id <- TRUE
    }

    # remove chromosomes with < 2 markers
    n.mar <- nmar(cross)
    cross <- subset(cross,chr=names(n.mar)[n.mar >= 2])
    n.chr <- nchr(cross)

    errlod <- NULL
    for(i in 1:n.chr) {
        if(!("errorlod" %in% names(cross$geno[[i]]))) { # need to run calc.errorlod
            warning("First running calc.errorlod.")
            cross <- calc.errorlod(cross,error.prob=0.01,map.function="haldane")
        }
        errlod <- cbind(errlod,cross$geno[[i]]$errorlod)
    }

    errlod <- t(errlod)

    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd=TRUE,las=1)
    on.exit(par(xpd=old.xpd,las=old.las))

    # plot grid
    breaks[breaks==Inf] <- max(errlod)
    breaks[breaks==-Inf] <- min(errlod)
    image(1:nrow(errlod),1:ncol(errlod),errlod,
          ylab="Individuals",xlab="Markers",col=col,
          breaks=breaks, yaxt="n")
    if(use.id)
        axis(side=2, at=1:nind(cross), labels=getid(cross))
    else axis(side=2)

    # plot lines at the chromosome boundaries
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)
    chr.names <- names(cross$geno)
    a <- c(0.5,cumsum(n.mar)+0.5)

    # the following makes the lines go slightly above the plotting region
    b <- par("usr")
    segments(a,b[3],a,b[4]+diff(b[3:4])*0.02)

    # this line adds a line above and below the image
    #     (the image function seems to leave these out)
    abline(h=0.5+c(0,ncol(errlod)),xpd=FALSE)

    # add chromosome numbers
    a <- par("usr")
    wh <- cumsum(c(0.5,n.mar))
    chrpos <- (wh[-1]+wh[-length(wh)])/2

    if(!alternate.chrid || length(chrpos) < 2) {
        for(i in seq(along=chrpos))
            axis(side=3, at=chrpos[i], labels=chr.names[i])
    }
    else {
        odd <- seq(1, length(chrpos), by=2)
        even <- seq(2, length(chrpos), by=2)
        for(i in odd) {
            axis(side=3, at=chrpos[i], labels="")
            axis(side=3, at=chrpos[i], labels=chr.names[i], line=-0.4, tick=FALSE)
        }
        for(i in even) {
            axis(side=3, at=chrpos[i], labels="")
            axis(side=3, at=chrpos[i], labels=chr.names[i], line=+0.4, tick=FALSE)
        }
    }

    title(main="Genotyping error LOD scores")

    invisible()
}


######################################################################
#
# top.errorlod
#
# Picks out the genotypes having errorlod values above some cutoff
#
######################################################################

top.errorlod <-
    function(cross, chr, cutoff=4, msg=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(!missing(chr)) cross <- subset(cross,chr=chr)

    id <- getid(cross)
    if(is.null(id)) id <- 1:nind(cross)

    mar <- ind <- lod <- chr <- NULL

    # remove chromosomes with < 2 markers
    n.mar <- nmar(cross)
    cross <- subset(cross,chr=names(n.mar)[n.mar >= 2])

    flag <- 0
    for(i in 1:nchr(cross)) {

        if(!("errorlod" %in% names(cross$geno[[i]])))
            stop("You first need to run calc.errorlod.")

        el <- cross$geno[[i]]$errorlod

        if(any(el > cutoff)) {
            o <- (el > cutoff)
            mar <- c(mar,colnames(el)[col(el)[o]])
            ind <- c(ind,as.character(id[row(el)][o]))
            lod <- c(lod,el[o])
            chr <- c(chr,rep(names(cross$geno)[i],sum(o)))
            flag <- 1
        }
    }
    if(!flag) {
        if(msg) cat("\tNo errorlods above cutoff.\n")
        return(invisible(NULL))
    }

    suppressWarnings(asnum <- as.numeric(ind))
    if(!any(is.na(asnum)) && all(ind==asnum)) ind <- asnum

    o <- data.frame(chr=chr,id=ind,marker=mar,errorlod=lod,stringsAsFactors=FALSE)[order(-lod,ind),]
    rownames(o) <- 1:nrow(o)
    o
}



# end of errorlod.R
