######################################################################
#
# est.rf.R
#
# copyright (c) 2001-2016, Karl W Broman
# last modified Jan, 2016
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
# Contains: est.rf, plotRF, checkAlleles, pull.rf, plot.rfmatrix
#
######################################################################

######################################################################
#
# est.rf: Estimate sex-averaged recombination fractions between
#         all pairs of markers
#
######################################################################

est.rf <-
    function(cross, maxit=10000, tol=1e-6)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    n.chr <- nchr(cross)
    n.mar <- totmar(cross)
    n.ind <- nind(cross)
    mar.names <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

    type <- class(cross)[1]
    chrtype <- sapply(cross$geno,class)

    is.bcsft <- (type == "bcsft")
    if(is.bcsft) {
        cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
        is.bcsft <- cross.scheme[2] > 0 ## used for fixX only
    }

    xchrcol <- NULL
    fixX <- FALSE
    Geno <- NULL
    # create full genotype matrix
    for(i in 1:n.chr) {
        temp <- cross$geno[[i]]$data

        # treat X chromosome specially in an intercross or BCsFt with t>0.
        if((type=="f2" || is.bcsft) && chrtype[i]=="X") {
            fixX <- TRUE
            if(i != 1) xchrcol <- c(xchrcol,ncol(Geno)+(1:ncol(cross$geno[[i]]$data)))
            else xchrcol <- 1:ncol(cross$geno[[i]]$data)
            xchr <- temp
            xchr[is.na(xchr)] <- 0
            temp <- reviseXdata("f2","simple",getsex(cross),geno=temp,
                                cross.attr=attributes(cross))
        }
        Geno <- cbind(Geno,temp)
    }

    # which type of cross is this?
    if(type == "f2")
        cfunc <- "est_rf_f2"
    else if(type == "bc" || type=="risib" || type=="riself" || type=="dh" || type=="haploid")
        cfunc <- "est_rf_bc"
    else if(type == "4way")
        cfunc <- "est_rf_4way"
    else if(type=="ri8sib" || type=="ri8self" || type=="ri4sib" || type=="ri4self") {
        cfunc <- paste("est_rf_", type, sep="")
        if(any(chrtype == "X"))
            warning("est.rf not working properly for the X chromosome for 4- or 8-way RIL.")
    }
    else if(type == "bcsft")
        cfunc <- "est_rf_bcsft"
    else
        stop("est.rf not available for cross type ", type, ".")

    Geno[is.na(Geno)] <- 0

    if(type=="bc" || type=="risib" || type=="riself" || type=="dh" || type=="haploid")
        z <- .C(cfunc,
                as.integer(n.ind),         # number of individuals
                as.integer(n.mar),         # number of markers
                as.integer(Geno),
                rf = as.double(rep(0,n.mar*n.mar)),
                PACKAGE="qtl")
    else {
        ## Hide cross scheme in genoprob to pass to routine. BY
        temp <- as.double(rep(0,n.mar*n.mar))
        if(type == "bcsft")
            temp[1:2] <- cross.scheme

        z <- .C(cfunc,
                as.integer(n.ind),         # number of individuals
                as.integer(n.mar),         # number of markers
                as.integer(Geno),
                rf = as.double(temp),
                as.integer(maxit),
                as.double(tol),
                PACKAGE="qtl")
    }

    cross$rf <- matrix(z$rf,ncol=n.mar)
    dimnames(cross$rf) <- list(mar.names,mar.names)

    if(fixX) {
        temp <- as.double(rep(0, ncol(xchr) ^ 2))
        if(type == "bcsft") {
            temp[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)

            zz <- .C(cfunc,
                     as.integer(n.ind),         # number of individuals
                     as.integer(ncol(xchr)),   # number of markers on X chr.
                     as.integer(xchr),
                     rf = as.double(temp),
                     as.integer(maxit),
                     as.double(tol),
                     PACKAGE="qtl")
        }
        else {
            zz <- .C("est_rf_bc",
                     as.integer(n.ind),
                     as.integer(ncol(xchr)),
                     as.integer(xchr),
                     rf=as.double(temp),
                     PACKAGE="qtl")
        }
        zz <- matrix(zz$rf,ncol=ncol(xchr))
        cross$rf[xchrcol,xchrcol] <- zz
    }

    # check for alleles switches
    if(type == "risib" || type=="riself" || type=="f2" || type=="bc" || type=="dh" || type=="haploid") {
        out <- checkAlleles(cross, 5, FALSE)
        if(!is.null(out)) {
            out <- as.character(out[,1])
            warning("Alleles potentially switched at markers \n  ",
                    paste(out, collapse=" "))
        }
    }

    cross
}



plotRF <- plot.rf <-
    function(x, chr, what=c("both","lod","rf"),
             alternate.chrid=FALSE, zmax=12,
             mark.diagonal=FALSE,
             col.scheme=c("viridis", "redblue"),
             ...)
{
    if(!any(class(x) == "cross"))
        stop("Input should have class \"cross\".")

    what <- match.arg(what)
    if("onlylod" %in% names(attributes(x$rf)) && attr(x$rf, "onlylod")) {
        onlylod <- TRUE
        what <- "lod"
    }
    else onlylod <- FALSE

    if(!missing(chr)) x <- subset(x,chr=chr)

    if(!("rf" %in% names(x))) {
        warning("Running est.rf.")
        x <- est.rf(x)
    }
    g <- x$rf

    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd=TRUE,las=1)
    on.exit(par(xpd=old.xpd,las=old.las))

    if(!onlylod) {
        # if any of the rf's are NA (ie no data), put NAs in corresponding LODs
        if(any(is.na(g))) g[is.na(t(g))] <- NA

        # convert rf to -2*(log2(rf)+1); place zmax's on the diagonal;
        #    anything above zmax replaced by zmax;
        #    NA's replaced by -1
        g[row(g) > col(g) & g > 0.5] <- 0.5
        g[row(g) > col(g)] <- -4*(log2(g[row(g) > col(g)])+1)/12*zmax
    }
    diag(g) <- zmax
    g[!is.na(g) & g>zmax] <- zmax

    g[is.na(g)] <- -1

    if(what=="lod" && !onlylod) { # plot LOD scores
        # copy upper triangle (LODs) to lower triangle (rec fracs)
        g[row(g) > col(g)] <- t(g)[row(g) > col(g)]
    }
    else if(what=="rf") { # plot recombination fractions
        # copy lower triangle (rec fracs) to upper triangle (LODs)
        g[row(g) < col(g)] <- t(g)[row(g) < col(g)]
    }
    br <- c(-1, seq(-1e-6, zmax, length=257))

    col.scheme <- match.arg(col.scheme)
    if(col.scheme=="redblue") {
        # convert colors using gamma=0.6 (which will no longer be available in R)
        thecol <- rev(rainbow(256, start=0, end=2/3))
        rgbval <- (col2rgb(thecol)/255)^0.6
        thecol <- rgb(rgbval[1,], rgbval[2,], rgbval[3,])
    } else {
        thecol <- viridis_qtl(256) # the new default
    }

    image(1:ncol(g),1:nrow(g),t(g),ylab="Markers",xlab="Markers",breaks=br,
          col=c("lightgray",thecol))

    if(mark.diagonal) {
        for(i in 1:ncol(g))
            segments(i+c(-0.5, -0.5, -0.5, +0.5), i+c(-0.5, +0.5, -0.5, -0.5),
                     i+c(-0.5, +0.5, +0.5, +0.5), i+c(+0.5, +0.5, -0.5, +0.5))
    }

    # plot lines at the chromosome boundaries
    n.mar <- nmar(x)
    n.chr <- nchr(x)
    a <- c(0.5,cumsum(n.mar)+0.5)
    abline(v=a,xpd=FALSE,col="white")
    abline(h=a,xpd=FALSE,col="white")

    # this line adds a line above the image
    #     (the image function leaves it out)
    abline(h=0.5+c(0,nrow(g)),xpd=FALSE)
    abline(v=0.5+c(0,nrow(g)),xpd=FALSE)

    # add chromosome numbers
    a <- par("usr")
    wh <- cumsum(c(0.5,n.mar))
    chrnam <- names(x$geno)
    chrpos <- (wh[-1] + wh[-length(wh)])/2
    if(!alternate.chrid || length(chrnam) < 2) {
        for(i in seq(along=chrpos)) {
            axis(side=3, at=chrpos[i], labels=chrnam[i], tick=FALSE, line=-0.8)
            axis(side=4, at=chrpos[i], labels=chrnam[i], tick=FALSE, line=-0.8)
        }
    }
    else {
        odd <- seq(1, length(chrpos), by=2)
        even <- seq(2, length(chrpos), by=2)
        for(i in odd) {
            axis(side=3, at=chrpos[i], labels=chrnam[i], line=-0.8, tick=FALSE)
            axis(side=4, at=chrpos[i], labels=chrnam[i], line=-0.8, tick=FALSE)
        }
        for(i in even) {
            axis(side=3, at=chrpos[i], labels=chrnam[i], line=0, tick=FALSE)
            axis(side=4, at=chrpos[i], labels=chrnam[i], line=0, tick=FALSE)
        }

    }


    dots <- list(...)
    if("main" %in% names(dots))
        title(main=dots$main)
    else {
        if(what=="lod") title(main="Pairwise LOD scores")
        else if(what=="rf") title(main="Recombination fractions")
        else title("Pairwise recombination fractions and LOD scores")
    }

    invisible()
}

######################################################################
# check for apparent errors in the recombination fractions
######################################################################
#checkrf <-
#function(cross, threshold=5)
#{
#  rf <- cross$rf
#  n.mar <- nmar(cross)
#  map <- pull.map(cross)
#  n <- ncol(rf)
#  mnam <- colnames(rf)
#  whpos <- unlist(lapply(map,function(a) 1:length(a)))
#  whchr <- rep(names(map),sapply(map,length))
#
#  # first check whether a locus has "significant" pairwise recombination
#  #     with rf > 0.5
#  for(i in 1:n) {
#    if(i == 1) {
#      lod <- rf[1,-1]
#      r <- rf[-1,1]
#    }
#    else if(i == n) {
#      lod <- rf[-n,n]
#      r <- rf[n,-n]
#    }
#    else {
#      lod <- c(rf[1:(i-1),i],rf[i,(i+1):n])
#      r <- c(rf[i,1:(i-1)],rf[(i+1):n,i])
#    }
#
#    # if rf > 1/2 and LOD > threshold for more than two other markers
#    if(sum(!is.na(lod) & !is.na(r) & lod > threshold & r > 0.5) >= 2)
#      warning("Genotypes potentially switched for marker ", mnam[i],
#          paste(" (",whpos[i],")",sep=""), " on chr ", whchr[i], "\n")
#
#  }
#
#}


######################################################################
# checkAlleles()
#
# Function to find markers that may have alleles miscoded;
# we go through each marker, one at a time, swap alleles and
# then see what it does to pairwise linkage against all other
# markers
######################################################################

checkAlleles <-
    function(cross, threshold=3, verbose=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("checkAlleles() only works for cross objects.")

    type <- class(cross)[1]
    if(type != "f2" && type != "bc" &&
       type != "risib" && type != "riself" && type != "dh" && type!="haploid")
        stop("checkAlleles not available for cross type ", type, ".")

    # drop X chromosome
    chrtype <- sapply(cross$geno,class)
    if(all(chrtype=="X")) {
        if(verbose) cat("checkAlleles() only works for autosomal data.\n")
        return(NULL)
    }

    cross <- subset(cross, chr = (chrtype != "X"))

    n.mar <- nmar(cross)
    mar.names <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

    if(!("rf" %in% names(cross))) {
        warning("First running est.rf.")
        cross <- est.rf(cross)
    }
    diag(cross$rf) <- 0
    lod <- rf <- cross$rf
    lod[lower.tri(lod)] <- t(lod)[lower.tri(lod)]
    rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]

    orig.lod <- rev.lod <- lod
    orig.lod[rf > 0.5] <- 0
    rev.lod[rf < 0.5] <- 0

    dif <- apply(rev.lod, 2, max, na.rm=TRUE) -
        apply(orig.lod, 2, max, na.rm=TRUE)

    results <- data.frame(marker=mar.names,
                          chr=rep(names(cross$geno), n.mar),
                          index=unlist(lapply(n.mar, function(a) 1:a)),
                          "diff in max LOD" = dif, stringsAsFactors=TRUE)
    rownames(results) <- 1:nrow(results)

    if(all(results[,4] < threshold)) {
        if(verbose) cat("No apparent problems.\n")
        return(invisible(NULL))
    }

    results[results[,4] >= threshold,]
}

######################################################################
# pull.rf
#
# pull out the pairwise marker recombination fraction estimates,
# with class "rfmatrix"
######################################################################

pull.rf <-
    function(cross, what=c("rf", "lod"), chr)
{
    if(!("cross" %in% class(cross)))
        stop("Input must have class \"cross\".")

    if(!missing(chr)) cross <- subset(cross, chr=chr)

    if(!("rf" %in% names(cross))) {
        warning(" -Running est.rf")
        cross <- est.rf(cross)
    }
    rf <- cross$rf
    if(nrow(rf) != totmar(cross) || ncol(rf) != totmar(cross) ||
       any(rownames(rf) != markernames(cross)) ||
       any(colnames(rf) != markernames(cross))) {
        warning(" -Rec. frac. estimates seem corrupted; re-running est.rf")
        cross <- est.rf(cross)
    }
    rf <- cross$rf

    diag(rf) <- NA
    what <- match.arg(what)
    if(what=="rf")
        rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
    else
        rf[lower.tri(rf)] <- t(rf)[lower.tri(rf)]

    attr(rf, "map") <- pull.map(cross, as.table=TRUE)
    attr(rf, "what") <- what
    class(rf) <- c("rfmatrix", "matrix")
    rf
}

######################################################################
# plot.rfmatrix:
#
# plot a slice through the matrix (that is, for one marker)
######################################################################

plot.rfmatrix <-
    function(x, marker, ...)
{
    if(!("rfmatrix" %in% class(x)))
        stop("Input must have class \"rfmatrix\".")

    if(missing(marker))
        stop("You must provide a marker name.")

    if(length(marker) > 1) {
        warning("Ignoring all but the first marker, ", marker[1])
        marker <- marker[1]
    }

    if(!(marker %in% rownames(x)))
        stop("Marker ", marker, " not found.")

    what <- attr(x, "what")
    x <- cbind(attr(x, "map"), x[marker,])
    x$chr <- factor(x$chr, levels=unique(x$chr))

    # fill in hole
    wh <- which(rownames(x)==marker)
    if(wh > 1 && wh < nrow(x) && x[wh-1,1] == x[wh,1] && x[wh+1,1] == x[wh,1]) {
        xl <- x[wh-1,2]
        xm <- x[wh,2]
        xr <- x[wh+1,2]
        yl <- x[wh-1,3]
        yr <- x[wh+1,3]
        x[wh,3] <- yl + (xm-xl)*(yr-yl)/(xr-xl)
    }

    colnames(x)[3] <- what
    class(x) <- c("scanone", "data.frame")

    dots <- list(...)
    if("main" %in% names(dots))
        plot(x, ...)
    else
        plot(x, main=marker, ...)

    invisible(x)
}


# end of est.rf.R
