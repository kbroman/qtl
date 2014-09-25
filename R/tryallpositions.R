######################################################################
#
# tryallpositions.R
#
# copyright (c) 2007-2013, Karl W Broman
# last modified Sep, 2013
# first written Oct, 2007
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
# Contains: tryallpositions, markerloglik, allchrsplits
#
######################################################################

######################################################################
# tryallpositions
#
# Place a given marker in all possible positions on a selected set of
# chromosomes, keeping the positions of all other markers fixed, and
# evaluate the likelihood and estimate the chromosome length
######################################################################

tryallpositions <-
    function(cross, marker, chr, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE, verbose=TRUE)
{
    map.function <- match.arg(map.function)
    if(missing(chr)) chr <- names(cross$geno)
    else chr <- matchchr(chr, names(cross$geno))

    thechr <- find.markerpos(cross, marker)[1,1]
    if(is.na(thechr))
        stop("Marker ", marker, " not found.")

    markerll <- markerloglik(cross, marker, error.prob)

    allchr <- names(subset(cross, chr=chr)$geno)

    results <- NULL
    for(i in allchr) {
        if(i == thechr) { # marker already on this chromosome
            pos <- cross$geno[[i]]$map
            if(is.matrix(pos)) {
                pos <- pos[1,]
                matrixmap <- TRUE
            }
            else matrixmap <- FALSE

            n.mar <- ncol(cross$geno[[i]]$data)
            if(n.mar == 1) { # this is the only marker
                pos <- 0
                length <- length.male <- 0
                llik <- 0
                if(verbose) cat(i, pos, llik/log(10), "\n")
                interval <- "---"
            } # just this marker
            else if(n.mar == 2) { # just two markers
                pos <- 0
                themar <- colnames(cross$geno[[i]]$data)
                initialloglik <- sum(sapply(themar, function(a, b, c) markerloglik(b, a, c), cross, error.prob))

                nm <- est.map(cross, chr=i, error.prob=error.prob,
                              map.function=map.function, m=m, p=p, maxit=maxit,
                              tol=tol, sex.sp=sex.sp, verbose=FALSE,
                              omit.noninformative=FALSE)[[1]]
                llik <- attr(nm, "loglik") - initialloglik
                if(verbose) cat(i, pos, llik/log(10), "\n")

                if(is.matrix(nm)) {
                    length <- diff(range(nm[1,]))
                    length.male <- diff(range(nm[2,]))
                }
                else length <- diff(range(nm))
                interval <- "---"
            } # just two markers
            else { # >2 markers
                temp <- drop.markers(cross, marker)

                pos <- temp$geno[[i]]$map
                if(is.matrix(pos)) {
                    pos <- pos[1,]
                    matrixmap <- TRUE
                }
                else matrixmap <- FALSE

                pos <- c(min(pos)-10, pos, max(pos)+10)
                pos <- (pos[-1] + pos[-length(pos)])/2

                themarkers <- colnames(temp$geno[[i]]$data)
                int2 <- c("*", match(themarkers, colnames(cross$geno[[i]]$data)), "*")
                themarkers <- c("pter", themarkers, "qter")
                interval <- paste(themarkers[-length(themarkers)], themarkers[-1],
                                  sep="-")
                int2 <- paste("(", int2[-length(int2)], "-", int2[-1], ")", sep="")
                interval <- paste(interval, int2)

                initialmap <- est.map(temp, chr=i,
                                      error.prob=error.prob,
                                      map.function=map.function, m=m, p=p, maxit=maxit,
                                      tol=tol, sex.sp=sex.sp, verbose=FALSE,
                                      omit.noninformative=FALSE)[[1]]
                initialloglik <- attr(initialmap, "loglik") + markerll

                llik <- length <- length.male <- rep(NA, length(pos))
                for(j in seq(along=pos)) {
                    temp <- movemarker(cross, marker, i, pos[j])
                    nm <- est.map(temp, chr=i, error.prob=error.prob,
                                  map.function=map.function, m=m, p=p, maxit=maxit,
                                  tol=tol, sex.sp=sex.sp, verbose=FALSE,
                                  omit.noninformative=FALSE)[[1]]
                    llik[j] <- attr(nm, "loglik") - initialloglik
                    if(verbose) cat(i, pos[j], llik[j]/log(10), "\n")

                    if(is.matrix(nm)) {
                        length[j] <- diff(range(nm[1,]))
                        length.male[j] <- diff(range(nm[2,]))
                    }
                    else length[j] <- diff(range(nm))
                }
            } # >2 markers

        }

        else { # marker not on this chromosome

            pos <- cross$geno[[i]]$map
            if(is.matrix(pos)) {
                pos <- pos[1,]
                matrixmap <- TRUE
            }
            else matrixmap <- FALSE

            n.mar <- ncol(cross$geno[[i]]$data)

            if(n.mar > 1) {
                pos <- c(min(pos)-10, pos, max(pos)+10)
                pos <- (pos[-1] + pos[-length(pos)])/2

                themarkers <- colnames(cross$geno[[i]]$data)
                int2 <- c("*", 1:length(themarkers), "*")
                themarkers <- c("pter", themarkers, "qter")
                interval <- paste(themarkers[-length(themarkers)], themarkers[-1],
                                  sep="-")
                int2 <- paste("(", int2[-length(int2)], "-", int2[-1], ")", sep="")
                interval <- paste(interval, int2)

                initialmap <- est.map(cross, chr=i, error.prob=error.prob,
                                      map.function=map.function, m=m, p=p, maxit=maxit,
                                      tol=tol, sex.sp=sex.sp, verbose=FALSE,
                                      omit.noninformative=FALSE)[[1]]
                initialloglik <- attr(initialmap, "loglik") + markerll

                llik <- length <- length.male <- rep(NA, length(pos))
                for(j in seq(along=pos)) {
                    temp <- movemarker(cross, marker, i, pos[j])

                    nm <- est.map(temp, chr=i, error.prob=error.prob,
                                  map.function=map.function, m=m, p=p, maxit=maxit,
                                  tol=tol, sex.sp=sex.sp, verbose=FALSE,
                                  omit.noninformative=FALSE)[[1]]
                    llik[j] <- attr(nm, "loglik") - initialloglik
                    if(verbose) cat(i, pos[j], llik[j]/log(10), "\n")

                    if(is.matrix(nm)) {
                        length[j] <- diff(range(nm[1,]))
                        length.male[j] <- diff(range(nm[2,]))
                    }
                    else length[j] <- diff(range(nm))
                }
            } # >1 marker on chromosome
            else {
                initialloglik <- markerloglik(cross, colnames(cross$geno[[i]]$data), error.prob) + markerll

                pos <- pos+10
                interval <- "---"

                temp <- movemarker(cross, marker, i, pos)

                nm <- est.map(temp, chr=i, error.prob=error.prob,
                              map.function=map.function, m=m, p=p, maxit=maxit,
                              tol=tol, sex.sp=sex.sp, verbose=FALSE,
                              omit.noninformative=FALSE)[[1]]
                llik <- attr(nm, "loglik") - initialloglik
                if(verbose) cat(i, pos, llik/log(10), "\n")

                if(is.matrix(nm)) {
                    length <- diff(range(nm[1,]))
                    length.male <- diff(range(nm[2,]))
                }
                else length <- diff(range(nm))

            } # one marker on chromosome
        } # marker not on this chr

        if(matrixmap && sex.sp)
            tempres <- data.frame(chr=rep(i, length(pos)),
                                  pos=pos,
                                  lod=llik/log(10),
                                  length.female=length,
                                  length.male=length.male,
                                  interval=interval,
                                  stringsAsFactors=FALSE)
        else
            tempres <- data.frame(chr=rep(i, length(pos)),
                                  pos=pos,
                                  lod=llik/log(10),
                                  length=length,
                                  interval=interval,
                                  stringsAsFactors=FALSE)


        results <- rbind(results, tempres)
    } # loop over chromosomes

    rownames(results) <- results$interval
    results <- results[,-ncol(results)]
    class(results) <- c("scanone", "data.frame")

    results
}

######################################################################
#
# markerloglik: Calculate log likelihood for a given marker
#
######################################################################

markerloglik <-
    function(cross, marker, error.prob=0.0001)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    type <- class(cross)[1]

    if(length(marker) > 1) {
        ll <- sapply(marker, function(a,b,d) markerloglik(b, a, d), cross, error.prob)
        names(ll) <- marker
        return(ll)
    }

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }

    n.ind <- nind(cross)

    thechr <- find.markerpos(cross, marker)[1,1]
    if(is.na(thechr))
        stop("Marker ", marker, " not found.")

    chrtype <- class(cross$geno[[thechr]])

    g <- pull.geno(cross, chr=thechr)
    m <- match(marker, colnames(g))
    if(is.na(m))
        stop("Marker ", marker, " not found.")
    g <- g[,m]
    g[is.na(g)] <- 0

    # which type of cross is this?
    if(type == "f2") {
        if(chrtype == "A") # autosomal
            cfunc <- "marker_loglik_f2"
        else                  # X chromsome
            cfunc <- "marker_loglik_bc"
    }
    else if(type == "bc" || type=="riself" || type=="risib" || type=="dh" || type=="haploid") {
        cfunc <- "marker_loglik_bc"
    }
    else if(type == "4way") {
        cfunc <- "marker_loglik_4way"
    }
    else if(type=="ri4sib" || type=="ri4self" || type=="ri8sib" || type=="ri8self" || type=="bgmagic16") {
        cfunc <- paste("marker_loglik_", type, sep="")
        if(chrtype=="X")
            warning("markerloglik not working properly for the X chromosome for 4- or 8-way RIL.")
    }
    else if(type == "bcsft") {
        cfunc <- "marker_loglik_bcsft"
        cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
        if(chrtype != "A") { ## X chromosome
            cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
            cross.scheme[2] <- 0
        }
    }
    else
        stop("markerloglik not available for cross type ", type, ".")

    ## Hide cross scheme in genoprob to pass to routine. BY
    temp <- 0
    if(type == "bcsft")
        temp[1] <- cross.scheme[1] * 1000 + cross.scheme[2]

    # call the C function
    z <- .C(cfunc,
            as.integer(n.ind),       # number of individuals
            as.integer(g),           # genotype data
            as.double(error.prob),
            loglik=as.double(temp),     # log likelihood
            PACKAGE="qtl")

    z$loglik
}


######################################################################
# allchrsplits
#
# get LOD scores for each possible split of each chromosome into
# two pieces
#
######################################################################

allchrsplits <-
    function(cross, chr, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE, verbose=TRUE)
{
    map.function <- match.arg(map.function)
    if(!missing(chr)) cross <- subset(cross, chr=chr)

    biggap <- imf.h(0.5 - 1e-14)

    n.mar <- nmar(cross)
    chrnam <- names(cross$geno)

    result <- NULL

    for(i in seq(along=cross$geno)) {
        if(n.mar[i] == 1) {
            #      temp <- data.frame(chr=chrnam[i], pos=pos, lod=NA, gap=0)
            #      rownames(temp) <- themarkers
            #      result <- cbind(result, temp)
            next
        }
        thischr <- subset(cross, chr=chrnam[i])
        if(verbose) cat("Chr ", chrnam[i], " (", n.mar[i], " markers)\n", sep="")
        pos <- cross$geno[[i]]$map
        themarkers <- colnames(cross$geno[[i]]$data)
        if(is.matrix(pos)) {
            pos <- pos[1,]
            matrixmap <- TRUE
        }
        else matrixmap <- FALSE


        gap <- (pos[-1] - pos[-length(pos)])
        pos <- (pos[-1] + pos[-length(pos)])/2
        int2 <- match(themarkers, colnames(cross$geno[[i]]$data))
        interval <- paste(themarkers[-length(themarkers)], themarkers[-1],
                          sep="-")
        int2 <- paste("(", int2[-length(int2)], "-", int2[-1], ")", sep="")
        interval <- paste(interval, int2)

        initialmap <- est.map(thischr, error.prob=error.prob,
                              map.function=map.function, m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp)
        thischr$geno[[1]]$map <- initialmap[[1]]
        initialloglik <- attr(initialmap[[1]], "loglik")

        if(n.mar[i] == 2) { # 2 markers
            mmll <- markerloglik(thischr, markernames(thischr), error.prob=error.prob)
            temp <- data.frame(chr=chrnam[i], pos=pos,
                               lod=(initialloglik - sum(mmll))/log(10), gap=gap, stringsAsFactors=TRUE)
            rownames(temp) <- interval
        }
        else { # >2 markers
            mn <- markernames(thischr)
            mmll <- markerloglik(thischr, mn[c(1,length(mn))], error.prob=error.prob)
            lod <- rep(NA, length(mn)-1)

            if(verbose) cat("  interval 1\n")
            # first interval
            lod[1] <- initialloglik - mmll[1] -
                attr(est.map(drop.markers(thischr, mn[1]), error.prob=error.prob,
                             map.function=map.function, m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp)[[1]], "loglik")

            if(n.mar[i] > 3) {
                for(j in 2:(n.mar[i]-2)) {
                    if(verbose) cat("  interval", j, "\n")
                    temp1 <- est.map(pull.markers(thischr, mn[1:j]), error.prob=error.prob,
                                     map.function=map.function, m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp)[[1]]
                    temp2 <- est.map(drop.markers(thischr, mn[1:j]), error.prob=error.prob,
                                     map.function=map.function, m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp)[[1]]
                    #          lod[j] <- initialloglik - attr(temp1, "loglik") - attr(temp2, "loglik")

                    if(any(is.na(temp1)) || any(is.na(temp2)))
                        stop("Missing values in estimated map on chr ", chrnam[i], " with split at interval ", j, "\n")

                    # the likelihoods aren't adding properly, so I'll use the following kluge:
                    temp3 <- thischr
                    if(is.matrix(temp1))
                        temp3$geno[[1]]$map <- cbind(temp1, biggap+temp2)
                    else
                        temp3$geno[[1]]$map <- c(temp1, biggap+temp2)

                    temp3 <- est.map(temp3, error.prob=error.prob, map.function=map.function, m=m, p=p, maxit=0, tol=tol,
                                     sex.sp=sex.sp)[[1]]
                    lod[j] <- initialloglik - attr(temp3, "loglik")
                }
            }

            if(verbose) cat("  interval", n.mar[i]-1, "\n")
            # last interval
            lod[length(mn)-1] <- initialloglik - mmll[2] -
                attr(est.map(drop.markers(thischr, mn[length(mn)]), error.prob=error.prob,
                             map.function=map.function, m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp)[[1]], "loglik")

            temp <- data.frame(chr=rep(chrnam[i], length(interval)), pos=pos, lod=lod/log(10), gap=gap, stringsAsFactors=TRUE)
            rownames(temp) <- interval
        }
        result <- rbind(result, temp)
    }
    class(result) <- c("scanone", "data.frame")
    result
}

# end of tryallpositions.R
