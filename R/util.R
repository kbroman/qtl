#####################################################################
#
# util.R
#
# copyright (c) 2001-2015, Karl W Broman
#     [find.pheno, find.flanking, and a modification to create.map
#      from Brian Yandell]
# last modified Aug, 2015
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
# Contains: markernames, c.cross, create.map, reduce2grid,
#           clean, clean.cross, drop.nullmarkers, nullmarkers
#           drop.markers, pull.markers, drop.dupmarkers
#           geno.table, genotab.em
#           mf.k, mf.h, imf.k, imf.h, mf.cf, imf.cf, mf.m, imf.m,
#           mf.stahl, imf.stahl
#           switch.order, flip.order, makeSSmap,
#           subset.cross, fill.geno, checkcovar, find.marker,
#           find.pseudomarker,
#           adjust.rf.ri, lodint, bayesint,
#           comparecrosses, movemarker, summary.map (aka summaryMap),
#           print.summary.map, find.pheno,
#           convert, convert.scanone, convert.scantwo
#           find.flanking, strip.partials, comparegeno
#           qtlversion, locateXO, jittermap, getid,
#           find.markerpos, geno.crosstab, LikePheVector,
#           matchchr, convert2sa, charround, testchr,
#           scantwoperm2scanoneperm, subset.map, [.map, [.cross,
#           findDupMarkers, convert2riself, convert2risib,
#           switchAlleles, nqrank, cleanGeno, typingGap,
#           calcPermPval, phenames, updateParallelRNG
#
######################################################################

######################################################################
#
# markernames
#
# pull out the marker names for selected chromosomes as one big vector
#
######################################################################

markernames <-
    function(cross, chr)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(!missing(chr)) cross <- subset(cross, chr=chr)
    temp <- unlist(lapply(cross$geno, function(a) colnames(a$data)))
    names(temp) <- NULL
    temp
}

######################################################################
#
# chrnames
#
# pull out the chrnames for a cross
#
######################################################################

chrnames <-
    function(cross)
{
    names(cross$geno)
}

######################################################################
#
# create.map
#
# create a new map with inserted inter-marker locations
#
# Note: map is a vector or a matrix with 2 rows
#
# stepwidth = "fixed" is the original R/qtl version;
# stepwidth="variable" is for Brian Yandell and the qtlbim package
# stepwidth="max" creates the minimal number of inserted pseudomarkers
#                 to have the maximum stepwidth = step
######################################################################
create.map <-
    function(map, step, off.end, stepwidth = c("fixed", "variable", "max"))
{
    stepwidth <- match.arg(stepwidth)
    if(step<0 || off.end<0) stop("step and off.end must be > 0.")

    if(!is.matrix(map)) { # sex-ave map
        if(stepwidth == "variable") {
            if(off.end > 0) {
                tmp <- names(map)
                ## Append and prepend by off.end value (exact here, no roundoff).
                map <- c(map[1] - off.end, map, map[length(map)] + off.end)
                names(map) <- c("loc000", tmp, "loc999")
            }
            if(step == 0)
                return(unclass(map))

            ## Determine differences and expansion vector.
            dif <- diff(map)
            expand <- pmax(1, floor(dif / step))

            ## Create pseudomarker map.
            a <- min(map) + cumsum(c(0, rep(dif / expand, expand)))

            ## Names are marker names or locNNN.
            namesa <- paste("loc", seq(length(a)), sep = "")
            namesa[cumsum(c(1, expand))] <- names(map)
            names(a) <- namesa

            return(unclass(a))
        }
        if(stepwidth == "max") {
            if(off.end > 0) {
                toadd <- c(map[1] - off.end, map[length(map)]+off.end)

                if(step==0) {
                    names(toadd) <- paste("loc", 1:2, sep="")
                    map <- sort(c(map, toadd))
                    return(unclass(map))
                }

                nmap <- c(map[1] - off.end, map, map[length(map)]+off.end)
            }
            else {
                nmap <- map
                toadd <- NULL
            }

            if(step==0 || (length(map)==1 && off.end==0)) return(unclass(map))

            d <- diff(nmap)
            nadd <- ceiling(d/step)-1
            if(sum(nadd) > 0) {
                for(j in 1:(length(nmap)-1)) {
                    if(nadd[j]>0)
                        toadd <- c(toadd, seq(nmap[j], nmap[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
                }
            }
            if(length(toadd) > 0)  {
                names(toadd) <- paste("loc", 1:length(toadd), sep="")
                map <- sort(c(map, toadd))
            }
            return(unclass(map))
        }

        if(length(map) == 1) { # just one marker!
            if(off.end==0) {
                if(step == 0) step <- 1
                nam <- names(map)
                map <- c(map,map+step)
                names(map) <- c(nam,paste("loc",step,sep=""))
            }
            else {
                if(step==0) m <- c(-off.end,off.end)
                else m <- seq(-off.end,off.end,by=step)
                m <- m[m!=0]
                names(m) <- paste("loc",m,sep="")
                map <- sort(c(m+map,map))
            }
            return(map)
        }

        minloc <- min(map)
        map <- map-minloc

        if(step==0 && off.end==0) return(map+minloc)
        else if(step==0 && off.end > 0) {
            a <- c(floor(min(map)-off.end),ceiling(max(map)+off.end))
            names(a) <- paste("loc", a, sep="")
            return(sort(c(a,map))+minloc)
        }
        else if(step>0 && off.end == 0) {
            a <- seq(floor(min(map)),max(map),
                     by = step)
            if(any(is.na(match(a, map)))) {
                a <- a[is.na(match(a,map))]
                names(a) <- paste("loc",a,sep="")
                return(sort(c(a,map))+minloc)
            }
            else return(map+minloc)
        }
        else {
            a <- seq(floor(min(map)-off.end),ceiling(max(map)+off.end+step),
                     by = step)
            a <- a[is.na(match(a,map))]

            # no more than one point above max(map)+off.end
            z <- (seq(along=a))[a >= max(map)+off.end]
            if(length(z) > 1) a <- a[-z[-1]]

            names(a) <- paste("loc",a,sep="")
            return(sort(c(a,map))+minloc)
        }
    } # end sex-ave map
    else { # sex-specific map
        if(stepwidth == "variable") {
            if(off.end > 0) {
                tmp <- colnames(map)
                map <- cbind(map[, 1] - off.end, map, map[, ncol(map)] + off.end)
                dimnames(map) <- list(NULL, c("loc000", tmp, "loc999"))
            }
            if(step == 0)
                return(unclass(map))

            ## Determine differences and expansion vector.
            dif <- diff(map[1, ])
            expand <- pmax(1, floor(dif / step))

            ## Create pseudomarker map.
            a <- min(map[1, ]) + cumsum(c(0, rep(dif / expand, expand)))
            b <- min(map[2, ]) + cumsum(c(0, rep(diff(map[2, ]) / expand, expand)))

            namesa <- paste("loc", seq(length(a)), sep = "")
            namesa[cumsum(c(1, expand))] <- dimnames(map)[[2]]
            map <- rbind(a,b)
            dimnames(map) <- list(NULL, namesa)

            return(unclass(map))
        }
        if(stepwidth == "max") {
            if(step==0 && off.end==0) return(unclass(map))
            if(step==0 && off.end>0) {
                if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
                    L1 <- L2 <- 1
                }
                else {
                    L1 <- diff(range(map[1,]))
                    L2 <- diff(range(map[2,]))
                }

                nam <- colnames(map)
                nmap1 <- c(map[1,1]-off.end, map[1,], map[1,ncol(map)]+off.end)
                nmap2 <- c(map[2,1]-off.end*L2/L1, map[2,], map[2,ncol(map)]+off.end*L2/L1)
                map <- rbind(nmap1, nmap2)
                colnames(map) <- c("loc1", nam, "loc2")
                return(unclass(map))
            }

            if(ncol(map)==1) L1 <- L2 <- 1
            else {
                L1 <- diff(range(map[1,]))
                L2 <- diff(range(map[2,]))
            }

            nam <- colnames(map)

            if(off.end > 0) {
                toadd1 <- c(map[1,1] - off.end, map[1,ncol(map)]+off.end)
                toadd2 <- c(map[2,1] + off.end*L2/L1, map[2,ncol(map)]+off.end*L2/L1)

                neword <- order(c(map[1,], toadd1))
                nmap1 <- c(map[1,], toadd1)[neword]
                nmap2 <- c(map[2,], toadd2)[neword]
            }
            else {
                nmap1 <- map[1,]
                nmap2 <- map[2,]
                toadd1 <- toadd2 <- NULL
            }

            d <- diff(nmap1)
            nadd <- ceiling(d/step)-1
            if(sum(nadd) > 0) {
                for(j in 1:(length(nmap1)-1)) {
                    if(nadd[j]>0) {
                        toadd1 <- c(toadd1, seq(nmap1[j], nmap1[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
                        toadd2 <- c(toadd2, seq(nmap2[j], nmap2[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
                    }
                }
            }
            newnam <- paste("loc", 1:length(toadd1), sep="")

            toadd1 <- sort(toadd1)
            toadd2 <- sort(toadd2)
            neword <- order(c(map[1,], toadd1))
            nmap1 <- c(map[1,], toadd1)[neword]
            nmap2 <- c(map[2,], toadd2)[neword]
            map <- rbind(nmap1, nmap2)
            colnames(map) <- c(nam, newnam)[neword]

            return(unclass(map))
        }

        minloc <- c(min(map[1,]),min(map[2,]))
        map <- unclass(map-minloc)
        markernames <- colnames(map)

        if(step==0 && off.end==0) return(map+minloc)
        else if(step==0 && off.end > 0) {
            map <- map+minloc
            if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
                L1 <- L2 <- 1
            }
            else {
                L1 <- diff(range(map[1,]))
                L2 <- diff(range(map[2,]))
            }

            nam <- colnames(map)
            nmap1 <- c(map[1,1]-off.end, map[1,], map[1,ncol(map)]+off.end)
            nmap2 <- c(map[2,1]-off.end*L2/L1, map[2,], map[2,ncol(map)]+off.end*L2/L1)
            map <- rbind(nmap1, nmap2)
            colnames(map) <- c("loc1", nam, "loc2")
            return(map)
        }
        else if(step>0 && off.end == 0) {

            if(ncol(map)==1) return(map+minloc)

            a <- seq(floor(min(map[1,])),max(map[1,]),
                     by = step)
            a <- a[is.na(match(a,map[1,]))]

            if(length(a)==0) return(map+minloc)

            b <- sapply(a,function(x,y,z) {
                ZZ <- min((seq(along=y))[y > x])
                (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1] }, map[1,],map[2,])

            m1 <- c(a,map[1,])
            m2 <- c(b,map[2,])

            names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
            return(rbind(sort(m1),sort(m2))+minloc)
        }
        else {
            a <- seq(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end+step),
                     by = step)
            a <- a[is.na(match(a,map[1,]))]
            # no more than one point above max(map)+off.end
            z <- (seq(along=a))[a >= max(map[1,])+off.end]
            if(length(z) > 1) a <- a[-z[-1]]

            b <- sapply(a,function(x,y,z,ml) {
                if(x < min(y)) {
                    return(min(z) - (min(y)-x)/diff(range(y))*diff(range(z)) - ml)
                }
                else if(x > max(y)) {
                    return(max(z) + (x - max(y))/diff(range(y))*diff(range(z)) - ml)
                }
                else {
                    ZZ <- min((seq(along=y))[y > x])
                    (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1]
                }
            }, map[1,],map[2,], minloc[2])
            m1 <- c(a,map[1,])
            m2 <- c(b,map[2,])
            names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
            return(rbind(sort(m1),sort(m2))+minloc)
        }
    }
}

######################################################################
# reduce2grid
#
# for high-density marker data, rather than run scanone at both the
# markers and at a set of pseudomarkers, we could reduce to just
# a set of evenly-spaced pseudomarkers
#
# first run calc.genoprob (or sim.geno) and then use this.
######################################################################
reduce2grid <-
    function(cross)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    # sample one element from a vector
    sampleone <- function(x) ifelse(length(x)==1, x, sample(x, 1))

    # for a map containing a grid with a given step size,
    # find the grid min(map), min(map)+step, min(map)+2step, ...
    gridindex <- function(map, step) {
        if(is.matrix(map)) stop("reduce2grid isn't working for sex-specific maps")
        grid <- seq(min(map), max(map), by=step)
        index <- match(grid, map)
        if(any(is.na(index)))
            index <- sapply(grid, function(a,b) { d <- abs(a-b); sampleone(which(d == min(d))) }, map)
        index
    }

    attr2fix <- c("error.prob", "step", "off.end", "map.function", "stepwidth")

    reduced <- FALSE
    if("prob" %in% names(cross$geno[[1]])) {
        stepwidth <- attr(cross$geno[[1]]$prob, "stepwidth")
        if(stepwidth != "fixed") {
            warning("You need to have run calc.genoprob with stepwidth=\"fixed\".")
            break
        }

        step <- attr(cross$geno[[1]]$prob, "step")

        for(i in 1:nchr(cross)) {
            pr <- cross$geno[[i]]$prob
            map <- attr(pr, "map")
            butes <- attributes(pr)

            reduced <- gridindex(map, step)
            pr <- pr[,reduced,,drop=FALSE]
            attr(pr, "map") <- map[reduced]
            for(a in attr2fix)
                attr(pr, a) <- butes[[a]]

            attr(pr, "reduced2grid") <- TRUE
            cross$geno[[i]]$prob <- pr
        }


        reduced <- TRUE
    }
    if("draws" %in% names(cross$geno[[1]])) {
        stepwidth <- attr(cross$geno[[1]]$draws, "stepwidth")
        if(stepwidth != "fixed") {
            warning("You need to have run sim.geno with stepwidth=\"fixed\".")
            break
        }

        step <- attr(cross$geno[[1]]$draws, "step")

        for(i in 1:nchr(cross)) {
            dr <- cross$geno[[i]]$draws
            map <- attr(dr, "map")
            butes <- attributes(dr)

            reduced <- gridindex(map, step)
            dr <- dr[,reduced,,drop=FALSE]
            attr(dr, "map") <- map[reduced]
            for(a in attr2fix)
                attr(dr, a) <- butes[[a]]

            attr(dr, "reduced2grid") <- TRUE
            cross$geno[[i]]$draws <- dr
        }

        reduced <- TRUE
    }

    if(!reduced)
        warning("You first need to run calc.genoprob or sim.geno with stepwidth=\"fixed\".")

    cross
}


######################################################################
# clean functions
######################################################################
clean <-
    function(object, ...)
    UseMethod("clean")

######################################################################
#
# clean.cross
#
# remove all of the extraneous stuff from a cross object, to get back
# to just the data
#
######################################################################

clean.cross <-
    function(object, ...)
{
    if(!any(class(object) == "cross"))
        stop("Input should have class \"cross\".")

    cross2 <- list(geno=object$geno,pheno=object$pheno)

    if("cross" %in% names(object))
        cross2$cross <- object$cross

    if("founderGeno" %in% names(object))
        cross2$founderGeno <- object$founderGeno

    if(!is.null(attr(object, "alleles")))
        attr(cross2, "alleles") <- attr(object, "alleles")

    if(!is.null(attr(object, "scheme")))
        attr(cross2, "scheme") <- attr(object, "scheme")

    for(i in 1:length(object$geno)) {
        cross2$geno[[i]] <- list(data=object$geno[[i]]$data,
                                 map=object$geno[[i]]$map)
        class(cross2$geno[[i]]) <- class(object$geno[[i]])
    }

    class(cross2) <- class(object)
    cross2
}


######################################################################
#
# drop.qtlgeno
#
# remove any QTLs from the genotype data and the genetic maps
# from data simulated via sim.cross. (They all have names "QTL*")
#
######################################################################

#drop.qtlgeno <-
#function(cross)
#{
#  n.chr <- nchr(cross)
#  mar.names <- lapply(cross$geno, function(a) {
#    m <- a$map
#    if(is.matrix(m)) return(colnames(m))
#    else return(names(m)) } )
#
#  for(i in 1:n.chr) {
#    o <- grep("^QTL[0-9]+",mar.names[[i]])
#    if(length(o) != 0) {
#      cross$geno[[i]]$data <- cross$geno[[i]]$data[,-o,drop=FALSE]
#      if(is.matrix(cross$geno[[i]]$map))
#        cross$geno[[i]]$map <- cross$geno[[i]]$map[,-o,drop=FALSE]
#      else
#        cross$geno[[i]]$map <- cross$geno[[i]]$map[-o]
#    }
#  }
#  cross
#}

######################################################################
#
# drop.nullmarkers
#
# remove markers that have no genotype data from the data matrix and
# genetic maps
#
######################################################################

drop.nullmarkers <-
    function(cross)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    n.chr <- nchr(cross)

    keep.chr <- rep(TRUE,n.chr)
    for(i in 1:n.chr) {
        o <- !apply(cross$geno[[i]]$data,2,function(a) sum(!is.na(a)))
        if(any(o)) { # remove from genotype data and map
            mn.drop <- colnames(cross$geno[[i]]$data)[o]
            if(length(mn.drop) == ncol(cross$geno[[i]]$data))
                keep.chr[i] <- FALSE # removing all markers from this chromosome

            cross$geno[[i]]$data <- cross$geno[[i]]$data[,!o,drop=FALSE]

            if(is.matrix(cross$geno[[i]]$map))
                cross$geno[[i]]$map <- cross$geno[[i]]$map[,!o,drop=FALSE]
            else
                cross$geno[[i]]$map <- cross$geno[[i]]$map[!o]

            # results of calc.genoprob
            if("prob" %in% names(cross$geno[[i]])) {
                o <- match(mn.drop,colnames(cross$geno[[i]]$prob))
                cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,-o,,drop=FALSE]
            }

            # results of argmax.geno
            if("argmax" %in% names(cross$geno[[i]])) {
                o <- match(mn.drop,colnames(cross$geno[[i]]$argmax))
                cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-o,drop=FALSE]
            }

            # results of sim.geno
            if("draws" %in% names(cross$geno[[i]])) {
                o <- match(mn.drop,colnames(cross$geno[[i]]$draws))
                cross$geno[[i]]$draws <- cross$geno[[i]]$draws[,-o,,drop=FALSE]
            }

            # results of est.rf
            if("rf" %in% names(cross)) {
                attrib <- attributes(cross$rf)

                o <- match(mn.drop,colnames(cross$rf))
                cross$rf <- cross$rf[-o,-o]

                if("onlylod" %in% names(attrib)) # save the onlylod attribute if its there
                    attr(cross$rf, "onlylod") <- attrib$onlylod
            }
        }
    }
    cross$geno <- cross$geno[keep.chr]

    if("founderGeno" %in% names(cross))
        cross$founderGeno <- cross$founderGeno[,markernames(cross)]

    cross
}


######################################################################
#
# nullmarkers
#
# identify markers that have no genotype data from the data matrix and
# genetic maps
#
######################################################################

nullmarkers <-
    function(cross)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    n.chr <- nchr(cross)

    keep.chr <- rep(TRUE,n.chr)
    all2drop <- NULL
    for(i in 1:n.chr) {
        o <- !apply(cross$geno[[i]]$data,2,function(a) sum(!is.na(a)))
        if(any(o)) { # remove from genotype data and map
            mn.drop <- colnames(cross$geno[[i]]$data)[o]
            all2drop <- c(all2drop, mn.drop)
        }
    }

    all2drop
}

######################################################################
#
# drop.markers
#
# remove a vector of markers from the data matrix and genetic maps
#
######################################################################

drop.markers <-
    function(cross, markers)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    n.chr <- nchr(cross)

    keep.chr <- rep(TRUE,n.chr)
    found <- rep(FALSE, length(markers))
    for(i in 1:n.chr) {
        # find markers on this chromosome
        o <- match(markers,colnames(cross$geno[[i]]$data))
        found[!is.na(o)] <- TRUE
        o <- o[!is.na(o)]
        a <- rep(FALSE,ncol(cross$geno[[i]]$data))
        a[o] <- TRUE
        o <- a

        if(any(o)) { # remove from genotype data and map
            mn.drop <- colnames(cross$geno[[i]]$data)[o]
            if(length(mn.drop) == ncol(cross$geno[[i]]$data))
                keep.chr[i] <- FALSE # removing all markers from this chromosome

            cross$geno[[i]]$data <- cross$geno[[i]]$data[,!o,drop=FALSE]

            if(is.matrix(cross$geno[[i]]$map))
                cross$geno[[i]]$map <- cross$geno[[i]]$map[,!o,drop=FALSE]
            else
                cross$geno[[i]]$map <- cross$geno[[i]]$map[!o]

            # results of calc.genoprob
            if("prob" %in% names(cross$geno[[i]])) {
                o <- match(mn.drop,colnames(cross$geno[[i]]$prob))
                cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,-o,,drop=FALSE]
            }

            # results of argmax.geno
            if("argmax" %in% names(cross$geno[[i]])) {
                o <- match(mn.drop,colnames(cross$geno[[i]]$argmax))
                cross$geno[[i]]$argmax <- cross$geno[[i]]$argmax[,-o,drop=FALSE]
            }

            # results of sim.geno
            if("draws" %in% names(cross$geno[[i]])) {
                o <- match(mn.drop,colnames(cross$geno[[i]]$draws))
                cross$geno[[i]]$draws <- cross$geno[[i]]$draws[,-o,,drop=FALSE]
            }

            # results of est.rf
            if("rf" %in% names(cross)) {
                attrib <- attributes(cross$rf)

                o <- match(mn.drop,colnames(cross$rf))
                cross$rf <- cross$rf[-o,-o]

                if("onlylod" %in% names(attrib))
                    attr(cross$rf, "onlylod") <- attrib$onlylod
            }
        }
    }

    if(sum(keep.chr) == 0)
        stop("You're attempting to drop *all* of the markers, which isn't allowed.")

    if(any(!found))
        warning("Markers not found: ", paste(markers[!found],collapse=" "))

    cross$geno <- cross$geno[keep.chr]

    if("founderGeno" %in% names(cross))
        cross$founderGeno <- cross$founderGeno[,markernames(cross)]

    cross
}

######################################################################
# pull.markers
#
# like drop.markers, but retain just those indicated
######################################################################
pull.markers <-
    function(cross, markers)
{
    mn <- markernames(cross)
    m <- match(markers, mn)
    if(any(is.na(m)))
        warning("Some markers couldn't be found: ", paste(markers[is.na(m)], collapse=" "))
    drop.markers(cross, mn[is.na(match(mn, markers))])
}

######################################################################
# drop.dupmarkers
#
# drop duplicate markers, retaining the consensus genotypes
######################################################################
drop.dupmarkers <-
    function(cross, verbose=TRUE)
{
    mn <- markernames(cross)
    tab <- table(mn)
    if(all(tab==1)) {
        if(verbose) cat("No duplicate markers.\n")
        return(cross)
    }

    dup <- names(tab[tab > 1])
    if(verbose) cat(" ", length(dup), "duplicate markers\n")

    # get consensus genotypes
    g <- pull.geno(cross)[,!is.na(match(mn, dup))]
    ng.omitted <- rep(NA, length(dup))
    tot.omitted <- 0
    nmar.omitted <- 0
    for(i in seq(along=dup)) {
        gg <- g[,colnames(g)==dup[i],drop=FALSE]
        res <- apply(gg, 1, function(a)
                 {
                     if(all(is.na(a))) return(c(NA, 0))
                     a <- unique(a[!is.na(a)])
                     if(length(a)==1) return(c(a, 0))
                     return(c(NA, 1))
                 } )
        if(verbose>1) {
            cat("  ", dup[i], ":\t", sum(res[2,]), " genotypes omitted\n", sep="")
            if(sum(res[2,]) > 1) cat("      ", paste(which(res[2,]>0), collapse=" "), "\n")
        }
        tot.omitted <- tot.omitted + sum(res[2,])

        flag <- FALSE
        for(j in seq(along=cross$geno)) {
            mn <- colnames(cross$geno[[j]]$data)
            wh <- mn==dup[i]
            if(!any(wh)) next
            wh <- which(wh)
            if(!flag) {
                flag <- TRUE
                cross$geno[[j]]$data[,wh[1]] <- res[1,]
                if(length(wh)==1) next
                else wh <- wh[-1]
            }
            if(length(wh) > 0) {
                nmar.omitted <- nmar.omitted + length(wh)
                cross$geno[[j]]$data <- cross$geno[[j]]$data[,-wh,drop=FALSE]
                if(is.matrix(cross$geno[[j]]$map))
                    cross$geno[[j]]$map <- cross$geno[[j]]$map[,-wh,drop=FALSE]
                else
                    cross$geno[[j]]$map <- cross$geno[[j]]$map[-wh]
            }
        }
    }
    if(verbose) {
        cat("  Total genotypes omitted:", tot.omitted, "\n")
        cat("  Total markers omitted:  ", nmar.omitted, "\n")
    }

    if("founderGeno" %in% names(cross))
        cross$founderGeno <- cross$founderGeno[,markernames(cross)]

    clean(cross)
}



######################################################################
#
# geno.table
#
# create table showing observed numbers of individuals with each
# genotype at each marker
#
######################################################################

geno.table <-
    function(cross, chr, scanone.output=FALSE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(!missing(chr))
        cross <- subset(cross, chr=chr)

    n.chr <- nchr(cross)

    type <- class(cross)[1]
    is.bcs <- type == "bcsft"
    cross.scheme <- attr(cross, "scheme")
    if(is.bcs)
        is.bcs <- (cross.scheme[2] == 0)

    chrtype <- sapply(cross$geno, class)
    allchrtype <- rep(chrtype, nmar(cross))
    chrname <- names(cross$geno)
    allchrname <- rep(chrname, nmar(cross))

    ## Actually plan to have our own geno.table.bcsft routine.
    if(type == "f2" || (type == "bcsft" && !is.bcs)) {
        n.gen <- 5
        temp <- getgenonames("f2", "A", cross.attr=attributes(cross))
        gen.names <- c(temp, paste("not", temp[c(3,1)]))
    }
    else if(type %in% c("bc", "riself", "risib", "dh", "haploid", "bcsft")) {
        n.gen <- 2
        gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
    }
    else if(type == "4way") {
        n.gen <- 14
        temp <- getgenonames("4way", "A", cross.attr=attributes(cross))
        gen.names <- c(temp,
                       paste(temp[c(1,3)], collapse="/"),
                       paste(temp[c(2,4)], collapse="/"),
                       paste(temp[c(1,2)], collapse="/"),
                       paste(temp[c(3,4)], collapse="/"),
                       paste(temp[c(1,4)], collapse="/"),
                       paste(temp[c(2,3)], collapse="/"),
                       paste("not", temp[1]),
                       paste("not", temp[2]),
                       paste("not", temp[3]),
                       paste("not", temp[4]))

        gen.names[5:8] <- substr(temp[c(1,2,1,3)], c(1,1,2,2), c(1,1,2,2))
    }
    else stop("Unknown cross type: ",type)

    res <- lapply(cross$geno, function(a,ngen) {
        a <- a$data; a[is.na(a)] <- 0
        apply(a,2,function(b,ngen) table(factor(b,levels=0:ngen)),ngen)
    },n.gen)

    results <- NULL
    for(i in 1:length(res))
        results <- rbind(results,t(res[[i]]))

    colnames(results) <- c("missing",gen.names)
    rownames(results) <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

    pval <- rep(NA,nrow(results))
    if(type %in% c("bc","risib","riself","dh","haploid") || (type=="bcsft" & is.bcs)) {
        sexpgm <- getsex(cross)
        if((type == "bc" || type=="bcsft") && any(chrtype == "X") && !is.null(sexpgm$sex) && any(sexpgm$sex==1)) {
            for(i in which(allchrtype=="A")) {
                x <- results[i,2:3]
                if(sum(x) > 0)
                    pval[i] <- chisq.test(x,p=c(0.5,0.5))$p.value
                else pval[i] <- 1
            }

            gn <- getgenonames("bc","X", "full", sexpgm, attributes(cross))
            wh <- which(is.na(match(gn, colnames(results))))
            temp <- matrix(0, nrow=nrow(results), ncol=length(wh))
            colnames(temp) <- gn[wh]
            results <- cbind(results, temp)

            for(i in which(chrtype=="X")) {
                dat <- reviseXdata("bc", "full", sexpgm, geno=cross$geno[[i]]$data,
                                   cross.attr=attributes(cross))
                dat[is.na(dat)] <- 0
                tab <- t(apply(dat, 2, function(x) table(factor(x, levels=0:length(gn)))))
                colnames(tab) <- c("missing", gn)
                results[allchrname==chrname[i],] <- 0
                results[allchrname==chrname[i],colnames(tab)] <- tab

                for(j in 1:ncol(dat)) {
                    stat <- apply(table(sexpgm$sex, cross$geno[[i]]$data[,j]),1,
                                  function(a) if(length(a) > 1 && sum(a)>0) return(chisq.test(a,p=c(0.5,0.5))$stat)
                                  else return(0))
                    pval[allchrname==chrname[i]][j] <- 1-pchisq(sum(stat),length(stat))
                }
            }
            results <- cbind(results, P.value=pval)
        }
        else {
            for(i in 1:length(pval)) {
                x <- results[i,2:3]
                if(sum(x) > 0)
                    pval[i] <- chisq.test(x,p=c(0.5,0.5))$p.value
                else
                    pval[i] <- 1
            }
            results <- cbind(results, P.value=pval)
        }
    }
    else  if(type == "f2" || (type == "bcsft" && !is.bcs)) {
        sexpgm <- getsex(cross)

        ## F2 with set initial genotype probabilities.
        init.geno <- c(0.25,0.5,0.25,0.75,0.75)
        ## BCsFt initial genotype probabilities need to be computed.
        if(type == "bcsft") {
            ret <- .C("genotab_em_bcsft",
                      as.integer(cross.scheme),
                      init.geno = as.double(init.geno))
            init.geno <- ret$init.geno
        }

        for(i in which(allchrtype=="A")) {
            dat <- results[i,2:6]
            if(sum(dat)==0) pval[i] <- 1
            else if(dat[4]==0 && dat[5]==0)
                pval[i] <- chisq.test(dat[1:3],   p=init.geno[1:3]   )$p.value
            else if(all(dat[2:4]==0))
                pval[i] <- chisq.test(dat[c(1,5)],p=init.geno[c(1,5)])$p.value
            else if(all(dat[c(1,2,5)]==0))
                pval[i] <- chisq.test(dat[3:4],   p=init.geno[3:4]   )$p.value
            else { # this is harder: some dominant and some not
                pval[i] <- genotab.em(dat, init.geno)
            }
        }

        for(i in which(chrtype=="X")) {
            gn <- getgenonames("f2","X","full", getsex(cross), attributes(cross))
            wh <- which(is.na(match(gn, colnames(results))))
            temp <- matrix(0, nrow=nrow(results), ncol=length(wh))
            colnames(temp) <- gn[wh]
            results <- cbind(results, temp)

            dat <- reviseXdata("f2", "full", sexpgm, geno=cross$geno[[i]]$data,
                               cross.attr=attributes(cross))
            dat[is.na(dat)] <- 0
            tab <- t(apply(dat, 2, function(x) table(factor(x, levels=0:length(gn)))))
            colnames(tab) <- c("missing", gn)
            results[allchrname==chrname[i],] <- 0
            results[allchrname==chrname[i],colnames(tab)] <- tab

            cn <- colnames(results)
            f <- grep("f$", cn)
            r <- grep("r$", cn)
            if(length(f)>0 && length(r)>0) {
                results <- results[,c(1:2,f,3,r,(1:ncol(results))[-c(1:3,f,r)])]
                colnames(results) <- cn[c(1:2,f,3,r,(1:ncol(results))[-c(1:3,f,r)])]
            }

            sex <- sexpgm$sex
            pgm <- sexpgm$pgm
            for(j in 1:ncol(dat)) {
                g <- cross$geno[[i]]$data[,j]
                if(!is.null(sex)) {
                    if(!is.null(pgm))
                        g <- matrix(as.numeric(table(g,sex,pgm)),ncol=2,byrow=TRUE)
                    else
                        g <- matrix(as.numeric(table(g,sex)),ncol=2,byrow=TRUE)
                }
                else {
                    if(!is.null(pgm)) {
                        g <- matrix(as.numeric(table(g,pgm)),ncol=2,byrow=TRUE)
                    }
                    else g <- matrix(table(g),ncol=2,byrow=TRUE)
                }
                stat <- apply(g, 1,
                              function(a) if(sum(a)>0) return(chisq.test(a,p=c(0.5,0.5))$stat)
                              else return(0))
                pval[allchrname==chrname[i]][j] <- 1-pchisq(sum(stat),length(stat))
            }
        }

        results <- cbind(results, P.value=pval)
    }
    else if(type == "4way") {
        for(i in 1:length(pval)) {
            x <- results[i,2:5]
            y <- results[i,-(1:5)]
            if(sum(x) > 0 && sum(y)==0)
                pval[i] <- chisq.test(x,p=c(0.25,0.25,0.25,0.25))$p.value
            else {
                if(allchrtype[i] == "A") {
                    res <- results[i,-1]
                    if(all(res==0)) pval[i] <- 1 # entirely missing
                    else if(all(res[-c(1,11)]==0))      # AC/not AC
                        pval[i] <- chisq.test(res[c(1,11)], p=c(0.25, 0.75))$p.value
                    else if(all(res[-c(2,12)]==0)) # BC/not BC
                        pval[i] <- chisq.test(res[c(2,12)], p=c(0.25, 0.75))$p.value
                    else if(all(res[-c(3,13)]==0)) # AD/not AD
                        pval[i] <- chisq.test(res[c(3,13)], p=c(0.25, 0.75))$p.value
                    else if(all(res[-c(4,14)]==0)) # BD/not BD
                        pval[i] <- chisq.test(res[c(4,14)], p=c(0.25, 0.75))$p.value
                    else if(all(res[-c(5,6)]==0)) # A/B
                        pval[i] <- chisq.test(res[c(5,6)], p=c(0.5, 0.5))$p.value
                    else if(all(res[-c(7,8)]==0)) # C/D
                        pval[i] <- chisq.test(res[c(7,8)], p=c(0.5, 0.5))$p.value
                    else if(all(res[-c(9,10)]==0)) # AC/BD or AD/BC
                        pval[i] <- chisq.test(res[c(9,10)], p=c(0.5, 0.5))$p.value
                    else if(all(res[-c(2,4,5)]==0)) # BC/BD/A
                        pval[i] <- chisq.test(res[c(2,4,5)], p=c(0.25, 0.25, 0.5))$p.value
                    else if(all(res[-c(1,3,6)]==0)) # AC/AD/B
                        pval[i] <- chisq.test(res[c(1,3,6)], p=c(0.25, 0.25, 0.5))$p.value
                    else if(all(res[-c(3,4,7)]==0)) # AD/BD/C
                        pval[i] <- chisq.test(res[c(3,4,7)], p=c(0.25, 0.25, 0.5))$p.value
                    else if(all(res[-c(1,2,8)]==0)) # AC/BC/D
                        pval[i] <- chisq.test(res[c(1,2,8)], p=c(0.25, 0.25, 0.5))$p.value
                    else if(all(res[-c(2,3,9)]==0)) # AC/BD or AD or BC
                        pval[i] <- chisq.test(res[c(2,3,9)], p=c(0.25, 0.25, 0.5))$p.value
                    else if(all(res[-c(1,4,10)]==0)) # AD/BC or AC or BD
                        pval[i] <- chisq.test(res[c(1,4,10)], p=c(0.25, 0.25, 0.5))$p.value
                }
            }
        }
        results <- cbind(results, P.value=pval)
    }

    if(!scanone.output)
        return(data.frame(chr=rep(names(cross$geno),nmar(cross)),results, stringsAsFactors=TRUE))

    themap <- pull.map(cross)
    if(is.matrix(themap[[1]]))
        thepos <- unlist(lapply(themap, function(a) a[1,]))
    else thepos <- unlist(themap)

    temp <- results[,1:(ncol(results)-1),drop=FALSE]
    res <- data.frame(chr=rep(names(cross$geno),nmar(cross)),
                      pos=thepos,
                      neglog10P=-log10(results[,ncol(results)]),
                      missing=temp[,1]/apply(temp, 1, sum),
                      temp[,-1]/apply(temp[,-1], 1, sum), stringsAsFactors=TRUE)
    class(res) <- c("scanone", "data.frame")
    rownames(res) <- rownames(results)
    res[,1] <- factor(as.character(res[,1]), levels=unique(as.character(res[,1])))
    res
}



genotab.em <-
    function(dat, init.geno, tol=1e-6, maxit=10000, verbose=FALSE)
{

    genotab.ll <-
        function(dat, gam, init.geno)
        {
            p <- c(init.geno[1]*(1-gam[2]), init.geno[2]*gam[1], init.geno[3]*(1-gam[3]),
                   gam[2]*init.geno[4], gam[3]*init.geno[5])
            if(any(p==0 & dat > 0)) return(-Inf)
            return( sum((dat*log(p))[dat>0 & p>0]) )
        }

    n <- sum(dat)

    gam <- c(sum(dat[1:3]), dat[4], dat[5])/n

    curll <- genotab.ll(dat, gam, init.geno)

    flag <- 0
    if(verbose) cat(0, gam, curll, "\n")
    for(i in 1:maxit) {
        # estep
        zAA <- dat[1]*gam[3]/(gam[1]+gam[3])
        zBB <- dat[3]*gam[2]/(gam[1]+gam[2])

        # mstep
        gamnew <- c(sum(dat[1:3])-zAA-zBB, dat[4]+zBB, dat[5]+zAA)/n

        newll <- genotab.ll(dat, gamnew, init.geno)

        if(verbose) cat(i, gamnew, newll, "\n")

        if(abs(curll-newll) < tol) {
            flag <- 1
            break
        }
        gam <- gamnew
        curll <- newll
    }
    if(!flag) warning("Didn't converge.")

    gam <- gamnew
    p <- c(init.geno[1]*(1-gam[2]), init.geno[2]*gam[1], init.geno[3]*(1-gam[3]),
           gam[2]*init.geno[4], gam[3]*init.geno[5])
    ex <- p*n
    1-pchisq(sum(((dat-ex)^2/ex)[ex>0]), 2)
}

######################################################################
# geno.crosstab
#
# Get a cross-tabulation of the genotypes at two markers,
# with the markers specified by name
######################################################################
geno.crosstab <-
    function(cross, mname1, mname2, eliminate.zeros=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(missing(mname2) && length(mname1)>1) {
        mname2 <- mname1[2]
        mname1 <- mname1[1]
    }

    if(length(mname1) > 1 || length(mname2) > 1)
        stop("mname1 and mname2 should both have lenght 1, or mname1 should have length 2 and mname1 should be missing.")

    if(mname1==mname2)
        stop("You must give two distinct marker names.")

    pos <- find.markerpos(cross, c(mname1, mname2))
    if(any(is.na(pos$chr))) {
        if(all(is.na(pos$chr)))
            stop("Markers ", mname1, " and ", mname2, " not found.")
        else
            stop("Marker ", rownames(pos)[is.na(pos$chr)], " not found.")
    }

    chrtype <- sapply(cross$geno[pos$chr], class)
    crosstype <- class(cross)[1]

    g1 <- pull.geno(cross, pos$chr[1])[,mname1, drop=FALSE]
    g2 <- pull.geno(cross, pos$chr[2])[,mname2, drop=FALSE]

    if(chrtype[1] == "X")
        g1 <- reviseXdata(crosstype, "full", getsex(cross), geno=g1, cross.attr=attributes(cross))

    if(chrtype[2] == "X")
        g2 <- reviseXdata(crosstype, "full", getsex(cross), geno=g2, cross.attr=attributes(cross))

    g1[is.na(g1)] <- 0
    g2[is.na(g2)] <- 0

    g1names <- getgenonames(crosstype, chrtype[1], "full", getsex(cross), attributes(cross))
    g2names <- getgenonames(crosstype, chrtype[2], "full", getsex(cross), attributes(cross))

    if(chrtype[1] != "X") {
        if(crosstype == "f2")
            g1names <- c(g1names, paste("not", g1names[c(3,1)]))
        else if(crosstype == "bc" || crosstype == "risib" || crosstype=="riself" || crosstype=="dh" || crosstype=="haploid") {
        }
        else if(crosstype == "4way") {
            temp <- g1names
            g1names <- c(temp,
                         paste(temp[c(1,3)], collapse="/"),
                         paste(temp[c(2,4)], collapse="/"),
                         paste(temp[c(1,2)], collapse="/"),
                         paste(temp[c(3,4)], collapse="/"),
                         paste(temp[c(1,4)], collapse="/"),
                         paste(temp[c(2,3)], collapse="/"),
                         paste("not", temp[1]),
                         paste("not", temp[2]),
                         paste("not", temp[3]),
                         paste("not", temp[4]))

            g1names[5:8] <- substr(temp[c(1,2,1,3)], c(1,1,2,2), c(1,1,2,2))
        }
        else stop("Unknown cross type: ",crosstype)
    }
    if(chrtype[2] != "X") {
        if(crosstype == "f2")
            g2names <- c(g2names, paste("not", g2names[c(3,1)]))
        else if(crosstype == "bc" || crosstype == "risib" || crosstype=="riself" || crosstype=="dh" || crosstype=="haploid") {
        }
        else if(crosstype == "4way") {
            temp <- g2names
            g2names <- c(temp,
                         paste(temp[c(1,3)], collapse="/"),
                         paste(temp[c(2,4)], collapse="/"),
                         paste(temp[c(1,2)], collapse="/"),
                         paste(temp[c(3,4)], collapse="/"),
                         paste(temp[c(1,4)], collapse="/"),
                         paste(temp[c(2,3)], collapse="/"),
                         paste("not", temp[1]),
                         paste("not", temp[2]),
                         paste("not", temp[3]),
                         paste("not", temp[4]))

            g2names[5:8] <- substr(temp[c(1,2,1,3)], c(1,1,2,2), c(1,1,2,2))
        }
        else stop("Unknown cross type: ",crosstype)
    }
    g1names <- c("-", g1names)
    g2names <- c("-", g2names)

    g1 <- as.character(g1)
    g2 <- as.character(g2)

    for(i in 1:length(g1names)) {
        j <- as.character(i-1)
        g1[g1==j] <- g1names[i]
    }
    g1 <- factor(g1, levels=g1names)

    for(i in 1:length(g2names)) {
        j <- as.character(i-1)
        g2[g2==j] <- g2names[i]
    }
    g2 <- factor(g2, levels=g2names)

    tab <- table(g1, g2)
    names(attributes(tab)$dimnames) <- c(mname1, mname2)

    if(eliminate.zeros) { # eliminate rows and columns with no data (but always include missing data row and column)
        rs <- apply(tab, 1, sum); rs[1] <- 1
        tab <- tab[rs>0,,drop=FALSE]
        cs <- apply(tab, 2, sum); cs[1] <- 1
        tab <- tab[,cs>0,drop=FALSE]
    }

    tab
}


# map functions
mf.k <- function(d) 0.5*tanh(d/50)
mf.h <- function(d) 0.5*(1-exp(-d/50))
imf.k <- function(r) { r[r >= 0.5] <- 0.5-1e-14; 50*atanh(2*r) }
imf.h <- function(r) { r[r >= 0.5] <- 0.5-1e-14; -50*log(1-2*r) }
mf.m <- function(d) sapply(d,function(a) min(a/100,0.5))
imf.m <- function(r) sapply(r,function(a) min(a*100,50))

# carter-falconer: mf.cf, imf.cf
imf.cf <- function(r) { r[r >= 0.5] <- 0.5-1e-14; 12.5*(log(1+2*r)-log(1-2*r))+25*atan(2*r) }

mf.cf <-
    function(d)
{
    d[d >= 300] <- 300

    icf <- function(r,d)
        imf.cf(r)-d

    sapply(d,function(a) {
        if(a==0) return(0)
        uniroot(icf, c(0,0.5-1e-14),d=a,tol=1e-12)$root })
}


######################################################################
#
# switch.order: change the marker order on a given chromosome to some
#               specified order
#
######################################################################

switch.order <-
    function(cross, chr, order, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             maxit=4000, tol=1e-6, sex.sp=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    map.function <- match.arg(map.function)

    # check chromosome argument
    if(missing(chr)) {
        chr <- names(cross$geno)[1]
        warning("Assuming you mean chromosome ", chr)
    }
    else {
        if(length(chr) > 1) {
            chr <- chr[1]
            warning("switch.order can deal with just one chromosome at a time; assuming you want chr ", chr)
        }
        if(!testchr(chr, names(cross$geno)))
            stop("Chr ", chr, " not found.")
    }
    chr <- matchchr(chr, names(cross$geno))

    # check order argument
    n.mar <- nmar(cross)
    if(n.mar[chr] == length(order)-2 || n.mar[chr]==length(order)-1)
        order <- order[1:n.mar[chr]]     # useful for output from ripple()
    if(n.mar[chr] != length(order))
        stop("Incorrect number of markers.")

    # save recombination fractions
    flag <- 0
    if("rf" %in% names(cross)) {
        attrib <- attributes(cross$rf)

        rf <- cross$rf
        # determine column within rec fracs
        whchr <- which(names(cross$geno)==chr)
        oldcols <- cumsum(c(0,n.mar))[whchr]+seq(along=order)
        newcols <- cumsum(c(0,n.mar))[whchr]+order
        rf[oldcols,] <- rf[newcols,]
        rf[,oldcols] <- rf[,newcols]
        colnames(rf)[oldcols] <- colnames(rf)[newcols]
        rownames(rf)[oldcols] <- rownames(rf)[newcols]
        flag <- 1
    }

    # remove any intermediate calculations (except rec fracs),
    #   as they will no longer be meaningful
    cross <- clean(cross)

    if(!is.matrix(cross$geno[[chr]]$map))
        first <- min(cross$geno[[chr]]$map)
    else
        first <- apply(cross$geno[[chr]]$map,1,min)

    # re-order markers
    cross$geno[[chr]]$data <- cross$geno[[chr]]$data[,order,drop=FALSE]
    m <- seq(0,by=5,length=ncol(cross$geno[[chr]]$data))
    names(m) <- colnames(cross$geno[[chr]]$data)
    if(is.matrix(cross$geno[[chr]]$map))
        cross$geno[[chr]]$map <- rbind(m,m)
    else
        cross$geno[[chr]]$map <- m

    # re-estimate rec fracs for re-ordered chromosome
    if(flag==1) {
        temp <- est.rf(subset(cross, chr=chr))$rf
        rf[oldcols,oldcols] <- temp
        cross$rf <- rf
        if("onlylod" %in% names(attrib))
            attr(cross$rf, "onlylod") <- attrib$onlylod
    }

    # re-estimate map
    newmap <- est.map(cross, chr=chr,
                      error.prob=error.prob, map.function=map.function,
                      maxit=maxit, tol=tol, sex.sp=sex.sp)

    if(!is.matrix(newmap[[1]]))
        cross$geno[[chr]]$map <- newmap[[1]] + first
    else {
        cross$geno[[chr]]$map[1,] <- newmap[[1]][1,] + first[1]
        cross$geno[[chr]]$map[2,] <- newmap[[1]][2,] + first[2]
        rownames(cross$geno[[chr]]$map) <- NULL
    }

    cross
}


######################################################################
# flip.order: flip the order of markers on a set of chromosomes
######################################################################

flip.order <-
    function(cross, chr)
{
    # utility function to flip a map
    flipChrOnMap <-
        function(map)
        {
            if(is.matrix(map)) {
                map <- map[,ncol(map):1,drop=FALSE]
                for(j in 1:nrow(map))
                    map[j,] <- max(map[j,]) - map[j,]
            } else {
                map <- max(map) - map[length(map):1]
            }
            map
        }

    # utility function to flip intermediate calc
    flipChrInterCalc <-
        function(object, attr2consider=c("error.prob", "step", "off.end", "map.function", "stepwidth"))
        {
            the_attr <- attributes(object)
            ndim <- length(dim(object))
            if(ndim == 3)
                object <- object[,ncol(object):1,,drop=FALSE]
            else
                object <- object[,ncol(object):1,drop=FALSE]
            for(a in attr2consider) {
                if(a %in% names(the_attr))
                    attr(object, a) <- the_attr[[a]]
            }
            if("map" %in% names(the_attr))
                attr(object, "map") <- flipChrOnMap(the_attr$map)

            object
        }

    chr <- matchchr(chr, names(cross$geno))

    for(i in chr) {
        nc <- ncol(cross$geno[[i]]$data)
        cross$geno[[i]]$data <- cross$geno[[i]]$data[,nc:1,drop=FALSE]
        cross$geno[[i]]$map <- flipChrOnMap(cross$geno[[i]]$map)

        if("prob" %in% names(cross$geno[[i]]))
            cross$geno[[i]]$prob <- flipChrInterCalc(cross$geno[[i]]$prob)
        if("argmax" %in% names(cross$geno[[i]]))
            cross$geno[[i]]$argmax <- flipChrInterCalc(cross$geno[[i]]$argmax)
        if("draws" %in% names(cross$geno[[i]]))
            cross$geno[[i]]$draws <- flipChrInterCalc(cross$geno[[i]]$draws)
        if("errorlod" %in% names(cross$geno[[i]]))
            cross$geno[[i]]$errorlod <- flipChrInterCalc(cross$geno[[i]]$errorlod,
                                                          c("error.prob", "map.function"))
    }

    cross$rf <- NULL

    cross
}


######################################################################
#
# subset.cross: General subsetting function for a cross object
#
######################################################################

subset.cross <-
    function(x, chr, ind, ...)
{
    if(!any(class(x) == "cross"))
        stop("Input should have class \"cross\".")

    if(missing(chr) && missing(ind))
        stop("You must specify either chr or ind.")

    n.chr <- nchr(x)
    n.ind <- nind(x)

    # pull out relevant chromosomes
    if(!missing(chr)) {
        chr <- matchchr(chr, names(x$geno))

        if("rf" %in% names(x)) { # pull out part of rec fracs
            if(totmar(x) != ncol(x$rf))
                x <- clean(x)
            else {
                attrib <- attributes(x$rf)

                n.mar <- nmar(x)
                n.chr <- n.chr
                wh <- rbind(c(0,cumsum(n.mar)[-n.chr])+1,cumsum(n.mar))
                dimnames(wh) <- list(NULL, names(n.mar))
                wh <- wh[,chr,drop=FALSE]
                wh <- unlist(apply(wh,2,function(a) a[1]:a[2]))
                x$rf <- x$rf[wh,wh]

                if("onlylod" %in% names(attrib)) # save the onlylod attribute if its there
                    attr(x$rf, "onlylod") <- attrib$onlylod
            }
        }

        x$geno <- x$geno[chr]

        if("founderGeno" %in% names(x))
            x$founderGeno <- x$founderGeno[,unlist(lapply(x$geno, function(a) colnames(a$data)))]
    }

    if(!missing(ind)) {
        theid <- getid(x)

        if(is.logical(ind)) {
            ind[is.na(ind)] <- FALSE
            if(length(ind) != n.ind)
                stop("ind argument has wrong length (", length(ind), "; should be ", n.ind, ")")
            ind <- (1:n.ind)[ind]
        }

        else if(is.numeric(ind)) { # treat as numeric indices; don't match against individual identifiers
            if(all(ind < 0)) { # drop all but these
                ind <- -ind
                if(any(ind > nind(x))) {
                    a <- -ind[ind > nind(x)]
                    if(length(a) > 1) a <- sample(a, 1)
                    stop("Invalid ind values (e.g., ", a, ")")
                }
                ind <- (1:n.ind)[-ind]
            }
            else if(all(ind > 0)) { # keep these
                if(any(ind > nind(x))) {
                    a <- ind[ind > nind(x)]
                    if(length(a) > 1) a <- sample(a, 1)
                    stop("Invalid ind values (e.g., ", a, ")")
                }
            }
            else {
                stop("Need ind to be all > 0 or all < 0.")
            }
        }


        else if(!is.null(theid)) { # cross has individual IDs
            ind <- as.character(ind)
            if(all(substr(ind, 1,1) == "-")) {
                ind <- substr(ind, 2, nchar(ind))
                m <- match(ind, theid)
                if(all(is.na(m)))
                    stop("No matching individuals.")
                if(any(is.na(m)))
                    warning("Individuals not found: ", paste(ind[is.na(m)]))
                ind <- (1:n.ind)[-m[!is.na(m)]]
            }
            else  {
                m <- match(ind, theid)
                if(any(is.na(m)))
                    warning("Individuals not found: ", paste(ind[is.na(m)], collapse=" "))
                ind <- m[!is.na(m)]
            }
            ind <- ind[!is.na(ind)]
        }
        else { # no individual IDs
            stop("In the absense of individual IDs, ind should be logical or numeric.")
        }
        # Note: ind should now be a numeric vector

        if(length(ind) == 0)
            stop("Must retain at least one individual.")
        if(length(ind) == 1)
            warning("Retained only one individual!")

        x$pheno <- x$pheno[ind,,drop=FALSE]

        if("cross" %in% names(x))
            x$cross <- x$cross[ind,,drop=FALSE]

        for(i in 1:nchr(x)) {
            x$geno[[i]]$data <- x$geno[[i]]$data[ind,,drop=FALSE]

            if("prob" %in% names(x$geno[[i]])) {
                temp <- attributes(x$geno[[i]]$prob) # all attributes but dim and dimnames
                x$geno[[i]]$prob <- x$geno[[i]]$prob[ind,,,drop=FALSE]
                # put attributes back in
                for(k in seq(along=temp)) {
                    if(names(temp)[k] != "dim" && names(temp)[k] != "dimnames")
                        attr(x$geno[[i]]$prob, names(temp)[k]) <- temp[[k]]
                }
            }

            if("errorlod" %in% names(x$geno[[i]])) {
                temp <- attributes(x$geno[[i]]$prob) # all attributes but dim and dimnames
                x$geno[[i]]$errorlod <- x$geno[[i]]$errorlod[ind,,drop=FALSE]
                # put attributes back in
                for(k in seq(along=temp)) {
                    if(names(temp)[k] != "dim" && names(temp)[k] != "dimnames")
                        attr(x$geno[[i]]$errorlod, names(temp)[k]) <- temp[[k]]
                }
            }

            if("argmax" %in% names(x$geno[[i]])) {
                temp <- attributes(x$geno[[i]]$argmax) # all attributes but dim and dimnames
                x$geno[[i]]$argmax <- x$geno[[i]]$argmax[ind,,drop=FALSE]
                # put attributes back in
                for(k in seq(along=temp)) {
                    if(names(temp)[k] != "dim" && names(temp)[k] != "dimnames")
                        attr(x$geno[[i]]$argmax, names(temp)[k]) <- temp[[k]]
                }
            }

            if("draws" %in% names(x$geno[[i]])) {
                temp <- attributes(x$geno[[i]]$draws) # all attributes but dim and dimnames
                x$geno[[i]]$draws <- x$geno[[i]]$draws[ind,,,drop=FALSE]
                # put attributes back in
                for(k in seq(along=temp)) {
                    if(names(temp)[k] != "dim" && names(temp)[k] != "dimnames")
                        attr(x$geno[[i]]$draws, names(temp)[k]) <- temp[[k]]
                }
            }
        }

        if("qtlgeno" %in% names(x))
            x$qtlgeno <- x$qtlgeno[ind,,drop=FALSE]
    }
    x
}

#pull.chr <-
#function(cross, chr) {
#  warning("pull.chr is deprecated; use subset.cross.")
#  subset.cross(cross, chr)
#}


######################################################################
#
# c.cross: Combine crosses
#
######################################################################

c.cross <-
    function(...)
{
    args <- list(...)
    n.args <- length(args)

    for(i in seq(along=args)) {
        if(!any(class(args[[i]]) == "cross"))
            stop("Input should have class \"cross\".")
    }

    # if only one cross, just return it
    if(n.args==1) return(args[[1]])

    if(any(sapply(args, function(a) !any(class(a) == "cross"))))
        stop("All arguments must be cross objects.")

    # crosses must be all the same, or must be combination of F2 and BC
    classes <- sapply(args,function(a) class(a)[1])
    if(length(unique(classes))==1) {
        allsame <- TRUE
        type <- classes[1]
    }
    else {
        if(any(classes != "bc" & classes != "f2"))
            stop("Experiments must be either the same type or be bc/f2.")
        allsame <- FALSE
        type <- "f2"
    }

    if(length(unique(sapply(args, nchr))) > 1)
        stop("All arguments must have the same number of chromosomes.")

    x <- args[[1]]
    chr <- names(x$geno)
    n.mar <- nmar(x)
    marnam <- unlist(lapply(x$geno,function(b) colnames(b$data)))
    marpos <- unlist(lapply(x$geno,function(b) b$map))

    map.mismatch <- 0
    for(i in 2:n.args) {
        y <- args[[i]]
        y.marnam <- unlist(lapply(y$geno, function(b) colnames(b$data)))
        y.marpos <- unlist(lapply(y$geno, function(b) b$map))
        if(chr != names(y$geno) || any(n.mar != nmar(y)) ||
           any(marnam != y.marnam) || any(marpos != y.marpos)) {
            map.mismatch <- 1
            break
        }
    }
    if(map.mismatch) { # get the maps to line up
        for(i in 1:nchr(args[[1]])) {
            themap <- NULL

            themaps <- vector("list", n.args)

            for(j in 1:n.args)  {
                if(is.matrix(args[[j]]$map))
                    stop("c.cross() won't work with sex-specific maps.")

                themaps[[j]] <- args[[j]]$geno[[i]]$map
                themap <- c(themap, themaps[[j]])
            }
            themap <- sort(themap)
            mn <- unique(names(themap))

            newmap <- rep(0,length(mn))
            names(newmap) <- mn
            for(j in 1:length(newmap))
                newmap[j] <- mean(themap[names(themap) == mn[j]])

            for(j in 1:n.args) {
                m <- match(names(themaps[[j]]), mn)
                m <- m[!is.na(m)]
                if(any(diff(m)) < 0)
                    stop(" Markers must all be in the same order.")

                if(!all(mn %in% names(themaps[[j]]))) {
                    temp <- matrix(ncol=length(mn), nrow=nind(args[[j]]))
                    colnames(temp) <- mn
                    temp[,names(themaps[[j]])] <- args[[j]]$geno[[i]]$data
                    args[[j]]$geno[[i]]$data <- temp
                }
                args[[j]]$geno[[i]]$map <- newmap
            }

        }

    } # end of map mismatch fix

    x <- args[[1]]
    chr <- names(x$geno)
    n.mar <- nmar(x)
    marnam <- unlist(lapply(x$geno,function(b) colnames(b$data)))


    # create genotype information
    geno <- x$geno
    for(j in 1:nchr(x)) { # drop extraneous stuff
        geno[[j]] <- list(data=geno[[j]]$data, map=geno[[j]]$map)
        class(geno[[j]]) <- class(x$geno[[j]])
    }
    for(i in 2:n.args)
        for(j in 1:nchr(x))
            geno[[j]]$data <- rbind(geno[[j]]$data,args[[i]]$geno[[j]]$data)

    # get all phenotype names
    phenam <- names(x$pheno)
    for(i in 2:n.args)
        phenam <- c(phenam, names(args[[i]]$pheno))
    phenam <- unique(phenam)

    # form big phenotype matrix
    n.ind <- sapply(args,nind)
    pheno <- matrix(nrow=sum(n.ind),ncol=length(phenam))
    colnames(pheno) <- phenam
    pheno <- as.data.frame(pheno, stringsAsFactors=TRUE)

    if(!allsame) {
        crosstype <- factor(rep(c("bc","f2")[match(classes,c("bc","f2"))],n.ind),
                            levels=c("bc","f2"))
        pheno <- cbind(pheno,crosstype=crosstype)
    }

    for(i in 1:length(phenam)) {
        phe <- vector("list",n.args)
        for(j in 1:n.args) {
            o <- match(phenam[i],names(args[[j]]$pheno))
            if(is.na(o)) phe[[j]] <- rep(NA,n.ind[j])
            else phe[[j]] <- args[[j]]$pheno[,o]
        }
        pheno[,i] <- unlist(phe)
    }

    # indicator of which cross
    whichcross <- matrix(0,ncol=n.args,nrow=sum(n.ind))
    colnames(whichcross) <- paste("cross",1:n.args,sep="")
    thecross <- rep(NA, sum(n.ind))
    prev <- 0
    for(i in 1:n.args) {
        wh <- prev + 1:n.ind[i]
        prev <- prev + n.ind[i]
        whichcross[wh,i] <- 1
        thecross[wh] <- i
    }
    pheno <- cbind(pheno,cross=thecross,whichcross)

    x$pheno <- pheno

    if(!map.mismatch) { # keep probs and draws only if we've not re-aligned the maps
        # if probs exist in each and all have the same
        #     set up values, keep them
        wh <- sapply(args, function(a) match("prob",names(a$geno[[1]])))
        step <- sapply(args,function(a) attr(a$geno[[1]]$prob,"step"))
        error.prob <- sapply(args,function(a) attr(a$geno[[1]]$prob,"error.prob"))
        off.end <- sapply(args,function(a) attr(a$geno[[1]]$prob,"off.end"))
        map.function <- sapply(args,function(a) attr(a$geno[[1]]$prob,"map.function"))
        map <- lapply(args,function(a) attr(a$geno[[1]]$prob,"map"))
        if(!any(is.na(wh)) && length(unique(step))==1 &&
           length(unique(error.prob))==1 && length(unique(off.end))==1 &&
           length(unique(map.function))==1) {
            if(allsame) { # all same cross type
                for(j in 1:nchr(x)) {
                    geno[[j]]$prob <- array(dim=c(sum(n.ind),dim(x$geno[[j]]$prob)[-1]))
                    dimnames(geno[[j]]$prob) <- dimnames(x$geno[[j]]$prob)
                    prev <- 0
                    for(i in 1:n.args) {
                        wh <- prev + 1:n.ind[i]
                        prev <- prev + n.ind[i]
                        geno[[j]]$prob[wh,,] <- args[[i]]$geno[[j]]$prob
                    }
                }
            }
            else { # mixed F2 and BC
                for(j in 1:nchr(x)) {
                    wh <- match("f2",classes)
                    geno[[j]]$prob <- array(0,dim=c(sum(n.ind),dim(args[[wh]]$geno[[j]]$prob)[-1]))
                    dimnames(geno[[j]]$prob) <- dimnames(args[[wh]]$geno[[j]]$prob)
                    prev <- 0
                    for(i in 1:n.args) {
                        wh <- prev + 1:n.ind[i]
                        prev <- prev + n.ind[i]

                        if(classes[i]=="f2")
                            geno[[j]]$prob[wh,,] <- args[[i]]$geno[[j]]$prob
                        else # backcross
                            geno[[j]]$prob[wh,,1:2] <- args[[i]]$geno[[j]]$prob
                    }
                }
            }

            for(j in 1:nchr(x)) {
                wh <- sapply(args, function(a) match("prob",names(a$geno[[j]])))
                step <- sapply(args,function(a) attr(a$geno[[j]]$prob,"step"))
                error.prob <- sapply(args,function(a) attr(a$geno[[j]]$prob,"error.prob"))
                off.end <- sapply(args,function(a) attr(a$geno[[j]]$prob,"off.end"))
                map.function <- sapply(args,function(a) attr(a$geno[[j]]$prob,"map.function"))
                map <- lapply(args,function(a) attr(a$geno[[j]]$prob,"map"))

                attr(geno[[j]]$prob,"step") <- step[1]
                attr(geno[[j]]$prob,"error.prob") <- error.prob[1]
                attr(geno[[j]]$prob,"off.end") <- off.end[1]
                attr(geno[[j]]$prob,"map.function") <- map.function[1]
                attr(geno[[j]]$prob,"map") <- map[[1]]
            }
        }

        # if draws exist in each and all have the same
        #     set up values, keep them
        wh <- sapply(args, function(a) match("draws",names(a$geno[[1]])))
        step <- sapply(args,function(a) attr(a$geno[[1]]$draws,"step"))
        error.prob <- sapply(args,function(a) attr(a$geno[[1]]$draws,"error.prob"))
        off.end <- sapply(args,function(a) attr(a$geno[[1]]$draws,"off.end"))
        map.function <- sapply(args,function(a) attr(a$geno[[1]]$draws,"map.function"))
        map <- lapply(args,function(a) attr(a$geno[[1]]$draws,"map"))
        ndraws <- sapply(args,function(a) dim(a$geno[[1]]$draws)[3])
        if(!any(is.na(wh)) && length(unique(step))==1 &&
           length(unique(error.prob))==1 && length(unique(off.end))==1 &&
           length(unique(map.function))==1 && length(unique(ndraws))==1) {
            for(j in 1:nchr(x)) {
                geno[[j]]$draws <- array(0,dim=c(sum(n.ind),dim(x$geno[[j]]$draws)[-1]))
                dimnames(geno[[j]]$draws) <- dimnames(x$geno[[j]]$draws)
                prev <- 0
                for(i in 1:n.args) {
                    wh <- prev + 1:n.ind[i]
                    prev <- prev + n.ind[i]
                    geno[[j]]$draws[wh,,] <- args[[i]]$geno[[j]]$draws
                }

                wh <- sapply(args, function(a) match("draws",names(a$geno[[j]])))
                step <- sapply(args,function(a) attr(a$geno[[j]]$draws,"step"))
                error.prob <- sapply(args,function(a) attr(a$geno[[j]]$draws,"error.prob"))
                off.end <- sapply(args,function(a) attr(a$geno[[j]]$draws,"off.end"))
                map.function <- sapply(args,function(a) attr(a$geno[[j]]$draws,"map.function"))
                map <- lapply(args,function(a) attr(a$geno[[j]]$draws,"map"))

                attr(geno[[j]]$draws,"step") <- step[1]
                attr(geno[[j]]$draws,"error.prob") <- error.prob[1]
                attr(geno[[j]]$draws,"off.end") <- off.end[1]
                attr(geno[[j]]$draws,"map.function") <- map.function[1]
                attr(geno[[j]]$draws,"map") <- map[[1]]
            }
        }
    }

    x <- list(geno=geno, pheno=pheno)
    class(x) <- c(type,"cross")
    x
}

######################################################################
#
# fill.geno: Run argmax.geno or sim.geno and then fill in the
#            genotype data with the results.  This will allow
#            rough genome scans by marker regression without
#            holes.  WE WOULD NOT PLACE ANY TRUST IN THE RESULTS!
#
#   the newer method, "no_dbl_XO", fills in missing genotypes between
#   markers with matching genotypes
#
######################################################################

fill.geno <-
    function(cross, method=c("imp","argmax", "no_dbl_XO", "maxmarginal"),
             error.prob=0.0001, map.function=c("haldane","kosambi","c-f","morgan"),
             min.prob=0.95)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    method <- match.arg(method)

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }

    # remove any extraneous material
    cross <- clean(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    if(method=="imp") {
        # do one imputation
        temp <- sim.geno(cross,n.draws=1,step=0,off.end=0,
                         error.prob=error.prob,map.function=map.function)
        # replace the genotype data with the results,
        #     stripping off any attributes
        for(i in 1:n.chr) {
            nam <- colnames(cross$geno[[i]]$data)
            if(n.mar[i] == 1)
                cross$geno[[i]]$data <-
                    matrix(as.numeric(temp$geno[[i]]$draws[,2,1]),ncol=n.mar[i])
            else
                cross$geno[[i]]$data <-
                    matrix(as.numeric(temp$geno[[i]]$draws[,,1]),ncol=n.mar[i])
            colnames(cross$geno[[i]]$data) <- nam
        }
    }
    else if(method=="argmax") {
        # run the Viterbi algorithm
        temp <- argmax.geno(cross,step=0,off.end=0,error.prob=error.prob,
                            map.function=map.function)
        # replace the genotype data with the results,
        #     stripping off any attributes
        for(i in 1:n.chr) {
            nam <- colnames(cross$geno[[i]]$data)
            if(n.mar[i] == 1)
                cross$geno[[i]]$data <-
                    matrix(as.numeric(temp$geno[[i]]$argmax[,2]),ncol=n.mar[i])
            else
                cross$geno[[i]]$data <-
                    matrix(as.numeric(temp$geno[[i]]$argmax),ncol=n.mar[i])
            colnames(cross$geno[[i]]$data) <- nam
        }
    }
    else if(method=="maxmarginal") {
        temp <- calc.genoprob(cross, step=0, off.end=0, error.prob=error.prob,
                              map.function=map.function)

        for(i in 1:n.chr) {
            p <- temp$geno[[i]]$prob
            whmax <- apply(p, 1:2, which.max)
            maxpr <- apply(p, 1:2, max)
            g <- cross$geno[[i]]$data
            g[maxpr >= min.prob] <- whmax[maxpr > min.prob]
            g[maxpr < min.prob] <- NA
            cross$geno[[i]]$data <- g
        }
    }
    else {
        for(i in 1:n.chr) {
            dat <- cross$geno[[i]]$data
            dat[is.na(dat)] <- 0
            nr <- nrow(dat)
            nc <- ncol(dat)
            dn <- dimnames(dat)
            dat <- .C("R_fill_geno_nodblXO",
                      as.integer(nr),
                      as.integer(nc),
                      dat=as.integer(dat),
                      PACKAGE="qtl")$dat
            dat[dat==0] <- NA
            dat <- matrix(dat, ncol=nc, nrow=nr)
            dimnames(dat) <- dn
            cross$geno[[i]]$data <- dat
        }
    }
    cross
}

######################################################################
#
# checkcovar
#
# This is a utility function for scanone and scantwo.  We remove
# individuals with missing phenotypes or covariates and check
# that the covariates are of the appropriate form.
#
######################################################################

checkcovar <-
    function(cross, pheno.col, addcovar, intcovar, perm.strata, ind.noqtl=NULL, weights=NULL, verbose=TRUE)
{
    chrtype <- sapply(cross$geno, class)

    # drop individuals whose sex or pgm is missing if X chr is included
    if(any(chrtype=="X")) {
        sexpgm <- getsex(cross)
        keep <- rep(TRUE,nind(cross))
        flag <- 0
        if(!is.null(sexpgm$sex)) {
            if(any(is.na(sexpgm$sex))) {
                keep[is.na(sexpgm$sex)] <- FALSE
                flag <- 1
            }
        }
        if(!is.null(sexpgm$pgm)) {
            if(any(is.na(sexpgm$pgm))) {
                keep[is.na(sexpgm$pgm)] <- FALSE
                flag <- 1
            }
        }
        if(flag) {
            if(verbose) warning("Dropping ", sum(!keep), " individuals with missing sex or pgm.\n")
            cross <- subset(cross, ind=keep)
            if(!is.null(addcovar)) {
                if(!is.matrix(addcovar)) addcovar <- addcovar[keep]
                else addcovar <- addcovar[keep,]
            }
            if(!is.null(intcovar)) {
                if(!is.matrix(intcovar)) intcovar <- intcovar[keep]
                else intcovar <- intcovar[keep,]
            }
        }
    }

    # check phenotypes - we allow multiple phenotypes here
    if(pheno.col < 1 || pheno.col > nphe(cross))
        stop("Specified phenotype column is invalid.")
    # check if all phenotypes are numeric
    pheno <- cross$pheno[,pheno.col,drop=FALSE]
    idx.nonnum <- which(!apply(pheno,2, is.numeric))
    if(length(idx.nonnum) > 0)
        stop("Following phenotypes are not numeric: Column ",
             paste(idx.nonnum, collapse=","))

    orig.n.ind <- nind(cross)

    # drop individuals with missing phenotypes
    if(any(!is.finite(unlist(pheno)))) {
        keep.ind <- as.numeric(which(apply(pheno, 1, function(x) !any(!is.finite(x)))))
        #    keep.ind <- (1:length(pheno))[!is.na(pheno)]
        n.drop <- nind(cross) - length(keep.ind)
        keep.ind.boolean <- rep(FALSE, nind(cross))
        keep.ind.boolean[keep.ind] <- TRUE
        cross <- subset.cross(cross,ind=keep.ind.boolean)
        pheno <- pheno[keep.ind,,drop=FALSE]
        if(verbose) warning("Dropping ", n.drop, " individuals with missing phenotypes.\n")
    }
    else keep.ind <- 1:nind(cross)
    n.ind <- nind(cross)
    n.chr <- nchr(cross)      # number of chromosomes
    type <- class(cross)[1]   # type of cross

    n.addcovar <- n.intcovar <- 0
    if(!is.null(addcovar)) { # for additive covariates
        if(!is.matrix(addcovar)) {
            if(!is.numeric(as.matrix(addcovar)))
                stop("addcovar should be numeric")
            if(is.vector(addcovar) || is.data.frame(addcovar))
                addcovar <- as.matrix(addcovar)
            else stop("addcovar should be a matrix")
        }
        if(!all(apply(addcovar,2,is.numeric)))
            stop("All columns of addcovar must be numeric")
        if( nrow(addcovar) != orig.n.ind ) {
            # the length of additive covariates is incorrect
            stop("Number of rows in additive covariates is incorrect")
        }
        addcovar <- addcovar[keep.ind,,drop=FALSE]
        n.addcovar <- ncol(addcovar)

    }
    if(!is.null(intcovar)) { # interacting covariates
        if(!is.matrix(intcovar)) {
            if(!is.numeric(as.matrix(intcovar)))
                stop("intcovar should be a numeric")
            if(is.vector(intcovar) || is.data.frame(intcovar))
                intcovar <- as.matrix(intcovar)
            else stop("intcovar should be a matrix")
        }
        if(!all(apply(intcovar,2,is.numeric)))
            stop("All columns of intcovar must be numeric")
        if(nrow(intcovar)[1] != orig.n.ind) {
            # the length of interacting covariates is incorrect
            stop("The length of interacting covariates is incorrect!")
        }
        intcovar <- intcovar[keep.ind,,drop=FALSE]
        n.intcovar <- ncol(intcovar)
    }
    if(!is.null(perm.strata)) perm.strata <- perm.strata[keep.ind]
    if(!is.null(ind.noqtl)) ind.noqtl <- ind.noqtl[keep.ind]
    if(!is.null(weights)) weights <- weights[keep.ind]

    # drop individuals missing any covariates
    if(!is.null(addcovar)) { # note that intcovar is contained in addcovar
        wh <- apply(cbind(addcovar,intcovar),1,function(a) any(!is.finite(a)))
        if(any(wh)) {
            cross <- subset.cross(cross,ind=(!wh))
            pheno <- pheno[!wh,,drop=FALSE]
            addcovar <- addcovar[!wh,,drop=FALSE]
            if(!is.null(intcovar)) intcovar <- intcovar[!wh,,drop=FALSE]
            n.ind <- nind(cross)
            if(!is.null(perm.strata)) perm.strata <- perm.strata[!wh]
            if(!is.null(ind.noqtl)) ind.noqtl <- ind.noqtl[!wh]
            if(!is.null(weights)) weights <- weights[!wh]
            if(verbose) warning("Dropping ", sum(wh), " individuals with missing covariates.\n")
        }
    }

    # make sure columns of intcovar are contained in addcovar
    if(!is.null(intcovar)) {
        if(is.null(addcovar)) {
            addcovar <- intcovar
            n.addcovar <- n.intcovar
            if(verbose) warning("addcovar forced to contain all columns of intcovar\n")
        }
        else {
            wh <- 1:n.intcovar
            for(i in 1:n.intcovar) {
                o <- (apply(addcovar,2,function(a,b) max(abs(a-b)),intcovar[,i])<1e-14)
                if(any(o)) wh[i] <- (1:n.addcovar)[o]
                else wh[i] <- NA
            }
            if(any(!is.finite(wh))) {
                addcovar <- cbind(addcovar,intcovar[,!is.finite(wh)])
                n.addcovar <- ncol(addcovar)
                if(verbose) warning("addcovar forced to contain all columns of intcovar")
            }
        }
    }

    if(n.addcovar > 0) { # check rank
        if(qr(cbind(1,addcovar))$rank < n.addcovar+1)
            if(verbose) warning("addcovar appears to be over-specified; consider dropping columns.\n")
    }
    if(n.intcovar > 0) { # check rank
        if(qr(cbind(1,intcovar))$rank < n.intcovar+1)
            if(verbose) warning("intcovar appears to be over-specified; consider dropping columns.\n")
    }

    pheno <- as.matrix(pheno)

    list(cross=cross, pheno=pheno, addcovar=addcovar, intcovar=intcovar,
         n.addcovar=n.addcovar, n.intcovar=n.intcovar, perm.strata=perm.strata,
         ind.noqtl=ind.noqtl, weights=weights)
}

# Find the nearest marker to a particular position
find.marker <-
    function(cross, chr, pos, index)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(missing(pos) && missing(index))
        stop("Give either pos or index.")
    if(!missing(pos) && !missing(index))
        stop("Give just one of pos or index.")

    # if chr has length 1, expand if necessary
    if(length(chr) == 1) {
        if(!missing(pos))
            chr <- rep(chr,length(pos))
        else
            chr <- rep(chr, length(index))
    }
    # otherwise, chr and pos should have same length
    else {
        if(!missing(pos) && length(chr) != length(pos))
            stop("chr and pos must be the same length.")
        if(!missing(index) && length(chr) != length(index))
            stop("chr and index must be the same length.")
    }

    markers <- rep("",length(chr))
    chrnotfound <- NULL
    for(i in 1:length(chr)) {
        # find chromosome
        o <- match(chr[i], names(cross$geno))
        if(is.na(o)) {
            markers[i] <- NA  # chr not matched
            chrnotfound <- c(chrnotfound, chr[i])
        }
        else {
            thismap <- cross$geno[[o]]$map # genetic map
            # sex-specific map; look at female positions
            if(is.matrix(thismap)) thismap <- thismap[1,]

            if(!missing(pos)) {
                # find closest marker
                d <- abs(thismap-pos[i])
                o2 <- (1:length(d))[d==min(d)]
                if(length(o2)==1) markers[i] <- names(thismap)[o2]
                # if multiple markers are equidistant,
                #     choose the one with the most data
                #     or choose among them at random
                else {
                    x <- names(thismap)[o2]
                    n.geno <- apply(cross$geno[[o]]$data[,o2],2,function(a) sum(!is.na(a)))

                    o2 <- o2[n.geno==max(n.geno)]
                    if(length(o2) == 1)
                        markers[i] <- names(thismap)[o2]
                    else
                        markers[i] <- names(thismap)[sample(o2,1)]

                }
            }
            else { # by index
                if(index[i] < 1 || index[i] > length(thismap))
                    stop("Misspecified index ", index[i], " on chr ", chr[i])
                markers[i] <- names(thismap)[index[i]]
            }
        }
    }
    if(length(chrnotfound) > 0) {
        chrnotfound <- sort(unique(chrnotfound))
        if(length(chrnotfound) == 1)
            warning("Chromosome ", paste("\"", chrnotfound, "\"", sep=""), " not found")
        else
            warning("Chromosomes ", paste("\"", chrnotfound, "\"", sep="", collapse=", "), " not found")
    }

    markers
}


### Find the nearest pseudomarker to a particular position
find.pseudomarker <-
    function(cross, chr, pos, where=c("draws","prob"), addchr=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    # if chr has length 1, expand if necessary
    if(length(chr) == 1)
        chr <- rep(chr,length(pos))
    # otherwise, chr and pos should have same length
    else if(length(chr) != length(pos))
        stop("chr and pos must be the same length.")

    markers <- rep("",length(chr))

    where <- match.arg(where)
    if(where=="draws" && !("draws" %in% names(cross$geno[[1]])))
        stop("You'll need to first run sim.geno")
    if(where=="prob" && !("prob" %in% names(cross$geno[[1]])))
        stop("You'll need to first run calc.genoprob")

    for(i in 1:length(chr)) {
        # find chromosome
        o <- match(chr[i], names(cross$geno))
        if(is.na(o)) markers[i] <- NA  # chr not matched
        else {
            if(where=="draws") {
                if("map" %in% names(attributes(cross$geno[[o]]$draws)))
                    thismap <- attr(cross$geno[[o]]$draws, "map")
                else {
                    stp <- attr(cross$geno[[o]]$draws, "step")
                    oe <- attr(cross$geno[[o]]$draws, "off.end")

                    if("stepwidth" %in% names(attributes(cross$geno[[o]]$draws)))
                        stpw <- attr(cross$geno[[o]]$draws, "stepwidth")
                    else stpw <- "fixed"
                    thismap <- create.map(cross$geno[[o]]$map,stp,oe,stpw)
                    attr(cross$geno[[o]]$draws, "map") <- thismap
                }
            }
            else { # prob
                if("map" %in% names(attributes(cross$geno[[o]]$prob)))
                    thismap <- attr(cross$geno[[o]]$prob, "map")
                else {
                    stp <- attr(cross$geno[[o]]$prob, "step")
                    oe <- attr(cross$geno[[o]]$prob, "off.end")

                    if("stepwidth" %in% names(attributes(cross$geno[[o]]$prob)))
                        stpw <- attr(cross$geno[[o]]$prob, "stepwidth")
                    else stpw <- "fixed"
                    thismap <- create.map(cross$geno[[o]]$map,stp,oe,stpw)
                    attr(cross$geno[[o]]$prob, "map") <- thismap
                }
            }

            # sex-specific map; look at female positions
            if(is.matrix(thismap)) thismap <- thismap[1,]

            # find closest marker
            d <- abs(thismap-pos[i])
            o2 <- (1:length(d))[d==min(d)]
            if(length(o2)==1) themarker <- names(thismap)[o2]
            else themarker <- names(thismap)[sample(o2, 1)]

            if(addchr && length(grep("^loc[0-9]+\\.*[0-9]*(\\.[0-9]+)*$", themarker)) > 0)
                themarker <- paste("c", chr[i], ".", themarker, sep="")
            markers[i] <- themarker
        }
    }

    markers
}



# expand recombination fractions for RI lines
adjust.rf.ri <-
    function(r, type=c("self","sib"), chrtype=c("A","X"), expand=TRUE)
{
    # type of RI lines
    type <- match.arg(type)
    chrtype <- match.arg(chrtype)

    if(type=="self") {
        if(expand) return(r*2/(1+2*r))
        else return(r/2/(1-r))
    }
    else {
        if(chrtype=="A") { # autosome / sib mating
            if(expand) return(r*4/(1+6*r))
            else return(r/(4-6*r))
        }
        else { # X chromosome/ sib mating
            if(expand) return(8/3*r/(1+4*r))
            else return(3/8*r/(1-1.5*r))
        }
    }
}

######################################################################
# lodint: function to get lod support interval
######################################################################
lodint <-
    function(results, chr, qtl.index, drop=1.5, lodcolumn=1,
             expandtomarkers=FALSE)
{
    if(!("scanone" %in% class(results))) {
        if(!("qtl" %in% class(results)))
            stop("Input must have class \"scanone\" or \"qtl\".")
        else {
            if(!("lodprofile" %in% names(attributes(results))))
                stop("qtl object needs to be produced by refineqtl with keeplodprofile=TRUE.")
            else { # qtl object
                if(lodcolumn != 1) {
                    warning("lod column ignored if input is a qtl object.")
                    lodcolumn <- 1
                }
                results <- attr(results, "lodprofile")
                if(missing(qtl.index)) {
                    if(length(results)==1)
                        results <- results[[1]]
                    else
                        stop("You must specify qtl.index.")
                }
                else {
                    if(length(qtl.index)>1) stop("qtl.index should have length 1")
                    if(qtl.index < 1 || qtl.index > length(results))
                        stop("qtl.index misspecified.")
                    results <- results[[qtl.index]]
                }
                chr <- results[1,1]
            }
        }
    }
    else {
        if(lodcolumn < 1 || lodcolumn +2 > ncol(results))
            stop("Argument lodcolumn misspecified.")

        if(missing(chr)) {
            if(length(unique(results[,1]))>1)
                stop("Give a chromosome ID.")
        }
        else {
            if(length(chr) > 1) stop("chr should have length 1")
            if(is.na(match(chr, results[,1])))
                stop("Chromosome misspecified.")
            results <- results[results[,1]==chr,]
        }
    }

    if(all(is.na(results[,lodcolumn+2]))) return(NULL)

    maxlod <- max(results[,lodcolumn+2],na.rm=TRUE)
    w <- which(!is.na(results[,lodcolumn+2]) & results[,lodcolumn+2] == maxlod)
    o <- range(which(!is.na(results[,lodcolumn+2]) & results[,lodcolumn+2] > maxlod-drop))

    if(length(o)==0) o <- c(1,nrow(results))

    else {
        if(o[1] > 1) o[1] <- o[1]-1
        if(o[2] < nrow(results)) o[2] <- o[2]+1
    }

    if(expandtomarkers) {
        markerpos <- (1:nrow(results))[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", rownames(results))]
        if(any(markerpos <= o[1]))
            o[1] <- max(markerpos[markerpos <= o[1]])
        if(any(markerpos >= o[2]))
            o[2] <- min(markerpos[markerpos >= o[2]])
    }

    rn <- rownames(results)[c(o[1],w,o[2])]

    # look for duplicate rows
    if(length(w)>1 && rn[length(rn)]==rn[length(rn)-1]) w <- w[-length(w)]
    else if(length(w)>1 && rn[2]==rn[1]) w <- w[-1]
    rn <- rownames(results)[c(o[1],w,o[2])]

    # look for more duplicate rows
    if(any(table(rn)> 1)) {
        tab <- table(rn)
        temp <- which(tab>1)
        for(j in temp) {
            z <- which(rn==names(tab)[j])
            for(k in 2:length(z))
                rn[z[k:length(z)]] <- paste(rn[z[k:length(z)]], " ", sep="")
        }
    }

    results <- results[c(o[1],w,o[2]),]
    rownames(results) <- rn
    class(results) <- c("scanone","data.frame")

    results
}


######################################################################
# bayesint: function to get Bayesian probability interval
######################################################################
bayesint <-
    function(results, chr, qtl.index, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
{
    if(!("scanone" %in% class(results))) {
        if(!("qtl" %in% class(results)))
            stop("Input must have class \"scanone\" or \"qtl\".")
        else {
            if(!("lodprofile" %in% names(attributes(results))))
                stop("qtl object needs to be produced by refineqtl with keeplodprofile=TRUE.")
            else { # qtl object
                if(lodcolumn != 1) {
                    warning("lod column ignored if input is a qtl object.")
                    lodcolumn <- 1
                }
                results <- attr(results, "lodprofile")
                if(missing(qtl.index)) {
                    if(length(results)==1)
                        results <- results[[1]]
                    else
                        stop("You must specify qtl.index.")
                }
                else {
                    if(length(qtl.index)>1) stop("qtl.index should have length 1")
                    if(qtl.index < 1 || qtl.index > length(results))
                        stop("qtl.index misspecified.")
                    results <- results[[qtl.index]]
                }
                chr <- results[1,1]
            }
        }
    }
    else {
        if(lodcolumn < 1 || lodcolumn +2 > ncol(results))
            stop("Argument lodcolumn misspecified.")

        if(missing(chr)) {
            if(length(unique(results[,1]))>1)
                stop("Give a chromosome ID.")
        }
        else {
            if(length(chr) > 1) stop("chr should have length 1")
            if(is.na(match(chr, results[,1])))
                stop("Chromosome misspecified.")
            results <- results[results[,1]==chr,]
        }
    }

    if(all(is.na(results[,lodcolumn+2]))) return(NULL)

    loc <- results[,2]
    width <- diff(( c(loc[1],loc) + c(loc, loc[length(loc)]) )/ 2)

    area <- 10^results[,lodcolumn+2]*width
    area <- area/sum(area)

    o <- order(results[,lodcolumn+2], decreasing=TRUE)

    cs <- cumsum(area[o])

    wh <- min((1:length(loc))[cs >= prob])
    int <- range(o[1:wh])

    if(expandtomarkers) {
        markerpos <- (1:nrow(results))[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", rownames(results))]
        if(any(markerpos <= int[1]))
            int[1] <- max(markerpos[markerpos <= int[1]])
        if(any(markerpos >= int[2]))
            int[2] <- min(markerpos[markerpos >= int[2]])
    }

    rn <- rownames(results)[c(int[1],o[1],int[2])]
    # look for duplicate rows
    if(any(table(rn)> 1)) {
        rn[2] <- paste(rn[2], "")
        if(rn[1] == rn[3]) rn[3] <- paste(rn[3], " ")
    }

    results <- results[c(int[1],o[1],int[2]),]
    rownames(results) <- rn
    class(results) <- c("scanone", "data.frame")
    results
}


######################################################################
# makeSSmap: convert a genetic map, or the genetic maps in a cross
#            object, to be sex-specific (i.e., 2-row matrices)
######################################################################
makeSSmap <-
    function(cross)
{
    if(!any(class(cross) == "map")) {
        # input object is a genetic map
        for(i in 1:length(cross)) {
            if(!is.matrix(cross[[i]]))
                cross[[i]] <- rbind(cross[[i]], cross[[i]])
        }
    }
    else { # input object is assumed to be a "cross" object
        n.chr <- nchr(cross)
        for(i in 1:n.chr) {
            if(!is.matrix(cross$geno[[i]]$map))
                cross$geno[[i]]$map <-
                    rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
        }
    }

    cross
}

######################################################################
# comparecrosses: verify that two cross objects have identical
#                 classes, chromosomes, markers, genotypes, maps,
#                 and phenotypes
######################################################################
comparecrosses <-
    function(cross1, cross2, tol=1e-5)
{
    if(missing(cross1) || missing(cross2))
        stop("Two crosses must be input.")

    # both are of class "cross"
    if(!any(class(cross1) == "cross") ||
       !any(class(cross2) == "cross"))
        stop("Input should have class \"cross\".")

    # classes are the same
    if(any(class(cross1) != class(cross2)))
        stop("crosses are not the same type.")

    if(nchr(cross1) != nchr(cross2))
        stop("crosses do not have the same number of chromosomes.")

    if(any(names(cross1$geno) != names(cross2$geno)))
        stop("Chromosome names do not match.")

    if(any(nmar(cross1) != nmar(cross2)))
        stop("Number of markers per chromosome do not match.")

    mnames1 <- unlist(lapply(cross1$geno, function(a) colnames(a$data)))
    mnames2 <- unlist(lapply(cross2$geno, function(a) colnames(a$data)))
    if(any(mnames1 != mnames2)) {
        #    stop("Markers names do not match.")
        for(i in 1:nchr(cross1))
            if(any(colnames(cross1$geno[[i]]$data) != colnames(cross2$geno[[i]]$data)))
                stop("Marker names on chr ", names(cross1$geno)[i], " don't match.")
    }


    chrtype1 <- sapply(cross1$geno, class)
    chrtype2 <- sapply(cross2$geno, class)
    if(any(chrtype1 != chrtype2))
        stop("Chromosome types (autosomal vs X) do not match.")

    for(i in 1:nchr(cross1)) {
        if(any(abs(diff(cross1$geno[[i]]$map) - diff(cross2$geno[[i]]$map)) > tol))
            stop("Genetic maps for chromosome ", names(cross1$geno)[i],
                 " do not match.")

        if(abs(cross1$geno[[i]]$map[1] - cross2$geno[[i]]$map[1]) > tol)
            warning("Initial marker positions for chromosome ", names(cross1$geno)[i],
                    " do not match.")
    }

    if(nind(cross1) != nind(cross2))
        stop("Number of individuals do not match.")

    for(i in 1:nchr(cross1)) {
        g1 <- cross1$geno[[i]]$data
        g2 <- cross2$geno[[i]]$data
        if(any((is.na(g1) & !is.na(g2)) | (!is.na(g1) & is.na(g2)) |
               (!is.na(g1) & !is.na(g2) & g1!=g2)))
            stop("Genotype data for chromosome ", names(cross1$geno)[i],
                 " do not match.")
    }

    if(nphe(cross1) != nphe(cross2))
        stop("Number of phenotypes do not match.")

    if(any(names(cross1$pheno) != names(cross2$pheno)))
        stop("Phenotype names do not match.")

    for(i in 1:nphe(cross1)) {
        phe1 <- cross1$pheno[,i]
        phe2 <- cross2$pheno[,i]
        if(is.numeric(phe1) & is.numeric(phe2)) {
            phe1[phe1 == Inf] <- max(phe1[phe1 < Inf], na.rm=TRUE)+5
            phe2[phe2 == Inf] <- max(phe2[phe2 < Inf], na.rm=TRUE)+5
            phe1[phe1 == -Inf] <- min(phe1[phe1 > -Inf], na.rm=TRUE)-5
            phe2[phe2 == -Inf] <- min(phe2[phe2 > -Inf], na.rm=TRUE)-5
            if(any((is.na(phe1) & !is.na(phe2)) | (!is.na(phe1) & is.na(phe2)) |
                   (!is.na(phe1) & !is.na(phe2) & abs(phe1-phe2) > tol))) {
                stop("Data for phenotype ", names(cross1$pheno)[i],
                     " do not match.")
            }
        }
        else {
            if(any((is.na(phe1) & !is.na(phe2)) | (!is.na(phe1) & is.na(phe2)) |
                   (!is.na(phe1) & !is.na(phe2) &
                    as.character(phe1) != as.character(phe2)))) {
                stop("Data for phenotype ", names(cross1$pheno)[i], " do not match.")
            }
        }
    }

    message("\tCrosses are identical.")
}


######################################################################
# move marker
# Move a marker to a new chromosome...placed at the end
######################################################################
movemarker <-
    function(cross, marker, newchr, newpos)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    mnames <- unlist(lapply(cross$geno,function(a) colnames(a$data)))
    chr <- rep(names(cross$geno),nmar(cross))
    pos <- unlist(lapply(cross$geno,function(a) 1:ncol(a$data)))
    oldindex <- match(marker, mnames)

    # Marker found precisely once?
    if(is.na(oldindex)) stop(marker, " not found.\n")
    if(length(oldindex) > 1) stop(marker, " found multiple times.\n")

    chr <- chr[oldindex]
    pos <- pos[oldindex]

    # pull out genotype data
    g <- cross$geno[[chr]]$data[,pos]

    chrtype <- class(cross$geno[[chr]])
    mapmatrix <- is.matrix(cross$geno[[chr]]$map)

    # delete marker
    if(nmar(cross)[chr] == 1)  { # only marker on that chromosome, so drop the chromosome
        cross$geno <- cross$geno[-match(chr,names(cross$geno))]
        delchr <- TRUE
    }
    else {
        delchr <- FALSE
        cross$geno[[chr]]$data <- cross$geno[[chr]]$data[,-pos,drop=FALSE]
        if(is.matrix(cross$geno[[chr]]$map))
            cross$geno[[chr]]$map <- cross$geno[[chr]]$map[,-pos,drop=FALSE]
        else
            cross$geno[[chr]]$map <- cross$geno[[chr]]$map[-pos]
    }

    if(is.numeric(newchr)) newchr <- as.character(newchr)

    if(!(newchr %in% names(cross$geno))) { # create a new chromosome
        n <- length(cross$geno)
        cross$geno[[n+1]] <- list("data"=as.matrix(g),
                                  "map"=as.numeric(0))
        names(cross$geno)[n+1] <- newchr
        class(cross$geno[[n+1]]) <- chrtype
        colnames(cross$geno[[n+1]]$data) <- marker
        if(mapmatrix) {
            if(missing(newpos)) newpos <- 0
            cross$geno[[n+1]]$map <- matrix(newpos, ncol=1, nrow=2)
            colnames(cross$geno[[n+1]]$map) <- marker
        }
        else {
            if(missing(newpos)) newpos <- 0
            cross$geno[[n+1]]$map <- newpos
            names(cross$geno[[n+1]]$map) <- marker
        }
        return(cross)
    }

    if(missing(newpos)) {
        # add marker to end of new chromosome
        n.mar <- nmar(cross)[newchr]

        cross$geno[[newchr]]$data <- cbind(cross$geno[[newchr]]$data,g)
        colnames(cross$geno[[newchr]]$data)[n.mar+1] <- marker

        if(is.matrix(cross$geno[[newchr]]$map)) {
            cross$geno[[newchr]]$map <- cbind(cross$geno[[newchr]]$map,
                                              cross$geno[[newchr]]$map[,n.mar]+10)
            colnames(cross$geno[[newchr]]$map)[n.mar+1] <- marker
        }
        else {
            cross$geno[[newchr]]$map <- c(cross$geno[[newchr]]$map,
                                          cross$geno[[newchr]]$map[n.mar]+10)
            names(cross$geno[[newchr]]$map)[n.mar+1] <- marker
        }
    }
    else {
        # add marker to the specified position
        dat <- cross$geno[[newchr]]$data
        map <- cross$geno[[newchr]]$map

        if(length(newpos) != 1)
            stop("newpos should be a single number.")

        if(is.matrix(map)) { # sex-specific maps
            wh <- which(map[1,] < newpos)
            if(length(wh) == 0) { # place in first spot
                map <- cbind(c(newpos,map[2,1]-(map[1,1]-newpos)),map)
                colnames(map)[1] <- marker
            }
            else {
                wh <- max(wh)
                if(wh == ncol(map)) { # place at end of chromosome
                    map <- cbind(map,c(newpos,map[2,ncol(map)]+(newpos-map[1,ncol(map)])))
                    colnames(map)[ncol(map)] <- marker
                }
                else {
                    left <- map[,1:wh,drop=FALSE]
                    right <- map[,-(1:wh),drop=FALSE]
                    marleft <- colnames(left)
                    marright <- colnames(right)
                    left <- left[,ncol(left)]
                    right <- right[,1]

                    newpos2 <- (newpos-left[1])/(right[1]-left[1])*(right[2]-left[2])+left[2]
                    map <- cbind(map[,1:wh], c(newpos,newpos2), map[,-(1:wh)])
                    colnames(map) <- c(marleft, marker, marright)
                }
            }
        }
        else {
            wh <- which(map < newpos)
            if(length(wh) == 0) { # place in first position
                map <- c(newpos,map)
                names(map)[1] <- marker
            }
            else {
                wh <- max(wh)
                if(wh == length(map)) { # place in last position
                    map <- c(map,newpos)
                    names(map)[length(map)] <- marker
                }
                else {
                    map <- c(map[1:wh],newpos,map[-(1:wh)])
                    names(map)[wh+1] <- marker
                }
            }
        }
        cross$geno[[newchr]]$map <- map

        if(length(wh)==0) { # place marker in first position
            dat <- cbind(g, dat)
            colnames(dat)[1] <- marker
        }
        else if(wh == ncol(dat)) { # place marker in last position
            dat <- cbind(dat, g)
            colnames(dat)[ncol(dat)] <- marker
        }
        else { # place marker in the middle
            dat <- cbind(dat[,1:wh],g,dat[,-(1:wh)])
            colnames(dat)[wh+1] <- marker
        }
        cross$geno[[newchr]]$data <- dat

        # make sure the marker names for the data and the genetic map match
        if(is.matrix(cross$geno[[newchr]]$map))
            colnames(cross$geno[[newchr]]$data) <- colnames(cross$geno[[newchr]]$map)
        else
            colnames(cross$geno[[newchr]]$data) <- names(cross$geno[[newchr]]$map)
    }

    # update genoprob, errorlod, argmax, draws, rf
    if("rf" %in% names(cross)) {
        # reorder the recombination fractions
        #  -- a bit of pain, 'cause we need LODs in upper triangle
        #     and rec fracs in lower triangle
        newmar <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

        attrib <- attributes(cross$rf)

        rf <- cross$rf
        lods <- rf;lods[lower.tri(rf)] <- t(rf)[lower.tri(rf)]
        rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
        lods <- lods[newmar,newmar]
        rf <- rf[newmar,newmar]
        rf[upper.tri(rf)] <- lods[upper.tri(rf)]

        cross$rf <- rf

        if("onlylod" %in% names(attrib)) # save the onlylod attribute if its there
            attr(cross$rf, "onlylod") <- attrib$onlylod
    }


    if(!delchr) thechr <- c(chr,newchr)
    else thechr <- newchr
    for(i in thechr) {
        tempg <- cross$geno[[i]]
        tempx <- subset(cross, chr=i)
        if("prob" %in% names(tempg))
            atp <- attributes(tempg$prob)
        if("draws" %in% names(tempg)) {
            at <- attributes(cross$geno[[1]]$draws)
            tempg$draws <- sim.geno(tempx,
                                    n.draws=at$dim[3],
                                    step=at$step,
                                    off.end=at$off.end,
                                    map.function=at$map.function,
                                    error.prob=at$error.prob)$geno[[1]]$draws
        }
        if("argmax" %in% names(tempg)) {
            at <- attributes(cross$geno[[1]]$argmax)
            tempg$argmax <- argmax.geno(tempx,
                                        step=at$step,
                                        off.end=at$off.end,
                                        map.function=at$map.function,
                                        error.prob=at$error.prob)$geno[[1]]$argmax
        }
        if("errorlod" %in% names(tempg)) {
            at <- attributes(cross$geno[[1]]$errorlod)
            tempg$errorlod <- argmax.geno(tempx,
                                          map.function=at$map.function,
                                          error.prob=at$error.prob)$geno[[1]]$errorlod
        }
        if("prob" %in% names(tempg))
            tempg$prob <- calc.genoprob(tempx,
                                        step=atp$step,
                                        off.end=atp$off.end,
                                        map.function=atp$map.function,
                                        error.prob=atp$error.prob)$geno[[1]]$prob

        cross$geno[[i]] <- tempg
    }

    cross
}

######################################################################
#
# summary.map
#
# Give a short summary of a genetic map object.
#
######################################################################
summaryMap <- summary.map <-
    function(object, ...)
{
    map <- object
    if(any(class(map) == "cross")) # a cross object
        map <- pull.map(map)
    if(!any(class(map) == "map"))
        stop("Input should have class \"cross\" or \"map\".")

    n.chr <- length(map)
    chrnames <- names(map)
    if(is.null(chrnames)) chrnames <- 1:length(map)
    if(is.matrix(map[[1]])) { # sex-specific map
        sexsp <- TRUE
        n.mar <- sapply(map,ncol)
        tot.mar <- sum(n.mar)
        fmap <- lapply(map,function(a) a[1,])
        mmap <- lapply(map,function(a) a[2,])

        len.f <- sapply(fmap,function(a) diff(range(a)))
        len.m <- sapply(mmap,function(a) diff(range(a)))
        avesp.f <- sapply(fmap,function(a) {if(length(a)<2) return(NA); mean(diff(a))})
        avesp.m <- sapply(mmap,function(a) {if(length(a)<2) return(NA); mean(diff(a))})
        maxsp.f <- sapply(fmap,function(a) {if(length(a)<2) return(NA); max(diff(a))})
        maxsp.m <- sapply(mmap,function(a) {if(length(a)<2) return(NA); max(diff(a))})
        totlen.f <- sum(len.f)
        totlen.m <- sum(len.m)

        tot.avesp.f <- mean(unlist(lapply(fmap,diff)), na.rm=TRUE)
        tot.avesp.m <- mean(unlist(lapply(mmap,diff)), na.rm=TRUE)
        tot.maxsp.f <- max(maxsp.f,na.rm=TRUE)
        tot.maxsp.m <- max(maxsp.m,na.rm=TRUE)

        output <- rbind(cbind(n.mar,len.f,len.m,avesp.f,avesp.m, maxsp.f, maxsp.m),
                        c(tot.mar,totlen.f,totlen.m,tot.avesp.f,tot.avesp.m,
                          tot.maxsp.f, tot.maxsp.m))
        dimnames(output) <- list(c(chrnames,"overall"),
                                 c("n.mar","length.female","length.male",
                                   "ave.spacing.female","ave.spacing.male",
                                   "max.spacing.female", "max.spacing.male"))
    }
    else {
        sexsp=FALSE
        n.mar <- sapply(map,length)
        len <- sapply(map,function(a) diff(range(a)))
        tot.mar <- sum(n.mar)

        len <- sapply(map,function(a) diff(range(a)))
        avesp <- sapply(map,function(a) {if(length(a)<2) return(NA); mean(diff(a))})
        maxsp <- sapply(map,function(a) {if(length(a)<2) return(NA); max(diff(a))})
        totlen <- sum(len)
        tot.avesp <- mean(unlist(lapply(map,diff)), na.rm=TRUE)
        tot.maxsp <- max(maxsp, na.rm=TRUE)

        output <- rbind(cbind(n.mar,len,avesp, maxsp),
                        c(tot.mar,totlen,tot.avesp, tot.maxsp))
        dimnames(output) <- list(c(chrnames,"overall"),
                                 c("n.mar","length","ave.spacing", "max.spacing"))
    }

    output <- as.data.frame(output)
    attr(output, "sexsp") <- sexsp
    class(output) <- c("summary.map", "data.frame")
    output
}


######################################################################
#
# print.summary.map
#
# Print out the result of summary.map()
#
######################################################################
print.summary.map <-
    function(x, ...)
{
    sexsp <- attr(x, "sexsp")
    if(sexsp) cat("Sex-specific map\n\n")

    x <- apply(x,2,round,1)
    print(x)
}

######################################################################
# convert functions
######################################################################
convert <-
    function(object, ...)
    UseMethod("convert")


######################################################################
#
# convert.scanone
#
# Convert scanone output from the format for R/qtl ver 0.97 to
# that for R/qtl ver 0.98
# (previously, inter-maker locations named loc*.c*; now c*.loc*)
#
######################################################################
convert.scanone <-
    function(object, ...)
{
    if(!any(class(object) == "scanone"))
        stop("Input should have class \"scanone\".")

    rn <- rownames(object)
    o <- grep("^loc-*[0-9]+(\\.[0-9]+)*\\.c[0-9A-Za-z]+$", rn)
    if(length(o) > 0) {
        temp <- rn[o]
        temp <- strsplit(temp,"\\.")
        temp <- sapply(temp, function(a)
                       paste(a[c(length(a),1:(length(a)-1))],collapse="."))
        rownames(object)[o] <- temp
    }
    object
}

######################################################################
# convert.scantwo
#
# convert scantwo output from the format used in R/qtl version 1.03
# and earlier to that used in R/qtl version 1.04 and later.
#
# 1.03 and earlier: contained joint and interactive LOD scores
# 1.04 and later:   contains joint LOD scores and LOD scores from the
#                   additive QTL model
######################################################################
convert.scantwo <-
    function(object, ...)
{
    if(!any(class(object) == "scantwo"))
        stop("Input should have class \"scantwo\".")

    lod <- object$lod
    if(length(dim(lod)) == 2) {
        u <- upper.tri(lod)
        lod[u] <- t(lod)[u] - lod[u]
    }
    else { # multiple phenotypes
        u <- upper.tri(lod[,,1])
        for(i in 1:dim(lod)[3])
            lod[,,i][u] <- t(lod[,,i])[u] - lod[,,i][u]
    }
    object$lod <- lod
    object
}


######################################################################
# convert.map
#
# convert a genetic map from one map function to another
######################################################################
convert.map <-
    function(object, old.map.function=c("haldane", "kosambi", "c-f", "morgan"),
             new.map.function=c("haldane", "kosambi", "c-f", "morgan"), ...)
{
    old.map.function <- match.arg(old.map.function)
    new.map.function <- match.arg(new.map.function)
    if(!("map" %in% class(object)))
        stop("Input should have class \"map\".")

    if(old.map.function==new.map.function) {
        warning("old and new map functions are the same; no change.")
        return(object)
    }

    mf <- switch(old.map.function,
                 "haldane"=mf.h,
                 "kosambi"=mf.k,
                 "c-f"=mf.cf,
                 "morgan"=mf.m)
    imf <- switch(new.map.function,
                  "haldane"=imf.h,
                  "kosambi"=imf.k,
                  "c-f"=imf.cf,
                  "morgan"=imf.m)

    if(is.matrix(object[[1]])) { # sex-specific map
        for(i in seq(along=object)) {
            theclass <- class(object[[i]])
            thenames <- colnames(object[[i]])
            for(j in 1:2)
                object[[i]][j,] <- cumsum(c(object[[i]][j,1], imf(mf(diff(object[[i]][j,])))))

            class(object[[i]]) <- theclass
            colnames(object[[i]]) <- thenames
        }
    }
    else {
        for(i in seq(along=object)) {
            theclass <- class(object[[i]])
            thenames <- names(object[[i]])
            object[[i]] <- cumsum(c(object[[i]][1], imf(mf(diff(object[[i]])))))
            class(object[[i]]) <- theclass
            names(object[[i]]) <- thenames
        }
    }

    object
}


######################################################################
# find.pheno
#
# utility to get pheno number given pheno name
######################################################################
find.pheno <-
    function( cross,  pheno )
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    seq( ncol( cross$pheno ))[match(pheno,names(cross$pheno))]
}

######################################################################
# find.flanking
#
# utility to get flanking and/or closest marker to chr and pos
######################################################################
find.flanking <-
    function( cross, chr, pos)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    map = pull.map(cross)

    if(is.matrix(map[[1]]) && nrow(map[[1]]) > 1)
        stop("This function works only for crosses with sex-averaged maps.")

    if(length(chr) == 1 && length(pos) > 1) {
        chr <- rep(chr,length(pos))
    }
    wh <- match(chr, names(cross$geno))
    if(any(is.na(wh))) {
        stop("Chr ", paste(chr[is.na(wh)], collapse=", "), " not found.")
        wh <- wh[!is.na(wh)]
    }
    chr <- names(cross$geno)[wh]

    marker = NULL
    for (i in seq(length(chr))) {
        tmp = map[[chr[i]]]-pos[i]
        m = names(map[[chr[i]]])
        left = sum(tmp < 0)
        at = sum(tmp == 0)
        right = sum(tmp > 0)
        f <- if (at > 0)
            left+at[c(1,length(at))]
        else {
            if (right > 0)
                c(left,left+at+1)
            else
                c(left,left+at)
        }
        marker = rbind(marker,m[f[c(1:2,order(abs(tmp[f]))[1])]])
    }
    dimnames(marker) <- list(paste(chr,":",pos,sep=""),
                             c("left","right","close"))
    as.data.frame(marker, stringsAsFactors=TRUE)
}

######################################################################
# strip.partials
#
# Removes partially informative genotypes in an intercross.
#
# Input: a cross object; if verbose=TRUE, a statement regarding the
#        number of genotypes removed is printed.
######################################################################
strip.partials <-
    function(cross, verbose=TRUE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    type <- class(cross)[1]
    if(type != "f2")
        stop("This is for intercrosses only")

    n.removed <- 0
    for(i in 1:nchr(cross)) {
        g <- cross$geno[[i]]$data
        wh <- !is.na(g) & g>3
        if(any(wh)) {
            g[wh] <- NA
            cross$geno[[i]]$data <- g
            n.removed <- n.removed + sum(wh)
        }
    }
    if(verbose) {
        if(n.removed == 0) cat(" --Didn't remove any genotypes.\n")
        else cat("Removed", n.removed, "genotypes.\n")
    }
    cross
}

######################################################################
# comparegeno
######################################################################
comparegeno <-
    function(cross, what=c("proportion","number", "both"))
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    what <- match.arg(what)
    g <- pull.geno(cross)
    g[is.na(g)] <- 0
    n.ind <- nrow(g)
    n.mar <- ncol(g)
    z <- .C("R_comparegeno",
            as.integer(g),
            as.integer(n.ind),
            as.integer(n.mar),
            n.match=as.integer(rep(0,n.ind^2)),
            n.missing=as.integer(rep(0,n.ind^2)),
            PACKAGE="qtl")

    if(what=="number") {
        z <- matrix(z$n.match,n.ind,n.ind)
    }
    else {
        if(what=="proportion") {
            z <- matrix(z$n.match/(n.mar-z$n.missing),n.ind,n.ind)
            diag(z) <- NA
        }
        else {
            prop <- matrix(z$n.match/(n.mar-z$n.missing),n.ind,n.ind)
            z <- matrix(z$n.match,n.ind,n.ind)
            z[lower.tri(z)] <- prop[lower.tri(z)]
        }
    }
    z
}


######################################################################
# print the installed version of R/qtl
######################################################################
qtlversion <-
    function()
{
    version <- unlist(packageVersion("qtl"))

    # make it like #.#-#
    paste(c(version,".","-")[c(1,4,2,5,3)], collapse="")
}


######################################################################
#
# locateXO
#
# Locate crossovers on a single chromosome in each individual
# Look at just the first chromosome
#
######################################################################

locateXO <-
    function(cross, chr, full.info=FALSE)
{
    if(!missing(chr)) {
        cross <- subset(cross, chr=chr)
        if(nchr(cross) != 1)
            warning("locateXO works on just one chr; considering chr ", names(cross$geno)[1])
    }

    # individual IDs
    id <- getid(cross)
    if(is.null(id)) id <- as.character(1:nind(cross))

    if(nmar(cross)[1] == 1) { # just one marker; don't need to do anything
        warning("Just one marker.")
        res <- vector("list", nind(cross))
        names(res) <- id
        for(i in seq(along=res)) res[[i]] <- numeric(0)
        return(res)
    }

    geno <- cross$geno[[1]]$data
    geno[is.na(geno)] <- 0
    type <- class(cross)[1]

    if(type != "bc" && type != "f2" && type != "riself" && type != "risib" && type!="dh" && type!="haploid")
        stop("locateXO only working for backcross, intercross or RI strains.")

    map <- cross$geno[[1]]$map
    if(is.matrix(map)) map <- map[1,]
    #  map <- map - map[1] # shift first marker to 0

    # bc or intercross?  thetype==0 for BC and ==1 for intercross
    if(type=="f2") {
        if(class(geno) == "X") thetype <- 0
        else thetype <- 1
    }
    else thetype <- 0

    n.ind <- nrow(geno)
    n.mar <- ncol(geno)

    z <- .C("R_locate_xo",
            as.integer(n.ind),
            as.integer(n.mar),
            as.integer(thetype),
            as.integer(geno),
            as.double(map),
            location=as.double(rep(0,n.ind*2*(n.mar-1))),
            nseen=as.integer(rep(0,n.ind)),
            ileft=as.integer(rep(0,n.ind*2*(n.mar-1))),
            iright=as.integer(rep(0,n.ind*2*(n.mar-1))),
            left=as.double(rep(0,n.ind*2*(n.mar-1))),
            right=as.double(rep(0,n.ind*2*(n.mar-1))),
            gleft=as.integer(rep(0, n.ind*2*(n.mar-1))),
            gright=as.integer(rep(0, n.ind*2*(n.mar-1))),
            ntype=as.integer(rep(0,n.ind*2*(n.mar-1))),
            as.integer(full.info),
            PACKAGE="qtl")
    location <- t(matrix(z$location, nrow=n.ind))
    nseen <- z$nseen
    if(full.info) {
        ileft <- t(matrix(z$ileft, nrow=n.ind))
        iright <- t(matrix(z$iright, nrow=n.ind))
        left <- t(matrix(z$left, nrow=n.ind))
        right <- t(matrix(z$right, nrow=n.ind))
        gleft <- t(matrix(z$gleft, nrow=n.ind))
        gright <- t(matrix(z$gright, nrow=n.ind))
        ntype <- t(matrix(z$ntype, nrow=n.ind))
    }

    if(!full.info)
        res <- lapply(as.data.frame(rbind(nseen, location), stringsAsFactors=TRUE),
                      function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })
    else {
        location <- lapply(as.data.frame(rbind(nseen, location), stringsAsFactors=TRUE),
                           function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        ileft <- lapply(as.data.frame(rbind(nseen, ileft), stringsAsFactors=TRUE),
                        function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        iright <- lapply(as.data.frame(rbind(nseen, iright), stringsAsFactors=TRUE),
                         function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        left <- lapply(as.data.frame(rbind(nseen, left), stringsAsFactors=TRUE),
                       function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        right <- lapply(as.data.frame(rbind(nseen, right), stringsAsFactors=TRUE),
                        function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        gleft <- lapply(as.data.frame(rbind(nseen, gleft), stringsAsFactors=TRUE),
                        function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        gright <- lapply(as.data.frame(rbind(nseen, gright), stringsAsFactors=TRUE),
                         function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        ntype <- lapply(as.data.frame(rbind(nseen, ntype), stringsAsFactors=TRUE),
                        function(a) { if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })

        res <- location
        for(i in seq(along=res)) {
            if(length(res[[i]])>0) {
                ntype[[i]][length(ntype[[i]])] <- NA
                res[[i]] <- cbind(location=location[[i]],
                                  left=left[[i]],
                                  right=right[[i]],
                                  ileft=ileft[[i]],
                                  iright=iright[[i]],
                                  gleft=gleft[[i]],
                                  gright=gright[[i]],
                                  nTypedBetween=ntype[[i]])
            }
        }
    }

    names(res) <- id
    res
}

# jittermap: make sure no two markers are at precisely the same position
jittermap <-
    function(object, amount=1e-6) # x is either a cross object or a map
{
    if(any(class(object) == "cross")) {
        themap <- pull.map(object)
        return.cross <- TRUE
    }
    else {
        if(!any(class(object) == "map"))
            stop("Input must be a cross or a map")
        return.cross <- FALSE
        themap <- object
    }

    for(i in 1:length(themap)) {
        if(is.matrix(themap[[i]])) { # sex-specific maps
            n <- ncol(themap[[i]])
            if(n > 1) {
                for(j in 1:2)
                    themap[[i]][j,] <- themap[[i]][j,] + c(0,cumsum(rep(amount, n-1)))
            }
        }
        else {
            n <- length(themap[[i]])
            if(n > 1)
                themap[[i]] <- themap[[i]] + c(0,cumsum(rep(amount, n-1)))
        }
    }

    if(return.cross)
        return(clean(replace.map(object, themap)))

    themap
}


print.map <-
    function(x, ...)
{
    for(i in seq(along=x))
        attr(x[[i]], "loglik") <- NULL

    if(length(x) == 1)
        print(unclass(x[[1]]))
    else
        print(unclass(lapply(x, unclass)))
}


# map function for Stahl model
mf.stahl <-
    function(d, m=0, p=0)
{
    if(any(d < 0))
        stop("d must be >= 0\n")
    if(m < 0)
        stop("Must have m >= 0\n")
    if(p < 0 || p > 1)
        stop("Must have 0 <= p <= 1\n")

    # handle missing values
    if(any(!is.finite(d))) {
        which.finite <- is.finite(d)
        dsub <- d[which.finite]
        dropped.some <- TRUE
    } else {
        dsub <- d
        dropped.some <- FALSE
    }

    result <- .C("R_mf_stahl",
                 as.integer(length(dsub)),
                 as.double(dsub/100), # convert to Morgans
                 as.integer(m),
                 as.double(p),
                 result=as.double(rep(0,length(dsub))),
                 PACKAGE="qtl")$result

    if(dropped.some) {
        d[!which.finite] <- NA
        d[which.finite] <- result
        result <- d
    }
    result
}

# inverse map function for Stahl model
imf.stahl <-
    function(r, m=0, p=0, tol=1e-12, maxit=1000)
{
    if(any(r < 0 | r > 0.5))
        stop("r must be >= 0 or <= 0.5\n")
    if(m < 0)
        stop("Must have m >= 0\n")
    if(p < 0 || p > 1)
        stop("Must have 0 <= p <= 1\n")

    # handle missing values
    if(any(!is.finite(r))) {
        which.finite <- is.finite(r)
        rsub <- r[which.finite]
        dropped.some <- TRUE
    } else {
        rsub <- r
        dropped.some <- FALSE
    }

    result <- .C("R_imf_stahl",
                 as.integer(length(rsub)),
                 as.double(rsub),
                 as.integer(m),
                 as.double(p),
                 result=as.double(rep(0,length(rsub))),
                 as.double(tol),
                 as.integer(maxit),
                 PACKAGE="qtl")$result*100 # convert to cM

    if(dropped.some) {
        r[!which.finite] <- NA
        r[which.finite] <- result
        result <- r
    }
    result
}

######################################################################
# getid: internal function to pull out the "ID" column from the
#        phenotype data, if there is one
######################################################################
getid <-
    function(cross)
{
    phe <- cross$pheno
    nam <- names(phe)
    if("id" %in% nam) {
        id <- phe$id
        phenam <- "id"
    }
    else if("ID" %in% nam) {
        id <- phe$ID
        phenam <- "ID"
    }
    else if("Id" %in% nam) {
        id <- phe$Id
        phenam <- "Id"
    }
    else if("iD" %in% nam) {
        id <- phe$iD
        phenam <- "iD"
    }
    else {
        id <- NULL
        phenam <- NULL
    }

    if(is.factor(id))
        id <- as.character(id)

    attr(id, "phenam") <- phenam

    id
}

######################################################################
# find the chromosome and position of a vector of markers
######################################################################
find.markerpos <-
    function(cross, marker)
{
    if(length(marker) != length(unique(marker))) {
        warning("Dropping duplicate marker names.")
        marker <- unique(marker)
    }

    output <- data.frame(chr=rep("", length(marker)),
                         pos=rep(NA, length(marker)), stringsAsFactors=TRUE)
    output$chr <- as.character(output$chr)
    rownames(output) <- marker

    map <- pull.map(cross)
    n.mar <- nmar(cross)
    chr <- rep(names(map), n.mar)

    if(!is.matrix(map[[1]])) {
        pos <- unlist(map)
        onemap <- TRUE
    }
    else {
        pos <- unlist(lapply(map, function(a) a[1,]))
        pos2 <- unlist(lapply(map, function(a) a[2,]))
        onemap <- FALSE
        output <- cbind(output, pos2=rep(NA, length(marker)))
        colnames(output)[2:3] <- c("pos.female","pos.male")
    }

    mnam <- markernames(cross)

    for(i in seq(along=marker)) {

        wh <- match(marker[i], mnam)
        if(is.na(wh)) {
            output[i,1] <- NA
            output[i,2] <- NA
        }
        else {
            if(length(wh) > 1) {
                warning("Marker ", marker[i], " found multiple times.")
                wh <- sample(wh, 1)
            }
            output[i,1] <- chr[wh]
            output[i,2] <- pos[wh]
            if(!onemap) output[i,3] <- pos2[wh]
        }
    }

    output
}


######################################################################
# find the chromosome and position of a vector of pseudomarkers
######################################################################
find.pseudomarkerpos <-
    function(cross, marker, where=c("draws","prob"))
{
    if(length(marker) != length(unique(marker))) {
        warning("Dropping duplicate pseudomarker names.")
        marker <- unique(marker)
    }

    output <- data.frame(chr=rep("", length(marker)),
                         pos=rep(NA, length(marker)), stringsAsFactors=TRUE)
    output$chr <- as.character(output$chr)
    rownames(output) <- marker

    where <- match.arg(where)
    if(where=="draws" && !("draws" %in% names(cross$geno[[1]])))
        stop("You'll need to first run sim.geno")
    if(where=="prob" && !("prob" %in% names(cross$geno[[1]])))
        stop("You'll need to first run calc.genoprob")

    themap <- vector("list", nchr(cross))
    names(themap) <- names(cross$geno)
    for(i in 1:nchr(cross)) {
        if(where=="draws")
            themap[[i]] <- attr(cross$geno[[i]]$draws, "map")
        else
            themap[[i]] <- attr(cross$geno[[i]]$prob, "map")
    }

    if(!is.matrix(themap[[1]])) {
        pmar <- unlist(lapply(themap, names))
        pos <- unlist(themap)
        onemap <- TRUE
        chr <- rep(names(themap), sapply(themap, length))
    }
    else {
        pos <- unlist(lapply(themap, function(a) a[1,]))
        pos2 <- unlist(lapply(themap, function(a) a[2,]))
        onemap <- FALSE
        pmar <- unlist(lapply(themap, colnames))
        output <- cbind(output, pos2=rep(NA, length(marker)))
        colnames(output)[2:3] <- c("pos.female","pos.male")
        chr <- rep(names(themap), sapply(themap, ncol))
    }
    whnotmarker <- grep("^loc-*[0-9]*", pmar)
    pmar[whnotmarker] <- paste("c", chr[whnotmarker], ".", pmar[whnotmarker], sep="")

    for(i in seq(along=marker)) {
        wh <- match(marker[i], pmar)
        if(is.na(wh)) {
            output[i,1] <- NA
            output[i,2] <- NA
        }
        else {
            if(length(wh) > 1) {
                warning("Pseudomarker ", marker[i], " found multiple times.")
                wh <- sample(wh, 1)
            }
            output[i,1] <- chr[wh]
            output[i,2] <- pos[wh]
            if(!onemap) output[i,3] <- pos2[wh]
        }
    }

    output
}


######################################################################
# utility function for determining whether pheno.col (as argument
# to scanone etc) can be interpreted as a vector of phenotypes,
# versus a vector of phenotype columns
######################################################################
LikePheVector <-
    function(pheno, n.ind, n.phe)
{
    if(is.numeric(pheno) && length(pheno)==n.ind &&
       any(pheno < 1 | pheno > n.phe | pheno!=round(pheno)))
        return(TRUE)
    FALSE
}

######################################################################
# for matching chromosome names
######################################################################
matchchr <-
    function(selection, thechr)
{
    if(is.factor(thechr)) thechr <- as.character(thechr)
    if(length(thechr) > length(unique(thechr)))
        stop("Duplicate chromosome names.")

    if(is.logical(selection)) {
        if(length(selection) != length(thechr))
            stop("Logical vector to select chromosomes is the wrong length")
        return(thechr[selection])
    }

    if(is.numeric(selection)) selection <- as.character(selection)

    if(length(selection) > length(unique(selection))) {
        warning("Dropping duplicate chromosomes")
        selection <- unique(selection)
    }

    g <- grep("^-", selection)
    if(length(g) > 0 && length(g) < length(selection))
        stop("In selecting chromosomes, all must start with '-' or none should.")
    if(length(g) > 0) {
        selectomit <- TRUE
        selection <- substr(selection, 2, nchar(selection))
    }
    else selectomit <- FALSE

    wh <- match(selection, thechr)
    if(any(is.na(wh))) {
        warning("Chr ", paste(selection[is.na(wh)], collapse=", "), " not found")
        wh <- wh[!is.na(wh)]
        if(length(wh) == 0) return(thechr)
    }

    if(selectomit) return(thechr[-wh])
    thechr[sort(wh)]
}

######################################################################
# check that chromosomes match appropriately
# TRUE = chr okay
# FALSE = problem
######################################################################
testchr <-
    function(selection, thechr)
{
    if(is.factor(thechr)) thechr <- as.character(thechr)
    if(length(thechr) > length(unique(thechr))) {
        #    warning("Duplicate chromosome names.")
        return(FALSE)
    }

    if(is.logical(selection)) {
        if(length(selection) != length(thechr)) {
            #      warning("Logical vector to select chromosomes is the wrong length")
            return(FALSE)
        }
        return(TRUE)
    }

    if(is.numeric(selection)) selection <- as.character(selection)

    if(length(selection) > length(unique(selection))) {
        #    warning("Dropping duplicate chromosomes")
        selection <- unique(selection)
    }

    g <- grep("^-", selection)
    if(length(g) > 0 && length(g) < length(selection)) {
        #    stop("In selecting chromosomes, all must start with '-' or none should.")
        return(FALSE)
    }
    if(length(g) > 0) {
        selectomit <- TRUE
        selection <- substr(selection, 2, nchar(selection))
    }
    else selectomit <- FALSE

    wh <- match(selection, thechr)
    if(any(is.na(wh))) return(FALSE)

    TRUE
}

######################################################################
# convert2sa
#
# convert a sex-specific maps to a sex-averaged one.
# We pull out just the female map, and give a warning if the male and
# female maps are too different
######################################################################
convert2sa <-
    function(map, tol=1e-4)
{
    if(!("map" %in% class(map)))
        stop("Input should have class \"map\".")

    if(!is.matrix(map[[1]]))
        stop("Input map doesn't seem to be a sex-specific map.")

    theclass <- sapply(map, class)

    fem <- lapply(map, function(a) a[1,])

    dif <- sapply(map, function(a) { if(ncol(a)==1) return(diff(a))
                                     a <- apply(a, 1, diff);
                                     if(is.matrix(a)) return(max(abs(apply(a, 1, diff))))
                                     abs(diff(a)) })

    if(max(dif) > tol)
        warning("Female and male inter-marker distances differ by as much as ", max(dif), ".")

    for(i in seq(along=theclass)) class(fem[[i]]) <- theclass[i]

    class(fem) <- "map"
    fem
}

# round as character string, ensuring ending 0's are kept.
charround <-
    function(x, digits=1)
{
    if(digits < 1)
        stop("This is intended for the case digits >= 1.")

    y <- as.character(round(x, digits))

    z <- strsplit(y, "\\.")
    sapply(z, function(a, digits)
       {
           if(length(a) == 1)
               b <- paste(a[1], ".", paste(rep("0", digits),collapse=""), sep="")
           else {
               if(nchar(a[2]) == digits)
                   b <- paste(a, collapse=".")
               else
                   b <- paste(a[1], ".", a[2],
                              paste(rep("0", digits - nchar(a[2])), collapse=""),
                              sep="")
           }
       }, digits)
}

######################################################################
# scantwoperm2scanoneperm
#
# pull out the scanone permutations from scantwo permutation results,
# so that one may use the scantwo perms in calls to summary.scanone
######################################################################
scantwoperm2scanoneperm <-
    function(scantwoperms)
{
    if(!("scantwoperm" %in% class(scantwoperms)))
        stop("Input must have class \"scantwoperm\".")
    if(!("one" %in% names(scantwoperms)))
        stop("Input doesn't contain scanone permutation results.")
    scanoneperms <- scantwoperms$one
    class(scanoneperms) <- c("scanoneperm")
    scanoneperms
}


######################################################################
# subset.map
######################################################################
subset.map <-
    function(x, ...)
{
    cl <- class(x)
    x <- x[...]
    class(x) <- cl
    x
}

`[.map` <-
    function(x, ...)
{
    x <- unclass(x)[...]
    class(x) <- "map"
    x
}


######################################################################
# subset.cross with [ ]
######################################################################
`[.cross` <-
    function(x, chr, ind)
    subset(x, chr, ind)


######################################################################
# findDupMarkers
#
# find markers whose genotype data is identical to another marker
# (which thus might be dropped, as extraneous)
#
# chr          (Optional) set of chromosomes to consider
#
# exact.only   If TRUE, look for markers with the same genotypes
#                       and the same pattern of missing data
#              If FALSE, also identify markers whose observed
#                        genotypes match another marker, with no
#                        observed genotypes for which the other is
#                        missing
#
# adjacent.only   If TRUE, only consider adjacent markers
######################################################################

findDupMarkers <-
    function(cross, chr, exact.only=TRUE, adjacent.only=FALSE)
{
    if(!missing(chr)) cross <- subset(cross, chr=chr)

    g <- pull.geno(cross)
    markers <- colnames(g)
    markerloc <- lapply(nmar(cross), function(a) 1:a)
    if(length(markerloc) > 1) {
        for(j in 2:length(markerloc))
            markerloc[[j]] <- markerloc[[j]] + max(markerloc[[j-1]]) + 10
    }
    markerloc <- unlist(markerloc)

    if(exact.only) {
        g[is.na(g)] <- 0

        # genotype data patterns
        pat <- apply(g, 2, paste, collapse="")

        # table of unique values
        tab <- table(pat)

        # no duplicates; return
        if(all(tab == 1)) return(NULL)

        wh <- which(tab > 1)
        theloc <- themar <- vector("list", length(wh))
        for(i in seq(along=wh)) {
            themar[[i]] <- names(pat)[pat==names(tab)[wh[i]]]
            theloc[[i]] <- markerloc[pat==names(tab)[wh[i]]]
        }

        if(adjacent.only) {
            extraloc <- list()
            extramar <- list()
            for(i in seq(along=theloc)) {
                d <- diff(theloc[[i]]) # look for adjacency
                if(any(d>1)) { # split into adjacent groups
                    temp <- which(d>1)
                    first <- c(1, temp+1)
                    last <- c(temp, length(theloc[[i]]))
                    for(j in 2:length(first)) {
                        extraloc[[length(extraloc)+1]] <- theloc[[i]][first[j]:last[j]]
                        extramar[[length(extramar)+1]] <- themar[[i]][first[j]:last[j]]
                    }
                    themar[[i]] <- themar[[i]][first[1]:last[1]]
                    theloc[[i]] <- theloc[[i]][first[1]:last[1]]
                }
            }
            themar <- c(themar, extramar)
            theloc <- c(theloc, extraloc)

            nm <- sapply(themar, length)
            if(all(nm==1)) return(NULL) # nothing left
            themar <- themar[nm>1]
            theloc <- theloc[nm>1]
        }

        # order by first locus
        o <- order(sapply(theloc, min))
        themar <- themar[o]

        randompics <- sapply(themar, function(a) sample(length(a), 1))
        for(i in seq(along=themar)) {
            names(themar)[i] <- themar[[i]][randompics[i]]
            themar[[i]] <- themar[[i]][-randompics[i]]
        }

    }
    else {
        themar <- NULL

        ntyp <- ntyped(cross, "mar")
        o <- order(ntyp, decreasing=TRUE)

        g[is.na(g)] <- 0
        z <- .C("R_findDupMarkers_notexact",
                as.integer(nrow(g)),
                as.integer(ncol(g)),
                as.integer(g),
                as.integer(o),
                as.integer(markerloc),
                as.integer(adjacent.only),
                result=as.integer(rep(0,length(o))),
                PACKAGE="qtl")
        if(all(z$result==0)) return(NULL)
        u <- unique(z$result[z$result != 0])
        themar <- vector("list", length(u))
        names(themar) <- colnames(g)[u]
        for(i in seq(along=themar))
            themar[[i]] <- colnames(g)[z$result==u[i]]
    }

    themar
}

######################################################################
# convert2riself
######################################################################
convert2riself <-
    function(cross)
{
    if(class(cross)[2] != "cross")
        stop("input must be a cross object.")
    curtype <- class(cross)[1]
    chrtype <- sapply(cross$geno, class)
    whX <- NULL
    if(any(chrtype != "A")) { # there's an X chromosome
        whX <- names(cross$geno)[chrtype != "A"]
        if(length(whX) > 1)
            warning("Converting chromosomes ", paste(whX, collapse=" "), " to autosomal.")
        else
            warning("Converting chromosome ", whX, " to autosomal.")
        for(i in whX) class(cross$geno[[i]]) <- "A"
    }

    gtab <- table(pull.geno(cross))
    usethree <- FALSE
    if(!is.na(gtab["3"])) {
        if(!is.na(gtab["2"]) && gtab["3"] < gtab["2"])
            usethree <- FALSE
        else usethree <- TRUE
    }

    g2omit <- g3omit <- g4omit <- 0
    for(i in 1:nchr(cross)) {
        dat <- cross$geno[[i]]$data
        g1 <- sum(!is.na(dat) & dat==1)
        g2 <- sum(!is.na(dat) & dat==2)
        g3 <- sum(!is.na(dat) & dat==3)
        g4 <- sum(!is.na(dat) & dat>3)
        if(usethree && chrtype[i] == "A") {
            dat[!is.na(dat) & dat!=1 & dat!=3] <- NA
            dat[!is.na(dat) & dat==3] <- 2
            g2omit <- g2omit + g2
            g4omit <- g4omit + g4
        }
        else {
            dat[!is.na(dat) & dat!=1 & dat!=2] <- NA
            g3omit <- g3omit + g3
            g4omit <- g4omit + g4
        }
        cross$geno[[i]]$data <- dat
    }
    if(g2omit > 0)
        warning("Omitting ", g2omit, " genotypes with code==2.")
    if(g3omit > 0)
        warning("Omitting ", g3omit, " genotypes with code==3.")
    if(g4omit > 0)
        warning("Omitting ", g4omit, " genotypes with code>3.")

    class(cross)[1] <- "riself"

    cross
}


######################################################################
# convert2risib
######################################################################
convert2risib <-
    function(cross)
{
    if(class(cross)[2] != "cross")
        stop("input must be a cross object.")
    curtype <- class(cross)[1]
    chrtype <- sapply(cross$geno, class)

    gtab <- table(pull.geno(cross))
    usethree <- FALSE
    if(!is.na(gtab["3"])) {
        if(!is.na(gtab["2"]) && gtab["3"] < gtab["2"])
            usethree <- FALSE
        else usethree <- TRUE
    }

    g2omit <- g3omit <- g4omit <- 0
    for(i in 1:nchr(cross)) {
        dat <- cross$geno[[i]]$data
        g1 <- sum(!is.na(dat) & dat==1)
        g2 <- sum(!is.na(dat) & dat==2)
        g3 <- sum(!is.na(dat) & dat==3)
        g4 <- sum(!is.na(dat) & dat>3)
        if(usethree) {
            if(chrtype[i] == "A") {
                dat[!is.na(dat) & dat!=1 & dat!=3] <- NA
                dat[!is.na(dat) & dat==3] <- 2
                g2omit <- g2omit + g2
                g4omit <- g4omit + g4
            }
            else { # X chromosome
                if(g2 >= g3) {
                    dat[!is.na(dat) & dat!=1 & dat!=2] <- NA
                    g3omit <- g3omit + g3
                    g4omit <- g4omit + g4
                }
                else {
                    dat[!is.na(dat) & dat!=1 & dat!=3] <- NA
                    dat[!is.na(dat) & dat==3] <- 2
                    g2omit <- g2omit + g2
                    g4omit <- g4omit + g4
                }
            }
        }
        else {
            dat[!is.na(dat) & dat!=1 & dat!=2] <- NA
            g3omit <- g3omit + g3
            g4omit <- g4omit + g4
        }
        cross$geno[[i]]$data <- dat
    }
    if(g2omit > 0)
        warning("Omitting ", g2omit, " genotypes with code==2.")
    if(g3omit > 0)
        warning("Omitting ", g3omit, " genotypes with code==3.")
    if(g4omit > 0)
        warning("Omitting ", g4omit, " genotypes with code>3.")

    class(cross)[1] <- "risib"

    cross
}


rescalemap <-
    function(object, scale=1e-6)
{
    if("cross" %in% class(object)) {
        for(i in 1:nchr(object)) {
            object$geno[[i]]$map <-
                object$geno[[i]]$map * scale
        }
        if(abs(scale - 1) > 1e-6)
            object <- clean(object) # strip off intermediate calculations
    } else if("map" %in% class(object)) {
        for(i in seq(along=object)) {
            object[[i]] <- object[[i]] * scale
        }
    } else
        stop("rescalemap works only for objects of class \"cross\" or \"map\".")

    object
}


shiftmap <-
    function(object, offset=0)
{
    if("cross" %in% class(object)) {
        if(length(offset) == 1) offset <- rep(offset, nchr(object))
        else if(length(offset) != nchr(object))
            stop("offset must have length 1 or n.chr (", nchr(object), ")")

        for(i in 1:nchr(object)) {
            if(is.matrix(object$geno[[i]]$map)) {
                for(j in 1:2)
                    object$geno[[i]]$map[j,] <- object$geno[[i]]$map[j,] - object$geno[[i]]$map[j,1] + offset[i]
            } else {
                object$geno[[i]]$map <- object$geno[[i]]$map - object$geno[[i]]$map[1] + offset[i]
            }
        }
    } else if("map" %in% class(object)) {
        if(length(offset) == 1) offset <- rep(offset, length(object))
        else if(length(offset) != length(object))
            stop("offset must have length 1 or n.chr (", length(object), ")")
        for(i in seq(along=object)) {
            if(is.matrix(object[[i]])) {
                for(j in 1:2)
                    object[[i]][j,] <- object[[i]][j,] - object[[i]][j,1] + offset[i]
            } else {
                object[[i]] <- object[[i]] - object[[i]][1] + offset[i]
            }
        }
    } else
        stop("shiftmap works only for objects of class \"cross\" or \"map\".")

    object
}

######################################################################
# switch alleles in a cross
######################################################################
switchAlleles <-
    function(cross, markers, switch=c("AB","CD","ABCD", "parents"))
{
    type <- class(cross)[1]
    switch <- match.arg(switch)

    if(type %in% c("bc", "risib", "riself", "dh", "haploid")) {
        if(switch != "AB")
            warning("Using switch = \"AB\".")

        found <- rep(FALSE, length(markers))
        for(i in 1:nchr(cross)) {
            cn <- colnames(cross$geno[[i]]$data)
            m <- match(markers, cn)
            if(any(!is.na(m))) {
                found[!is.na(m)] <- TRUE
                for(j in m[!is.na(m)]) {
                    g <- cross$geno[[i]]$data[,j]
                    cross$geno[[i]]$data[!is.na(g) & g==1,j] <- 2
                    cross$geno[[i]]$data[!is.na(g) & g==2,j] <- 1
                }
            }
        }

    }
    else if(type=="f2") {
        if(switch != "AB")
            warning("Using switch = \"AB\".")

        found <- rep(FALSE, length(markers))
        for(i in 1:nchr(cross)) {
            cn <- colnames(cross$geno[[i]]$data)
            m <- match(markers, cn)
            if(any(!is.na(m))) {
                found[!is.na(m)] <- TRUE
                for(j in m[!is.na(m)]) {
                    g <- cross$geno[[i]]$data[,j]
                    cross$geno[[i]]$data[!is.na(g) & g==1,j] <- 3
                    cross$geno[[i]]$data[!is.na(g) & g==3,j] <- 1
                    cross$geno[[i]]$data[!is.na(g) & g==4,j] <- 5
                    cross$geno[[i]]$data[!is.na(g) & g==5,j] <- 4
                }
            }
        }
    }
    else if(type=="4way") {

        if(switch=="AB")
            newg <- c(2,1,4,3,6,5,7,8,10,9,12,11,14,13)
        else if(switch=="CD")
            newg <- c(3,4,1,2,5,6,8,7,10,9,13,14,11,12)
        else if(switch=="ABCD")
            newg <- c(4,3,2,1,6,5,8,7,9,10,14,13,12,11)
        else # switch parents
            newg <- c(1,3,2,4,7,8,5,6,9,10,11,13,12,14)

        found <- rep(FALSE, length(markers))
        for(i in 1:nchr(cross)) {
            cn <- colnames(cross$geno[[i]]$data)
            m <- match(markers, cn)
            if(any(!is.na(m))) {
                found[!is.na(m)] <- TRUE
                for(j in m[!is.na(m)]) {
                    g <- cross$geno[[i]]$data[,j]
                    for(k in 1:14)
                        cross$geno[[i]]$data[!is.na(g) & g==k,j] <- newg[k]
                }
            }
        }

    }

    if(any(!found))
        warning("Some markers not found: ", paste(markers[!found], collapse=" "))

    clean(cross)
}

######################################################################
#
# nqrank: Convert a set of quantitative values to the corresponding
#         normal quantiles (preserving mean and SD)
#
######################################################################

nqrank <-
    function(x, jitter=FALSE)
{
    y <- x[!is.na(x)]
    themean <- mean(y, na.rm=TRUE)
    thesd <- sd(y, na.rm=TRUE)

    y[y == Inf] <- max(y[y<Inf])+10
    y[y == -Inf] <- min(y[y > -Inf]) - 10
    if(jitter)
        y <- rank(y+runif(length(y))/(sd(y)*10^8))
    else y <- rank(y)

    x[!is.na(x)] <- qnorm((y-0.5)/length(y))

    x*thesd/sd(x, na.rm=TRUE)-mean(x,na.rm=TRUE)+themean
}

######################################################################
#
# cleanGeno: omit genotypes that are possibly in error, as indicated
#            by apparent double-crossovers separated by a distance of
#            no more than maxdist and having no more than maxmark
#            interior typed markers
#
######################################################################

cleanGeno <-
    function(cross, chr, maxdist=2.5, maxmark=2, verbose=TRUE)
{
    if(!(class(cross)[1] %in% c("bc", "riself", "risib", "dh", "haploid")))
        stop("This function currently only works for crosses with two genotypes")

    if(!missing(chr)) cleaned <- subset(cross, chr=chr)
    else cleaned <- cross

    thechr <- names(cleaned$geno)
    totdrop <- 0
    maxmaxdist <- max(maxdist)
    for(i in thechr) {
        xoloc <- locateXO(cleaned, chr=i, full.info=TRUE)
        nxo <- sapply(xoloc, function(a) if(is.matrix(a)) return(nrow(a)) else return(0))
        g <- pull.geno(cleaned, chr=i)

        ndrop <- 0
        for(j in which(nxo > 1)) {
            maxd <- xoloc[[j]][-1,"right"] - xoloc[[j]][-nrow(xoloc[[j]]),"left"]
            wh <- maxd <= maxmaxdist
            if(any(wh)) {
                for(k in which(wh)) {
                    nt <- sum(!is.na(g[j,(xoloc[[j]][k,"ileft"]+1):(xoloc[[j]][k+1,"iright"]-1)]))
                    if(nt > 0 && any(nt <= maxmark & maxd[k] < maxdist)) {
                        cleaned$geno[[i]]$data[j,(xoloc[[j]][k,"ileft"]+1):(xoloc[[j]][k+1,"iright"]-1)] <- NA
                        ndrop <- ndrop + nt
                        totdrop <- totdrop + nt
                    }
                }
            }
        }
        if(verbose && ndrop > 0) {
            totgen <- sum(ntyped(subset(cross, chr=i)))
            cat(" ---Dropping ", ndrop, " genotypes (out of ", totgen, ") on chr ", i, "\n", sep="")
        }
    }

    if(verbose && nchr(cleaned)>1 && totdrop > 0) {
        totgen <- sum(ntyped(subset(cross, chr=thechr)))
        cat(" ---Dropped ", totdrop, " genotypes (out of ", totgen, ") in total\n", sep="")
    }

    for(i in names(cleaned$geno))
        cross$geno[[i]] <- cleaned$geno[[i]]

    cross
}

######################################################################
# typingGap: calculate gaps between typed markers
######################################################################

typingGap <-
    function(cross, chr, terminal=FALSE)
{
    if(!missing(chr))
        cross <- subset(cross, chr)

    n.ind <- nind(cross)
    n.chr <- nchr(cross)

    gaps <- matrix(nrow=n.ind, ncol=n.chr)
    colnames(gaps) <- names(cross$geno)

    for(i in 1:n.chr) {
        map <- cross$geno[[i]]$map
        map <- c(map[1], map, map[length(map)])
        if(is.matrix(map)) stop("This function can't currently handle sex-specific maps.")

        if(terminal) # just look at terminal gaps
            gaps[,i] <- apply(cbind(1,cross$geno[[i]]$data,1), 1,
                              function(a,b) {d <- diff(b[!is.na(a)]); max(d[c(1,length(d))]) }, map)
        else
            gaps[,i] <- apply(cbind(1,cross$geno[[i]]$data,1), 1,
                              function(a,b) max(diff(b[!is.na(a)])), map)
    }
    if(n.chr==1) gaps <- as.numeric(gaps)
    gaps
}

######################################################################
# calcPermPval
#
# calculate permutation pvalues for summary.scanone()
######################################################################
calcPermPval <-
    function(peaks, perms)
{
    if(!is.matrix(peaks))
        peaks <- as.matrix(peaks)
    if(!is.matrix(perms))
        perms <- as.matrix(perms)

    ncol.peaks <- ncol(peaks)
    nrow.peaks <- nrow(peaks)
    n.perms <- nrow(perms)

    if(ncol.peaks != ncol(perms))
        stop("ncol(peaks) != ncol(perms)")

    pval <- .C("R_calcPermPval",
               as.double(peaks),
               as.integer(ncol.peaks),
               as.integer(nrow.peaks),
               as.double(perms),
               as.integer(n.perms),
               pval=as.double(rep(0,ncol.peaks*nrow.peaks)),
               PACKAGE="qtl")$pval

    matrix(pval, ncol=ncol.peaks, nrow=nrow.peaks)
}

######################################################################
# phenames: pull out phenotype names
######################################################################
phenames <-
    function(cross)
    colnames(cross$pheno)


######################################################################
# updateParallelRNG
#
# set RNGkind
# advance RNGstream by no. clusters
######################################################################
updateParallelRNG <-
    function(n.cluster=1)
{
    kind <- RNGkind()[1]
    if(kind != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")

    s <- .Random.seed
    if(n.cluster < 1) n.cluster <- 1
    for(i in 1:n.cluster) s <- nextRNGStream(s)
    .Random.seed <<- s ## global assign new .Random.seed
}

######################################################################
# formMarkerCovar
#
# cross: cross object
#
# markers: marker names or pseudomarker names (like "c5loc25.1" or "5@25.1")
#
# method: use genotype probabilities or imputated genotypes (imputed with imp or argmax)
#
# ...: passed to fill.geno, if necessary
#
######################################################################
formMarkerCovar <-
    function(cross, markers, method=c("prob", "imp", "argmax"), ...)
{
    method <- match.arg(method)

    # check if the marker names are all like "5@25.1"
    grepresult <- grep("@", markers)
    if(length(grepresult) == length(markers) && all(grepresult == seq(along=markers))) {
        spl <- strsplit(markers, "@")
        chr <- sapply(spl, "[", 1)
        pos <- as.numeric(sapply(spl, "[", 2))
        m <- match(chr, chrnames(cross))
        if(any(is.na(m)))
            stop("Some chr not found: ", paste(unique(chr[m]), collapse=" "))

        if(method=="prob")
            markers <- find.pseudomarker(cross, chr, pos, where="prob")
        else
            markers <- find.marker(cross, chr, pos)
    }

    chr <- unique(find.markerpos(cross, markers)[,1])

    cross <- subset(cross, chr=chr)
    isXchr <- (sapply(cross$geno, class) == "X")
    crosstype <- class(cross)[1]
    sexpgm <- getsex(cross)
    crossattr <- attributes(cross)

    if(method=="prob") {
        if(any(isXchr) && crosstype %in% c("f2", "bc", "bcsft")) {
            for(i in which(isXchr))
                cross$geno[[i]]$prob <- reviseXdata(crosstype, "full", sexpgm=sexpgm, prob=cross$geno[[i]]$prob,
                                                    cross.attr=crossattr)

        }

        if(any(isXchr) && any(!isXchr)) # some X, some not
            prob <- cbind(pull.genoprob(cross[!isXchr,], omit.first.prob=TRUE),
                          pull.genoprob(cross[isXchr,], omit.first.prob=TRUE))
        else # all X or all not X
            prob <- pull.genoprob(cross, omit.first.prob=TRUE)

        markercols <- sapply(strsplit(colnames(prob), ":"), "[", 1)
        m <- match(markers, markercols)
        if(any(is.na(m)))
            warning("Some markers/pseudomarkers not found: ", paste(unique(markers[is.na(m)]), collapse=" "))
        return(prob[,!is.na(match(markercols, markers)), drop=FALSE])
    }
    else {
        cross <- fill.geno(cross, method=method, ...)

        if(any(isXchr) && crosstype %in% c("f2", "bc", "bcsft")) {
            for(i in which(isXchr))
                cross$geno[[i]]$data <- reviseXdata(crosstype, "full", sexpgm=sexpgm, geno=cross$geno[[i]]$data,
                                                    cross.attr=crossattr)
        }

        geno <- pull.geno(cross)
        markercols <- colnames(geno)
        m <- match(markers, markercols)
        if(any(is.na(m)))
            warning("Some markers not found: ", paste(unique(markers[is.na(m)]), collapse=" "))

        geno <- geno[,!is.na(match(markercols, markers)), drop=FALSE]

        # expand each column
        nalle <- apply(geno, 2, function(a) length(unique(a)))
        g <- matrix(ncol=sum(nalle-1), nrow=nrow(geno))
        colnames(g) <- as.character(1:ncol(g))
        cur <- 0
        for(i in 1:ncol(geno)) {
            if(nalle[i] <= 1) next
            for(j in 2:nalle[i])
                g[,cur+j-1] <- as.numeric(geno[,i] == j)
            colnames(g)[cur - 1 + (2:nalle[i])] <- paste(colnames(geno)[i], 2:nalle[i], sep=".")
            cur <- cur + nalle[i] - 1
        }
        return(g)
    }
}

# end of util.R
