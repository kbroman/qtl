#####################################################################
#
# replacemap.R
#
# copyright (c) 2001-2011, Karl W Broman
# last modified May, 2011
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
# Contains: replace.map, replacemap, replacemap.cross,
#           replacemap.scanone, replacemap.scantwo
######################################################################

######################################################################
#
# replace.map
#
# replace the map portion of a cross object with a list defining a map
#
######################################################################

replace.map <- function(cross, map) replacemap.cross(cross, map)

replacemap.cross <-
    function(object, map)
{
    cross <- object
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    chr.names <- names(cross$geno)
    chr.names2 <- names(map)

    m <- match(chr.names2, chr.names)
    if(any(is.na(m))) { # some extra chr in map
        extra <- chr.names2[is.na(m)]
        map <- map[!is.na(m)]
        chr.names2 <- names(map)
        warning("Extra chr in map: ", paste(extra, collapse=" "))
    }

    cross.sexsp <- sapply(cross$geno, function(a) is.matrix(a$map))
    map.sexsp <- sapply(map, is.matrix)
    if(all(cross.sexsp)) cross.sexsp <- TRUE
    else if(all(!cross.sexsp)) cross.sexsp <- FALSE
    else stop("In cross, some maps sex-specific, some not.")

    if(all(map.sexsp)) map.sexsp <- TRUE
    else if(all(!map.sexsp)) map.sexsp <- FALSE
    else stop("In map, some chrsex-specific, some not.")

    m <- match(chr.names, chr.names2)
    if(any(is.na(m))) {
        for(i in chr.names2) {
            if(cross.sexsp) mnames <- colnames(cross$geno[[i]]$map)
            else mnames <- names(cross$geno[[i]]$map)

            if(map.sexsp) mnames2 <- colnames(map[[i]])
            else mnames2 <- names(map[[i]])

            if(length(mnames) != length(mnames2))
                stop("Different numbers of markers on chr ", i)
            if(!all(mnames==mnames2))
                stop("Different marker names on chr ", i)

            cross$geno[[i]]$map <- map[[i]]
        }
    }
    else { # the same chromosomes
        n.chr <- nchr(cross)
        n.mar <- nmar(cross)

        n.chr2 <- length(map)

        if(cross.sexsp) mnames <- unlist(lapply(cross$geno, function(a) colnames(a$map)))
        else mnames <- unlist(lapply(cross$geno, function(a) names(a$map)))

        if(map.sexsp) {
            mnames2 <- unlist(lapply(map, colnames))
            n.mar2 <- sapply(map, ncol)
        }
        else {
            mnames2 <- unlist(lapply(map, names))
            n.mar2 <- sapply(map, length)
        }

        # check that things line up properly
        if(n.chr != n.chr2)
            stop("Numbers of chromosomes don't match.")
        if(any(names(cross$geno) != names(map)))
            stop("Chromosome names don't match.")
        if(any(n.mar != n.mar2))
            stop("Number of markers don't match.")
        if(any(mnames != mnames2))
            stop("Marker names don't match.")

        # proceed if no errors
        for(i in 1:length(cross$geno))
            cross$geno[[i]]$map <- map[[i]]
    }

    #### < BUG > ####
    # the next two things don't work for sex-sp maps (4-way cross)
    #################

    # maps in geno prob
    if("prob" %in% names(cross$geno[[1]])) {
        for(i in names(cross$geno)) {
            if("map" %in% names(attributes(cross$geno[[i]]$prob))) {
                temp <- attr(cross$geno[[i]]$prob, "map")
                tempr <- interpmap(data.frame(chr=rep(i, length(temp)), pos=temp, stringsAsFactors=TRUE), map)[,2]
                names(tempr) <- names(temp)
                attr(cross$geno[[i]]$prob, "map") <- tempr
            }
        }
    }

    # maps in draws
    if("draws" %in% names(cross$geno[[1]])) {
        for(i in names(cross$geno)) {
            if("map" %in% names(attributes(cross$geno[[i]]$draws))) {
                temp <- attr(cross$geno[[i]]$draws, "map")
                tempr <- interpmap(data.frame(chr=rep(i, length(temp)), pos=temp, stringsAsFactors=TRUE), map)[,2]
                names(tempr) <- names(temp)
                attr(cross$geno[[i]]$draws, "map") <- tempr
            }
        }
    }

    # maps in argmax
    if("argmax" %in% names(cross$geno[[1]])) {
        for(i in names(cross$geno)) {
            if("map" %in% names(attributes(cross$geno[[i]]$argmax))) {
                temp <- attr(cross$geno[[i]]$argmax, "map")
                tempr <- interpmap(data.frame(chr=rep(i, length(temp)), pos=temp, stringsAsFactors=TRUE), map)[,2]
                names(tempr) <- names(temp)
                attr(cross$geno[[i]]$argmax, "map") <- tempr
            }
        }
    }

    cross
}

# generic function
replacemap <- function(object, map) UseMethod("replacemap")

replacemap.scanone <-
    function(object, map)
{
    object[,2] <- interpmap(object[,1:2], map)[,2]
    object
}


replacemap.scantwo <-
    function(object, map)
{
    object$map[,2] <- interpmap4scantwo(object, map)[,2]
    object
}


######################################################################
# interpolate map positions from one to another
######################################################################
interpmap <-
    function(oldmap, newmap)
{
    returnasmap <- TRUE
    if(is.data.frame(oldmap)) {
        origmap <- oldmap
        pos <- oldmap[,2]
        names(pos) <- rownames(oldmap)
        oldmap <- split(pos, oldmap[,1])
        returnasmap <- FALSE
    }

    ochrnam <- names(oldmap)
    nchrnam <- names(newmap)
    m <- match(ochrnam, nchrnam)
    if(any(is.na(m))) {
        u <- ochrnam[is.na(m)]
        stop("Chr ", paste(u, collapse=" "), " not found in new map")
    }
    for(i in seq(along=ochrnam)) {
        omap <- oldmap[[i]]
        nmap <- newmap[[m[i]]]

        mm <- match(names(omap), names(nmap))
        wh <- is.na(mm)

        if(sum(!wh) > 1 && any(diff(mm[!wh])<0))
            stop("Need old and new maps to have markers in the same order.")

        if(sum(wh)==0) {
            oldmap[[i]] <- nmap
            next
        }

        if(sum(!wh) < 2)
            stop("Need at least two markers per chromosome in both old and new maps")

        nnmap <- omap
        nnmap[!wh] <- nmap[mm[!wh]]

        nL <- diff(range(nmap[mm[!wh]]))
        oL <- diff(range(omap[!wh]))

        onmap <- which(!wh)
        first <- onmap[1]
        last <- max(onmap)
        notonmap <- which(wh)

        for(j in notonmap) {
            if(!any(onmap < j))  # before first marker
                nnmap[j] <- nnmap[first] - (omap[first] - omap[j])*nL/oL
            else if(!any(onmap > j))  # after last marker
                nnmap[j] <- nnmap[last] + (omap[j] - omap[last])*nL/oL
            else {
                left <- max(onmap[onmap < j])
                right <- min(onmap[onmap > j])
                nnmap[j] <- nnmap[left] + (omap[j]-omap[left])*(nnmap[right]-nnmap[left])/
                    (omap[right]-omap[left])
            }
        }
        oldmap[[i]] <- nnmap
    }

    if(!returnasmap) {
        origmap[,2] <- unlist(oldmap)
        return(origmap)
    }
    oldmap
}

# like interpmap, but special for scantwo to deal with
# the case that scantwo was run without incl.markers=TRUE,
# so we need to add all marker positions to the map thing
interpmap4scantwo <-
    function(output, newmap)
{
    themap <- output$map[,1:2]

    pos <- themap[,2]
    names(pos) <- rownames(themap)
    themapalt <- split(pos, themap[,1])

    markermap <- attr(output, "fullmap")

    omapnam <- names(themapalt)
    nmapnam <- names(markermap)

    m <- match(omapnam, nmapnam)
    if(any(is.na(m))) {
        u <- omapnam[is.na(m)]
        stop("Chr ", paste(u, collapse=" "), " not found in new map")
    }

    flag <- FALSE
    for(i in seq(along=m)) {
        omap <- themapalt[[i]]
        nmap <- markermap[[m[i]]]

        mm <- match(names(nmap), names(omap))
        if(any(is.na(mm))) {
            flag <- TRUE
            themapalt[[i]] <- sort(c(omap, nmap[is.na(mm)]))
        }
    }

    if(flag) {
        revmap <- data.frame(chr=factor(rep(names(themapalt), sapply(themapalt, length)),
                             levels=names(themapalt)),
                             pos=unlist(themapalt), stringsAsFactors=TRUE)
        rownames(revmap) <- unlist(lapply(themapalt, names))
        revmap[,2] <- interpmap(revmap, newmap)[,2]
        themap <- revmap[rownames(themap),]
    }
    else
        themap[,2] <- interpmap(themap, newmap)[,2]

    themap
}

# end of replacemap.R
