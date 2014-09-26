#####################################################################
#
# addmarker.R
#
# copyright (c) 2013, Karl W Broman
# last modified Sep, 2013
# first written May, 2013
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
# Contains: addmarker.R
######################################################################

######################################################################
# add a marker to a cross object
#
# cross: the cross object
# genotypes: vector of numeric genotypes
# markername: character string with the marker name
# chr: character string with the chromosome ID
# pos: position of the marker
######################################################################
addmarker <-
    function(cross, genotypes, markername, chr, pos)
{
    if(!("cross" %in% class(cross)))
        stop("cross must have class \"cross\"")

    chr <- matchchr(chr, names(cross$geno))
    if(length(chr) != 1)
        stop("Only 1 chromosome should be specified.")
    if(length(pos) != 1)
        stop("length(pos) != 1: length(pos) = ", length(pos))
    if(length(markername) != 1)
        stop("length(markername) != 1: length(markername) = ", length(markername))
    if(length(genotypes) != nind(cross))
        stop("length(genotypes) != nind(cross): length(genotypes) = ", length(genotypes),
             ", nind(cross) = ", nind(cross))

    map <- cross$geno[[chr]]$map
    where <- sum(map < pos)
    g <- cross$geno[[chr]]$data
    if(where == 0) {
        map <- c("tmp"=pos, map)
        names(map)[1] <- markername
        cross$geno[[chr]]$map <- map
        g <- cbind("tmp"=genotypes, g)
        colnames(g)[1] <- markername
        cross$geno[[chr]]$data <- g
    }
    else if(where == length(map)) {
        map <- c(map,"tmp"=pos)
        names(map)[length(map)] <- markername
        cross$geno[[chr]]$map <- map
        g <- cbind(g, "tmp"=genotypes)
        colnames(g)[length(map)] <- markername
        cross$geno[[chr]]$data <- g
    }
    else {
        mnames <- names(map)
        map <- c(map[1:where], "tmp"=pos, map[-(1:where)])
        names(map) <- c(mnames[1:where], markername, mnames[-(1:where)])
        cross$geno[[chr]]$map <- map
        g <- cbind(g[,1:where, drop=FALSE], genotypes, g[,-(1:where), drop=FALSE])
        colnames(g) <- names(map)
        cross$geno[[chr]]$data <- g
    }
    cross
}

# end of addmarker.R
