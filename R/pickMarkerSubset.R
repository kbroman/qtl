#####################################################################
#
# pickMarkerSubset.R
#
# copyright (c) 2011-2013, Karl W Broman
# last modified Apr, 2013
# first written Nov, 2011
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
######################################################################

pickMarkerSubset <-
    function(locations, min.distance, weights)
{
    n.loc <- length(locations)
    if(n.loc==1) return(names(locations)) # just one marker

    if(missing(weights))
        weights <- rep(1, n.loc)
    else {
        if(n.loc != length(weights))
            stop("length(locations) != length(weights) [", n.loc, " != ", length(weights), "]")
    }
    if(is.null(names(locations)))
        names(locations) <- 1:n.loc

    if(any(diff(locations) < 0)) {
        o <- order(locations)
        weights <- weights[o]
        locations <- locations[o]
        warning("Markers are not in order; sorting them.")
    }

    z <- .C("R_pickMarkerSubset",
            as.double(locations),
            as.integer(n.loc),
            as.double(weights),
            as.double(min.distance),
            path=as.integer(rep(0, n.loc)),
            n.path=as.integer(0),
            PACKAGE="qtl")

    path <- rev(z$path[1:z$n.path]+1) # reverse and add 1

    return(names(locations)[path])
}

# end of pickMarkerSubset.R
