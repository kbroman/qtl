#####################################################################
#
# interpPositions.R
#
# copyright (c) 2011-2012, Karl W Broman
# last modified Mar, 2012
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
# Takes an "oldmap" (e.g., a physical map in bp or Mbp)
#     and a "newmap" (e.g., a genetic map in cM)
# take additional positions in the "oldmap" scale and estimate
# the corresponding positions (by interpolation or extrapolation)
# in the "newmap" scale
######################################################################

# oldmap and newmap in "map" format (list of vectors of positions)
# oldpositions as dataframe with $chr and $pos
interpPositions <-
    function(oldpositions, oldmap, newmap)
{
    orig.rownames <- rownames(oldpositions)
    # not sure why this is necessary, but it avoids a bug
    if(is.null(rownames(oldpositions)))
        rownames(oldpositions) <- paste("temprn", 1:nrow(oldpositions), sep="")
    else
        rownames(oldpositions) <- paste("temprn", rownames(oldpositions), sep="")

    oldchrnum <- match(oldpositions$chr, names(oldmap))
    newchrnum <- match(oldpositions$chr, names(newmap))

    missingchr <- is.na(oldchrnum) | is.na(newchrnum)
    if(any(missingchr))
        warning("Chromosomes ", paste(sort(unique(oldpositions$chr[missingchr])), collapse=", "), " not found")

    newpositions <- cbind(oldpositions, newpos=rep(NA, nrow(oldpositions)))
    u <- unique(oldchrnum)

    for(i in seq(along=u)) { # loop over chromosomes
        chrnam <- names(oldmap)[u[i]] # name of chromosome

        # the positions to be interpolated
        wholdpositions <- !missingchr & oldchrnum==u[i]
        theposnames <- rownames(oldpositions)[wholdpositions]
        if(!any(wholdpositions)) next

        # data frame with oldmap positions for this chromosome
        tempoldmap <- oldmap[[u[i]]]
        tempoldmap.df <- data.frame(chr=rep(chrnam, length(tempoldmap)),
                                    pos=as.numeric(tempoldmap))
        rownames(tempoldmap.df) <- names(tempoldmap)

        # add the positions to be interpolated
        tempoldmap.df <- rbind(tempoldmap.df, oldpositions[wholdpositions,,drop=FALSE])
        tempoldmap.df <- tempoldmap.df[order(tempoldmap.df$pos),,drop=FALSE]
        tempoldmap.df$chr <- as.character(tempoldmap.df$chr)

        # do the interpolation
        result <- interpmap(tempoldmap.df, newmap[chrnam])

        # paste in the interpolated positions
        newpositions[theposnames, "newpos"] <- result[theposnames, "pos"]
    }

    rownames(newpositions) <- orig.rownames
    newpositions
}

# end of interpPositions.R
