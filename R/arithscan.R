#####################################################################
#
# arithscan.R
#
# copyright (c) 2005-2015, Karl W Broman
# last modified Oct, 2015
# first written Mar, 2005
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
# Contains: +.scanone, -.scanone, +.scanoneperm, -.scanoneperm
#           +.scantwo, -.scantwo, +.scantwoperm, -.scantwoperm
#
######################################################################


"-.scanone" <-
    function(e1,e2)
{
    if(!any(class(e1) == "scanone"))
        stop("Input should have class \"scanone\".")

    if(missing(e2)) {
        class(e1) <- "data.frame"
        e1[,-(1:2)] <- -e1[,-(1:2)]
        class(e1) <- c("scanone","data.frame")

        return(e1)
    }
    if(!any(class(e2) == "scanone"))
        stop("Input should have class \"scanone\".")

    class(e1) <- class(e2) <- "data.frame"
    if(nrow(e1) != nrow(e2)) {
        u1 <- levels(e1[,1])
        u2 <- levels(e2[,2])
        u <- unique(c(u1,u2))
        if(length(u) == 0) stop("Can't subtract; no chromosomes in common.")
        e1 <- e1[!is.na(match(e1[,1],u)),]
        e2 <- e2[!is.na(match(e2[,1],u)),]
        if(nrow(e1) != nrow(e2) || any(e1[,1] != e2[,1]) || max(abs(e1[,2]-e2[,2])) < 0.01)
            stop("Can't subtract; arguments not compatible")
    }
    nc1 <- ncol(e1)
    nc2 <- ncol(e2)
    nc <- min(c(nc1,nc2))
    e1 <- e1[,1:nc]
    e1[,3:nc] <- e1[,3:nc] - e2[,3:nc]

    # zero out small stuff
    temp <- e1[,3:nc]
    temp[!is.na(temp) & abs(temp) < 1e-6] <- 0
    e1[,3:nc] <- temp

    class(e1) <- c("scanone","data.frame")

    e1
}

"+.scanone" <-
    function(e1,e2)
{
    if(!any(class(e1) == "scanone"))
        stop("Input should have class \"scanone\".")

    if(missing(e2)) return(e1)
    if(!any(class(e2) == "scanone"))
        stop("Input should have class \"scanone\".")

    class(e1) <- class(e2) <- "data.frame"
    if(nrow(e1) != nrow(e2)) {
        u1 <- levels(e1[,1])
        u2 <- levels(e2[,2])
        u <- unique(c(u1,u2))
        if(length(u) == 0) stop("Can't add; no chromosomes in common.")
        e1 <- e1[!is.na(match(e1[,1],u)),]
        e2 <- e2[!is.na(match(e2[,1],u)),]
        if(nrow(e1) != nrow(e2) || any(e1[,1] != e2[,1]) || max(abs(e1[,2]-e2[,2])) < 0.01)
            stop("Can't add; arguments not compatible")
    }
    nc1 <- ncol(e1)
    nc2 <- ncol(e2)
    nc <- min(c(nc1,nc2))
    e1 <- e1[,1:nc]
    e1[,3:nc] <- e1[,3:nc] + e2[,3:nc]
    class(e1) <- c("scanone","data.frame")

    e1
}

"-.scanoneperm" <-
    function(e1, e2)
{
    if(!any(class(e1) == "scanoneperm"))
        stop("Input should have class \"scanoneperm\".")

    if(missing(e2)) {
        e1.x <- ("xchr" %in% names(attributes(e1)))
        if(e1.x) {
            e1$A <- -e1$A
            e1$X <- -e1$X
        }
        else {
            theclass <- class(e1)
            e1 <- -unclass(e1)
            class(e1) <- theclass
        }

        return(e1)
    }
    if(!any(class(e2) == "scanoneperm"))
        stop("Input should have class \"scanoneperm\".")

    # check input
    e1.x <- ("xchr" %in% names(attributes(e1)))
    e2.x <- ("xchr" %in% names(attributes(e2)))

    if(e1.x != e2.x)
        stop("Need both or neither input to be X-chr specific.\n")

    if((e1.x && (any(dim(e1$A)!=dim(e2$A)) || any(dim(e1$X)!=dim(e2$X)))) ||
       (!e1.x && any(dim(e1) != dim(e2))))
        stop("Need input to concern the same phenotypes and no. permutations.\n")

    if(e1.x) {
        e1$A <- e1$A - e2$A
        e1$X <- e1$X - e2$X

        # zero out small stuff
        e1$A[!is.na(e1$A) & abs(e1$A) < 1e-6] <- 0
        e1$X[!is.na(e1$X) & abs(e1$X) < 1e-6] <- 0
    }
    else {
        theclass <- class(e1)
        e1 <- unclass(e1) - unclass(e2)

        # zero out small stuff
        e1[!is.na(e1) & abs(e1) < 1e-6] <- 0

        class(e1) <- theclass
    }

    e1
}

"+.scanoneperm" <-
    function(e1, e2)
{
    if(!any(class(e1) == "scanoneperm"))
        stop("Input should have class \"scanoneperm\".")
    if(missing(e2)) return(e1)
    if(!any(class(e2) == "scanoneperm"))
        stop("Input should have class \"scanoneperm\".")

    # check input
    e1.x <- ("xchr" %in% names(attributes(e1)))
    e2.x <- ("xchr" %in% names(attributes(e2)))

    if(e1.x != e2.x)
        stop("Need both or neither input to be X-chr specific.\n")

    if((e1.x && (any(dim(e1$A)!=dim(e2$A)) || any(dim(e1$X)!=dim(e2$X)))) ||
       (!e1.x && any(dim(e1) != dim(e2))))
        stop("Need input to concern the same phenotypes and no. permutations.\n")

    if(e1.x) {
        e1$A <- e1$A + e2$A
        e1$X <- e1$X + e2$X
    }
    else {
        theclass <- class(e1)
        e1 <- unclass(e1) + unclass(e2)
        class(e1) <- theclass
    }

    e1
}

######################################################################
# -.scantwo: subtract LOD scores in two scantwo results
######################################################################
"-.scantwo" <-
    function(e1, e2)
{
    if(!any(class(e1) == "scantwo"))
        stop("Input should have class \"scantwo\".")

    if(missing(e2)) {
        e1$lod <- -e1$lod
        if("scanoneX" %in% names(e1))
            e1$scanoneX <- -e1$scanoneX

        return(e1)
    }
    if(!any(class(e2) == "scantwo"))
        stop("Input should have class \"scantwo\".")

    e1x <- "scanoneX" %in% names(e1)
    e2x <- "scanoneX" %in% names(e2)

    if(any(dim(e1$map) != dim(e2$map)) ||
       length(dim(e1$lod)) != length(dim(e2$lod)) ||
       any(dim(e1$lod) != dim(e2$lod)) || e1x != e2x)
        stop("input arguments do not conform.")

    e1$lod <- e1$lod - e2$lod
    if(e1x) {
        if(!is.null(e1$scanoneX) && !is.null(e2$scanoneX))
            e1$scanoneX <- e1$scanoneX - e2$scanoneX
    }

    e1
}


######################################################################
# +.scantwo: add LOD scores in two scantwo results
######################################################################
"+.scantwo" <-
    function(e1, e2)
{
    if(!any(class(e1) == "scantwo"))
        stop("Input should have class \"scantwo\".")
    if(missing(e2)) return(e1)
    if(!any(class(e2) == "scantwo"))
        stop("Input should have class \"scantwo\".")

    e1x <- "scanoneX" %in% names(e1)
    e2x <- "scanoneX" %in% names(e2)

    if(any(dim(e1$map) != dim(e2$map)) ||
       length(dim(e1$lod)) != length(dim(e2$lod)) ||
       any(dim(e1$lod) != dim(e2$lod)) || e1x != e2x)
        stop("input arguments do not conform.")

    e1$lod <- e1$lod + e2$lod
    if(e1x) {
        if(!is.null(e1$scanoneX) && !is.null(e2$scanoneX))
            e1$scanoneX <- e1$scanoneX + e2$scanoneX
    }

    e1
}


######################################################################
# -.scantwoperm: subtract LOD scores in two scantwo permutation results
######################################################################
"-.scantwoperm" <-
    function(e1, e2)
{
    if(!any(class(e1) == "scantwoperm"))
        stop("Input should have class \"scantwoperm\".")

    if(missing(e2)) {
        # x-chr-specific
        if("AA" %in% names(e1)) {
            for(i in seq(along=e1))
                for(j in seq(along=e1[[i]]))
                    e1[[i]][[j]] <- -e1[[i]][[j]]
            return(e1)
        }

        for(i in 1:length(e1))
            e1[[i]] <- -e1[[i]]

        return(e1)
    }

    if(!any(class(e2) == "scantwoperm"))
        stop("Input should have class \"scantwoperm\".")

    # x-chr-specific
    if("AA" %in% names(e1) || "AA" %in% names(e2)) {
        if(!("AA" %in% names(e1) && "AA" %in% names(e2)))
            stop("Input must both be Xchr-specific, or neither")
        for(i in seq(along=e1)) {
            for(j in seq(along=e1[[i]])) {
                if(any(dim(e1[[i]][[j]]) != dim(e1[[i]][[j]])))
                    stop("dimensions do not match")
                e1[[i]][[j]] <- e1[[i]][[j]] - e2[[i]][[j]]
            }
        }
        return(e1)
    }

    dim1 <- sapply(e1, dim)
    dim2 <- sapply(e2, dim)
    if(any(dim1 != dim2))
        stop("Need input to concern the same phenotypes and no. permutations.\n")

    for(i in 1:length(e1))
        e1[[i]] <- e1[[i]] - e2[[i]]

    e1
}


######################################################################
# +.scantwoperm: add LOD scores in two scantwo permutation results
######################################################################
"+.scantwoperm" <-
    function(e1, e2)
{
    if(!any(class(e1) == "scantwoperm"))
        stop("Input should have class \"scantwoperm\".")
    if(missing(e2)) return(e1)
    if(!any(class(e2) == "scantwoperm"))
        stop("Input should have class \"scantwoperm\".")

    # x-chr-specific
    if("AA" %in% names(e1) || "AA" %in% names(e2)) {
        if(!("AA" %in% names(e1) && "AA" %in% names(e2)))
            stop("Input must both be Xchr-specific, or neither")
        for(i in seq(along=e1)) {
            for(j in seq(along=e1[[i]])) {
                if(any(dim(e1[[i]][[j]]) != dim(e1[[i]][[j]])))
                    stop("dimensions do not match")
                e1[[i]][[j]] <- e1[[i]][[j]] + e2[[i]][[j]]
            }
        }
        return(e1)
    }

    dim1 <- sapply(e1, dim)
    dim2 <- sapply(e2, dim)
    if(any(dim1 != dim2))
        stop("Need input to concern the same phenotypes and no. permutations.\n")

    for(i in 1:length(e1))
        e1[[i]] <- e1[[i]] + e2[[i]]

    e1
}


# end of arithscan.R
