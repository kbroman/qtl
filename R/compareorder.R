######################################################################
#
# compareorder.R
#
# copyright (c) 2007-2011, Karl W Broman
# last modified May, 2011
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
# Contains: compareorder
#
######################################################################

######################################################################
# Calculate likelihood for a fixed order of markers on a given
# chromosome versus the current one
######################################################################
compareorder <-
    function(cross, chr, order, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             maxit=4000, tol=1e-6, sex.sp=TRUE)
{
    if(missing(chr)) chr <- names(cross$geno)[1]
    if(length(chr) > 1) {
        chr <- chr[1]
        warning("compareorder works on a single chromosome.")
    }

    map.function <- match.arg(map.function)
    cross <- subset(cross, chr)

    if(length(order) != totmar(cross)) {
        if(length(order) == totmar(cross)+1 || length(order) == totmar(cross)+2)
            order <- order[1:totmar(cross)]
        else
            stop("Argument 'order' should have length ", totmar(cross))
    }
    if(any(is.na(match(1:totmar(cross), order))))
        stop("order should be a permutation of the numbers 1, 2, ..., ", totmar(cross))

    orig <- est.map(cross, error.prob=error.prob, map.function=map.function,
                    maxit=maxit, tol=tol, sex.sp=sex.sp, verbose=FALSE)

    cross$geno[[1]]$data <- cross$geno[[1]]$data[,order]

    new <- est.map(cross, error.prob=error.prob, map.function=map.function,
                   maxit=maxit, tol=tol, sex.sp=sex.sp, verbose=FALSE)

    result <- matrix(0, ncol=2, nrow=2)
    dimnames(result) <- list(c("orig","new"), c("LOD", "length"))
    result[2,1] <- (attr(new[[1]], "loglik") - attr(orig[[1]], "loglik"))/log(10)

    if(is.matrix(orig[[1]])) {
        result[,2] <- c(diff(range(orig[[1]][1,])),
                        diff(range(new[[1]][1,])))
        if(sex.sp) {
            result <- cbind(result, c(diff(range(orig[[1]][2,])),
                                      diff(range(new[[1]][2,]))))
            colnames(result)[2:3] <- c("femaleLength","maleLength")
        }
    }
    else
        result[,2] <- c(diff(range(orig[[1]])),
                        diff(range(new[[1]])))

    as.data.frame(result, stringsAsFactors=TRUE)
}

# end of compareorder.R
