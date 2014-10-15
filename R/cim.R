######################################################################
#
# cim.R
#
# copyright (c) 2009, Karl W Broman
# last modified Apr, 2009
# first written Jan, 2007
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
# Contains: cim, markerforwsel, markerforwself2, expandf2covar
#
######################################################################

######################################################################
# CIM by first doing forward selection at the markers (with filled-in
# data) to a fixed number of markers, followed by interval mapping with
# the selected markers as covariates, dropping marker covariates if
# they are within some fixed window of the location under test.
######################################################################
cim <-
    function(cross, pheno.col=1, n.marcovar=3, window=10,
             method=c("em", "imp", "hk", "ehk"),
             imp.method=c("imp", "argmax"), error.prob=0.0001,
             map.function=c("haldane", "kosambi", "c-v", "morgan"),
             n.perm)
{
    method <- match.arg(method)
    imp.method <- match.arg(imp.method)
    map.function <- match.arg(map.function)

    type <- class(cross)[1]

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(any(is.na(num))) {
            if(sum(is.na(num)) > 1)
                stop("Couldn't identify phenotypes ", paste(paste("\"", pheno.col[is.na(num)], "\"", sep=""),
                                                            collapse=" "))
            else
                stop("Couldn't identify phenotype \"", pheno.col[is.na(num)], "\"")
        }
        pheno.col <- num
    }

    if(any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")

    y <- cross$pheno[,pheno.col]
    if(any(is.na(y))) {
        cross <- subset(cross, ind=(!is.na(y)))
        y <- y[!is.na(y)]
    }

    if(!missing(n.perm) && n.perm > 0) {
        results <- matrix(ncol=1, nrow=n.perm)
        for(i in 1:n.perm) {
            o <- sample(length(y))
            y <- y[o]
            cross$pheno[,pheno.col] <- y
            temp <- cim(cross, pheno.col=pheno.col, n.marcovar=n.marcovar,
                        window=window, method=method, imp.method=imp.method,
                        error.prob=error.prob, map.function=map.function)
            results[i,1] <- max(temp[,3], na.rm=TRUE)
        }
        class(results) <- c("scanoneperm", "matrix")
        return(results)
    }

    window <- window/2 # window specifies twice the distance between marker and test position

    g <- pull.geno(cross)

    if(any(is.na(g)))
        g <- pull.geno(fill.geno(cross, method=imp.method, error.prob=error.prob,
                                 map.function=map.function))

    if(type=="f2") out.forw <- markerforwself2(g, y, n.marcovar)
    else out.forw <- markerforwsel(g, y, n.marcovar)

    mar <- colnames(g)[out.forw]
    chrpos <- find.markerpos(cross, mar)

    ac <- g[,mar,drop=FALSE]

    if(type=="f2") useac <- expandf2covar(ac)
    else useac <- ac

    firstscan <- scanone(cross, pheno.col=pheno.col, addcovar=useac,
                         method=method)

    # scan again, dropping one marker covariate at a time
    for(i in seq(along=mar)) {
        if(type=="f2") useac <- expandf2covar(ac[,-i])
        else useac <- ac[,-i]

        temp <- scanone(cross, pheno.col=pheno.col, addcovar=useac,
                        method=method, chr=chrpos[i,1])
        wh1 <- (firstscan[,1]==chrpos[i,1] & firstscan[,2] >= chrpos[i,2]-window &
                firstscan[,2] <= chrpos[i,2]+window)
        wh2 <- (temp[,2] >= chrpos[i,2]-window & temp[,2] <= chrpos[i,2] + window)

        firstscan[wh1,3] <- temp[wh2,3]
    }

    attr(firstscan, "marker.covar") <- mar
    attr(firstscan, "marker.covar.pos") <- chrpos

    u <- table(chrpos[,1])
    if(any(u>1)) { # chromosomes with multiple marker covariates
        u <- names(u)[u>1]

        for(j in u) {
            wh <- which(chrpos[,1]==j)
            pos <- chrpos[wh,2] # positions of the covariates on this chromosome

            scanpos <- firstscan[firstscan[,1]==j,2] # positions at which the genome scan is performed

            # matrix indicating, for each position, whether marker covariates need to be dropped
            need2drop <- t(sapply(scanpos, function(a,b,d) as.numeric(abs(a-b) <= d), pos, window))
            n2drop <- apply(need2drop, 1, sum)
            if(any(n2drop > 1)) {
                pat2drop <- apply(need2drop, 1, paste, collapse="")
                thepat <- unique(pat2drop[n2drop > 1])

                for(k in thepat) {
                    whpos <- which(pat2drop==k)
                    whpos2 <- (firstscan[,1]==j & !is.na(match(firstscan[,2], scanpos[whpos])))

                    tempac <- ac[,-wh[as.logical(need2drop[whpos[1],])]]

                    if(type=="f2") useac <- expandf2covar(tempac)
                    else useac <- tempac
                    temp <- scanone(cross, pheno.col=pheno.col, addcovar=useac,
                                    method=method, chr=j)

                    firstscan[whpos2,3] <- temp[whpos,3]
                }
            }

        }
    }
    firstscan
}

######################################################################
# Simple forward selection to a fixed number of covariates
#
# x = matrix of covariates
# y = outcome
# maxsize = maximum size of model
#
# output: indices of chosen covariates [1, 2, ..., ncol(x)]
######################################################################
markerforwsel <-
    function(x, y, maxsize=7)
{
    if(length(y) != nrow(x))
        stop("Need length(y) == nrow(x).")

    if(maxsize < 0 || maxsize > ncol(x))
        stop("Need maxsize between 1 and ncol(x).")

    out <- .C("R_markerforwsel",
              as.integer(nrow(x)),
              as.integer(ncol(x)),
              as.double(x),
              as.double(y),
              as.integer(maxsize),
              chosen=as.integer(rep(0,maxsize)),
              rss=as.double(rep(0,maxsize)),
              PACKAGE="qtl")

    #  out$chosen+1
    temp <- out$chosen+1
    attr(temp, "rss") <- out$rss
    temp
}

######################################################################
# The same as markerforwsel, but for an intercross, in which
# we need to expand each marker to two columns and do selection
# with those pairs of columns
######################################################################
markerforwself2 <-
    function(x, y, maxsize=7)
{
    if(length(y) != nrow(x))
        stop("Need length(y) == nrow(x).")

    if(maxsize < 0 || maxsize > ncol(x))
        stop("Need maxsize between 1 and ncol(x).")

    out <- .C("R_markerforwself2",
              as.integer(nrow(x)),
              as.integer(ncol(x)),
              as.integer(x),
              as.double(y),
              as.integer(maxsize),
              chosen=as.integer(rep(0,maxsize)),
              rss=as.double(rep(0,maxsize)),
              PACKAGE="qtl")

    #  out$chosen+1
    temp <- out$chosen+1
    attr(temp, "rss") <- out$rss
    temp
}

######################################################################
# expand covariates for F2
######################################################################
expandf2covar <-
    function(thecovar)
{
    if(is.null(thecovar) || (is.matrix(thecovar) && ncol(thecovar)==0)) return(NULL)
    if(!is.matrix(thecovar))
        return(cbind((thecovar==3) - (thecovar==1),
                     as.numeric(thecovar==2)))

    revcovar <- matrix(ncol=ncol(thecovar)*2, nrow=nrow(thecovar))
    for(i in 1:ncol(thecovar)) {
        revcovar[,i*2-1] <- (thecovar[,i]==3) - (thecovar[,i]==1)
        revcovar[,i*2] <- as.numeric(thecovar[,i]==2)
    }

    revcovar
}

# end of cim.R
