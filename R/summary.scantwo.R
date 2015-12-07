######################################################################
#
# summary.scantwo.R
#
# copyright (c) 2001-2015, Karl W Broman, Hao Wu, and Brian Yandell
# last modified Oct, 2015
# first written Nov, 2001
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
# Contains: summary.scantwo, print.summary.scantwo,
#           max.scantwo, clean.scantwo, print.scantwo, subset.scantwo
#           summary.scantwoperm, print.summary.scantwoperm
#           condense.scantwo, summary.scantwocondensed
#           max.scantwocondensed, print.summary.addpair
#           rbind.scantwoperm, c.scantwoperm, subset.scantwoperm
#           [.scantwoperm
#
######################################################################
#
# summarize the result from scantwo
#
######################################################################
summary.scantwo <-
    function(object, thresholds,
             what=c("best", "full", "add", "int"),
             perms, alphas, lodcolumn=1, pvalues=FALSE,
             allpairs=TRUE, ...)
{
    if(!any(class(object) == "scantwo") &&
       !any(class(object) == "scantwocondensed"))
        stop("Input should have class \"scantwo\".")

    addpair <- attr(object, "addpair")
    if(!is.null(addpair) && addpair) { # results from addpair() that need special treatment
        if("lod.minus1" %in% names(attributes(object))) { # asymmetric formula
            attr(object, "addpair") <- NULL
            x <- summary.scantwo(object, allpairs=allpairs)
            class(x) <- "data.frame"
            x <- x[,c(1,2,8,9,10,3,4,5),drop=FALSE]

            mlod.minus1 <- tapply(attr(object, "lod.minus1"), object$map$chr, max, na.rm=TRUE)
            mlod.minus2 <- tapply(attr(object, "lod.minus2"), object$map$chr, max, na.rm=TRUE)

            w <- which(x[,1] == x[,2])
            for(i in seq(along=w))
                if(x[w[i],5] < x[w[i],8]) x[w[i],3:5] <- x[w[i],c(7,6,8)]

            if(any(x[,1] != x[,2])) {
                w <- which(x[,1] != x[,2])
                y <- x[w,,drop=FALSE]
                for(i in 1:nrow(y))
                    y[,1:5] <- y[,c(2,1,7,6,8)]

                neworder <- cbind(1:nrow(x), rep(NA, nrow(x)))
                neworder[w,2] <- nrow(x) + 1:nrow(y)
                neworder <- as.numeric(t(neworder))

                x <- rbind(x, y)[neworder[!is.na(neworder)],,drop=FALSE]
            }

            x <- cbind(x[,1:5], lod.2v1b=rep(NA,nrow(x)), lod.2v1a=rep(NA,nrow(x)))
            names(x)[5] <- "lod.2v0"

            x[,6] <- x[,5] - mlod.minus1[as.character(x[,2])]
            x[,7] <- x[,5] - mlod.minus2[as.character(x[,1])]

            if(!missing(thresholds)) {
                if(length(thresholds) > 2)
                    warning("Only the first two values in thresholds are used.")

                if(length(thresholds) == 1) thresholds <- c(thresholds, 0)
                x <- x[!is.na(x[,5]) & x[,5] >= thresholds[1] &
                       ((!is.na(x[,6]) & x[,6] >= thresholds[2]) |
                        (!is.na(x[,7]) & x[,7] >= thresholds[2])),,drop=FALSE]
            }

            class(x) <- c("summary.addpair", "data.frame")
            return(x)
        }
        else { # symmetric formula
            attr(object, "addpair") <- NULL
            if(missing(thresholds))
                x <- summary.scantwo(object, allpairs=allpairs)
            else
                x <- summary.scantwo(object, thresholds=thresholds, allpairs=allpairs)

            x <- x[,1:6]
            colnames(x)[5:6] <- c("lod.2v0", "lod.2v1")

            if(!missing(thresholds)) {
                if(length(thresholds) > 2)
                    warning("Only the first two values in thresholds are used.")

                if(length(thresholds) == 1) thresholds <- c(thresholds, 0)
                x <- x[!is.na(x[,5]) & x[,5] >= thresholds[1] &
                       !is.na(x[,6]) & x[,6] >= thresholds[2],,drop=FALSE]
            }

            class(x) <- c("summary.addpair", "data.frame")
            return(x)
        }
    }

    what <- match.arg(what)

    if(!missing(thresholds)) {
        if(length(thresholds) != 5)
            stop("If thresholds are given, there must be 5 of them.")
        else if(!is.numeric(thresholds))
            stop("thresholds should be a numeric vector")
    }

    if(!missing(alphas)) {
        if(length(alphas) != 5) {
            if(length(alphas)==1)
                alphas <- rep(alphas, 5)
            else
                stop("If alphas are given, there must be 5 of them.")
        }
        else if(!is.numeric(alphas))
            stop("alphas should be a numeric vector")
    }

    if(!missing(perms) && !any(class(perms) == "scantwoperm"))
        stop("perms must be in scantwoperm format.")

    # subset object and permutations, if necessary
    if(any(class(object) == "scantwo")) {
        d <- dim(object$lod)
        if(length(d)==3) {
            if(!missing(perms)) {
                ncp <- sapply(perms, ncol)
                if(all(ncp==1)) onepermcol <- TRUE
                else onepermcol <- FALSE
                if(any(ncp != d[3])) {
                    if(onepermcol)  {
                        if(lodcolumn > 1)
                            warning("Just one column of permutation results; assuming they apply to all LOD score columns.")
                    }
                    else
                        stop("perms have different numbers of columns as object input.\n")
                }
            }

            if(lodcolumn < 1 || lodcolumn > d[3])
                stop("lodcolumn must be between 1 and ", d[3])

            object$lod <- object$lod[,,lodcolumn]
            if(!missing(perms) && !onepermcol)
                perms <- lapply(perms, function(a, b) a[,b,drop=FALSE], lodcolumn)
        }
    }
    else { # condensed version
        if(is.matrix(object$pos1.jnt)) {
            d <- ncol(object$pos1.jnt)
            if(!missing(perms)) {
                ncp <- sapply(perms, ncol)
                if(all(ncp==1)) onepermcol <- TRUE
                else onepermcol <- FALSE
                if(any(ncp != d[3])) {
                    if(onepermcol)
                        warning("Just one column of permutation results; reusing for all LOD score columns.")
                    else
                        stop("perms have different numbers of columns as object input.\n")
                }
            }

            if(lodcolumn < 1 || lodcolumn > d)
                stop("lodcolumn must be between 1 and ", d)

            for(i in 3:length(object))
                object[[i]] <- object[[i]][,lodcolumn]

            if(!missing(perms) && !onepermcol)
                perms <- lapply(perms, function(a, b) a[,b,drop=FALSE], lodcolumn)
        }
    }

    # check input
    if(missing(perms) && !missing(alphas))
        stop("If alphas are to be used, permutation results must be provided.")
    if(!missing(thresholds) && !missing(alphas))
        stop("Only one of threshold and alpha should be specified.")
    if(pvalues && what != "best") {
        pvalues <- FALSE
        warning("pvalues shown only with what=\"best\".")
    }
    if(pvalues && missing(perms)) {
        pvalues <- FALSE
        warning("p-values may be calculated only if perms are provided.")
    }

    if(any(class(object) == "scantwo"))
        out <- subrousummaryscantwo(object, for.perm=FALSE)
    else
        out <- as.data.frame(object)

    if(!allpairs) # only look at self-self cases
        out <- out[out$chr1==out$chr2,]

    if(what=="best") {
        p1.f <- out$pos1.jnt
        p2.f <- out$pos2.jnt
        p1.a <- out$pos1.add
        p2.a <- out$pos2.add

        lf <- out$jnt.lod.full
        li <- out$jnt.lod.full - out$add.lod.add
        lfv1 <- out$jnt.lod.full - out$lod.1qtl
        la <- out$add.lod.add
        lav1 <- out$add.lod.add - out$lod.1qtl
    }
    else if(what=="full") {
        p1.f <- p1.a <- out$pos1.jnt
        p2.f <- p2.a <- out$pos2.jnt

        lf <- out$jnt.lod.full
        li <- out$jnt.lod.full - out$jnt.lod.add
        lfv1 <- out$jnt.lod.full - out$lod.1qtl
        la <- out$jnt.lod.add
        lav1 <- out$jnt.lod.add - out$lod.1qtl
    }
    else if(what=="add") {
        p1.f <- p1.a <- out$pos1.add
        p2.f <- p2.a <- out$pos2.add

        lf <- out$add.lod.full
        li <- out$add.lod.full - out$add.lod.add
        lfv1 <- out$add.lod.full - out$lod.1qtl
        la <- out$add.lod.add
        lav1 <- out$add.lod.add - out$lod.1qtl
    }
    else { # what == "int"
        p1.f <- p1.a <- out$pos1.int
        p2.f <- p2.a <- out$pos2.int

        lf <- out$int.lod.full
        li <- out$int.lod.full - out$int.lod.add
        lfv1 <- out$int.lod.full - out$lod.1qtl
        la <- out$int.lod.add
        lav1 <- out$int.lod.add - out$lod.1qtl
    }

    out <- data.frame(chr1=out$chr1, chr2=out$chr2,
                      pos1f=p1.f,
                      pos2f=p2.f,
                      lod.full=lf, lod.fv1=lfv1, lod.int=li,
                      pos1a=p1.a,
                      pos2a=p2.a,
                      lod.add=la, lod.av1=lav1)

    if(what != "best") {
        out <- out[,-(8:9)]
        names(out)[3:4] <- c("pos1", "pos2")
    }

    if(!missing(perms) && "AA" %in% names(perms)) {
        xchr_specific <- TRUE
        chrtype <- attr(perms, "chrtype")
        chrpair_type <- paste0(chrtype[out$chr1], chrtype[out$chr2])
        if(all(chrpair_type=="AA")) {
            perms <- perms$AA
            xchr_specific <- FALSE
        }
        else if(all(chrpair_type=="XX")) {
            perms <- perms$XX
            xchr_specific <- FALSE
        }
        else if(all(chrpair_type=="AX")) {
            perms <- perms$AX
            xchr_specific <- FALSE
        }
    }
    else xchr_specific <- FALSE

    if(!missing(alphas)) { # get thresholds
        if(!xchr_specific) {
            thresholds <- rep(0,5)
            for(i in 1:5)
                thresholds[i] <- quantile(perms[[i]], 1-alphas[i])
            thresholds[alphas==1] <- 0
            thresholds[alphas==0] <- Inf
        }
        else { # x-chr-specific
            thresholds <- matrix(nrow=3, ncol=5)
            for(j in 1:3) {
                for(i in 1:5)
                    thresholds[j,i] <- quantile(perms[[j]][[i]], 1-alphas[i])
            }
            thresholds[alphas==1] <- 0
            thresholds[alphas==0] <- Inf
        }
    }

    if(!missing(thresholds) || !missing(alphas)) { # apply thresholds
        if(!xchr_specific) {
            out <- out[(out$lod.full >= thresholds[1] & (out$lod.fv1 >= thresholds[2] |
                        out$lod.int >= thresholds[3])) |
                       (out$lod.add >= thresholds[4] & out$lod.av1 >= thresholds[5]),,drop=FALSE]
        }
        else {
            tokeep <- NULL
            for(i in 1:nrow(out)) {
                this_chrpair_type <- match(chrpair_type[i], c("AA", "AX", "XX"))
                if((out$lod.full[i] >= thresholds[this_chrpair_type, 1] &
                    (out$lod.fv1[i] >= thresholds[this_chrpair_type, 2] |
                     out$lod.int[i] >= thresholds[this_chrpair_type, 3])) |
                   (out$lod.add[i] >= thresholds[this_chrpair_type, 4] &
                    out$lod.av1[i] >= thresholds[this_chrpair_type, 5]))
                    tokeep <- c(tokeep, i)
            }
            out <- out[tokeep, , drop=FALSE]
        }
    }

    if(pvalues && nrow(out) > 0) {
        result <- as.data.frame(matrix(ncol=11+5, nrow=nrow(out)))
        wh <- c(1,2,3,4,5,7,9,11,12,13,15)
        wh2 <- (1:16)[-wh]
        result[,wh] <- out
        names(result)[wh] <- names(out)
        colnames(result)[wh2] <- rep("pval",5)

        if(!xchr_specific) {
            for(i in 1:5) {
                for(j in 1:nrow(out))
                    result[j,wh2[i]] <- mean(perms[[i]] >= result[j,wh2[i]-1], na.rm=TRUE)
            }
        }
        else { # X-chr-specific
            L <- attr(perms, "L")
            sumL <- sum(L)
            for(j in 1:nrow(out)) {
                this_chrpair_type <- match(chrpair_type[j], c("AA", "AX", "XX"))
                pow <- sumL/L[this_chrpair_type]
                for(i in 1:5) {
                    nominal_p <- mean(perms[[this_chrpair_type]][[i]] >= result[j,wh2[i]-1], na.rm=TRUE)
                    result[j,wh2[i]] <- 1 - (1-nominal_p)^pow # adjusted P-value
                }
            }
        }
        out <- result
    }

    class(out) <- c("summary.scantwo", "data.frame")

    out
}




# subroutine for summary.scantwo; pulls out the key info
#     on each pair of chromosomes
subrousummaryscantwo <-
    function(object, for.perm=FALSE)
{
    lod <- object$lod
    lod[is.na(lod) | lod == Inf | lod == -Inf] <- 0
    map <- object$map

    pos <- map[,2]
    chr <- map[,1]

    tchr <- as.numeric(chr)
    n.chr <- max(tchr)
    xchr <- tapply(map[,4], map[,1], function(a) a[1])
    xchr <- xchr[!is.na(xchr)]
    n.phe <- 1
    if(length(dim(lod)) == 3) n.phe <- dim(lod)[3]

    if(!("scanoneX" %in% names(object)) ||
       is.null(object$scanoneX) || length(object$scanoneX)==0) {
        if(n.phe==1) scanoneX <- diag(lod)
        else {
            if(nrow(lod)==1) scanoneX <- lod[1,1,1]
            else scanoneX <- diag(lod[,,1])
            for(i in 2:n.phe) {
                if(nrow(lod)==1) scanoneX <- cbind(scanoneX, lod[1,1,i])
                else scanoneX <- cbind(scanoneX, diag(lod[,,i]))
            }
        }
    }
    else scanoneX <- object$scanoneX

    if((is.matrix(scanoneX) && nrow(scanoneX) != nrow(lod)) ||
       (!is.matrix(scanoneX) && length(scanoneX) != nrow(lod)))
        stop("scanoneX component has length ", length(scanoneX), " but should have length ", nrow(lod))

    n.chrpair <- n.chr*(n.chr+1)/2

    fill <- matrix(0, nrow=n.chrpair, ncol=n.phe)

    out <- .C("R_summary_scantwo",
              as.integer(nrow(map)),
              as.integer(n.phe),
              as.double(lod),
              as.integer(n.chr),
              as.integer(tchr),
              as.double(pos),
              as.integer(xchr),
              as.double(scanoneX),
              as.integer(n.chrpair),
              chr1=as.integer(rep(0,n.chrpair)),
              chr2=as.integer(rep(0,n.chrpair)),
              as.integer(rep(0,n.chr*n.chr)),
              pos1.jnt=as.double(fill),
              pos2.jnt=as.double(fill),
              pos1.add=as.double(fill),
              pos2.add=as.double(fill),
              pos1.int=as.double(fill),
              pos2.int=as.double(fill),
              jnt.lod.full=as.double(fill),
              jnt.lod.add=as.double(fill),
              add.lod.full=as.double(fill),
              add.lod.add=as.double(fill),
              int.lod.full=as.double(fill),
              int.lod.add=as.double(fill),
              lod.1qtl=as.double(fill),
              PACKAGE="qtl")

    chr1 <- factor(levels(chr)[out$chr1+1], levels=levels(chr))
    chr2 <- factor(levels(chr)[out$chr2+1], levels=levels(chr))

    if(n.phe == 1)
        out <- data.frame(chr1=chr1, chr2=chr2,
                          pos1.jnt=out$pos1.jnt,
                          pos2.jnt=out$pos2.jnt,
                          jnt.lod.full=out$jnt.lod.full,
                          jnt.lod.add=out$jnt.lod.add,
                          pos1.add=out$pos1.add,
                          pos2.add=out$pos2.add,
                          add.lod.full=out$add.lod.full,
                          add.lod.add=out$add.lod.add,
                          pos1.int=out$pos1.int,
                          pos2.int=out$pos2.int,
                          int.lod.full=out$int.lod.full,
                          int.lod.add=out$int.lod.add,
                          lod.1qtl=out$lod.1qtl)
    else
        out <- list(chr1=chr1, chr2=chr2,
                    pos1.jnt=matrix(out$pos1.jnt, ncol=n.phe),
                    pos2.jnt=matrix(out$pos2.jnt, ncol=n.phe),
                    jnt.lod.full=matrix(out$jnt.lod.full, ncol=n.phe),
                    jnt.lod.add=matrix(out$jnt.lod.add, ncol=n.phe),
                    pos1.add=matrix(out$pos1.add, ncol=n.phe),
                    pos2.add=matrix(out$pos2.add, ncol=n.phe),
                    add.lod.full=matrix(out$add.lod.full, ncol=n.phe),
                    add.lod.add=matrix(out$add.lod.add, ncol=n.phe),
                    pos1.int=matrix(out$pos1.int, ncol=n.phe),
                    pos2.int=matrix(out$pos2.int, ncol=n.phe),
                    int.lod.full=matrix(out$int.lod.full, ncol=n.phe),
                    int.lod.add=matrix(out$int.lod.add, ncol=n.phe),
                    lod.1qtl=matrix(out$lod.1qtl, ncol=n.phe))

    if(for.perm) {
        if(n.phe==1) {
            out <- c("full"=max(out$jnt.lod.full,na.rm=TRUE),
                     "fv1"=max(out$jnt.lod.full - out$lod.1qtl, na.rm=TRUE),
                     "int"=max(out$jnt.lod.full - out$add.lod.add, na.rm=TRUE),
                     "add"=max(out$add.lod.add, na.rm=TRUE),
                     "av1"=max(out$add.lod.add - out$lod.1qtl, na.rm=TRUE),
                     "one"=max(out$lod.1qtl, na.rm=TRUE))
        }
        else {
            out <- list("full"=apply(out$jnt.lod.full, 2, max, na.rm=TRUE),
                        "fv1"=apply(out$jnt.lod.full - out$lod.1qtl, 2, max, na.rm=TRUE),
                        "int"=apply(out$jnt.lod.full - out$add.lod.add, 2, max, na.rm=TRUE),
                        "add"=apply(out$add.lod.add, 2, max, na.rm=TRUE),
                        "av1"=apply(out$add.lod.add - out$lod.1qtl, 2, max, na.rm=TRUE),
                        "one"=apply(out$lod.1qtl, 2, max, na.rm=TRUE))
        }
    }

    out
}


print.summary.scantwo <-
    function(x, ...)
{
    if(nrow(x)==0) {
        cat("    There were no pairs of loci meeting the criteria.\n")
        return(invisible(NULL))
    }

    z <- as.character(unlist(x[,1]))

    if(max(nchar(z)) == 1)
        rownames(x) <- apply(x[,1:2], 1, function(a)
                             paste("c", a, collapse=":", sep=""))
    else
        rownames(x) <- apply(x[,1:2], 1, function(a)
                             paste(sprintf("c%-2s", a), collapse=":"))

    x <- x[,-(1:2)]

    cn <- colnames(x)
    if(any(cn=="pos1a"))
        cn[cn=="pos1a"] <- "    pos1a"

    wh <- grep("^pval", cn)
    if(length(wh) > 0)
        cn[wh] <- "pval"
    colnames(x) <- cn

    print.data.frame(x, digits=3)
}


print.summary.addpair <-
    function(x, ...)
{
    if(nrow(x)==0) {
        cat("    There were no pairs of loci meeting the criteria.\n")
        return(invisible(NULL))
    }

    z <- as.character(unlist(x[,1]))

    if(max(nchar(z)) == 1)
        rownames(x) <- apply(x[,1:2], 1, function(a)
                             paste("c", a, collapse=":", sep=""))
    else
        rownames(x) <- apply(x[,1:2], 1, function(a)
                             paste(sprintf("c%-2s", a), collapse=":"))

    x <- x[,-(1:2), drop=FALSE]

    print.data.frame(x, digits=3)
}



print.scantwo <-
    function(x, ...)
{
    d <- dim(x$lod)
    dn <- dimnames(x$lod)[[3]]
    if(is.null(dn)) dn <- attr(x, "phenotypes")


    if(nrow(x$lod) == 0)
        cat("Empty scantwo object.\n")
    else {

        if(length(d)==2)
            print(summary(x))
        else {
            for(i in 1:d[3]) {
                if(is.null(dn)) cat("Phenotype", i, "\n")
                else cat(dn[i], "\n")
                print(summary(x, lod=i))
            }
        }
    }
}


######################################################################
#
# max.scantwo:  Give pair of chromosome with maximum 2-locus LOD score
#
######################################################################

max.scantwo <-
    function(object, lodcolumn=1,
             what=c("best", "full", "add", "int"),
             na.rm=TRUE, ...)
{
    if(class(object)[1] != "scantwo" &&
       class(object)[1] != "scantwocondensed")
        stop("Input must have class \"scantwo\".")

    addpair <- attr(object, "addpair")
    if(!is.null(addpair) && addpair) { # special treatment for output for addpair
        temp <- summary(object)
        mx <- max(temp[,5],na.rm=TRUE)
        return(temp[!is.na(temp[,5]) & temp[,5]==mx,])
    }

    what <- match.arg(what)
    if(class(object)[1] == "scantwo") {
        d <- dim(object$lod)
        if(length(d)==3) {
            if(lodcolumn < 1 || lodcolumn > d[3])
                stop("lodcolumn must be between 1 and ", d[3])
            object$lod <- object$lod[,,lodcolumn]
        }

        out <- subrousummaryscantwo(object, for.perm=FALSE)
    }
    else { # condensed version
        if(is.matrix(object$pos1.jnt)) {
            d <- ncol(object$pos1.jnt)
            if(lodcolumn < 1 || lodcolumn > d)
                stop("lodcolumn must be between 1 and ", d)

            for(i in 3:length(object))
                object[[i]] <- object[[i]][,lodcolumn]
        }
        out <- as.data.frame(object)
    }

    if(what=="best") {
        wh <- which(!is.na(out$jnt.lod.full) & out$jnt.lod.full==max(out$jnt.lod.full, na.rm=TRUE))
        p1.f <- out$pos1.jnt[wh]
        p2.f <- out$pos2.jnt[wh]
        p1.a <- out$pos1.add[wh]
        p2.a <- out$pos2.add[wh]

        lf <- out$jnt.lod.full[wh]
        li <- out$jnt.lod.full[wh] - out$add.lod.add[wh]
        lfv1 <- out$jnt.lod.full[wh] - out$lod.1qtl[wh]
        la <- out$add.lod.add[wh]
        lav1 <- out$add.lod.add[wh] - out$lod.1qtl[wh]
    }
    else if(what=="full") {
        wh <- which(!is.na(out$jnt.lod.full) & out$jnt.lod.full==max(out$jnt.lod.full, na.rm=TRUE))
        p1.f <- p1.a <- out$pos1.jnt[wh]
        p2.f <- p2.a <- out$pos2.jnt[wh]

        lf <- out$jnt.lod.full[wh]
        li <- out$jnt.lod.full[wh] - out$jnt.lod.add[wh]
        lfv1 <- out$jnt.lod.full[wh] - out$lod.1qtl[wh]
        la <- out$jnt.lod.add[wh]
        lav1 <- out$jnt.lod.add[wh] - out$lod.1qtl[wh]
    }
    else if(what=="add") {
        wh <- which(!is.na(out$add.lod.add) & out$add.lod.add==max(out$add.lod.add, na.rm=TRUE))
        p1.f <- p1.a <- out$pos1.add[wh]
        p2.f <- p2.a <- out$pos2.add[wh]

        lf <- out$add.lod.full[wh]
        li <- out$add.lod.full[wh] - out$add.lod.add[wh]
        lfv1 <- out$add.lod.full[wh] - out$lod.1qtl[wh]
        la <- out$add.lod.add[wh]
        lav1 <- out$add.lod.add[wh] - out$lod.1qtl[wh]
    }
    else { # what == "int"
        lod.int <- out$int.lod.full - out$int.lod.add
        wh <- which(!is.na(lod.int) & lod.int == max(lod.int, na.rm=TRUE))
        p1.f <- p1.a <- out$pos1.int[wh]
        p2.f <- p2.a <- out$pos2.int[wh]

        lf <- out$int.lod.full[wh]
        li <- out$int.lod.full[wh] - out$int.lod.add[wh]
        lfv1 <- out$int.lod.full[wh] - out$lod.1qtl[wh]
        la <- out$int.lod.add[wh]
        lav1 <- out$int.lod.add[wh] - out$lod.1qtl[wh]
    }

    out <- data.frame(chr1=out$chr1[wh], chr2=out$chr2[wh],
                      pos1f=p1.f,
                      pos2f=p2.f,
                      lod.full=lf, lod.fv1=lfv1, lod.int=li,
                      pos1a=p1.a,
                      pos2a=p2.a,
                      lod.add=la, lod.av1=lav1)

    if(what != "best") {
        out <- out[,-(8:9)]
        names(out)[3:4] <- c("pos1", "pos2")
    }

    class(out) <- c("summary.scantwo", "data.frame")
    rownames(out) <- what

    out
}



######################################################################
# clean.scantwo
#
# sets LOD scores that are missing or < 0 to 0
# If full LOD < add've LOD, set full = add've
# sets LOD scores, for pairs of positions that are not separated
#     by n.mar markers and distance cM, to 0
######################################################################

clean.scantwo <-
    function(object, n.mar=1, distance=0, ...)
{
    if(class(object)[1] != "scantwo")
        stop("Input should have class \"scantwo\".")

    addpair <- attr(object, "addpair")
    if(is.null(addpair)) addpair <- FALSE

    lod <- object$lod
    map <- object$map
    if(!("fullmap" %in% names(attributes(object))))
        stop("clean.scantwo only works on scantwo objects created with R/qtl ver >= 1.04-38.\n")
    fmap <- attr(object, "fullmap")

    lod[is.na(lod) | lod < 0] <- 0

    subrou <- function(x,y,z)
    {
        out <- x
        for(i in seq(along=x))
            out[i] <- sum((z > x[i] & z < y[i]) | (z < x[i] & z > y[i])) +
                (x[i] == y[i])
        out
    }

    last <- 0
    for(i in seq(along=fmap)) {
        m <- map[map[,1]==names(fmap)[i],2]
        idx <- 1:length(m)+last
        last <- last + length(m)

        toclean <- (outer(m, m, subrou, fmap[[i]]) < n.mar) | (abs(outer(m, m, "-")) <distance)
        diag(toclean) <- FALSE

        if(length(dim(lod))==3) # case of multiple phenotypes
            lod[idx,idx,][toclean] <- 0
        else
            lod[idx,idx][toclean] <- 0
    }

    # if full LOD < add've LOD, set full = add've
    if(!addpair) {
        if(length(dim(lod)) == 2) {
            lo <- lower.tri(lod)
            wh <- (lod[lo] < t(lod)[lo])
            if(any(wh))
                lod[lo][wh] <- t(lod)[lo][wh]
        }
        else {
            lo <- lower.tri(lod[,,1])
            for(i in 1:dim(lod)[3]) {
                wh <- (lod[,,i][lo] < t(lod[,,i])[lo])
                if(any(wh))
                    lod[,,i][lo][wh] <- t(lod[,,i])[lo][wh]
            }
        }
    }

    object$lod <- lod

    object
}

######################################################################
#
# subset.scantwo
#
######################################################################
subset.scantwo <-
    function(x, chr, lodcolumn, ...)
{
    if(class(x)[1] != "scantwo")
        stop("Input should have class \"scantwo\".")

    if((missing(chr) || length(chr)==0) && missing(lodcolumn)) return(x)

    if(!missing(lodcolumn)) {
        if(any(lodcolumn>0) && any(lodcolumn<0))
            stop("lodcolumn values can't be both >0 and <0.")

        if(any(lodcolumn<0) || is.logical(lodcolumn))
            lodcolumn <- (1:(dim(x$lod)[3]))[lodcolumn]
        if(length(lodcolumn)==0)
            stop("You must retain at least one LOD column.")
        if(any(lodcolumn < 1 || lodcolumn > dim(x$lod)[3]))
            stop("lodcolumn values must be >=1 and <=",dim(x$lod)[3])

        x$lod <- x$lod[,,lodcolumn]
        if("scanoneX" %in% names(x))
            x$scanoneX <- x$scanoneX[,lodcolumn]
    }

    if(!missing(chr)) {
        chr <- matchchr(chr, unique(x$map[,1]))

        wh <- x$map[,1] %in% chr

        x$map <- x$map[wh,,drop=FALSE]
        x$map[,1] <- droplevels(x$map[,1])

        if(length(dim(x$lod))==2)
            x$lod <- x$lod[wh,wh,drop=FALSE]
        else
            x$lod <- x$lod[wh,wh,,drop=FALSE]

        if(!is.null(x$scanoneX))
            x$scanoneX <- x$scanoneX[wh]

        if("fullmap" %in% names(attributes(x))) {
            fmap <- attr(x, "fullmap")
            fmap <- fmap[chr]
            attr(x, "fullmap") <- fmap
        }

    }

    x
}

######################################################################
# summary.scantwoperm
#
# Give genome-wide LOD thresholds on the basis of results of
# scantwo permutation test (from scantwo with n.perm > 0)
######################################################################
summary.scantwoperm <-
    function(object, alpha=c(0.05, 0.10), ...)
{
    if(!any(class(object) == "scantwoperm"))
        stop("Input should have class \"scantwoperm\".")

    if("AA" %in% names(object)) { # X-chr-specific version
        # get region-specific (nominal) signif levels
        L <- attr(object, "L")
        LL <- attr(object, "LL")
        if(is.null(LL)) stop("LL attribute not found in input object")
        if(length(alpha)==1) {
            tmp <- L/sum(L); tmp <- c(tmp[1], 1, tmp[2])
            one_minus_alpha_onechr <- cbind((1-alpha)^tmp)
            one_minus_alpha <- cbind((1-alpha)^(LL/sum(LL)))
        }
        else {
            tmp <- L/sum(L); tmp <- c(tmp[1], 1, tmp[2])
            one_minus_alpha_onechr <- vapply(alpha, function(a,b) (1-a)^b, rep(0, length(tmp)), tmp)
            one_minus_alpha <- vapply(alpha, function(a,b) (1-a)^b, rep(0, length(LL)), LL/sum(LL))
        }
        # get quantiles
        out <- vector("list", length(object))
        names(out) <- names(object)
        for(i in seq(along=out)) {
            f <- function(a, qu) {
                b <- apply(a, 2, quantile, qu)
                if(!is.matrix(b)) {
                    nam <- names(b)
                    b <- matrix(b, nrow=1)
                    colnames(b) <- nam
                }
                rownames(b) <- paste0(alpha*100, "%")
                b
            }
            out[[i]] <- lapply(unclass(object[[i]])[1:5], f, one_minus_alpha[i,,drop=FALSE])
            qu_one <- f(object[[i]][[6]], one_minus_alpha_onechr[i,,drop=FALSE])
            out[[i]] <- c(out[[i]], "one"=list(qu_one))
        }

        attr(out, "n.perm") <- vapply(object, function(a) nrow(a[[1]]), 0)
        class(out) <- c("summary.scantwoperm", "list")

        return(out)
    }

    out <- lapply(object, apply, 2, quantile, 1-alpha)

    for(i in 1:length(out)) {
        if(!is.matrix(out[[i]]))
            out[[i]] <- matrix(out[[i]], nrow=length(alpha))
        rownames(out[[i]]) <- paste0(alpha*100, "%")
    }

    attr(out, "n.perm") <- nrow(object[[1]])
    class(out) <- c("summary.scantwoperm", "list")
    out
}

print.summary.scantwoperm <-
    function(x, ...)
{
    if("AA" %in% names(x)) { # X-sp perms
        nperm <- attr(x, "n.perm")
        rn <- rownames(x[[1]][[1]])
        nc <- rapply(x, ncol)
        if(length(unique(nc)) != 1)
            stop("The components shouldn't have varying numbers of columns.\n")
        nc <- nc[1]
        phe <- colnames(x[[1]][[1]])

        convert2matrix <-
            function(a, which_col=1) {
                a <- lapply(a, "[", , which_col, drop=FALSE)
                b <- matrix(unlist(a), ncol=6)
                rownames(b) <- rn
                colnames(b) <- names(a)
                b
            }

        for(i in seq(along=phe)) {
            cat(phe[i], ":\n", sep="")
            y <- lapply(x, convert2matrix, i)
            lab <- c(AA="A:A", AX="A:X", XX="X:X")
            for(j in c("AA", "AX", "XX")) {
                cat("  ", lab[j], " (", nperm[j], " permutations)\n", sep="")
                print(y[[j]], digits=3)
                cat("\n")
            }
            if(i != length(phe)) cat("\n")
        }
        invisible(return(x))
    }

    nam <- names(x)
    rn <- rownames(x[[1]])
    n.perm <- attr(x, "n.perm")

    nc <- sapply(x, ncol)
    if(length(unique(nc)) != 1)
        stop("The components shouldn't have varying numbers of columns.\n")
    nc <- nc[1]

    if(nc==1) {
        phe <- colnames(x[[1]])
        if(is.null(phe)) phe <- ""

        x <- matrix(unlist(x), nrow=length(rn))
        dimnames(x) <- list(rn, nam)

        if(is.null(phe))
            cat("(", n.perm, " permutations)\n", sep="")
        else
            cat(phe, " (", n.perm, " permutations)\n", sep="")
        print(x, digits=3)
    }
    else {
        phe <- colnames(x[[1]])
        if(is.null(phe)) phe <- paste("pheno", 1:nc)

        for(i in 1:nc) {
            y <- matrix(ncol=length(x), nrow=nrow(x[[1]]))
            for(j in 1:length(x))
                y[,j] <- x[[j]][,i]
            dimnames(y) <- list(rn, nam)

            cat(phe[i], " (", n.perm, " permutations)\n", sep="")
            print(y, digits=3)
            if(i != nc) cat("\n")
        }
    }
}

######################################################################
# combine scantwo results ... paste the phenotype columns together
######################################################################
cbind.scantwo <- c.scantwo <-
    function(...)
{
    dots <- list(...)

    cl1 <- class(dots[[1]])
    if(length(dots)==1 && length(cl1)==1 && cl1=="list") dots <- dots[[1]]

    if(length(dots)==1) return(dots[[1]])
    for(i in seq(along=dots)) {
        if(!any(class(dots[[i]]) == "scantwo"))
            stop("Input should have class \"scantwo\".")
    }

    # check dimensions of LODs
    nr <- sapply(dots, function(a) nrow(a$lod))
    nc <- sapply(dots, function(a) ncol(a$lod))
    nd3 <- sapply(dots, function(a) dim(a$lod)[3])
    if(any(nr[-1]!=nr[1]) || any(nc[-1] != nc[1]))
        stop("Input objects are not the same dimensions.")

    # check maps
    map1 <- dots$map[[1]]
    for(i in 2:length(dots)) {
        if(!all(dots$map[[i]] != map1))
            stop("Maps are not all the same.")
    }

    output <- dots[[1]]
    nd3[is.na(nd3)] <- 1
    output$lod <- array(dim=c(nr[1], nc[1], sum(nd3)))
    start <- cumsum(c(0,nd3))[-length(nd3)-1]
    end <- start+nd3

    for(i in seq(along=dots))
        output$lod[,,(start[i]+1):end[i]] <- dots[[i]]$lod
    attr(output, "phenotypes") <- unlist(lapply(dots, attr, "phenotypes"))
    dimnames(output$lod) <- list(NULL, NULL, attr(output, "phenotypes"))

    # check scanoneX
    if(is.null(dots[[1]]$scanoneX)) {
        for(i in 2:length(dots)) {
            if(!is.null(dots[[i]]$scanoneX))
                stop("Some but not all input objects have null scanoneX.")
        }
    }
    else {
        nrx <- sapply(dots, function(a) nrow(a$scanoneX))
        if(any(nrx[-1]!=nrx[1]))
            stop("Mismatch in scanoneX dimensions.")
        for(i in 2:length(dots))
            output$scanoneX <- cbind(output$scanoneX, dots[[i]]$scanoneX)
        colnames(output$scanoneX) <- attr(output, "phenotypes")
    }

    methods <- sapply(dots, attr, "method")
    if(length(unique(methods)) != 1)
        attr(output, "method") <- rep(sapply(dots, attr, "method"), nd3)

    fm <- lapply(dots, attr, "fullmap")
    for(i in 2:length(fm)) {
        if(length(fm[[1]]) != length(fm[[i]]) || !all(names(fm[[1]]) == names(fm[[i]])))
            stop("Mismatch in \"fullmap\" attributes (1).")
        for(j in 1:length(fm[[1]])) {
            if(length(fm[[1]][[j]]) != length(fm[[i]][[j]]) ||
               !all(names(fm[[1]][[j]]) == names(fm[[i]][[j]])) ||
               max(abs(fm[[1]][[j]] - fm[[i]][[j]])) > 0.001)
                stop("Mismatch in \"fullmap\" attributes (2).")
        }
    }

    output
}


######################################################################
# combine scantwoperm results ... paste the rows together
######################################################################
rbind.scantwoperm <- c.scantwoperm <-
    function(...)
{
    dots <- list(...)

    cl1 <- class(dots[[1]])
    if(length(dots)==1 && length(cl1)==1 && cl1=="list") dots <- dots[[1]]
    if(length(dots)==1) return(dots[[1]])

    xchrsp <- vapply(dots, function(a) "AA" %in% names(a), TRUE)
    if(any(xchrsp)) {
        if(!all(xchrsp)) stop("Some but not all inputs are X-chr specific")

        for(i in 2:length(dots)) {
            for(j in seq(along=dots[[1]])) {
                for(k in seq(along=dots[[1]][[j]])) {
                    if(ncol(dots[[1]][[j]][[k]]) != ncol(dots[[i]][[j]][[k]]))
                        stop("Mismatch in no. columns")
                    if(any(colnames(dots[[1]][[j]][[k]]) != colnames(dots[[1]][[j]][[k]])))
                        warning("Mismatch in column names")
                    dots[[1]][[j]][[k]] <- rbind(dots[[1]][[j]][[k]], dots[[i]][[j]][[k]])
                }
            }
        }
        return(dots[[1]])
    }

    for(i in seq(along=dots)) {
        if(!any(class(dots[[i]]) == "scantwoperm"))
            stop("Input should have class \"scantwoperm\".")
    }

    nc <- sapply(dots, function(a) ncol(a[[1]]))
    if(length(unique(nc)) != 1)
        stop("Number of LOD columns in the input objects must be constant.\n")

    flag <- 0
    for(j in 1:length(dots[[1]])) {
        for(i in 2:length(dots)) {
            dots[[1]][[j]] <- rbind(dots[[1]][[j]], dots[[i]][[j]])
            if(any(colnames(dots[[i]][[j]]) != colnames(dots[[1]][[j]]))) flag <- 1
        }
    }
    if(flag) warning("Mismatch in column names; input may not be consistent.\n")

    dots[[1]]
}

# paste columns together
cbind.scantwoperm <-
function(...)
{
    dots <- list(...)
    cl1 <- class(dots[[1]])
    if(length(dots)==1 && length(cl1)==1 && cl1=="list") dots <- dots[[1]]

    if(length(dots)==1) return(dots)

    xchrsp <- vapply(dots, function(a) "AA" %in% names(a), TRUE)
    if(any(xchrsp)) {
        if(!all(xchrsp)) stop("Some but not all inputs are X-chr specific")

        for(i in 2:length(dots)) {
            for(j in seq(along=dots[[1]])) {
                for(k in seq(along=dots[[1]][[j]])) {
                    if(nrow(dots[[1]][[j]][[k]]) != nrow(dots[[i]][[j]][[k]]))
                        stop("Mismatch in no. permutations")
                    dots[[1]][[j]][[k]] <- cbind(dots[[1]][[j]][[k]], dots[[i]][[j]][[k]])
                }
            }
        }
        return(dots[[1]])
    }

    for(i in 2:length(dots)) {
        for(j in seq(along=dots[[1]])) {
            if(nrow(dots[[1]][[j]]) != nrow(dots[[i]][[j]]))
                stop("Mismatch in no. permutations")
            dots[[1]][[j]] <- cbind(dots[[1]][[j]], dots[[i]][[j]])
        }
    }

    dots[[1]]
}



######################################################################
# condensed scantwo output

condense <-
    function(object)
    UseMethod("condense")

condense.scantwo <-
    function(object)
{
    out <- subrousummaryscantwo(object, for.perm=FALSE)
    class(out) <- c("scantwocondensed", "list")
    out
}

summary.scantwocondensed <- summary.scantwo
max.scantwocondensed <- max.scantwo


##############################
# subset.scantwoperm: pull out a set of lodcolumns
##############################
subset.scantwoperm <-
    function(x, repl, lodcolumn, ...)
{
    if("AA" %in% names(x)) { # x chr specific
        for(j in seq(along=x)) {
            if(missing(lodcolumn)) lodcolumn <- 1:ncol(x[[j]][[1]])
            else if(!is.null(attr(try(x[[j]][[1]][,lodcolumn], silent=TRUE),"try-error")))
                stop("lodcolumn misspecified.")

            repl <- 1:nrow(x[[j]][[1]])

            x[[j]] <- lapply(x[[j]], function(a,b,d) unclass(a)[b,d,drop=FALSE], repl, lodcolumn)
        }
        return(x)
    }

    att <- attributes(x)

    if(any(!sapply(x, is.matrix)))
        x <- lapply(x, as.matrix)

    if(missing(lodcolumn)) lodcolumn <- 1:ncol(x[[1]])
    else if(!is.null(attr(try(x[[1]][,lodcolumn], silent=TRUE),"try-error")))
        stop("lodcolumn misspecified.")

    if(missing(repl)) repl <- 1:nrow(x[[1]])
    else if(!is.null(attr(try(x[[1]][repl,], silent=TRUE),"try-error")))
        stop("repl misspecified.")

    cl <- class(x)
    x <- lapply(x, function(a,b,d) unclass(a)[b,d,drop=FALSE], repl, lodcolumn)
    class(x) <- cl

    for(i in seq(along=att)) {
        if(names(att)[i] == "dim" || length(grep("names", names(att)[i]))>0) next
        attr(x, names(att)[i]) <- att[[i]]
    }

    x
}

# subset.scantwoperm using [,]
`[.scantwoperm` <-
    function(x, repl, lodcolumn)
    subset.scantwoperm(x, repl, lodcolumn)

# end of summary.scantwo.R
