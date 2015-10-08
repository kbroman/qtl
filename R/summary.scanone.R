######################################################################
#
# summary.scanone.R
#
# copyright (c) 2001-2015, Karl W Broman
# last modified Oct, 2015
# first written Sep, 2001
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
# Contains: summary.scanone, print.summary.scanone,
#           max.scanone, c.scanone, subset.scanone,
#           summary.scanoneperm, print.summary.scanoneperm
#           c.scanoneperm, rbind.scanoneperm, cbind.scanoneperm
#           grab.arg.names, subset.scanoneperm, [.scanoneperm
#
######################################################################

##################################################################
# summarize scanone results
##################################################################
summary.scanone <-
    function(object, threshold, format=c("onepheno", "allpheno", "allpeaks", "tabByCol", "tabByChr"),
             perms, alpha, lodcolumn=1, pvalues=FALSE,
             ci.function=c("lodint", "bayesint"), ...)
{
    if(!any(class(object) == "scanone"))
        stop("Input should have class \"scanone\".")

    format <- match.arg(format)
    ncol.object <- ncol(object)-2
    cn.object <- colnames(object)[-(1:2)]

    if(ncol.object==1 && (format == "allpeaks" || format == "allpheno")) {
        warning("With just one LOD column, format=\"onepheno\" used.")
        format <- "onepheno"
    }

    if(format != "onepheno" && !missing(lodcolumn))
        warning("lodcolumn ignored except when format=\"onepheno\".")

    if(!missing(perms)) {
        if("scantwoperm" %in% class(perms))
            perms <- scantwoperm2scanoneperm(perms)
        else if(!("scanoneperm" %in% class(perms)))
            warning("perms need to be in scanoneperm format.")
    }

    # check input
    if(missing(perms) && !missing(alpha))
        stop("If alpha is to be used, permutation results must be provided.")
    if(!missing(threshold) && !missing(alpha))
        stop("Only one of threshold and alpha should be specified.")

    if(format == "onepheno") {
        if(!missing(lodcolumn) && length(lodcolumn) > 1) {
            warning("With format=\"onepheno\", lodcolumn should have length 1.")
            lodcolumn <- lodcolumn[1]
        }
        if(lodcolumn < 1 || lodcolumn > ncol.object)
            stop("lodcolumn should be between 1 and no. LOD columns.")
    }

    if(!missing(alpha) && length(alpha) > 1) {
        warning("alpha should have length 1.")
        alpha <- alpha[1]
    }

    if(!missing(perms)) {
        if("xchr" %in% names(attributes(perms))) {
            ncol.perms <- ncol(perms$A)
            cn.perms <- colnames(perms$A)
        }
        else {
            ncol.perms <- ncol(perms)
            cn.perms <- colnames(perms)
        }

        if(ncol.object != ncol.perms) {
            if(ncol.perms==1) { # reuse the multiple columns
                origperms <- perms
                if("xchr" %in% names(attributes(perms))) {
                    for(j in 2:ncol.object) {
                        perms$A <- cbind(perms$A, origperms$A)
                        perms$X <- cbind(perms$X, origperms$X)
                    }
                    cn.perms <- colnames(perms$A) <- colnames(perms$X) <- cn.object
                }
                else {
                    for(j in 2:ncol.object)
                        perms <- cbind(perms, origperms)
                    cn.perms <- colnames(perms) <- cn.object
                }
                warning("Just one column of permutation results; reusing for all LOD score columns.")
            }
            else {
                if(ncol.object == 1) {
                    warning("Using just the first column in the perms input")
                    if("xchr" %in% names(attributes(perms))) {
                        perms$A <- perms$A[,1,drop=FALSE]
                        perms$X <- perms$X[,1,drop=FALSE]
                    }
                    else {
                        clp <- class(perms)
                        perms <- perms[,1,drop=FALSE]
                        class(perms) <- clp
                    }
                }
                else
                    stop("scanone input has different number of LOD columns as perms input.")
            }
        }
        if(!all(cn.object == cn.perms))
            warning("Column names in scanone input do not match those in perms input.")
    }

    if(format != "onepheno") {
        if(!missing(threshold)) {
            if(length(threshold)==1)
                threshold <- rep(threshold, ncol.object)
            else if(length(threshold) != ncol.object)
                stop("threshold should have length 1 or match number LOD scores in scanone input.")
        }
    }

    if(missing(perms) && pvalues) {
        warning("Can show p-values only if perms are provided.")
        pvalues <- FALSE
    }
    # end of check of input

    # chromosome IDs as a character string
    chr <- as.character(object[,1])

    if(format=="onepheno") {
        lodcolumn <- lodcolumn+2

        # pull out max on each chromosome
        wh <- NULL
        for(i in unique(chr)) {
            if(any(!is.na(object[chr==i,lodcolumn]))) {
                mx <- max(object[chr==i,lodcolumn],na.rm=TRUE)
                tmp <- which(chr==i & object[,lodcolumn]==mx)
                if(length(tmp) > 1) tmp <- sample(tmp, 1) # if multiple, pick at random
                wh <- c(wh, tmp)
            }
        }
        thechr <- as.character(object[wh,1])

        if(!missing(threshold))
            wh <- wh[object[wh,lodcolumn] > threshold]
        else if(!missing(alpha)) {
            thr <- summary(perms, alpha)
            if("xchr" %in% names(attributes(perms))) {
                thr <- sapply(thr, function(a,b) a[,b], lodcolumn-2)
                xchr <- attr(perms, "xchr")
                xchr <- names(xchr)[xchr]
                xchr <- thechr %in% xchr
                wh <- wh[(!xchr & object[wh,lodcolumn] > thr[1]) |
                         (xchr & object[wh,lodcolumn] > thr[2])]
            }
            else {
                thr <- thr[,lodcolumn-2]
                wh <- wh[object[wh,lodcolumn] > thr]
            }
        }

        result <- object[wh,]
    } # end of "onepheno" format

    else if(format=="allpheno") {
        # pull out max on each chromosome
        wh <- vector("list", ncol.object)
        for(lodcolumn in 1:ncol.object+2) {
            for(i in unique(chr)) {
                if(any(!is.na(object[chr==i,lodcolumn]))) {
                    mx <- max(object[chr==i,lodcolumn],na.rm=TRUE)
                    tmp <- which(chr==i & object[,lodcolumn]==mx)
                    if(length(tmp) > 1) tmp <- sample(tmp, 1)
                    wh[[lodcolumn-2]] <- c(wh[[lodcolumn-2]], tmp)
                }
            }
        }

        if(!missing(threshold)) { # rows with at least one LOD > threshold
            for(lodcolumn in 1:ncol.object) {
                temp <- wh[[lodcolumn]]
                wh[[lodcolumn]] <- temp[object[temp,lodcolumn+2] > threshold[lodcolumn]]
            }
        }
        else if(!missing(alpha)) {
            thr <- summary(perms, alpha)
            if("xchr" %in% names(attributes(perms))) {
                xchr <- attr(perms, "xchr")
                xchr <- names(xchr)[xchr]

                for(lodcolumn in 1:ncol.object) {
                    temp <- wh[[lodcolumn]]
                    thechr <- as.character(object[temp,1])
                    xchr <- thechr %in% xchr
                    wh[[lodcolumn]] <- temp[(!xchr & object[temp,lodcolumn+2] > thr$A[lodcolumn]) |
                                            (xchr & object[temp,lodcolumn+2] > thr$X[lodcolumn])]
                }
            }
            else {
                for(lodcolumn in 1:ncol.object) {
                    temp <- wh[[lodcolumn]]
                    wh[[lodcolumn]] <- temp[object[temp,lodcolumn+2] > thr[lodcolumn]]
                }
            }
        }

        wh <- sort(unique(unlist(wh)))
        result <- object[wh,]
    }  # end of format=="allpheno"

    else if(format=="allpeaks") {
        # pull out max on each chromosome
        wh <- vector("list", ncol.object)

        for(lodcolumn in (1:ncol.object)+2) {
            for(i in unique(chr)) {
                if(any(!is.na(object[chr==i,lodcolumn]))) {
                    mx <- max(object[chr==i,lodcolumn],na.rm=TRUE)
                    temp <- which(chr==i & object[,lodcolumn]==mx)
                    if(length(temp)>1) temp <- sample(temp, 1)
                    wh[[lodcolumn-2]] <- c(wh[[lodcolumn-2]], temp)
                }
                else
                    wh[[lodcolumn-2]] <- c(wh, NA)
            }
        }

        pos <- sapply(wh, function(a,b) b[a], object[,2])
        if(!is.matrix(pos)) pos <- as.matrix(pos)

        lod <- pos

        for(i in 1:ncol(pos))
            lod[,i] <- object[wh[[i]],i+2]
        thechr <- as.character(unique(object[,1]))

        if(!missing(threshold)) { # rows with at least one LOD > threshold
            keep <- NULL
            for(i in seq(along=thechr))
                if(any(lod[i,] > threshold)) keep <- c(keep, i)
        }
        else if(!missing(alpha)) {
            keep <- NULL
            thr <- summary(perms, alpha)
            if("xchr" %in% names(attributes(perms))) {
                xchr <- attr(perms, "xchr")
                xchr <- names(xchr)[xchr]
                xchr <- thechr %in% xchr

                for(i in seq(along=thechr)) {
                    if((xchr[i] && any(lod[i,] > thr$X)) ||
                       (!xchr[i] && any(lod[i,] > thr$A)))
                        keep <- c(keep, i)
                }
            }
            else {
                for(i in seq(along=thechr)) {
                    if(any(lod[i,] > thr))
                        keep <- c(keep, i)
                }
            }
        }
        else keep <- seq(along=thechr)

        if(is.null(keep))
            result <- object[NULL,,drop=FALSE]
        else {
            pos <- pos[keep,,drop=FALSE]
            lod <- lod[keep,,drop=FALSE]
            thechr <- thechr[keep]
            result <- as.data.frame(matrix(ncol=ncol.object*2+1,nrow=length(keep)), stringsAsFactors=TRUE)
            names(result)[1] <- "chr"
            names(result)[(1:ncol.object)*2] <- "pos"
            names(result)[(1:ncol.object)*2+1] <- names(object)[-(1:2)]
            result[,1] <- thechr
            result[,(1:ncol.object)*2] <- pos
            result[,(1:ncol.object)*2+1] <- lod
        }
    }
    else { # format=="tabByChr" or =="tabByCol"
        result <- vector("list", ncol.object)
        names(result) <- names(object)[-(1:2)]

        # pull out max on each chromosome
        wh <- vector("list", ncol.object)

        for(lodcolumn in (1:ncol.object)+2) {
            for(i in unique(chr)) {
                if(any(!is.na(object[chr==i,lodcolumn]))) {
                    mx <- max(object[chr==i,lodcolumn],na.rm=TRUE)
                    temp <- which(chr==i & object[,lodcolumn]==mx)
                    if(length(temp)>1) temp <- sample(temp, 1)
                    wh[[lodcolumn-2]] <- c(wh[[lodcolumn-2]], temp)
                }
                else
                    wh[[lodcolumn-2]] <- c(wh, NA)
            }
        }

        pos <- sapply(wh, function(a,b) b[a], object[,2])
        if(!is.matrix(pos)) pos <- as.matrix(pos)

        lod <- pos

        for(i in 1:ncol(pos))
            lod[,i] <- object[wh[[i]],i+2]
        thechr <- as.character(unique(object[,1]))

        for(i in 1:ncol.object)
            result[[i]] <- object[wh[[i]],c(1,2,i+2)]
    }


    if(pvalues) {
        if(format != "tabByCol" && format != "tabByChr") {
            if(nrow(result) > 0) { # get p-values and add to the results
                rn <- rownames(result)

                if("xchr" %in% names(attributes(perms))) {
                    xchr <- attr(perms, "xchr")
                    xchr <- names(xchr)[xchr]
                    xchr <- as.character(result[,1]) %in% xchr
                    L <- attr(perms, "L")
                    Lt <- sum(L)

                    if(format=="allpeaks")
                        thecol <- (1:ncol.object)*2+1
                    else thecol <- (1:ncol.object)+2

                    if(any(xchr)) {
                        tempX <- calcPermPval(result[xchr,thecol,drop=FALSE], perms$X)
                        tempX <- as.data.frame(1-(1-tempX)^(Lt/L[2]), stringsAsFactors=TRUE)
                    }
                    else tempX <- NULL
                    if(any(!xchr)) {
                        tempA <- calcPermPval(result[!xchr,thecol,drop=FALSE], perms$A)
                        tempA <- as.data.frame(1-(1-tempA)^(Lt/L[1]), stringsAsFactors=TRUE)
                    }
                    else tempA <- NULL
                    pval <- rbind(tempA, tempX)
                    if(any(xchr)) pval[xchr,] <- tempX
                    if(any(!xchr)) pval[!xchr,] <- tempA
                }
                else {
                    if(format=="allpeaks") thecol <- (1:ncol.object)*2+1
                    else thecol <- (1:ncol.object)+2

                    pval <- as.data.frame(calcPermPval(result[,thecol,drop=FALSE], perms), stringsAsFactors=TRUE)
                }

                if(format == "allpeaks") {
                    temp <- as.data.frame(matrix(nrow=nrow(result), ncol=ncol.object*3+1), stringsAsFactors=TRUE)

                    names(temp)[1] <- names(result)[1]
                    temp[,1] <- result[,1]

                    for(i in 1:ncol.object) {
                        names(temp)[i*3+(-1:1)] <- c(names(result)[i*2+(0:1)], "pval")
                        temp[,i*3-1:0] <- result[,i*2+(0:1)]
                        temp[,i*3+1] <- pval[[i]]
                    }
                }
                else if(format != "tabByCol" && format != "tabByChr") {
                    temp <- as.data.frame(matrix(nrow=nrow(result), ncol=ncol.object*2+2), stringsAsFactors=TRUE)

                    names(temp)[1:2] <- names(result)[1:2]
                    temp[,1:2] <- result[,1:2]
                    for(i in 1:ncol.object) {
                        names(temp)[i*2+1:2] <- c(names(result)[i+2], "pval")
                        temp[,i*2+1] <- result[,i+2]
                        temp[,i*2+2] <- pval[[i]]
                    }
                }
                result <- temp
                rownames(result) <- rn

            }
        }
        else { # format=="tabByCol" || format=="tabByChr"
            peaks <- as.data.frame(lapply(result, function(a) a[,3]), stringsAsFactors=TRUE)

            if("xchr" %in% names(attributes(perms))) {
                xchr <- attr(perms, "xchr")
                xchr <- names(xchr)[xchr]
                xchr <- as.character(result[[1]][,1]) %in% xchr
                L <- attr(perms, "L")
                Lt <- sum(L)

                if(any(xchr)) {
                    tempX <- as.data.frame(calcPermPval(peaks[xchr,,drop=FALSE], perms$X), stringsAsFactors=TRUE)
                    tempX <- 1-(1-tempX)^(Lt/L[2])
                }
                else tempX <- NULL
                if(any(!xchr)) {
                    tempA <- as.data.frame(calcPermPval(peaks[!xchr,,drop=FALSE], perms$A), stringsAsFactors=TRUE)
                    tempA <- 1-(1-tempA)^(Lt/L[1])
                }
                else tempA <- NULL

                pval <- rbind(tempA, tempX)
                if(any(xchr)) pval[xchr,] <- tempX
                if(any(!xchr)) pval[!xchr,] <- tempA
            }
            else
                pval <- as.data.frame(calcPermPval(peaks, perms), stringsAsFactors=TRUE)

            for(i in seq(along=result))
                result[[i]] <- cbind(as.data.frame(result[[i]]), pval=pval[,i], stringsAsFactors=TRUE)
        }
    }

    if(format=="tabByCol" || format=="tabByChr") {

        # drop insignificant peaks
        if(!missing(threshold)) { # rows with at least one LOD > threshold
            for(i in seq(along=result))
                result[[i]] <- result[[i]][lod[,i] > threshold[i],,drop=FALSE]
        }
        else if(!missing(alpha)) {
            keep <- NULL
            thr <- summary(perms, alpha)

            if("xchr" %in% names(attributes(perms))) {
                xchr <- attr(perms, "xchr")
                xchr <- names(xchr)[xchr]
                xchr <- thechr %in% xchr

                for(i in seq(along=result))
                    result[[i]] <- result[[i]][(lod[,i] > thr$A[i] & !xchr) |
                                               (lod[,i] > thr$X[i] & xchr), , drop=FALSE]
            }
            else {
                for(i in seq(along=result))
                    result[[i]] <- result[[i]][lod[,i] > thr[i],,drop=FALSE]
            }
        }

        # add intervals
        ci.function <- match.arg(ci.function)
        if(ci.function=="lodint") cif <- lodint
        else cif <- bayesint

        for(i in seq(along=result)) {
            if(nrow(result[[i]]) == 0) next

            lo <- hi <- rep(NA, nrow(result[[i]]))
            for(j in 1:nrow(result[[i]])) {
                temp <- cif(object, chr=as.character(result[[i]][j,1]), lodcolumn=i, ...)
                lo[j] <- temp[1,2]
                hi[j] <- temp[nrow(temp),2]
            }
            result[[i]] <- cbind(as.data.frame(result[[i]]), ci.low=lo, ci.high=hi, stringsAsFactors=TRUE)
            colnames(result[[i]])[3] <- "lod"
        }

        if(format=="tabByChr" && length(result)==1)
            format <- "tabByCol"     # no need to do by chr in this case

        if(format=="tabByChr") {
            temp <- vector("list", length(thechr))
            names(temp) <- thechr

            for(i in seq(along=result)) {
                if(nrow(result[[i]])==0) next
                rownames(result[[i]]) <- paste(names(result)[i], rownames(result[[i]]), sep=" : ")
                for(j in 1:nrow(result[[i]])) {
                    thischr <- match(result[[i]][j,1], thechr)
                    if(length(temp[[thischr]])==0)
                        temp[[thischr]] <- result[[i]][j,,drop=FALSE]
                    else
                        temp[[thischr]] <- rbind(temp[[thischr]], result[[i]][j,,drop=FALSE])
                }
            }
            result <- temp
        }

        # move CI to before the lod score
        for(i in seq(along=result)) {
            if(is.null(result[[i]]) || nrow(result[[i]])==0) next
            nc <- ncol(result[[i]])
            result[[i]] <- result[[i]][,c(1,2,nc-1,nc,3:(nc-2)),drop=FALSE]
        }

        attr(result, "tab") <- format
    }

    if(format=="allpeaks") rownames(result) <- as.character(result$chr)

    if(format=="tabByCol" || format=="tabByChr")
        class(result) <- c("summary.scanone", "list")
    else
        class(result) <- c("summary.scanone", "data.frame")
    result
}

# print output of summary.scanone
print.summary.scanone <-
    function(x, ...)
{
    tab <- attr(x, "tab")

    if(is.null(tab) && nrow(x) == 0) {
        cat("    There were no LOD peaks above the threshold.\n")
        return(invisible(NULL))
    }

    flag <- FALSE
    if(is.null(tab)) {
        print.data.frame(x,digits=3)
        flag <- TRUE
    }
    else if(tab=="tabByChr") {
        for(i in seq(along=x)) {
            if(is.null(x[[i]])) next
            else {
                flag <- TRUE
                cat("Chr ", names(x)[i], ":\n", sep="")
                print(x[[i]], digits=3)
                cat("\n")
            }
        }
    }
    else if(tab=="tabByCol") {
        for(i in seq(along=x)) {
            if(nrow(x[[i]])==0) next
            else {
                flag <- TRUE
                if(length(x) > 1) cat(names(x)[i], ":\n", sep="")
                print(x[[i]], digits=3)
                if(length(x) > 1) cat("\n")
            }
        }
    }
    if(!flag)
        cat("    There were no LOD peaks above the threshold.\n")
}

# pull out maximum LOD peak, genome-wide
max.scanone <-
    function(object, chr, lodcolumn=1, na.rm=TRUE, ...)
{
    if(!any(class(object) == "scanone"))
        stop("Input must have class \"scanone\".")

    if(lodcolumn < 1 || lodcolumn+2 > ncol(object))
        stop("Argument lodcolumn misspecified.")

    if(!missing(chr)) object <- subset(object, chr=chr)

    maxlod <- max(object[,lodcolumn+2],na.rm=TRUE)
    wh <- which(!is.na(object[,lodcolumn+2]) & object[,lodcolumn+2]==maxlod)
    if(length(wh) > 1) wh <- sample(wh, 1)
    object <- object[wh,]

    object[,1] <- factor(as.character(unique(object[,1])))

    summary.scanone(object,threshold=0,lodcolumn=lodcolumn)
}

######################################################################
# subset.scanone
######################################################################
subset.scanone <-
    function(x, chr, lodcolumn, ...)
{
    if(!any(class(x) == "scanone"))
        stop("Input should have class \"scanone\".")

    if(missing(chr) && missing(lodcolumn))
        stop("You must specify either chr or lodcolumn.")

    y <- x

    if(!missing(chr)) {
        chr <- matchchr(chr, unique(x[,1]))
        x <- x[!is.na(match(x[,1],chr)), ,drop=FALSE]
        thechr <- as.character(x[,1])
        x[,1] <- factor(thechr, levels=unique(thechr))
    }

    if(!missing(lodcolumn)) {
        if(any(lodcolumn>0) && any(lodcolumn<0))
            stop("lodcolumn values can't be both >0 and <0.")
        if(any(lodcolumn<0) || is.logical(lodcolumn))
            lodcolumn <- (1:(ncol(x)-2))[lodcolumn]
        if(length(lodcolumn)==0)
            stop("You must retain at least one LOD column.")
        if(any(lodcolumn < 1 || lodcolumn > ncol(x)-2))
            stop("lodcolumn values must be >=1 and <=",ncol(x)-2)
        x <- x[,c(1,2,lodcolumn+2)]
    }

    class(x) <- class(y)
    nam <- names(attributes(y))
    if("method" %in% nam)
        attr(x, "method") <- attr(y,"method")
    if("type" %in% nam)
        attr(x, "type") <- attr(y,"type")
    if("model" %in% nam)
        attr(x, "model") <- attr(y,"model")

    x
}

######################################################################
# c.scanone
#
# Combine the results of multiple runs of scanone into single object
# (pasting the columns together).
######################################################################
c.scanone <-
    function(..., labels)
{
    dots <- list(...)
    cl1 <- class(dots[[1]])
    if(length(dots)==1 && length(cl1)==1 && cl1=="list") dots <- dots[[1]]

    if(length(dots)==1) return(dots[[1]])
    for(i in seq(along=dots)) {
        if(!any(class(dots[[i]]) == "scanone"))
            stop("Input should have class \"scanone\".")
    }

    if(!missing(labels)) {
        if(length(labels)==1)
            labels <- rep(labels, length(dots))
        if(length(labels) != length(dots))
            stop("labels needs to be the same length as the number of objects input.")
        gavelabels <- TRUE
    }
    else {
        labels <- grab.arg.names(...)
        gavelabels <- FALSE
    }

    nr <- sapply(dots, nrow)
    if(length(unique(nr)) != 1)
        stop("The input must all have the same number of rows.")

    chr <- lapply(dots, function(a) a$chr)
    pos <- lapply(dots, function(a) a$pos)
    for(i in 2:length(dots)) {
        if(any(chr[[1]] != chr[[i]]) || any(pos[[1]] != pos[[i]])) {
            cat("The input must conform exactly (same chr and positions\n")
            stop("(That is, calc.genoprob and/or sim.geno must have used the same step and off.end)\n")
        }
    }
    cl <- class(dots[[1]])

    thenam <- unlist(lapply(dots, function(a) colnames(a)[-(1:2)]))
    if(length(unique(thenam)) == length(thenam))
        repeats <- FALSE
    else repeats <- TRUE

    if(repeats || gavelabels) {
        for(i in 1:length(dots)) {
            colnames(dots[[i]])[-(1:2)] <- paste(colnames(dots[[i]])[-(1:2)], labels[i], sep=".")
            dots[[i]] <- as.data.frame(dots[[i]], stringsAsFactors=TRUE)
        }
    }

    result <- dots[[1]]

    for(i in 2:length(dots))
        result <- cbind(as.data.frame(result, stringsAsFactors=TRUE),
                        as.data.frame(dots[[i]][,-(1:2),drop=FALSE], stringsAsFactors=TRUE))

    class(result) <- cl
    result
}

cbind.scanone <- c.scanone

grab.arg.names <-
    function(...)
{
    # pull out the names from the input
    temp <- deparse(substitute(c(...)))
    temp <- unlist(strsplit(temp, ","))
    for(i in seq(along=temp))
        temp[i] <-  paste(unlist(strsplit(temp[i]," ")),collapse="")
    temp[1] <- substr(temp[1], 3, nchar(temp[1]))
    temp[length(temp)] <- substr(temp[length(temp)], 1, nchar(temp[length(temp)])-1)

    temp
}


######################################################################
# summary.scanoneperm
#
# Give genome-wide LOD thresholds on the basis of the results of
# scanone permutation test (from scanone with n.perm > 0)
######################################################################
summary.scanoneperm <-
    function(object, alpha=c(0.05, 0.10), controlAcrossCol=FALSE, ...)
{
    if(!any(class(object) == "scanoneperm"))
        stop("Input should have class \"scanoneperm\".")

    if(any(alpha < 0 | alpha > 1))
        stop("alpha should be between 0 and 1.")

    if("xchr" %in% names(attributes(object))) { # X-chromosome-specific results
        L <- attr(object, "L")
        thealpha <- cbind(1 - (1-alpha)^(L[1]/sum(L)), 1 - (1-alpha)^(L[2]/sum(L)))
        v <- c("A","X")

        quant <- vector("list", 2)
        names(quant) <- c("A","X")
        for(k in 1:2) {
            if(!is.matrix(object[[v[k]]])) object[[v[k]]] <- as.matrix(object[[v[k]]])

            if(controlAcrossCol) {
                if(any(is.na(object[[v[k]]])))
                    object[[v[k]]] <- object[[v[k]]][apply(object[[v[k]]],1,function(a) !any(is.na(a))),,drop=FALSE]

                r <- apply(object[[v[k]]], 2, rank, ties.method="random", na.last=FALSE)
                print(is.matrix(r))
                rmax <- apply(r, 1, max)
                rqu <- quantile(rmax, 1-thealpha[,k], na.rm=TRUE)
                qu <- matrix(nrow=length(thealpha[,k]), ncol=ncol(object[[v[k]]]))
                object.sort <- apply(object[[v[k]]], 2, sort, na.last=FALSE)

                for(i in seq(along=rqu)) {
                    if(fl==ce) # exact
                        qu[i,] <- object.sort[rqu[i],]
                    else # need to interpolate
                        qu[i,] <- object.sort[fl,]*(1-(ce-fl)) + object.sort[ce,]*(ce-fl)
                }
                colnames(qu) <- colnames(object[[v[k]]])
            }
            else
                qu <- apply(object[[v[k]]], 2, quantile, 1-thealpha[,k], na.rm=TRUE)

            if(!is.matrix(qu)) {
                nam <- names(qu)
                qu <- matrix(qu, nrow=length(alpha))
                dimnames(qu) <- list(paste(100*alpha,"%", sep=""), nam)
            }
            else rownames(qu) <- paste(100*alpha, "%", sep="")

            quant[[k]] <- qu
        }

        attr(quant, "n.perm") <- c("A"=nrow(object$A), "X"=nrow(object$X))
        class(quant) <- "summary.scanoneperm"
    }
    else {
        if(!is.matrix(object)) object <- as.matrix(object)

        if(controlAcrossCol) {
            if(any(is.na(object)))
                object <- object[apply(object,1,function(a) !any(is.na(a))),,drop=FALSE]

            r <- apply(object, 2, rank, ties.method="random", na.last=FALSE)
            rmax <- apply(r, 1, max)
            rqu <- quantile(rmax, 1-alpha, na.rm=TRUE)
            quant <- matrix(nrow=length(alpha), ncol=ncol(object))
            object.sort <- apply(object, 2, sort, na.last=FALSE)
            for(i in seq(along=rqu)) {
                fl <- floor(rqu[i])
                ce <- ceiling(rqu[i])
                if(fl==ce) # exact
                    quant[i,] <- object.sort[rqu[i],]
                else # need to interpolate
                    quant[i,] <- object.sort[fl,]*(1-(ce-fl)) + object.sort[ce,]*(ce-fl)
            }
            colnames(quant) <- colnames(object)
        }
        else
            quant <- apply(object, 2, quantile, 1-alpha, na.rm=TRUE)

        if(!is.matrix(quant)) {
            nam <- names(quant)
            quant <- matrix(quant, nrow=length(alpha))
            dimnames(quant) <- list(paste(100*alpha,"%", sep=""), nam)
        }
        else rownames(quant) <- paste(100*alpha, "%", sep="")

        attr(quant, "n.perm") <- nrow(object)
        class(quant) <- "summary.scanoneperm"
    }
    quant
}

print.summary.scanoneperm <-
    function(x, ...)
{
    n.perm <- attr(x, "n.perm")
    if(length(n.perm)==2) {
        cat("Autosome LOD thresholds (", n.perm[1], " permutations)\n", sep="")
        x$A <- x$A[1:nrow(x$A),,drop=FALSE]
        print(x$A, digits=3)

        cat("\nX chromosome LOD thresholds (", n.perm[2], " permutations)\n", sep="")
        x$X <- x$X[1:nrow(x$X),,drop=FALSE]
        print(x$X, digits=3)
    }
    else {
        cat("LOD thresholds (", n.perm, " permutations)\n", sep="")
        x <- x[1:nrow(x),,drop=FALSE]
        print(x, digits=3)
    }
}

######################################################################
# combine scanoneperm results ... paste the rows together
######################################################################
rbind.scanoneperm <- c.scanoneperm <-
    function(...)
{
    dots <- list(...)

    cl1 <- class(dots[[1]])
    if(length(dots)==1 && length(cl1)==1 && cl1=="list") dots <- dots[[1]]

    if(length(dots)==1) return(dots[[1]])
    for(i in seq(along=dots)) {
        if(!any(class(dots[[i]]) == "scanoneperm"))
            stop("Input should have class \"scanoneperm\".")
    }

    if("xchr" %in% names(attributes(dots[[1]]))) {
        xchr <- lapply(dots, attr, "xchr")
        L <- lapply(dots, attr, "L")
        for(i in 2:length(dots)) {
            if(length(xchr[[1]]) != length(xchr[[i]]) ||
               any(xchr[[1]] != xchr[[i]]))
                stop("xchr attributes in the input must be consistent.")
            if(length(L[[1]]) != length(L[[i]]) ||
               any(L[[1]] != L[[i]]))
                stop("L attributes in the input must be consistent.")
        }
        for(i in 1:length(dots)) dots[[i]] <- unclass(dots[[i]])

        ncA <- sapply(dots, function(a) ncol(a$A))
        ncX <- sapply(dots, function(a) ncol(a$X))
        if(length(unique(ncA)) != 1 || length(unique(ncX)) != 1)
            stop("The input must all have the same number of columns.")
        result <- dots[[1]]
        for(i in 2:length(dots)) {
            result$A <- rbind(result$A, dots[[i]]$A)
            result$X <- rbind(result$X, dots[[i]]$X)
        }
        class(result) <- "scanoneperm"
        attr(result, "xchr") <- xchr[[1]]
        attr(result, "L") <- L[[1]]
    }
    else {
        nc <- sapply(dots, ncol)
        if(length(unique(nc)) != 1)
            stop("The input must all have the same number of columns.")
        for(i in 1:length(dots)) dots[[i]] <- unclass(dots[[i]])
        result <- dots[[1]]
        for(i in 2:length(dots))
            result <- rbind(result, dots[[i]])
        class(result) <- "scanoneperm"
    }
    result
}

######################################################################
# combine scanoneperm results ... paste the columns together
######################################################################
cbind.scanoneperm <-
    function(..., labels)
{
    dots <- list(...)
    if(length(dots)==1) return(dots[[1]])

    for(i in seq(along=dots)) {
        if(!any(class(dots[[i]]) == "scanoneperm"))
            stop("Input should have class \"scanoneperm\".")
    }

    if(!missing(labels)) {
        if(length(labels)==1)
            labels <- rep(labels, length(dots))
        if(length(labels) != length(dots))
            stop("labels needs to be the same length as the number of objects input.")
        gavelabels <- TRUE
    }
    else {
        labels <- grab.arg.names(...)
        gavelabels <- FALSE
    }


    if("xchr" %in% names(attributes(dots[[1]]))) {
        xchr <- lapply(dots, attr, "xchr")
        L <- lapply(dots, attr, "L")
        for(i in 2:length(dots)) {
            if(length(xchr[[1]]) != length(xchr[[i]]) ||
               any(xchr[[1]] != xchr[[i]]))
                stop("xchr attributes in the input must be consistent.")
            if(length(L[[1]]) != length(L[[i]]) ||
               any(L[[1]] != L[[i]]))
                stop("L attributes in the input must be consistent.")
        }
        for(i in 1:length(dots)) dots[[i]] <- unclass(dots[[i]])

        nr <- sapply(dots, function(a) nrow(a$A))
        mnr <- max(nr)
        if(any(nr < mnr)) { # pad with NAs
            for(i in which(nr < mnr))
                dots[[i]]$A <- rbind(dots[[i]]$A, matrix(NA, ncol=ncol(dots[[i]]$A),
                                                         nrow=mnr-nr[i]))
        }

        nr <- sapply(dots, function(a) nrow(a$X))
        mnr <- max(nr)
        if(any(nr < mnr)) { # pad with NAs
            for(i in which(nr < mnr))
                dots[[i]]$X <- rbind(dots[[i]]$X, matrix(NA, ncol=ncol(dots[[i]]$X),
                                                         nrow=mnr-nr[i]))
        }

        thenamA <- unlist(lapply(dots, function(a) colnames(a$A)))
        thenamX <- unlist(lapply(dots, function(a) colnames(a$X)))
        if(length(unique(thenamA)) == length(thenamA) &&
           length(unique(thenamX)) == length(thenamX))
            repeats <- FALSE
        else repeats <- TRUE

        if(repeats || gavelabels) {
            colnames(dots[[1]]$A) <- paste(colnames(dots[[1]]$A),labels[1],sep=".")
            colnames(dots[[1]]$X) <- paste(colnames(dots[[1]]$X),labels[1],sep=".")
            for(i in 2:length(dots)) {
                colnames(dots[[i]]$A) <- paste(colnames(dots[[i]]$A),labels[i],sep=".")
                colnames(dots[[i]]$X) <- paste(colnames(dots[[i]]$X),labels[i],sep=".")
            }
        }

        result <- dots[[1]]
        for(i in 2:length(dots)) {
            result$A <- cbind(result$A, dots[[i]]$A)
            result$X <- cbind(result$X, dots[[i]]$X)
        }
        class(result) <- "scanoneperm"
        attr(result, "xchr") <- xchr[[1]]
        attr(result, "L") <- L[[1]]
    }
    else {
        for(i in 1:length(dots)) dots[[i]] <- unclass(dots[[i]])

        nr <- sapply(dots, nrow)
        mnr <- max(nr)
        if(any(nr < mnr)) { # pad with NAs
            for(i in which(nr < mnr))
                dots[[i]] <- rbind(dots[[i]], matrix(NA, ncol=ncol(dots[[i]]),
                                                     nrow=mnr-nr[i]))
        }

        thenam <- unlist(lapply(dots, colnames))
        if(length(unique(thenam)) == length(thenam))
            repeats <- FALSE
        else repeats <- TRUE

        if(repeats || gavelabels) {
            colnames(dots[[1]]) <- paste(colnames(dots[[1]]),labels[1],sep=".")
            for(i in 2:length(dots))
                colnames(dots[[i]]) <- paste(colnames(dots[[i]]), labels[i], sep=".")
        }
        result <- dots[[1]]
        for(i in 2:length(dots))
            result <- cbind(result, dots[[i]])

        class(result) <- "scanoneperm"
    }
    result
}


##############################
# subset.scanoneperm: pull out a set of lodcolumns
##############################
subset.scanoneperm <-
    function(x, repl, lodcolumn, ...)
{
    att <- attributes(x)

    if(is.list(x)) {
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
    }
    else {
        if(!is.matrix(x)) x <- as.matrix(x)

        if(missing(lodcolumn)) lodcolumn <- 1:ncol(x)
        else if(!is.null(attr(try(x[,lodcolumn], silent=TRUE),"try-error")))
            stop("lodcolumn misspecified.")

        if(missing(repl)) repl <- 1:nrow(x)
        else if(!is.null(attr(try(x[repl,], silent=TRUE),"try-error")))
            stop("repl misspecified.")

        cl <- class(x)
        x <- unclass(x)[repl,lodcolumn,drop=FALSE]
        class(x) <- cl
    }

    for(i in seq(along=att)) {
        if(names(att)[i] == "dim" || length(grep("names", names(att)[i]))>0) next
        attr(x, names(att)[i]) <- att[[i]]
    }

    x
}

# subset.scanoneperm using [,]
`[.scanoneperm` <-
    function(x, repl, lodcolumn)
    subset.scanoneperm(x, repl, lodcolumn)

# end of summary.scanone.R
