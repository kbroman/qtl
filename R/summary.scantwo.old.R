######################################################################
#
# summaryScantwoOld.R
#
# copyright (c) 2001-2011, Karl W Broman, Hao Wu, and Brian Yandell
# last modified May, 2011
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
# Hao Wu (The Jackson Lab) wrote the initial code for summary.scantwo
# function.  Brian Yandell made further modifications/enhancements to
# summary.scantwo, but Karl re-wrote most of it later.
#
# Part of the R/qtl package
# Contains: summaryScantwoOld, print.summary.scantwo.old
#
######################################################################

summaryScantwoOld <-
    function (object, thresholds = c(0, 0, 0), lodcolumn=1,
              type = c("joint","interaction"), ...)
{
    warning("This function is provided solely for continuity of the software;\n",
            "it is not recommended.\n")
    if(!any(class(object) == "scantwo"))
        stop("Input should have class \"scantwo\".")

    type <- match.arg(type)

    if(length(dim(object$lod)) > 2) { # results from multiple phenotypes
        if(length(lodcolumn) > 1) {
            warning("Argument lodcolumn should be of length 1.")
            lodcolumn <- lodcolumn[1]
        }

        if(lodcolumn < 0 || lodcolumn > dim(object$lod)[3])
            stop("Argument lodcolumn misspecified.")
        object$lod <- object$lod[,,lodcolumn]
    }

    if(length(thresholds) < 3) {
        if(length(thresholds) == 1) thresholds <- c(thresholds, 0, 0)
        else stop("You must give three thresholds: full, interaction and main\n")
    }

    thrfull <- thresholds[1]
    thrint <- thresholds[2]
    thrcond <- thresholds[3]

    lod <- object$lod
    map <- object$map

    # backward compatibility for previous version of R/qtl
    if(is.na(match("scanoneX",names(object)))) {
        warning("It would be best to re-run scantwo() with the R/qtl version 0.98 or later.")
        scanoneX <- NULL
    }
    else scanoneX <- object$scanoneX

    # deal with bad LOD score values
    if(any(is.na(lod) | lod < -1e-06 | lod == Inf))
        warning("Some LOD scores NA, Inf or < 0; set to 0")
    lod[is.na(lod) | lod < 0 | lod == Inf] <- 0

    # if there's no mainscan result, ignore the thresholds
    #     and don't include the 4 conditional LOD columns
    if(all(is.na(diag(lod)) | diag(lod) < 1e-10))
        includes.scanone <- FALSE
    else includes.scanone <- TRUE

    # change lod scores to old version
    u <- upper.tri(lod)
    lod[u] <- t(lod)[u] - lod[u]

    # If scanone results available, calculate conditional LOD scores
    if(includes.scanone) {
        d <- diag(lod)
        q1 <- matrix(rep(d,length(d)),ncol=length(d))
        q2 <- matrix(rep(d,length(d)),ncol=length(d),byrow=TRUE)

        if(!is.null(scanoneX) && any(map[,4])) {
            d <- scanoneX
            q1X <- matrix(rep(d,length(d)),ncol=length(d))
            q2X <- matrix(rep(d,length(d)),ncol=length(d),byrow=TRUE)
            q1[map[,4],] <- q1X[map[,4],]
            q2[,map[,4]] <- q2X[,map[,4]]
        }

        q1[lower.tri(q1)] <- t(q2)[lower.tri(q2)]
        condlod <- abs(lod - t(lod)) - q1
        diag(condlod) <- 0
    }
    else condlod <- NULL

    # Negative thresholds are interpreted relative to the maximum LOD score
    if(thrfull < 0)
        thrfull <- max(0,max(lod[lower.tri(lod)]) + thrfull)
    if(thrint < 0)
        thrint <- max(0,max(lod[upper.tri(lod)]) + thrint)
    if(thrcond < 0 && includes.scanone)
        thrcond <- max(0,max(condlod) + thrcond)

    crosstype <- attr(object, "type")
    if(is.null(crosstype)) {
        warning("No type attribute in input data; assuming backcross.")
        crosstype <- "bc"
    }

    # calculate the degree of freedom
    if(crosstype == "bc" || crosstype == "riself" || crosstype ==
       "risib" || crosstype=="dh") {
        df.int <- 1
        df.add <- 1
    }
    else if(crosstype == "f2") {
        df.int <- 4
        df.add <- 2
    }
    else if(crosstype == "4way") {
        df.int <- 9
        df.add <- 3
    }
    else {
        stop("Don't know what to do with cross type ", crosstype)
    }

    # chromsomes in the result
    chr <- unique(map[, 1])
    n.chr <- length(chr)

    # calculate the locations of each chromosome within the LOD matrix
    wh.index <- vector("list", n.chr)
    n <- nrow(map)
    for(i in 1:n.chr)
        wh.index[[i]] <- which(map[, 1] == chr[i])

    results <- NULL

    # go through each pair of chromosomes
    for(i in 1:n.chr) {
        for(j in i:n.chr) {
            tmplod1 <- lod[wh.index[[j]], wh.index[[i]],drop=FALSE]
            if(!is.null(condlod)) {
                if(i==j) tmpcondlod <- condlod[wh.index[[i]],wh.index[[i]],drop=FALSE]
                else {
                    tmpcondlod1 <- condlod[wh.index[[j]],wh.index[[i]],drop=FALSE]
                    tmpcondlod2 <- condlod[wh.index[[i]],wh.index[[j]],drop=FALSE]
                }
            }

            if(i != j) tmplod2 <- lod[wh.index[[i]], wh.index[[j]],drop=FALSE]
            else tmplod2 <- tmplod1


            if(type == "joint") {
                if(i == j) {
                    tri <- lower.tri(tmplod1)
                    lod.joint <- max(tmplod1[tri])
                    idx <- which(tmplod1 == lod.joint & tri, arr.ind=TRUE)
                    if(!is.matrix(idx)) {
                        cat("problem\n")
                        return(tmplod1)
                    }
                }
                else {
                    lod.joint <- max(tmplod1)
                    idx <- which(tmplod1 == lod.joint, arr.ind=TRUE)
                    if(!is.matrix(idx)) {
                        cat("problem\n")
                        return(tmplod1)
                    }
                }
                if(nrow(idx)>1) idx <- idx[sample(nrow(idx),1),,drop=FALSE]
                idx.row <- idx[1]
                idx.col <- idx[2]

                lod.int <- tmplod2[idx.col, idx.row,drop=FALSE]
            }
            else { # interaction lod
                if(i == j) {
                    tri <- upper.tri(tmplod2)
                    lod.int <- max(tmplod2[tri])
                    idx <- which(tmplod2 == lod.int & tri, arr.ind=TRUE)
                }
                else {
                    lod.int <- max(tmplod2)
                    idx <- which(tmplod2 == lod.int)
                }
                if(nrow(idx)>1) idx <- idx[sample(nrow(idx),1),,drop=FALSE]
                idx.row <- idx[2]
                idx.col <- idx[1]

                lod.joint <- tmplod1[idx.row, idx.col,drop=FALSE]
            }

            full.idx.row <- idx.row + wh.index[[j]][1] - 1
            full.idx.col <- idx.col + wh.index[[i]][1] - 1

            flag <- FALSE # a flag to indicate whether there's any peak on this pair
            if(lod.joint >= thrfull) {
                if(includes.scanone) {
                    if(i==j) {
                        lod.q1 <- tmpcondlod[idx.row,idx.col,drop=FALSE]
                        lod.q2 <- tmpcondlod[idx.col,idx.row,drop=FALSE]
                    }
                    else {
                        lod.q1 <- tmpcondlod1[idx.row,idx.col,drop=FALSE]
                        lod.q2 <- tmpcondlod2[idx.col,idx.row,drop=FALSE]
                    }

                    if(lod.int >= thrint || min(c(lod.q1, lod.q2)) >= thrcond) {
                        flag <- TRUE
                        i.pos <- map[full.idx.col, 2]
                        j.pos <- map[full.idx.row, 2]
                        results <- rbind(results,
                                         data.frame(chr[i], chr[j], i.pos, j.pos,
                                                    lod.joint, 1 - pchisq(2 * log(10) * lod.joint,
                                                                          df.int + 2 * df.add),
                                                    lod.int, 1 - pchisq(2 * log(10) * lod.int, df.int),
                                                    lod.q1, 1 - pchisq(2 * log(10) * lod.q1, df.add),
                                                    lod.q2, 1 - pchisq(2 * log(10) * lod.q2, df.add), stringsAsFactors=TRUE)
                                         )
                    }
                }
                else { # no scanone output
                    flag <- TRUE
                    i.pos <- map[full.idx.col, 2]
                    j.pos <- map[full.idx.row, 2]
                    results <- rbind(results,
                                     data.frame(chr[i], chr[j], i.pos, j.pos,
                                                lod.joint, 1 - pchisq(2 * log(10) * lod.joint,
                                                                      df.int + 2 * df.add),
                                                lod.int, 1 - pchisq(2 * log(10) * lod.int, df.int), stringsAsFactors=TRUE)
                                     )
                }
                # give the new row (if any) a name
                if(flag) {
                    mname <- rownames(map)
                    rownames(results)[nrow(results)] <- paste(mname[full.idx.col], ":",
                                                              mname[full.idx.row], sep="")
                }
            } # lod joint above threshold

        } # end loop over chromosomes
    }

    if(is.null(results)) {
        results <- numeric(0)
    }
    else {
        if(includes.scanone)
            colnames(results) <- c("chr1", "chr2", "pos1", "pos2",
                                   "lod.joint", "p.joint", "lod.int", "p.int", "lod.q1",
                                   "p.q1", "lod.q2", "p.q2")
        else colnames(results) <- c("chr1", "chr2", "pos1", "pos2",
                                    "lod.joint", "p.joint", "lod.int", "p.int")
        results <- as.data.frame(results, stringsAsFactors=TRUE)
    }
    class(results) <- c("summary.scantwo.old", "data.frame")
    results
}


print.summary.scantwo.old <-
    function(x,...)
{
    if(length(x)==0) {
        cat("    There were no pairs of loci meeting the criteria.\n")
        return(invisible(NULL))
    }

    # column names
    cnames <- c("pos1", "pos2", "  LODjnt", "-logP",
                "  LODint", "-logP", "  LODq1", "-logP",
                "  LODq2", "-logP")

    # chr names
    chr1 <- paste("c",x[,1],sep="")
    chr2 <- paste("c",x[,2],sep="")

    # pad chr names with spaces; this isn't really necessary
    nchar.c1 <- nchar(chr1); max.nchar.c1 <- max(nchar.c1)
    nchar.c2 <- nchar(chr2); max.nchar.c2 <- max(nchar.c2)
    if(any(nchar.c1 < max.nchar.c1 | nchar.c2 < max.nchar.c2)) {
        for(i in 1:length(nchar.c2)) {
            if(nchar.c1[i] < max.nchar.c1)
                chr1[i] <- paste(paste(rep(" ", max.nchar.c1-nchar.c1[i]),collapse=""),
                                 chr1[i],sep="")
            if(nchar.c2[i] < max.nchar.c2)
                chr2[i] <- paste(paste(rep(" ", max.nchar.c2-nchar.c2[i]),collapse=""),
                                 chr2[i],sep="")
        }
    }
    chr <- paste(chr1,chr2,sep=":")

    # round the rest; take -log10(P-values)
    for(j in 3:ncol(x)) {
        if(j<5)
            x[,j] <- round(x[,j])
        else if(j %% 2)  # odd
            x[,j] <- round(x[,j],2)
        else
            x[,j] <- -round(log10(x[,j]),1)
    }

    res <- as.data.frame(x[,-(1:2)], stringsAsFactors=TRUE)
    names(res) <- cnames[1:ncol(res)]
    rownames(res) <- chr

    cat("\n")
    print.data.frame(res)
    cat("\n")
}

# end of summary.scantwo.old.R
