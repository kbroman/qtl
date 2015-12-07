######################################################################
#
# refineqtl.R
#
# copyright (c) 2006-2015, Karl W. Broman
# last modified Oct, 2015
# first written Jun, 2006
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
# Contains: refineqtl, plotLodProfile
#
######################################################################

######################################################################
# this is like scanqtl, though the positions are all fixed loci;
# it calls scanqtl iteratively, trying to find the best positions for
# each QTL.
#
# the method is like that of that described by Zeng and colleagues
# regarding MIM; each QTL is slid from one end of the chromosome to
# the other, or to the next flanking QTLs, if there are linked ones
#
# maxit is the maximum number of iterations (going through all QTLs
# in each iteration) to perform
######################################################################
refineqtl <-
    function(cross, pheno.col=1, qtl, chr, pos, qtl.name, covar=NULL, formula,
             method=c("imp", "hk"), model=c("normal", "binary"), verbose=TRUE, maxit=10,
             incl.markers=TRUE, keeplodprofile=TRUE, tol=1e-4, maxit.fitqtl=1000,
             forceXcovar=FALSE)
{
    method <- match.arg(method)
    model <- match.arg(model)

    if( !("cross" %in% class(cross)) )
        stop("The cross argument must be an object of class \"cross\".")

    # allow formula to be a character string
    if(!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

    if(!is.null(covar) && !is.data.frame(covar)) {
        if(is.matrix(covar) && is.numeric(covar))
            covar <- as.data.frame(covar, stringsAsFactors=TRUE)
        else stop("covar should be a data.frame")
    }

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    if(!missing(qtl) && (!missing(chr) || !missing(pos) || !missing(qtl.name)))
        warning("qtl argument is provided, and so chr, pos and qtl.name are ignored.")

    if(missing(qtl) && (missing(chr) || missing(pos)))
        stop("Provide either qtl or both chr and pos.")

    if(!missing(qtl)) {
        chr <- qtl$chr
        pos <- qtl$pos
    }
    else { # chr and pos provided
        if(missing(qtl.name)) {
            if(method=="imp")
                qtl <- makeqtl(cross, chr=chr, pos=pos, what="draws")
            else
                qtl <- makeqtl(cross, chr=chr, pos=pos, what="prob")
        }
        else {
            if(method=="imp")
                qtl <- makeqtl(cross, chr=chr, pos=pos, qtl.name=qtl.name, what="draws")
            else
                qtl <- makeqtl(cross, chr=chr, pos=pos, qtl.name=qtl.name, what="prob")
        }
    }

    if(method=="imp") {
        if(!("geno" %in% names(qtl))) {
            if("prob" %in% names(qtl)) {
                warning("The qtl object doesn't contain imputations; using method=\"hk\".")
                method <- "hk"
            }
            else
                stop("The qtl object needs to be created with makeqtl with what=\"draws\".")
        }
    }
    else {
        if(!("prob" %in% names(qtl))) {
            if("geno" %in% names(qtl)) {
                warning("The qtl object doesn't contain QTL genotype probabilities; using method=\"imp\".")
                method <- "imp"
            }
            else
                stop("The qtl object needs to be created with makeqtl with what=\"prob\".")
        }
    }


    if(!all(chr %in% names(cross$geno)))
        stop("Chr ", paste(unique(chr[!(chr %in% cross$geno)]), sep=" "),
             " not found in cross.")

    if(verbose > 1) scanqtl.verbose <- TRUE
    else scanqtl.verbose <- FALSE

    cross <- subset(cross, chr=as.character(unique(chr))) # pull out just those chromosomes

    # save map in the object
    map <- attr(qtl, "map")

    if(qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if(method=="imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        else
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
        attr(qtl, "map") <- map
    }
    if(method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3])  {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        attr(qtl, "map") <- map
    }

    # minimum distance between pseudomarkers
    if(is.null(map))
        stop("Input qtl object should contain the genetic map.")
    mind <- min(sapply(map, function(a) { if(is.matrix(a)) a <- a[1,]; min(diff(a)) }))/2
    if(mind <= 0) mind <- 1e-6

    # check phenotypes and covariates; drop ind'ls with missing values
    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("refineqtl can take just one phenotype; only the first will be used")
    }

    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }

    if(pheno.col < 1 || pheno.col > nphe(cross))
        stop("pheno.col should be between 1 and ", nphe(cross))
    pheno <- cross$pheno[,pheno.col]
    if(!is.null(covar) && nrow(covar) != length(pheno))
        stop("nrow(covar) != no. individuals in cross.")
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors=TRUE)
    hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
    if(all(hasmissing))
        stop("All individuals are missing phenotypes or covariates.")

    if(any(hasmissing)) {
        origcross <- cross
        origqtl <- qtl

        cross <- subset(cross, ind=!hasmissing)
        pheno <- pheno[!hasmissing]
        if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]

        if(method=="imp")
            qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
        else
            qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,,drop=FALSE])

        qtl$n.ind <- sum(!hasmissing)
    }

    # if missing formula, include the additive QTL plus all covariates
    if(missing(formula)) {
        formula <- paste("y ~", paste(qtl$altname, collapse="+"))
        if(!is.null(covar))
            formula <- paste(formula, "+", paste(colnames(covar), collapse="+"))
        formula <- as.formula(formula)
    }

    # drop covariates that are not in the formula
    if(!is.null(covar)) {
        theterms <- rownames(attr(terms(formula), "factors"))
        m <- match(colnames(covar), theterms)
        if(all(is.na(m))) covar <- NULL
        else covar <- covar[,!is.na(m),drop=FALSE]
    }

    formula <- checkformula(formula, qtl$altname, colnames(covar))

    # identify which QTL are in the model formula
    tovary <- sort(parseformula(formula, qtl$altname, colnames(covar))$idx.qtl)
    if(length(tovary) != qtl$n.qtl)
        reducedqtl <- dropfromqtl(qtl, index=(1:qtl$n.qtl)[-tovary])
    else reducedqtl <- qtl

    # if a QTL is missing from the formula, we need to revise the formula, moving
    # everything over, for use in scanqtl
    if(any(1:length(tovary) != tovary)) {
        tempform <- strsplit(deparseQTLformula(formula), " *~ *")[[1]][2]
        terms <- strsplit(tempform, " *\\+ *")[[1]]
        for(j in seq(along=terms)) {
            if(length(grep(":", terms[j])) > 0) { # interaction
                temp <- strsplit(terms[j], " *: *")[[1]]

                for(k in seq(along=temp)) {
                    g <- grep("^[Qq][0-9]+$", temp[k])
                    if(length(g) > 0) {
                        num <- as.numeric(substr(temp[k], 2, nchar(temp[k])))
                        temp[k] <- paste("Q", which(tovary == num), sep="")
                    }
                }
                terms[j] <- paste(temp, collapse=":")
            }
            else {
                g <- grep("^[Qq][0-9]+$", terms[j])
                if(length(g) > 0) {
                    num <- as.numeric(substr(terms[j], 2, nchar(terms[j])))
                    terms[j] <- paste("Q", which(tovary == num), sep="")
                }
            }
        }
        formula <- as.formula(paste("y ~", paste(terms, collapse=" + ")))
    }

    curpos <- pos[tovary]
    chrnam <- chr[tovary]
    if(verbose) cat("pos:", curpos, "\n")
    converged <- FALSE

    oldo <- NULL
    lc <- length(chrnam)

    lastout <- vector("list", length(curpos))
    names(lastout) <- qtl$name[tovary]

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

    for(i in 1:maxit) {
        if(keeplodprofile) # do drop-one analysis
            basefit <- fitqtlengine(pheno=pheno, qtl=reducedqtl, covar=covar, formula=formula,
                                    method=method, model=model, dropone=TRUE, get.ests=FALSE,
                                    run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                                    tol=tol, maxit=maxit.fitqtl, forceXcovar=forceXcovar)
        else
            basefit <- fitqtlengine(pheno=pheno, qtl=reducedqtl, covar=covar, formula=formula,
                                    method=method, model=model, dropone=FALSE, get.ests=FALSE,
                                    run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                                    tol=tol, maxit=maxit.fitqtl, forceXcovar=forceXcovar)

        if(i==1) {
            origlod <- curlod <- thisitlod <- basefit$result.full[1,4]
            origpos <- curpos
        }

        if(verbose) cat("Iteration", i, "\n")
        o <- sample(lc)

        # make sure the first here was not the last before
        if(!is.null(oldo))
            while(o[1] != oldo[lc])
                o <- sample(lc)
        oldo <- o

        newpos <- curpos
        for(j in o) {

            otherchr <- chrnam[-j]
            otherpos <- newpos[-j]

            thispos <- as.list(newpos)
            if(any(otherchr == chrnam[j])) { # linked QTLs
                linkedpos <- otherpos[otherchr==chr[j]]
                if(any(linkedpos < newpos[j]))
                    low <- max(linkedpos[linkedpos < newpos[j]])
                else low <- -Inf
                if(any(linkedpos > newpos[j]))
                    high <- min(linkedpos[linkedpos > newpos[j]])
                else high <- Inf

                thispos[[j]] <- c(low, high)
            }
            else
                thispos[[j]] <- c(-Inf, Inf)

            out <- scanqtl(cross=cross, pheno.col=pheno.col, chr=chrnam, pos=thispos,
                           covar=covar, formula=formula, method=method, model=model,
                           incl.markers=incl.markers,
                           verbose=scanqtl.verbose, tol=tol, maxit=maxit.fitqtl,
                           forceXcovar=forceXcovar)

            lastout[[j]] <- out

            newpos[j] <- as.numeric(strsplit(names(out)[out==max(out)],"@")[[1]][2])

            if(verbose) {
                cat(" Q", j, " pos: ", curpos[j], " -> ", newpos[j], "\n", sep="")
                cat("    LOD increase: ", round(max(out) - curlod, 3), "\n")
            }
            curlod <- max(out)
        }

        if(verbose) {
            cat("all pos:", curpos, "->", newpos, "\n")
            cat("LOD increase at this iteration: ", round(curlod - thisitlod, 3), "\n")
        }
        thisitlod <- curlod

        if(max(abs(curpos - newpos)) < mind) {
            converged <- TRUE
            break
        }
        curpos <- newpos

        reducedqtl <- replaceqtl(cross, reducedqtl, seq(length(curpos)),
                                 reducedqtl$chr, curpos, reducedqtl$name)
    }

    if(verbose) {
        cat("overall pos:", origpos, "->", newpos, "\n")
        cat("LOD increase overall: ", round(curlod - origlod, 3), "\n")
    }

    if(!converged) warning("Didn't converge.")

    # do the qtl have custom names?
    g <- grep("^.+@[0-9\\.]+$", qtl$name)
    if(length(g) == length(qtl$name)) thenames <- NULL
    else thenames <- qtl$name

    if(any(hasmissing)) {
        qtl <- origqtl
        cross <- origcross
    }
    for(j in seq(along=tovary))
        qtl <- replaceqtl(cross, qtl, tovary[j], chrnam[j], newpos[j])

    if(!is.null(thenames)) qtl$name <- thenames

    if(keeplodprofile) {
        # subtract off the results from the drop-one analysis from the LOD profiles
        dropresult <- basefit$result.drop
        if(is.null(dropresult)) {
            if(length(lastout)==1) {
                dropresult <- rbind(c(NA,NA, basefit$result.full[1,4]))
                rownames(dropresult) <- names(lastout)
            }
            else
                stop("There's a problem: need dropresult, but didn't obtain one.")
        }

        rn <- rownames(dropresult)
        qn <- names(lastout)

        for(i in seq(along=lastout)) {
            if(sum(rn==qn[i])>1) # ack! multiple QTL at same position
                warning("Multiple QTL at the same location.")
            lastout[[i]] <- lastout[[i]] - (max(lastout[[i]]) - max(dropresult[rn==qn[i],3]))
            pos <- as.numeric(matrix(unlist(strsplit(names(lastout[[i]]), "@")),byrow=TRUE,ncol=2)[,2])
            chr <- rep(qtl$chr[tovary][i], length(pos))
            lastout[[i]] <- data.frame(chr=chr, pos=pos, lod=as.numeric(lastout[[i]]), stringsAsFactors=TRUE)
        }
        names(lastout) <- qtl$name[tovary]

        # make the profiles scanone objects
        for(i in seq(along=lastout)) {
            class(lastout[[i]]) <- c("scanone", "data.frame")
            thechr <- qtl$chr[i]
            if(method=="imp")
                detailedmap <- attr(cross$geno[[thechr]]$draws,"map")
            else
                detailedmap <- attr(cross$geno[[thechr]]$prob,"map")

            if(is.matrix(detailedmap)) detailedmap <- detailedmap[1,]

            r <- range(lastout[[i]][,2])+c(-1e-5, 1e-5)
            rn <- names(detailedmap)[detailedmap>=r[1] & detailedmap<=r[2]]
            o <- grep("^loc-*[0-9]+",rn)
            if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
                rn[o] <- paste("c",thechr,".",rn[o],sep="")
            #      if(length(rn) != nrow(lastout[[i]])) return(list(lastout[[i]], rn, detailedmap))
            if(length(rn) == nrow(lastout[[i]])) rownames(lastout[[i]]) <- rn
        }

        attr(qtl, "lodprofile") <- lastout
    }

    # if there's a pLOD attribute, revise it
    if("pLOD" %in% names(attributes(qtl)) && curlod > origlod)
        attr(qtl,"pLOD") <- attr(qtl,"pLOD") + curlod - origlod

    qtl
}

######################################################################
# plotLodProfile
#
# This is for creating a plot of 1-d LOD profiles calculated within
# refineqtl.
######################################################################

plotLodProfile <-
    function(qtl, chr, incl.markers=TRUE, gap=25, lwd=2, lty=1, col="black",
             qtl.labels=TRUE, mtick=c("line", "triangle"),
             show.marker.names=FALSE, alternate.chrid=FALSE, add=FALSE,
             showallchr=FALSE, labelsep=5, ...)
{
    if(!("qtl" %in% class(qtl)))
        stop("Input qtl is not a qtl object")

    lodprof <- attr(qtl, "lodprofile")
    if(is.null(lodprof))
        stop("You must first run refineqtl, using keeplodprofile=TRUE")

    m <- match(qtl$name, names(lodprof))
    if(any(is.na(m)))
        qtl <- dropfromqtl(qtl, index=which(is.na(m)), drop.lod.profile=FALSE)

    # reorder qtl by position
    if(qtl$n.qtl > 1) {
        chrindex <- match(qtl$chr, names(attr(qtl, "map")))
        if(any(is.na(chrindex)))
            stop("Cannot find chr ", qtl$chr[is.na(chrindex)])
        neworder <- order(chrindex, qtl$pos)
        if(any(neworder != seq(qtl$n.qtl))) {
            qtl <- reorderqtl(qtl, neworder)
            lodprof <- attr(qtl, "lodprofile")
        }
    }

    n.qtl <- length(lodprof)
    if(length(lwd) == 1) lwd <- rep(lwd, n.qtl)
    else {
        if(length(lwd) != n.qtl) {
            warning("lwd should have length 1 or ", n.qtl)
            lwd <- rep(lwd[1], n.qtl)
        }
        else lwd <- lwd[neworder]
    }


    if(length(lty) == 1) lty <- rep(lty, n.qtl)
    else {
        if(length(lty) != n.qtl) {
            warning("lty should have length 1 or ", n.qtl)
            lty <- rep(lty[1], n.qtl)
        }
        else lty <- lty[neworder]
    }
    if(length(col) == 1) col <- rep(col, n.qtl)
    else {
        if(length(col) != n.qtl) {
            warning("col should have length 1 or ", n.qtl)
            if(length(col) < n.qtl)
                col <- rep(col, n.qtl)[1:n.qtl]
            else
                col <- col[1:n.qtl]
        }
        else col <- col[neworder]
    }

    map <- attr(qtl, "map")
    if(is.null(map))
        stop("Input qtl object should contain the genetic map.")

    thechr <- unique(qtl$chr)
    orderedchr <- names(map)
    if(showallchr) chr2keep <- seq(along=orderedchr)
    else chr2keep <- which(!is.na(match(orderedchr, thechr)))

    tempscan <- NULL
    for(i in chr2keep) {
        temp <- data.frame(chr=orderedchr[i],
                           pos=as.numeric(map[[i]]),
                           lod=NA, stringsAsFactors=TRUE)
        rownames(temp) <- names(map[[i]])
        tempscan <- rbind(tempscan, temp)
    }
    class(tempscan) <- c("scanone", "data.frame")

    if(missing(chr)) {
        if(showallchr) chr <- orderedchr
        else chr <- thechr
    }

    dontskip <- which(!is.na(match(qtl$chr, chr)))
    if(length(dontskip)==0)
        stop("Nothing to plot.")

    ymax <- max(sapply(lodprof[dontskip], function(a) max(a[,3], na.rm=TRUE)))
    begend <- matrix(unlist(tapply(tempscan[,2],tempscan[,1],range)),ncol=2,byrow=TRUE)
    rownames(begend) <- unique(tempscan[,1])
    begend <- begend[as.character(chr),,drop=FALSE]
    len <- begend[,2]-begend[,1]
    if(length(chr)==1) start <- 0
    else start <- c(0,cumsum(len+gap))-c(begend[,1],0)
    names(start) <- chr

    dots <- list(...)

    if(!add) {
        if("ylim" %in% names(dots)) {
            plot.scanone(tempscan, chr=chr, incl.markers=incl.markers, gap=gap,
                         mtick=mtick, show.marker.names=show.marker.names,
                         alternate.chrid=alternate.chrid, type="n", ...)
        }
        else {
            if(qtl.labels)
                ylim <- c(0, ymax+1)
            else
                ylim <- c(0, ymax)

            plot.scanone(tempscan, chr=chr, incl.markers=incl.markers, gap=gap,
                         mtick=mtick, show.marker.names=show.marker.names,
                         alternate.chrid=alternate.chrid, type="n", ylim=ylim,
                         ...)
        }
    }

    for(i in dontskip) {
        temp <- rbind(tempscan[tempscan[,1] != qtl$chr[i] |
                               (tempscan[,1] == qtl$chr[i] & (tempscan[,2] < min(lodprof[[i]][,2]) |
                                        tempscan[,2] > max(lodprof[[i]][,2]))),],
                      lodprof[[i]])

        temp <- temp[order(match(temp[,1], orderedchr), temp[,2]),]
        class(temp) <- c("scanone", "data.frame")

        plot.scanone(temp, chr=chr, incl.markers=FALSE, gap=gap, add=TRUE,
                     col=col[i], lwd=lwd[i], lty=lty[i], ...)

        if(qtl.labels) {
            maxlod <- max(temp[,3], na.rm=TRUE)
            maxpos <- median(temp[!is.na(temp[,3]) & temp[,3]==maxlod,2] + start[qtl$chr[i]])
            d <- diff(par("usr")[3:4]*labelsep/100)

            text(maxpos, maxlod + d, names(lodprof)[i], col=col[i], font=(lwd[i]>1)+1, ...)
        }
    }

    invisible()
}

# end of refineqtl.R
