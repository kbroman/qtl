######################################################################
#
# scanqtl.R
#
# copyright (c) 2002-2014, Hao Wu and Karl W. Broman
# last modified Jan, 2014
# first written Apr, 2002
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
# Contains: scanqtl
#
######################################################################

scanqtl <-
    function(cross, pheno.col=1, chr, pos, covar=NULL, formula,
             method=c("imp", "hk"), model=c("normal", "binary"),
             incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
             forceXcovar=FALSE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(!is.null(covar) && !is.data.frame(covar)) {
        if(is.matrix(covar) && is.numeric(covar))
            covar <- as.data.frame(covar, stringsAsFactors=TRUE)
        else stop("covar should be a data.frame")
    }

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("scanqtl can take just one phenotype; only the first will be used")
    }

    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }

    if(pheno.col < 1 | pheno.col > nphe(cross))
        stop("pheno.col values should be between 1 and the no. phenotypes")

    pheno <- cross$pheno[,pheno.col]
    if(!is.null(covar) && nrow(covar) != length(pheno))
        stop("nrow(covar) != no. individuals in cross.")

    method <- match.arg(method)
    model <- match.arg(model)

    # allow formula to be a character string
    if(!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

    if(method=="imp") {
        if(!("draws" %in% names(cross$geno[[1]]))) {
            if("prob" %in% names(cross$geno[[1]])) {
                warning("The cross doesn't contain imputations; using method=\"hk\".")
                method <- "hk"
            }
            else
                stop("You need to first run sim.geno.")
        }
    }
    else {
        if(!("prob" %in% names(cross$geno[[1]]))) {
            if("draws" %in% names(cross$geno[[1]])) {
                warning("The cross doesn't contain QTL genotype probabilities; using method=\"imp\".")
                method <- "imp"
            }
            else
                stop("You need to first run calc.genoprob.")
        }
    }

    if(method=="imp") {
        if("stepwidth" %in% names(attributes(cross$geno[[1]]$draws)) &&
           attr(cross$geno[[1]]$draws, "stepwidth") != "fixed") {
            stepwidth.var <- TRUE
            incl.markers <- TRUE
        }
        else stepwidth.var <- FALSE
    }
    else {
        if("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
           attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
            stepwidth.var <- TRUE
            incl.markers <- TRUE
        }
        else stepwidth.var <- FALSE
    }

    type <- class(cross)[1]
    chrtype <- sapply(cross$geno,class)

    # input data checking
    if( length(chr) != length(pos))
        stop("Input chr and pos must have the same length")
    # note that input chr is a vector and pos is a list

    method <- match.arg(method)

    ichr <- match(chr, names(cross$geno))
    if(any(is.na(ichr)))
        stop("There's no chromosome number ", chr[is.na(ichr)],
             " in input cross object")


    # if formula is missing, make one.
    # All QTLs and covariates will be additive by default
    n.qtl <- length(chr)
    n.covar <- length(covar)
    if(missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
        formula <- "y~Q1"
        if(n.qtl > 1)
            for (i in 2:n.qtl)
                formula <- paste(formula, tmp.Q[i], sep="+")
        if (n.covar) { # if covariate is not empty
            tmp.C <- names(covar) # covariate term names
            for(i in 1:n.covar)
                formula <- paste(formula, tmp.C[i], sep="+")
        }
        formula <- as.formula(formula)
    }
    else {
        # include all input QTLs and covariates in the formula additively
        formula.str <- deparseQTLformula(formula) # deparse formula as a string
        for(i in 1:n.qtl) { # loop thru the QTLs
            qtl.term <- paste("Q", i, sep="")
            if( length(grep(qtl.term, formula.str, ignore.case=TRUE))==0 )
                # this term is not in the formula
                # add it to the formula
                formula.str <- paste(formula.str, qtl.term, sep="+")
        }
        if(n.covar) { # covariate is not empty
            for(i in 1:n.covar) {
                covar.term <- names(covar)[i]
                if( length(grep(covar.term, formula.str, ignore.case=TRUE))==0 )
                    # this term is not in the formula
                    # add it to the formula
                    formula.str <- paste(formula.str, covar.term, sep="+")
            }
        }
        formula <- as.formula(formula.str)
    }

    # check the formula
    formula <- checkformula(formula, paste("Q", 1:length(chr), sep=""),
                            colnames(covar))

    # drop covariates that are not in the formula
    if(!is.null(covar)) {
        theterms <- rownames(attr(terms(formula), "factors"))
        m <- match(colnames(covar), theterms)
        if(all(is.na(m))) covar <- NULL
        else covar <- covar[,!is.na(m),drop=FALSE]
    }

    # check phenotypes and covariates; drop ind'ls with missing values
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- cbind(pheno)
    if(any(is.na(phcovar))) {
        if(ncol(phcovar)==1) hasmissing <- is.na(phcovar)
        else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
        if(all(hasmissing))
            stop("All individuals are missing phenotypes or covariates.")
        if(any(hasmissing)) {
            warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
            cross <- subset(cross, ind=!hasmissing)
            pheno <- pheno[!hasmissing]
            if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]
        }
    }
    sexpgm <- getsex(cross)

    # find the chromosome with multiple QTLs
    # indices for chromosomes with multiple QTLs
    idx.varied <- NULL
    indices <- pos  ## added by Karl 8/23/05
    for(i in 1:length(pos)) {
        l <- length(pos[[i]] )
        if( l >= 2 ) {
            # if there're more than two elements in pos, issue warning message
            if(l > 2) {
                msg <- "There are more than two elements in "
                msg <- paste(msg, i, "th input pos.")
                msg <- paste(msg, "The first two are taken as starting and ending position.")
                warning(msg)
            }

            # user specified a range
            # find all markers in this range
            idx.varied <- c(idx.varied, i)
            # make the genetic map on this chromosome

            # make genetic map
            if(method=="imp") {
                if("map" %in% names(attributes(cross$geno[[ichr[i]]]$draws)))
                    map <- attr(cross$geno[[ichr[i]]]$draws,"map")
                else {
                    stp <- attr(cross$geno[[ichr[i]]]$draws, "step")
                    oe <- attr(cross$geno[[ichr[i]]]$draws, "off.end")

                    if("stepwidth" %in% names(attributes(cross$geno[[ichr[i]]]$draws)))
                        stpw <- attr(cross$geno[[ichr[i]]]$draws, "stepwidth")
                    else stpw <- "fixed"
                    map <- create.map(cross$geno[[ichr[i]]]$map,stp,oe,stpw)
                }
            }
            else {
                if("map" %in% names(attributes(cross$geno[[ichr[i]]]$prob)))
                    map <- attr(cross$geno[[ichr[i]]]$prob,"map")
                else {
                    stp <- attr(cross$geno[[ichr[i]]]$prob, "step")
                    oe <- attr(cross$geno[[ichr[i]]]$prob, "off.end")

                    if("stepwidth" %in% names(attributes(cross$geno[[ichr[i]]]$prob)))
                        stpw <- attr(cross$geno[[ichr[i]]]$prob, "stepwidth")
                    else stpw <- "fixed"
                    map <- create.map(cross$geno[[ichr[i]]]$map,stp,oe,stpw)
                }

            }
            # pull out the female map if there are sex-specific maps
            if(is.matrix(map)) map <- map[1,]

            indices[[i]] <- seq(along=map)
            if(method=="imp")
                step <- attr(cross$geno[[ichr[i]]]$draws,"step")
            else
                step <- attr(cross$geno[[ichr[i]]]$prob,"step")

            if(!incl.markers && step>0) { # equally spaced positions
                eq.sp.pos <- seq(min(map), max(map), by=step)
                wh.eq.pos <- match(eq.sp.pos, map)
                map <- map[wh.eq.pos]
                indices[[i]] <- indices[[i]][wh.eq.pos]
            }

            # locate the markers given starting and ending postion
            # we should do this before or after incl.markers?
            start <- pos[[i]][1]
            end <- pos[[i]][2]
            # replace pos[[i]] (a range) by the marker positions within the range
            # extend the position to the nearest markers outside the ranges
            tmp <- which( (map - start)<=0 )
            if(length(tmp) != 0) # starting position is after the first marker
                start <- map[max(tmp)]
            tmp <- which( (end-map) <= 0 )
            if(length(tmp) != 0) # ending position is before the last marker
                end <- map[min(tmp)]
            pos[[i]] <- as.vector( map[(map>=start)&(map<=end)] )
            indices[[i]] <- indices[[i]][(map>=start)&(map<=end)]
        }
    }
    # Now, pos contains all the marker positions for all chromosomes

    #########################
    # Now start general scan
    #########################
    # There might be several chromosomes with multiple QTLs
    # Use one loop

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

    # number of chromosomes with multiple positions to be scanned
    n.idx.varied <- length(idx.varied)
    n.loop <- 1 # total number of loops
    if(n.idx.varied != 0) { # there IS some chromosomes with multiple QTL
        # vector to indicate the positions indices for those chromosomes
        idx.pos <- rep(0, n.idx.varied)
        l.varied <- NULL
        for(i in 1:n.idx.varied) {
            l.varied[i] <- length(pos[[idx.varied[i]]])
            n.loop <- n.loop * l.varied[i]
        }
        # initialize output variable
        result <- array(rep(0, n.loop), rev(l.varied))
        matrix.rank <- matrix.ncol <- array(rep(0, n.loop), rev(l.varied))
    }
    else { # fixed QTL model (no scanning)
        if(method=="imp")
            qtl <- makeqtl(cross, chr=chr, pos=unlist(pos), what="draws")
        else
            qtl <- makeqtl(cross, chr=chr, pos=unlist(pos), what="prob")

        result <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar,
                               formula=formula, method=method, model=model, dropone=FALSE,
                               get.ests=FALSE, run.checks=FALSE, cross.attr=cross.attr,
                               sexpgm=sexpgm, tol=tol, maxit=maxit, forceXcovar=forceXcovar)
        matrix.rank <- attr(result, "matrix.rank")
        matrix.ncol <- attr(result, "matrix.ncol")

        result <- result[[1]][1,4]
        names(result) <- "LOD"
        class(result) <- "scanqtl"
        attr(result, "method") <- method
        attr(result, "formula") <- deparseQTLformula(formula)
        attr(result, "matrix.rank") <- matrix.rank
        attr(result, "matrix.ncol") <- matrix.ncol
        return(result)
    }

    # loop thru all varied QTLs
    if(verbose) {
        cat(" ",n.loop, "models to fit\n")
        n.prnt <- floor(n.loop/20)
        if(n.prnt < 1) n.prnt <- 1
    }
    current.pos <- NULL ## added by Karl 8/23/05
    for(i in 1:n.loop) {
        # find the indices for positions
        remain <- i
        if(n.idx.varied > 1) {
            for(j in 1:(n.idx.varied-1)) {
                ns <- 1
                for( k in (j+1):n.idx.varied )
                    ns <- ns * length(pos[[idx.varied[k]]])
                idx.pos[j] <- floor(remain / ns) + 1
                remain <- remain - (idx.pos[j]-1) * ns
                # remain cannot be zero
                if(remain == 0) {
                    idx.pos[j] <- idx.pos[j] - 1
                    remain <- remain + ns
                }
            }
        }
        idx.pos[n.idx.varied] <- remain

        # make an QTL object
        pos.tmp <- NULL
        for(j in 1:length(pos)) {
            if(j %in% idx.varied) {
                idx.tmp <- which(j==idx.varied)
                pos.tmp <- c(pos.tmp, pos[[j]][idx.pos[idx.tmp]])
            }
            else
                pos.tmp <- c(pos.tmp, pos[[j]])
        }

        # this bit revised by Karl 8/23/05; now we make the qtl object
        #     once, and copy stuff over otherwise
        if(is.null(current.pos)) {
            if(method=="imp")
                qtl.obj <- makeqtl(cross, chr, pos.tmp, what="draws")
            else
                qtl.obj <- makeqtl(cross, chr, pos.tmp, what="prob")
            current.pos <- pos.tmp
        }
        else {
            thew <- rep(NA, length(pos.tmp)) ####
            for(kk in seq(along=pos.tmp)) {
                if(pos.tmp[kk] != current.pos[kk]) {
                    u <- abs(pos.tmp[kk]-pos[[kk]])
                    w <- indices[[kk]][u==min(u)]
                    if(length(w) > 1) {
                        warning("Confused about QTL positions.  You should probably run jittermap to ensure that no two markers conincide.")
                        w <- sample(w, 1)
                    }
                    if(method=="imp")
                        qtl.obj$geno[,kk,] <- cross$geno[[ichr[kk]]]$draws[,w,]
                    else
                        qtl.obj$prob[[kk]] <- cross$geno[[ichr[kk]]]$prob[,w,]
                    thew[kk] <- w ####

                    if(chrtype[ichr[kk]]=="X" && (type=="bc" || type=="f2")) {
                        if(method=="imp")
                            qtl.obj$geno[,kk,] <-
                                reviseXdata(type,"full",sexpgm,draws=qtl.obj$geno[,kk,,drop=FALSE],
                                            cross.attr=attributes(cross))
                        else {
                            temp <- qtl.obj$prob[[kk]]
                            temp <- array(temp, dim=c(nrow(temp),1,ncol(temp)))
                            dimnames(temp) <- list(NULL,"loc", 1:ncol(qtl.obj$prob[[kk]]))
                            qtl.obj$prob[[kk]] <- reviseXdata(type,"full",sexpgm,prob=temp,
                                                              cross.attr=attributes(cross))[,1,]
                        }
                    }

                    current.pos[kk] <- pos.tmp[kk]
                }
            }
        }
        # end of Karl's 8/23/05 addition

        # fit QTL, don't do drop one at a time
        fit <- fitqtlengine(pheno=pheno, qtl=qtl.obj, covar=covar,
                            formula=formula, method=method, model=model, dropone=FALSE,
                            get.ests=FALSE, run.checks=FALSE,
                            cross.attr=cross.attr, sexpgm=sexpgm, tol=tol, maxit=maxit,
                            forceXcovar=forceXcovar)

        matrix.rank[i] <- attr(fit, "matrix.rank")
        matrix.ncol[i] <- attr(fit, "matrix.ncol")

        if(verbose && ((i-1) %% n.prnt) == 0)
            cat("    ", i,"/", n.loop, "\n")

        # assign to result matrix
        #     Note: [[1]][1,4] picks out the LOD score
        result[i] <- fit[[1]][1,4]
    }

    # make the row and column names for the result matrix
    dnames <- list(NULL)
    for(i in 1:n.idx.varied) {
        i.chr <- chr[idx.varied[n.idx.varied-i+1]]
        i.pos <- pos[[idx.varied[n.idx.varied-i+1]]]
        dnames[[i]] <- paste( paste("Chr", i.chr,sep=""),
                             i.pos, sep="@")
    }
    dimnames(result) <- dnames

    class(result) <- "scanqtl"
    attr(result, "method") <- method
    attr(result, "formula") <- deparseQTLformula(formula)
    attr(result, "matrix.rank") <- matrix.rank
    attr(result, "matrix.ncol") <- matrix.ncol

    result
}


#summary.scanqtl <- function(object, ...)
#{
#}

#print.summary.qtl <- function(x, ...)
#{
#}

# end of scanqtl.R
