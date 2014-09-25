######################################################################
#
# addqtl.R
#
# copyright (c) 2007-2012, Karl W. Broman
# last modified Aug, 2012
# first written Nov, 2007
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
# Contains: addint, print.addint, addqtl, addpair,
#           reviseqtlnuminformula, qtlformulasymmetric
#           dropfromqtlformula
#           addcovarint, print.addcovarint, summary.addcovarint
#
######################################################################

######################################################################
# addint
#
# Try adding each possible QTL:QTL interaction (that is not
# already in the formula), and give results similar to the drop-one
# analysis.
######################################################################
addint <-
    function(cross, pheno.col=1, qtl, covar=NULL, formula,
             method=c("imp","hk"), model=c("normal", "binary"),
             qtl.only=FALSE, verbose=TRUE, pvalues=TRUE, simple=FALSE,
             tol=1e-4, maxit=1000, require.fullrank=FALSE)
{
    if( !("cross" %in% class(cross)) )
        stop("The cross argument must be an object of class \"cross\".")

    if( !("qtl" %in% class(qtl)) )
        stop("The qtl argument must be an object of class \"qtl\".")

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
        warning("addint can take just one phenotype; only the first will be used")
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

    if(qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if(method=="imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        else
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
    }
    if(method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3])  {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
    }

    # check phenotypes and covariates; drop ind'ls with missing values
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors=TRUE)
    if(any(is.na(phcovar))) {
        if(ncol(phcovar)==1) hasmissing <- is.na(phcovar)
        else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
        if(all(hasmissing))
            stop("All individuals are missing phenotypes or covariates.")
        if(any(hasmissing)) {
            warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
            pheno <- pheno[!hasmissing]
            qtl$n.ind <- sum(!hasmissing)
            if(method=="imp")
                qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
            else
                qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,,drop=FALSE])

            if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]
        }
    }

    # number of covariates
    if( is.null(covar) ) n.covar <- 0
    else n.covar <- ncol(covar)

    # if formula is missing, build one
    # all QTLs and covarariates will be additive by default
    n.qtl <- qtl$n.qtl
    if(missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
        formula <- "y~Q1"
        if(n.qtl > 1)
            for (i in 2:n.qtl)
                formula <- paste(formula, tmp.Q[i], sep="+")
        if (n.covar) { # if covarariate is not empty
            tmp.C <- colnames(covar) # covarariate term names
            for(i in 1:n.covar)
                formula <- paste(formula, tmp.C[i], sep="+")
        }
        formula <- as.formula(formula)
    }

    # check input formula
    formula <- checkformula(formula, qtl$altname, colnames(covar))

    # look for interactions that haven't been added
    factors <- attr(terms(formula), "factors")
    if(sum(factors[1,])==0) factors <- factors[-1,]

    # replace QTL altnames (Q1 etc) with real names (chr1@20 etc)
    fn <- fn.alt <- rownames(factors)
    qan <- qtl$altname
    qn <- qtl$name
    m <- match(fn, qan)
    fn.alt[!is.na(m)] <- qn[m[!is.na(m)]]

    # all possible interactions
    int2test <- int2test.alt <- NULL
    for(i in 1:(nrow(factors)-1)) {
        for(j in (i+1):nrow(factors)) {
            temp <- rep(0, nrow(factors))
            temp[c(i,j)] <- 1
            if(!any(apply(factors, 2, function(a, b) all(a==b), temp))) {
                int2test <- c(int2test, paste(fn[i], fn[j], sep=":"))
                int2test.alt <- c(int2test.alt, paste(fn.alt[i], fn.alt[j], sep=":"))
            }
        }
    }

    if(qtl.only && length(int2test) > 0) {
        z <- matrix(unlist(strsplit(int2test, ":")), ncol=2, byrow=TRUE)
        wh <- apply(z, 1, function(a)
                    length(grep("^[Qq][0-9]+$", a)) )
        int2test <- int2test[wh==2]
        int2test.alt <- int2test.alt[wh==2]
    }

    n2test <- length(int2test)

    if(n2test == 0) {
        if(verbose) cat("No pairwise interactions to add.\n")
        return(NULL)
    }

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

    # fit base model
    thefit0 <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar, formula=formula,
                            method=method, model=model, dropone=FALSE, get.ests=FALSE,
                            run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                            tol=tol, maxit=maxit)

    matrix0.rank <- attr(thefit0, "matrix.rank")
    matrix0.ncol <- attr(thefit0, "matrix.ncol")

    results <- matrix(ncol=7, nrow=n2test)
    dimnames(results) <- list(int2test.alt, c("df", "Type III SS", "LOD", "%var",
                                              "F value", "Pvalue(Chi2)", "Pvalue(F)"))

    matrix1.rank <- matrix1.ncol <- rep(0, n2test)

    for(k in seq(along=int2test)) {
        thefit1 <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar,
                                formula=as.formula(paste(deparseQTLformula(formula), int2test[k], sep="+")),
                                method=method, model=model, dropone=FALSE, get.ests=FALSE,
                                run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                                tol=tol, maxit=maxit)

        results[k,1] <- thefit1$result.full[1,1] - thefit0$result.full[1,1]
        results[k,2] <- thefit1$result.full[1,2] - thefit0$result.full[1,2]
        results[k,3] <- thefit1$result.full[1,4] - thefit0$result.full[1,4]

        results[k,4] <- 100*(1-10^(-2*thefit1$result.full[1,4]/qtl$n.ind)) -
            100*(1-10^(-2*thefit0$result.full[1,4]/qtl$n.ind))

        results[k,5] <- (results[k,2]/results[k,1])/thefit1$result.full[2,3]
        results[k,6] <- pchisq(results[k,3]*2*log(10), results[k,1], lower.tail=FALSE)
        results[k,7] <- pf(results[k,5], results[k,1], thefit1$result.full[3,1], lower.tail=FALSE)

        matrix1.rank[k] <- attr(thefit1, "matrix.rank")
        matrix1.ncol[k] <- attr(thefit1, "matrix.ncol")
    }
    matrix.fullrank <- (matrix1.rank - matrix0.rank == matrix1.ncol - matrix0.ncol)

    results <- as.data.frame(results, stringsAsFactors=TRUE)
    class(results) <- c("addint", "data.frame")
    attr(results, "method") <- method
    attr(results, "model") <- model
    attr(results, "formula") <- deparseQTLformula(formula)
    if(simple) pvalues <- FALSE
    attr(results, "pvalues") <- pvalues
    attr(results, "simple") <- simple
    attr(results, "matrix.fullrank") <- matrix.fullrank

    if(require.fullrank) results[!matrix.fullrank,3] <- 0

    results
}

print.addint <-
    function(x, ...)
{
    meth <- attr(x, "method")
    mod <- attr(x, "model")
    simp <- attr(x, "simple")
    if(is.null(mod)) mod <- "normal"
    if(is.null(meth)) meth <- "unknown"
    if(mod=="binary" || simp) attr(x, "pvalues") <- FALSE
    if(meth=="imp") meth <- "multiple imputation"
    else if(meth=="hk") meth <- "Haley-Knott regression"
    cat("Method:", meth, "\n")
    cat("Model: ", mod, "phenotype\n")

    cat("Model formula:")
    w <- options("width")[[1]]
    printQTLformulanicely(attr(x, "formula"), "                   ", w+5, w)
    cat("\n")

    cat("Add one pairwise interaction at a time table:\n")
    cat("--------------------------------------------\n")
    pval <- attr(x, "pvalues")
    if(!is.null(pval) && !pval)
        x <- x[,-ncol(x)+(0:1)]

    if(mod == "binary" || simp) x <- x[,c(1,3,4), drop=FALSE]

    printCoefmat(x, digits=4, cs.ind=1, P.values=pval, has.Pvalue=pval)

    cat("\n")
}

summary.addint <- function(object, ...) object


######################################################################
# addqtl
#
# scan for an additional QTL in the context of a multiple-QTL model
#
# If the formula includes one more QTL than in the QTL object, we
# use it as given; otherwise, a main effect for the additional QTL
# is added
#
# the output is like scanone
######################################################################
addqtl <-
    function(cross, chr, pheno.col=1, qtl, covar=NULL, formula,
             method=c("imp","hk"), model=c("normal", "binary"),
             incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
             forceXcovar=FALSE, require.fullrank=FALSE)
{
    method <- match.arg(method)
    model <- match.arg(model)

    if( !("cross" %in% class(cross)) )
        stop("The cross argument must be an object of class \"cross\".")

    if( !("qtl" %in% class(qtl)) )
        stop("The qtl argument must be an object of class \"qtl\".")

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

    if(verbose > 1) {
        verbose <- TRUE
        verbose.scanqtl <- TRUE
    }
    else verbose.scanqtl <- FALSE

    n.qtl <- qtl$n.qtl
    qtlchr <- qtl$chr
    qtlpos <- qtl$pos

    if(qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if(method=="imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        else
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
    }
    if(method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3])  {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
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

    # look for the chr
    if(missing(chr)) chr <- names(cross$geno)
    else chr <- matchchr(chr, names(cross$geno))

    # if formula is missing, make one.
    # All QTLs and covariates will be additive by default
    if(is.null(covar)) n.covar <- 0
    else n.covar <- ncol(covar)

    if(missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
        formula <- "y~Q1"
        if(n.qtl > 1)
            for (i in 2:n.qtl)
                formula <- paste(formula, tmp.Q[i], sep="+")
        if(n.covar) { # if covariate is not empty
            tmp.C <- names(covar) # covariate term names
            for(i in 1:n.covar)
                formula <- paste(formula, tmp.C[i], sep="+")
        }
        newformula <- as.formula(paste(formula, "+Q", n.qtl+1, sep=""))
        formula <- as.formula(formula)
    }
    else { # formula given
        newqtl <- paste("Q", n.qtl+1, sep="")

        # check the formula
        formula <- checkformula(formula, c(qtl$altname, newqtl), colnames(covar))

        theterms <- rownames(attr(terms(formula), "factors"))

        # is new QTL in the formula?
        g <- grep(paste("^[Qq]", n.qtl+1, "$", sep=""), theterms)
        if(length(g) == 0)  { # no; add to formula
            newformula <- as.formula(paste(deparseQTLformula(formula), "+ Q", n.qtl+1, sep=""))
        }
        else { # need a version without it
            newformula <- formula

            theterms <- colnames(attr(terms(formula), "factors"))
            g <- unique(c(grep(paste("^[Qq]", n.qtl+1, "$", sep=""), theterms),
                          grep(paste("^[Qq]", n.qtl+1, " *:", sep=""), theterms),
                          grep(paste(": +[Qq]", n.qtl+1, " *:", sep=""), theterms),
                          grep(paste(": +[Qq]", n.qtl+1, "$", sep=""), theterms)))

            if(length(g) > 0) {
                theterms <- theterms[-g]
                formula <- as.formula(paste("y ~ ", paste(theterms, collapse=" + "), sep=""))
            }
        }
    }

    # drop qtl that are not in the formula
    thefactors <- rownames(attr(terms(formula), "factors"))
    todrop <- NULL
    for(i in 1:n.qtl) {
        if(length(grep(paste("^[Qq]", i, "$", sep=""), thefactors))==0)
            todrop <- c(todrop, i)
    }

    if(length(todrop) > 0) {
        newqtlnum <- n.qtl+1

        notdropped <- (1:n.qtl)[-todrop]
        newnum <- 1:length(notdropped)

        qtl <- dropfromqtl(qtl, index=todrop)
        qtlchr <- qtlchr[-todrop]
        qtlpos <- qtlpos[-todrop]
        n.qtl <- n.qtl - length(todrop)
        revnewqtlnum <- n.qtl+1

        formula <- reviseqtlnuminformula(formula, notdropped, newnum)
        newformula <- reviseqtlnuminformula(newformula, c(notdropped,newqtlnum),
                                            c(newnum, revnewqtlnum))
    }

    # drop covariates that are not in the formula
    if(!is.null(covar)) {
        theterms <- rownames(attr(terms(formula), "factors"))
        m <- match(colnames(covar), theterms)
        if(all(is.na(m))) covar <- NULL
        else covar <- covar[,!is.na(m),drop=FALSE]
    }

    # phenotype column
    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("addqtl can take just one phenotype; only the first will be used")
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

    # check phenotypes and covariates; drop ind'ls with missing values
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors=TRUE)
    if(any(is.na(phcovar))) {
        if(ncol(phcovar)==1) hasmissing <- is.na(phcovar)
        else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
        if(all(hasmissing))
            stop("All individuals are missing phenotypes or covariates.")
        if(any(hasmissing)) {
            warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
            pheno <- pheno[!hasmissing]
            qtl$n.ind <- sum(!hasmissing)
            if(method=="imp")
                qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
            else
                qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,,drop=FALSE])

            if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]
            cross <- subset(cross, ind=!hasmissing)
        }
    }

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

    # fit the base model
    fit0 <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar, formula=formula,
                         method=method, model=model, dropone=FALSE, get.ests=FALSE,
                         run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                         tol=tol, maxit=maxit, forceXcovar=forceXcovar)
    lod0 <- fit0$result.full[1,4]
    matrix0.rank <- attr(fit0, "matrix.rank")
    matrix0.ncol <- attr(fit0, "matrix.ncol")

    results <- matrix1.rank <- matrix1.ncol <- NULL
    for(i in chr) {
        if(verbose) cat("Scanning chr", i, "\n")
        thechr <- c(qtlchr, i)
        thepos <- c(as.list(qtlpos), list(c(-Inf,Inf)))

        sqout <- scanqtl(cross, pheno.col=pheno.col, chr=thechr, pos=thepos,
                         covar=covar, formula=newformula, method=method, model=model,
                         incl.markers=incl.markers, verbose=verbose.scanqtl,
                         tol=tol, maxit=maxit)

        matrix1.rank <- c(matrix1.rank, attr(sqout, "matrix.rank"))
        matrix1.ncol <- c(matrix1.ncol, attr(sqout, "matrix.ncol"))


        # get map of positions
        if(method=="imp") {
            if("map" %in% names(attributes(cross$geno[[i]]$draws)))
                map <- attr(cross$geno[[i]]$draws,"map")
            else {
                stp <- attr(cross$geno[[i]]$draws, "step")
                oe <- attr(cross$geno[[i]]$draws, "off.end")

                if("stepwidth" %in% names(attributes(cross$geno[[i]]$draws)))
                    stpw <- attr(cross$geno[[i]]$draws, "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[i]]$map,stp,oe,stpw)
            }
        }
        else {
            if("map" %in% names(attributes(cross$geno[[i]]$prob)))
                map <- attr(cross$geno[[i]]$prob,"map")
            else {
                stp <- attr(cross$geno[[i]]$prob, "step")
                oe <- attr(cross$geno[[i]]$prob, "off.end")

                if("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
                    stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[i]]$map,stp,oe,stpw)
            }
        }
        # pull out the female map if there are sex-specific maps
        if(is.matrix(map)) map <- map[1,]

        if(method=="imp")
            step <- attr(cross$geno[[i]]$draws,"step")
        else
            step <- attr(cross$geno[[i]]$prob,"step")

        if(!incl.markers && step>0) { # equally spaced positions
            eq.sp.pos <- seq(min(map), max(map), by=step)
            wh.eq.pos <- match(eq.sp.pos, map)
            map <- map[wh.eq.pos]
        }

        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",i,".",w[o],sep="")

        z <- data.frame(lod=as.numeric(sqout)-lod0, stringsAsFactors=TRUE)
        z <- cbind(chr=rep(i,length(map)),
                   pos=as.numeric(map), z)
        rownames(z) <- w

        results <- rbind(results, z)
    }
    matrix.fullrank <- (matrix1.rank - matrix0.rank == matrix1.ncol - matrix0.ncol)

    class(results) <- c("scanone","data.frame")
    attr(results,"method") <- method
    attr(results,"formula") <- deparseQTLformula(newformula)

    attr(results, "matrix.fullrank") <- matrix.fullrank
    if(require.fullrank) results[!matrix.fullrank,3] <- 0

    results
}


######################################################################
# addpair
#
# scan for an additional pair of QTL in the context of a multiple-QTL
# model
#
# If the formula includes one more QTL than in the QTL object, we
# use it as given (perhaps adding the second, additively);
# otherwise, we do as in scantwo, performing both additive and
# interactive models, plus a single-QTL scan
#
# The output is like scantwo.  If we didn't do the scantwo type format,
# the results are placed where the full LODs usually are, and everything
# else is NA
######################################################################
addpair <-
    function(cross, chr, pheno.col=1, qtl, covar=NULL, formula,
             method=c("imp","hk"), model=c("normal", "binary"),
             incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
             forceXcovar=FALSE)
{
    method <- match.arg(method)
    model <- match.arg(model)

    if( !("cross" %in% class(cross)) )
        stop("The cross argument must be an object of class \"cross\".")

    if( !("qtl" %in% class(qtl)) )
        stop("The qtl argument must be an object of class \"qtl\".")

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

    if(verbose > 1) {
        verbose <- TRUE
        verbose.scanqtl <- TRUE
    }
    else verbose.scanqtl <- FALSE

    n.qtl <- qtl$n.qtl
    qtlchr <- qtl$chr
    qtlpos <- qtl$pos
    if(qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if(method=="imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        else
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
    }
    if(method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3]) {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
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

    # look for the chr
    if(missing(chr)) chr <- names(cross$geno)
    else chr <- matchchr(chr, names(cross$geno))

    fullmap <- pull.map(cross, chr)

    # if formula is missing, make one.
    # All QTLs and covariates will be additive by default
    if(is.null(covar)) n.covar <- 0
    else n.covar <- ncol(covar)

    if(missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
        formula <- "y~Q1"
        if(n.qtl > 1)
            for (i in 2:n.qtl)
                formula <- paste(formula, tmp.Q[i], sep="+")
        if(n.covar) { # if covariate is not empty
            tmp.C <- names(covar) # covariate term names
            for(i in 1:n.covar)
                formula <- paste(formula, tmp.C[i], sep="+")
        }
        newformula1 <- as.formula(paste(formula, " + Q", n.qtl+1,
                                        " + Q", n.qtl+2, " + Q", n.qtl+1,
                                        ":Q", n.qtl+2, sep=""))
        newformula2 <- as.formula(paste(formula, " + Q", n.qtl+1,
                                        " + Q", n.qtl+2, sep=""))
        formula <- as.formula(formula)
        scanbothways <- FALSE
    }
    else { # formula given
        newqtl <- paste("Q", n.qtl+1:2, sep="")

        # check the formula
        formula <- checkformula(formula, c(qtl$altname, newqtl), colnames(covar))

        theterms <- rownames(attr(terms(formula), "factors"))

        # are either of the new QTL in the formula?
        g <- c(grep(paste("^[Qq]", n.qtl+1, "$", sep=""), theterms),
               grep(paste("^[Qq]", n.qtl+2, "$", sep=""), theterms))
        g1 <- grep(paste("^[Qq]", n.qtl+1, "$", sep=""), theterms)
        g2 <- grep(paste("^[Qq]", n.qtl+2, "$", sep=""), theterms)
        if(length(g) == 0)  { # no; add to formula
            newformula1 <- as.formula(paste(deparseQTLformula(formula), "+ Q", n.qtl+1,
                                            " + Q", n.qtl+2, " + Q", n.qtl+1,
                                            ":Q", n.qtl+2, sep=""))
            newformula2 <- as.formula(paste(deparseQTLformula(formula), "+ Q", n.qtl+1,
                                            " + Q", n.qtl+2, sep=""))
            scanbothways <- FALSE
        }
        else { # need a version without them
            # first make sure that *both* terms are in the formula
            if(length(g1)==0) # add Q1
                formula <- as.formula(paste(deparseQTLformula(formula), "+ Q", n.qtl+1, sep=""))
            if(length(g2)==0) # add Q2
                formula <- as.formula(paste(deparseQTLformula(formula), "+ Q", n.qtl+2, sep=""))

            newformula1 <- formula
            newformula2 <- NULL

            theterms <- colnames(attr(terms(formula), "factors"))
            g <- unique(c(grep(paste("^[Qq]", n.qtl+1, "$", sep=""), theterms),
                          grep(paste("^[Qq]", n.qtl+1, " *:", sep=""), theterms),
                          grep(paste(": *[Qq]", n.qtl+1, " *:", sep=""), theterms),
                          grep(paste(": *[Qq]", n.qtl+1, "$", sep=""), theterms),
                          grep(paste("^[Qq]", n.qtl+2, "$", sep=""), theterms),
                          grep(paste("^[Qq]", n.qtl+2, " *:", sep=""), theterms),
                          grep(paste(": *[Qq]", n.qtl+2, " *:", sep=""), theterms),
                          grep(paste(": *[Qq]", n.qtl+2, "$", sep=""), theterms)))
            if(length(g) > 0) {
                theterms <- theterms[-g]
                formula <- as.formula(paste("y ~ ", paste(theterms, collapse=" + "), sep=""))
            }

            # if the QTL formula is symmetric in the two new QTL, need scan only for i<j
            # if not symmetric, we also need to scan with j<i
            scanbothways <- !qtlformulasymmetric(newformula1, n.qtl+1, n.qtl+2)

            if(scanbothways) {
                newformula1.minus1 <- dropfromqtlformula(newformula1, n.qtl+1)
                newformula1.minus2 <- dropfromqtlformula(newformula1, n.qtl+2)
                newformula1.minus1 <- reviseqtlnuminformula(newformula1.minus1, n.qtl+2, n.qtl+1)
            }
        }
    }

    # drop qtl that are not in the formula
    thefactors <- rownames(attr(terms(formula), "factors"))
    todrop <- NULL
    for(i in 1:n.qtl) {
        if(length(grep(paste("^[Qq]", i, "$", sep=""), thefactors))==0)
            todrop <- c(todrop, i)
    }
    if(length(todrop) > 0) {
        newqtlnum1 <- n.qtl+1
        newqtlnum2 <- n.qtl+2

        notdropped <- (1:n.qtl)[-todrop]
        newnum <- 1:length(notdropped)

        qtl <- dropfromqtl(qtl, index=todrop)
        qtlchr <- qtlchr[-todrop]
        qtlpos <- qtlpos[-todrop]
        n.qtl <- n.qtl - length(todrop)
        revnewqtlnum1 <- n.qtl+1
        revnewqtlnum2 <- n.qtl+2

        formula <- reviseqtlnuminformula(formula, notdropped, newnum)
        newformula1 <- reviseqtlnuminformula(newformula1, c(notdropped, newqtlnum1, newqtlnum2),
                                             c(newnum, revnewqtlnum1, revnewqtlnum2))
        if(!is.null(newformula2))
            newformula2 <- reviseqtlnuminformula(newformula2, c(notdropped, newqtlnum1, newqtlnum2),
                                                 c(newnum, revnewqtlnum1, revnewqtlnum2))
        if(scanbothways) {
            newformula1.minus1 <- reviseqtlnuminformula(newformula1.minus1,
                                                        c(notdropped, newqtlnum1),
                                                        c(newnum, revnewqtlnum1))

            newformula1.minus2 <- reviseqtlnuminformula(newformula1.minus2,
                                                        c(notdropped, newqtlnum1),
                                                        c(newnum, revnewqtlnum1))
        }
    }

    # drop covariates that are not in the formula
    if(!is.null(covar)) {
        theterms <- rownames(attr(terms(formula), "factors"))
        m <- match(colnames(covar), theterms)
        if(all(is.na(m))) covar <- NULL
        else covar <- covar[,!is.na(m),drop=FALSE]
    }

    # phenotype column
    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("addpair can take just one phenotype; only the first will be used")
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

    # check phenotypes and covariates; drop ind'ls with missing values
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors=TRUE)
    if(any(is.na(phcovar))) {
        if(ncol(phcovar)==1) hasmissing <- is.na(phcovar)
        else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
        if(all(hasmissing))
            stop("All individuals are missing phenotypes or covariates.")
        if(any(hasmissing)) {
            warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
            pheno <- pheno[!hasmissing]
            qtl$n.ind <- sum(!hasmissing)
            if(method=="imp")
                qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
            else
                qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,,drop=FALSE])

            if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]
            cross <- subset(cross, ind=!hasmissing)
        }
    }

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

    # fit the base model
    fit0 <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar, formula=formula,
                         method=method, model=model, dropone=FALSE, get.ests=FALSE,
                         run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                         tol=tol, maxit=maxit, forceXcovar=forceXcovar)
    lod0 <- fit0$result.full[1,4]

    gmap <- NULL

    # form map
    for(i in 1:length(chr)) {
        ci <- chr[i]
        # get map of positions
        if(method=="imp") {
            if("map" %in% names(attributes(cross$geno[[ci]]$draws)))
                map <- attr(cross$geno[[ci]]$draws,"map")
            else {
                stp <- attr(cross$geno[[ci]]$draws, "step")
                oe <- attr(cross$geno[[ci]]$draws, "off.end")

                if("stepwidth" %in% names(attributes(cross$geno[[ci]]$draws)))
                    stpw <- attr(cross$geno[[ci]]$draws, "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[ci]]$map,stp,oe,stpw)
            }
        }
        else {
            if("map" %in% names(attributes(cross$geno[[ci]]$prob)))
                map <- attr(cross$geno[[ci]]$prob,"map")
            else {
                stp <- attr(cross$geno[[ci]]$prob, "step")
                oe <- attr(cross$geno[[ci]]$prob, "off.end")

                if("stepwidth" %in% names(attributes(cross$geno[[ci]]$prob)))
                    stpw <- attr(cross$geno[[ci]]$prob, "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[ci]]$map,stp,oe,stpw)
            }
        }
        # pull out the female map if there are sex-specific maps
        if(is.matrix(map)) map <- map[1,]

        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",ci,".",w[o],sep="")
        map <- cbind(chr=rep(ci,length(map)),
                     pos=as.data.frame(map, stringsAsFactors=TRUE) )
        rownames(map) <- w

        if(method=="imp")
            step <- attr(cross$geno[[ci]]$draws,"step")
        else
            step <- attr(cross$geno[[ci]]$prob,"step")

        if(step==0 || stepwidth.var)  # just use markers
            eq.sp.pos <- rep(1,nrow(map))
        else {
            eq.sp.pos <- seq(min(map[,2]),max(map[,2]),by=step)
            wh.eq.sp <- match(eq.sp.pos,map[,2])
            if(any(is.na(wh.eq.sp))) { # this shouldn't happen
                warning("Possible error in determining the equally spaced positions.")
                wh.eq.sp <- wh.eq.sp[!is.na(wh.eq.sp)]
            }
            eq.sp.pos <- rep(0,nrow(map))
            eq.sp.pos[wh.eq.sp] <- 1
        }
        if(!incl.markers && any(eq.sp.pos==0)) {
            map <- map[eq.sp.pos==1,]
            eq.sp.pos <- eq.sp.pos[eq.sp.pos==1]
        }
        gmap <- rbind(gmap, cbind(map,
                                  eq.spacing=eq.sp.pos,
                                  xchr=(class(cross$geno[[i]])=="X")))
    }

    lod <- matrix(ncol=nrow(gmap), nrow=nrow(gmap))
    if(scanbothways) lod.m1 <- lod.m2 <- rep(NA, nrow(gmap))

    for(i in 1:length(chr)) {
        ci <- chr[i]
        whi <- which(gmap[,1]==ci)
        for(j in i:length(chr)) {
            cj <- chr[j]
            whj <- which(gmap[,1]==cj)

            if(verbose) {
                if(is.null(newformula2))
                    cat("Scanning chr", ci, "and", cj, "\n")
                else
                    cat("Scanning full model for chr", ci, "and", cj, "\n")
            }

            thechr <- c(qtlchr, ci, cj)
            thepos <- c(as.list(qtlpos), list(c(-Inf, Inf)), list(c(-Inf, Inf)))

            temp1 <- scanqtl(cross, pheno.col=pheno.col, chr=thechr, pos=thepos,
                             covar=covar, formula=newformula1, method=method, model=model,
                             incl.markers=incl.markers, verbose=verbose.scanqtl,
                             tol=tol, maxit=maxit) - lod0

            if(!is.null(newformula2)) {
                if(verbose)
                    cat("Scanning add've model for chr", ci, "and", cj, "\n")

                temp2 <- scanqtl(cross, pheno.col=pheno.col, chr=thechr, pos=thepos,
                                 covar=covar, formula=newformula2, method=method, model=model,
                                 incl.markers=incl.markers, verbose=verbose.scanqtl,
                                 tol=tol, maxit=maxit) - lod0
            }
            else {
                if(i != j && scanbothways) {
                    if(verbose) cat("Scanning chr", cj, "and", ci, "\n")
                    thechr <- c(qtlchr, cj, ci)
                    temp1r <- scanqtl(cross, pheno.col=pheno.col, chr=thechr, pos=thepos,
                                      covar=covar, formula=newformula1, method=method, model=model,
                                      incl.markers=incl.markers, verbose=verbose.scanqtl,
                                      tol=tol, maxit=maxit) - lod0
                }
            }

            if(i==j) {
                if(!is.null(newformula2))
                    temp1[upper.tri(temp1)] <- temp2[upper.tri(temp1)]
                else temp1 <- t(temp1)

                lod[whi,whi] <- temp1
            }
            else {
                if(!is.null(newformula2)) {
                    lod[whi,whj] <- t(temp2)
                    lod[whj,whi] <- temp1
                }
                else {
                    lod[whi,whj] <- t(temp1)
                    if(scanbothways)
                        lod[whj,whi] <- t(temp1r)
                    else
                        lod[whj,whi] <- temp1
                }
            }

        }

        if(scanbothways) {
            if(verbose) cat("Scanning chr", ci, "alone\n")

            thechr <- c(qtlchr, ci)
            thepos <- c(as.list(qtlpos), list(c(-Inf, Inf)))

            lod.m1[whi] <- scanqtl(cross, pheno.col=pheno.col, chr=thechr, pos=thepos,
                                   covar=covar, formula=newformula1.minus1, method=method,
                                   model=model, incl.markers=incl.markers,
                                   verbose=verbose.scanqtl, tol=tol, maxit=maxit) - lod0
            lod.m2[whi] <- scanqtl(cross, pheno.col=pheno.col, chr=thechr, pos=thepos,
                                   covar=covar, formula=newformula1.minus2, method=method,
                                   model=model, incl.markers=incl.markers,
                                   verbose=verbose.scanqtl, tol=tol, maxit=maxit) - lod0

        }

    }

    result <- list(lod=lod, map=gmap, scanoneX=NULL)
    class(result) <- "scantwo"
    attr(result, "fullmap") <- fullmap
    attr(result,"method") <- method
    attr(result,"formula") <- deparseQTLformula(newformula1)

    if(scanbothways) {
        attr(result, "lod.minus1") <- lod.m1
        attr(result, "lod.minus2") <- lod.m2
    }

    if(is.null(newformula2)) attr(result, "addpair") <- TRUE

    result
}



######################################################################
# consider, e.g., oldnum=5 and newnum = 3
# This functions replaces any instances of "Q5" or "q5" in the
# formula with "Q3".
######################################################################
reviseqtlnuminformula <-
    function(formula, oldnum, newnum)
{
    if(is.character(formula))
        formula <- as.formula(formula)

    if(length(oldnum) != length(newnum))
        stop("oldnum and newnum must be the same length.")

    newterms <- theterms <- colnames(attr(terms(formula), "factors"))

    for(i in seq(along=oldnum)) {
        g <- grep(paste("^[Qq]", oldnum[i], "$", sep=""), theterms)
        if(length(g) > 0)
            newterms[g] <- paste("Q", newnum[i], sep="")
    }

    intxn <- grep(":", theterms)
    if(length(intxn) > 0) {
        temp <- strsplit(theterms[intxn], ":")
        for(i in seq(along=temp)) {
            for(j in seq(along=oldnum)) {
                g <- grep(paste("^[Qq]", oldnum[j], "$", sep=""), temp[[i]])
                if(length(g) > 0)
                    temp[[i]][g] <- paste("Q", newnum[j], sep="")
            }
        }
        temp <- sapply(temp, paste, collapse=":")
        newterms[intxn] <- temp
    }
    as.formula(paste("y ~ ", paste(newterms, collapse=" + "), sep=""))
}

######################################################################
# return TRUE or FALSE according to whether the formula is symmetric
# in QTL qtlnum1 and qtlnum2
######################################################################
qtlformulasymmetric <-
    function(formula, qtlnum1, qtlnum2)
{
    theterms <- attr(terms(formula), "factors")
    rn <- rownames(theterms)
    wh1 <- grep(paste("^[Qq]", qtlnum1, "$", sep=""), rn)
    wh2 <- grep(paste("^[Qq]", qtlnum2, "$", sep=""), rn)
    cn <- colnames(theterms)

    if(length(wh1)==0 && length(wh2)==0) return(TRUE)
    if(length(wh1)==0 || length(wh2)==0) return(FALSE)

    revterms <- theterms
    revterms[c(wh1,wh2),] <- revterms[c(wh2, wh1),]

    theterms <- sort(apply(theterms, 2, paste, collapse=""))
    revterms <- sort(apply(revterms, 2, paste, collapse=""))

    all(theterms==revterms)
}


######################################################################
# If qtlnum=5, drop any terms containing Q5 or q5
######################################################################
dropfromqtlformula <-
    function(formula, qtlnum)
{
    theterms <- colnames(attr(terms(formula), "factors"))

    todrop <- NULL

    for(i in seq(along=qtlnum)) {
        g <- grep(paste("^[Qq]", qtlnum[i], "$", sep=""), theterms)
        if(length(g) > 0) todrop <- c(todrop, g)
    }

    intxn <- grep(":", theterms)
    if(length(intxn) > 0) {
        temp <- strsplit(theterms[intxn], ":")
        for(i in seq(along=temp)) {
            for(j in seq(along=qtlnum)) {
                g <- grep(paste("^[Qq]", qtlnum[j], "$", sep=""), temp[[i]])
                if(length(g) > 0) todrop <- c(todrop, intxn[i])
            }
        }
    }

    todrop <- unique(todrop)

    as.formula(paste("y ~ ", paste(theterms[-todrop], collapse=" + "), sep=""))
}



######################################################################
# addcovarint
#
# Try adding each QTL x covariate interaction (that is not
# already in the formula), and give results similar to the drop-one
# analysis.
######################################################################
addcovarint <-
    function(cross, pheno.col=1, qtl, covar=NULL, icovar, formula,
             method=c("imp","hk"), model=c("normal", "binary"),
             verbose=TRUE, pvalues=TRUE, simple=FALSE, tol=1e-4, maxit=1000,
             require.fullrank=FALSE)
{
    if( !("cross" %in% class(cross)) )
        stop("The cross argument must be an object of class \"cross\".")

    if( !("qtl" %in% class(qtl)) )
        stop("The qtl argument must be an object of class \"qtl\".")

    if(missing(covar) || is.null(covar))
        stop("Must include covariate data frame.")
    if(!is.data.frame(covar)) {
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
        warning("addcovarint can take just one phenotype; only the first will be used")
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
    if(nrow(covar) != length(pheno))
        stop("nrow(covar) != no. individuals in cross.")

    if(missing(icovar))
        stop("Must include icovar (the covariate to consider in interactions)")

    if(!is.character(icovar) || any(is.na(match(icovar, colnames(covar)))))
        stop("icovar must be a vector of character strings corresonding to columns in covar.")

    method <- match.arg(method)
    model <- match.arg(model)

    # allow formula to be a character string
    if(!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

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

    if(qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if(method=="imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        else
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
    }
    if(method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3])  {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
    }

    # check phenotypes and covariates; drop ind'ls with missing values
    phcovar <- cbind(pheno, covar)
    if(any(is.na(phcovar))) {
        if(ncol(phcovar)==1) hasmissing <- is.na(phcovar)
        else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
        if(all(hasmissing))
            stop("All individuals are missing phenotypes or covariates.")
        if(any(hasmissing)) {
            warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
            pheno <- pheno[!hasmissing]
            qtl$n.ind <- sum(!hasmissing)
            if(method=="imp")
                qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
            else
                qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,,drop=FALSE])

            covar <- covar[!hasmissing,,drop=FALSE]
        }
    }

    # number of covariates
    n.covar <- ncol(covar)

    # if formula is missing, build one
    # all QTLs and covarariates will be additive by default
    n.qtl <- qtl$n.qtl
    if(missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
        formula <- "y~Q1"
        if(n.qtl > 1)
            for (i in 2:n.qtl)
                formula <- paste(formula, tmp.Q[i], sep="+")
        if (n.covar) { # if covarariate is not empty
            tmp.C <- colnames(covar) # covarariate term names
            for(i in 1:n.covar)
                formula <- paste(formula, tmp.C[i], sep="+")
        }
        formula <- as.formula(formula)
    }

    # check input formula
    formula <- checkformula(formula, qtl$altname, colnames(covar))

    # make sure icovar is in the formula
    m <- is.na(match(icovar, rownames(attr(terms(formula), "factors"))))
    if(any(m))
        formula <- as.formula(paste(deparseQTLformula(formula), "+",
                                    paste(icovar[m], collapse="+"), sep=""))

    # look for interactions that haven't been added
    factors <- attr(terms(formula), "factors")
    if(sum(factors[1,])==0) factors <- factors[-1,]

    # replace QTL altnames (Q1 etc) with real names (chr1@20 etc)
    fn <- fn.alt <- rownames(factors)
    qan <- qtl$altname
    qn <- qtl$name
    m <- match(fn, qan)
    fn.alt[!is.na(m)] <- qn[m[!is.na(m)]]

    theqtl <- fn[fn != fn.alt]
    theqtl.alt <- fn.alt[fn != fn.alt]

    theint <- theint.alt <- NULL
    for(i in icovar) {
        theint <- c(theint, paste(theqtl, ":", i, sep=""))
        theint.alt <- c(theint.alt, paste(theqtl.alt, ":", i, sep=""))
    }

    wh <- match(theint, colnames(factors))
    theint <- theint[is.na(wh)]
    theint.alt <- theint.alt[is.na(wh)]

    n2test <- length(theint)

    if(n2test == 0) {
        if(verbose) cat("No QTL x covariate interactions to add.\n")
        return(NULL)
    }

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

    # fit base model
    thefit0 <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar, formula=formula,
                            method=method, model=model, dropone=FALSE, get.ests=FALSE,
                            run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                            tol=tol, maxit=maxit)
    matrix0.rank <- attr(thefit0, "matrix.rank")
    matrix0.ncol <- attr(thefit0, "matrix.ncol")

    results <- matrix(ncol=7, nrow=n2test)
    dimnames(results) <- list(theint.alt, c("df", "Type III SS", "LOD", "%var",
                                            "F value", "Pvalue(Chi2)", "Pvalue(F)"))

    matrix1.rank <- matrix1.ncol <- rep(0, n2test)
    for(k in seq(along=theint)) {
        thefit1 <- fitqtlengine(pheno=pheno, qtl=qtl, covar=covar,
                                formula=as.formula(paste(deparseQTLformula(formula), theint[k], sep="+")),
                                method=method, model=model, dropone=FALSE, get.ests=FALSE,
                                run.checks=FALSE, cross.attr=cross.attr, sexpgm=sexpgm,
                                tol=tol, maxit=maxit)

        results[k,1] <- thefit1$result.full[1,1] - thefit0$result.full[1,1]
        results[k,2] <- thefit1$result.full[1,2] - thefit0$result.full[1,2]
        results[k,3] <- thefit1$result.full[1,4] - thefit0$result.full[1,4]

        results[k,4] <- 100*(1-10^(-2*thefit1$result.full[1,4]/qtl$n.ind)) -
            100*(1-10^(-2*thefit0$result.full[1,4]/qtl$n.ind))

        results[k,5] <- (results[k,2]/results[k,1])/thefit1$result.full[2,3]
        results[k,6] <- pchisq(results[k,3]*2*log(10), results[k,1], lower.tail=FALSE)
        results[k,7] <- pf(results[k,5], results[k,1], thefit1$result.full[3,1], lower.tail=FALSE)
        matrix1.rank[k] <- attr(thefit1, "matrix.rank")
        matrix1.ncol[k] <- attr(thefit1, "matrix.ncol")
    }
    matrix.fullrank <- (matrix1.rank - matrix0.rank == matrix1.ncol - matrix0.ncol)

    results <- as.data.frame(results, stringsAsFactors=TRUE)
    class(results) <- c("addcovarint", "data.frame")
    attr(results, "model") <- model
    attr(results, "method") <- method
    attr(results, "formula") <- deparseQTLformula(formula)
    if(simple) pvalues <- FALSE
    attr(results, "pvalues") <- pvalues
    attr(results, "simple") <- simple

    attr(results, "matrix.fullrank") <- matrix.fullrank
    if(require.fullrank) results[!matrix.fullrank,3] <- 0

    results
}

print.addcovarint <-
    function(x, ...)
{
    meth <- attr(x, "method")
    mod <- attr(x, "model")
    simp <- attr(x, "simple")
    if(is.null(mod)) mod <- "normal"
    if(is.null(meth)) meth <- "unknown"
    if(mod=="binary" || simp) attr(x, "pvalues") <- FALSE
    if(meth=="imp") meth <- "multiple imputation"
    else if(meth=="hk") meth <- "Haley-Knott regression"
    cat("Method:", meth, "\n")
    cat("Model: ", mod, "phenotype\n")

    cat("Model formula:")
    w <- options("width")[[1]]
    printQTLformulanicely(attr(x, "formula"), "                   ", w+5, w)
    cat("\n")

    cat("Add one QTL x covar interaction at a time table:\n")
    cat("--------------------------------------------\n")
    pval <- attr(x, "pvalues")
    if(!is.null(pval) && !pval)
        x <- x[,-ncol(x)+(0:1)]

    if(mod == "binary" || simp) x <- x[,c(1,3,4), drop=FALSE]

    printCoefmat(x, digits=4, cs.ind=1, P.values=pval, has.Pvalue=pval)

    cat("\n")
}

summary.addcovarint <- function(object, ...) object

# end of addqtl.R
