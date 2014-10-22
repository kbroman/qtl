######################################################################
# stepwiseqtlX.R
#
# copyright (c) 2013-2014, Karl W Broman and Quoc Tran
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
#
# Functions related to stepwise QTL analysis but with separate
# penalties for autosome and X chromosome
#
# Contains: countqtltermsX, calc.penalties.X, stepwiseqtlX, calc.plod.X
######################################################################

######################################################################
# count terms in a model, for use by plod, Quoc editted
######################################################################
countqtltermsX <-
    function(formula, qtl, ignore.covar=TRUE) # qtl is used to extract a logical vector: on "X" chr == 1
{
    if(is.character(formula)) formula <- as.formula(formula)
    factors <- attr(terms(formula), "factors")[-1,,drop=FALSE]
    if(any(factors > 1))  {
        warning("some formula terms > 1; may be a problem with the formula:\n    ", deparseQTLformula(formula))
        factors[factors > 1] <- 1
    }
    nterm <- apply(factors, 2, sum)

    if(any(nterm>2))
        stop("Can't deal with higher-order interactions\n")

    # need to check for QTL x covariate interactions in here!

    if(ignore.covar) {
        cn <- colnames(factors)
        wh <- c(grep("^[Qq][0-9]+$", cn),
                grep("^[Qq][0-9]+:[Qq][0-9]+$", cn))
        rn <- rownames(factors)
        wh2 <- c(grep("^[Qq][0-9]+$", rn),
                 grep("^[Qq][0-9]+:[Qq][0-9]+$", rn))
        factors <- factors[wh2,wh, drop=FALSE]
    }
    nterm <- apply(factors, 2, sum)

    nmain <- sum(nterm==1)

    qtltype <- qtl$chrtype
    names(qtltype) <- qtl$altname
    factors.main <- factors[,nterm==1, drop=FALSE] # Get first order factors
    nmainA <- sum(qtltype[colnames(factors.main)]=="A")
    nmainX <- sum(qtltype[colnames(factors.main)]=="X")


    if(all(nterm==1))
        return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=0, intAX=0, intXX=0, inttot=0))

    n.int <- sum(nterm==2)

    if(n.int <=0) # 0 interactions, so no need to figure them out
        return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=n.int, intAX=0, intXX=0, inttot=n.int))

    factors <- factors[, nterm==2, drop=FALSE] # Get second order factors

    wh <- apply(factors, 2, function(a) which(a==1)) # Get the location of interaction, check for X here.

    if(n.int ==1){ # 1 interaction, check if it is AA AX or XX
        int.type <- sum(qtltype[rownames(factors)[wh[,1]]]=="X") #0 is AA, 1 is AX, 2 is XX
        if (int.type == 0) return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=1, intAX=0, intXX=0, inttot=n.int))
        else
            if (int.type == 1) return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=0, intAX=1, intXX=0, inttot=n.int))
            else return(c(mainA=nmainA, mainX=nmainX, intH=0, intL=0, intAX=0, intXX=1, inttot=n.int))
    }

    u <- sort(unique(as.numeric(wh))) # Unique nodes with interactions, in increasing order; does not contain isolated nodes
    grp <- rep(NA, length(u)) # Group member, length = number of node
    names(grp) <- u

    ngrp <- 0 # Number of group
    nintAA <- NULL # Number of AA interaction of specific group
    nintAX <- NULL # Number of AX interaction of specific group
    nintXX <- NULL # Number of XX interaction of specific group


    for(i in 1:ncol(wh)) {
        thegrp <- grp[as.character(wh[,i])]
        int.type <- sum(qtltype[rownames(factors)[wh[,i]]]=="X") #0 is AA, 1 is AX, 2 is XX
        if(all(!is.na(thegrp))) { # Merge 2 groups
            nintAA[as.character(thegrp[1])] <-
                sum(nintAA[unique(as.character(thegrp))])
            nintAX[as.character(thegrp[1])] <-
                sum(nintAX[unique(as.character(thegrp))])
            nintXX[as.character(thegrp[1])] <-
                sum(nintXX[unique(as.character(thegrp))])
            if (int.type==0) nintAA[as.character(thegrp[1])] <- nintAA[as.character(thegrp[1])] + 1
            if (int.type==1) nintAX[as.character(thegrp[1])] <- nintAX[as.character(thegrp[1])] + 1
            if (int.type==2) nintXX[as.character(thegrp[1])] <- nintXX[as.character(thegrp[1])] + 1
            grp[grp==thegrp[1] | grp==thegrp[2]] <- thegrp[1]
        } # two connected group become one group and has the name of the first group, number of interaction is updated
        else if(any(!is.na(thegrp))) { # add 1 more interaction to current group
            grp[as.character(wh[,i])] <- thegrp[!is.na(thegrp)]
            if (int.type==0) nintAA[as.character(thegrp[!is.na(thegrp)])] <-
                nintAA[as.character(thegrp[!is.na(thegrp)])] + 1
            if (int.type==1) nintAX[as.character(thegrp[!is.na(thegrp)])] <-
                nintAX[as.character(thegrp[!is.na(thegrp)])] + 1
            if (int.type==2) nintXX[as.character(thegrp[!is.na(thegrp)])] <-
                nintXX[as.character(thegrp[!is.na(thegrp)])] + 1
        }
        else { # introduce new group
            ngrp <- ngrp+1
            grp[as.character(wh[,i])] <- ngrp
            # Initialize nint
            nintAA[as.character(ngrp)] <- 0
            nintAX[as.character(ngrp)] <- 0
            nintXX[as.character(ngrp)] <- 0
            if (int.type==0) nintAA[as.character(ngrp)] <- 1
            if (int.type==1) nintAX[as.character(ngrp)] <- 1
            if (int.type==2) nintXX[as.character(ngrp)] <- 1
        }
    }

    nintAA <- nintAA[as.character(unique(grp))]
    nintAX <- nintAX[as.character(unique(grp))]
    nintXX <- nintXX[as.character(unique(grp))]
    nL <- sum(nintAA>0)
    nH <- sum(nintAA)-nL
    c(mainA=nmainA, mainX=nmainX, intH=nH, intL=nL,
      intAX=sum(nintAX), intXX=sum(nintXX),
      inttot=nH+nL+sum(nintAX)+sum(nintXX))
}

######################################################################
# calculate penalties for pLOD using scantwo permutation results,  Quoc editted
######################################################################
calc.penalties.X <-
    function(perms, alpha=0.05, lodcolumn)
{
    if(missing(perms) || !("scantwoperm" %in% class(perms)))
        stop("You must include permutation results from scantwo.")

    if(!("AA" %in% names(perms)))
        stop("perms need to be X-chr-specific")

    if(length(alpha) > 1) {
        alpha <- alpha[1]
        warning("alpha needs to be a single value; only the first will be used.")
    }

    if(missing(lodcolumn) || is.null(lodcolumn)) {
        if(is.matrix(perms[[1]][[1]]) && ncol(perms[[1]][[1]]) > 1)
            lodcolumn <- 1:ncol(perms[[1]][[1]])
        else lodcolumn <- 1
    }

    if(length(lodcolumn)>1) {
        result <- NULL
        for(i in seq(along=lodcolumn)) {
            temp <- calc.penalties.X(perms, alpha, lodcolumn[i])
            result <- rbind(result, temp)
        }
        dimnames(result) <- list(colnames(perms[[1]][[1]])[lodcolumn], names(temp))
        return(result)
    }

    if(is.matrix(perms[[1]][[1]]) && ncol(perms[[1]][[1]]) >1) {
        if(lodcolumn < 1 || lodcolumn > ncol(perms[[1]][[1]]))
            stop("lodcolumn misspecified")
        for(j in 1:3) {
            for(i in seq(along=perms[[j]]))
                perms[[j]][[i]] <- perms[[j]][[i]][,lodcolumn,drop=FALSE]
        }
    }

    qu <- summary(perms, alpha=alpha)

    c(mainA=qu$AA$one,
      mainX=qu$XX$one,
      intH=qu$AA$int,
      intL=qu$AA$fv1 - qu$AA$one,
      intAX=qu$AX$int,
      intXX=qu$XX$int)
}

######################################################################
# stepwiseqtlX, work for data with both A and X chromosomes,
# Need a check for X chromosome to switch back to standard stepwiseqtl
# This is derived from stepwiseqtl (R/qtl v1.25.14)
# 2 more variables than stepwiseqtl:
#   stop.rule (default at 0), 0 - no early stopping rule,
#                             1 - early stoping rule; plod<-k_f*mainA,
#                             2 - early stoping rule; plod<bestplod-k_f*mainA.
#   k_f (default at 3) number of main penalty plod can go below.
######################################################################
stepwiseqtlX <-
    function(cross, chr, pheno.col=1, qtl, formula, max.qtl=10, k_f=3, stop.rule=0, covar=NULL,
             method=c("imp", "hk"), model=c("normal", "binary"), incl.markers=TRUE, refine.locations=TRUE,
             additive.only=FALSE, penalties,
             keeplodprofile=TRUE, keeptrace=FALSE, verbose=TRUE,
             tol=1e-4, maxit=1000, require.fullrank=TRUE)
{
    if(verbose) message("Running stepwiseqtlX")

    if(!("cross" %in% class(cross)))
        stop("Input should have class \"cross\".")

    if(missing(penalties))
        stop("penalties must be provided")

    if(missing(qtl)) qtl <- NULL
    if(missing(formula)) formula <- NULL

    if(!missing(chr)) cross <- subset(cross, chr)

    method <- match.arg(method)
    model <- match.arg(model)

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    chrtype <- sapply(cross$geno, class)
    if(any(chrtype=="X")) {
        Xadjustment <- scanoneXnull(class(cross)[1], getsex(cross))
        forceXcovar <- Xadjustment$adjustX
        Xcovar <- Xadjustment$sexpgmcovar
    } else forceXcovar <- FALSE

    if(!is.null(qtl)) { # start f.s. at somewhere other than the null
        if( !("qtl" %in% class(qtl)) )
            stop("The qtl argument must be an object of class \"qtl\".")

        # check that chromosomes were retained, otherwise give error
        m <- is.na(match(qtl$chr, names(cross$geno)))
        if(any(m)) {
            wh <- qtl$chr[m]
            if(length(wh) > 1)
                stop("Chromosomes ", paste(wh, collapse=", "), " (in QTL object) not in cross object.")
            else
                stop("Chromosome ", wh, " (in QTL object) not in cross object.")
        }
        if(is.null(formula)) { # create a formula with all covariates and all QTL add've
            if(!is.null(covar))
                formula <- paste("y ~ ", paste(names(covar), collapse="+"), "+")
            else
                formula <- "y ~ "
            formula <- paste(formula, paste(paste("Q", 1:length(qtl$chr), sep=""), collapse="+"))
        }
        else {
            temp <- checkStepwiseqtlStart(qtl, formula, covar)
            qtl <- temp$qtl
            formula <- temp$formula
        }
        startatnull <- FALSE
    }
    else {
        if(!is.null(formula))
            warning("formula ignored if qtl is not provided.")
        startatnull <- TRUE
    }

    # revise names in qtl object
    if(!startatnull)
        qtl$name <- qtl$altname

    # check that we have the right stuff for the selected method
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
    if(method=="imp") qtlmethod <- "draws"
    else qtlmethod <- "prob"

    if(!is.null(qtl) && qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if(method=="imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
        else
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="prob")
    }

    if(!is.null(qtl) && method=="imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3])  {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what="draws")
    }

    # check that qtl object matches the method
    if(!startatnull) {
        if(method=="imp" && !("geno" %in% names(qtl)))
            stop("The qtl object doesn't contain imputations; re-run makeqtl with what=\"draws\".")
        else if(method=="hk" && !("prob" %in% names(qtl)))
            stop("The qtl object doesn't contain QTL genotype probabilities; re-run makeqtl with what=\"prob\".")
    }

    # check phenotypes and covariates; drop ind'ls with missing values
    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("stepwiseqtl can take just one phenotype; only the first will be used")
    }
    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }
    if(any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")
    pheno <- cross$pheno[,pheno.col]
    if(!is.null(covar)) phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors=TRUE)
    hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
    if(all(hasmissing))
        stop("All individuals are missing phenotypes or covariates.")
    if(any(hasmissing)) {
        pheno <- pheno[!hasmissing]
        cross <- subset(cross, ind=!hasmissing)
        if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]

        if(!startatnull) {
            if(method=="imp")
                qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
            else {
                for(i in seq(along=qtl$prob))
                    qtl$prob[[i]] <- qtl$prob[[i]][!hasmissing,,drop=FALSE]
            }
            qtl$n.ind <- sum(!hasmissing)
        }
    }

    if(max.qtl < 1)
        stop("Need max.qtl > 0 if we are to scan for qtl")

    if(is.null(covar)) {
        lod0 <- 0
        if(startatnull)
            firstformula <- y~Q1
        else firstformula <- formula
    }
    else {
        lod0 <- length(pheno)/2 * log10(sum((pheno-mean(pheno))^2) / sum(lm(pheno ~ as.matrix(covar))$resid^2))
        if(startatnull)
            firstformula <- as.formula(paste("y~", paste(names(covar), collapse="+"), "+", "Q1"))
        else firstformula <- formula
    }


    cross.type <- class(cross)[1]
    if(verbose > 2) verbose.scan <- TRUE
    else verbose.scan <- FALSE

    curbest <- NULL
    curbestplod <- 0

    # initial scan : either 1d or 2d
    if(verbose) cat(" -Initial scan\n")
    if(startatnull) { # Quoc: changed

        if(forceXcovar) {
            if(is.null(covar)) covar.w.X <- Xcovar
            else covar.w.X <- cbind(covar, Xcovar)
        }  else covar.w.X <- covar

        suppressWarnings(out.scanone <- scanone(cross, pheno.col=pheno.col, method=method, model=model,
                                                addcovar=covar.w.X))
        out.A <- subset(out.scanone, chr='-X')
        out.X <- subset(out.scanone, chr='X')
        formula <- firstformula
        n.qtl <- 1
        # Calculate max plod in Autosome
        lod.A <- max(out.A[,3], na.rm=TRUE)

        wh <- which(!is.na(out.A[,3]) & out.A[,3]==lod.A)
        if(length(wh) > 1) wh <- sample(wh, 1)

        qtl.A <- makeqtl(cross, as.character(out.A[wh,1]), out.A[wh,2], "Q1",
                         what=qtlmethod)

        # update Autosome plod after choosing formula and qtl
        curplod.A <- calc.plod.X(lod.A, formula=firstformula, qtl=qtl.A, penalties=penalties)
        # Calculate max plod in X chr
        lod.X <- max(out.X[,3], na.rm=TRUE)

        wh <- which(!is.na(out.X[,3]) & out.X[,3]==lod.X)
        if(length(wh) > 1) wh <- sample(wh, 1)

        qtl.X <- makeqtl(cross, as.character(out.X[wh,1]), out.X[wh,2], "Q1",
                         what=qtlmethod)

        # update X chr plod after choosing formula and qtl
        curplod.X <- calc.plod.X(lod.X, formula=firstformula, qtl=qtl.X, penalties=penalties)
        # Compare and choose plod between Autosome and X Chr
        if (curplod.X > curplod.A) {
            curplod <- curplod.X
            qtl <- qtl.X
            lod <- lod.X
        } else {
            curplod <- curplod.A
            qtl <- qtl.A
            lod <- lod.A
        }
        if(verbose) cat("initial lod: ", lod, "\n")
    } # start at null, Quoc: changed
    else {
        if(verbose) cat(" ---Starting at a model with", length(qtl$chr), "QTL\n")
        if(refine.locations) {
            if(verbose) cat(" ---Refining positions\n")
            rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                              covar=covar, formula=formula, method=method,
                              verbose=verbose.scan, incl.markers=incl.markers,
                              keeplodprofile=FALSE, forceXcovar=forceXcovar)
            if(any(rqtl$pos != qtl$pos)) { # updated positions
                if(verbose) cat(" ---  Moved a bit\n")
            }
            qtl <- rqtl
        }
        fit <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                      method=method, model=model, dropone=FALSE, get.ests=FALSE,
                      run.checks=FALSE, tol=tol, maxit=maxit, forceXcovar=forceXcovar)
        lod <- fit$result.full[1,4] - lod0
        if(require.fullrank && attr(fit, "matrix.rank") < attr(fit, "matrix.ncol")) lod <- 0
        curplod <- calc.plod.X(lod, formula=formula, qtl=qtl,
                               penalties=penalties) # Quoc changed
        attr(qtl, "pLOD") <- curplod
        n.qtl <- length(qtl$chr)
    }

    attr(qtl, "formula") <- deparseQTLformula(formula)
    attr(qtl, "pLOD") <- curplod

    if(curplod > 0) {
        curbest <- qtl
        curbestplod <- curplod

        if(verbose)
            cat("** new best ** (pLOD increased by ", round(curplod, 4), ")\n", sep="")
    }

    if(keeptrace) {
        temp <- list(chr=qtl$chr, pos=qtl$pos)
        attr(temp, "formula") <- deparseQTLformula(formula)
        attr(temp, "pLOD") <- curplod
        attr(temp, "LOD") <- lod # check LOD
        class(temp) <- c("compactqtl", "list")
        thetrace <- list("0"=temp)
    }

    if(verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
            deparseQTLformula(formula), "\n")
    if(verbose > 1)
        cat("         qtl:", paste(qtl$chr, round(qtl$pos,1), sep="@"), "\n")

    # start stepwise search
    i <- 0
    while(n.qtl < max.qtl) {
        i <- i+1

        if(verbose) {
            cat(" -Step", i, "\n")
            cat(" ---Scanning for additive qtl\n")
        }
        curformula <- as.formula(paste(deparseQTLformula(formula), "+Q", n.qtl+1, sep=""))
        # Add QTL in X chr
        out.X <- addqtl(cross, chr="X", pheno.col=pheno.col, qtl=qtl, covar=covar,
                        formula=formula, method=method, incl.markers=incl.markers,
                        verbose=verbose.scan, forceXcovar=forceXcovar,
                        require.fullrank=require.fullrank)

        curlod.X <- max(out.X[,3], na.rm=TRUE)
        wh <- which(!is.na(out.X[,3]) & out.X[,3]==curlod.X)
        if(length(wh) > 1) wh <- sample(wh,1)
        qtl.X <- addtoqtl(cross, qtl, as.character(out.X[wh,1]), out.X[wh,2],
                          paste("Q", n.qtl+1, sep=""))
        curlod.X <- curlod.X+lod
        plod.X <- calc.plod.X(curlod.X, formula=curformula, qtl=qtl.X,
                              penalties=penalties)
        # Add QTL in Autosome chr
        out.A <- addqtl(cross, chr="-X", pheno.col=pheno.col, qtl=qtl, covar=covar,
                        formula=formula, method=method, incl.markers=incl.markers,
                        verbose=verbose.scan, forceXcovar=forceXcovar,
                        require.fullrank=require.fullrank)

        curlod.A <- max(out.A[,3], na.rm=TRUE)
        wh <- which(!is.na(out.A[,3]) & out.A[,3]==curlod.A)
        if(length(wh) > 1) wh <- sample(wh,1)
        qtl.A <- addtoqtl(cross, qtl, as.character(out.A[wh,1]), out.A[wh,2],
                          paste("Q", n.qtl+1, sep=""))
        curlod.A <-curlod.A+lod
        plod.A <- calc.plod.X(curlod.A, formula=curformula, qtl=qtl.A,
                              penalties=penalties)
        if (plod.X > plod.A) {
            curplod <- plod.X
            curqtl <- qtl.X
            curlod <- curlod.X
        } else {
            curplod <- plod.A
            curqtl <- qtl.A
            curlod <- curlod.A
        }
        if(verbose) cat("        plod =", curplod, "\n")

        curnqtl <- n.qtl+1
        if(!additive.only) {
            for(j in 1:n.qtl) {

                if(verbose)
                    cat(" ---Scanning for QTL interacting with Q", j, "\n", sep="")

                thisformula <- as.formula(paste(deparseQTLformula(formula), "+Q", n.qtl+1,
                                                "+Q", j, ":Q", n.qtl+1, sep=""))
                # Scanning for QTL interacting with Qj in Autosome
                if(verbose)
                    cat(" ---Scanning for QTL in Autosome interacting with Q", j, "\n", sep="")

                out <- addqtl(cross, chr="-X", pheno.col=pheno.col, qtl=qtl, covar=covar,
                              formula=thisformula, method=method, incl.markers=incl.markers,
                              verbose=verbose.scan, forceXcovar=forceXcovar,
                              require.fullrank=require.fullrank)
                thislod <- max(out[,3], na.rm=TRUE)

                wh <- which(!is.na(out[,3]) & out[,3]==thislod)
                if(length(wh) > 1) wh <- sample(wh,1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,1]), out[wh,2],
                                    paste("Q", n.qtl+1, sep=""))

                thislod <- thislod + lod
                thisplod <- calc.plod.X(thislod, formula=thisformula, qtl=thisqtl,
                                        penalties=penalties)
                if(verbose) cat("        plod =", thisplod, "\n")

                if(thisplod > curplod) {
                    curformula <- thisformula
                    curplod <- thisplod
                    curlod <- thislod
                    curqtl <- thisqtl

                    curnqtl <- n.qtl+1
                }
                # End: Scanning for QTL interacting with Qj in Autosome

                # Scanning for QTL interacting with Qj in X chromosome
                if(verbose)
                    cat(" ---Scanning for QTL in X chromosome interacting with Q", j, "\n", sep="")

                out <- addqtl(cross, chr="X", pheno.col=pheno.col, qtl=qtl, covar=covar,
                              formula=thisformula, method=method, incl.markers=incl.markers,
                              verbose=verbose.scan, forceXcovar=forceXcovar,
                              require.fullrank=require.fullrank)
                thislod <- max(out[,3], na.rm=TRUE)

                wh <- which(!is.na(out[,3]) & out[,3]==thislod)
                if(length(wh) > 1) wh <- sample(wh,1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,1]), out[wh,2],
                                    paste("Q", n.qtl+1, sep=""))

                thislod <- thislod + lod
                thisplod <- calc.plod.X(thislod, formula=thisformula, qtl=thisqtl,
                                        penalties=penalties)
                if(verbose) cat("        plod =", thisplod, "\n")

                if(thisplod > curplod) {
                    curformula <- thisformula
                    curplod <- thisplod
                    curlod <- thislod
                    curqtl <- thisqtl

                    curnqtl <- n.qtl+1
                }
                # End: Scanning for QTL interacting with Qj in X chromosome
            }

            if(n.qtl > 1) {
                if(verbose)
                    cat(" ---Look for additional interactions\n")
                temp <- addint(cross, pheno.col, qtl, covar=covar, formula=formula,
                               method=method, qtl.only=TRUE, verbose=verbose.scan,
                               require.fullrank=require.fullrank)
                if(!is.null(temp)) {

                    formula.addint <- paste(deparseQTLformula(formula), "+", rownames(temp)) # formula of addint model
                    thelod <- temp[,3]
                    thelod <- thelod + lod
                    #             browser()
                    plod.addint <- mapply(FUN=calc.plod.X, lod=thelod, formula=formula.addint, MoreArgs=list(qtl=qtl,penalties=penalties))
                    thisplod <- max(plod.addint, na.rm=TRUE)
                    wh <- which(!is.na(thelod) & plod.addint==thisplod)
                    if(length(wh) > 1) wh <- sample(wh, 1)
                    thisformula <- as.formula(formula.addint[wh])
                    thislod <- thelod[wh]

                    if(verbose) cat("        plod =", thisplod, "\n")
                    if(thisplod > curplod) {
                        curformula <- thisformula
                        curplod <- thisplod
                        curlod <- thislod
                        curqtl <- qtl
                        curnqtl <- n.qtl
                    }
                }
            }
        }

        qtl <- curqtl
        n.qtl <- curnqtl
        attr(qtl, "formula") <- deparseQTLformula(curformula)
        attr(qtl, "pLOD") <- curplod
        formula <- curformula
        lod <- curlod

        if(refine.locations) {
            if(verbose) cat(" ---Refining positions\n")
            rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                              covar=covar, formula=formula, method=method,
                              verbose=verbose.scan, incl.markers=incl.markers,
                              keeplodprofile=FALSE, forceXcovar=forceXcovar)
            if(any(rqtl$pos != qtl$pos)) { # updated positions
                if(verbose) cat(" ---  Moved a bit\n")
                qtl <- rqtl
                fit <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                              method=method, model=model, dropone=FALSE, get.ests=FALSE,
                              run.checks=FALSE, tol=tol, maxit=maxit, forceXcovar=forceXcovar)
                lod <- fit$result.full[1,4] - lod0
                if(require.fullrank && attr(fit, "matrix.rank") < attr(fit, "matrix.ncol")) lod <- 0
                curplod <- calc.plod.X(lod, formula=formula, qtl=qtl,
                                       penalties=penalties)
                attr(qtl, "pLOD") <- curplod
            }

        }

        if(verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
                deparseQTLformula(formula), "\n")
        if(verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,1), sep="@"), "\n")


        if(curplod > curbestplod) {
            if(verbose)
                cat("** new best ** (pLOD increased by ", round(curplod - curbestplod, 4),
                    ")\n", sep="")

            curbest <- qtl
            curbestplod <- curplod
        }

        if(keeptrace) {
            temp <- list(chr=qtl$chr, pos=qtl$pos)
            attr(temp, "formula") <- deparseQTLformula(formula)
            attr(temp, "pLOD") <- curplod
            attr(temp, "LOD") <- lod # check LOD
            class(temp) <- c("compactqtl", "list")
            temp <- list(temp)
            names(temp) <- i
            thetrace <- c(thetrace, temp)
        }

        if(n.qtl >= max.qtl) break
        if (stop.rule==1) if(curplod < -k_f*penalties[1]) break # Add 1 Early stoping rule; plod<-k_f*mainA; k_f default at 3.
        if (stop.rule==2) if(curplod < (curbestplod-k_f*penalties[1])) break # Add 1 Early stoping rule; plod<bestplod-k_f*mainA k_f default at 3.
        # stop.rule==0 mean no early stop.
    }

    if(verbose) cat(" -Starting backward deletion\n")

    while(n.qtl > 1) {
        i <- i+1
        out <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                      method=method, model=model, dropone=TRUE, get.ests=FALSE,
                      run.checks=FALSE, tol=tol, maxit=maxit, forceXcovar=forceXcovar)$result.drop

        formulas <- attr(out, "formulas") # Extract formula of dropone model
        #       lods <- attr(out, "lods")

        rn <- rownames(out)
        # ignore things with covariates
        wh <- c(grep("^[Qq][0-9]+$", rn),
                grep("^[Qq][0-9]+:[Qq][0-9]+$", rn))
        out <- out[wh,,drop=FALSE]
        formulas <- formulas[wh,drop=FALSE]
        #       lods <- lods[wh,drop=FALSE]
        thelod <- out[,3]
        plod.dropone <- mapply(FUN=calc.plod.X, lod=lod-thelod, formula=formulas, MoreArgs=list(qtl=qtl,penalties=penalties))
        maxplod <- max(plod.dropone, na.rm=TRUE)

        wh <- which(!is.na(thelod) & plod.dropone==maxplod)
        if(length(wh) > 1) wh <- sample(wh, 1)

        lod <- lod - thelod[wh] # End of Quoc change, below code is unchanged.
        todrop <- rownames(out)[wh]

        if(verbose) cat(" ---Dropping", todrop, "\n")

        if(length(grep(":", todrop)) > 0) { # dropping an interaction
            theterms <- attr(terms(formula), "factors")
            wh <- colnames(theterms)==todrop
            if(!any(wh)) stop("Confusion about what interation to drop!")
            theterms <- colnames(theterms)[!wh]
            formula <- as.formula(paste("y~", paste(theterms, collapse="+")))
        }
        else {
            numtodrop <- as.numeric(substr(todrop, 2, nchar(todrop)))

            theterms <- attr(terms(formula), "factors")
            cn <- colnames(theterms)
            g <- c(grep(paste("^[Qq]", numtodrop, "$", sep=""), cn),
                   grep(paste("^[Qq]", numtodrop, ":", sep=""), cn),
                   grep(paste(":[Qq]", numtodrop, "$", sep=""), cn))
            cn <- cn[-g]
            formula <- as.formula(paste("y~", paste(cn, collapse="+")))

            if(n.qtl > numtodrop) {
                for(j in (numtodrop+1):n.qtl)
                    formula <- reviseqtlnuminformula(formula, j, j-1)
            }

            qtl <- dropfromqtl(qtl, index=numtodrop)
            qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl, sep="")
            n.qtl <- n.qtl - 1
        }

        curplod <- calc.plod.X(lod, formula=formula, qtl=qtl,
                               penalties=penalties)
        #       browser()
        if(verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
                deparseQTLformula(formula), "\n")
        if(verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,1), sep=":"), "\n")

        attr(qtl, "formula") <- deparseQTLformula(formula)
        attr(qtl, "pLOD") <- curplod

        if(refine.locations) {
            if(verbose) cat(" ---Refining positions\n")
            if(!is.null(qtl)) {
                rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                                  covar=covar, formula=formula, method=method,
                                  verbose=verbose.scan, incl.markers=incl.markers,
                                  keeplodprofile=FALSE, forceXcovar=forceXcovar)
                if(any(rqtl$pos != qtl$pos)) { # updated positions
                    if(verbose) cat(" ---  Moved a bit\n")
                    qtl <- rqtl
                    fit <- fitqtl(cross, pheno.col, qtl, covar=covar, formula=formula,
                                  method=method, model=model, dropone=FALSE, get.ests=FALSE,
                                  run.checks=FALSE, tol=tol, maxit=maxit, forceXcovar=forceXcovar)
                    lod <- fit$result.full[1,4] - lod0
                    if(require.fullrank && attr(fit, "matrix.rank") < attr(fit, "matrix.ncol")) lod <- 0
                    curplod <- calc.plod.X(lod, formula=formula, qtl=qtl,
                                           penalties=penalties)
                    attr(qtl, "pLOD") <- curplod
                }
            }
        }

        if(curplod > curbestplod) {
            if(verbose)
                cat("** new best ** (pLOD increased by ", round(curplod - curbestplod, 4),
                    ")\n", sep="")
            #         browser()
            curbestplod <- curplod
            curbest <- qtl
        }

        if(keeptrace) {
            temp <- list(chr=qtl$chr, pos=qtl$pos)
            attr(temp, "formula") <- deparseQTLformula(formula)
            attr(temp, "pLOD") <- curplod
            attr(temp, "LOD") <- lod # check LOD
            class(temp) <- c("compactqtl", "list")
            temp <- list(temp)
            names(temp) <- i
            thetrace <- c(thetrace, temp)
        }
    }

    # re-form the qtl
    if(!is.null(curbest)) {
        chr <- curbest$chr
        pos <- curbest$pos
        o <- order(factor(chr, levels=names(cross$geno)), pos)

        qtl <- makeqtl(cross, chr[o], pos[o], what=qtlmethod)

        # need to redo numbering in formula
        formula <- as.formula(attr(curbest, "formula"))

        if(length(chr) > 1) {
            n.qtl <- length(chr)
            for(i in 1:n.qtl)
                formula <- reviseqtlnuminformula(formula, i, n.qtl+i)
            for(i in 1:n.qtl)
                formula <- reviseqtlnuminformula(formula, n.qtl+o[i], i)
        }

        if(keeplodprofile) {
            if(verbose) cat(" ---One last pass through refineqtl\n")
            qtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl,
                             covar=covar, formula=formula, method=method,
                             verbose=verbose.scan, incl.markers=incl.markers,
                             keeplodprofile=TRUE, forceXcovar=forceXcovar)
        }
        attr(qtl, "formula") <- deparseQTLformula(formula)
        attr(qtl, "pLOD") <- attr(curbest, "pLOD")
        curbest <- qtl
    }
    else {
        curbest <- numeric(0)
        class(curbest) <- "qtl"
        attr(curbest,"pLOD") <- 0
    }

    if(keeptrace)
        attr(curbest, "trace") <- thetrace

    attr(curbest, "formula") <- deparseQTLformula(attr(curbest, "formula"), TRUE)
    attr(curbest, "penalties") <- penalties

    curbest
}


######################################################################
# penalized LOD score, Quoc changed, now central function for pLOD
######################################################################
calc.plod.X <-
    function(lod, nterms, type=c("f2","bc"), penalties, formula, qtl)
{
    if (missing(nterms) && (!missing(formula) && !is.null(formula)) && !missing(qtl))
        nterms <- countqtltermsX(formula, qtl)
    nterms <- nterms[1:6]
    if(any(penalties==Inf & nterms > 0)) return(-Inf)

    as.numeric(lod - sum((nterms*penalties)[nterms > 0]))
}

# end of stepwiseqtlX.R
