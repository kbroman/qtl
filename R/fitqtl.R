######################################################################
#
# fitqtl.R
#
# copyright (c) 2002-2013, Hao Wu and Karl W. Broman
# last modified Sep, 2013
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
# Contains: fitqtl, fitqtlengine, parseformula, summary.fitqtl,
#           print.summary.fitqtl, mybinaryrep, deparseQTLformula
#           printQTLformulanicely
#
######################################################################

######################################################################
#
# This is the function to fit a model and generate some tables
#
# Only Haley-Knott regression and multiple imputation are implemented.
#
######################################################################

fitqtl <-
    function(cross, pheno.col=1, qtl, covar=NULL, formula, method=c("imp", "hk"),
             model=c("normal", "binary"), dropone=TRUE, get.ests=FALSE,
             run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
{
    # some input checking stuff in here
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
        warning("fitqtl can take just one phenotype; only the first will be used")
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

    if(model=="binary" && any(!is.na(pheno) & pheno != 0 & pheno != 1))
        stop("For model=\"binary\", phenotypes must by 0 or 1.")

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

    fitqtlengine(pheno=pheno, qtl=qtl, covar=covar, formula=formula,
                 method=method, model=model, dropone=dropone, get.ests=get.ests,
                 run.checks=run.checks, cross.attr=attributes(cross),
                 sexpgm=getsex(cross), tol=tol, maxit=maxit, forceXcovar=forceXcovar)
}


fitqtlengine <-
    function(pheno, qtl, covar=NULL, formula, method=c("imp", "hk"),
             model=c("normal", "binary"),
             dropone=TRUE, get.ests=FALSE, run.checks=TRUE, cross.attr,
             sexpgm, tol, maxit, forceXcovar=FALSE)
{
    model <- match.arg(model)
    method <- match.arg(method)

    # local variables
    n.ind <- qtl$n.ind # number of individuals
    n.qtl <- qtl$n.qtl # number of selected markers
    n.gen <- qtl$n.gen # number of genotypes
    if(method=="imp")
        n.draws <- dim(qtl$geno)[3] # number of draws

    if( is.null(covar) )  # number of covariates
        n.covar <- 0
    else
        n.covar <- ncol(covar)

    # if formula is missing, build one
    # all QTLs and covariates will be additive by default
    if(missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep="") # QTL term names
        formula <- "y~Q1"
        if(n.qtl > 1)
            for (i in 2:n.qtl)
                formula <- paste(formula, tmp.Q[i], sep="+")
        if (n.covar) { # if covariate is not empty
            tmp.C <- colnames(covar) # covariate term names
            for(i in 1:n.covar)
                formula <- paste(formula, tmp.C[i], sep="+")
        }
        formula <- as.formula(formula)
    }

    # check input formula
    if(run.checks) {
        formula <- checkformula(formula, qtl$altname, colnames(covar))

        if(qtl$n.ind != length(pheno))
            stop("No. individuals in qtl object doesn't match length of input phenotype.")

        # drop covariates that are not in the formula
        if(!is.null(covar)) {
            theterms <- rownames(attr(terms(formula), "factors"))
            m <- match(colnames(covar), theterms)
            if(all(is.na(m))) covar <- NULL
            else covar <- covar[,!is.na(m),drop=FALSE]
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
                n.ind <- qtl$n.ind <- sum(!hasmissing)
                if(method=="imp")
                    qtl$geno <- qtl$geno[!hasmissing,,,drop=FALSE]
                else
                    qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,,drop=FALSE])

                if(!is.null(covar)) covar <- covar[!hasmissing,,drop=FALSE]

                # subset sexpgm
                for(i in seq(along=sexpgm))
                    if(!is.null(sexpgm[[i]]))
                        sexpgm[[i]] <- sexpgm[[i]][!hasmissing]
            }
        }
    }

    # parse the input formula
    p <- parseformula(formula, qtl$altname, colnames(covar))
    # make an array n.gen.QC to represent the genotype numbers
    # for all input QTLs and covariates. For covariates the
    # number of genotyps is 1. This makes programming easier
    n.gen.QC <- c(n.gen[p$idx.qtl]-1, rep(1, p$n.covar))

    # covariates to be passed to C function
    # This is done in case of that user input covar but has no covar in formula
    covar.C <- NULL
    if(!is.null(p$idx.covar))
        covar.C <- as.matrix(covar[,p$idx.covar,drop=FALSE])

    sizefull <- 1+sum(n.gen.QC)
    if(p$n.int > 0) {
        form <- p$formula.intmtx*n.gen.QC
        if(!is.matrix(form)) {
            sizefull <- sizefull + prod(form[form!=0])
        }
        else {
            form <- apply(form,2,function(a) prod(a[a != 0]))
            sizefull <- sizefull + sum(form)
        }
    }

    if(method != "imp") {
        # form genotype probabilities as a matrix
        prob <- matrix(ncol=sum(qtl$n.gen[p$idx.qtl]), nrow=n.ind)

        curcol <- 0
        for(i in p$idx.qtl) {
            prob[,curcol+1:n.gen[i]] <- qtl$prob[[i]]
            curcol <- curcol + n.gen[i]
        }
    }

    Xadjustment <- scanoneXnull(cross.attr$class[1], sexpgm, cross.attr)
    n.origcovar <- p$n.covar
    if((sum(qtl$chrtype[p$idx.qtl]=="X") >= 1 || forceXcovar) && Xadjustment$adjustX) { # need to include X chromosome covariates
        adjustX <- TRUE

        n.newcovar <- ncol(Xadjustment$sexpgmcovar)
        n.gen.QC <- c(n.gen.QC, rep(1, n.newcovar))
        p$n.covar <- p$n.covar + n.newcovar
        covar.C <- cbind(covar.C, Xadjustment$sexpgmcovar)
        sizefull <- sizefull + n.newcovar

        if(p$n.int==1)
            p$formula.intmtx <- c(p$formula.intmtx, rep(0,n.newcovar))
        if(p$n.int>1) {
            for(i in 1:n.newcovar)
                p$formula.intmtx <- rbind(p$formula.intmtx, rep(0,p$n.int))
        }
    }
    else adjustX <- FALSE

    # call C function to do the genome scan
    if(model=="normal") {
        if(method == "imp") {
            z <- .C("R_fitqtl_imp",
                    as.integer(n.ind), # number of individuals
                    as.integer(p$n.qtl), # number of qtls
                    as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                    as.integer(n.draws), # number of draws
                    as.integer(qtl$geno[,p$idx.qtl,]), # genotypes for selected marker
                    as.integer(p$n.covar), # number of covariate
                    as.double(covar.C), # covariate
                    as.integer(p$formula.intmtx),  # formula matrix for interactive terms
                    as.integer(p$n.int), # number of interactions in the formula
                    as.double(pheno), # phenotype
                    as.integer(get.ests), # get estimates?
                    # return variables
                    lod=as.double(0), # LOD score
                    df=as.integer(0), # degree of freedom
                    ests=as.double(rep(0,sizefull)),
                    ests.cov=as.double(rep(0,sizefull*sizefull)),
                    design.mat=as.double(rep(0,sizefull*n.ind)),
                    matrix.rank=as.integer(0), # on return, minimum of matrix rank across imputations
                    resid=as.double(rep(0,n.ind)), # on return, the average residuals across imputations
                    PACKAGE="qtl")
        }
        else {
            z <- .C("R_fitqtl_hk",
                    as.integer(n.ind), # number of individuals
                    as.integer(p$n.qtl), # number of qtls
                    as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                    as.double(prob),      # QTL genotype probabilities
                    as.integer(p$n.covar), # number of covariate
                    as.double(covar.C), # covariates
                    as.integer(p$formula.intmtx),  # formula matrix for interactive terms
                    as.integer(p$n.int), # number of interactions in the formula
                    as.double(pheno), # phenotype
                    as.integer(get.ests), # get estimates?
                    # return variables
                    lod=as.double(0), # LOD score
                    df=as.integer(0), # degree of freedom
                    ests=as.double(rep(0,sizefull)),
                    ests.cov=as.double(rep(0,sizefull*sizefull)),
                    design.mat=as.double(rep(0,sizefull*n.ind)),
                    matrix.rank=as.integer(0), # on return, rank of matrix
                    resid=as.double(rep(0,n.ind)), # on return, residuals from the fit
                    PACKAGE="qtl")
        }
    }
    else {
        if(method=="imp") {
            z <- .C("R_fitqtl_imp_binary",
                    as.integer(n.ind), # number of individuals
                    as.integer(p$n.qtl), # number of qtls
                    as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                    as.integer(n.draws), # number of draws
                    as.integer(qtl$geno[,p$idx.qtl,]), # genotypes for selected marker
                    as.integer(p$n.covar), # number of covariate
                    as.double(covar.C), # covariate
                    as.integer(p$formula.intmtx),  # formula matrix for interactive terms
                    as.integer(p$n.int), # number of interactions in the formula
                    as.double(pheno), # phenotype
                    as.integer(get.ests), # get estimates?
                    # return variables
                    lod=as.double(0), # LOD score
                    df=as.integer(0), # degree of freedom
                    ests=as.double(rep(0,sizefull)),
                    ests.cov=as.double(rep(0,sizefull*sizefull)),
                    design.mat=as.double(rep(0,sizefull*n.ind)),
                    as.double(tol),
                    as.integer(maxit),
                    matrix.rank=as.integer(0), # on return, minimum of matrix rank across imputations
                    PACKAGE="qtl")
        }
        else {
            z <- .C("R_fitqtl_hk_binary",
                    as.integer(n.ind), # number of individuals
                    as.integer(p$n.qtl), # number of qtls
                    as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                    as.double(prob),      # QTL genotype probabilities
                    as.integer(p$n.covar), # number of covariate
                    as.double(covar.C), # covariates
                    as.integer(p$formula.intmtx),  # formula matrix for interactive terms
                    as.integer(p$n.int), # number of interactions in the formula
                    as.double(pheno), # phenotype
                    as.integer(get.ests), # get estimates?
                    # return variables
                    lod=as.double(0), # LOD score
                    df=as.integer(0), # degree of freedom
                    ests=as.double(rep(0,sizefull)),
                    ests.cov=as.double(rep(0,sizefull*sizefull)),
                    design.mat=as.double(rep(0,sizefull*n.ind)),
                    # convergence
                    as.double(tol),
                    as.integer(maxit),
                    matrix.rank=as.integer(0), # on return, rank of matrix
                    PACKAGE="qtl")
        }

    }
    matrix.rank <- z$matrix.rank
    matrix.ncol <- sizefull
    residuals <- z$resid

    if(get.ests) {
        # first, construct the new design matrix
        #  X = the matrix used in the C coe
        #  Z = the matrix we want
        thenames <- qtl$name[p$idx.qtl]
        if(n.covar > 0) thenames <- c(thenames,names(covar)[p$idx.covar])

        ests <- z$ests
        ests.cov <- matrix(z$ests.cov,ncol=sizefull)

        if(adjustX) {
            keep <- 1:(length(ests)-n.newcovar)
            ests <- ests[keep]
            ests.cov <- ests.cov[keep,keep]
        }

        if(any(qtl$n.gen[p$idx.qtl]>=4)) {
            type <- cross.attr$class[1]
            if(type == "4way")
                genotypes <- c("AC","BC","AD","BD")
            else {
                genotypes <- cross.attr$genotypes
                if(is.null(genotypes))
                    genotypes <- as.character(1:max(qtl$n.gen))
            }

            # just attach rownames for this case
            thenames <- "Intercept"

            if(length(p$idx.qtl) > 0) { # qtl names
                qtlnames <- vector("list", length(p$idx.qtl))
                names(qtlnames) <- paste("Q", 1:length(p$idx.qtl), sep="")
                for(i in seq(along=p$idx.qtl)) {
                    qtlnames[[i]] <- paste(qtl$name[p$idx.qtl[i]], genotypes[2:qtl$n.gen[p$idx.qtl[i]]], sep=".")
                    thenames <- c(thenames, qtlnames[[i]])
                }
            }

            if(p$n.covar > 0) {
                covnames <- colnames(covar)[p$idx.covar]
                thenames <- c(thenames, covnames)
            }

            if(p$n.int > 0) { # interactions
                if(!is.matrix(p$formula.intmtx)) {
                    nam <- names(p$formula.intmtx)
                    p$formula.intmtx <- matrix(p$formula.intmtx, nrow=p$n.int)
                    colnames(p$formula.intmtx) <- nam
                }

                for(i in 1:p$n.int) {
                    wh <- which(p$formula.intmtx[i,]==1)
                    nam <- colnames(p$formula.intmtx)[wh]

                    if(length(grep("^Q[0-9]+$", nam[1])) > 0)
                        curnam <- qtlnames[[nam[1]]]
                    else
                        curnam <- nam[1]

                    for(j in 2:length(nam)) {
                        if(length(grep("Q[0-9]+", nam[j])) > 0)
                            curnam <- paste(curnam, qtlnames[[nam[j]]], sep=".")
                        else
                            curnam <- paste(curnam, nam[j], sep=".")
                    }
                    thenames <- c(thenames, curnam)
                }
            }

            if(length(thenames) ==length(ests)) {
                names(ests) <- thenames
                dimnames(ests.cov) <- list(thenames, thenames)
            }
            else
                warning("Estimated QTL effects not yet made meaningful for this case.\n   ")

        }
        else {
            X <- matrix(z$design.mat,ncol=sizefull)
            Z <- matrix(0,nrow=n.ind,ncol=sizefull)
            colnames(Z) <- rep("",sizefull)

            # mean column
            Z[,1] <- 1
            colnames(Z)[1] <- "Intercept"

            # ZZ stores the main effects matrices, for creating the interactions
            ZZ <- vector("list",p$n.qtl+p$n.covar)
            curcol <- 1
            # covariates
            if(p$n.covar > 0) {
                for(j in 1:p$n.covar) {
                    Z[,curcol+j] <- ZZ[[p$n.qtl+j]] <- as.matrix(covar[,p$idx.covar[j],drop=FALSE])
                    colnames(Z)[curcol+j] <- colnames(ZZ[[p$n.qtl+j]]) <- names(covar)[p$idx.covar[j]]
                }
                curcol <- curcol + p$n.covar
            }
            # QTL main effects
            for(i in seq(along=p$idx.qtl)) {
                if(n.gen[i]==2) {
                    if(method=="imp") {
                        if(cross.attr$class[1] == "bc") {
                            Z[qtl$geno[,p$idx.qtl[i],1]==1,curcol+1] <- -0.5
                            Z[qtl$geno[,p$idx.qtl[i],1]==2,curcol+1] <- 0.5
                        } else {
                            Z[qtl$geno[,p$idx.qtl[i],1]==1,curcol+1] <- -1
                            Z[qtl$geno[,p$idx.qtl[i],1]==2,curcol+1] <- 1
                        }
                    }
                    else
                        if(cross.attr$class[1] == "bc") {
                            Z[,curcol+1] <- (qtl$prob[[p$idx.qtl[i]]][,2] - qtl$prob[[p$idx.qtl[i]]][,1])/2
                        } else {
                            Z[,curcol+1] <- (qtl$prob[[p$idx.qtl[i]]][,2] - qtl$prob[[p$idx.qtl[i]]][,1])
                        }
                    colnames(Z)[curcol+1] <- thenames[i]
                }
                else { # 3 genotypes
                    if(method=="imp") {
                        Z[qtl$geno[,p$idx.qtl[i],1]==1,curcol+1] <- -1
                        Z[qtl$geno[,p$idx.qtl[i],1]==3,curcol+1] <- 1
                        Z[qtl$geno[,p$idx.qtl[i],1]==2,curcol+2] <- 0.5
                        Z[qtl$geno[,p$idx.qtl[i],1]!=2,curcol+2] <- -0.5
                    }
                    else {
                        Z[,curcol+1] <- qtl$prob[[p$idx.qtl[i]]][,3] - qtl$prob[[p$idx.qtl[i]]][,1]
                        Z[,curcol+2] <- (qtl$prob[[p$idx.qtl[i]]][,2] - qtl$prob[[p$idx.qtl[i]]][,1] -
                                         qtl$prob[[p$idx.qtl[i]]][,3])/2
                    }
                    colnames(Z)[curcol+1:2] <- paste(thenames[i],c("a","d"),sep="")
                }
                ZZ[[i]] <- Z[,curcol+1:(n.gen[i]-1),drop=FALSE]
                curcol <- curcol + n.gen[i]-1
            }

            if(p$n.int>0) {
                for(i in 1:p$n.int) {
                    if(p$n.int==1)
                        intform <- p$formula.intmtx
                    else
                        intform <- p$formula.intmtx[,i]
                    tempZ <- matrix(1,ncol=1,nrow=nrow(Z))
                    colnames(tempZ) <- ""
                    for(j in seq(along=intform)) {
                        if(intform[j]==1) {
                            tZ <- NULL
                            for(k in 1:ncol(ZZ[[j]])) {
                                tZZ <- tempZ * ZZ[[j]][,k]
                                if(all(colnames(tempZ) == ""))
                                    colnames(tZZ) <- colnames(ZZ[[j]])[k]
                                else
                                    colnames(tZZ) <- paste(colnames(tempZ),colnames(ZZ[[j]])[k],sep=":")
                                tZ <- cbind(tZ,tZZ)
                            }
                            tempZ <- tZ
                        }
                    }
                    Z[,curcol+1:ncol(tempZ)] <- tempZ
                    colnames(Z)[curcol+1:ncol(tempZ)] <- colnames(tempZ)
                    curcol <- curcol + ncol(tempZ)
                }
            }

            b <- solve(t(Z) %*% Z, t(Z) %*% X)
            ests <- as.numeric(b %*% ests)
            ests.cov <- b %*% ests.cov %*% t(b)
            names(ests) <- colnames(Z)
            dimnames(ests.cov) <- list(colnames(Z),colnames(Z))
        }
    }

    ##### output ANOVA table for full model #####
    result.full <- matrix(NA, 3, 7)
    colnames(result.full) <- c("df", "SS", "MS", "LOD", "%var", "Pvalue(Chi2)",
                               "Pvalue(F)")
    rownames(result.full) <- c("Model", "Error", "Total")
    result.full[1,1] <- z$df # model degree of freedom

    if(model=="normal") {
        # compute the SS for total

        if(adjustX) {
            mpheno <- mean(pheno)
            Rss0 <- sum( (pheno-mpheno)^2 )

            Rss0adj <- Rss0x <- sum( lm(pheno ~ Xadjustment$sexpgmcovar)$resid^2 )
            OrigModellod <- z$lod
            Modellod <- z$lod + length(pheno)/2 * (log10(Rss0x) - log10(Rss0))
        } else {
            mpheno <- mean(pheno)
            Rss0adj <- Rss0 <- sum( (pheno-mpheno)^2 )
            OrigModellod <- Modellod <- z$lod
        }

        # third row, for Total
        result.full[3,1] <- length(pheno) - 1 # total degree of freedom
        result.full[3,2] <- Rss0adj # total sum of squares

        # first row, for Model
        result.full[1,1] <- z$df # df for Model
        # Variance explained by model
        result.full[1,5] <- 100 * (1 - exp(-2*Modellod*log(10)/n.ind))
        result.full[1,2] <- Rss0adj * result.full[1,5]/100  # SS for model
        result.full[1,3] <- result.full[1,2]/z$df # MS for model
        result.full[1,4] <- Modellod # Model LOD score

        # Second row, for Error
        # df
        result.full[2,1] <- result.full[3,1] - result.full[1,1]
        # SS
        result.full[2,2] <- result.full[3,2] - result.full[1,2]
        # MS
        result.full[2,3] <- result.full[2,2] / result.full[2,1]

        # first row, P values
        # P value (chi2) for model
        result.full[1,6] <- 1 - pchisq(2*log(10)*Modellod, z$df)
        # P value (F statistics) for model
        df0 <- result.full[3,1]; df1 <- result.full[2,1];
        Rss1 <- result.full[2,2]
        Fstat <- ((Rss0adj-Rss1)/(df0-df1)) / (Rss1/df1)
        result.full[1,7] <- 1 - pf(Fstat, df0-df1, df1)
    }
    else {
        # third row, for Total
        result.full[3,1] <- length(pheno) - 1 # total degree of freedom
        result.full[3,2] <- NA # total sum of squares

        # first row, for Model
        result.full[1,1] <- z$df # df for Model
        # Variance explained by model
        result.full[1,5] <- 100 * (1 - exp(-2*z$lod*log(10)/n.ind))
        result.full[1,2] <- NA  # SS for model
        result.full[1,3] <- NA # MS for model
        result.full[1,4] <- z$lod # Model LOD score
        OrigModellod <- Modellod <- z$lod

        # Second row, for Error
        # df
        result.full[2,1] <- result.full[3,1] - result.full[1,1]
        # SS
        result.full[2,2] <- NA
        # MS
        result.full[2,3] <- NA

        # first row, P values
        # P value (chi2) for model
        result.full[1,6] <- 1 - pchisq(2*log(10)*z$lod, z$df)
        # P value (F statistics) for model
        result.full[1,7] <- NA
    }
    ############# Finish ANOVA table for full model


    # initialize output object
    output <- NULL
    output$result.full <- result.full

    # drop one at a time?
    if(dropone && (p$n.qtl+n.origcovar)>1) {
        # user wants to do drop one term at a time and output anova table

        # get the terms etc. for input formula
        f.terms <- terms(formula)
        f.order <- attr(f.terms, "order")
        f.label <- attr(f.terms, "term.labels")

        # initialize output matrix
        # ANOVA table will have five columns, e.g., df,Type III SS,
        # LOD, %var, Pvalue for each dropping term
        # Full model result will not be in this table
        result <- matrix(0, length(f.order), 7)
        colnames(result) <- c("df", "Type III SS", "LOD", "%var", "F value",
                              "Pvalue(Chi2)", "Pvalue(F)")
        rownames(result) <- rep("",length(f.order))

        drop.term.name <- NULL
        formulas <- rep("", length(f.order))
        lods <- rep(NA, length(f.order))

        for( i in (1:length(f.order)) ) {
            # loop thru all terms in formula, from the highest order
            # the label of the term to be droped
            label.term.drop <- f.label[i]

            ### find the corresponding QTL name for this term ###
            # This is used for output ANOVA table
            if(f.order[i] == 1) {
                # this is a first order term
                # if the term label is like Q(q)1, Q(q)2, etc., then it's a QTL
                if( length(grep("Q[0-9]", label.term.drop, ignore.case=TRUE)) != 0) {
                    idx.qtlname <- as.integer(substr(label.term.drop, 2, nchar(label.term.drop)))
                    drop.term.name[i] <- qtl$name[idx.qtlname]
                }
                else { # this is a covariate
                    drop.term.name[i] <- label.term.drop
                }
            }
            else {
                # this is a 2nd (or higher)order and the term is a string like "Q2:Q3:C1"
                # I use strsplit to split it to a string vector "Q2" "Q3" "C1".
                # then take out 2 and 3 as integer. Then find out the
                # QTL name from the input QTL object and concatenate them
                tmp.str <- strsplit(label.term.drop,":")[[1]]
                for(j in 1:length(tmp.str)) {
                    if( length(grep("Q[0-9]", tmp.str[j], ignore.case=TRUE)) != 0 ) {
                        # this is a QTL
                        idx.qtlname <- as.integer(substr(tmp.str[j], 2, nchar(tmp.str[j])))
                        tmp.str[j] <- qtl$name[idx.qtlname]
                    }
                    if(j == 1) # first term
                        drop.term.name[i] <- tmp.str[j]
                    else # not the first term
                        drop.term.name[i] <- paste(drop.term.name[i], tmp.str[j], sep=":")
                }
            }
            ### Finish QTL name ###

            # find the indices of the term(s) to be dropped
            # All terms contain label.term.drop will be dropped
            idx.term.drop <- NULL
            tmp.str.drop <- tolower(strsplit(label.term.drop,":")[[1]])
            for(j in 1:length(f.label)) {
                tmp.str.label <- tolower(strsplit(f.label[j], ":")[[1]])
                if(all(tmp.str.drop %in% tmp.str.label))
                    idx.term.drop <- c(idx.term.drop, j)
            }

            # the indices of term(s) to be kept
            idx.term.kept <- setdiff(1:length(f.order), idx.term.drop)

            #### regenerate a formula with the kept terms additive ###
            if(length(idx.term.kept) == 0) # nothing left after drop label.term.drop
                stop("There will be nothing left if drop ", drop.term.name[i])
            else {
                # All terms for idx.term.kept will be additive
                formula.new <- as.formula(paste("y~", paste(f.label[idx.term.kept], collapse="+"), sep=""))
            }

            ### Start fitting model again
            # parse the input formula
            p.new <- parseformula(formula.new, qtl$altname, colnames(covar))
            n.gen.QC <- c(n.gen[p.new$idx.qtl]-1, rep(1, p.new$n.covar))

            formulas[i] <- deparseQTLformula(formula.new)

            # covariate to be passed to C function
            covar.C <- NULL
            if(!is.null(p.new$idx.covar))
                covar.C <- as.matrix(covar[,p.new$idx.covar,drop=FALSE])

            if(method != "imp") {
                # form genotype probabilities as a matrix
                prob <- matrix(ncol=sum(qtl$n.gen[p.new$idx.qtl]), nrow=n.ind)

                curcol <- 0
                for(z in p.new$idx.qtl) {
                    prob[,curcol+1:n.gen[z]] <- qtl$prob[[z]]
                    curcol <- curcol + n.gen[z]
                }
            }

            if(adjustX)  { # need to include X chromosome covariates
                n.newcovar <- ncol(Xadjustment$sexpgmcovar)
                n.gen.QC <- c(n.gen.QC, rep(1, n.newcovar))
                p.new$n.covar <- p.new$n.covar + n.newcovar
                covar.C <- cbind(covar.C, Xadjustment$sexpgmcovar)
                sizefull <- sizefull + n.newcovar

                if(p.new$n.int==1)
                    p.new$formula.intmtx <- c(p.new$formula.intmtx, rep(0,n.newcovar))
                if(p.new$n.int>1) {
                    for(i2 in 1:n.newcovar)
                        p.new$formula.intmtx <- rbind(p.new$formula.intmtx, rep(0,p.new$n.int))
                }
            }

            # call C function fit model
            if(model=="normal") {
                if(method == "imp") {
                    z <- .C("R_fitqtl_imp",
                            as.integer(n.ind), # number of individuals
                            as.integer(p.new$n.qtl), # number of qtls
                            as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                            as.integer(n.draws), # number of draws
                            as.integer(qtl$geno[,p.new$idx.qtl,]), # genotypes for selected marker
                            as.integer(p.new$n.covar), # number of covariate
                            as.double(covar.C), # covariate
                            as.integer(p.new$formula.intmtx),  # formula matrix for interactive terms
                            as.integer(p.new$n.int), # number of interactions in the formula
                            as.double(pheno), # phenotype
                            as.integer(0),
                            # return variables
                            lod=as.double(0), # LOD score
                            df=as.integer(0), # degree of freedom
                            as.double(rep(0,sizefull)),
                            as.double(rep(0,sizefull*sizefull)),
                            as.double(rep(0,n.ind*sizefull)),
                            matrix.rank=as.integer(0),
                            resid=as.double(rep(0,n.ind)), # on return, the average residuals across imputations
                            PACKAGE="qtl")
                }

                else {
                    z <- .C("R_fitqtl_hk",
                            as.integer(n.ind), # number of individuals
                            as.integer(p.new$n.qtl), # number of qtls
                            as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                            as.double(prob),
                            as.integer(p.new$n.covar), # number of covariate
                            as.double(covar.C), # covariate
                            as.integer(p.new$formula.intmtx),  # formula matrix for interactive terms
                            as.integer(p.new$n.int), # number of interactions in the formula
                            as.double(pheno), # phenotype
                            as.integer(0),
                            # return variables
                            lod=as.double(0), # LOD score
                            df=as.integer(0), # degree of freedom
                            as.double(rep(0,sizefull)),
                            as.double(rep(0,sizefull*sizefull)),
                            as.double(rep(0,n.ind*sizefull)),
                            matrix.rank=as.integer(0),
                            resid=as.double(rep(0,n.ind)), # on return, the residuals
                            PACKAGE="qtl")
                }
            }
            else { # binary trait
                if(method=="imp") {
                    z <- .C("R_fitqtl_imp_binary",
                            as.integer(n.ind), # number of individuals
                            as.integer(p.new$n.qtl), # number of qtls
                            as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                            as.integer(n.draws), # number of draws
                            as.integer(qtl$geno[,p.new$idx.qtl,]), # genotypes for selected marker
                            as.integer(p.new$n.covar), # number of covariate
                            as.double(covar.C), # covariate
                            as.integer(p.new$formula.intmtx),  # formula matrix for interactive terms
                            as.integer(p.new$n.int), # number of interactions in the formula
                            as.double(pheno), # phenotype
                            as.integer(0),
                            # return variables
                            lod=as.double(0), # LOD score
                            df=as.integer(0), # degree of freedom
                            as.double(rep(0,sizefull)),
                            as.double(rep(0,sizefull*sizefull)),
                            as.double(rep(0,n.ind*sizefull)),
                            as.double(tol),
                            as.integer(maxit),
                            matrix.rank=as.integer(0),
                            PACKAGE="qtl")
                }
                else {
                    z <- .C("R_fitqtl_hk_binary",
                            as.integer(n.ind), # number of individuals
                            as.integer(p.new$n.qtl), # number of qtls
                            as.integer(n.gen.QC), # number of genotypes QTLs and covariates
                            as.double(prob),
                            as.integer(p.new$n.covar), # number of covariate
                            as.double(covar.C), # covariate
                            as.integer(p.new$formula.intmtx),  # formula matrix for interactive terms
                            as.integer(p.new$n.int), # number of interactions in the formula
                            as.double(pheno), # phenotype
                            as.integer(0),
                            # return variables
                            lod=as.double(0), # LOD score
                            df=as.integer(0), # degree of freedom
                            as.double(rep(0,sizefull)),
                            as.double(rep(0,sizefull*sizefull)),
                            as.double(rep(0,n.ind*sizefull)),
                            # convergence
                            as.double(tol),
                            as.integer(maxit),
                            matrix.rank = as.integer(0),
                            PACKAGE="qtl")
                }
            }

            if(model=="normal" && adjustX) # adjust for X chromosome covariates
                z$lod <- z$lod + length(pheno)/2 * (log10(Rss0x) - log10(Rss0))

            # record the result for dropping this term
            # df
            result[i,1] <- result.full[1,1] - z$df
            # LOD score
            result[i,3] <- Modellod - z$lod
            # % variance explained
            result[i,4] <- result.full[1,5] - 100*(1 - 10^(-2*z$lod/n.ind))

            # lod score for reduced model
            lods[i] <- z$lod

            # Type III SS for this term - computed from %var
            if(model=="normal")
                result[i,2] <- result.full[3,2] * result[i,4] / 100
            else
                result[i,2] <- NA
            # F value
            if(model=="normal") {
                df0 <- length(pheno) - z$df - 1; df1 <- result.full[2,1];
                Rss0p <- result.full[2,2] + result[i,2];
                Rss1p <- result.full[2,2]
                Fstat <- ((Rss0p-Rss1p)/(df0-df1)) / (Rss1/df1)
                result[i,5] <- Fstat
                # P value (F)
                result[i,7] <- 1 - pf(Fstat, df0-df1, df1)
            }
            else # ignore F stat for binary trait
                result[i,c(5,7)] <- NA
            # P value (chi2)
            result[i,6] <- 1 - pchisq(2*log(10)*result[i,3], result[i,1])
            # assign row name
            rownames(result)[i] <- drop.term.name[i]
        } # finish dropping terms loop

        attr(result, "formulas") <- formulas
        attr(result, "lods") <- lods

        # assign output object
        output$result.drop <- result

    }  ## if(dropone)

    if(get.ests)
        output$ests <- list(ests=ests, covar=ests.cov)

    output$lod <- output$result.full[1,4]

    class(output) <- "fitqtl"
    attr(output, "method") <- method
    attr(output, "model") <- model
    attr(output, "formula") <- deparseQTLformula(formula)
    attr(output, "type") <- qtl$type
    attr(output, "nind") <- length(pheno)

    attr(output, "matrix.rank") <- matrix.rank
    attr(output, "matrix.ncol") <- matrix.ncol

    if(!is.null(residuals))
        attr(output, "residuals") <- residuals

    output
}


######################################################################
# checkformula
#
# check that input formula satisfies our imposed hiearchy: that
# main effects for terms in any interactions are also included
######################################################################
checkformula <-
    function(formula, qtl.name, covar.name)
{
    factors <- attr(terms(formula), "factors")
    altform <- deparseQTLformula(formula)

    if(sum(factors[1,])==0) factors <- factors[-1,,drop=FALSE]

    rn <- rownames(factors)
    # if mentions of "q1" or such, convert to "Q1" and such
    g <- grep("^[Qq][0-9]+$", rn)
    todrop <- NULL
    if(length(g) >= 1) {
        rownames(factors)[g] <- rn[g] <- toupper(rn[g])
        if(any(table(rn) > 1)) { # now there are some duplicates
            urn <- unique(rn)
            for(i in urn) {
                wh <- which(rn == i)
                if(length(wh) > 1) {
                    factors[wh[1],] <- apply(factors[wh,], 2, sum)
                    todrop <- c(todrop, wh[-1])
                    rownames(factors)[wh[1]] <- rn[wh[1]]
                }
            }
        }
    }
    if(length(todrop) > 0) factors <- factors[-todrop,]
    rn <- rownames(factors)

    if(!missing(qtl.name) || !missing(covar.name)) {
        m <- match(rn, c(qtl.name, covar.name))
        if(any(is.na(m)))
            warning("Terms ", paste(rn[is.na(m)], collapse=" "), " not understood.")
    }

    # paste rows together
    zo <- factors
    zo[zo>1] <- 1
    pzo <- apply(zo, 2, paste, collapse="")
    nt <- apply(zo, 2, sum)

    # form binary representations
    maxnt <- max(nt)
    v <- vector("list", maxnt)
    for(i in 2:maxnt) {
        v[[i]] <- mybinaryrep(i)
        v[[i]] <- v[[i]][,-ncol(v[[i]])+c(0,1)]
    }

    # for each higher-order column, form all lower-order terms
    for(i in which(nt > 1)) {
        cur <- zo[,i]
        wh <- which(cur==1)
        z <- v[[nt[i]]]
        zz <- matrix(0, ncol=ncol(z), nrow=length(cur))
        for(j in seq(along=wh))
            zz[wh[j],] <- z[j,]
        pzo <- unique(c(pzo, apply(zz, 2, paste, collapse="")))
    }

    zo <- matrix(as.numeric(unlist(strsplit(pzo, ""))), ncol=length(pzo))
    nt <- apply(zo, 2, sum)
    zo <- zo[,order(nt, apply(1-zo, 2, paste, collapse="")), drop=FALSE]
    rownames(zo) <- rn

    # form column names
    theterms <- apply(zo, 2, function(a, b) paste(b[as.logical(a)], collapse=":"),
                      rownames(zo))

    as.formula(paste("y ~ ", paste(theterms, collapse=" + ")))
}

#####################################################################
#
# parseformula
#
# Function to be called by fitqtl. It's used to
# parse the input formula
#
# This is the internal function and not supposed to be used by user
#
#####################################################################

parseformula <- function(formula, qtl.dimname, covar.dimname)
{
    # The terms for input formula
    f.formula <- terms(formula)
    order.term <- attr(f.formula, "order") # get the order of the terms
    idx.term <- which(order.term==1) # get the first order terms
    label.term <- attr(f.formula, "term.labels")[idx.term]
    formula.mtx <- attr(f.formula, "factors") # formula matrix

    idx.qtl <- NULL
    idx.covar <- NULL

    # loop thru all terms and find out how many QTLs and covariates
    # are there in the formula. Construct idx.qtl and idx.covar at the same time
    termisqtl <- rep(0, length(idx.term))
    for (i in 1:length(idx.term)) {
        # find out if there term is a QTL or a covariate
        # ignore the case for QTLs, e.g., Q1 is equivalent to q1
        idx.tmp <- grep(paste(label.term[i],"$", sep=""),
                        qtl.dimname, ignore.case=TRUE)
        if( length(idx.tmp) ) { # it's a QTL
            idx.qtl <- c(idx.qtl, idx.tmp)
            termisqtl[i] <- 1
        }
        else if(label.term[i] %in% covar.dimname) # it's a covariate
            idx.covar <- c(idx.covar, which(label.term[i]==covar.dimname))
        else
            stop("Unrecognized term ", label.term[i], " in formula")
    }
    n.qtl <- length(idx.qtl) # number of QTLs in formula
    n.covar <- length(idx.covar) # number of covariates in formula
    # now idx.qtl and idx.covar are the indices for genotype
    # and covariate matrices according to input formula

    # loop thru all terms again and reorganize formula.mtx
    formula.idx <- NULL
    ii <- 1
    jj <- 1
    for (i in 1:length(idx.term)) {
        #    if(label.term[i] %in% qtl.dimname) {  # it's a QTL
        if(termisqtl[i]) {
            formula.idx <- c(formula.idx, ii)
            ii <- ii+1
        }
        else { # it's a covariate
            formula.idx <- c(formula.idx, jj+n.qtl)
            jj <- jj+1
        }
    }

    # reorganize formula.mtx according to formula.idx
    # remove the first row (for y)
    formula.mtx <- formula.mtx[2:nrow(formula.mtx),]
    # rearrange the rows according to formula.idx if there's more than one row
    if(length(formula.idx) > 1)
        formula.mtx <- formula.mtx[order(formula.idx),]
    # take out only part of the matrix for interactions and pass to C function
    # all the input QTLs and covariates for C function will be additive
    n.int <- length(order.term) - length(idx.term) # number of interactions
    if(n.int != 0)
        formula.intmtx <- formula.mtx[,(length(idx.term)+1):length(order.term)]
    else # no interaction terms
        formula.intmtx <- NULL

    # return object
    result <- NULL
    result$idx.qtl <- idx.qtl
    result$n.qtl <- n.qtl
    result$idx.covar <- idx.covar
    result$n.covar <- n.covar
    result$formula.intmtx <- formula.intmtx
    result$n.int <- n.int

    result

}


#####################################################################
#
# summary.fitqtl
#
#####################################################################
summary.fitqtl <-
    function(object, pvalues=TRUE, simple=FALSE, ...)
{
    if(!any(class(object) == "fitqtl"))
        stop("Input should have class \"fitqtl\".")

    # this is just an interface.
    if("ests" %in% names(object)) {
        ests <- object$ests$ests
        se <- sqrt(diag(object$ests$covar))
        object$ests <- cbind(est=ests, SE=se, t=ests/se)
    }
    class(object) <- "summary.fitqtl"
    if(simple) pvalues <- FALSE
    attr(object, "pvalues") <- pvalues
    attr(object, "simple") <- simple
    object
}


#####################################################################
#
# print.summary.fitqtl
#
#####################################################################
print.summary.fitqtl <- function(x, ...)
{
    cat("\n")
    cat("\t\tfitqtl summary\n\n")
    meth <- attr(x, "method")
    mod <- attr(x, "model")
    if(is.null(mod)) mod <- "normal"
    if(meth=="imp") meth <- "multiple imputation"
    else if(meth=="hk") meth <- "Haley-Knott regression"
    cat("Method:", meth, "\n")
    cat("Model: ", mod, "phenotype\n")
    cat("Number of observations :", attr(x, "nind"), "\n\n")

    # print ANOVA table for full model
    cat("Full model result\n")
    cat("----------------------------------  \n")
    cat("Model formula:")
    w <- options("width")[[1]]
    printQTLformulanicely(attr(x, "formula"), "                   ", w+5, w)
    cat("\n")

    pval <- attr(x, "pvalues")
    simple <- attr(x, "simple")
    if(!is.null(pval) && !pval)
        x$result.full <- x$result.full[,-ncol(x$result.full)+(0:1)]
    if(mod=="binary" || (!is.null(simple) && simple))
        x$result.full <- x$result.full[1,-c(2:3,7),drop=FALSE]

    print(x$result.full, quote=FALSE, na.print="")
    cat("\n")

    # print ANOVA table for dropping one at a time analysis (if any)
    if("result.drop" %in% names(x)) {
        cat("\n")
        cat("Drop one QTL at a time ANOVA table: \n")
        cat("----------------------------------  \n")
        # use printCoefmat instead of print.data.frame
        # make sure the last column is P value
        if(!is.null(pval) && !pval)
            x$result.drop <- x$result.drop[,-ncol(x$result.drop)+(0:1)]
        if(mod=="binary" || (!is.null(simple) && simple))
            x$result.drop <- x$result.drop[,-c(2,5,7),drop=FALSE]

        printCoefmat(x$result.drop, digits=4, cs.ind=1, P.values=TRUE, has.Pvalue=TRUE)
        cat("\n")
    }

    if("ests" %in% names(x)) {
        cat("\n")
        cat("Estimated effects:\n")
        cat("-----------------\n")
        printCoefmat(x$ests,digits=4)
        cat("\n")
    }

}

######################################################################
# binary repreentation of the numbers 1...2^n;
# used in checkformula
######################################################################
mybinaryrep <-
    function(n)
{
    lx <- 2^n
    x <- 1:lx
    ans <- 0:(n-1)
    x <- matrix(rep(x,rep(n, lx)), ncol=lx)
    (x %/% 2^ans) %% 2
}

######################################################################
# deparseQTLformula: turn QTL formula into a string
######################################################################
deparseQTLformula <-
    function(formula, reorderterms=FALSE)
{
    if(is.null(formula)) return(NULL)
    if(reorderterms) {
        if(is.character(formula)) formula <- as.formula(formula)
        factors <- colnames(attr(terms(formula), "factors"))
        wh <- grep("^[Qq][0-9]+$", factors)
        if(length(wh)>0)
            factors[wh] <- paste("Q", sort(as.numeric(substr(factors[wh], 2, nchar(factors[wh])))), sep="")
        wh <- grep(":", factors)
        if(length(wh)>0) {
            temp <- strsplit(factors[wh], ":")
            temp <- sapply(temp, function(a) {
                wh <- grep("^[Qq][0-9]+$", a)
                if(any(wh)) a[wh] <- paste("Q", sort(as.numeric(substr(a[wh], 2, nchar(a[wh])))), sep="")
                paste(a[order(is.na(match(seq(along=a),wh)))], collapse=":")
            })
            factors[wh] <- temp
        }
        return(paste("y ~ ", paste(factors, collapse=" + "), sep=""))
    }


    if(is.character(formula)) return(formula)
    paste(as.character(formula)[c(2,1,3)], collapse=" ")
}

printQTLformulanicely <-
    function(formula, header, width, width2, sep=" ")
{
    if(!is.character(formula)) formula <- deparseQTLformula(formula)
    thetext <- unlist(strsplit(formula, " "))

    if(missing(width2)) width2 <- width
    nleft <- width - nchar(header)
    nsep <- nchar(sep)
    if(length(thetext) < 2) cat("", thetext, "\n", sep=sep)
    else {
        z <- paste("", thetext[1], sep=sep, collapse=sep)
        for(j in 2:length(thetext)) {
            if(nchar(z) + nsep + nchar(thetext[j]) > nleft) {
                cat(z, "\n")
                nleft <- width2
                z <- paste(header, thetext[j], sep=sep)
            }
            else {
                z <- paste(z, thetext[j], sep=sep)
            }
        }
        cat(z, "\n")
    }
}

# end of fitqtl.R
