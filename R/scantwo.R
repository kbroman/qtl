######################################################################
#
# scantwo.R
#
# copyright (c) 2001-2015, Karl W Broman and Hao Wu
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
# Hao Wu (The Jackson Lab) wrote the initial code for the imputation
# method.
#
# Part of the R/qtl package
# Contains: scantwo, scantwo.perm, scantwo.perm.engine
#
######################################################################

######################################################################
#
# scantwo: Do 2-dimensional genome scan with a two-QTL model,
#          calculating joint LOD scores and LOD scores testing
#          epistasis.
#
######################################################################

scantwo <-
    function(cross, chr, pheno.col=1,
             model=c("normal","binary"),
             method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"),
             incl.markers=FALSE, clean.output=FALSE, clean.nmar=1,
             clean.distance=0,
             maxit=4000, tol=1e-4, verbose=TRUE, n.perm,
             perm.Xsp=FALSE, perm.strata=NULL, assumeCondIndep=FALSE,
             batchsize=250, n.cluster=1)
{
    if(batchsize < 1) stop("batchsize must be >= 1.")
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    method <- match.arg(method)
    model <- match.arg(model)
    use <- match.arg(use)

    # pull out chromosomes to be scanned
    if(missing(chr)) chr1 <- chr2 <- chr <- names(cross$geno)
    else {
        thechr <- names(cross$geno)
        if(is.list(chr)) {
            # special case: do just specific pairs (each of chr1 vs each of chr2, except when chr2 < chr1)
            chr1 <- matchchr(chr[[1]], thechr)
            chr2 <- matchchr(chr[[2]], thechr)
        }
        else
            chr1 <- chr2 <- matchchr(chr, thechr)
    }

    cross <- subset(cross, unique(c(chr1, chr2)))
    thechr <- names(cross$geno)
    nchr1 <- match(chr1, thechr)
    nchr2 <- match(chr2, thechr)
    if(!any(sapply(nchr1, function(a,b) any(a <= b), nchr2)))
        stop("Need some of first chr to be <= some of second chr")

    if(missing(n.perm)) n.perm <- 0

    if((method=="hk" || method=="em") && !assumeCondIndep) { # if reduce2grid was used, for assumeCondIndep
        # if reduced2grid, force assumeCondIndep=TRUE
        reduced2grid <- attr(cross$geno[[1]]$prob, "reduced2grid")
        if(!is.null(reduced2grid) && reduced2grid) {
            assumeCondIndep <- TRUE
            warning("Using assumeCondIndep=TRUE, since probabilities reduced to grid")
        }
    }

    # in RIL, treat X chromomse like an autosome
    chrtype <- sapply(cross$geno, class)
    if(any(chrtype=="X") && (class(cross)[1] == "risib" || class(cross)[1] == "riself"))
        for(i in which(chrtype=="X")) class(cross$geno[[i]]) <- "A"

    # X-chr-specific perms (actually A:A, A:X, and X:X specific)
    if(perm.Xsp && n.perm > 0) {
        if(all(chrtype=="X") || all(chrtype=="A")) break # no need for x-chr-specific perms
        if(length(chr1) != length(chr2) || any(chr1 != chr2))
            stop("With perm.Xsp=TRUE, chr can't be a list")

        chr.names <- chrnames(cross)
        chrL <- sapply(cross$geno, function(a) diff(range(a$map)))
        AL <- sum(chrL[chrtype=="A"])
        XL <- sum(chrL[chrtype=="X"])
        AAL <- AL*AL/2
        XXL <- XL*XL/2
        AXL <- AL*XL
        n.permAA <- n.perm
        n.permXX <- ceiling(n.perm * AAL/XXL)
        n.permAX <- ceiling(n.perm * AAL/AXL)

        # names of autosomes and X chr
        Achr <- chr.names[chrtype=="A"]
        Xchr <- chr.names[chrtype=="X"]

        # X chr covariates
        Xnull <- scanoneXnull(class(cross)[1], getsex(cross), attributes(cross))
        Xcovar <- Xnull$sexpgmcovar

        if(verbose) message("Running ", n.permAA, " A:A permutations")
        AAresult <- scantwo(cross, chr=Achr, pheno.col=pheno.col,
                            model=model, method=method,
                            addcovar=cbind(addcovar, Xcovar), intcovar=intcovar,
                            weights=weights, use=use,
                            incl.markers=incl.markers, clean.output=clean.output,
                            clean.nmar=clean.nmar, clean.distance=clean.distance,
                            maxit=maxit, tol=tol, verbose=verbose, n.perm=n.permAA,
                            perm.Xsp=FALSE, perm.strata=perm.strata,
                            assumeCondIndep=assumeCondIndep,
                            batchsize=batchsize, n.cluster=n.cluster)

        if(verbose) message("Running ", n.permXX, " X:X permutations")
        XXresult <- scantwo(cross, chr=Xchr, pheno.col=pheno.col,
                            model=model, method=method,
                            addcovar=addcovar, intcovar=intcovar,
                            weights=weights, use=use,
                            incl.markers=incl.markers, clean.output=clean.output,
                            clean.nmar=clean.nmar, clean.distance=clean.distance,
                            maxit=maxit, tol=tol, verbose=verbose, n.perm=n.permXX,
                            perm.Xsp=FALSE, perm.strata=perm.strata,
                            assumeCondIndep=assumeCondIndep,
                            batchsize=batchsize, n.cluster=n.cluster)

        if(verbose) message("Running ", n.permAX, " A:X permutations")
        AXresult <- scantwo(cross, chr=list(Achr, Xchr), pheno.col=pheno.col,
                            model=model, method=method,
                            addcovar=addcovar, intcovar=intcovar,
                            weights=weights, use=use,
                            incl.markers=incl.markers, clean.output=clean.output,
                            clean.nmar=clean.nmar, clean.distance=clean.distance,
                            maxit=maxit, tol=tol, verbose=verbose, n.perm=n.permAX,
                            perm.Xsp=FALSE, perm.strata=perm.strata,
                            assumeCondIndep=assumeCondIndep,
                            batchsize=batchsize, n.cluster=n.cluster)

        result <- list(AA=AAresult, AX=AXresult, XX=XXresult)
        attr(result, "L") <- c(A=AL, X=XL)
        attr(result, "LL") <- c(AA=AAL, AX=AXL, XX=XXL)
        names(chrtype) <- chr.names
        attr(result, "chrtype") <- chrtype
        class(result) <- c("scantwoperm", "list")
        return(result)
    }

    if(!missing(n.perm) && n.perm > 0 && n.cluster > 1) {
        cat(" -Running permutations via a cluster of", n.cluster, "nodes.\n")
        updateParallelRNG(n.cluster)
        n.perm <- ceiling(n.perm/n.cluster)
        scantwoPermInParallel <- function(n.perm, cross, chr, pheno.col, model, method, addcovar, intcovar, weights,
                                          incl.markers, clean.output, clean.nmar, clean.distance, maxit, tol,
                                          perm.strata, assumeCondIndep, batchsize)
            scantwo(cross=cross, chr=chr, pheno.col=pheno.col, model=model, method=method, addcovar=addcovar,
                    intcovar=intcovar, weights=weights, incl.markers=incl.markers, clean.output=clean.output,
                    clean.nmar=clean.nmar, clean.distance=clean.distance, maxit=maxit, tol=tol, perm.strata=perm.strata,
                    assumeCondIndep=assumeCondIndep, batchsize=batchsize, n.cluster=0, verbose=FALSE, n.perm=n.perm)

        if(Sys.info()[1] == "Windows") { # Windows doesn't support mclapply, but it's faster if available
            cl <- makeCluster(n.cluster)
            on.exit(stopCluster(cl))
            operm <- clusterApply(cl, rep(n.perm, n.cluster), scantwoPermInParallel, cross=cross, chr=chr, pheno.col=pheno.col,
                                  model=model, method=method, addcovar=addcovar, intcovar=intcovar,
                                  weights=weights, incl.markers=incl.markers, clean.output=clean.output,
                                  clean.nmar=clean.nmar, clean.distance=clean.distance,
                                  maxit=maxit, tol=tol, perm.strata=perm.strata,
                                  assumeCondIndep=assumeCondIndep, batchsize=batchsize)
        }
        else {
            operm <- mclapply(rep(n.perm, n.cluster), scantwoPermInParallel, cross=cross, chr=chr, pheno.col=pheno.col,
                              model=model, method=method, addcovar=addcovar, intcovar=intcovar,
                              weights=weights, incl.markers=incl.markers, clean.output=clean.output,
                              clean.nmar=clean.nmar, clean.distance=clean.distance,
                              maxit=maxit, tol=tol, perm.strata=perm.strata,
                              assumeCondIndep=assumeCondIndep, batchsize=batchsize, mc.cores=n.cluster)
        }
        for(j in 2:length(operm))
            operm[[1]] <- c(operm[[1]], operm[[j]])
        return(operm[[1]])
    }

    # check perm.strat
    if(!missing(perm.strata) && !is.null(perm.strata)) {
        if(length(perm.strata) != nind(cross))
            stop("perm.strata, if given, must have length = nind(cross) [", nind(cross), "]")
    }

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    origcross <- cross

    fullmap <- pull.map(cross)

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

    # if stepwidth="variable" or stepwidth=="max" when calling calc.genoprob or sim.geno,
    # we force incl.markers=TRUE; I assume it is the same for all chromosomes
    stepwidth.var <- FALSE
    if(method=="em" || method=="hk") {
        if("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
           attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
            stepwidth.var <- TRUE
            incl.markers <- TRUE
        }
    }
    else if(method=="imp") {
        if("stepwidth" %in% names(attributes(cross$geno[[1]]$draws)) &&
           attr(cross$geno[[1]]$draws, "stepwidth") != "fixed") {
            stepwidth.var <- TRUE
            incl.markers <- TRUE
        }
    }

    if(length(pheno.col)==1 && n.perm>=0) use <- "complete.obs"

    if(n.perm >= 0) {    # not in the midst of a permutation test
        # If use="all.obs", check whether there are individuals missing some
        # phenotypes but not others.  If not, act like "complete.obs".
        if(use=="all.obs" && length(pheno.col) > 1) {
            n.phe <- length(pheno.col)
            temp <- apply(cross$pheno[,pheno.col], 1, function(a) sum(is.na(a)))
            if(all(temp==0 | temp==n.phe)) use <- "complete.obs"
        }

        # If use="complete.obs", drop individuals with any missing phenotypes
        if(use=="complete.obs") {
            temp <- checkcovar(cross, pheno.col, addcovar, intcovar,
                               perm.strata, ind.noqtl=NULL, weights, TRUE)
            cross <- temp[[1]]
            pheno <- temp[[2]]
            addcovar <- temp[[3]]
            intcovar <- temp[[4]]
            n.addcovar <- temp[[5]]
            n.intcovar <- temp[[6]]
            perm.strata <- temp[[7]]
            weights <- temp[[9]]
        }
    }

    # use all observations; not in a permutation test; different phenotypes have different sets of missing values
    #   -> want to do in batches, but need to define batches by the pattern of missing data
    if(n.perm <= 0 && use=="all.obs" && length(pheno.col) > 1 && (method=="hk" || method=="imp")) {
        # drop individuals with missing covariates
        cross$pheno <- cbind(cross$pheno, rep(1, nind(cross)))
        temp <- checkcovar(cross, nphe(cross), addcovar, intcovar,
                           perm.strata, ind.noqtl=NULL, weights, TRUE)
        cross <- temp[[1]]
        pheno <- cross$pheno[,pheno.col, drop=FALSE]
        addcovar <- temp[[3]]
        intcovar <- temp[[4]]
        n.addcovar <- temp[[5]]
        n.intcovar <- temp[[6]]
        perm.strata <- temp[[7]]
        weights <- temp[[9]]

        # determine the batches (defined by the pattern of missing data)
        patterns <- apply(pheno, 2, function(a) paste(!is.na(a), collapse=":"))
        upat <- unique(patterns)
        m <- match(patterns, upat)
        batches <- vector("list", length(upat))
        upat <- lapply(strsplit(upat, ":"), function(a) as.logical(a))
        for(i in seq(along=batches)) batches[[i]] <- pheno.col[m==i]

        # run scanone for one batch at a time
        out <- NULL
        for(i in seq(along=batches)) {
            if(!is.null(addcovar)) {
                if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
                tempac <- addcovar[upat[[i]],,drop=FALSE]
            }
            else tempac <- addcovar
            if(!is.null(intcovar)) {
                if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)
                tempic <- intcovar[upat[[i]],,drop=FALSE]
            }
            else tempic <- intcovar

            temp <- scantwo(subset(cross, ind=upat[[i]]), chr=chr, pheno.col=batches[[i]], model=model,
                            method=method, addcovar=tempac, intcovar=tempic,
                            weights=weights[upat[[i]]], use=use, incl.markers=incl.markers, clean.output=clean.output,
                            clean.nmar=clean.nmar, clean.distance=clean.distance,
                            maxit=maxit, tol=tol, verbose=verbose, n.perm=n.perm,
                            perm.strata=perm.strata[upat[[i]]], assumeCondIndep=assumeCondIndep,
                            batchsize=batchsize, n.cluster=n.cluster)
            if(is.null(out)) out <- temp
            else out <- cbind(out, temp)
        }

        # reorder LOD score columns and make sure that the names are correct
        dimnames(out$lod) <- list(NULL, NULL, colnames(cross$pheno)[unlist(batches)])
        out$lod <- out$lod[,,colnames(cross$pheno)[pheno.col]]
        dimnames(out$lod)[[3]] <- colnames(cross$pheno)[pheno.col]
        attr(out,"phenotypes") <- colnames(cross$pheno)[pheno.col]

        return(out)
    }

    # multiple phenotype for methods other than imp and hk
    if(length(pheno.col)>1 && n.perm <= 0 && (model != "normal" ||
                              (method!="imp" && method != "hk"))) {
        n.phe <- length(pheno.col)
        if(verbose) cat(" -Phenotype 1\n")
        output <- scantwo(cross, pheno.col=pheno.col[1], model=model,
                          method=method, addcovar=addcovar,
                          intcovar=intcovar, weights=weights, use=use,
                          incl.markers=incl.markers,
                          clean.output=clean.output, clean.nmar=clean.nmar,
                          clean.distance=clean.distance,
                          maxit=maxit, tol=tol, verbose=verbose, n.perm=n.perm,
                          perm.strata=perm.strata,
                          assumeCondIndep=assumeCondIndep,
                          batchsize=batchsize, n.cluster=0)
        temp <- array(dim=c(nrow(output$lod), ncol(output$lod), n.phe))
        temp[,,1] <- output$lod
        output$lod <- temp
        for(i in 2:n.phe) {
            if(verbose) cat(" -Phenotype ", i, "\n")
            temp <- scantwo(cross, pheno.col=pheno.col[i], model=model, method=method,
                            addcovar=addcovar, intcovar=intcovar, weights=weights,
                            use=use, incl.markers=incl.markers,
                            clean.output=clean.output, clean.nmar=clean.nmar,
                            clean.distance=clean.distance,
                            maxit=maxit, tol=tol, verbose=verbose, n.perm=n.perm,
                            perm.strata=perm.strata,
                            assumeCondIndep=assumeCondIndep,
                            batchsize=batchsize, n.cluster=0)
            output$lod[,,i] <- temp$lod
            if(!is.null(output$scanoneX))
                output$scanoneX <- cbind(output$scanoneX, temp$scanoneX)
        }
        attr(output,"fullmap") <- fullmap
        attr(output,"phenotypes") <- colnames(cross$pheno)[pheno.col]
        names(output$map)[2] <- "pos"
        dimnames(output$lod) <- list(NULL, NULL, colnames(cross$pheno)[pheno.col])
        return(output)
    }

    # if n.perm specified, do a permutation test
    if(n.perm>0) {
        return(scantwo.perm(cross, pheno.col=pheno.col, model=model, method=method,
                            addcovar=addcovar, intcovar=intcovar, weights=weights,
                            use=use, incl.markers=incl.markers, clean.output=clean.output,
                            clean.nmar=clean.nmar,
                            clean.distance=clean.distance,
                            maxit=maxit,
                            tol=tol, verbose=verbose, n.perm=n.perm,
                            perm.strata=perm.strata, assumeCondIndep=assumeCondIndep,
                            batchsize=batchsize, chr=chr))
    }


    if(n.perm < 0) { # in the midst of permutations
        if(use=="all.obs") {
            temp <- checkcovar(cross, pheno.col, addcovar, intcovar,
                               perm.strata, ind.noqtl=NULL, weights, n.perm==-1)
            cross <- temp[[1]]
            pheno <- temp[[2]]
            addcovar <- temp[[3]]
            intcovar <- temp[[4]]
            n.addcovar <- temp[[5]]
            n.intcovar <- temp[[6]]
            perm.strata <- temp[[7]]
            weights <- temp[[9]]
        }
        else {
            pheno <- as.matrix(cross$pheno[,pheno.col])
            if(is.null(addcovar)) n.addcovar <- 0
            else n.addcovar <- ncol(addcovar)
            if(is.null(intcovar)) n.intcovar <- 0
            else n.intcovar <- ncol(intcovar)
        }
    }

    n.chr <- nchr(cross)
    n.ind <- nind(cross)
    n.phe <- length(pheno.col)
    type <- class(cross)[1]
    chrtype <- sapply(cross$geno,class)
    is.bcs <- FALSE
    if(type == "bcsft") {
        cross.scheme <- attr(cross, "scheme")
        is.bcs <- (cross.scheme[2] == 0)
    }

    if(any(chrtype=="X")) {
        sexpgm <- getsex(cross)

        addcovarX <- revisecovar(sexpgm,addcovar)
        if(!is.null(addcovar) && (nd <- attr(addcovarX, "n.dropped")) > 0 && n.perm > -2)
            warning("Dropped ", nd, " additive covariates on X chromosome.")
        if(length(addcovarX)==0) {
            n.acX <- 0
            addcovarX <- NULL
        }
        else n.acX <- ncol(addcovarX)

        intcovarX <- revisecovar(sexpgm,intcovar)
        if(!is.null(intcovar) && (nd <- attr(intcovarX, "n.dropped")) > 0 && n.perm > -2)
            warning("Dropped ", nd, " interactive covariates on X chromosome.")
        if(length(intcovarX)==0) {
            n.icX <- 0
            intcovarX <- NULL
        }
        else n.icX <- ncol(intcovarX)
    }

    if(model=="binary") {
        if(method != "em" && method != "hk") {
            method <- "em"
            if(n.perm > -2) warning("Only EM algorithm and Haley-Knott regression coded for binary traits; using EM")
        }
        if(!is.null(weights)) {
            weights <- NULL
            if(n.perm > -2) warning("weights ignored for binary traits.")
        }

        u <- unique(pheno)
        if(any(u!=0 & u!=1))
            stop("Phenotypes must be either 0 or 1.")
    }

    if(n.perm == 0) { # not in the midst of permutations
        if(method=="mr-argmax")
            cross <- fill.geno(cross,method="argmax")
        if(method=="mr-imp")
            cross <- fill.geno(cross,method="imp")
    }

    # weights of individuals
    if(model == "normal") {
        if(is.null(weights))
            weights <- rep(1, nind(cross))
        if(length(weights) != nind(cross))
            stop("weights should either be NULL or a vector of length n.ind")
        if(any(weights <= 0))
            stop("weights should be entirely positive")
        weights <- sqrt(weights)
    }

    if(verbose) cat(" --Running scanone\n")
    temp <- scanone(cross, pheno.col=pheno.col, model=model, method=method,
                    addcovar=addcovar, intcovar=intcovar, weights=weights,
                    use=use, maxit=maxit, tol=tol, verbose=FALSE)
    out.scanone <- temp[,-(1:2),drop=FALSE]
    if(verbose) cat(" --Running scantwo\n")

    if(method=="mr" || method=="mr-imp" || method=="mr-argmax") { # marker regression
        # number of genotypes on each chromosome,
        #     combine the genetic maps for all chromosomes
        map <- unlist(pull.map(cross))
        names(map) <- unlist(lapply(pull.map(cross),names))
        n.pos <- nmar(cross)
        gmap <- data.frame(chr=rep(names(cross$geno),n.pos),
                           pos=map,
                           eq.spacing=rep(1,sum(n.pos)),
                           xchr=rep(sapply(cross$geno,class)=="X",nmar(cross)), stringsAsFactors=TRUE)

        # number of possible genotypes for each chromosome
        n.gen <- 1:n.chr
        for(i in 1:n.chr) {
            gen.names <- getgenonames(type, chrtype[i], "full", sexpgm, attributes(cross))
            n.gen[i] <- length(gen.names)
        }
    } # end of if(method=="mr")

    else { # all methods except "mr"
        # check for genotype probabilities or simulated genotypes
        steps <- rep(0,n.chr) # step length on each chromosome
        if(method=="imp") {
            for(i in 1:n.chr) {
                if(!("draws" %in% names(cross$geno[[i]]))) {
                    # need to run sim.geno
                    if(n.perm > -2) warning("First running sim.geno.")
                    cross <- sim.geno(cross)
                }
                steps[i] <- attr(cross$geno[[i]]$draws,"step")
            }

            # make sure all chromosomes have the same number of imputations
            n.draws <- sapply(cross$geno, function(a) dim(a$draws)[3])
            if(length(unique(n.draws)) > 1) {
                if(n.perm > -2) warning("Re-running sim.geno to have a fixed number of imputations.")
                cross <- sim.geno(cross, n.draws=max(n.draws),
                                  step=attr(cross$geno[[1]]$draws,"step"),
                                  off.end=attr(cross$geno[[1]]$draws,"off.end"))
            }
            n.draws <- max(n.draws)
        }
        else { # H-K or EM
            for(i in 1:n.chr) {
                if(!("prob" %in% names(cross$geno[[i]]))) {
                    # need to run calc.genoprob
                    if(n.perm > -2) warning("First running calc.genoprob.")
                    cross <- calc.genoprob(cross)
                }
                steps[i] <- attr(cross$geno[[i]]$prob,"step")
            }
        }

        # number of genotypes on each chromosome,
        #     construct the genetic map for all chromosomes
        #     and possibly drop marker positions
        gmap <- NULL
        n.pos <- n.gen <- rep(0,n.chr)
        keep.pos <- vector("list",n.chr)
        some.dropped <- rep(FALSE,n.chr)

        for(i in 1:n.chr) {
            gen.names <- getgenonames(type, chrtype[i], "full", sexpgm, attributes(cross))
            n.gen[i] <- length(gen.names)

            # construct the genetic map for this chromesome
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

            if(is.matrix(map)) map <- map[1,] # in case of sex-specific map

            w <- names(map)
            o <- grep("^loc-*[0-9]+",w)

            if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
                w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
            map <- cbind(chr=rep(names(cross$geno)[i],length(map)),
                         pos=as.data.frame(as.numeric(map), stringsAsFactors=TRUE) )
            rownames(map) <- w

            # equally spaced positions
            if(steps[i]==0 || stepwidth.var)  # just use markers
                eq.sp.pos <- rep(1,nrow(map))
            else {
                eq.sp.pos <- seq(min(map[,2]),max(map[,2]),by=steps[i])
                wh.eq.sp <- match(eq.sp.pos,map[,2])
                if(any(is.na(wh.eq.sp))) { # this shouldn't happen
                    if(n.perm > -2) warning("Possible error in determining the equally spaced positions.")
                    wh.eq.sp <- wh.eq.sp[!is.na(wh.eq.sp)]
                }
                eq.sp.pos <- rep(0,nrow(map))
                eq.sp.pos[wh.eq.sp] <- 1
            }
            if(!incl.markers && any(eq.sp.pos==0)) {
                keep.pos[[i]] <- (seq(along=eq.sp.pos))[eq.sp.pos==1]
                map <- map[eq.sp.pos==1,]
                eq.sp.pos <- eq.sp.pos[eq.sp.pos==1]
                some.dropped[i] <- TRUE # indicates some positions were dropped
            }
            else keep.pos[[i]] <- seq(along=eq.sp.pos)
            gmap <- rbind(gmap, cbind(map,eq.spacing=eq.sp.pos,
                                      xchr=(class(cross$geno[[i]])=="X")))
            n.pos[i] <- length(keep.pos[[i]])

            # Revise X chromosome genotype probabilities or imputations
            if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft"))) {
                if(method=="imp")
                    cross$geno[[i]]$draws <-
                        reviseXdata(type, "full", sexpgm, draws=cross$geno[[i]]$draws,
                                    cross.attr=attributes(cross))
                else if(method=="hk" || method=="em") {
                    oldXchr <- subset(cross, chr=thechr[i])
                    cross$geno[[i]]$prob <-
                        reviseXdata(type, "full", sexpgm, prob=cross$geno[[i]]$prob,
                                    cross.attr=attributes(cross))
                }
                else
                    cross$geno[[i]]$data <-
                        reviseXdata(type, "full", sexpgm, geno=cross$geno[[i]]$data,
                                    cross.attr=attributes(cross))
            }

        } # end loop over chromosomes
    } # end of if/else for method="mr" vs other

    # columns in result matrix for each chromosome
    wh.col <- vector("list",n.chr)
    first.pos <- cumsum(c(1,n.pos))
    for(i in 1:n.chr)
        wh.col[[i]] <- seq(first.pos[i],by=1,length=n.pos[i])

    # initialize the results matrix
    if(n.phe > 1)
        results <- array(NA,dim=c(sum(n.pos),sum(n.pos), n.phe))
    else
        results <- matrix(NA,sum(n.pos),sum(n.pos))

    # do the 2-dimensional genome scan
    do.nllik0 <- TRUE
    for(i in nchr1) { # loop over the 1st chromosome
        for(j in nchr2) { # loop over the 2nd chromosome
            if(j < i) next

            if(chrtype[i]=="X" || chrtype[j]=="X") {
                ac <- addcovarX
                n.ac <- n.acX
                ic <- intcovarX
                n.ic <- n.icX
            }
            else {
                ac <- addcovar
                n.ac <- n.addcovar
                ic <- intcovar
                n.ic <- n.intcovar
            }
            if(i==j && chrtype[i]=="X") {
                col2drop <- dropXcol(type, sexpgm, attributes(cross))
                n.col2drop <- sum(col2drop)
                n.col2drop.addmodel <- sum(col2drop[1:(2*n.gen[i]-1)])
            }
            else {
                col2drop <- rep(0,n.gen[i]*n.gen[j])
                n.col2drop <- 0
            }

            # print the current working pair
            if(verbose) cat(" (", names(cross$geno)[i], ",",
                            names(cross$geno)[j],")\n",sep="")

            if(method=="imp") {
                if(n.phe > batchsize) {
                    firstcol <- 1
                    z <- NULL
                    while(firstcol <= n.phe) {
                        thiscol <- firstcol + 0:(batchsize-1)
                        thiscol <- thiscol[thiscol <= n.phe]
                        thisz <- .C("R_scantwo_imp",
                                    as.integer(n.ind),
                                    as.integer(i==j),
                                    as.integer(n.pos[i]),
                                    as.integer(n.pos[j]),
                                    as.integer(n.gen[i]),
                                    as.integer(n.gen[j]),
                                    as.integer(n.draws),
                                    as.integer(cross$geno[[i]]$draws[,keep.pos[[i]],]),
                                    as.integer(cross$geno[[j]]$draws[,keep.pos[[j]],]),
                                    as.double(ac),
                                    as.integer(n.ac),
                                    as.double(ic),
                                    as.integer(n.ic),
                                    as.double(pheno[,thiscol]),
                                    as.integer(length(thiscol)),
                                    as.double(weights),
                                    result=as.double(rep(0,2*n.pos[i]*n.pos[j]*length(thiscol))),
                                    as.integer(n.col2drop),
                                    as.integer(col2drop),
                                    PACKAGE="qtl")
                        firstcol <- firstcol + batchsize
                        if(is.null(z))
                            z <- array(NA, dim=c(n.pos[i], n.pos[j], 2*n.phe))
                        thisz$result <- array(thisz$result, dim=c(n.pos[i], n.pos[j], 2*length(thiscol)))
                        z[,,thiscol] <- thisz$result[,,1:length(thiscol)]
                        z[,,n.phe+thiscol] <- thisz$result[,,length(thiscol)+1:length(thiscol)]
                    }
                }
                else {
                    z <- .C("R_scantwo_imp",
                            as.integer(n.ind),
                            as.integer(i==j),
                            as.integer(n.pos[i]),
                            as.integer(n.pos[j]),
                            as.integer(n.gen[i]),
                            as.integer(n.gen[j]),
                            as.integer(n.draws),
                            as.integer(cross$geno[[i]]$draws[,keep.pos[[i]],]),
                            as.integer(cross$geno[[j]]$draws[,keep.pos[[j]],]),
                            as.double(ac),
                            as.integer(n.ac),
                            as.double(ic),
                            as.integer(n.ic),
                            as.double(pheno),
                            as.integer(n.phe),
                            as.double(weights),
                            result=as.double(rep(0,2*n.pos[i]*n.pos[j]*n.phe)),
                            as.integer(n.col2drop),
                            as.integer(col2drop),
                            PACKAGE="qtl")
                    z <- array(z$result,dim=c(n.pos[i], n.pos[j], 2*n.phe)) # rearrange the result
                }

                # do this just once: do null model and get neg log10 likelihood
                if(do.nllik0) {
                    do.nllik0 <- FALSE
                    if(n.ac > 0)
                        resid0 <- lm(pheno ~ ac, weights=weights^2)$resid
                    else
                        resid0 <- lm(pheno ~ 1, weights=weights^2)$resid
                    sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
                    nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
                }

                # update the final result matrix
                if(i == j) { # on same chromosome
                    if(n.phe > 1)
                        results[wh.col[[i]],wh.col[[j]],] <- z[,,1:n.phe]
                    else
                        results[wh.col[[i]],wh.col[[j]]] <- z[,,1]

                }
                else { # on different chromosomes
                    if(n.phe > 1) {
                        # full lod
                        results[wh.col[[i]],wh.col[[j]],] <- z[,,1:n.phe]
                        # epistasis lod - need to reshape the matrix
                        results[wh.col[[j]],wh.col[[i]],] <- array(z[,,1:n.phe+n.phe],
                                                                   c(n.pos[j],n.pos[i], n.phe))
                    }
                    else { # only one phenotype, result is a matrix
                        # full lod
                        results[wh.col[[i]],wh.col[[j]]] <- z[,,1]
                        # epistasis - need to reshape the matrix
                        results[wh.col[[j]],wh.col[[i]]] <- matrix(z[,,2],n.pos[j],n.pos[i])
                    }
                }
            }
            else if(model=="normal" && (method=="hk" || method=="em")) {
                if(do.nllik0) { # first time! do null model and get neg log10 likelihood
                    do.nllik0 <- FALSE
                    if(n.ac > 0)
                        resid0 <- lm(pheno ~ ac, weights=weights^2)$resid
                    else
                        resid0 <- lm(pheno ~ 1, weights=weights^2)$resid
                    if(method=="hk") {
                        if(n.phe == 1)
                            nllik0 <- (n.ind/2)*log10(sum((resid0*weights)^2))
                        else # multiple phenotypes
                            nllik0 <- apply(resid0, 2, function(x)
                                            (n.ind/2)*log10(sum((x*weights)^2)))
                    }
                    else {
                        sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
                        nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
                    }
                }

                if(i==j) { # same chromosome


                    if(verbose>1) cat("  --Calculating joint probs.\n")

                    if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft"))) {
                        # calculate joint genotype probabilities for all pairs of positions
                        stp <- attr(oldXchr$geno[[1]]$prob, "step")
                        oe <- attr(oldXchr$geno[[1]]$prob, "off.end")
                        err <- attr(oldXchr$geno[[1]]$prob, "error.prob")
                        mf <- attr(oldXchr$geno[[1]]$prob, "map.function")

                        if("stepwidth" %in% names(attributes(oldXchr$geno[[1]]$prob)))
                            stpw <- attr(oldXchr$geno[[1]]$prob, "stepwidth")
                        else stpw <- "fixed"
                        if("map" %in% names(attributes(oldXchr$geno[[1]]$prob)))
                            tmap <- attr(oldXchr$geno[[1]]$prob,"map")
                        else
                            tmap <- create.map(oldXchr$geno[[1]]$map, stp, oe, stpw)

                        temp <- calc.pairprob(oldXchr,stp,oe,err,mf,tmap,
                                              assumeCondIndep=assumeCondIndep)
                    }
                    else {
                        # calculate joint genotype probabilities for all pairs of positions
                        stp <- attr(cross$geno[[i]]$prob, "step")
                        oe <- attr(cross$geno[[i]]$prob, "off.end")
                        err <- attr(cross$geno[[i]]$prob, "error.prob")
                        mf <- attr(cross$geno[[i]]$prob, "map.function")

                        if("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
                            stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
                        else stpw <- "fixed"
                        if("map" %in% names(attributes(cross$geno[[i]]$prob)))
                            tmap <- attr(cross$geno[[i]]$prob,"map")
                        else
                            tmap <- create.map(cross$geno[[i]]$map, stp, oe, stpw)

                        temp <- calc.pairprob(subset(cross,chr=thechr[i]),stp,oe,err,mf,tmap,
                                              assumeCondIndep=assumeCondIndep)
                    }

                    # pull out positions from genotype probs
                    if(some.dropped[i]) {
                        # figure out pos'ns corresponding to columns of temp
                        nc <- ncol(cross$geno[[i]]$prob)
                        ind <- matrix(rep(1:nc,nc),ncol=nc)
                        w <- lower.tri(ind)
                        ind <- cbind(first=t(ind)[w],second=ind[w])

                        # which part to keep
                        keep <- apply(ind,1,function(a,b) all(a %in% b),
                                      keep.pos[[i]])
                        temp <- temp[,keep,,]
                    }

                    # revise pair probilities for X chromosome
                    if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft")))
                        temp <- reviseXdata(type, "full", sexpgm, pairprob=temp,
                                            cross.attr=attributes(cross))

                    if(verbose>1) cat("  --Done.\n")

                    if(method=="hk") {
                        if(n.phe > batchsize) {
                            firstcol <- 1
                            z <- NULL
                            while(firstcol <= n.phe) {
                                thiscol <- firstcol + 0:(batchsize-1)
                                thiscol <- thiscol[thiscol <= n.phe]

                                thisz <- .C("R_scantwo_1chr_hk",
                                            as.integer(n.ind),
                                            as.integer(n.pos[i]),
                                            as.integer(n.gen[i]),
                                            as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                            as.double(temp),
                                            as.double(ac),
                                            as.integer(n.ac),
                                            as.double(ic),
                                            as.integer(n.ic),
                                            as.double(pheno[,thiscol]),
                                            as.integer(length(thiscol)),
                                            as.double(weights),
                                            result=as.double(rep(0,n.pos[i]^2*length(thiscol))),
                                            as.integer(n.col2drop),
                                            as.integer(col2drop),
                                            PACKAGE="qtl")

                                firstcol <- firstcol + batchsize
                                if(is.null(z)) {
                                    z <- thisz
                                    z$result <- array(NA, dim=c(n.pos[i], n.pos[i], n.phe))
                                }
                                z$result[,,thiscol] <- array(thisz$result, dim=c(n.pos[i], n.pos[i], length(thiscol)))
                            }
                        }
                        else {
                            z <- .C("R_scantwo_1chr_hk",
                                    as.integer(n.ind),
                                    as.integer(n.pos[i]),
                                    as.integer(n.gen[i]),
                                    as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                    as.double(temp),
                                    as.double(ac),
                                    as.integer(n.ac),
                                    as.double(ic),
                                    as.integer(n.ic),
                                    as.double(pheno),
                                    as.integer(n.phe),
                                    as.double(weights),
                                    result=as.double(rep(0,n.pos[i]^2*n.phe)),
                                    as.integer(n.col2drop),
                                    as.integer(col2drop),
                                    PACKAGE="qtl")
                        }

                        ## fill results matrix
                        if(n.phe == 1)
                            results[wh.col[[i]],wh.col[[i]]] <-
                                matrix(z$result,ncol=n.pos[i])
                        else # multiple phenotypes
                            results[wh.col[[i]],wh.col[[i]],] <-
                                array(z$result,c(n.pos[i],n.pos[i], n.phe))
                        z <- 0
                    }
                    else {
                        z <- .C("R_scantwo_1chr_em",
                                as.integer(n.ind),
                                as.integer(n.pos[i]),
                                as.integer(n.gen[i]),
                                as.double(temp),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.double(pheno),
                                as.double(weights),
                                result=as.double(rep(0,n.pos[i]^2)),
                                as.integer(maxit),
                                as.double(tol),
                                as.integer(verbose),
                                as.integer(n.col2drop),
                                as.integer(col2drop),
                                PACKAGE="qtl")
                        # re-organize results
                        results[wh.col[[i]],wh.col[[i]]] <-
                            matrix(z$result,ncol=n.pos[i])
                    }

                    z <- 0
                    temp <- 0 # remove the joint genotype probabilities
                } # end same chromosome

                else {
                    if(method=="hk") {
                        if(n.phe > batchsize) {
                            firstcol <- 1
                            z <- NULL
                            while(firstcol <= n.phe) {
                                thiscol <- firstcol + 0:(batchsize-1)
                                thiscol <- thiscol[thiscol <= n.phe]

                                thisz <- .C("R_scantwo_2chr_hk",
                                            as.integer(n.ind),
                                            as.integer(n.pos[i]),
                                            as.integer(n.pos[j]),
                                            as.integer(n.gen[i]),
                                            as.integer(n.gen[j]),
                                            as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                            as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                                            as.double(ac),
                                            as.integer(n.ac),
                                            as.double(ic),
                                            as.integer(n.ic),
                                            as.double(pheno[,thiscol]),
                                            as.integer(length(thiscol)),
                                            as.double(weights),
                                            full=as.double(rep(0,n.pos[i]*n.pos[j]*length(thiscol))),
                                            int=as.double(rep(0,n.pos[i]*n.pos[j]*length(thiscol))),
                                            PACKAGE="qtl")

                                firstcol <- firstcol + batchsize
                                if(is.null(z)) {
                                    z <- thisz
                                    z$full <- z$int <- array(NA, dim=c(n.pos[j], n.pos[i], n.phe))
                                }
                                z$full[,,thiscol] <- array(thisz$full, dim=c(n.pos[j], n.pos[i], length(thiscol)))
                                z$int[,,thiscol] <- array(thisz$int, dim=c(n.pos[j], n.pos[i], length(thiscol)))
                            }
                        }
                        else {
                            z <- .C("R_scantwo_2chr_hk",
                                    as.integer(n.ind),
                                    as.integer(n.pos[i]),
                                    as.integer(n.pos[j]),
                                    as.integer(n.gen[i]),
                                    as.integer(n.gen[j]),
                                    as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                    as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                                    as.double(ac),
                                    as.integer(n.ac),
                                    as.double(ic),
                                    as.integer(n.ic),
                                    as.double(pheno),
                                    as.integer(n.phe),
                                    as.double(weights),
                                    full=as.double(rep(0,n.pos[i]*n.pos[j]*n.phe)),
                                    int=as.double(rep(0,n.pos[i]*n.pos[j]*n.phe)),
                                    PACKAGE="qtl")
                        }
                        ## reorgnize results
                        if(n.phe == 1) {
                            results[wh.col[[j]],wh.col[[i]]] <-
                                matrix(z$full,ncol=n.pos[j])
                            results[wh.col[[i]],wh.col[[j]]] <-
                                matrix(z$int,ncol=n.pos[j])
                        }
                        else { # multiple phenotypes
                            results[wh.col[[j]],wh.col[[i]],] <-
                                array(z$full,c(n.pos[j], n.pos[i], n.phe))
                            results[wh.col[[i]],wh.col[[j]],] <-
                                array(z$int,c(n.pos[j], n.pos[i], n.phe))
                        }
                        z <- 0
                    }
                    else  {
                        z <- .C("R_scantwo_2chr_em",
                                as.integer(n.ind),
                                as.integer(n.pos[i]),
                                as.integer(n.pos[j]),
                                as.integer(n.gen[i]),
                                as.integer(n.gen[j]),
                                as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.double(pheno),
                                as.double(weights),
                                full=as.double(rep(0,n.pos[i]*n.pos[j])),
                                int=as.double(rep(0,n.pos[i]*n.pos[j])),
                                as.integer(maxit),
                                as.double(tol),
                                as.integer(verbose),
                                PACKAGE="qtl")

                        results[wh.col[[j]],wh.col[[i]]] <-
                            t(matrix(z$full,ncol=n.pos[j]))
                        results[wh.col[[i]],wh.col[[j]]] <-
                            matrix(z$int,ncol=n.pos[j])
                        z <- 0
                    }
                } # end different chromosome
            }
            else if(model=="binary") {
                if(do.nllik0) { # first time! do null model and get neg log10 likelihood
                    do.nllik0 <- FALSE
                    if(n.ac > 0)
                        nullfit <- glm(pheno ~ ac, family=binomial(link="logit"))
                    else
                        nullfit <- glm(pheno ~ 1, family=binomial(link="logit"))
                    fitted <- nullfit$fitted
                    nullcoef <- nullfit$coef
                    nllik0 <- -sum(pheno*log10(fitted) + (1-pheno)*log10(1-fitted))
                    if(verbose > 1) cat("null log lik: ", nllik0, "\n")
                }

                if(i==j) { # same chromosome

                    start <- c(rep(nullcoef[1],n.gen[i]),rep(0,n.gen[i]-1),
                               nullcoef[-1],rep(0,n.gen[i]*n.gen[i]+
                                                (n.gen[i]*(n.gen[i]-1)*n.ic)))

                    if(n.col2drop)
                        start <- c(start[!col2drop], rep(0,sum(col2drop)))

                    if(verbose>1) cat("  --Calculating joint probs.\n")

                    if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft"))) {
                        # calculate joint genotype probabilities for all pairs of positions
                        stp <- attr(oldXchr$geno[[1]]$prob, "step")
                        oe <- attr(oldXchr$geno[[1]]$prob, "off.end")
                        err <- attr(oldXchr$geno[[1]]$prob, "error.prob")
                        mf <- attr(oldXchr$geno[[1]]$prob, "map.function")

                        if("stepwidth" %in% names(attributes(oldXchr$geno[[1]]$prob)))
                            stpw <- attr(oldXchr$geno[[1]]$prob, "stepwidth")
                        else stpw <- "fixed"
                        if("map" %in% names(attributes(oldXchr$geno[[1]]$prob)))
                            tmap <- attr(oldXchr$geno[[1]]$prob,"map")
                        else
                            tmap <- create.map(oldXchr$geno[[1]]$map, stp, oe, stpw)

                        temp <- calc.pairprob(oldXchr,stp,oe,err,mf,tmap,
                                              assumeCondIndep=assumeCondIndep)
                    }
                    else {
                        # calculate joint genotype probabilities for all pairs of positions
                        stp <- attr(cross$geno[[i]]$prob, "step")
                        oe <- attr(cross$geno[[i]]$prob, "off.end")
                        err <- attr(cross$geno[[i]]$prob, "error.prob")
                        mf <- attr(cross$geno[[i]]$prob, "map.function")

                        if("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
                            stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
                        else stpw <- "fixed"
                        if("map" %in% names(attributes(cross$geno[[i]]$prob)))
                            tmap <- attr(cross$geno[[i]]$prob,"map")
                        else
                            tmap <- create.map(cross$geno[[i]]$map, stp, oe, stpw)

                        temp <- calc.pairprob(subset(cross,chr=thechr[i]),stp,oe,err,mf,tmap,
                                              assumeCondIndep=assumeCondIndep)
                    }

                    # pull out positions from genotype probs
                    if(some.dropped[i]) {
                        # figure out pos'ns corresponding to columns of temp
                        nc <- ncol(cross$geno[[i]]$prob)
                        ind <- matrix(rep(1:nc,nc),ncol=nc)
                        w <- lower.tri(ind)
                        ind <- cbind(first=t(ind)[w],second=ind[w])

                        # which part to keep
                        keep <- apply(ind,1,function(a,b) all(a %in% b),
                                      keep.pos[[i]])
                        temp <- temp[,keep,,]
                    }

                    # revise pair probilities for X chromosome
                    if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft")))
                        temp <- reviseXdata(type, "full", sexpgm, pairprob=temp,
                                            cross.attr=attributes(cross))

                    if(verbose>1) cat("  --Done.\n")

                    if(method=="em")
                        z <- .C("R_scantwo_1chr_binary_em",
                                as.integer(n.ind),
                                as.integer(n.pos[i]),
                                as.integer(n.gen[i]),
                                as.double(temp),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.integer(pheno),
                                as.double(start),
                                result=as.double(rep(0,n.pos[i]^2)),
                                as.integer(maxit),
                                as.double(tol),
                                as.integer(verbose),
                                as.integer(n.col2drop),
                                as.integer(col2drop),
                                PACKAGE="qtl")
                    else # h-k regression
                        z <- .C("R_scantwo_1chr_binary_hk",
                                as.integer(n.ind),
                                as.integer(n.pos[i]),
                                as.integer(n.gen[i]),
                                as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                as.double(temp),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.double(pheno),
                                result=as.double(rep(0,n.pos[i]^2)),
                                as.integer(n.col2drop),
                                as.integer(col2drop),
                                as.double(tol),
                                as.integer(maxit),
                                as.integer(verbose),
                                PACKAGE="qtl")


                    # re-organize results
                    results[wh.col[[i]],wh.col[[i]]] <-
                        matrix(z$result,ncol=n.pos[i])

                    z <- 0
                    temp <- 0 # remove the joint genotype probabilities
                } # end same chromosome
                else {
                    start <- c(rep(nullcoef[1],n.gen[i]),rep(0,n.gen[j]-1),
                               nullcoef[-1],rep(0,n.gen[i]*n.gen[j]+
                                                (n.gen[i]*(n.gen[j]-1)*n.intcovar)));

                    if(method=="em")
                        z <- .C("R_scantwo_2chr_binary_em",
                                as.integer(n.ind),
                                as.integer(n.pos[i]),
                                as.integer(n.pos[j]),
                                as.integer(n.gen[i]),
                                as.integer(n.gen[j]),
                                as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.integer(pheno),
                                as.double(start),
                                full=as.double(rep(0,n.pos[i]*n.pos[j])),
                                int=as.double(rep(0,n.pos[i]*n.pos[j])),
                                as.integer(maxit),
                                as.double(tol),
                                as.integer(verbose),
                                PACKAGE="qtl")
                    else # h-k regression
                        z <- .C("R_scantwo_2chr_binary_hk",
                                as.integer(n.ind),
                                as.integer(n.pos[i]),
                                as.integer(n.pos[j]),
                                as.integer(n.gen[i]),
                                as.integer(n.gen[j]),
                                as.double(cross$geno[[i]]$prob[,keep.pos[[i]],]),
                                as.double(cross$geno[[j]]$prob[,keep.pos[[j]],]),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.double(pheno),
                                full=as.double(rep(0,n.pos[i]*n.pos[j])),
                                int=as.double(rep(0,n.pos[i]*n.pos[j])),
                                as.double(tol),
                                as.integer(maxit),
                                as.integer(verbose),
                                PACKAGE="qtl")

                    results[wh.col[[j]],wh.col[[i]]] <-
                        t(matrix(z$full,ncol=n.pos[j]))
                    results[wh.col[[i]],wh.col[[j]]] <-
                        matrix(z$int,ncol=n.pos[j])
                    z <- 0
                } # end same chromosome
            }
            else { # marker regression
                # replace missing and partially informative genotypes with 0's
                datai <- cross$geno[[i]]$data
                datai[is.na(datai)] <- 0

                if(type=="f2" || (type=="bcsft" & !is.bcs)) datai[datai>3] <- 0
                else if(type=="4way") datai[datai>4] <- 0

                if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft")))
                    datai <- reviseXdata(type, "full", sexpgm, geno=datai,
                                         cross.attr=attributes(cross))

                if(i==j) { # same chromosome

                    z <- .C("R_scantwo_1chr_mr",
                            as.integer(n.ind),
                            as.integer(n.pos[i]),
                            as.integer(n.gen[i]),
                            as.integer(datai),
                            as.double(ac),
                            as.integer(n.ac),
                            as.double(ic),
                            as.integer(n.ic),
                            as.double(pheno),
                            as.double(weights),
                            result=as.double(rep(0,n.pos[i]^2)),
                            as.integer(n.col2drop),
                            as.integer(col2drop),
                            PACKAGE="qtl")

                    # re-organize results
                    results[wh.col[[i]],wh.col[[i]]] <-
                        matrix(z$result,ncol=n.pos[i])
                    z <- 0
                } # end same chromosome
                else {

                    # replace missing and partially informative genotypes with 0's
                    dataj <- cross$geno[[j]]$data
                    dataj[is.na(dataj)] <- 0
                    if(type=="f2" || (type=="bcsft" && !is.bcs)) dataj[dataj>3] <- 0
                    else if(type=="4way") dataj[dataj>4] <- 0

                    if(chrtype[j]=="X" && (type %in% c("bc","f2","bcsft")))
                        dataj <- reviseXdata(type, "full", sexpgm, geno=dataj,
                                             cross.attr=attributes(cross))

                    z <- .C("R_scantwo_2chr_mr",
                            as.integer(n.ind),
                            as.integer(n.pos[i]),
                            as.integer(n.pos[j]),
                            as.integer(n.gen[i]),
                            as.integer(n.gen[j]),
                            as.integer(datai),
                            as.integer(dataj),
                            as.double(ac),
                            as.integer(n.ac),
                            as.double(ic),
                            as.integer(n.ic),
                            as.double(pheno),
                            as.double(weights),
                            full=as.double(rep(0,n.pos[i]*n.pos[j])),
                            int=as.double(rep(0,n.pos[i]*n.pos[j])),
                            PACKAGE="qtl")

                    results[wh.col[[j]],wh.col[[i]]] <-
                        t(matrix(z$full,ncol=n.pos[j]))
                    results[wh.col[[i]],wh.col[[j]]] <-
                        matrix(z$int,ncol=n.pos[j])

                    z <- 0
                }
            }

        } # end loop over second chr
    } # end loop over first chromosome

    # subtract null neg log lik from lower tri
    if(method=="em") {
        offdiag <- lower.tri(results) | upper.tri(results)
        results[offdiag] <- nllik0 - results[offdiag]
    }
    else if(method=="hk") {
        if(n.phe == 1) {
            offdiag <- lower.tri(results) | upper.tri(results)
            results[offdiag] <- nllik0 - results[offdiag]
        }
        else { # multiple phenotypes
            offdiag <- lower.tri(results[,,1]) | upper.tri(results[,,1])
            for(itmp in 1:n.phe) {
                # I'm doing a loop here. I should put null model back to C function
                results[,,itmp][offdiag] <- nllik0[itmp] - results[,,itmp][offdiag]
            }
        }
    }

    # If the X chromosome was included, need to do an adjustment...
    scanoneX <- NULL
    if(any(gmap[,4])) { # the X chromosome was included

        # determine which covariates belong in null hypothesis
        temp <- scanoneXnull(type, sexpgm, cross.attr=attributes(cross))
        adjustX <- temp$adjustX
        parX0 <- temp$parX0
        sexpgmcovar <- temp$sexpgmcovar

        if(adjustX) {
            if(method=="mr" && any(is.na(pull.geno(cross))))
                if(n.perm > -2) warning("Scantwo with the X chr doesn't work quite right when method=\"mr\"\n",
                                        "    when there is missing genotype data.")

            if(model=="binary") {
                if(n.ac > 0) {
                    nullfitX <- glm(pheno ~ ac+sexpgmcovar, family=binomial(link="logit"))
                    parX0 <- lm(pheno ~ ac+sexpgmcovar)$rank
                }
                else {
                    nullfitX <- glm(pheno ~ sexpgmcovar, family=binomial(link="logit"))
                    parX0 <- ncol(sexpgmcovar)
                }
                fittedX <- nullfitX$fitted
                nullcoefX <- nullfitX$coef
                nllikX <- -sum(pheno*log10(fittedX) + (1-pheno)*log10(1-fittedX))
                if(verbose > 1) cat("X chr null log lik: ", nllikX, "\n")
            }
            else {
                if(n.ac > 0) {
                    outX <- lm(pheno ~ ac+sexpgmcovar, weights=weights^2)
                    residX <- outX$resid
                    parX0 <- outX$rank
                }
                else
                    residX <- lm(pheno ~ sexpgmcovar, weights=weights^2)$resid

                if(method=="hk") {
                    if(n.phe==1)
                        nllikX <- (n.ind/2)*log10(sum((residX*weights)^2))
                    else
                        nllikX <- (n.ind/2)*apply(residX, 2, function(a,b)
                                                  log10(sum((a*b)^2)), weights)
                }
                else {
                    if(method=="imp" || method=="mr") {
                        if(n.ac > 0) {
                            out0 <- lm(pheno ~ ac, weights=weights^2)
                            resid0 <- out0$resid
                        }
                        else {
                            out0 <- lm(pheno ~ 1, weights=weights^2)
                            resid0 <- out0$resid
                        }

                        if(n.phe > 1) {
                            sig0 <- sqrt(apply(resid0, 2, function(a,b) sum((a*b)^2),weights)/n.ind)
                            nllik0 <- sig0
                            for(i in seq(along=nllik0))
                                nllik0[i] <- -sum(dnorm(resid0[,i],0,sig0[i]/weights,log=TRUE))/log(10)
                        }
                        else {
                            sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
                            nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
                        }
                    }

                    if(n.phe > 1) {
                        sigX <- sqrt(apply(residX, 2, function(a,b) sum((a*b)^2),weights)/n.ind)
                        nllikX <- sigX
                        for(i in seq(along=nllikX))
                            nllikX[i] <- -sum(dnorm(residX[,i],0,sigX[i]/weights,log=TRUE))/log(10)
                    }
                    else {
                        sigX <- sqrt(sum((residX*weights)^2)/n.ind)
                        nllikX <- -sum(dnorm(residX,0,sigX/weights,log=TRUE))/log(10)
                    }
                }
            }

            if(n.phe > 1)
                wh <- ((gmap[row(results[,,1]),4] | gmap[col(results[,,1]),4]) &
                       (lower.tri(results[,,1]) | upper.tri(results[,,1])))
            else
                wh <- ((gmap[row(results),4] | gmap[col(results),4]) &
                       (lower.tri(results) | upper.tri(results)))

            if(n.phe > 1) {
                for(i in 1:n.phe)
                    results[,,i][wh] <- results[,,i][wh] + nllikX[i] - nllik0[i]
            }
            else
                results[wh] <- results[wh] + nllikX - nllik0

            notxchr <- names(cross$geno)[sapply(cross$geno,class)!="X"]
            if(length(notxchr) > 0) {
                if(verbose) cat(" --Running scanone with special X chr covariates\n")
                temp <- scanone(subset(cross,chr=notxchr),
                                pheno.col=pheno.col,
                                model=model, method=method,
                                addcovar=cbind(ac,sexpgmcovar),
                                intcovar=ic, weights=weights, use=use,
                                maxit=maxit, tol=tol, verbose=FALSE)

                scanoneX <- temp[,-(1:2),drop=FALSE]
                scanoneX <- rbind(scanoneX,
                                  out.scanone[rownames(gmap),,drop=FALSE][gmap[,4],,drop=FALSE])

                scanoneX <- scanoneX[rownames(gmap),,drop=FALSE]
            }
            else {
                scanoneX <- out.scanone[rownames(gmap),,drop=FALSE][gmap[,4],,drop=FALSE]
                scanoneX <- scanoneX[rownames(gmap),,drop=FALSE]
            }
        }

    }

    #  if(any(is.na(results)) && n.perm > -2)
    #    warning("Some LOD scores NA")
    #  if(any(!is.na(results) & results < 0) && n.perm > -2)
    #    warning("Some LOD scores < 0")
    #  if(any(!is.na(results) & (results == Inf | results == -Inf)) && n.perm > -2)
    #    warning("Some LOD scores = Inf or -Inf")

    if(!is.null(scanoneX)) scanoneX <- as.matrix(scanoneX)

    # output has 2 fields, lod and map
    out <- list(lod=results,map=gmap,scanoneX=scanoneX)

    # fill in scanone result
    if(n.phe == 1)
        diag(out$lod) <- out.scanone[rownames(out$map),]
    else {
        for(iphe in 1:n.phe) {
            if(nrow(out$lod)==1)
                out$lod[1,1,iphe] <- out.scanone[rownames(out$map),iphe]
            else
                diag(out$lod[,,iphe]) <- out.scanone[rownames(out$map),iphe]
        }
    }

    attr(out,"method") <- method
    attr(out,"type") <- type
    attr(out, "fullmap") <- fullmap
    class(out) <- "scantwo"

    if(clean.output) # remove NA, 0 out positions between markers
        out <- clean(out, clean.nmar, clean.distance)

    attr(out, "phenotypes") <- colnames(pheno)
    if(length(colnames(pheno)) > 1)
        dimnames(out$lod) <- list(NULL, NULL, colnames(pheno))
    names(out$map)[2] <- "pos"
    out
}

######################################################################
#
# scantwo.perm: Permutation test of scantwo
#
######################################################################

scantwo.perm <-
    function(cross, pheno.col=1, model=c("normal","binary"),
             method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"),
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=4000, tol=1e-4, verbose=FALSE,
             n.perm=1000, perm.strata,
             assumeCondIndep=FALSE, batchsize=250, chr)
{
    method <- match.arg(method)
    model <- match.arg(model)
    use <- match.arg(use)
    if(missing(chr)) chr <- names(chr$geno)

    scantwo.perm.engine(n.perm, cross=cross, pheno.col=pheno.col,
                        model=model, method=method, addcovar=addcovar,
                        intcovar=intcovar, weights=weights, use=use,
                        incl.markers=incl.markers, clean.output=clean.output,
                        clean.nmar=clean.nmar,
                        clean.distance=clean.distance,
                        maxit=maxit, tol=tol, verbose=verbose,
                        perm.strata=perm.strata,
                        assumeCondIndep=assumeCondIndep, batchsize=batchsize,
                        chr=chr)

}


######################################################################
#
# Engine function to scantwo permutation
#
######################################################################
scantwo.perm.engine <-
    function(n.perm, cross, pheno.col, model,
             method, addcovar, intcovar, weights, use,
             incl.markers, clean.output, clean.nmar=1, clean.distance=0,
             maxit, tol, verbose, perm.strata,
             assumeCondIndep=FALSE, batchsize=250, chr)
{
    if(missing(chr)) chr <- names(chr$geno)

    ## local variables
    n.phe <- length(pheno.col)
    n.ind <- dim(cross$pheno)[1]
    pn <- colnames(cross$pheno)[pheno.col]

    ## if there's only one phenotype, no covariate, and method is imp or hk,
    ## generate permuted phenotype as a matrix and do permutation
    ## as multiple phenotypes
    ## we also need one sex and one direction, or that the
    ##     stratification is within those groups
    batch.mode <- FALSE
    if( (n.phe==1) && ((method=="imp") || (method=="hk")) &&
       model=="normal" &&
       is.null(addcovar) && is.null(intcovar) ) {
        chrtype <- sapply(cross$geno, class)
        sexpgm <- getsex(cross)
        sex <- sexpgm$sex
        pgm <- sexpgm$pgm
        if(all(chrtype=="A"))
            batch.mode <- TRUE
        else if((is.null(sex) || length(unique(sex))==1) &&
                (is.null(pgm) || length(unique(pgm))==1))
            batch.mode <- TRUE
        else if(!is.null(perm.strata)) {
            sp <- paste(sex, pgm, sep=":")
            tab <- table(sp, perm.strata)
            if(all(apply(tab, 2, function(a) sum(a != 0))==1))
                batch.mode <- TRUE
        }
    }

    if(batch.mode) {
        if(verbose)
            cat("Doing permutation in batch mode ...\n")
        ord <- matrix(0, n.ind, n.perm)
        if(!is.null(perm.strata)) {  # stratified permutation test
            if(length(perm.strata) != n.ind)
                stop("perm.strata must be NULL or have length nind(cross).")
            u <- unique(perm.strata)
            theindex <- 1:n.ind
            if(length(u)==n.ind)
                stop("All elements of perm.strata are unique, so there will be no real permutation.")
            if(length(u)==1)
                warning("Just one unique element in perm.strata, so the perms are not stratified.")

            for(iperm in 1:n.perm) {
                for(j in u) {
                    wh <- perm.strata==j
                    if(sum(wh)==1) ord[wh,iperm] <- theindex[wh]
                    else ord[wh,iperm] <- sample(theindex[wh])
                }
            }
        }
        else {
            for(iperm in 1:n.perm)
                ord[,iperm] <- sample(n.ind)
        }

        cross$pheno <- cbind(matrix(cross$pheno[,pheno.col][ord], nrow=n.ind), cross$pheno)

        pheno.col <- 1:n.perm

        if(is.list(chr)) {
            chr1 <- chr[[1]]
            chr2 <- chr[[2]]
        }
        else chr1 <- chr2 <- chr

        thechr <- names(cross$geno)
        nchr1 <- match(chr1, thechr)
        nchr2 <- match(chr2, thechr)

        perm.result <- NULL
        for(i in nchr1) {
            for(j in nchr2) {
                if(j < i) next

                tem <- scantwo(cross, pheno.col=pheno.col, model=model, method=method,
                               addcovar=addcovar, intcovar=intcovar, weights=weights,
                               use=use, incl.markers=incl.markers,
                               clean.output=clean.output, clean.nmar=clean.nmar,
                               clean.distance=clean.distance,
                               maxit=maxit, tol=tol,verbose=FALSE, n.perm=-1,
                               perm.strata=perm.strata,
                               assumeCondIndep=assumeCondIndep,
                               batchsize=batchsize, n.cluster=0, chr=list(thechr[i],thechr[j]))

                if(clean.output) tem <- clean(tem, clean.nmar, clean.distance)

                ## find the maximum LOD on each permutation
                if(is.null(perm.result)) {
                    perm.result <- lapply(subrousummaryscantwo(tem,for.perm=TRUE), as.matrix)
                }
                else {
                    tem <- lapply(subrousummaryscantwo(tem,for.perm=TRUE), as.matrix)
                    for(k in seq(along=perm.result))
                        perm.result[[k]] <- as.matrix(apply(cbind(perm.result[[k]], tem[[k]]), 1, max, na.rm=TRUE))
                }


            }
        }
    }
    else { ## all other cases, do one permutation at a time
        if(method=="mr-imp") # save version with missing genotypes
            tempcross <- cross
        if(method=="mr-argmax") # impute genotypes
            cross <- fill.geno(cross,method="argmax")
        if(!is.null(addcovar)) addcovar <- as.matrix(addcovar)
        if(!is.null(intcovar)) intcovar <- as.matrix(intcovar)
        addcovarp <- addcovar
        intcovarp <- intcovar


        ## initialize result
        temp <- matrix(ncol=n.phe, nrow=n.perm)
        perm.result <- list("full"=temp,
                            "fv1"=temp,
                            "int"=temp,
                            "add"=temp,
                            "av1"=temp,
                            "one"=temp)

        if(is.list(chr)) {
            chr1 <- chr[[1]]
            chr2 <- chr[[2]]
        }
        else chr1 <- chr2 <- chr

        thechr <- names(cross$geno)
        nchr1 <- match(chr1, thechr)
        nchr2 <- match(chr2, thechr)

        ## permutation loop
        for(i in 1:n.perm) {
            if(verbose) cat("Permutation", i, "\n")
            ## impute genotypes for method "mr-imp"
            if(method=="mr-imp") cross <- fill.geno(tempcross)

            if(!is.null(perm.strata)) {  # stratified permutation test
                if(length(perm.strata) != n.ind)
                    stop("perm.strata must be NULL or have length nind(cross).")
                u <- unique(perm.strata)
                theindex <- 1:n.ind
                if(length(u)==n.ind)
                    stop("All elements of perm.strata are unique, so no real permutations.")
                if(length(u)==1 && i==1)
                    warning("Just one unique element in perm.strata, so the perms are not stratified.")

                o <- 1:n.ind
                for(j in u) {
                    wh <- perm.strata==j
                    if(sum(wh)>1) o[wh] <- sample(o[wh])
                }
            }
            else
                o <- sample(1:n.ind)

            cross$pheno <- cross$pheno[o,,drop=FALSE]
            if(!is.null(addcovar)) addcovarp <- addcovarp[o,,drop=FALSE]
            if(!is.null(intcovar)) intcovarp <- intcovarp[o,,drop=FALSE]

            temp <- NULL
            for(ii in nchr1) {
                for(jj in nchr2) {
                    if(jj < ii) next

                    tem <- scantwo(cross, pheno.col=pheno.col, model=model, method=method,
                                   addcovar=addcovar, intcovar=intcovar, weights=weights,
                                   use=use, incl.markers=incl.markers,
                                   clean.output=clean.output, clean.nmar=clean.nmar,
                                   clean.distance=clean.distance,
                                   maxit=maxit, tol=tol,verbose=FALSE, n.perm=-i,
                                   perm.strata=perm.strata,
                                   assumeCondIndep=assumeCondIndep,
                                   batchsize=batchsize, n.cluster=0, chr=list(thechr[ii],thechr[jj]))

                    if(clean.output) tem <- clean(tem, clean.nmar, clean.distance)

                    ## find the maximum LOD on each permutation
                    if(is.null(temp)) {
                        temp <- lapply(subrousummaryscantwo(tem,for.perm=TRUE), as.matrix)
                    }
                    else {
                        tem <- lapply(subrousummaryscantwo(tem,for.perm=TRUE), as.matrix)
                        for(k in seq(along=temp))
                            temp[[k]] <- as.matrix(apply(cbind(temp[[k]], tem[[k]]), 1, max, na.rm=TRUE))
                    }

                }
            }


            # maxima
            for(j in 1:6) perm.result[[j]][i,] <- temp[[j]]

        }
    }

    ## make result
    attr(perm.result,"method") <- method
    class(perm.result) <- c("scantwoperm", "list")

    ## add column names
    for(i in 1:length(perm.result))
        colnames(perm.result[[i]]) <- pn

    perm.result
}


# end of scantwo.R
