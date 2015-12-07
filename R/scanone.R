#####################################################################
#
# scanone.R
#
# copyright (c) 2001-2015, Karl W Broman
# last modified Oct, 2015
# first written Feb, 2001
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
# Hao Wu (The Jackson Lab) wrote the imputation method
#
# Part of the R/qtl package
# Contains: scanone, scanone.perm, scanone.perm.engine
#
######################################################################

######################################################################
#
# scanone: scan genome, calculating LOD scores with single QTL model
#          (covariates are not allowed for models other than "normal"
#           and "binary")
#
######################################################################

scanone <-
    function(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"), upper=FALSE,
             ties.random=FALSE, start=NULL, maxit=4000, tol=1e-4,
             n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose, batchsize=250,
             n.cluster=1, ind.noqtl)
{
    if(batchsize < 1) stop("batchsize must be >= 1.")
    model <- match.arg(model)
    method <- match.arg(method)
    use <- match.arg(use)

    # in RIL, treat X chromomse like an autosome
    chrtype <- sapply(cross$geno, class)
    if(any(chrtype=="X") && (class(cross)[1] == "risib" || class(cross)[1] == "riself"))
        for(i in which(chrtype=="X")) class(cross$geno[[i]]) <- "A"

    if(!missing(n.perm) && n.perm > 0 && n.cluster > 1) {
        cat(" -Running permutations via a cluster of", n.cluster, "nodes.\n")
        updateParallelRNG(n.cluster)

        scanonePermInParallel <- function(n.perm, cross, chr, pheno.col, model, method, addcovar,
                                          intcovar, weights, use, upper, ties.random, start, maxit, tol,
                                          perm.Xsp, perm.strata, batchsize)
            scanone(cross=cross, chr=chr, pheno.col=pheno.col, model=model, method=method, addcovar=addcovar,
                    intcovar=intcovar, weights=weights, use=use, upper=upper, ties.random=ties.random, start=start,
                    maxit=maxit, tol=tol, n.perm=n.perm, perm.Xsp=perm.Xsp, perm.strata=perm.strata, batchsize=batchsize,
                    n.cluster=0, verbose=FALSE)

        n.perm <- ceiling(n.perm/n.cluster)
        if(missing(chr)) chr <- names(cross$geno)
        if(Sys.info()[1] == "Windows") { # Windows doesn't support mclapply, but it's faster if available
            cl <- makeCluster(n.cluster)
            on.exit(stopCluster(cl))
            operm <- clusterApply(cl, rep(n.perm, n.cluster), scanonePermInParallel, cross=cross, chr=chr, pheno.col=pheno.col,
                                  model=model, method=method, addcovar=addcovar, intcovar=intcovar, weights=weights, use=use,
                                  upper=upper, ties.random=ties.random, start=start, maxit=maxit, tol=tol, perm.Xsp=perm.Xsp,
                                  perm.strata=perm.strata, batchsize=batchsize)
        }
        else {
            operm <- mclapply(rep(n.perm, n.cluster), scanonePermInParallel, cross=cross, chr=chr, pheno.col=pheno.col,
                              model=model, method=method, addcovar=addcovar, intcovar=intcovar, weights=weights, use=use,
                              upper=upper, ties.random=ties.random, start=start, maxit=maxit, tol=tol, perm.Xsp=perm.Xsp,
                              perm.strata=perm.strata, batchsize=batchsize, mc.cores=n.cluster)
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

    # individuals with no QTL effect
    if(missing(ind.noqtl)) ind.noqtl <- rep(FALSE, nind(cross))
    else {
        if(method %in% c("mr", "mr-imp", "mr-argmax","ehk") ||
           model %in% c("2part", "np")) {
            ind.noqtl <- rep(FALSE, nind(cross))
            warning("ind.noqtl ignored for model=", model, ", method=", method)
        }
        else if(is.null(addcovar) && (!is.logical(ind.noqtl) || any(ind.noqtl))) {
            ind.noqtl <- rep(FALSE, nind(cross))
            warning("ind.noqtl ignored when no additive covariates")
        }
        else if(!is.logical(ind.noqtl) || length(ind.noqtl) != nind(cross))
            stop("ind.noqtl be a logical vector of length n.ind (", nind(cross), ")")
    }

    if(!missing(chr)) cross <- subset(cross, chr)
    if(missing(n.perm)) n.perm <- 0

    if(missing(verbose)) {
        if(!missing(n.perm) && n.perm > 0) verbose <- TRUE
        else verbose <- FALSE
    }

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
                               perm.strata, ind.noqtl, weights, TRUE)
            cross <- temp[[1]]
            pheno <- temp[[2]]
            addcovar <- temp[[3]]
            intcovar <- temp[[4]]
            n.addcovar <- temp[[5]]
            n.intcovar <- temp[[6]]
            perm.strata <- temp[[7]]
            ind.noqtl <- temp[[8]]
            weights <- temp[[9]]
        }
    }

    # use all observations; not in a permutation test; different phenotypes have different sets of missing values
    #   -> want to do in batches, but need to define batches by the pattern of missing data
    if(n.perm <= 0 && use=="all.obs" && length(pheno.col) > 1 && (method=="hk" || method=="imp") && model == "normal") {
        # drop individuals with missing covariates
        cross$pheno <- cbind(cross$pheno, rep(1, nind(cross)))
        temp <- checkcovar(cross, nphe(cross), addcovar, intcovar,
                           perm.strata, ind.noqtl, weights, TRUE)
        cross <- temp[[1]]
        pheno <- cross$pheno[,pheno.col, drop=FALSE]
        addcovar <- temp[[3]]
        intcovar <- temp[[4]]
        n.addcovar <- temp[[5]]
        n.intcovar <- temp[[6]]
        perm.strata <- temp[[7]]
        ind.noqtl <- temp[[8]]
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

            temp <- scanone(subset(cross, ind=upat[[i]]), chr=chr, pheno.col=batches[[i]], model=model,
                            method=method, addcovar=tempac, intcovar=tempic,
                            weights=weights[upat[[i]]], use="complete.obs", upper=upper, ties.random=ties.random,
                            start=start, maxit=maxit, tol=tol, n.perm=n.perm, perm.Xsp=perm.Xsp,
                            perm.strata=perm.strata[upat[[i]]], verbose=verbose, batchsize=batchsize,
                            n.cluster=n.cluster)
            if(is.null(out)) out <- temp
            else out <- cbind(out, temp)
        }

        # reorder LOD score columns and make sure that the names are correct
        colnames(out)[-(1:2)] <- colnames(cross$pheno)[unlist(batches)]
        out[,-(1:2)] <- out[,colnames(cross$pheno)[pheno.col]]
        colnames(out)[-(1:2)] <- colnames(cross$pheno)[pheno.col]

        return(out)
    }

    # multiple phenotype for methods except imp or hk with normal model
    if(length(pheno.col)>1 && n.perm <= 0 && (model != "normal" ||
                              (method!="imp" && method != "hk"))) {
        out <- scanone(cross, chr, pheno.col[1], model, method,
                       addcovar, intcovar, weights, use, upper, ties.random,
                       start, maxit, tol, n.perm, perm.Xsp, perm.strata,
                       verbose, batchsize, n.cluster=0)
        nc <- ncol(out)-2
        cn <- colnames(out)[-(1:2)]
        for(i in 2:length(pheno.col))
            out[,ncol(out)+1:nc] <- scanone(cross, chr, pheno.col[i], model,
                                            method, addcovar, intcovar, weights,
                                            use, upper, ties.random, start,
                                            maxit, tol, n.perm, perm.Xsp,
                                            perm.strata, verbose, batchsize, n.cluster=0)[,-(1:2)]

        if(length(cn) > 1)
            colnames(out)[-(1:2)] <- paste(rep(cn,length(pheno.col)),
                                           rep(colnames(cross$pheno)[pheno.col],
                                               rep(nc,length(pheno.col))),
                                           sep=".")
        else
            colnames(out)[-(1:2)] <- colnames(cross$pheno)[pheno.col]

        return(out)
    }

    if(n.perm>0) {
        return(scanone.perm(cross, pheno.col, model, method, addcovar,
                            intcovar, weights, use, upper, ties.random,
                            start, maxit, tol, n.perm, perm.Xsp, perm.strata,
                            verbose, batchsize))
    }

    if(n.perm < 0) { # in the midst of permutations
        if(use=="all.obs") {
            temp <- checkcovar(cross, pheno.col, addcovar, intcovar,
                               perm.strata, ind.noqtl, weights, n.perm==-1)
            cross <- temp[[1]]
            pheno <- temp[[2]]
            addcovar <- temp[[3]]
            intcovar <- temp[[4]]
            n.addcovar <- temp[[5]]
            n.intcovar <- temp[[6]]
            perm.strata <- temp[[7]]
            ind.noqtl <- temp[[8]]
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
    n.phe <- ncol(pheno)
    type <- class(cross)[1]
    is.bcs <- FALSE
    if(type == "bcsft") {
        cross.scheme <- attr(cross, "scheme")
        is.bcs <- (cross.scheme[2] == 0)
    }

    # fill in missing genotypes with imputed values
    if(n.perm==0) { # not in the midst of permutations
        if(method=="mr-argmax")
            cross <- fill.geno(cross,method="argmax")
        if(method=="mr-imp")
            cross <- fill.geno(cross,method="imp")
    }

    # weights for model="normal"
    if(model != "normal") {
        if(!is.null(weights) && !all(weights==1) && n.perm > -2)
            warning("weights used only for normal model.")
    }
    else {
        if(is.null(weights))
            weights <- rep(1, nind(cross))
        else
            if(length(weights) != nind(cross))
                stop("weights should either be NULL or a vector of length n.ind")
        if(any(weights <= 0))
            stop("weights should be entirely positive")
        weights <- sqrt(weights)
    }

    if(model=="binary") {
        if(method=="imp") {
            if(n.perm > -2) warning("Method imp not available for binary model; using em")
            method <- "em"
        }
        return(discan(cross, pheno, method, addcovar, intcovar, maxit, tol, verbose, n.perm > -2, ind.noqtl))
    }
    else if(model=="2part") {
        if((n.addcovar > 0 || n.intcovar > 0) && n.perm > -2)
            warning("Covariates ignored for the two-part model.")
        if(method!="em") {
            if(n.perm > -2) warning("Only em method is available for the two-part model")
            method <- "em"
        }
        return(vbscan(cross, pheno.col, upper, method, maxit, tol))
    }
    else if(model=="np") {
        if((n.addcovar > 0 || n.intcovar > 0) && n.perm > -2)
            warning("Covariates ignored for non-parametric interval mapping.")
        if(method!="em") {
            if(n.perm > -2) warning("Method argument ignored for non-parametric interval mapping.")
            method <- "em"
        }
    }

    # if non-parametric, convert phenotypes to ranks
    if(model=="np") {
        if(ties.random) {
            y <- pheno[!is.na(pheno)]
            y <- rank(y+runif(length(y))/(sd(y)*10^8))
            pheno[!is.na(pheno)] <- y
            correct <- 1
        }
        else {
            ties <- table(pheno)
            if(any(ties > 1)) {
                ties <- ties[ties>1]
                correct <- 1-sum(ties^3-ties)/(n.ind^3-n.ind)
            }
            else correct <- 1
            pheno <- rank(pheno)
        }
    }

    results <- NULL

    # starting points for interval mapping
    if(method=="em" && model=="normal") {
        if(is.null(start)) std.start <- 1
        else if(length(start)==1) std.start <- -1
        else std.start <- 0
    }

    # scan genome one chromosome at a time
    for(i in 1:n.chr) {
        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") {
            sexpgm <- getsex(cross)
            ac <- revisecovar(sexpgm,addcovar)

            if(!is.null(addcovar) && (nd <- attr(ac, "n.dropped")) > 0 && n.perm > -2)
                warning("Dropped ", nd, " additive covariates on X chromosome.")

            if(length(ac)==0) {
                n.ac <- 0
                ac <- NULL
            }
            else n.ac <- ncol(ac)
            ic <- revisecovar(sexpgm,intcovar)
            if(!is.null(intcovar) && (nd <- attr(ic, "n.dropped")) > 0 && n.perm > -2)
                warning("Dropped ", nd, " interactive covariates on X chromosome.")
            if(length(ic)==0) {
                n.ic <- 0
                ic <- NULL
            }
            else n.ic <- ncol(ic)
        }
        else {
            sexpgm <- NULL
            ac <- addcovar
            n.ac <- n.addcovar
            ic <- intcovar
            n.ic <- n.intcovar
        }

        # get genotype names
        gen.names <- getgenonames(type,chrtype,"full",sexpgm,attributes(cross))
        n.gen <- length(gen.names)

        # starting values for interval mapping
        if(method=="em" && model=="normal") {
            this.start <- rep(0,n.gen+1)
            if(std.start == 0) {
                if(length(start) < n.gen+1)
                    stop("Length of start argument should be 0, 1 or ", n.gen+1)
                this.start <- c(start[1:n.gen],start[length(start)])
            }
        }

        # pull out reconstructed genotypes (mr)
        # or imputations (imp)
        # or genotype probabilities (em or hk)
        if(method=="mr" || method=="mr-imp" || method=="mr-argmax") {
            newgeno <- cross$geno[[i]]$data
            newgeno[is.na(newgeno)] <- 0

            # discard partially informative genotypes
            if(type=="f2" || (type=="bcsft" && !is.bcs)) newgeno[newgeno>3] <- 0
            if(type=="4way") newgeno[newgeno>4] <- 0

            # revise X chromosome genotypes
            if(chrtype=="X" && (type %in% c("bc","f2","bcsft")))
                newgeno <- reviseXdata(type, "full", sexpgm, geno=newgeno,
                                       cross.attr=attributes(cross))

            n.pos <- ncol(newgeno)
            map <- cross$geno[[i]]$map
            if(is.matrix(map)) {
                marnam <- colnames(map)
                map <- map[1,]
            }
            else marnam <- names(map)
        }
        else if(method == "imp") {
            if(!("draws" %in% names(cross$geno[[i]]))) {
                # need to run sim.geno
                if(n.perm > -2) warning("First running sim.geno.")
                cross <- sim.geno(cross)
            }

            draws <- cross$geno[[i]]$draws
            n.pos <- ncol(draws)
            n.draws <- dim(draws)[3]

            # revise X chromosome genotypes
            if(chrtype=="X" && (type %in% c("bc","f2","bcsft")))
                draws <- reviseXdata(type, "full", sexpgm, draws=draws,
                                     cross.attr=attributes(cross))

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

            if(is.matrix(map)) {
                marnam <- colnames(map)
                map <- map[1,]
            }
            else marnam <- names(map)
        }
        else {
            if(!("prob" %in% names(cross$geno[[i]]))) {
                # need to run calc.genoprob
                if(n.perm > -2) warning("First running calc.genoprob.")
                cross <- calc.genoprob(cross)
            }
            genoprob <- cross$geno[[i]]$prob
            n.pos <- ncol(genoprob)

            # revise X chromosome genotypes
            if(chrtype=="X" && (type %in% c("bc","f2","bcsft")))
                genoprob <- reviseXdata(type, "full", sexpgm, prob=genoprob,
                                        cross.attr=attributes(cross))

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

            if(is.matrix(map)) {
                marnam <- colnames(map)
                map <- map[1,]
            }
            else marnam <- names(map)
        }

        # call the C function
        if(method == "mr" || method=="mr-imp" || method=="mr-argmax")
            z <- .C("R_scanone_mr",
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.gen),         # number of possible genotypes
                    as.integer(newgeno),       # genotype data
                    as.double(ac),       # additive covariates
                    as.integer(n.ac),
                    as.double(ic),       # interactive covariates
                    as.integer(n.ic),
                    as.double(pheno),          # phenotype data
                    as.double(weights),        # weights
                    result=as.double(rep(0,n.pos)),
                    PACKAGE="qtl")
        else if(method=="imp") {
            if(n.phe > batchsize) {
                firstcol <- 1
                z <- NULL
                while(firstcol <= n.phe) {
                    if(verbose > 2) cat("chr", names(cross$geno)[i], "phe", firstcol, "...\n")
                    thiscol <- firstcol + 0:(batchsize-1)
                    thiscol <- thiscol[thiscol <= n.phe]
                    thisz <- .C("R_scanone_imp",
                                as.integer(n.ind),
                                as.integer(n.pos),
                                as.integer(n.gen),
                                as.integer(n.draws),
                                as.integer(draws),
                                as.double(ac),
                                as.integer(n.ac),
                                as.double(ic),
                                as.integer(n.ic),
                                as.double(pheno[,thiscol]),
                                as.integer(length(thiscol)), # number of phenotypes
                                as.double(weights),
                                result=as.double(rep(0,length(thiscol)*n.pos)),
                                as.integer(ind.noqtl),       # indicators of ind'l w/o QTL effects
                                PACKAGE="qtl")
                    firstcol <- firstcol + batchsize
                    if(is.null(z)) {
                        z <- thisz
                        z$result <- matrix(ncol=n.phe, nrow=n.pos)
                    }
                    z$result[,thiscol] <- matrix(thisz$result, nrow=n.pos)
                }
            }
            else {
                z <- .C("R_scanone_imp",
                        as.integer(n.ind),
                        as.integer(n.pos),
                        as.integer(n.gen),
                        as.integer(n.draws),
                        as.integer(draws),
                        as.double(ac),
                        as.integer(n.ac),
                        as.double(ic),
                        as.integer(n.ic),
                        as.double(pheno),
                        as.integer(n.phe), # number of phenotypes
                        as.double(weights),
                        result=as.double(rep(0,n.phe*n.pos)),
                        as.integer(ind.noqtl),       # indicators of ind'l w/o QTL effects
                        PACKAGE="qtl")
            }
        }
        else if(method=="hk") { # Haley-Knott regression
            if(n.phe > batchsize) {
                firstcol <- 1
                z <- NULL
                while(firstcol <= n.phe) {
                    if(verbose > 2) cat("chr", names(cross$geno)[i], "phe", firstcol, "...\n")
                    thiscol <- firstcol + 0:(batchsize-1)
                    thiscol <- thiscol[thiscol <= n.phe]
                    thisz <- .C("R_scanone_hk",
                                as.integer(n.ind),         # number of individuals
                                as.integer(n.pos),         # number of markers
                                as.integer(n.gen),         # number of possible genotypes
                                as.double(genoprob),       # genotype probabilities
                                as.double(ac),         # additive covariates
                                as.integer(n.ac),
                                as.double(ic),         # interactive covariates
                                as.integer(n.ic),
                                as.double(pheno[,thiscol]),          # phenotype data
                                as.integer(length(thiscol)), # number of phenotypes
                                as.double(weights),
                                result=as.double(rep(0,length(thiscol)*n.pos)),
                                as.integer(ind.noqtl),       # indicators of ind'l w/o QTL effects
                                PACKAGE="qtl")
                    firstcol <- firstcol + batchsize
                    if(is.null(z)) {
                        z <- thisz
                        z$result <- matrix(ncol=n.phe, nrow=n.pos)
                    }
                    z$result[,thiscol] <- matrix(thisz$result, nrow=n.pos)
                }
            }
            else {
                z <- .C("R_scanone_hk",
                        as.integer(n.ind),         # number of individuals
                        as.integer(n.pos),         # number of markers
                        as.integer(n.gen),         # number of possible genotypes
                        as.double(genoprob),       # genotype probabilities
                        as.double(ac),         # additive covariates
                        as.integer(n.ac),
                        as.double(ic),         # interactive covariates
                        as.integer(n.ic),
                        as.double(pheno),          # phenotype data
                        as.integer(n.phe), # number of phenotypes
                        as.double(weights),
                        result=as.double(rep(0,n.phe*n.pos)),
                        as.integer(ind.noqtl),       # indicators of ind'l w/o QTL effects
                        PACKAGE="qtl")
            }
        }
        else if(method=="ehk")  { # extended Haley-Knott method
            z <- .C("R_scanone_ehk",
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.gen),         # number of possible genotypes
                    as.double(genoprob),       # genotype probabilities
                    as.double(ac),         # additive covariates
                    as.integer(n.ac),
                    as.double(ic),         # interactive covariates
                    as.integer(n.ic),
                    as.double(pheno),          # phenotype data
                    as.double(weights),
                    result=as.double(rep(0,n.pos)),
                    as.integer(maxit),
                    as.double(tol),
                    PACKAGE="qtl")
        }

        else if(method=="em" && model=="normal")  # interval mapping
            z <- .C("R_scanone_em",
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.gen),         # number of possible genotypes
                    as.double(genoprob),       # genotype probabilities
                    as.double(ac),
                    as.integer(n.ac),
                    as.double(ic),
                    as.integer(n.ic),
                    as.double(pheno),          # phenotype data
                    as.double(weights),
                    result=as.double(rep(0,n.pos)),
                    as.integer(std.start),
                    as.double(this.start),
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(0), # debugging verbose off
                    as.integer(ind.noqtl),       # indicators of ind'l w/o QTL effects
                    PACKAGE="qtl")

        else if(model=="np")  # non-parametric interval mapping
            z <- .C("R_scanone_np",
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.gen),         # number of possible genotypes
                    as.double(genoprob),       # genotype probabilities
                    as.double(pheno) ,         # phenotype data
                    result=as.double(rep(0,n.pos)),
                    PACKAGE="qtl")

        else
            stop("Model ", model, " with method ", method, " not available")

        z <- matrix(z$result,nrow=n.pos)

        # interval mapping without covariates:
        #   rescale log likelihood
        if(model == "np" && !ties.random)
            z <- z/correct  # correct for ties

        # setup column names for z
        if(length(pheno.col)==1)
            colnames(z) <- "lod"
        else {
            if(is.null(colnames(pheno)))
                colnames(z) <- paste("lod", pheno.col, sep="")
            else
                colnames(z) <- colnames(pheno)
        }

        # get null log10 likelihood
        if(i==1 & model != "np") {
            if(n.ac > 0)
                resid0 <- lm(pheno ~ ac, weights=weights^2)$resid
            else
                resid0 <- lm(pheno ~ 1, weights=weights^2)$resid

            if(method=="hk")  {
                if(n.phe > 1) {
                    nllik0 <- apply(resid0, 2, function(x)
                                    (n.ind/2)*log10(sum((x*weights)^2)))
                }
                else
                    nllik0 <- (n.ind/2)*log10(sum((resid0*weights)^2))
            }

            else {
                sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
                nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
            }
        }

        # re-scale with null log10 likel for methods em and ehk
        if(method=="em" && model=="normal")
            z <- nllik0 - z

        else if(method=="ehk")
            z <- z/log(10) + nllik0

        else if(method == "hk") {
            if(n.phe > 1) {
                z <- t(nllik0 - t(z))
            }
            else
                z <- nllik0 - z
        }

        # get null log10 likelihood for the X chromosome
        if(chrtype=="X") {

            # determine which covariates belong in null hypothesis
            temp <- scanoneXnull(type, sexpgm, cross.attr=attributes(cross))
            adjustX <- temp$adjustX
            parX0 <- temp$parX0+n.ac
            sexpgmcovar <- temp$sexpgmcovar

            if(adjustX) {
                if(model == "np") {
                    sexpgmcovar <- factor(apply(sexpgmcovar,1,paste,collapse=":"))
                    nllikX <- kruskal.test(pheno ~ sexpgmcovar)$stat/(2*log(10))
                    z <- z - nllikX
                }
                else if(method=="mr") {
                    for(s in 1:ncol(newgeno)) {
                        wh <- newgeno[,s] != 0

                        if(n.ac > 0) {
                            residX <- lm(pheno ~ ac+sexpgmcovar, weights=weights^2,subset=wh)$resid
                            resid0 <- lm(pheno ~ ac, weights=weights^2,subset=wh)$resid
                        }
                        else {
                            residX <- lm(pheno ~ sexpgmcovar, weights=weights^2,subset=wh)$resid
                            resid0 <- lm(pheno ~ 1, weights=weights^2,subset=wh)$resid
                        }
                        nllikX <- (sum(wh)/2)*log10(sum((residX*weights[wh])^2))
                        nllik0 <- (sum(wh)/2)*log10(sum((resid0*weights[wh])^2))

                        # rescale LOD score
                        z[s,] <- z[s,] + nllikX - nllik0
                    }
                }
                else {
                    if(n.ac > 0) {
                        outX <- lm(pheno ~ ac+sexpgmcovar, weights=weights^2)
                        residX <- outX$resid
                        # revise the parX0, if some columns got dropped
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
                        if(method=="imp") {
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
                                for(j in seq(along=nllik0))
                                    nllik0[j] <- -sum(dnorm(resid0[,j],0,sig0[j]/weights,log=TRUE))/log(10)
                            }
                            else {
                                sig0 <- sqrt(sum((resid0*weights)^2)/n.ind)
                                nllik0 <- -sum(dnorm(resid0,0,sig0/weights,log=TRUE))/log(10)
                            }
                        }

                        if(n.phe > 1) {
                            sigX <- sqrt(apply(residX, 2, function(a,b) sum((a*b)^2),weights)/n.ind)
                            nllikX <- sigX
                            for(j in seq(along=nllikX))
                                nllikX[j] <- -sum(dnorm(residX[,j],0,sigX[j]/weights,log=TRUE))/log(10)
                        }
                        else {
                            sigX <- sqrt(sum((residX*weights)^2)/n.ind)
                            nllikX <- -sum(dnorm(residX,0,sigX/weights,log=TRUE))/log(10)
                        }
                    }
                }

                if(method != "mr" && model != "np") {
                    # rescale LOD score
                    z <- t(t(z) + nllikX - nllik0)
                }
            }
        }

        # replace missing or negative LODs with 0
        z[is.na(z) | z<0] <- 0

        w <- marnam
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")

        z <- as.data.frame(z, stringsAsFactors=TRUE)
        z <- cbind(chr=rep(names(cross$geno)[i],length(map)),
                   pos=as.numeric(map), z)
        rownames(z) <- w

        results <- rbind(results, z)
    }

    class(results) <- c("scanone","data.frame")
    attr(results,"method") <- method
    attr(results,"type") <- type
    attr(results,"model") <- model

    results
}




######################################################################
#
# scanone.perm: Permutation test of scanone
#
######################################################################

scanone.perm <-
    function(cross, pheno.col=1, model=c("normal","binary","2part","np"),
             method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use=c("all.obs", "complete.obs"), upper=FALSE,
             ties.random=FALSE, start=NULL, maxit=4000, tol=1e-4,
             n.perm=1000, perm.Xsp=FALSE, perm.strata=NULL,
             verbose=TRUE, batchsize=250)
{
    method <- match.arg(method)
    model <- match.arg(model)
    use <- match.arg(use)

    if((model=="2part" || model=="np") && (!is.null(addcovar) || !is.null(intcovar))) {
        warning("Use of covariates not available for model ", model)
        addcovar <- intcovar <- NULL
    }

    chr.type <- sapply(cross$geno, class)
    if((all(chr.type=="X") || all(chr.type=="X")) && perm.Xsp==TRUE)
        warning("All chromosomes of the same type, so X-chr specific permutations not needed.\n")

    if(any(chr.type=="X") && any(chr.type=="A") && perm.Xsp) { # autosome and X-specific perms
        # X chr versus autosomes
        xchr <- chr.type=="X"
        names(xchr) <- names(cross$geno)

        # chromosome lengths
        L <- summary(pull.map(cross))[,2]
        L <- L[-length(L)]
        La <- sum(L[!xchr])
        Lx <- sum(L[xchr])

        n.perm.X <- ceiling(La/Lx*n.perm)

        if(verbose) cat("--Autosome permutations\n")
        resA <- scanone.perm.engine(n.perm, subset(cross, chr=!xchr),
                                    pheno.col, model,
                                    method, addcovar, intcovar, weights, use,
                                    upper, ties.random, start, maxit, tol,
                                    verbose, perm.strata, batchsize)

        if(verbose) cat("--X chromosome permutations\n")
        resX <- scanone.perm.engine(n.perm.X, subset(cross, chr=xchr),
                                    pheno.col, model,
                                    method, addcovar, intcovar, weights, use,
                                    upper, ties.random, start, maxit, tol,
                                    verbose, perm.strata, batchsize)
        res <- list("A"=resA, "X"=resX)
        attr(res, "xchr") <- xchr
        attr(res, "L") <- c("A"=La, "X"=Lx)
    }

    else {
        res <- scanone.perm.engine(n.perm, cross, pheno.col, model,
                                   method, addcovar, intcovar, weights, use,
                                   upper, ties.random, start, maxit, tol,
                                   verbose, perm.strata, batchsize)
    }

    attr(res,"method") <- method
    attr(res,"model") <- model
    attr(res,"type") <- class(cross)[1]

    if(any(chr.type=="X") && any(chr.type=="A") && perm.Xsp)
        class(res) <- c("scanoneperm", "list")
    else
        class(res) <- c("scanoneperm", "matrix")

    res
}

##################################################################
# engine function to permform permutation test
##################################################################
scanone.perm.engine <-
    function(n.perm, cross, pheno.col, model,
             method, addcovar, intcovar, weights, use,
             upper, ties.random, start, maxit, tol,
             verbose, perm.strata, batchsize=250)
{
    ## local variables
    n.phe <- length(pheno.col)
    n.addcov <- ncol(addcovar)
    n.intcovar <- ncol(intcovar)
    n.ind <- dim(cross$pheno)[1]

    if(method=="mr-imp") # save version with missing genotypes
        tempcross <- cross
    if(method=="mr-argmax") # impute genotypes
        cross <- fill.geno(cross,method="argmax")

    ## if there's only one phenotype, no covariate, and method is imp or hk,
    ## generate permuted phenotype as a matrix
    ## we also need one sex and one direction, or that the
    ##     stratification is within those groups
    batch.mode <- FALSE
    if( (n.phe==1) && ((method=="imp") || (method=="hk")) &&
       model == "normal" &&
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
        ## if there's only one phenotype, no covariate, and method is imp or hk,
        ## generate permuted phenotype as a matrix and do permutation
        ## as multiple phenotypes
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
        tem <- scanone(cross,,pheno.col,model,method,addcovar,
                       intcovar, weights, use, upper,ties.random,start,
                       maxit,tol,n.perm= -1, perm.Xsp=FALSE, perm.strata, verbose=FALSE, batchsize,
                       n.cluster=0)

        res <- matrix(apply(tem[,-(1:2),drop=FALSE], 2, max, na.rm=TRUE), ncol=1)
        colnames(res) <- "lod"
    }
    else { ## all other cases, do one permutation at a time
        ## rnd: how often to print tracing information
        if(verbose > 1) rnd <- 1
        else {
            if(n.perm >= 1000) rnd <- 20
            else if(n.perm >= 100) rnd <- 5
            else rnd <- 1
        }

        addcovarp <- addcovar
        intcovarp <- intcovar
        if(!is.null(addcovar)) addcovarp <- as.matrix(addcovarp)
        if(!is.null(intcovar)) intcovarp <- as.matrix(intcovarp)

        if(model=="2part") res <- matrix(ncol=3*n.phe,nrow=n.perm)
        else res <- matrix(0, n.perm, n.phe)

        for(i in 1:n.perm) {
            if(verbose && i/rnd == floor(i/rnd))
                cat("Permutation", i, "\n")
            # impute genotypes for method "mr-imp"
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
            if(!is.null(weights)) weights <- weights[o]

            tem <- scanone(cross,,pheno.col,model,method,addcovarp,
                           intcovarp,weights,use,upper,ties.random,start,
                           maxit,tol,n.perm= -i, perm.Xsp=FALSE, perm.strata, verbose=FALSE, batchsize,
                           n.cluster=0)

            res[i,] <- apply(tem[,-(1:2),drop=FALSE], 2, max, na.rm=TRUE)

        } # finish permutation

        colnames(res) <- colnames(tem)[-(1:2)]
    }

    ## set row and column names when n.phe>1
    rownames(res) <- 1:n.perm

    ## return
    res

}

# end of scanone.R
