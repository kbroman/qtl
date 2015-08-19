## scantwopermhk.R

scantwopermhk <-
    function(cross, chr, pheno.col=1, addcovar=NULL, weights=NULL, n.perm=1,
             batchsize=1000, perm.strata=NULL, perm.Xsp=NULL, verbose=FALSE,
             assumeCondIndep=FALSE)
{
    if(!missing(chr) && !is.null(chr))
        cross <- subset(cross, chr)

    # in RIL, treat X chromosome like an autosome
    chr.names <- chrnames(cross)
    chrtype <- sapply(cross$geno, class)
    type <- class(cross)[1]
    if(any(chrtype=="X") && (type=="risib" || type=="riself"))
        for(i in which(chrtype=="X"))
            class(cross$geno[[i]]) <- chrtype[i] <- "A"

    if(any(chrtype=="X") && (type=="bc" || type=="f2")) # force stratified permutation test
        perm.strata <- force_sexstrata(cross, perm.strata)

    if(!assumeCondIndep) { # if reduce2grid was used, for assumeCondIndep
        # if reduced2grid, force assumeCondIndep=TRUE
        reduced2grid <- attr(cross$geno[[1]]$prob, "reduced2grid")
        if(!is.null(reduced2grid) && reduced2grid) {
            assumeCondIndep <- TRUE
            warning("Using assumeCondIndep=TRUE, since probabilities reduced to grid")
        }
    }

    if(is.null(perm.Xsp) || !perm.Xsp || !any(chrtype=="X")) { # all autosomes
        result <- .scantwopermhk(cross, pheno.col=pheno.col,
                                 addcovar=addcovar, weights=weights, n.perm=n.perm,
                                 batchsize=batchsize,
                                 perm.strata=perm.strata, verbose=verbose,
                                 assumeCondIndep=assumeCondIndep)
    }
    else { # separate A:A, A:X, and X:X
        # covariates for X chr
        Xcovar <- scanoneXnull(class(cross)[1], getsex(cross), attributes(cross))$sexpgmcovar

        # lengths of autosomes and X chr
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


        if(verbose) message("Running ", n.permAA, " A:A permutations")
        AAresult <- .scantwopermhk(cross, chr=Achr, pheno.col=pheno.col,
                                   addcovar=cbind(addcovar, Xcovar), # include X covariates
                                   weights=weights, n.perm=n.permAA,
                                   batchsize=batchsize,
                                   perm.strata=perm.strata, verbose=verbose,
                                   assumeCondIndep=assumeCondIndep)

        if(verbose) message("Running ", n.permXX, " X:X permutations")
        XXresult <- .scantwopermhk(cross, chr=Xchr, pheno.col=pheno.col,
                                   addcovar=addcovar, weights=weights, n.perm=n.permXX,
                                   batchsize=batchsize,
                                   perm.strata=perm.strata, verbose=verbose,
                                   assumeCondIndep=assumeCondIndep)

        if(verbose) message("Running ", n.permAX, " A:X permutations")
        AXresult <- .scantwopermhk(cross, chr=list(Achr, Xchr),
                                   pheno.col=pheno.col,
                                   addcovar=addcovar, weights=weights, n.perm=n.permAX,
                                   batchsize=batchsize,
                                   perm.strata=perm.strata, verbose=verbose,
                                   assumeCondIndep=assumeCondIndep)

        result <- list(AA=AAresult, AX=AXresult, XX=XXresult)
        attr(result, "L") <- c(A=AL, X=XL)
        attr(result, "LL") <- c(AA=AAL, AX=AXL, XX=XXL)
        names(chrtype) <- chr.names
        attr(result, "chrtype") <- chrtype
    }

    class(result) <- "scantwoperm"
    result
}



.scantwopermhk <-
    function(cross, chr, pheno.col=1, addcovar=NULL,
             weights=NULL, n.perm=1, batchsize=1000,
             perm.strata=NULL, verbose=FALSE, assumeCondIndep=FALSE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(n.perm > batchsize) { # run in batches
        if(missing(chr)) chr <- chrnames(cross)

        # determine size of batches: as equal as possible
        n.batches <- ceiling(n.perm/batchsize)
        batchsizes <- rep(ceiling(n.perm/n.batches), n.batches)
        batchsizes[n.batches] <- n.perm - sum(batchsizes[-n.batches])

        result <- NULL
        if(verbose) message(" - perms in ", n.batches, " batches")
        for(i in seq(along=batchsizes)) {
            if(verbose) message("   -- batch ", i)
            thisresult <- .scantwopermhk(cross, chr, pheno.col=pheno.col, addcovar=addcovar,
                                         weights=weights, n.perm=batchsizes[i], batchsize=Inf,
                                         perm.strata=perm.strata, verbose=verbose, assumeCondIndep)
            if(is.null(result)) result <- thisresult
            else {
                for(j in seq(along=result))
                    result[[j]] <- rbind(result[[j]], thisresult[[j]])
            }
        }
        return(result)
    }


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

    # subset cross and check chr arguments
    cross <- subset(cross, unique(c(chr1, chr2)))
    thechr <- names(cross$geno)
    nchr1 <- match(chr1, thechr)
    nchr2 <- match(chr2, thechr)
    if(!any(sapply(nchr1, function(a,b) any(a <= b), nchr2)))
        stop("Need some of first chr to be <= some of second chr")

    # in RIL, treat X chromomse like an autosome
    chrtype <- sapply(cross$geno, class)
    type <- class(cross)[1]
    if(any(chrtype=="X") && (type == "risib" || type == "riself"))
        for(i in which(chrtype=="X")) class(cross$geno[[i]]) <- "A"

    # check perm.strat
    if(!missing(perm.strata) && !is.null(perm.strata)) {
        if(length(perm.strata) != nind(cross))
            stop("perm.strata, if given, must have length = nind(cross) [", nind(cross), "]")
        if(any(is.na(perm.strata)))
            stop("perm.strata cannot have missing values")
    }

    # grab phenotypes
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
    if(length(pheno.col) > 1)
        stop("Only one phenotype column allowed")
    if(pheno.col < 1 || pheno.col > nphe(cross))
        stop("pheno.col values should be between 1 and the no. phenotypes")

    # if stepwidth="variable" or stepwidth=="max" when calling calc.genoprob or sim.geno,
    # we force incl.markers=TRUE; I assume it is the same for all chromosomes
    stepwidth.var <- FALSE
    if("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
       attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
        stepwidth.var <- TRUE
        incl.markers <- TRUE
    }

    # omit individuals with missing phenotypes or covariates
    temp <- checkcovar(cross, pheno.col, addcovar, intcovar=NULL,
                       perm.strata, ind.noqtl=NULL, weights, TRUE)
    cross <- temp[[1]]
    pheno <- temp[[2]]
    addcovar <- temp[[3]]
    n.addcovar <- temp[[5]]
    perm.strata <- temp[[7]]
    weights <- temp[[9]]
    n.ind <- length(pheno)
    n.chr <- nchr(cross)
    if(is.null(weights)) weights <- rep(1, n.ind)

    # null log likelihood
    if(n.addcovar > 0)
        resid0 <- lm(pheno ~ addcovar, weights=weights^2)$resid
    else
        resid0 <- lm(pheno ~ 1, weights=weights^2)$resid
    nllik0X <- nllik0 <- (n.ind/2)*log10(sum((resid0*weights)^2))

    # reorganize perm.strata
    if(!missing(perm.strata) && !is.null(perm.strata)) {
        u <- unique(perm.strata)
        perm.strata <- match(perm.strata, u) # turn into integers in {1, ..., n.strata}
        n.strata <- length(u)
    }
    else {
        perm.strata <- rep(0, n.ind)
        n.strata <- 0
    }

    # X chromosome covariates
    if(any(chrtype=="X")) {
        sexpgm <- getsex(cross)

        # get null loglik for X chr
        Xnullcovar <- scanoneXnull(class(cross)[1], sexpgm, attributes(cross))[[3]]
        adjcovar <- cbind(Xnullcovar, addcovar)
        if(!is.null(adjcovar) && ncol(adjcovar) > 0) {
            resid0 <- lm(pheno ~ adjcovar, weights=weights^2)$resid
            nllik0X <- (n.ind/2)*log10(sum((resid0*weights)^2))
        }

        # covariates to use on X chr
        addcovarX <- revisecovar(sexpgm,addcovar)
        if(!is.null(addcovar) && (nd <- attr(addcovarX, "n.dropped")) > 0 && n.perm > -2)
            warning("Dropped ", nd, " additive covariates on X chromosome.")
        if(length(addcovarX)==0) {
            n.acX <- 0
            addcovarX <- NULL
        }
        else n.acX <- ncol(addcovarX)


        for(i in 1:n.chr) {
            if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft"))) {
                oldXchr <- subset(cross, chr=thechr[i])
                cross$geno[[i]]$prob <-
                    reviseXdata(type, "full", sexpgm, prob=cross$geno[[i]]$prob,
                                cross.attr=attributes(cross))
            }
        }
    }

    # no. genotypes and positions on each chromosome
    n.gen <- n.pos <- rep(0, n.chr)
    for(i in 1:n.chr) {
        n.pos[i] <- ncol(cross$geno[[i]]$prob)
        n.gen[i] <- dim(cross$geno[[i]]$prob)[3]
    }

    # create shuffled indices
    if(n.strata==0)
        permindex <- replicate(n.perm, sample(1:n.ind))
    else { # stratified permutation
        permindex <- replicate(n.perm, 1:n.ind)
        for(i in 1:n.strata)
            permindex[perm.strata==i,] <- apply(permindex[perm.strata==i,], 2, sample)
    }
    permindex <- permindex-1

    # begin loop over pairs of chromosomes
    result <- vector("list", length(nchr1)*length(nchr2))
    k <- 0
    for(i in nchr1) { # loop over the 1st chromosome
        for(j in nchr2) { # loop over the 2nd chromosome
            if(j < i) next
            k <- k + 1

            if(chrtype[i]=="X" || chrtype[j]=="X") {
                ac <- addcovarX
                n.ac <- n.acX
            }
            else {
                ac <- addcovar
                n.ac <- n.addcovar
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
            if(verbose>1) cat(" (", names(cross$geno)[i], ",",
                            names(cross$geno)[j],")\n",sep="")

            if(i==j) { # same chromosome

                if(verbose>2) cat("  --Calculating joint probs.\n")

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

                # revise pair probilities for X chromosome
                if(chrtype[i]=="X" && (type %in% c("bc","f2","bcsft")))
                    temp <- reviseXdata(type, "full", sexpgm, pairprob=temp,
                                        cross.attr=attributes(cross))

                if(verbose>2) cat("  --Done.\n")

                thisz <- .C("R_scantwopermhk_1chr",
                            as.integer(n.ind),
                            as.integer(n.pos[i]),
                            as.integer(n.gen[i]),
                            as.double(cross$geno[[i]]$prob),
                            as.double(temp),
                            as.double(ac),
                            as.integer(n.ac),
                            as.double(pheno),
                            as.integer(n.perm),
                            as.integer(permindex),
                            as.double(weights),
                            result=as.double(rep(0,n.perm*6)),
                            as.integer(n.col2drop),
                            as.integer(col2drop),
                            PACKAGE="qtl")
                result[[k]] <- -matrix(thisz$result, nrow=n.perm)
                result[[k]][,c(1,4,6)] <- result[[k]][,c(1,4,6)] +
                    ifelse(chrtype[i]=="X", nllik0X, nllik0)
            } # end same chromosome
            else {
                thisz <- .C("R_scantwopermhk_2chr",
                            as.integer(n.ind),
                            as.integer(n.pos[i]),
                            as.integer(n.pos[j]),
                            as.integer(n.gen[i]),
                            as.integer(n.gen[j]),
                            as.double(cross$geno[[i]]$prob),
                            as.double(cross$geno[[j]]$prob),
                            as.double(ac),
                            as.integer(n.ac),
                            as.double(pheno),
                            as.integer(n.perm),
                            as.integer(permindex),
                            as.double(weights),
                            result=as.double(rep(0,n.perm*6)),
                            PACKAGE="qtl")
                result[[k]] <- -matrix(thisz$result, nrow=n.perm)
                result[[k]][,c(1,4,6)] <- result[[k]][,c(1,4,6)] +
                    ifelse(chrtype[i]=="X" || chrtype[j]=="X", nllik0X, nllik0)
            } # end diff chr
        } # end loop chr 2
    } # end loop chr 1

    result <- apply( array(unlist(result[1:k]), dim=c(n.perm, 6, k)), 1:2, max, na.rm=TRUE)
    result <- as.list(as.data.frame(result))
    phename <- phenames(cross)[1]
    result <- lapply(result, function(a) { a <- as.matrix(a); colnames(a) <- phenames(cross)[1]; a })
    names(result) <- c("full", "fv1", "int", "add", "av1", "one")
    class(result) <- c("scantwoperm", "list")
    result
}

# force stratified permutation test for X chr
force_sexstrata <-
function(cross, perm.strata)
{
    type <- class(cross)[1]
    if(type!="bc" && type!="f2")
        return(perm.strata)

    Xcovar <- scanoneXnull(type, getsex(cross), attributes(cross))$sexpgmcovar

    if(!is.null(Xcovar)) {
        pastedcovar <- apply(Xcovar, 1, paste, collapse="")
        upastedcovar <- unique(pastedcovar)
        sexstrata <- match(pastedcovar, upastedcovar)
        if(is.null(perm.strata)) perm.strata <- sexstrata
        else {
            u <- unique(perm.strata)
            perm.strata <- match(perm.strata, u)
            m <- max(perm.strata)
            perm.strata <- perm.strata + m*(sexstrata-1)
        }
    }

    perm.strata
}
