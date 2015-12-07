######################################################################
#
# discan.R
#
# copyright (c) 2001-2015, Karl W Broman
# last modified Oct, 2015
# first written Oct, 2001
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
# Contains: discan
#
######################################################################

######################################################################
#
# discan: scan genome, calculating LOD scores with single QTL model
#         for a dichotomous trait
#
######################################################################

discan <-
    function(cross, pheno, method=c("em","hk","mr"),
             addcovar=NULL, intcovar=NULL, maxit=4000, tol=1e-4,
             verbose=FALSE, give.warnings=TRUE, ind.noqtl)
{
    method <- match.arg(method)

    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    type <- class(cross)[1]
    if(is.null(addcovar)) n.addcovar <- 0
    else n.addcovar <- ncol(addcovar)
    if(is.null(intcovar)) n.intcovar <- 0
    else n.intcovar <- ncol(intcovar)

    # individuals with no QTL effect
    if(missing(ind.noqtl)) ind.noqtl <- rep(FALSE, nind(cross))
    else {
        if(!is.logical(ind.noqtl) || length(ind.noqtl) != nind(cross))
            stop("ind.noqtl be a logical vector of length n.ind (", nind(cross), ")")

        if(sum(ind.noqtl) > 1) {
            if(method == "mr") {
                ind.noqtl <- rep(FALSE, nind(cross))
                warning("ind.noqtl ignored for method=", method, ", model=binary")
            }
            else if(is.null(addcovar) && (!is.logical(ind.noqtl) || any(ind.noqtl))) {
                ind.noqtl <- rep(FALSE, nind(cross))
                warning("ind.noqtl ignored when no additive covariates")
            }
        }
    }

    if(method=="mr" && n.addcovar+n.intcovar>0)  {
        if(give.warnings) warning("Covariates ignored with method=\"mr\"; use \"em\" instead")
        n.addcovar <- n.intcovar <- addcovar <- intcovar <- 0
    }

    u <- unique(pheno)
    if(any(u != 0 & u != 1))
        stop("Phenotypes must be either 0 or 1.")

    results <- NULL

    llik0 <- c("A"=NA,"X"=NA)
    nullcoef <- list("A"=NA,"X"=NA)

    for(i in 1:n.chr) {

        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") {
            sexpgm <- getsex(cross)

            ac <- revisecovar(sexpgm,addcovar)
            if(!is.null(addcovar) && (nd <- attr(ac, "n.dropped")) > 0 && give.warnings)
                warning("Dropped ", nd, " additive covariates on X chromosome.")
            if(length(ac)==0) {
                n.ac <- 0
                ac <- NULL
            }
            else n.ac <- ncol(ac)
            ic <- revisecovar(sexpgm,intcovar)
            if(!is.null(intcovar) && (nd <- attr(ic, "n.dropped")) > 0 && give.warnings)
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

        # get null log liklihood
        if(is.na(llik0[chrtype])) {
            if(n.ac > 0)
                nullfit <- glm(pheno ~ ac, family=binomial(link="logit"))
            else if(i==1)
                nullfit <- glm(pheno ~ 1, family=binomial(link="logit"))
            fitted <- nullfit$fitted
            nullcoef[[chrtype]] <- nullfit$coef
            llik0[chrtype] <- sum(pheno*log10(fitted) + (1-pheno)*log10(1-fitted))
        }

        # get genotype names
        gen.names <- getgenonames(type,chrtype,"full",sexpgm,attributes(cross))
        n.gen <- length(gen.names)

        # pull out genotype data (mr)
        # or genotype probabilities (em)
        if(method == "mr") {
            newgeno <- cross$geno[[i]]$data
            newgeno[is.na(newgeno)] <- 0

            # discard partially informative genotypes
            if(type=="f2") newgeno[newgeno>3] <- 0
            if(type=="4way") newgeno[newgeno>4] <- 0

            # revise X chromosome genotypes
            if(chrtype=="X" && (type=="bc" || type=="f2"))
                newgeno <- reviseXdata(type, "full", sexpgm, geno=newgeno,
                                       cross.attr=attributes(cross))

            n.pos <- ncol(newgeno)
            map <- cross$geno[[i]]$map
            if(is.matrix(map)) map <- map[1,]

            z <- .C("R_discan_mr",
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.gen),         # number of possible genotypes
                    as.integer(newgeno),       # genotype data
                    as.integer(pheno),          # phenotype data
                    result=as.double(rep(0,n.pos*(n.gen+1))),
                    PACKAGE="qtl")

        }
        else {
            if(!("prob" %in% names(cross$geno[[i]]))) { # need to run calc.genoprob
                if(give.warnings) warning("First running calc.genoprob.")
                cross <- calc.genoprob(cross)
            }
            genoprob <- cross$geno[[i]]$prob
            n.pos <- ncol(genoprob)

            # revise X chromosome genotypes
            if(chrtype=="X" && (type=="bc" || type=="f2"))
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

            if(is.matrix(map)) map <- map[1,]

            if(method=="hk") {
                z <- .C("R_scanone_hk_binary",
                        as.integer(n.ind),         # number of individuals
                        as.integer(n.pos),         # number of markers
                        as.integer(n.gen),         # number of possible genotypes
                        as.double(genoprob),       # genotype probabilities
                        as.double(ac),
                        as.integer(n.ac),
                        as.double(ic),
                        as.integer(n.ic),
                        as.double(pheno),          # phenotype data
                        result=as.double(rep(0,n.pos)),
                        as.double(tol),
                        as.integer(maxit),
                        as.integer(verbose),
                        as.integer(ind.noqtl),
                        PACKAGE="qtl")
            }
            else if(n.ac + n.ic > 0) {

                start <- rep(nullcoef[[chrtype]][1],n.gen)
                if(n.ac > 0)
                    start <- c(start, nullcoef[[chrtype]][-1])
                if(n.ic > 0)
                    start <- c(start, rep(0, n.ic*(n.gen-1)))

                z <- .C("R_discan_covar",
                        as.integer(n.ind),         # number of individuals
                        as.integer(n.pos),         # number of markers
                        as.integer(n.gen),         # number of possible genotypes
                        as.double(genoprob),       # genotype probabilities
                        as.double(ac),
                        as.integer(n.ac),
                        as.double(ic),
                        as.integer(n.ic),
                        as.integer(pheno),          # phenotype data
                        as.double(start),
                        result=as.double(rep(0,n.pos)),
                        as.integer(maxit),
                        as.double(tol),
                        as.integer(verbose),
                        as.integer(ind.noqtl),
                        PACKAGE="qtl")
            }
            else {
                z <- .C("R_discan_im",
                        as.integer(n.ind),         # number of individuals
                        as.integer(n.pos),         # number of markers
                        as.integer(n.gen),         # number of possible genotypes
                        as.double(genoprob),       # genotype probabilities
                        as.integer(pheno),          # phenotype data
                        result=as.double(rep(0,n.pos)),
                        as.integer(maxit),
                        as.double(tol),
                        PACKAGE="qtl")
            }

        }
        z <- matrix(z$result,nrow=n.pos)

        if(method != "mr") z[,1] <- z[,1] - llik0[chrtype]
        z[is.na(z[,1]),1] <- 0

        z <- z[,1,drop=FALSE]

        colnames(z)[1] <- "lod"

        # get null log10 likelihood for the X chromosome
        adjustX <- FALSE
        if(chrtype=="X") {
            # determine which covariates belong in null hypothesis
            temp <- scanoneXnull(type, sexpgm, cross.attr=attributes(cross))
            adjustX <- temp$adjustX
            parX0 <- temp$parX0+n.ac
            sexpgmcovar <- temp$sexpgmcovar
            sexpgmcovar.alt <- temp$sexpgmcovar.alt

            if(adjustX) { # get LOD-score adjustment
                if(n.ac > 0) {
                    nullfitX <- glm(pheno ~ ac+sexpgmcovar,
                                    family=binomial(link="logit"))
                    parX0 <- lm(pheno~ac+sexpgmcovar)$rank
                }
                else
                    nullfitX <- glm(pheno ~ sexpgmcovar,
                                    family=binomial(link="logit"))
                fittedX <- nullfitX$fitted
                llik0X <- sum(pheno*log10(fittedX) + (1-pheno)*log10(1-fittedX))

                # adjust LOD curve
                z <- z - (llik0X - llik0["X"])
            }

        }

        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")

        z <- data.frame(chr=rep(names(cross$geno)[i],length(map)),
                        pos=as.numeric(map), z, stringsAsFactors=TRUE)
        rownames(z) <- w

        results <- rbind(results, z)
    } # loop over chromosomes

    class(results) <- c("scanone","data.frame")
    attr(results,"method") <- method
    attr(results,"type") <- type
    attr(results,"model") <- "binary"
    attr(results,"null.log10.lik") <- llik0["A"]
    if(adjustX)
        attr(results,"null.log10.lik.X") <- llik0X

    results
}

# end of discan.R
