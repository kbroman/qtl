######################################################################
#
# vbscan.R
#
# copyright (c) 2001-2015, Karl W Broman
# last modified Oct, 2015
# first written May, 2001
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
# Contains: vbscan
#
######################################################################

######################################################################
#
# vbscan: scan genome for a quantitative phenotype for which some
# individuals' phenotype is undefined (for example, the size of a
# lesion, where some individuals have no lesion).
#
######################################################################

vbscan <-
    function(cross, pheno.col=1, upper=FALSE, method="em", maxit=4000,
             tol=1e-4)
{
    method <- match.arg(method)
    type <- class(cross)[1]

    # check arguments are okay
    if(length(pheno.col) > 1) pheno.col <- pheno.col[1]
    if(pheno.col > nphe(cross))
        stop("Specified phenotype column exceeds the number of phenotypes")
    y <- cross$pheno[,pheno.col]

    # modify phenotypes
    if(upper) {
        if(!any(y == Inf)) y[y==max(y)] <- Inf
    }
    else {
        if(!any(y == -Inf)) y[y==min(y)] <- -Inf
    }
    survived <- rep(0,length(y))
    survived[y == -Inf | y == Inf] <- 1

    # The following line is included since .C() doesn't accept Infs
    y[y == -Inf | y == Inf] <- 99999

    n.chr <- nchr(cross)
    results <- NULL

    for(i in 1:n.chr) {
        # make sure inferred genotypes or genotype probabilities are available
        if(!("prob" %in% names(cross$geno[[i]]))) {
            cat(" -Calculating genotype probabilities\n")
            cross <- calc.genoprob(cross)
        }

        genoprob <- cross$geno[[i]]$prob
        n.pos <- dim(genoprob)[2]
        n.ind <- length(y)

        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") sexpgm <- getsex(cross)
        else sexpgm <- NULL

        gen.names <- getgenonames(type,chrtype,"full", sexpgm, attributes(cross))
        n.gen <- length(gen.names)

        # revise X chromosome genotypes
        if(chrtype=="X" && (type=="f2" || type=="bc"))
            genoprob <- reviseXdata(type, "full", sexpgm, prob=genoprob,
                                    cross.attr=attributes(cross))

        z <- .C("R_vbscan",
                as.integer(n.pos),
                as.integer(n.ind),
                as.integer(n.gen),
                as.double(genoprob),
                as.double(y),
                as.integer(survived),
                lod=as.double(rep(0, n.pos*3)),
                as.integer(maxit),
                as.double(tol),
                PACKAGE="qtl")

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

        res <- data.frame(chr=rep(names(cross$geno)[i],length(map)),
                          pos = as.numeric(map),
                          matrix(z$lod,nrow=n.pos,byrow=TRUE), stringsAsFactors=TRUE)

        colnames(res)[-(1:2)] <- c("lod.p.mu","lod.p","lod.mu")

        w <- names(map)
        o <- grep("^loc-*[0-9]+",w)

        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
        rownames(res) <- w

        z <- res

        # get null log10 likelihood for the X chromosome
        if(chrtype=="X") {

            # determine which covariates belong in null hypothesis
            temp <- scanoneXnull(type, sexpgm, cross.attr=attributes(cross))
            adjustX <- temp$adjustX
            parX0 <- temp$parX0
            sexpgmcovar <- temp$sexpgmcovar
            sexpgmcovar.alt <- temp$sexpgmcovar.alt

            if(adjustX) { # get LOD-score adjustment
                n.gen <- ncol(sexpgmcovar)+1
                genoprob <- matrix(0,nrow=n.ind,ncol=n.gen)
                for(i in 1:n.gen)
                    genoprob[sexpgmcovar.alt==i,i] <- 1

                nullz <- .C("R_vbscan",
                            as.integer(1),
                            as.integer(n.ind),
                            as.integer(n.gen),
                            as.double(genoprob),
                            as.double(y),
                            as.integer(survived),
                            lod=as.double(rep(0,(4+2*n.gen))),
                            as.integer(maxit),
                            as.double(tol),
                            PACKAGE="qtl")

                # adjust LOD curve
                for(i in 1:3) z[,i+2] <- z[,i+2] - nullz$lod[i]
            }
        }

        results <- rbind(results, z)
    }

    class(results) <- c("scanone","data.frame")
    attr(results,"method") <- method
    attr(results,"type") <- class(cross)[1]
    attr(results,"model") <- "twopart"

    results
}

# end of vbscan.R
