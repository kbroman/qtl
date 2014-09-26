######################################################################
#
# calc.pairprob.R
#
# copyright (c) 2001-2013, Karl W Broman
# last modified Sep, 2013
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
# Part of the R/qtl package
# Contains: calc.pairprob
#
######################################################################

######################################################################
#
# calc.pairprob: calculate joint genotype probabilities for all pairs
#                of putative QTLs, conditional on the observed marker
#                data
#
# This is an *internal* function, not to be called by the user.
#
# The input argument cross is assumed to have just one chromosome.
#
######################################################################

calc.pairprob <-
    function(cross, step=0, off.end=0, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             map, assumeCondIndep=FALSE)
{
    # which type of cross is this?
    type <- class(cross)[1]

    if(assumeCondIndep) { # assume conditional independence of QTL given markers
        if(!("prob" %in% names(cross$geno[[1]]))) {
            cross <- calc.genoprob(subset(cross, chr=1), step=step, off.end=off.end,
                                   error.prob=error.prob, map.function=map.function)
        }
        prob <- cross$geno[[1]]$prob
        n.ind <- dim(prob)[1]
        n.pos <- dim(prob)[2]
        n.gen <- dim(prob)[3]

        if(n.pos < 2) return(NULL)

        z <- .C("R_calc_pairprob_condindep",
                as.integer(n.ind),
                as.integer(n.pos),
                as.integer(n.gen),
                as.double(prob),
                pairprob=as.double(rep(0,n.ind*choose(n.pos, 2)*n.gen*n.gen)),
                PACKAGE="qtl")

        pairprob <- array(z$pairprob, dim=c(n.ind,n.pos*(n.pos-1)/2,n.gen,n.gen))

        return(pairprob)
    }

    if(step==0 && off.end > 0) step <- off.end*2

    # map function
    map.function <- match.arg(map.function)
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }
    n.ind <- nind(cross)
    n.chr <- nchr(cross)

    # type of chromosome?
    chrtype <- class(cross$geno[[1]])
    if(chrtype=="X") xchr <- TRUE
    else xchr <- FALSE

    if(type == "f2") {
        one.map <- TRUE
        if(!xchr) { # autosome
            cfunc <- "calc_pairprob_f2"
            n.gen <- 3
            gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
        }
        else {                             # X chromsome
            cfunc <- "calc_pairprob_bc"
            n.gen <- 2
            gen.names <- c("g1","g2")
        }
    }
    else if(type == "bc") {
        cfunc <- "calc_pairprob_bc"
        n.gen <- 2
        if(!xchr) # autosome
            gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
        else gen.names <- c("g1","g2")
        one.map <- TRUE
    }
    else if(type == "riself" || type=="risib" || type=="dh" || type=="haploid") {
        cfunc <- "calc_pairprob_bc"
        n.gen <- 2
        gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
        one.map <- TRUE
    }
    else if(type == "4way") {
        cfunc <- "calc_pairprob_4way"
        n.gen <- 4
        one.map <- FALSE
        gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
    }
    else if(type == "ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="magic16") {
        cfunc <- paste("calc_pairprob_", type, sep="")
        if(type=="magic16") n.gen <- 16
        else n.gen <- as.numeric(substr(type, 3, 3))
        one.map <- TRUE
        gen.names <- LETTERS[1:n.gen]
        if(xchr)
            warning("calc.pairprob not working properly for the X chromosome for 4- or 8-way RIL.")
    }
    else if(type == "bcsft") {
        one.map <- TRUE
        cfunc <- "calc_pairprob_bcsft"
        cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
        if(!xchr) { # autosome
            gen.names <- getgenonames("bcsft", "A", cross.attr=attributes(cross))
            n.gen <- 2 + (cross.scheme[2] > 0)
        }
        else { ## X chromsome
            cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
            cross.scheme[2] <- 0
            gen.names <- c("g1","g2")
            n.gen <- 2
        }
    }
    else
        stop("calc.pairprob not available for cross type ", type, ".")

    # genotype data
    gen <- cross$geno[[1]]$data
    gen[is.na(gen)] <- 0

    # get recombination fractions
    if(one.map) {
        #    map <- create.map(cross$geno[[1]]$map,step,off.end)
        rf <- mf(diff(map))
        if(type=="risib" || type=="riself")
            rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),class(cross$geno[[1]]))
        rf[rf < 1e-14] <- 1e-14

        # new genotype matrix with pseudomarkers filled in
        newgen <- matrix(ncol=length(map),nrow=nrow(gen))
        colnames(newgen) <- names(map)
        newgen[,colnames(gen)] <- gen
        newgen[is.na(newgen)] <- 0
        n.pos <- ncol(newgen)
        marnames <- names(map)
    }
    else {
        #    map <- create.map(cross$geno[[1]]$map,step,off.end)
        rf <- mf(diff(map[1,]))
        rf[rf < 1e-14] <- 1e-14
        rf2 <- mf(diff(map[2,]))
        rf2[rf2 < 1e-14] <- 1e-14

        # new genotype matrix with pseudomarkers filled in
        newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
        colnames(newgen) <- colnames(map)
        newgen[,colnames(gen)] <- gen
        newgen[is.na(newgen)] <- 0
        n.pos <- ncol(newgen)
        marnames <- colnames(map)
    }

    if(n.pos < 2) return(NULL)

    # below: at least two positions
    # call the C function
    if(one.map) {
        ## Hide cross scheme in genoprob to pass to routine. BY
        temp <- as.double(rep(0,n.gen*n.ind*n.pos))
        if(type == "bcsft")
            temp[1:2] <- cross.scheme

        z <- .C(cfunc,
                as.integer(n.ind),         # number of individuals
                as.integer(n.pos),         # number of markers
                as.integer(newgen),        # genotype data
                as.double(rf),             # recombination fractions
                as.double(error.prob),     #
                as.double(temp),
                pairprob=as.double(rep(0,n.ind*n.pos*(n.pos-1)/2*n.gen^2)),
                PACKAGE="qtl")
    }
    else {
        z <- .C(cfunc,
                as.integer(n.ind),         # number of individuals
                as.integer(n.pos),         # number of markers
                as.integer(newgen),        # genotype data
                as.double(rf),             # recombination fractions
                as.double(rf2),            # recombination fractions
                as.double(error.prob),     #
                as.double(rep(0,n.gen*n.ind*n.pos)),
                pairprob=as.double(rep(0,n.ind*n.pos*(n.pos-1)/2*n.gen^2)),
                PACKAGE="qtl")
    }

    pairprob <- array(z$pairprob, dim=c(n.ind,n.pos*(n.pos-1)/2,n.gen,n.gen))

    # 4- and 8-way RIL: reorganize the results
    if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="bgmagic16")
        pairprob <- reorgRIpairprob(cross, pairprob)

    pairprob
}

# end of calc.pairprob.R
