######################################################################
#
# calc.genoprob.R
#
# copyright (c) 2001-2013, Karl W Broman
# last modified Sep, 2013
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
# Part of the R/qtl package
# Contains: calc.genoprob
#
######################################################################

######################################################################
#
# calc.genoprob: calculate genotype probabilities conditional on
#                observed marker genotypes
#
######################################################################

calc.genoprob <-
    function(cross, step=0, off.end=0, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             stepwidth=c("fixed", "variable", "max"))
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    # map function
    map.function <- match.arg(map.function)
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    stepwidth <- match.arg(stepwidth)

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }

    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    type <- class(cross)[1]

    # calculate genotype probabilities one chromosome at a time
    for(i in 1:n.chr) {
        if(n.mar[i]==1) temp.offend <- max(c(off.end,5))
        else temp.offend <- off.end

        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") xchr <- TRUE
        else xchr <- FALSE

        # which type of cross is this?
        if(type == "f2") {
            one.map <- TRUE
            if(!xchr) { # autosomal
                cfunc <- "calc_genoprob_f2"
                n.gen <- 3
                gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
            }
            else {                             # X chromsome
                cfunc <- "calc_genoprob_bc"
                n.gen <- 2
                gen.names <- c("g1","g2")
            }
        }
        else if(type == "bc") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            if(!xchr)
                gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
            else gen.names <- c("g1","g2")
            one.map <- TRUE
        }
        else if(type == "riself" || type=="risib" || type=="dh" || type=="haploid") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
            one.map <- TRUE
        }
        else if(type == "4way") {
            cfunc <- "calc_genoprob_4way"
            n.gen <- 4
            one.map <- FALSE
            gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
        }
        else if(type=="ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self" || type=="bgmagic16") {
            cfunc <- paste("calc_genoprob_", type, sep="")
            if(type=="bgmagic16") n.gen <- 16
            else n.gen <- as.numeric(substr(type, 3, 3))
            one.map <- TRUE
            gen.names <- LETTERS[1:n.gen]
            if(xchr)
                warning("calc.genoprob not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else if(type == "bcsft") {
            one.map <- TRUE
            cfunc <- "calc_genoprob_bcsft"
            cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
            if(!xchr) { # autosomal
                gen.names <- getgenonames("bcsft", "A", cross.attr=attributes(cross))
                n.gen <- 2 + (cross.scheme[2] > 0)
            }
            else { ## X chr
                cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
                cross.scheme[2] <- 0
                gen.names <- c("g1","g2")
                n.gen <- 2
            }
        }
        else if(type == "ri8selfIRIP1") {
            one.map <- TRUE
            cfunc <- "calc_genoprob_ri8selfIRIP1"
            n.gen <- 8
            gen.names <- LETTERS[1:n.gen]
        }
        else
            stop("calc.genoprob not available for cross type ", type, ".")

        # genotype data
        gen <- cross$geno[[i]]$data
        gen[is.na(gen)] <- 0

        # recombination fractions
        if(one.map) {
            # recombination fractions
            map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
            rf <- mf(diff(map))
            if(type=="risib" || type=="riself")
                rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),chrtype)
            rf[rf < 1e-14] <- 1e-14

            # new genotype matrix with pseudomarkers filled in
            newgen <- matrix(ncol=length(map),nrow=nrow(gen))
            dimnames(newgen) <- list(NULL,names(map))
            newgen[,colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- names(map)
        }
        else {
            map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
            rf <- mf(diff(map[1,]))
            rf[rf < 1e-14] <- 1e-14
            rf2 <- mf(diff(map[2,]))
            rf2[rf2 < 1e-14] <- 1e-14

            # new genotype matrix with pseudomarkers filled in
            newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
            dimnames(newgen) <- list(NULL,dimnames(map)[[2]])
            newgen[,colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- colnames(map)
        }
        ## Cross scheme being added for Ft and BCs.
        ## cross_scheme = c(BC = s, F = t).
        ## BC has cross_scheme = c(1,0)
        ## F2 has cross_scheme = c(0,2)
        ## BCsFt has cross_scheme = c(s,t)
        ## Other designs such as Ft with test cross don't quite fit in (yet).
        ## *** Need to change all the C routines!!

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
                    genoprob=as.double(temp),
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
                    genoprob=as.double(rep(0,n.gen*n.ind*n.pos)),
                    PACKAGE="qtl")
        }

        # re-arrange marginal probabilites
        cross$geno[[i]]$prob <- array(z$genoprob,dim=c(n.ind,n.pos,n.gen))
        dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, gen.names)
        # attribute set to the error.prob value used, for later
        #     reference, especially by calc.errorlod()

        attr(cross$geno[[i]]$prob, "map") <- map
        attr(cross$geno[[i]]$prob,"error.prob") <- error.prob
        attr(cross$geno[[i]]$prob,"step") <- step
        attr(cross$geno[[i]]$prob,"off.end") <- temp.offend
        attr(cross$geno[[i]]$prob,"map.function") <- map.function
        attr(cross$geno[[i]]$prob,"stepwidth") <- stepwidth
    } # end loop over chromosomes

    # 4- and 8-way RIL: reorganize the results
    if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="bgmagic16")
        cross <- reorgRIgenoprob(cross)

    cross
}

######################################################################
#
# calc.genoprob.special:
#    special version used by calc.errorlod
#    for each individual and marker, calculate probabilities allowing
#    that genotype to be in error but assuming that all other genotypes
#    are correct
#
######################################################################

calc.genoprob.special <-
    function(cross, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"))
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    step <- 0
    off.end <- 0
    stepwidth <- "fixed"

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
    n.mar <- nmar(cross)

    type <- class(cross)[1]

    # calculate genotype probabilities one chromosome at a time
    for(i in 1:n.chr) {
        if(n.mar[i]==1) temp.offend <- max(c(off.end,5))
        else temp.offend <- off.end

        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") xchr <- TRUE
        else xchr <- FALSE

        # which type of cross is this?
        if(type == "f2") {
            one.map <- TRUE
            if(!xchr) { # autosomal
                cfunc <- "calc_genoprob_special_f2"
                n.gen <- 3
                gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
            }
            else {                             # X chromsome
                cfunc <- "calc_genoprob_special_bc"
                n.gen <- 2
                gen.names <- c("g1","g2")
            }
        }
        else if(type == "bc") {
            cfunc <- "calc_genoprob_special_bc"
            n.gen <- 2
            if(!xchr)
                gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
            else gen.names <- c("g1","g2")
            one.map <- TRUE
        }
        else if(type == "riself" || type=="risib" || type=="dh" || type=="haploid") {
            cfunc <- "calc_genoprob_special_bc"
            n.gen <- 2
            gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
            one.map <- TRUE
        }
        else if(type == "4way") {
            cfunc <- "calc_genoprob_special_4way"
            n.gen <- 4
            one.map <- FALSE
            gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
        }
        else if(type=="ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self" || type=="bgmagic16") {
            cfunc <- paste("calc_genoprob_special_", type, sep="")
            if(type=="bgmagic16") n.gen <- 16
            else n.gen <- as.numeric(substr(type, 3, 3))
            one.map <- TRUE
            gen.names <- LETTERS[1:n.gen]
            if(xchr)
                warning("calc.genoprob.special not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else if(type == "bcsft") {
            one.map <- TRUE
            cfunc <- "calc_genoprob_bcsft"
            cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
            if(!xchr) { # autosomal
                if(cross.scheme[2] == 0) {
                    gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
                    n.gen <- 2
                }
                else {
                    gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
                    n.gen <- 3
                }
            }
            else { ## X chr
                cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
                cross.scheme[2] <- 0
                gen.names <- c("g1","g2")
                n.gen <- 2
            }
        }
        else
            stop("calc.genoprob.special not available for cross type ", type, ".")

        # genotype data
        gen <- cross$geno[[i]]$data
        gen[is.na(gen)] <- 0

        # recombination fractions
        if(one.map) {
            # recombination fractions
            map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
            rf <- mf(diff(map))
            if(type=="risib" || type=="riself")
                rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),chrtype)
            rf[rf < 1e-14] <- 1e-14

            # new genotype matrix with pseudomarkers filled in
            newgen <- matrix(ncol=length(map),nrow=nrow(gen))
            dimnames(newgen) <- list(NULL,names(map))
            newgen[,colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- names(map)
        }
        else {
            map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
            rf <- mf(diff(map[1,]))
            rf[rf < 1e-14] <- 1e-14
            rf2 <- mf(diff(map[2,]))
            rf2[rf2 < 1e-14] <- 1e-14

            # new genotype matrix with pseudomarkers filled in
            newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
            dimnames(newgen) <- list(NULL,dimnames(map)[[2]])
            newgen[,colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- colnames(map)
        }

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
                    genoprob=as.double(temp),
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
                    genoprob=as.double(rep(0,n.gen*n.ind*n.pos)),
                    PACKAGE="qtl")
        }

        # re-arrange marginal probabilites
        cross$geno[[i]]$prob <- array(z$genoprob,dim=c(n.ind,n.pos,n.gen))
        dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, gen.names)
        # attribute set to the error.prob value used, for later
        #     reference, especially by calc.errorlod()

        attr(cross$geno[[i]]$prob, "map") <- map
        attr(cross$geno[[i]]$prob,"error.prob") <- error.prob
        attr(cross$geno[[i]]$prob,"step") <- step
        attr(cross$geno[[i]]$prob,"off.end") <- temp.offend
        attr(cross$geno[[i]]$prob,"map.function") <- map.function
        attr(cross$geno[[i]]$prob,"stepwidth") <- stepwidth
    } # end loop over chromosomes

    cross
}

# end of calc.genoprob.R
