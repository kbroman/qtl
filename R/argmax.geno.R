######################################################################
#
# argmax.geno.R
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
# Contains: argmax.geno
#
######################################################################

######################################################################
#
# argmax.geno: Use Viterbi algorithm to find most likely sequence of
#              underlying genotypes, given observed marker data
#
######################################################################

argmax.geno <-
    function(cross, step=0, off.end=0, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             stepwidth=c("fixed", "variable", "max"))
{
    if(!any(class(cross) == "cross"))
        stop("cross should have class \"cross\".")

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

    # loop over chromosomes
    for(i in 1:n.chr) {
        if(n.mar[i]==1) temp.offend <- max(c(off.end,5))
        else temp.offend <- off.end

        chrtype <- class(cross$geno[[i]])
        if(chrtype=="X") xchr <- TRUE
        else xchr <- FALSE

        # which type of cross is this?
        if(type=="f2") {
            one.map <- TRUE
            if(!xchr) # autosomal
                cfunc <- "argmax_geno_f2"
            else                              # X chromsome
                cfunc <- "argmax_geno_bc"
        }
        else if(type=="bc" || type=="dh" || type=="riself" || type=="risib" || type=="haploid") {
            cfunc <- "argmax_geno_bc"
            one.map <- TRUE
        }
        else if(type == "4way") {
            cfunc <- "argmax_geno_4way"
            one.map <- FALSE
        }
        else if(type == "ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self" || type=="bgmagic16") {
            cfunc <- paste("argmax_geno_", type, sep="")
            one.map <- TRUE
            if(xchr)
                warning("argmax.geno not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else if(type == "bcsft") {
            one.map <- TRUE
            cfunc <- "argmax_geno_bcsft"
            cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
            if(xchr) { ## X chr
                cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
                cross.scheme[2] <- 0
            }
        }
        else
            stop("argmax.geno not available for cross type ", type, ".")

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
        }
        if(any(is.na(rf))) { # this occurs when there is only one marker
            rf <- rf2 <- 0
            warn <- paste("Only one marker on chr ", names(cross$geno)[i],
                          ": argmax results tenuous.", sep="")
            warning(warn)
        }


        # call the C function
        if(one.map) {
            ## Hide cross scheme in genoprob to pass to routine. BY
            temp <- newgen
            if(type == "bcsft")
                temp[1:2] <- cross.scheme

            z <- .C(cfunc,
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(newgen),        # genotype data
                    as.double(rf),             # recombination fractions
                    as.double(error.prob),
                    argmax=as.integer(temp), # the output
                    PACKAGE="qtl")

            cross$geno[[i]]$argmax <- matrix(z$argmax,ncol=n.pos)
            dimnames(cross$geno[[i]]$argmax) <- list(NULL, names(map))
        }
        else {
            z <- .C(cfunc,
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(newgen),        # genotype data
                    as.double(rf),             # recombination fractions
                    as.double(rf2),            # recombination fractions
                    as.double(error.prob),
                    argmax=as.integer(newgen), # the output
                    PACKAGE="qtl")

            cross$geno[[i]]$argmax <- matrix(z$argmax,ncol=n.pos)
            dimnames(cross$geno[[i]]$argmax) <- list(NULL, colnames(map))
        }

        # attribute set to the error.prob value used, for later
        #     reference
        attr(cross$geno[[i]]$argmax, "map") <- map
        attr(cross$geno[[i]]$argmax,"error.prob") <- error.prob
        attr(cross$geno[[i]]$argmax,"step") <- step
        attr(cross$geno[[i]]$argmax,"off.end") <- temp.offend
        attr(cross$geno[[i]]$argmax,"map.function") <- map.function
        attr(cross$geno[[i]]$argmax,"stepwidth") <- stepwidth
    }

    # store argmax values as integers
    for(i in 1:nchr(cross))
        storage.mode(cross$geno[[i]]$argmax) <- "integer"

    # 4- and 8-way RIL: reorganize the results
    if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="bgmagic16")
        cross <- reorgRIargmax(cross)

    cross
}

# end of argmax.geno.R
