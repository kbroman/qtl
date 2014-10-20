######################################################################
#
# sim.geno.R
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
# Contains: sim.geno
#
######################################################################

######################################################################
#
# sim.geno: simulate from the joint distribution Pr(g | O)
#
######################################################################

sim.geno <-
    function(cross, n.draws=16, step=0, off.end=0, error.prob=0.0001,
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

    # don't let error.prob be exactly zero, just in case
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
            n.gen <- 3
            one.map <- TRUE
            if(!xchr) # autosomal
                cfunc <- "sim_geno_f2"
            else                              # X chromsome
                cfunc <- "sim_geno_bc"
        }
        else if(type == "bc" || type=="riself" || type=="risib" || type=="dh" || type=="haploid") {
            cfunc <- "sim_geno_bc"
            n.gen <- 2
            one.map <- TRUE
        }
        else if(type == "4way") {
            n.gen <- 4
            cfunc <- "sim_geno_4way"
            one.map <- FALSE
        }
        else if(type=="ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self" || type=="bgmagic16") {
            cfunc <- paste("sim_geno_", type, sep="")
            if(type=="bgmagic16") n.gen <- 16
            else n.gen <- as.numeric(substr(type, 3, 3))
            one.map <- TRUE
            if(xchr)
                warning("sim.geno not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else if(type == "bcsft") {
            one.map <- TRUE
            cfunc <- "sim_geno_bcsft"
            cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
            if(!xchr) {# autosomal
                n.gen <- 3 - (cross.scheme[1] == 0)
            }
            else {                             # X chromsome
                cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
                cross.scheme[2] <- 0
                n.gen <- 2
            }
        }
        else
            stop("sim_geno not available for cross type ", type, ".")

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
            rf2 <- mf(diff(map[1,]))
            rf2[rf2 < 1e-14] <- 1e-14

            # new genotype matrix with pseudomarkers filled in
            newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
            dimnames(newgen) <- list(NULL,colnames(map))
            newgen[,colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
        }


        # call C function
        if(one.map) {
            ## Hide cross scheme in genoprob to pass to routine. BY
            temp <- as.integer(rep(0,n.draws*n.ind*n.pos))
            if(type == "bcsft")
                temp[1:2] <- cross.scheme

            z <- .C(cfunc,
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.draws),       # number of simulation replicates
                    as.integer(newgen),        # genotype data
                    as.double(rf),             # recombination fractions
                    as.double(error.prob),     #
                    draws=as.integer(temp),
                    PACKAGE="qtl")

            cross$geno[[i]]$draws <- array(z$draws,dim=c(n.ind,n.pos,n.draws))
            dimnames(cross$geno[[i]]$draws) <- list(NULL, names(map), NULL)
        }
        else {
            z <- .C(cfunc,
                    as.integer(n.ind),         # number of individuals
                    as.integer(n.pos),         # number of markers
                    as.integer(n.draws),       # number of simulation replicates
                    as.integer(newgen),        # genotype data
                    as.double(rf),             # recombination fractions
                    as.double(rf2),            # recombination fractions
                    as.double(error.prob),     #
                    draws=as.integer(rep(0,n.draws*n.ind*n.pos)),
                    PACKAGE="qtl")

            cross$geno[[i]]$draws <- array(z$draws,dim=c(n.ind,n.pos,n.draws))
            dimnames(cross$geno[[i]]$draws) <- list(NULL, colnames(map), NULL)

        }

        # attribute set to the error.prob value used, for later
        #     reference
        attr(cross$geno[[i]]$draws, "map") <- map
        attr(cross$geno[[i]]$draws,"error.prob") <- error.prob
        attr(cross$geno[[i]]$draws,"step") <- step
        attr(cross$geno[[i]]$draws,"off.end") <- temp.offend
        attr(cross$geno[[i]]$draws,"map.function") <- map.function
        attr(cross$geno[[i]]$draws,"stepwidth") <- stepwidth
    }

    # store simulated genotypes as integers
    for(i in 1:nchr(cross))
        storage.mode(cross$geno[[i]]$draws) <- "integer"

    return(cross)

    # 4- and 8-way RIL: reorganize the results
    if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib" || type=="bgmagic16")
        cross <- reorgRIdraws(cross)

    cross
}

# end of sim.geno.R
