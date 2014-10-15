#####################################################################
#
# ril48_reorg.R
#
# copyright (c) 2009-2012, Karl W Broman
# last modified Nov, 2012
# first written Apr, 2009
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
# Contains: reorgRIgenoprob, reorgRIdraws, reorgRIimp, reorgRIpairprob
#
######################################################################

######################################################################
# reorgRIgenoprob
#
# For 4- and 8-way RIL, reorganize the QTL genotype probabilities
# using the information on the order of the founder strains in each
# cross.
######################################################################

reorgRIgenoprob <-
    function(cross)
{
    crosses <- cross$cross
    flag <- 0
    for(i in 1:ncol(crosses)) {
        if(any(crosses[,i] != i)) {
            flag <- 1
            break
        }
    }
    if(!flag) return(cross) # no need to reorder

    crosstype <- class(cross)[1]
    if(crosstype != "ri4sib" && crosstype != "ri4self" &&
       crosstype != "ri8sib" && crosstype != "ri8self" && crosstype != "bgmagic16")
        stop("reorgRIgenoprob not appropriate for cross type ", crosstype)

    if(crosstype=="bgmagic16") n.str <- 16
    else n.str <- as.numeric(substr(crosstype, 3, 3))
    n.ind <- nind(cross)
    for(i in names(cross$geno)) { # loop over chromosomes
        chrtype <- class(cross$geno[[i]])
        if(chrtype == "X")
            warning("reorgRIgenoprob not working properly for the X chromosome.")

        if(!("prob" %in% names(cross$geno[[i]]))) {
            warning("No QTL genotype probabilities within cross.")
            return(cross)
        }

        prob <- cross$geno[[i]]$prob
        att <- attributes(prob)
        n.mar <- dim(prob)[2]
        if(dim(prob)[1] != n.ind)
            stop("Mismatch between no. individuals in cross and in genoprobs.")
        if(dim(crosses)[2] != n.str)
            stop("Invalid no. of founder strains specified")
        if(dim(prob)[3] != n.str) {
            warning("Odd no. columns in genoprobs for chromosome ", i)
            next
        }

        prob <- .C("R_reorgRIgenoprob",
                   as.integer(n.ind),
                   as.integer(n.mar),
                   as.integer(n.str),
                   prob=as.double(prob),
                   as.integer(crosses),
                   PACKAGE="qtl")$prob
        prob <- array(prob, dim=c(n.ind, n.mar, n.str))
        for(j in names(att))
            attr(prob, j) <- att[[j]]
        cross$geno[[i]]$prob <- prob
    }

    cross
}


######################################################################
# reorgRIdraws
#
# For 4- and 8-way RIL, reorganize the imputed QTL genotypes
# using the information on the order of the founder strains in each
# cross.
######################################################################

reorgRIdraws <-
    function(cross)
{
    crosses <- cross$cross
    flag <- 0
    for(i in 1:ncol(crosses)) {
        if(any(crosses[,i] != i)) {
            flag <- 1
            break
        }
    }
    if(!flag) return(cross) # no need to reorder

    crosstype <- class(cross)[1]
    if(crosstype != "ri4sib" && crosstype != "ri4self" &&
       crosstype != "ri8sib" && crosstype != "ri8self" && crosstype != "bgmagic16")
        stop("reorgRIdraws not appropriate for cross type ", crosstype)

    if(crosstype=="bgmagic16") n.str <- 16
    else n.str <- as.numeric(substr(crosstype, 3, 3))
    n.ind <- nind(cross)
    for(i in names(cross$geno)) { # loop over chromosomes
        chrtype <- class(cross$geno[[i]])
        if(chrtype == "X")
            warning("reorgRIdraws not working properly for the X chromosome.")

        if(!("draws" %in% names(cross$geno[[i]]))) {
            warning("No imputations within cross.")
            return(cross)
        }

        draws <- cross$geno[[i]]$draws
        att <- attributes(draws)
        n.mar <- dim(draws)[2]
        n.imp <- dim(draws)[3]
        if(dim(draws)[1] != n.ind)
            stop("Mismatch between no. individuals in cross and in draws.")
        if(dim(crosses)[2] != n.str)
            stop("Invalid no. of founder strains specified")

        draws <- .C("R_reorgRIdraws",
                    as.integer(n.ind),
                    as.integer(n.mar),
                    as.integer(n.str),
                    as.integer(n.imp),
                    draws=as.integer(draws),
                    as.integer(crosses),
                    PACKAGE="qtl")$draws
        draws <- array(draws, dim=c(n.ind, n.mar, n.imp))
        for(j in names(att))
            attr(draws, j) <- att[[j]]
        cross$geno[[i]]$draws <- draws
    }

    cross
}

######################################################################
# reorgRIargmax
#
# For 4- and 8-way RIL, reorganize the results of argmax.geno
# using the information on the order of the founder strains in each
# cross.
######################################################################

reorgRIargmax <-
    function(cross)
{
    crosses <- cross$cross
    flag <- 0
    for(i in 1:ncol(crosses)) {
        if(any(crosses[,i] != i)) {
            flag <- 1
            break
        }
    }
    if(!flag) return(cross) # no need to reorder

    crosstype <- class(cross)[1]
    if(crosstype != "ri4sib" && crosstype != "ri4self" &&
       crosstype != "ri8sib" && crosstype != "ri8self" && crosstype != "bgmagic16")
        stop("reorgRIargmax not appropriate for cross type ", crosstype)

    if(crosstype=="bgmagic16") n.str <- 16
    else n.str <- as.numeric(substr(crosstype, 3, 3))
    n.ind <- nind(cross)
    for(i in names(cross$geno)) { # loop over chromosomes
        chrtype <- class(cross$geno[[i]])
        if(chrtype == "X")
            warning("reorgRIargmax not working properly for the X chromosome.")

        if(!("argmax" %in% names(cross$geno[[i]]))) {
            warning("No argmax.geno results within cross.")
            return(cross)
        }

        argmax <- cross$geno[[i]]$argmax
        att <- attributes(argmax)
        n.mar <- dim(argmax)[2]
        if(dim(argmax)[1] != n.ind)
            stop("Mismatch between no. individuals in cross and in argmax.")
        if(dim(crosses)[2] != n.str)
            stop("Invalid no. of founder strains specified")

        argmax <- .C("R_reorgRIdraws",
                     as.integer(n.ind),
                     as.integer(n.mar),
                     as.integer(n.str),
                     as.integer(1),
                     argmax=as.integer(argmax),
                     as.integer(crosses),
                     PACKAGE="qtl")$argmax
        argmax <- matrix(argmax, nrow=n.ind, ncol=n.mar)
        for(j in names(att))
            attr(argmax, j) <- att[[j]]
        cross$geno[[i]]$argmax <- argmax
    }

    cross
}

######################################################################
# reorgRIpairprob
#
# For 4- and 8-way RIL, reorganize the results of calc.pairprob
# using the information on the order of the founder strains in each
# cross.
######################################################################

reorgRIpairprob <-
    function(cross, pairprob)
{
    crosses <- cross$cross
    flag <- 0
    for(i in 1:ncol(crosses)) {
        if(any(crosses[,i] != i)) {
            flag <- 1
            break
        }
    }
    if(!flag) return(pairprob) # no need to reorder

    crosstype <- class(cross)[1]
    if(crosstype != "ri4sib" && crosstype != "ri4self" &&
       crosstype != "ri8sib" && crosstype != "ri8self" && crosstype != "bgmagic16")
        stop("reorgRIargmax not appropriate for cross type ", crosstype)

    if(crosstype=="bgmagic16") n.str <- 16
    else n.str <- as.numeric(substr(crosstype, 3, 3))
    n.ind <- nind(cross)
    thedim <- dim(pairprob)

    chrtype <- class(cross$geno[[1]])
    if(chrtype == "X")
        warning("reorgRIpairprob not working properly for the X chromosome.")

    att <- attributes(pairprob)
    if(thedim[1] != n.ind)
        stop("Mismatch between no. individuals in cross and in pairprob.")

    if(thedim[3] != n.str || thedim[4] != n.str)
        stop("Mismatch between no. founder strains in cross and in pairprob.")

    if(dim(crosses)[2] != n.str)
        stop("Invalid no. of founder strains specified")
    n.mar <- nmar(cross)[1]
    if(n.mar*(n.mar-1)/2 != thedim[2])
        stop("Mismatch between no. markers in cross and in pairprob.")

    pairprob <- .C("R_reorgRIpairprob",
                   as.integer(n.ind),
                   as.integer(n.mar), # no. prob
                   as.integer(n.str),
                   pairprob=as.double(pairprob),
                   as.integer(crosses),
                   PACKAGE="qtl")$pairprob
    pairprob <- array(pairprob, dim=thedim)
    for(j in names(att))
        attr(pairprob, j) <- att[[j]]

    pairprob
}

# end of ril48_reorg.R
