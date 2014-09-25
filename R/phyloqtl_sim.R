######################################################################
# phyloqtl_sim.R
#
# copyright (c) 2009, Karl W Broman
# last modified May, 2009
# first written May, 2009
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
# Simulating multiple crosses with multiple taxa with diallelic QTL
# located on a phylogenetic tree
#
# Part of the R/qtl package
# Contains: simPhyloQTL
#
######################################################################

######################################################################
# simPhyloQTL
#
# function for simulating data on multiple intercrosses, for a single
# diallelic QTL
#
# n.taxa =
# partition = character string of the form "AB|CD" or "A|BCD"
#             if missing, simulate under the null
# crosses = set of two-character strings indicating the crosses to do
#             (of the form "AB", AC", etc.)
#            these will be sorted and then only unique ones used
#            if missing, we'll do all crosses
######################################################################

simPhyloQTL <-
    function(n.taxa=3, partition, crosses, map, n.ind=100, model,
             error.prob=0, missing.prob=0, partial.missing.prob=0,
             keep.qtlgeno=FALSE, keep.errorind=TRUE, m=0, p=0,
             map.function=c("haldane", "kosambi", "c-f", "morgan"))
{
    if(n.taxa < 3) stop("Should have 3 or more taxa.")
    if(n.taxa > 26) stop("We can only deal with 3-26 taxa.")

    map.function <- match.arg(map.function)
    if(!missing(model)) {
        if((is.matrix(model) && ncol(model)!=4) ||
           (!is.matrix(model) && length(model) != 4))
            stop("model should have QTL chr, pos, and additive and dominance effects.")
        if(is.matrix(model)) n.qtl <- nrow(model)
        else {
            n.qtl <- 1
            model <- rbind(model)
        }
    }
    else {
        model <- NULL
        n.qtl <- 0
    }

    taxa <- LETTERS[1:n.taxa]

    if(missing(partition) || is.null(partition) || partition=="") {
        partition <- NULL
        if(!is.null(model))
            warning("partition is NULL, so model is ignored and data are simulated with no QTL.")

    }
    else {
        if(is.null(model)) {
            warning("model is NULL, so partition is ignored and data are simulated with no QTL.")
            partition <- NULL
        }
        else {
            if(length(partition)==1 && n.qtl > 1)
                partition <- rep(partition, n.qtl)
            if(length(partition)>1 && n.qtl == 1)
                stop("Model indicates just one QTL, so partition should have length 1.")
            if(length(partition)>1 && n.qtl != length(partition))
                stop("No. QTL in model should match the length of partition.")

            for(i in unique(partition))
                checkPhyloPartition(i, taxa)
        }
    }

    if(missing(crosses) || is.null(crosses)) { # use all possible crosses
        crosses <- NULL
        for(i in 1:(n.taxa-1))
            for(j in (i+1):n.taxa)
                crosses <- c(crosses, paste(taxa[i], taxa[j], sep=""))
    }
    else  # check that the crosses are correct and sufficient
        crosses <- checkPhyloCrosses(crosses, taxa)

    if(length(n.ind) > 1 && length(n.ind) != length(crosses))
        stop("n.ind should have length 1 or have length equal to the number of crosses.")
    if(length(n.ind)==1) n.ind <- rep(n.ind, length(crosses))

    if(n.qtl > 0) {
        # for each partition, determine which crosses have the QTL
        crossmat <- qtlByPartition(crosses, partition)

        models <- vector("list", length(crosses))
        for(i in seq(along=crosses)) {
            models[[i]] <- model
            models[[i]][crossmat[i,] <0,3] <- -models[[i]][crossmat[i,] <0,3]
            models[[i]] <- models[[i]][crossmat[i,]!=0,,drop=FALSE]
        }
    }
    else models <- vector("list", length(crosses))

    thedata <- vector("list", length(crosses))
    names(thedata) <- crosses
    if(missing(map)) stop("Must provide a genetic map")

    for(i in seq(along=crosses))
        thedata[[i]] <- sim.cross(map, models[[i]], n.ind=n.ind[i], type="f2",
                                  error.prob=error.prob, missing.prob=missing.prob,
                                  partial.missing.prob=partial.missing.prob,
                                  keep.qtlgeno = keep.qtlgeno, keep.errorind=keep.errorind,
                                  m=m, p=p, map.function=map.function)
    thedata
}

# end of phyloqtl_sim.R
