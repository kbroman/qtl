######################################################################
# phyloqtl_util.R
#
# copyright (c) 2009-2010, Karl W Broman
# last modified Feb, 2010
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
# Utility functions for the phylo/qtl analyses
#
# Part of the R/qtl package
# Contains: checkPhyloPartition, checkPhyloCrosses, qtlByPartition,
#           flipcross, genAllPartitions, sortPhyloPartitions
#
######################################################################

# check that 'partition' is appropriate for a given set of taxa
checkPhyloPartition <-
    function(partition, taxa)
{
    n.taxa <- length(taxa)

    if(length(partition) > 1) {
        for(i in partition) checkPhyloPartition(i, taxa)
        return(TRUE)
    }

    # check that all the taxa are there and that there is just one "|"
    temp <- unlist(strsplit(partition, ""))
    if(length(temp) != n.taxa+1 || any(is.na(match(c(taxa, "|"), temp))))
        stop("partition is mis-specified; should be a character string with all taxa and one vertical bar (|).")

    # check that "|" is not at one end or the other
    ss <- unlist(strsplit(partition, "\\|"))
    if(length(ss) != 2 || any(nchar(ss)==0))
        stop("partition is mis-specified; vertical bar (|) should be somewhere in the middle.")

    TRUE
}

# check that the crosses are correct
checkPhyloCrosses <-
    function(crosses, taxa, nostop=FALSE)
{
    temp <- strsplit(crosses, "")
    if(any(is.na(match(unlist(temp), taxa)))) {
        temp <- unlist(temp)
        extras <- unique(temp[is.na(match(temp, taxa))])
        if(nostop) return(FALSE)
        stop("Crosses mis-specified; have extra taxa, ", paste(extras, collapse=" "))
    }
    crosses <- sapply(temp, paste, collapse="")
    if(length(crosses) != length(unique(crosses))) {
        warning("Some crosses given multiple times; these are omitted.")
        crosses <- unique(crosses)
        temp <- strsplit(crosses, "")
    }

    # all taxa included in the crosses?
    m <- is.na(match(taxa, unique(unlist(temp))))
    if(any(m)) {
        temp <- paste(taxa[m], collapse=" ")
        if(nostop) return(FALSE)
        stop("Some taxa missing from the crosses: ", temp)
    }

    # check that the crosses connect all n taxa
    thetaxa <- as.list(taxa)
    for(i in seq(along=temp)) {
        wh1 <- which(sapply(thetaxa, function(a,b) any(a==b), temp[[i]][1]))
        wh2 <- which(sapply(thetaxa, function(a,b) any(a==b), temp[[i]][2]))
        if(wh1 != wh2) {
            thetaxa[[wh1]] <- c(thetaxa[[wh1]], thetaxa[[wh2]])
            thetaxa <- thetaxa[-wh2]
        }
    }
    if(length(thetaxa) > 1) {
        if(nostop) return(FALSE)
        stop("The crosses are insufficient; they should connect all taxa.")
    }

    if(nostop) return(TRUE)
    crosses
}


######################################################################
# qtlByPartition
#
# for each cross and each partition, determine which crosses have a
# QTL and whether the alleles need to be swapped
######################################################################
qtlByPartition <-
    function(crosses, partition)
{
    # check for multiple crosses (of the form "AB:AC:AD")
    if(length(grep(":", crosses)) > 0) {
        result <- lapply(strsplit(crosses, ":"), qtlByPartition, partition)
        names(result) <- crosses
        return(result)
    }

    # for each partition, determine which crosses have the QTL
    crossmat <- matrix(NA, ncol=length(partition), nrow=length(crosses))
    crosssplit <- strsplit(crosses, "")
    partitionsplit <- lapply(strsplit(partition, "\\|"), strsplit, "")
    for(i in seq(along=crosses)) {
        for(j in seq(along=partition)) {
            v <- vector("list", 2)
            for(k in 1:2)
                v[[k]] <- sapply(partitionsplit[[j]], function(a,b) match(b, a), crosssplit[[i]][k])
            crossmat[i,j] <- diff(sapply(v, function(a) which(!is.na(a))))
        }
    }
    dimnames(crossmat) <- list(crosses, partition)
    crossmat
}


######################################################################
# flipcross
#
# The goal of this function is to take a QTL cross object (for R/qtl)
# and flip the alleles A <-> B.
#
# For the case of an intercross, the allele codes (and corresponding
# QTL genotype probabilities and/or imputated genotypes) are switched
# as follows:
#
# genotype   old code    new code
#   AA          1           3
#   AB          2           2
#   BB          3           1
#  not BB       4           5
#  not AA       5           4
######################################################################

flipcross <-
    function(cross)
{
    if(!("cross" %in% class(cross)))
        stop("The input should have class 'cross'")
    allowed_crosses <- c("f2", "riself", "risib", "dh", "haploid")
    crosstype <- class(cross)[1]
    if(!(crosstype %in% allowed_crosses))
        stop("The function is not working for cross type ", crosstype)

    chrtype <- sapply(cross$geno, "class")

    # omit X chr
    if(any(chrtype=="X")) {
        cross <- subset(cross, chr = (chrtype != "X"))
        warning("flipcross is not yet working for the X chromosome; X chr omitted from output.")
    }


    if(crosstype == "f2") {

        for(i in seq(along=cross$geno)) {
            nd <- d <- cross$geno[[i]]$data
            nd[!is.na(d) & d==1] <- 3
            nd[!is.na(d) & d==3] <- 1
            nd[!is.na(d) & d==4] <- 5
            nd[!is.na(d) & d==5] <- 4
            cross$geno[[i]]$data <- nd

            if("prob" %in% names(cross$geno[[i]])) {
                theattr <- attributes(cross$geno[[i]]$prob)
                cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,,3:1,drop=FALSE]
                attr(cross$geno[[i]]$prob,"map") <- theattr$map
                attr(cross$geno[[i]]$prob,"error.prob") <- theattr$error.prob
                attr(cross$geno[[i]]$prob,"step") <- theattr$step
                attr(cross$geno[[i]]$prob,"off.end") <- theattr$off.end
                attr(cross$geno[[i]]$prob,"map.function") <- theattr$map.function
                attr(cross$geno[[i]]$prob,"stepwidth") <- theattr$stepwidth
            }

            if("draws" %in% names(cross$geno[[i]])) {
                nd <- d <- cross$geno[[i]]$draws
                nd[d==3] <- 1
                nd[d==1] <- 3
                cross$geno[[i]]$draws <- nd
            }
        }
    }
    else { # riself/risib/dh/haploid

        for(i in seq(along=cross$geno)) {
            nd <- d <- cross$geno[[i]]$data
            nd[!is.na(d) & d==1] <- 2
            nd[!is.na(d) & d==2] <- 1

            if("prob" %in% names(cross$geno[[i]])) {
                theattr <- attributes(cross$geno[[i]]$prob)
                cross$geno[[i]]$prob <- cross$geno[[i]]$prob[,,2:1,drop=FALSE]
                attr(cross$geno[[i]]$prob,"map") <- theattr$map
                attr(cross$geno[[i]]$prob,"error.prob") <- theattr$error.prob
                attr(cross$geno[[i]]$prob,"step") <- theattr$step
                attr(cross$geno[[i]]$prob,"off.end") <- theattr$off.end
                attr(cross$geno[[i]]$prob,"map.function") <- theattr$map.function
                attr(cross$geno[[i]]$prob,"stepwidth") <- theattr$stepwidth
            }

            if("draws" %in% names(cross$geno[[i]])) {
                nd <- d <- cross$geno[[i]]$draws
                nd[d==2] <- 1
                nd[d==1] <- 2
                cross$geno[[i]]$draws <- nd
            }
        }
    }

    if("alleles" %in% names(attributes(cross)))
        attr(cross, "alleles") <- rev(attr(cross, "alleles"))

    cross
}

# generate all possible partitions (except the null)
genAllPartitions <-
    function(n.taxa, taxa)
{

    # Utility function
    #     returns binary representation of 1:(2^n)
    binary.v <-
        function(n)
        {
            x <- 1:(2^n)
            mx <- max(x)
            digits <- floor(log2(mx))
            ans <- 0:(digits-1); lx <- length(x)
            x <- matrix(rep(x,rep(digits, lx)),ncol=lx)
            (x %/% 2^ans) %% 2
        }

    if(missing(taxa))
        taxa <- LETTERS[1:n.taxa]
    else {
        if(missing(n.taxa)) n.taxa <- length(taxa)
        else {
            if(n.taxa != length(taxa))
                stop("n.taxa != length(taxa)")
        }
    }

    mat <- binary.v(n.taxa)
    colsum <- apply(mat, 2, sum)
    mat <- mat[,colsum > 0 & colsum <= floor(n.taxa/2)]

    result <- unique(apply(mat, 2, function(a,b) {
        x <- c(paste(b[a==1],collapse=""), paste(b[a==0],collapse=""))
        if(diff(nchar(x)) == 0) x <- sort(x)
        paste(x, collapse="|")}, taxa))

    sortPhyloPartitions(result, taxa)
}


sortPhyloPartitions <-
    function(partitions, taxa)
{
    if(missing(taxa)) taxa <- LETTERS[1:26]

    thesplit <- strsplit(partitions, "\\|")
    n1 <- sapply(thesplit, function(a) nchar(a[1]))
    c1 <- sapply(thesplit, function(a) a[1])
    c1 <- strsplit(c1, "")
    for(i in seq(along=c1))
        c1[[i]] <- as.numeric(paste(match(c1[[i]], taxa), collapse=""))
    c1 <- unlist(c1)
    partitions[order(n1, c1)]
}

# end of phyloqtl_util.R
