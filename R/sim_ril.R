#####################################################################
#
# sim_ril.R
#
# copyright (c) 2004-2011, Karl W Broman
# last modified May, 2011
# first written May, 2004
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
# Contains: sim.ril, simulateFounderSnps, convertMWril
#
######################################################################

######################################################################
#
# sim.ril
#
# Simulate RILs by selfing or sibling mating from 2, 4, or 8
# parental strains
# map = map in the usual R/qtl map format
# m = interference parameter (0 is no interference)
######################################################################
sim.ril <-
    function(map, n.ril=1, type=c("sibmating", "selfing"),
             n.str=c("2","4","8"), m=0, p=0,
             error.prob=0, missing.prob=0, random.cross=TRUE)
{
    type <- match.arg(type)
    if(type=="sibmating") selfing <- 0
    else selfing <- 1
    if(is.numeric(n.str)) n.str <- as.character(n.str)
    n.str <- as.numeric(match.arg(n.str))
    n.chr <- length(map)
    n.mar <- sapply(map,length)
    tot.mar <- sum(n.mar)

    if(m < 0) stop("Must have m >= 0.")
    if(p < 0 || p > 1) stop("Must have 0 <= p <= 1.")
    if(p == 1) {
        p <- 0
        m <- 0
    }

    omap <- map
    map <- lapply(map, function(a) a-min(a))

    if(!selfing && class(omap[[length(omap)]])=="X")
        include.x <- TRUE
    else {
        for(i in seq(along=omap)) class(omap[[i]]) <- "A"
        include.x <- FALSE
    }

    if(n.str==2) random.cross <- FALSE

    x <- .C("R_sim_ril",
            as.integer(n.chr),
            as.integer(n.mar),
            as.integer(n.ril),
            as.double(unlist(map)),
            as.integer(n.str),
            as.integer(m),
            as.double(p),
            as.integer(include.x),
            as.integer(random.cross),
            as.integer(selfing),
            cross=as.integer(rep(0,n.ril*n.str)),
            res=as.integer(rep(0,tot.mar*n.ril)),
            orig=as.integer(rep(0,tot.mar*n.ril)),
            as.double(error.prob),
            as.double(missing.prob),
            err=as.integer(rep(0,tot.mar*n.ril)),
            PACKAGE="qtl")

    cross <- t(matrix(x$cross,ncol=n.ril,nrow=n.str))
    err <- t(matrix(x$err,nrow=tot.mar,ncol=n.ril))
    truegeno <- t(matrix(x$orig, nrow=tot.mar, ncol=n.ril))
    x <- t(matrix(x$res,nrow=tot.mar,ncol=n.ril))
    x[x==0] <- NA

    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    cur <- 0
    for(i in 1:n.chr) {
        geno[[i]]$data <- x[,cur + 1:n.mar[i],drop=FALSE]
        colnames(geno[[i]]$data) <- names(map[[i]])
        geno[[i]]$map <- omap[[i]]
        if(missing.prob > 0 || (error.prob>0 && n.str==2))
            geno[[i]]$truegeno <- truegeno[,cur+1:n.mar[i],drop=FALSE]
        if(error.prob > 0 && n.str==2)
            geno[[i]]$errors <- err[,cur+1:n.mar[i],drop=FALSE]

        cur <- cur + n.mar[i]
        class(geno[[i]]) <- class(omap[[i]])

    }
    pheno <- data.frame(line=1:n.ril, stringsAsFactors=TRUE)
    x <- list(geno=geno,pheno=pheno,cross=cross)

    # ri[n][sib/self]un: un = genotypes not yet transformed
    if(type=="sibmating") {
        if(n.str=="2")
            class(x) <- c("risib","cross")
        else
            class(x) <- c(paste("ri", n.str, "sibun",sep=""),"cross")
    }
    else {
        if(n.str=="2")
            class(x) <- c("riself","cross")
        else
            class(x) <- c(paste("ri", n.str, "selfun",sep=""),"cross")
    }

    x
}

######################################################################
# simFounderSnps
#
# Simulate founder snp genotypes for a multiple-strain RIL
#
# map = genetic map of markers (used just to get no. markers per chr)
#
# n.str = number of founder strains (4 or 8)
#
# pat.freq = frequency of SNP genotype patterns (length n.str/2 + 1)
#            (monoallelic, snp unique to a founder,
#             snp present in 2 founder,
#             [for 8 founders: snp in 3/8, snp in 4/8] )
######################################################################
simFounderSnps <-
    function(map, n.str=c("4","8"), pat.freq)
{
    if(is.numeric(n.str)) n.str <- as.character(n.str)
    n.str <- as.numeric(match.arg(n.str))

    if(missing(pat.freq)) {
        if(n.str==8) pat.freq <- c(0, 0.4, 0.3, 0.2, 0.1)
        else pat.freq <- c(0, 0.7, 0.3)
    }

    if(length(pat.freq) < n.str/2+1)
        pat.freq <- c(pat.freq, rep(0, n.str/2+1 - length(pat.freq)))
    else pat.freq <- pat.freq[1:(n.str/2+1)]
    pat.freq <- pat.freq/sum(pat.freq)

    n.mar <- sapply(map, length)
    output <- vector("list", length(map))
    names(output) <- names(map)
    for(i in seq(along=map)) {
        thepat <- sample(seq(length(pat.freq))-1, n.mar[i], prob=pat.freq, replace=TRUE)
        output[[i]] <- matrix(0, ncol=n.str, nrow=n.mar[i])
        for(j in seq(along=thepat))
            output[[i]][j,sample(1:n.str, thepat[j])] <- 1
    }

    output
}

######################################################################
# convertMWril: Convert multiple-strain RIL genotypes using parental data
#
# parents = Parental genotype data, with genetic map
#           list with elements being chromosomes
#           each chromosome is a matrix n.mar x n.str,
######################################################################
convertMWril <-
    function(cross, parents, error.prob=0)
{
    crosstype <- class(cross)[1]
    n.str.by.crosstype <- as.numeric(substr(crosstype, 3, 3))
    un <- substr(crosstype, nchar(crosstype)-1, nchar(crosstype))
    if(un != "un")
        stop("cross appears to have already been converted.")
    class(cross)[1] <- substr(crosstype, 1, nchar(crosstype)-2)

    n.ril <- nind(cross)
    thecrosses <- cross$cross
    n.str <- ncol(thecrosses)
    if(n.str != ncol(parents[[1]]))
        stop("Different numbers of founders in cross and parents.")
    if(n.str != n.str.by.crosstype)
        stop("Confusion regarding no. founders within cross.")
    if(length(parents) != nchr(cross))
        stop("Different numbers of chromosomes in cross and parents.")

    n.mar <- nmar(cross)
    n.mar2 <- sapply(parents, nrow)
    if(any(n.mar != n.mar2))
        stop("Different numbers of markers in cross and parents.")

    pg <- unlist(parents)
    if(any(is.na(pg)))
        stop("Missing parental data not allowed.")

    # if positive error prob, check whether all parental data are snps
    if(error.prob > 0) {
        if(all(pg==0 | pg==1)) all.snps <- TRUE
        else {
            if(all(pg==1 | pg==2)) { # convert to 0/1
                parents <- lapply(parents, function(a) a - 1)
                all.snps <- TRUE
            }
            else all.snps <- FALSE
        }
    }
    else all.snps <- FALSE


    for(i in 1:nchr(cross)) {
        dat <- cross$geno[[i]]$data
        dat[is.na(dat)] <- 0

        results <-
            .C("R_convertMWril",
               as.integer(n.ril),        # no. ril
               as.integer(n.mar[i]),     # no. markers
               as.integer(n.str),        # no. founders
               as.integer(parents[[i]]), # SNP data on parents (n.mar x n.str)
               g=as.integer(dat), # SNP data on RIL (n.ril x n.mar)
               as.integer(thecrosses),   # the crosses (n.ril x n.str)
               as.integer(all.snps),
               as.double(error.prob),
               err=as.integer(rep(0,n.mar[i]*n.ril)),
               PACKAGE="qtl")

        # replace 0's with missing values
        newgeno <- results$g
        newgeno[newgeno==0] <- NA
        newgeno <- matrix(newgeno, n.ril, n.mar[i])
        colnames(newgeno) <- colnames(cross$geno[[i]]$data)
        cross$geno[[i]]$data <- newgeno
        if(error.prob > 0)
            cross$geno[[i]]$errors <- matrix(results$err, n.ril, n.mar[i])
    }

    cross
}

# end of sim_ril.R
