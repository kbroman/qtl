######################################################################
#
# read.cross.karl.R
#
# copyright (c) 2000-2011, Karl W Broman
# last modified May, 2011
# first written Aug, 2000
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
# Contains: read.cross.karl
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.karl
#
# read data in Karl's format
#
######################################################################

read.cross.karl <-
    function(dir,genfile,mapfile,phefile)
{
    # create file names
    if(missing(genfile)) genfile <- "gen.txt"
    if(missing(mapfile)) mapfile <- "map.txt"
    if(missing(phefile)) phefile <- "phe.txt"

    if(!missing(dir) && dir != "") {
        genfile <- file.path(dir, genfile)
        mapfile <- file.path(dir, mapfile)
        phefile <- file.path(dir, phefile)
    }

    # read data
    geno <- as.matrix(read.table(genfile,na.strings="0"))
    pheno <- as.matrix(read.table(phefile,na.strings="-",header=TRUE))
    tempmap <- scan(mapfile, what=character(),quiet=TRUE)

    # fix up map information
    # number of chromosomes
    n.chr <- as.numeric(tempmap[1])
    n.mar <- 1:n.chr
    g <- map <- geno.data <- vector("list", n.chr)
    cur <- 2
    min.mar <- 1
    names(g) <- as.character(1:n.chr)
    for(i in 1:n.chr) { # loop over chromosomes
        # number of markers
        n.mar[i] <- as.numeric(tempmap[cur])
        cur <- cur+1

        # pull out appropriate portion of genotype data
        geno.data[[i]] <- geno[,min.mar:(min.mar+n.mar[i]-1)]
        min.mar <- min.mar + n.mar[i]

        # recombination fractions
        r <- as.numeric(tempmap[cur:(cur+n.mar[i]-2)])

        # convert to cM distances (w/ Kosambi map function)
        d <- 0.25*log((1+2*r)/(1-2*r))*100

        # convert to locations
        map[[i]] <- round(c(0,cumsum(d)),2)
        cur <- cur+n.mar[i]-1

        # marker names
        names(map[[i]]) <- tempmap[cur:(cur+n.mar[i]-1)]
        dimnames(geno.data[[i]]) <- list(NULL, names(map[[i]]))
        cur <- cur+n.mar[i]

        g[[i]] <- list(data=geno.data[[i]],map=map[[i]])

        # attempt to pull out chromosome number
        mar.names <- names(map[[i]])
        twodig <- grep("[Dd][1-9][0-9][Mm]", mar.names)
        onedig <- grep("[Dd][1-9][Mm]", mar.names)
        xchr <- grep("[Dd][Xx][Mm]", mar.names)

        chr.num <- NULL
        if(length(twodig) > 0)
            chr.num <- c(chr.num,substr(mar.names[twodig],2,3))
        if(length(onedig) > 0)
            chr.num <- c(chr.num,substr(mar.names[onedig],2,2))
        if(length(xchr) > 0)
            chr.num <- c(chr.num,rep("X",length(xchr)))

        # no marker names of the form above
        if(is.null(chr.num)) {
            chr.num <- length(mar.names)
            names(chr.num) <- "1"
        }
        else {
            chr.num <- table(chr.num)
        }

        m <- max(chr.num)
        if(m > sum(chr.num)/2 && m > 1)
            names(g)[i] <- names(chr.num)[chr.num==m][1]

        if(names(g)[i] == "X" || names(g)[i] == "x") class(g[[i]]) <- "X"
        else class(g[[i]]) <- "A"
    }

    # check that data dimensions match
    n.mar1 <- sapply(g,function(a) ncol(a$data))
    n.mar2 <- sapply(g,function(a) length(a$map))
    n.phe <- ncol(pheno)
    n.ind1 <- nrow(pheno)
    n.ind2 <- sapply(g,function(a) nrow(a$data))
    if(any(n.ind1 != n.ind2)) {
        print(c(n.ind1,n.ind2))
        stop("Number of individuals in genotypes and phenotypes do not match.");
    }
    if(any(n.mar1 != n.mar2)) {
        print(c(n.mar,n.mar2))
        stop("Numbers of markers in genotypes and marker names files do not match.");
    }

    # print some information about the amount of data read
    cat(" --Read the following data:\n");
    cat("\t", n.ind1, " individuals\n");
    cat("\t", sum(n.mar1), " markers\n");
    cat("\t", n.phe, " phenotypes\n");

    # add phenotype names, if missing
    if(is.null(colnames(pheno)))
        dimnames(pheno) <- list(NULL, paste("phenotype", 1:n.phe,sep=""))

    # determine map type: f2 or bc or 4way?
    if(max(geno[!is.na(geno)])<=2) type <- "bc"
    else if(max(geno[!is.na(geno)])<=5) type <- "f2"
    else type <- "4way"
    cross <- list(geno=g,pheno=pheno)
    class(cross) <- c(type,"cross")

    # check that nothing is strange in the genotype data
    cross.type <- class(cross)[1]
    if(cross.type=="f2") max.gen <- 5
    else if(cross.type=="bc") max.gen <- 2
    else max.gen <- 14

    u <- unique(geno)
    if(any(!is.na(u) & (u > max.gen | u < 1)))
        stop("There are stange values in the genotype data : ",
             paste(u,collapse=":"), ".")

    cross$pheno <- as.data.frame(cross$pheno, stringsAsFactors=TRUE)

    # return cross + indicator of whether to run est.map
    list(cross,FALSE)
}

# end of read.cross.karl.R
