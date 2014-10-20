######################################################################
#
# read.cross.gary.R
#
# copyright (c) 2000-2014, Karl W Broman
# last modified Jun, 2014
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
# Contains: read.cross.gary
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.gary
#
# read data in Gary's format
#
######################################################################

read.cross.gary <-
    function(dir,genfile,mnamesfile,chridfile,phefile,pnamesfile,mapfile,
             estimate.map,na.strings)
{
    # create file names
    if(missing(genfile)) genfile <- "geno.dat"
    if(missing(mnamesfile)) mnamesfile <- "mnames.txt"
    if(missing(chridfile)) chridfile <- "chrid.dat"
    if(missing(phefile)) phefile <- "pheno.dat"
    if(missing(pnamesfile)) pnamesfile <- "pnames.txt"
    if(missing(mapfile)) mapfile <- "markerpos.txt"

    if(!missing(dir) && dir != "") {
        genfile <- file.path(dir, genfile)
        mnamesfile <- file.path(dir, mnamesfile)
        chridfile <- file.path(dir, chridfile)
        phefile <- file.path(dir, phefile)
        if(!is.null(pnamesfile)) pnamesfile <- file.path(dir, pnamesfile)
        if(!is.null(mapfile)) mapfile <- file.path(dir, mapfile)
    }

    # read data
    allgeno <- as.matrix(read.table(genfile,na.strings="9"))+1
    pheno <- as.matrix(read.table(phefile,na.strings=na.strings,header=FALSE))
    chr <- scan(chridfile,what=character(),quiet=TRUE)
    mnames <- scan(mnamesfile,what=character(),quiet=TRUE)

    if(!is.null(mapfile)) {
        map <- read.table(mapfile,row.names=1)
        map <- map[mnames,1]
        map.included <- TRUE
    }
    else {
        map <- seq(0,by=5,len=length(mnames))
        map.included <- FALSE
    }

    if(!is.null(pnamesfile)) pnames <- scan(pnamesfile,what=character(),quiet=TRUE)
    else pnames <- paste("pheno",1:ncol(pheno),sep="")


    # fix up map information
    # number of chromosomes
    uchr <- unique(chr)
    n.chr <- length(uchr)
    geno <- vector("list",n.chr)
    names(geno) <- uchr
    min.mar <- 1
    for(i in 1:n.chr) { # loop over chromosomes
        # create map
        temp.map <- map[chr==uchr[i]]

        # deal with any markers that didn't appear in the marker pos file
        if(any(is.na(temp.map))) {
            o <- (seq(along=temp.map))[is.na(temp.map)]
            for(j in o) {
                if(j==1 || all(is.na(temp.map[1:(j-1)]))) {
                    z <- min((seq(along=temp.map))[-o])
                    temp.map[j] <- min(temp.map,na.rm=TRUE)-(z-j+1)
                }
                else if(j==length(temp.map) || all(is.na(temp.map[-(1:j)]))) {
                    z <- max((seq(along=temp.map))[-o])
                    temp.map[j] <- max(temp.map,na.rm=TRUE)+(j-z+1)
                }
                else {
                    temp.map[j] <- (min(temp.map[-(1:j)],na.rm=TRUE)+
                                    max(temp.map[1:(j-1)],na.rm=TRUE))/2
                }
            }
        }

        names(temp.map) <- mnames[chr==uchr[i]]

        # pull out appropriate portion of genotype data
        data <- allgeno[,min.mar:(length(temp.map)+min.mar-1),drop=FALSE]
        min.mar <- min.mar + length(temp.map)
        colnames(data) <- names(temp.map)

        geno[[i]] <- list(data=data,map=temp.map)
        if(uchr[i] == "X" || uchr[i] == "x")
            class(geno[[i]]) <- "X"
        else class(geno[[i]]) <- "A"
    }
    colnames(pheno) <- pnames

    # fix up phenotype data: make things numeric that look numeric
    sw2numeric_gary <-
        function(x) {
            pattern <- "^[ \t]*-*[0-9]*[.]*[0-9]*[ \t]*$"
            n <- sum(!is.na(x))
            if(length(grep(pattern,as.character(x[!is.na(x)])))==n)
                return(as.numeric(as.character(x)))
            else return(x)
        }
    pheno <- data.frame(lapply(as.data.frame(pheno), sw2numeric_gary), stringsAsFactors=TRUE)

    # check that data dimensions match
    n.mar1 <- sapply(geno,function(a) ncol(a$data))
    n.mar2 <- sapply(geno,function(a) length(a$map))
    n.phe <- ncol(pheno)
    n.ind1 <- nrow(pheno)
    n.ind2 <- sapply(geno,function(a) nrow(a$data))
    if(any(n.ind1 != n.ind2)) {
        cat(n.ind1,n.ind2,"\n")
        stop("Number of individuals in genotypes and phenotypes do not match.");
    }
    if(any(n.mar1 != n.mar2)) {
        cat(n.mar1,n.mar2,"\n")
        stop("Numbers of markers in genotypes and marker names files do not match.");
    }

    # print some information about the amount of data read
    cat(" --Read the following data:\n");
    cat("\t", n.ind1, " individuals\n");
    cat("\t", sum(n.mar1), " markers\n");
    cat("\t", n.phe, " phenotypes\n");

    # determine map type: f2 or bc or 4way?
    if(max(allgeno[!is.na(allgeno)])<=2) type <- "bc"
    else type <- "f2"
    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c(type,"cross")

    # check that nothing is strange in the genotype data
    cross.type <- class(cross)[1]
    if(cross.type=="f2") max.gen <- 5
    else max.gen <- 2

    u <- unique(allgeno)
    if(any(!is.na(u) & (u > max.gen | u < 1)))
        stop("There are stange values in the genotype data : ",
             paste(sort(u),collapse=":"), ".")

    cross$pheno <- as.data.frame(cross$pheno, stringsAsFactors=TRUE)

    # if map wasn't included, go through each chromosome and
    # make first marker at 0 cM.
    if(!map.included) {
        for(i in 1:nchr(cross))
            cross$geno[[i]]$map <- cross$geno[[i]]$map - min(cross$geno[[i]]$map)
    }

    # return cross + indicator of whether to run est.map
    # [run est.map if map not included and estimate.map == TRUE]
    list(cross, (!map.included && estimate.map) )
}

# end of read.cross.gary.R
