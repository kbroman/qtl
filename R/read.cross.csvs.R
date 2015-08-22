######################################################################
#
# read.cross.csvs.R
#
# copyright (c) 2005-2015, Karl W Broman
# last modified Aug, 2015
# first written Oct, 2005
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
# Contains: read.cross.csvs
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.csvs
#
# read data in comma-delimited format, with separate files for phenotype
# and genotype data
#
######################################################################

read.cross.csvs <-
    function(dir, genfile, phefile, na.strings=c("-","NA"),
             genotypes=c("A","H","B","D","C"),
             estimate.map=TRUE, rotate=FALSE, ...)
{
    # create file names
    if(missing(genfile)) genfile <- "gen.csv"
    if(missing(phefile)) phefile <- "phe.csv"

    if(!missing(dir) && dir != "") {
        genfile <- file.path(dir, genfile)
        phefile <- file.path(dir, phefile)
    }

    args <- list(...)

    if("" %in% na.strings) {
        na.strings <- na.strings[na.strings != ""]
        warning("Including \"\" in na.strings will cause problems; omitted.")
    }

    # if user wants to use comma for decimal point, we need
    if(length(args) > 0 && "dec" %in% names(args)) {
        dec <- args[["dec"]]
    }
    else dec <- "."

    # read the data file
    if(length(args) < 1 || !("sep" %in% names(args))) {
        # "sep" not in the "..." argument and so take sep=","
        if(length(args) < 1 || !("comment.char" %in% names(args))) {
            gen <- read.table(genfile, sep=",", na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, comment.char="", ...)
            pheno <- read.table(phefile, sep=",", na.strings=na.strings,
                                colClasses="character", fill=TRUE,
                                blank.lines.skip=TRUE, comment.char="", ...)
        }
        else {
            gen <- read.table(genfile, sep=",", na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, ...)
            pheno <- read.table(phefile, sep=",", na.strings=na.strings,
                                colClasses="character", fill=TRUE,
                                blank.lines.skip=TRUE, ...)
        }
    }
    else {
        if(length(args) < 1 || !("comment.char" %in% names(args))) {
            gen <- read.table(genfile, na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, comment.char="", ...)
            pheno <- read.table(phefile, na.strings=na.strings,
                                colClasses="character", fill=TRUE,
                                blank.lines.skip=TRUE, comment.char="", ...)
        }
        else {
            gen <- read.table(genfile, na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, ...)
            pheno <- read.table(phefile, na.strings=na.strings,
                                colClasses="character", fill=TRUE,
                                blank.lines.skip=TRUE, ...)
        }
    }

    if(rotate) {
        gen <- as.data.frame(t(gen), stringsAsFactors=FALSE)
        pheno <- as.data.frame(t(pheno), stringsAsFactors=FALSE)
    }

    # We must make the first column have the individual IDs
    indname <- gen[1,1]
    if(gen[3,1] == "") {
        genind <- gen[-(1:3),1]
        map.included <- TRUE
        nondatrow <- 3 # last non-data row
    }
    else {
        genind <- gen[-(1:2),1]
        map.included <- FALSE
        nondatrow <- 2 # last non-data row
    }
    gen <- gen[,-1,drop=FALSE]

    wh <- which(pheno[1,] == indname)
    if(length(wh) < 1)
        stop("Can't find the individual ID column (expected '", indname,
             "') in the phenotype file.")

    pheind <- pheno[-1,wh[1]]

    if(length(genind) == length(pheind) && all(genind == pheind)) {
        if(length(genind) != length(unique(genind)))
            warning("Duplicate individual IDs")
    }
    else {

        if(any(is.na(genind)) && any(is.na(pheind)))
            stop("There are missing genotype and phenotype IDs")
        else if(any(is.na(genind)))
            stop("There are missing genotype IDs")
        else if(any(is.na(pheind)))
            stop("There are missing phenotype IDs")

        if(length(genind) != length(unique(genind)) && length(pheind) != length(unique(pheind)))
            stop("There are duplicate genotype and phenotype IDs, and they don't all line up.")
        else if(length(genind) != length(unique(genind)))
            stop("There are duplicate genotype IDs, and the genotype and phenotype IDs don't all line up.")
        else if(length(pheind) != length(unique(pheind)))
            stop("There are duplicate phenotype IDs, and the genotype and phenotype IDs don't all line up.")

        mgp <- match(genind, pheind)

        if(any(is.na(mgp))) { # individuals with genotypes but no phenotypes
            n.add <- sum(is.na(mgp))
            pheind <- c(pheind, genind[is.na(mgp)])
            pheno <- rbind(pheno, matrix(rep(NA, n.add*ncol(pheno)), ncol=ncol(pheno)))
            pheno[-1, wh[1]] <- pheind

            warning(n.add, " individuals with genotypes but no phenotypes\n    ", paste(genind[is.na(mgp)], collapse="|"), "\n")
        }

        mpg <- match(pheind, genind)
        if(any(is.na(mpg))) {
            n.add <- sum(is.na(mpg))
            genind <- c(genind, pheind[is.na(mpg)])

            genadd <- matrix(rep(NA, n.add*ncol(gen)), ncol=ncol(gen))
            colnames(genadd) <- colnames(gen)
            gen <- rbind(gen, genadd)

            warning(n.add, " individuals with phenotypes but no genotypes\n    ", paste(pheind[is.na(mpg)], collapse="|"), "\n")
        }

        mgp <- match(genind, pheind)
        pheind <- pheind[mgp]
        pheno <- pheno[1+c(0,mgp),,drop=FALSE]
    }
    n.phe <- ncol(pheno)

    if(map.included) {
        map <- asnumericwithdec(unlist(gen[3,]), dec=dec)
        if(any(is.na(map))) {
            temp <- unique(unlist(gen[3,])[is.na(map)])
            stop("There are missing marker positions.\n",
                 "   In particular, we see these value(s): ",
                 paste("\"",paste(temp,collapse="\",\"",sep=""),"\"",collapse=" ",sep=""),
                 " at position(s): ",
                 paste(which(is.na(map)),colapse=",",sep=""),sep="")
        }
    }
    else
        map <- rep(0,ncol(gen))

    colnames(pheno) <- unlist(pheno[1,])
    pheno <- apply(pheno, 2, function(a) { a[!is.na(a) & a==""] <- NA; a })
    pheno <- as.data.frame(pheno[-1,], stringsAsFactors=TRUE)

    # replace empty cells with NA
    gen <- sapply(gen,function(a) { a[!is.na(a) & a==""] <- NA; a })

    # pull apart phenotypes, genotypes and map
    mnames <- unlist(gen[1,])
    if(any(is.na(mnames)))
        stop("There are missing marker names. Check column(s) ",paste(which(is.na(mnames))+1+n.phe,collapse=","),sep="")
    chr <- unlist(gen[2,])
    if(any(is.na(chr)))
        stop("There are missing chromosome IDs. Check column(s) ",paste(which(is.na(chr))+1+n.phe,collapse=","),sep="")

    if(any(is.na(chr))) {
        na.positions <- which(is.na(chr))
        na.positions.str <- ""
        if (length(na.positions)<10) {
            na.positions.str <- paste(" at position(s) ",
                                      paste(na.positions,collapse=",",sep=""),sep="")
        }
        stop("There are ", length(na.positions), " missing chromosome IDs",
             na.positions.str, ".")
    }

    # look for strange entries in the genotype data
    if(length(genotypes) > 0) {
        temp <- unique(as.character(gen[-(1:nondatrow),,drop=FALSE]))
        temp <- temp[!is.na(temp)]
        wh <- !(temp %in% genotypes)
        if(any(wh)) {
            warn <- "The following unexpected genotype codes were treated as missing.\n    "
            ge <- paste("|", paste(temp[wh],collapse="|"),"|",sep="")
            warn <- paste(warn,ge,"\n",sep="")
            warning(warn)
        }

        # convert genotype data
        allgeno <- matrix(match(gen[-(1:nondatrow),,drop=FALSE],genotypes),
                          ncol=ncol(gen))
    }
    else
        allgeno <- matrix(as.numeric(gen[-(1:nondatrow),,drop=FALSE]),
                          ncol=ncol(gen))

    pheno <- data.frame(lapply(pheno, sw2numeric, dec=dec), stringsAsFactors=TRUE)

    # re-order the markers by chr and position
    # try to figure out the chr labels
    if(all(chr %in% c(1:999,"X","x"))) { # 1...19 + X
        tempchr <- chr
        tempchr[chr=="X" | chr=="x"] <- 1000
        tempchr <- as.numeric(tempchr)
        if(map.included) neworder <- order(tempchr, map)
        else neworder <- order(tempchr)
    }
    else {
        # don't let it reorder the chromosomes
        tempchr <- factor(chr, levels=unique(chr))

        if(map.included) neworder <- order(tempchr, map)
        else neworder <- order(tempchr)
    }
    chr <- chr[neworder]
    map <- map[neworder]
    allgeno <- allgeno[,neworder,drop=FALSE]
    mnames <- mnames[neworder]

    # fix up dummy map
    if(!map.included) {
        map <- split(rep(0,length(chr)),chr)[unique(chr)]
        map <- unlist(lapply(map,function(a) seq(0,length=length(a),by=5)))
        names(map) <- NULL
    }

    # fix up map information
    # number of chromosomes
    uchr <- unique(chr)
    n.chr <- length(uchr)
    geno <- vector("list",n.chr)
    names(geno) <- uchr
    min.mar <- 1
    allautogeno <- NULL
    for(i in 1:n.chr) { # loop over chromosomes
        # create map
        temp.map <- map[chr==uchr[i]]
        names(temp.map) <- mnames[chr==uchr[i]]

        # pull out appropriate portion of genotype data
        data <- allgeno[,min.mar:(length(temp.map)+min.mar-1),drop=FALSE]
        min.mar <- min.mar + length(temp.map)
        colnames(data) <- names(temp.map)

        geno[[i]] <- list(data=data,map=temp.map)
        if(uchr[i] == "X" || uchr[i] == "x")
            class(geno[[i]]) <- "X"
        else {
            class(geno[[i]]) <- "A"
            if(is.null(allautogeno)) allautogeno <- data
            else allautogeno <- cbind(allautogeno,data)
        }
    }

    if(is.null(allautogeno)) allautogeno <- allgeno

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

    if(all(is.na(allgeno)))
        warning("There is no genotype data!\n")

    # determine map type: f2 or bc or 4way?
    if(all(is.na(allgeno))) warning("There is no genotype data!\n")
    if(all(is.na(allautogeno)) || max(allautogeno,na.rm=TRUE)<=2) type <- "bc"
    else if(max(allautogeno,na.rm=TRUE)<=5) type <- "f2"
    else type <- "4way"
    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c(type,"cross")

    # check that nothing is strange in the genotype data
    cross.type <- class(cross)[1]
    if(cross.type=="f2") max.gen <- 5
    else if(cross.type=="bc") max.gen <- 2
    else max.gen <- 14

    # check that markers are in proper order
    #     if not, fix up the order
    for(i in 1:n.chr) {
        if(any(diff(cross$geno[[i]]$map)<0)) {
            o <- order(cross$geno[[i]]$map)
            cross$geno[[i]]$map <- cross$geno[[i]]$map[o]
            cross$geno[[i]]$data <- cross$geno[[i]]$data[,o,drop=FALSE]
        }
    }

    # if 4-way cross, make the maps matrices
    if(type=="4way") {
        for(i in 1:n.chr)
            cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
    }

    # estimate genetic map
    if(estimate.map && !map.included) estmap <- TRUE
    else estmap <- FALSE

    # return cross + indicator of whether to run est.map
    list(cross,estmap)
}

# end of read.cross.csvs.R
