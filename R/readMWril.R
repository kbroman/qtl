######################################################################
#
# readMWril.R
#
# copyright (c) 2009-2014, Karl W Broman
# last modified Jun, 2014
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
# Contains: readMWril
#
######################################################################

######################################################################
#
# readMWril
#
# read multi-way RIL data in comma-delimited format,
# with a separate file for the founder genotypes, and possible a
# separate file for the phenotype data.
#
######################################################################

readMWril <-
    function(dir, rilfile, founderfile,
             type=c("ri4self", "ri4sib", "ri8self", "ri8selfIRIP1", "ri8sib", "bgmagic16"),
             na.strings=c("-","NA"), rotate=FALSE,
             ...)
{
    # create file names
    if(missing(rilfile) || missing(founderfile))
        stop("Need to specify rilfile and founderfile.")

    if(!missing(dir) && dir != "") {
        rilfile <- file.path(dir, rilfile)
        founderfile <- file.path(dir, founderfile)
    }

    type <- match.arg(type)

    args <- list(...)

    # if user wants to use comma for decimal point, we need
    if(length(args) > 0 && "dec" %in% names(args)) {
        dec <- args[["dec"]]
    }
    else dec <- "."

    # read the data file
    if(length(args) < 1 || !("sep" %in% names(args))) {
        # "sep" not in the "..." argument and so take sep=","
        if(length(args) < 1 || !("comment.char" %in% names(args))) {
            gen <- read.table(rilfile, sep=",", na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, comment.char="", ...)
            founder <- read.table(founderfile, sep=",", na.strings=na.strings,
                                  colClasses="character", fill=TRUE,
                                  blank.lines.skip=TRUE, comment.char="", ...)
        }
        else {
            gen <- read.table(rilfile, sep=",", na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, ...)
            founder <- read.table(founderfile, sep=",", na.strings=na.strings,
                                  colClasses="character", fill=TRUE,
                                  blank.lines.skip=TRUE, ...)
        }
    }
    else {
        if(length(args) < 1 || !("comment.char" %in% names(args))) {
            gen <- read.table(rilfile, na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, comment.char="", ...)
            founder <- read.table(founderfile, na.strings=na.strings,
                                  colClasses="character", fill=TRUE,
                                  blank.lines.skip=TRUE, comment.char="", ...)
        }
        else {
            gen <- read.table(rilfile, na.strings=na.strings,
                              colClasses="character", fill=TRUE,
                              blank.lines.skip=TRUE, ...)
            founder <- read.table(founderfile, na.strings=na.strings,
                                  colClasses="character", fill=TRUE,
                                  blank.lines.skip=TRUE, ...)
        }
    }

    if(rotate) {
        gen <- as.data.frame(t(gen), stringsAsFactors=FALSE)
        founder <- as.data.frame(t(founder), stringsAsFactors=FALSE)
    }
    rn <- founder[-1,1]
    cn <- founder[1,-1]
    founder <- founder[-1,-1, drop=FALSE]
    founder <- matrix(as.numeric(unlist(founder)), ncol=length(cn))
    dimnames(founder) <- list(rn, cn)

    # determine number of phenotypes based on initial blanks in row 2
    n <- ncol(gen)
    temp <- rep(FALSE,n)
    for(i in 1:n) {
        temp[i] <- all(gen[2,1:i]=="")
        if(!temp[i]) break
    }

    if(!any(temp)) # no phenotypes!
        stop("You must include at least one phenotype (e.g., an index).")
    n.phe <- max((1:n)[temp])

    # Is map included?  yes if first n.phe columns in row 3 are all blank
    if(all(!is.na(gen[3,1:n.phe]) & gen[3,1:n.phe]=="")) {
        map.included <- TRUE
        map <- asnumericwithdec(unlist(gen[3,-(1:n.phe)]), dec=dec)
        if(any(is.na(map)))
            stop("There are missing marker positions.")
        nondatrow <- 3
    }
    else {
        map.included <- FALSE
        map <- rep(0,ncol(gen)-n.phe)
        nondatrow <- 2 # last non-data row
    }
    pheno <- as.data.frame(gen[-(1:nondatrow),1:n.phe,drop=FALSE], stringsAsFactors=FALSE)
    colnames(pheno) <- as.character(gen[1,1:n.phe])

    # replace empty cells with NA
    gen <- sapply(gen,function(a) { a[!is.na(a) & a==""] <- NA; a })

    # pull apart phenotypes, genotypes and map
    mnames <- gen[1,-(1:n.phe)]
    if(any(is.na(mnames)))  stop("There are missing marker names.")
    chr <- gen[2,-(1:n.phe)]
    if(any(is.na(chr))) stop("There are missing chromosome IDs.")

    gen <- matrix(as.numeric(gen[-(1:nondatrow),-(1:n.phe)]),
                  ncol=ncol(gen)-n.phe)

    pheno <- data.frame(lapply(pheno, sw2numeric, dec=dec), stringsAsFactors=TRUE)

    n.str <- nrow(founder)

    if(type == "ri8selfIRIP1"){
    }
    else if(!("cross" %in% names(pheno))) {
        warning("Need a phenotype named \"cross\"; assuming all come from the cross ",
                paste(LETTERS[1:n.str], collapse="x"))
        crosses <- matrix(1:n.str, ncol=n.str, nrow=nrow(gen), byrow=TRUE)
    }
    else {
        pheno$cross <- as.character(pheno$cross)
        if(any(nchar(pheno$cross) != n.str))
            stop("Mismatches in length of \"cross\" phenotype.")
        thecross <- matrix(unlist(strsplit(pheno$cross, "")), byrow=TRUE, ncol=n.str)
        crosses <- matrix(NA, ncol=n.str, nrow=nrow(gen))
        for(i in 1:nrow(thecross))
            crosses[i,] <- match(thecross[i,], LETTERS[1:n.str])
        if(any(is.na(crosses)))
            stop("Problems in the \"cross \" phenotype.")
    }

    # check founder data matches in dimension
    if(ncol(founder) != ncol(gen))
        stop("Different numbers of markers in RIL and founder files.")
    if(any(colnames(founder) != mnames)) {
        cnf <- colnames(founder)
        if(any(is.na(match(cnf, mnames))) || any(is.na(match(mnames, cnf))))
            stop("Mismatch in markers in RIL and founder files.")
        founder <- founder[,match(mnames, cnf),drop=FALSE]
    }

    wh <- which(is.na(gen))
    missingval <- min(as.numeric(c(gen, founder)), na.rm=TRUE)-1
    gen[wh] <- missingval
    founder[is.na(founder)] <- missingval
    d <- dim(gen)
    if(type == "ri8selfIRIP1")
    {
        gen <-
            .C("R_reviseMWrilNoCross",
               as.integer(d[1]),
               as.integer(d[2]),
               as.integer(n.str),
               as.integer(founder),
               gen=as.integer(gen),
               as.integer(missingval),
               PACKAGE="qtl")$gen
    }
    else
    {
        gen <-
            .C("R_reviseMWril",
               as.integer(d[1]),
               as.integer(d[2]),
               as.integer(n.str),
               as.integer(founder),
               gen=as.integer(gen),
               as.integer(crosses),
               as.integer(missingval),
               PACKAGE="qtl")$gen
    }
    gen[wh] <- NA
    gen <- matrix(gen, nrow=d[1])
    gen[gen==0 | gen==(2^n.str-1)] <- NA


    # re-order the markers by chr and position
    # try to figure out the chr labels
    if(all(chr %in% c(1:999,"X","x"))) { # 1...19 + X
        tempchr <- chr
        tempchr[chr=="X" | chr=="x"] <- 1000
        tempchr <- as.numeric(tempchr)
        if(map.included) neworder <- order(tempchr, map)
        else neworder <- order(tempchr)

        chr <- chr[neworder]
        map <- map[neworder]
        gen <- gen[,neworder,drop=FALSE]
        mnames <- mnames[neworder]
    }

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
    for(i in 1:n.chr) { # loop over chromosomes
        # create map
        temp.map <- map[chr==uchr[i]]
        names(temp.map) <- mnames[chr==uchr[i]]

        # pull out appropriate portion of genotype data
        data <- gen[,min.mar:(length(temp.map)+min.mar-1),drop=FALSE]
        min.mar <- min.mar + length(temp.map)
        colnames(data) <- names(temp.map)

        geno[[i]] <- list(data=data,map=temp.map)
        if(uchr[i] == "X" || uchr[i] == "x")
            class(geno[[i]]) <- "X"
        else
            class(geno[[i]]) <- "A"
    }

    cross <- list(geno=geno,pheno=pheno)

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

    if(all(is.na(gen)))
        warning("There is no genotype data!\n")

    if(type != "ri8selfIRIP1")
    {
        cross$cross <- crosses
    }

    # save founder genotypes in data
    founder[founder==missingval] <- NA
    ua <- apply(founder, 2, function(a) unique(a[!is.na(a)]))
    nua <- sapply(ua, length)
    if(all(nua <= 2)) { # re-code as 0/1
        for(i in 1:ncol(founder)) {
            if(nua[i]==1)
                founder[!is.na(founder[,i]),i] <- 0
            else if(nua[i]==2) {
                founder[!is.na(founder[,i]) & founder[,i]==ua[[i]][1],i] <- 0
                founder[!is.na(founder[,i]) & founder[,i]==ua[[i]][2],i] <- 1
            }
        }
    }
    cross$founderGeno <- founder

    class(cross) <- c(type, "cross")

    # check data
    summary(cross)

    cross
}


# end of readMWril.R
