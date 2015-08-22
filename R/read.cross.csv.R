######################################################################
#
# read.cross.csv.R
#
# copyright (c) 2000-2015, Karl W Broman
# last modified Aug, 2015
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
# Contains: read.cross.csv
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.csv
#
# read data in comma-delimited format
#
######################################################################

read.cross.csv <-
    function(dir, file, na.strings=c("-","NA"),
             genotypes=c("A","H","B","D","C"),
             estimate.map=TRUE, rotate=FALSE, ...)
{
    # create file names
    if(missing(file)) file <- "data.csv"

    if(!missing(dir) && dir != "") {
        file <- file.path(dir, file)
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
        if(length(args) < 1 || !("comment.char" %in% names(args)))
            data <- read.table(file, sep=",", na.strings=na.strings,
                               colClasses="character", fill=TRUE,
                               blank.lines.skip=TRUE, comment.char="", ...)
        else
            data <- read.table(file, sep=",", na.strings=na.strings,
                               colClasses="character", fill=TRUE,
                               blank.lines.skip=TRUE, ...)
    }
    else {
        if(length(args) < 1 || !("comment.char" %in% names(args)))
            data <- read.table(file, na.strings=na.strings,
                               colClasses="character", fill=TRUE,
                               blank.lines.skip=TRUE, comment.char="", ...)
        else
            data <- read.table(file, na.strings=na.strings,
                               colClasses="character", fill=TRUE,
                               blank.lines.skip=TRUE, ...)
    }

    if(rotate)
        data <- as.data.frame(t(data), stringsAsFactors=FALSE)

    empty <- grep("^\\s*$", data[2, ])

    if( ! 1 %in% empty)
        stop("You must include at least one phenotype (e.g., an index). ",
             "There was this value in the first column of the second row '",
             data[2,1],"' where was supposed to be nothing.",sep="")

    # determine number of phenotypes based on initial blanks in row 2
    if(length(empty)==ncol(data))
        stop("Second row has all blank cells; you need to include chromosome IDs for the markers.")
    n.phe <- min((1:ncol(data))[-empty])-1
    empty <- rep(FALSE, n.phe)
    empty[grep("^\\s*$", data[3,1:n.phe])] <- TRUE

    # Is map included?  yes if first n.phe columns in row 3 are all blank
    if(all(empty)) {
        map.included <- TRUE
        map <- asnumericwithdec(unlist(data[3,-(1:n.phe)]), dec=dec)
        if(any(is.na(map))) {
            temp <- unique(unlist(data[3,-(1:n.phe)])[is.na(map)])
            stop("There are missing marker positions.\n",
                 "   In particular, we see these value(s): ",
                 paste("\"",paste(temp,collapse="\",\"",sep=""),"\"",collapse=" ",sep=""),
                 " at position(s): ",
                 paste(which(is.na(map)),colapse=",",sep=""),sep="")
        }
        nondatrow <- 3
    }
    else {
        map.included <- FALSE
        map <- rep(0,ncol(data)-n.phe)
        nondatrow <- 2 # last non-data row
    }

    # replace empty cells with NA
    data <- sapply(data,function(a) { a[!is.na(a) & a==""] <- NA; a })

    pheno <- as.data.frame(data[-(1:nondatrow),1:n.phe,drop=FALSE], stringsAsFactors=TRUE)
    colnames(pheno) <- data[1,1:n.phe]

    # pull apart phenotypes, genotypes and map
    mnames <- data[1,-(1:n.phe)]
    if(any(is.na(mnames)))
        stop("There are missing marker names. Check column(s) ",paste(which(is.na(mnames))+1+n.phe,collapse=","),sep="")
    chr <- data[2,-(1:n.phe)]
    if(any(is.na(chr)))
        stop("There are missing chromosome IDs. Check column(s) ",paste(which(is.na(chr))+1+n.phe,collapse=","),sep="")

    if(length(genotypes) > 0) {
        # look for strange entries in the genotype data
        temp <- unique(as.character(data[-(1:nondatrow),-(1:n.phe),drop=FALSE]))

        temp <- temp[!is.na(temp)]
        wh <- !(temp %in% genotypes)
        if(any(wh)) {
            warn <- "The following unexpected genotype codes were treated as missing.\n    "
            ge <- paste("|", paste(temp[wh],collapse="|"),"|",sep="")
            warn <- paste(warn,ge,"\n",sep="")
            warning(warn)
        }

        # convert genotype data
        allgeno <- matrix(match(data[-(1:nondatrow),-(1:n.phe)],genotypes),
                          ncol=ncol(data)-n.phe)
    }
    else
        allgeno <- matrix(as.numeric(data[-(1:nondatrow),-(1:n.phe)]),
                          ncol=ncol(data)-n.phe)

    oldpheno <- pheno

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
        stop("Number of individuals in genotypes and phenotypes do not match.")
    }
    if(any(n.mar1 != n.mar2)) {
        cat(n.mar1,n.mar2,"\n")
        stop("Numbers of markers in genotypes and marker names files do not match.")
    }

    # print some information about the amount of data read
    cat(" --Read the following data:\n")
    cat("\t",n.ind1," individuals\n")
    cat("\t",sum(n.mar1)," markers\n")
    cat("\t",n.phe," phenotypes\n")

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

# end of read.cross.csv.R
