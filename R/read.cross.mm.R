######################################################################
#
# read.cross.mm.R
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
# Contains: read.cross.mm, read.maps.mm
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.mm: read data from an experimental cross in mapmaker
#                format.
#
# We need two files: a "raw" file containing the genotype and
# phenotype data and a "map" file containing the chromosomes
# assignments and (optionally) map positions.
#
# The map file contains two or three columns, separated by white
# space, with the chromosome number, marker name (with markers in
# order along the chromosomes) and (optionally) the map position.
#
######################################################################

read.cross.mm <-
    function(dir,rawfile,mapfile,estimate.map=TRUE)
{
    # create file names
    if(missing(mapfile)) stop("Missing mapfile.")
    if(missing(rawfile)) stop("Missing rawfile.")
    if(!missing(dir)  && dir != "") {
        mapfile <- file.path(dir, mapfile)
        rawfile <- file.path(dir, rawfile)
    }

    # count lines in rawfile
    n.lines <- length(scan(rawfile, what=character(), skip=0, nlines=0,
                           blank.lines.skip=FALSE,quiet=TRUE,sep="\n"))

    # read map file
    map <- read.table(mapfile,header=FALSE,colClasses="character",blank.lines.skip=FALSE)
    fixmap <- TRUE
    if(ncol(map) == 1)
        stop("Map file should contain the markers' chromosome IDs.")

    if(ncol(map) > 3) { # special maps format
        maps <- read.maps.mm(mapfile)
        chr <- rep(names(maps),sapply(maps,length))
        markers <- unlist(lapply(maps,names))
        includes.pos <- TRUE
        fixmap <- FALSE
    }

    if(fixmap) { # my map format: 2 or 3 column table
        # remove any rows lacking a chromosome ID
        o <- (1:nrow(map))[map[,1]==""]
        if(length(o) > 0) map <- map[-o,]

        # remove any leading *'s from the marker names
        g <- grep("^\\*",map[,2])
        if(length(g) > 0)
            map[g,2] <- substr(map[g,2],2,nchar(map[g,2]))
    }

    # begin reading/parsing the genotype data
    cur.mar <- 0
    cur.phe <- 0
    NEW.symb <- c("1","2","3","4","5","0")
    OLD.symb <- c("A","H","B","D","C","-")

    flag <- 0
    #  rawdata <- scan(rawfile,what=character(),sep="\n",
    #                  blank.lines.skip=TRUE,quiet=TRUE)
    #
    for(i in 1:n.lines) {
        a <- scan(rawfile,what=character(),skip=i-1,nlines=1,
                  blank.lines.skip=TRUE,quiet=TRUE)

        if(length(a) == 0) next
        if(length(grep("#", a[1])) != 0) next

        if(flag == 0) {
            flag <- 1
            if(!is.na(match("intercross", a))) type <- "f2"
            else if(!is.na(match("backcross", a)) || !is.na(match("bc", a))) type <- "bc"
            else
                stop("File indicates invalid cross type: ", a[length(a)], ".")
        }
        else if(flag == 1) {
            flag <- 2
            n.ind <- as.numeric(a[1])
            n.mar <- as.numeric(a[2])
            n.phe <- as.numeric(a[3])
            cat(" --Read the following data:\n")
            cat("\tType of cross:         ", type, "\n")
            cat("\tNumber of individuals: ", n.ind, "\n")
            cat("\tNumber of markers:     ", n.mar, "\n")
            cat("\tNumber of phenotypes:  ", n.phe, "\n")

            # if there's a set of "symbols" for non-standard symbols in
            #     the file, use them.
            if(length(a) > 3 && ("symbols" %in% a)) {
                o <- match("symbols",a)
                b <- a[-(1:o)]
                infile.symb <- substring(b,1,1)
                std.symb <- substring(b,3,3)

                wh <- rep(0,length(std.symb))
                fixed <- rep(0,length(OLD.symb))
                for(j in 1:length(std.symb))
                    if(std.symb[j] %in% OLD.symb)
                        wh[j] <- match(std.symb[j],OLD.symb)
                for(j in 1:length(std.symb))
                    if(wh[j] != 0) {
                        OLD.symb[wh[j]] <- infile.symb[j]
                        fixed[wh[j]] <- 1
                    }

                temp <- table(OLD.symb)
                if(any(temp>1)) {
                    for(j in names(temp)[temp>1]) {
                        o <- OLD.symb==j & fixed==0
                        if(any(o)) OLD.symb[o] <- paste(OLD.symb[o],"   ")
                    }
                }
            }

            marnames <- rep("", n.mar)
            geno <- matrix(0,ncol=n.mar,nrow=n.ind)
            if(n.phe == 0) {
                pheno <- matrix(1:n.ind,ncol=1)
                phenames <- c("number")
            }
            else {
                pheno <- matrix(0,ncol=n.phe,nrow=n.ind)
                phenames <- rep("", n.phe)
            }

        }
        else {
            if(substring(a[1],1,1) == "*") {
                cur.mar <- cur.mar+1
                cur.row <- 1

                if(cur.mar > n.mar) { # now reading phenotypes
                    cur.phe <- cur.phe+1
                    if(cur.phe > n.phe) next
                    phenames[cur.phe] <- substring(a[1],2)
                    if(length(a) > 1) {
                        p <- a[-1]
                        p[p=="-"] <- NA
                        n <- length(p)
                        oldna <- is.na(p)
                        numerp <- suppressWarnings(as.numeric(p))
                        newna <- is.na(numerp)
                        wh <- !oldna & newna
                        if(any(wh)) {
                            droppedasmissing <- unique(p[wh])
                            if(length(droppedasmissing) > 1) {
                                themessage <- paste("The values", paste("\"", droppedasmissing, "\"", sep="", collapse=" "))
                                themessage <- paste(themessage, " for phenotype \"", phenames[cur.phe], "\" were", sep="")
                            }
                            else {
                                themessage <- paste("The value \"", droppedasmissing, "\" ", sep="")
                                themessage <- paste(themessage, " for phenotype \"", phenames[cur.phe], "\" was", sep="")
                            }
                            themessage <- paste(themessage, "interpreted as missing.")
                            warning(themessage)
                        }

                        pheno[cur.row+(0:(n-1)),cur.phe] <- numerp
                    }
                    else n <- 0 ## ?
                    cur.row <- cur.row + n
                }

                else { # reading genotypes
                    marnames[cur.mar] <- substring(a[1],2)
                    if(length(a) > 1) {
                        g <- paste(a[-1],collapse="")
                        h <- g <- unlist(strsplit(g,""))
                        for(j in seq(along=NEW.symb)) {
                            if(any(h==OLD.symb[j]))
                                g[h==OLD.symb[j]] <- NEW.symb[j]
                        }

                        n <- length(g)

                        geno[cur.row+(0:(n-1)),cur.mar] <- as.numeric(g)
                    }
                    else n <- 0
                    cur.row <- cur.row + n
                }

            }
            else { # continuation lines
                if(cur.mar > n.mar) { # now reading phenotypes
                    a[a=="-"] <- NA
                    n <- length(a)
                    pheno[cur.row+(0:(n-1)),cur.phe] <- as.numeric(a)
                    cur.row <- cur.row + n
                }
                else {
                    g <- paste(a,collapse="")
                    h <- g <- unlist(strsplit(g,""))
                    for(j in seq(along=NEW.symb)) {
                        if(any(h==OLD.symb[j]))
                            g[h==OLD.symb[j]] <- NEW.symb[j]
                    }
                    n <- length(g)
                    geno[cur.row+(0:(n-1)),cur.mar] <- as.numeric(g)
                    cur.row <- cur.row + n
                }
            } # end continuation line
        } # end non-intro line
    }
    dimnames(pheno) <- list(NULL, phenames)
    # done reading the raw file

    if(fixmap) { # my map format: 2 or 3 column table
        # parse map file
        if(ncol(map) == 3) {
            includes.pos <- TRUE
            # make positions numeric
            pos <- as.numeric(map[,3])
        }
        else includes.pos <- FALSE

        chr <- as.character(map[,1])
        markers <- map[,2]

        # reorder markers?
        if(all(chr %in% c(1:999,"X","x"))) { # 1...19 + X
            tempchr <- chr
            tempchr[chr=="X" | chr=="x"] <- 1000
            tempchr <- as.numeric(tempchr)
            if(includes.pos) neworder <- order(tempchr, pos)
            else neworder <- order(tempchr)
        } else {
            # prevent reordering of chromosomes
            tempchr <- factor(chr, levels=unique(chr))

            if(includes.pos) neworder <- order(tempchr, pos)
            else neworder <- order(tempchr)
        }
        chr <- chr[neworder]
        if(includes.pos) pos <- pos[neworder]
        markers <- markers[neworder]
    }

    Geno <- vector("list",length(unique(chr)))
    names(Geno) <- unique(chr)

    for(i in unique(chr)) {
        mar <- markers[chr == i]

        if(fixmap) { # my map format: 2 or 3 column table
            # create map
            if(includes.pos) {
                map <- pos[chr == i]

                # reorder markers?
                if(any(diff(map)<0)) {
                    o <- order(map)
                    map <- map[o]
                    mar <- mar[o]
                }
            }
            else map <- seq(0,by=5,length=length(mar))
            names(map) <- mar
        }
        else map <- maps[[i]]

        # pull out genotype data
        o <- match(mar,marnames)
        if(any(is.na(o)))
            stop("Cannot find markers in genotype data: ",
                 paste(mar[is.na(o)],collapse=" "), ".",sep="")

        if(length(o)==1) data <- matrix(geno[,o],ncol=1)
        else data <- geno[,o]
        # add marker names to data
        colnames(data) <- mar
        # changes 0's to NA's
        data[!is.na(data) & data==0] <- NA

        Geno[[i]] <- list(data=data,map=map)
        if(i=="X" || i=="x") class(Geno[[i]]) <- "X"
        else class(Geno[[i]]) <- "A"
    }

    cross <- list(geno=Geno,pheno=pheno)
    class(cross) <- c(type,"cross")

    if(estimate.map && !includes.pos) estmap <- TRUE
    else estmap <- FALSE

    cross$pheno <- as.data.frame(cross$pheno, stringsAsFactors=TRUE)

    # return cross + indicator of whether to run est.map
    list(cross,estmap)
}

######################################################################
#
# read.maps.mm: Read genetic map for a special Mapmaker format
# Written by Brian S Yandell; modified by Karl W Broman
#
######################################################################
read.maps.mm <-
    function( mapsfile )
{
    if (missing(mapsfile)) stop("Missing mapsfile.")

    ## find where everything is
    f <- scan(mapsfile, what = "", blank.lines.skip = FALSE, sep = "\n",
              quiet = TRUE)
    start <- pmatch( paste( "*", c("OrderInfo","Classes","Chromosomes",
                                   "Assignments and Placements" ), ":", sep = "" ), f )

    ## marker names
    f <- scan( mapsfile, what = c("",rep(0,9)), skip = start[1],
              nlines = start[2] - start[1] - 1,
              blank.lines.skip = FALSE, quiet = TRUE)
    markers <- substring( f[ seq( 1, length( f ), by = 10 ) ], 2 )

    ## distances
    f <- scan( mapsfile, what = "", skip = start[3],
              nlines = start[4] - start[3] - 1,
              blank.lines.skip = FALSE, quiet = TRUE)
    chr <- grep( "^\\*", f)
    chrom <- substring( f[chr], 2 )
    nmark <- as.integer( f[ 1 + chr ] )
    chr <- c( chr[-1], 1 + length( f ))
    lo <- chr - 2 * nmark + 2
    hi <- chr - nmark
    map <- list()
    imark <- c( 0, cumsum( nmark ))
    for( i in seq( along = chrom )) {
        tmp <- cumsum( c(0,imf.h(as.numeric( f[ lo[i]:hi[i] ] ))))
        names( tmp ) <- markers[ imark[i] + seq( nmark[i] ) ]
        map[[ chrom[i] ]] <- tmp
    }
    map
}

# end of read.cross.mm.R
