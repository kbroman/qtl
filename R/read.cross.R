######################################################################
#
# read.cross.R
#
# copyright (c) 2000-2014, Karl W Broman
# last modified June, 2014
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
# Contains: read.cross, fixXgeno.bc, fixXgeno.f2, asnumericwithdec
#           [See read.cross.csv.R, read.cross.mm.R, read.cross.qtx.R,
#            qtlcart_io.R, read.cross.gary.R, and read.cross.karl.R
#            for the format-specific functions.]
#
######################################################################

######################################################################
#
# read.cross: read data from an experimental cross
#
######################################################################

read.cross <-
    function(format=c("csv", "csvr", "csvs", "csvsr", "mm", "qtx",
             "qtlcart", "gary", "karl", "mapqtl", "tidy"),
             dir="", file, genfile, mapfile, phefile, chridfile, mnamesfile, pnamesfile,
             na.strings=c("-","NA"), genotypes=c("A","H","B","D","C"),
             alleles=c("A","B"), estimate.map=TRUE, convertXdata=TRUE,
             error.prob=0.0001, map.function=c("haldane", "kosambi", "c-f", "morgan"),
             BC.gen = 0, F.gen = 0, crosstype, ...)
{
    if(format == "csvrs") {
        format <- "csvsr"
        warning("Assuming you mean 'csvsr' rather than 'csvrs'.\n")
    }
    format <- match.arg(format)

    if(format=="csv" || format=="csvr") { # comma-delimited format
        cross <- read.cross.csv(dir, file, na.strings, genotypes,
                                estimate.map, rotate=(format=="csvr"),
                                ...)
    }
    else if(format=="csvs" || format=="csvsr") { # comma-delimited format

        # allow easier input of filenames into function arguments
        if(missing(phefile) && !missing(file) && !missing(genfile)) {
            # read.cross("format", "dir", "genfile", "phefile")
            phefile <- genfile
            genfile <- file
        }
        else if(missing(genfile) && !missing(file) && !missing(phefile)) {
            # read.cross("format", "dir", "genfile", phefile="phefile")
            genfile <- file
        }

        cross <- read.cross.csvs(dir, genfile, phefile, na.strings, genotypes,
                                 estimate.map, rotate=(format=="csvsr"),
                                 ...)
    }
    else if(format=="qtx") { # Mapmanager QTX format
        cross <- read.cross.qtx(dir,file,estimate.map)
    }
    else if(format=="qtlcart") { # QTL Cartographer format
        # if missing mapfile but genfile is specified,
        #     use genfile as the map file.
        if(missing(mapfile) && !missing(genfile))
            mapfile <- genfile

        cross <- read.cross.qtlcart(dir, file, mapfile)
    }
    else if(format=="karl") { # karl's format
        # if missing file names, use standard ones
        if(missing(genfile)) genfile <- "gen.txt"
        if(missing(mapfile)) mapfile <- "map.txt"
        if(missing(phefile)) phefile <- "phe.txt"

        cross <- read.cross.karl(dir,genfile,mapfile,phefile)
    }
    else if(format=="mm") { # mapmaker format
        # if missing mapfile but genfile is specified,
        #     use genfile as the map file.
        if(missing(mapfile) && !missing(genfile))
            mapfile <- genfile

        cross <- read.cross.mm(dir,file,mapfile,estimate.map)
    }
    else if(format=="gary") { # gary's format
        # if missing file names, use the standard ones
        if(missing(genfile)) genfile <- "geno.dat"
        if(missing(mnamesfile)) mnamesfile <- "mnames.txt"
        if(missing(chridfile)) chridfile <- "chrid.dat"
        if(missing(phefile)) phefile <- "pheno.dat"
        if(missing(pnamesfile)) pnamesfile <- "pnames.txt"
        if(missing(mapfile)) mapfile <- "markerpos.txt"

        cross <- read.cross.gary(dir,genfile,mnamesfile,chridfile,
                                 phefile,pnamesfile,mapfile,estimate.map,na.strings)
    }
    else if(format == "mapqtl") { # MapQTL format (same as JoinMap)
        cross <- read.cross.mq(dir=dir, locfile=genfile, mapfile=mapfile,
                               quafile=phefile)
    }
    else if (format == "tidy") { # tidy format
        if(!missing(file) && !missing(genfile) && !missing(mapfile) && missing(phefile)) {
            phefile <- mapfile
            mapfile <- genfile
            genfile <- file
        }

        if(missing(genfile)) genfile <- "gen.csv"
        if(missing(phefile)) phefile <- "phe.csv"
        if(missing(mapfile)) mapfile <- "map.csv"

        cross <- read.cross.tidy(dir=dir,
                                 genfile=genfile, phefile=phefile, mapfile=mapfile,
                                 na.strings=na.strings, genotypes=genotypes)
    }

    estimate.map <- cross[[2]]
    cross <- cross[[1]]

    # if chr names all start with "chr" or "Chr", remove that part
    chrnam <- names(cross$geno)
    if(all(regexpr("^[Cc][Hh][Rr]",chrnam)>0)){
        chrnam <- substr(chrnam,4,nchar(chrnam))
        if(all(regexpr("^[Oo][Mm][Oo][Ss][Oo][Mm][Ee]",chrnam)>0))
            chrnam <- substr(chrnam,8,nchar(chrnam))
    }
    # if chr named "x" make it "X"
    if(sum(chrnam=="x")>0) chrnam[chrnam=="x"] <- "X"
    names(cross$geno) <- chrnam
    # make sure the class of chromosomes named "X" is "X"
    for(i in 1:length(cross$geno))
        if(names(cross$geno)[i] == "X")
            class(cross$geno[[i]]) <- "X"

    # Fix up the X chromosome data for a backcross or intercross
    chrtype <- sapply(cross$geno,class)
    if(any(chrtype=="X") && convertXdata) {
        if(class(cross)[1]=="bc")
            cross <- fixXgeno.bc(cross)
        if(class(cross)[1]=="f2") {
            if(missing(alleles)) alleles <- c("A","B")
            cross <- fixXgeno.f2(cross, alleles)
        }
    }

    ## Pass through read.cross.bcsft for BCsFt (convert if appropriate).
    cross <- read.cross.bcsft(cross = cross, BC.gen = BC.gen, F.gen = F.gen, ...)

    # re-estimate map?
    if(estimate.map) {
        cat(" --Estimating genetic map\n")
        map.function <- match.arg(map.function)
        newmap <- est.map(cross, error.prob=error.prob, map.function=map.function)
        cross <- replace.map(cross, newmap)
    }

    # store genotype data as integers
    for(i in 1:nchr(cross))
        storage.mode(cross$geno[[i]]$data) <- "integer"

    # check alleles
    if(class(cross)[1] != "4way") {
        if(length(alleles) > 2) {
            warning("length of arg alleles should be 2")
            alleles <- alleles[1:2]
        }
        if(length(alleles) < 2)
            stop("length of arg alleles should be 2")
    }
    else { # 4-way cross
        if(missing(alleles))
            alleles <- c("A","B","C","D")

        if(length(alleles) > 4) {
            warning("length of arg alleles should be 4 for a 4-way cross")
            alleles <- alleles[1:4]
        }
        if(length(alleles) < 4)
            stop("length of arg alleles should be 4 for a 4-way cross")
    }
    if(any(nchar(alleles)) != 1) {
        warning("Each item in arg alleles should be a single character")
        alleles <- substr(alleles, 1, 1)
    }

    attr(cross, "alleles") <- alleles

    type <- class(cross)[1]
    if(!missing(crosstype)) {
        if(crosstype=="risib") cross <- convert2risib(cross)
        else if(crosstype=="riself") cross <- convert2riself(cross)
        else class(cross)[1] <- crosstype
    }

    # run checks
    summary(cross)

    cat(" --Cross type:", class(cross)[1], "\n")

    cross
}




##############################
# fixXgeno.bc: fix up the X chromosome genotype data for backcross
##############################
fixXgeno.bc <-
    function(cross)
{
    omitX <- FALSE

    # pull out X chr genotype data
    chrtype <- sapply(cross$geno,class)
    xchr <- which(chrtype=="X")
    Xgeno <- cross$geno[[xchr]]$data

    # find "sex" and "pgm" in the phenotype data
    sexpgm <- getsex(cross)

    if(!is.null(sexpgm$sex)) {     # "sex" is provided
        malegeno <- Xgeno[sexpgm$sex==1,]
        if(any(!is.na(malegeno) & malegeno==2)) {
            n.omit <- sum(!is.na(malegeno) & malegeno==2)
            warning(" --Omitting ", n.omit, " male heterozygote genotypes on the X chromosome.")
            malegeno[!is.na(malegeno) & malegeno==2] <- NA
        }
        malegeno[!is.na(malegeno) & malegeno==3] <- 2

        femalegeno <- Xgeno[sexpgm$sex==0,]
        if(any(!is.na(femalegeno) & femalegeno==3)) {
            n.omit <- sum(!is.na(femalegeno) & femalegeno==3)
            warning(" --Omitting ", n.omit, " BB genotypes from females on the X chromosome.")
            femalegeno[!is.na(femalegeno) & femalegeno==3] <- NA
        }

        Xgeno[sexpgm$sex==1,] <- malegeno
        Xgeno[sexpgm$sex==0,] <- femalegeno

    }
    else {
        # "sex" not provided

        if(all(is.na(Xgeno) | Xgeno==1 | Xgeno==3)) { # look like all males
            warning(" --Assuming that all individuals are male.\n")
            Xgeno[!is.na(Xgeno) & Xgeno==3] <- 2
            cross$pheno$sex <- factor(rep("m",nind(cross)),levels=c("f","m"))
        }
        else if(all(is.na(Xgeno) | Xgeno==1 | Xgeno==2)) { # look like females A:H
            warning(" --Assuming that all individuals are female.\n")
            cross$pheno$sex <- factor(rep("f",nind(cross)),levels=c("f","m"))
        }
        else { # have some of each of the three genotypes
            warning(" --Can't figure out the X chromosome genotypes.\n   You need to provide phenotypes \"sex\"\n   See the help file for read.cross() for details.\n   Omitting the X chr for now.\n  ")
            omitX <- TRUE
        }
    }

    if(!omitX) {
        wh <- !is.na(Xgeno) & Xgeno!=1 & Xgeno!=2
        if(any(wh)) {
            Xgeno[wh] <- NA
            n.omit <- sum(wh)
            warning(" --Omitted ", n.omit, " additional X chr genotype(s).")
        }

        cross$geno[[xchr]]$data <- Xgeno
    }
    else cross <- subset(cross,chr= -xchr) # <- omit the X chr completely

    cross
}


##############################
# fixXgeno.f2: fix up the X chromosome genotype data for intercross
##############################
fixXgeno.f2 <-
    function(cross, alleles)
{
    omitX <- FALSE

    # pull out X chr genotype data
    chrtype <- sapply(cross$geno,class)
    xchr <- which(chrtype=="X")
    Xgeno <- cross$geno[[xchr]]$data

    # find "sex" and "pgm" in the phenotype data
    sexpgm <- getsex(cross)

    AA <- paste(rep(alleles[1], 2), collapse="")
    AB <- paste(alleles, collapse="")
    BB <- paste(rep(alleles[2], 2), collapse="")
    cross0 <- paste("(", alleles[1], "x", alleles[2], ")x(",
                    alleles[1], "x", alleles[2], ")", sep="")
    cross1 <- paste("(", alleles[2], "x", alleles[1], ")x(",
                    alleles[2], "x", alleles[1], ")", sep="")

    if(!is.null(sexpgm$sex) && !is.null(sexpgm$pgm)) {
        # both "sex" and "pgm" are provided

        if(any(sexpgm$sex == 1)) { # there are males
            malegeno <- Xgeno[sexpgm$sex==1,]
            if(any(!is.na(malegeno) & malegeno==2)) {
                n.omit <- sum(!is.na(malegeno) & malegeno==2)
                warning(" --Omitting ", n.omit, " male heterozygote genotypes on the X chromosome.")
                malegeno[!is.na(malegeno) & malegeno==2] <- NA
            }
            malegeno[!is.na(malegeno) & malegeno==3] <- 2
            Xgeno[sexpgm$sex==1,] <- malegeno
        }

        if(any(sexpgm$sex==0)) { # there are females
            femalegeno0 <- Xgeno[sexpgm$sex==0 & sexpgm$pgm==0,]
            femalegeno1 <- Xgeno[sexpgm$sex==0 & sexpgm$pgm==1,]

            if((any(!is.na(femalegeno0) & femalegeno0==3) ||
                any(!is.na(femalegeno1) & femalegeno1==1)) &&
               !any(!is.na(femalegeno0) & femalegeno0==1) &&
               !any(!is.na(femalegeno1) & femalegeno1==3)) {
                # appear to switched the "pgm" values
                warning(" --The 0/1 values for \"pgm\" appear to be switched; switching back.")
                sexpgm$pgm[sexpgm$pgm==1] <- 2
                sexpgm$pgm[sexpgm$pgm==0] <- 1
                sexpgm$pgm[sexpgm$pgm==2] <- 0
                cross$pheno$pgm[cross$pheno$pgm==1] <- 2
                cross$pheno$pgm[cross$pheno$pgm==0] <- 1
                cross$pheno$pgm[cross$pheno$pgm==2] <- 0

                temp <- femalegeno0
                femalegeno0 <- femalegeno1
                femalegeno1 <- temp
            }
            if(any(!is.na(femalegeno0) & femalegeno0==3)) {
                n.omit <- sum(!is.na(femalegeno0) & femalegeno0==3)
                warning(" --Omitting ", n.omit, " ", BB, " genotypes from females from cross ",
                        cross0, " on the X chr.\n")
                femalegeno0[!is.na(femalegeno0) & femalegeno0==3] <- NA
            }
            if(any(!is.na(femalegeno1) & femalegeno1==1)) {
                n.omit <- sum(!is.na(femalegeno1) & femalegeno1==1)
                warning(" --Omitting ", n.omit, " ", AA, " genotypes from females from cross ",
                        cross1, " on the X chr.\n")
                femalegeno1[!is.na(femalegeno1) & femalegeno1==1] <- NA
            }
            femalegeno1[!is.na(femalegeno1) & femalegeno1==3] <- 1
            Xgeno[sexpgm$sex==0 & sexpgm$pgm==0,] <- femalegeno0
            Xgeno[sexpgm$sex==0 & sexpgm$pgm==1,] <- femalegeno1
        }

    }

    else if(!is.null(sexpgm$sex) && is.null(sexpgm$pgm)) {
        # "sex" is provided but not "pgm"

        if(any(sexpgm$sex == 1)) { # there are males
            malegeno <- Xgeno[sexpgm$sex==1,]
            if(any(!is.na(malegeno) & malegeno==2)) {
                n.omit <- sum(!is.na(malegeno) & malegeno==2)
                warning(" --Omitting ", n.omit, " male heterozygote genotypes on the X chromosome.")
                malegeno[!is.na(malegeno) & malegeno==2] <- NA
            }
            malegeno[!is.na(malegeno) & malegeno==3] <- 2
            Xgeno[sexpgm$sex==1,] <- malegeno
        }

        if(any(sexpgm$sex==0)) { # there are females
            femalegeno <- Xgeno[sexpgm$sex==0,]

            if(any(!is.na(femalegeno) & femalegeno==3) &
               !any(!is.na(femalegeno) & femalegeno==1)) { # looks like (BxA)x(BxA)
                cross$pheno$pgm <- rep(1,nind(cross))
                femalegeno[!is.na(femalegeno) & femalegeno==3] <- 1
            }
            else if(any(!is.na(femalegeno) & femalegeno==1) &
                    !any(!is.na(femalegeno) & femalegeno==3)) { # looks like (AxB)x(AxB)
                cross$pheno$pgm <- rep(0,nind(cross))
            }
            else { # we have some 1's and some 3's
                warning(" --There appear to be some individuals of each cross direction, but \"pgm\" is not provided.\n   Check the X chr genotype data and include a \"pgm\" column in the phenotype data.\n   \"pgm\" was inferred (probably poorly).\n   ")

                cross$pheno$pgm <- rep(0,nind(cross))

                # females with no 3's -> assumed to be from (AxB)x(AxB)
                # females with both 3's and 1's -> assumed to be from (AxB)x(AxB); 3's tossed
                wh.have3 <- apply(femalegeno, 1, function(a) any(!is.na(a) & a==3))
                cross$pheno$pgm[sexpgm$sex==0][wh.have3] <- 1

                temp <- femalegeno[wh.have3,]
                temp[!is.na(temp) & temp==1] <- NA
                temp[!is.na(temp) & temp==3] <- 1
                femalegeno[wh.have3,] <- temp

            }

            Xgeno[sexpgm$sex==0,] <- femalegeno
        }
    }
    else if(is.null(sexpgm$sex) && !is.null(sexpgm$pgm)) {
        # "pgm" is provided but not "sex"

        if(all(is.na(Xgeno) | Xgeno==1 | Xgeno==3)) { # look like all males
            cross$pheno$sex <- factor(rep("m",nind(cross)),levels=c("f","m"))
            Xgeno[!is.na(Xgeno) & Xgeno==3] <- 2
        }
        else { # assume all females
            cross$pheno$sex <- factor(rep("f",nind(cross)),levels=c("f","m"))
            Xgeno.pgm0 <- Xgeno[sexpgm$pgm==0,]
            Xgeno.pgm1 <- Xgeno[sexpgm$pgm==1,]
            if(all(is.na(Xgeno.pgm0) | Xgeno.pgm0==2 | Xgeno.pgm0==3) &&
               all(is.na(Xgeno.pgm1) | Xgeno.pgm1==1 | Xgeno.pgm1==2)) {
                cross$pheno$pgm <- 1 - sexpgm$pgm
                temp <- Xgeno.pgm0
                Xgeno.pgm0 <- Xgeno.pgm1
                Xgeno.pgm1 <- temp
            }
            Xgeno.pgm1[!is.na(Xgeno.pgm1) & Xgeno.pgm1==1] <- NA
            Xgeno.pgm1[!is.na(Xgeno.pgm1) & Xgeno.pgm1==3] <- 1
            Xgeno.pgm0[!is.na(Xgeno.pgm0) & Xgeno.pgm0==3] <- NA
            Xgeno[sexpgm$pgm==0,] <- Xgeno.pgm0
            Xgeno[sexpgm$pgm==1,] <- Xgeno.pgm1
        }
    }


    else {
        # Neither "sex" and "pgm" provided

        if(all(is.na(Xgeno) | Xgeno==1 | Xgeno==3)) { # look like all males
            warning(" --Assuming that all individuals are male.\n")
            Xgeno[!is.na(Xgeno) & Xgeno==3] <- 2
            cross$pheno$sex <- factor(rep("m",nind(cross)),levels=c("f","m"))
            cross$pheno$pgm <- rep(0,nind(cross))
        }
        else if(all(is.na(Xgeno) | Xgeno==2 | Xgeno==3)) { # look like females H:B
            warning(" --Assuming that all individuals are female.\n")
            Xgeno[!is.na(Xgeno) & Xgeno==3] <- 1
            cross$pheno$sex <- factor(rep("f",nind(cross)),levels=c("f","m"))
            cross$pheno$pgm <- rep(1,nind(cross))
        }
        else if(all(is.na(Xgeno) | Xgeno==2 | Xgeno==1)) { # looks like females A:H
            warning(" --Assuming that all individuals are female.\n")
            cross$pheno$sex <- factor(rep("f",nind(cross)),levels=c("f","m"))
            cross$pheno$pgm <- rep(0,nind(cross))
        }
        else { # have some of each of the three genotypes
            warning(" --Can't figure out the X chromosome genotypes.\n   You need to provide phenotypes \"sex\" and/or \"pgm\"\n   See the help file for read.cross() for details.\n   Omitting the X chr for now.\n  ")

            omitX <- TRUE

        }
    }


    if(!omitX) {
        wh <- !is.na(Xgeno) & Xgeno!=1 & Xgeno!=2
        if(any(wh)) {
            Xgeno[wh] <- NA
            n.omit <- sum(wh)
            warning(" --Omitted ", n.omit, " additional X chr genotype(s).")
        }

        cross$geno[[xchr]]$data <- Xgeno
    }
    else cross <- subset(cross,chr= -xchr) # <- omit the X chr completely

    cross
}


######################################################################
# convert character to numeric, using dec as the decimal point
######################################################################
asnumericwithdec <-
    function(x, dec=".")
{
    if(dec!=".") x <- gsub(paste0("\\", dec), ".", x)
    as.numeric(x)
}

# Fix up phenotypes
sw2numeric <-
    function(x, dec) {
        x[x == ""] <- NA
        wh1 <- is.na(x)
        n <- sum(!is.na(x))
        y <- suppressWarnings(asnumericwithdec(as.character(x), dec))
        wh2 <- is.na(y)
        m <- sum(!is.na(y))
        if(n==m || (n-m) < 2 || (n-m) < n*0.05) {
            if(sum(!wh1 & wh2) > 0) {
                u <- unique(as.character(x[!wh1 & wh2]))
                if(length(u) > 1) {
                    themessage <- paste("The phenotype values", paste("\"", u, "\"", sep="", collapse=" "))
                    themessage <- paste(themessage, " were", sep="")
                }
                else {
                    themessage <- paste("The phenotype value \"", u, "\" ", sep="")
                    themessage <- paste(themessage, " was", sep="")
                }
                themessage <- paste(themessage, "interpreted as missing.")
                warning(themessage)

            }
            return(y)
        }
        else return(x)
    }

# end of read.cross.R
