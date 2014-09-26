######################################################################
#
# read.cross.qtx.R
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
# Contains: read.cross.qtx
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.qtx
#
# read data in Map Manager QTX format
#
######################################################################

read.cross.qtx <-
    function(dir, file, estimate.map=TRUE)
{
    if(!missing(dir) && dir != "") {
        file <- file.path(dir, file)
    }


    # This is a revised version of match which gives *all* matches
    # of x within the table
    mymatch <-
        function(x, table)
        {
            if(length(x) > 1) x <- x[1] # ignore any but the first element of x

            if(!any(x==table)) return(NA)
            seq(along=table)[x==table]
        }

    # read file into a big vector, each item one line
    cat(" --Read the following data:\n")
    x <- scan(file,what=character(0),sep="\n",quiet=TRUE)

    genoabbrev <- unlist(strsplit(x[9],""))
    if(length(genoabbrev) < 8)  # just in case, fill out to 8 chars
        genoabbrev <- c(genoabbrev,rep("H",8-length(genoabbrev)))
    myabbrev <- c(0,1,3,2,5,4,2,2)
    ugeno <- NULL

    # individuals
    ind.beg <- match("{pgy", x) # there should be just one
    ind.end <- match("}pgy", x)
    n.ind <- as.numeric(x[ind.beg+1])
    ind <- x[(ind.beg+2):(ind.end-1)]
    if(length(ind) != n.ind)
        stop("Problem with individual IDs ({pgy}).")
    cat("\t", n.ind, " individuals\n", sep="")

    # determine if individuals can be viewed as numbers
    g <- grep("^[0-9\\.]+$", ind)
    if(length(g) == n.ind)
        ind <- as.numeric(as.character(ind))

    # phenotypes
    phe.beg <- mymatch("{trt",x)
    phe.end <- mymatch("}trt",x)
    pheno <- NULL
    if(!is.na(phe.beg[1])) { # at least one phenotype
        pheno <- vector("list",length(phe.beg))
        names(pheno) <- paste(phe.beg)
        for(i in 1:length(phe.beg)) {
            z <- x[phe.beg[i]:phe.end[i]]
            names(pheno)[i] <- z[2]
            vals.beg <- match("{tvl", z)+1 # there should be just one match
            vals.end <- match("}tvl", z)-1

            # "X" or "x" is a missing phenotype
            temp <- unlist(strsplit(z[vals.beg[1]:vals.end[1]]," "))
            temp[temp=="X" | temp=="x"] <- NA
            pheno[[i]] <- as.numeric(temp)
        }
        pheno <- cbind(as.data.frame(pheno, stringsAsFactors=TRUE),ind=ind)
        cat("\t", length(pheno), " phenotypes\n",sep="")
    }
    else {
        pheno <- data.frame(ind=ind, stringsAsFactors=TRUE)
        cat("\t", 0, "  phenotypes\n",sep="")
    }

    # chromosomes
    chr.beg <- mymatch("{chx",x)
    chr.end <- mymatch("}chx",x)

    if(is.na(chr.beg[1])) # no genotype data
        stop("There appears to be no genotype data!")
    geno <- vector("list", length(chr.beg))
    names(geno) <- paste(chr.beg)
    has.loci <- rep(TRUE,length(chr.beg))
    map.offset <- rep(0,length(chr.beg))

    cat("\t", length(chr.beg), "  chromosomes\n",sep="")

    for(i in 1:length(chr.beg)) {
        z <- x[chr.beg[i]:chr.end[i]]
        names(geno)[i] <- z[2]
        map.offset <- as.numeric(z[5])

        # loci
        loc.beg <- mymatch("{lox",z)
        loc.end <- mymatch("}lox",z)
        if(all(is.na(loc.beg))) {
            has.loci[i] <- FALSE
            next
        }
        data <- matrix(ncol=length(loc.beg),nrow=n.ind)
        loctype <- rep(NA,length(loc.beg)) ####
        colnames(data) <- paste(loc.beg)
        has.geno <- rep(TRUE,length(loc.beg))
        for(j in 1:length(loc.beg)) {
            zz <- z[loc.beg[j]:loc.end[j]]
            colnames(data)[j] <- zz[2]
            loctype[j] <- zz[5] ####
            geno.beg <- match("{sdp",zz)+1 # should be just one match
            geno.end <- match("}sdp",zz)-1
            if(all(is.na(geno.beg))) { # no genotype data
                has.geno[j] <- FALSE
                next
            }
            dat <- unlist(strsplit(paste(zz[geno.beg[1]:geno.end[1]],collapse=""),""))

            data[,j] <- myabbrev[match(dat,genoabbrev)]
        } # end loop over loci

        # check that all loci have the same code
        if(all(loctype == loctype[1]))
            loctype <- loctype[1]
        # 0 = unknown
        # 1 = backcross codominant maternal unique
        # 2 = backcross codominant paternal unique
        # 3 = backcross maternal dominant
        # 4 = backcross paternal dominant
        # 5 = f2 codominant
        # 6 = f2 maternal dominant
        # 7 = f2 paternal dominant
        # 8 = doubled haploid
        # 9 = selfed RI
        # 10 = sib-mated RI
        # 11 = advanced backcross codominant maternal unique
        # 12 = advanced backcross codominant paternal udnique
        # 13 = advanced backcross maternal dominant
        # 14 = advanced backcross paternal dominant
        # 15 = AIL codominant
        # 16 = AIL maternal dominant
        # 17 = AIL paternal dominant
        # 18 = radiation hybrid data
        # 19 = radiation hybrid data
        # 20 = selfed RIX
        # 21 = sib-mated RIX

        # replace 0's with NA's
        data[!is.na(data) & data==0] <- NA

        # remove columns with no data
        data <- data[,has.geno,drop=FALSE]

        # temporary map
        map <- seq(0,length=ncol(data),by=5)+map.offset
        names(map) <- colnames(data)
        geno[[i]] <- list(data=data,map=map)
        if(length(grep("[Xx]", names(geno)[i]))>0) # X chromosome
            class(geno[[i]]) <- "X"
        else class(geno[[i]]) <- "A"
    } # end loop over chromosomes

    # unique genotypes
    for(i in 1:length(geno)) {
        ugeno <- unique(c(ugeno,unique(geno[[i]]$data)))
        ugeno <- ugeno[!is.na(ugeno)]
    }

    if(length(ugeno)==2) { # backcross
        # Fix if coded as A:B rather than A:H (RI lines)
        if(all(ugeno==1 | ugeno==3)) {
            for(i in 1:length(geno))
                geno[[i]]$data[geno[[i]]$data == 3] <- 2
        }
        # Fix if coded as H:B rather than A:H (other backcross)
        else if(all(ugeno==2 | ugeno==3)) {
            for(i in 1:length(geno))
                geno[[i]]$data[geno[[i]]$data == 3] <- 1
        }

        type <- "bc"
        for(i in 1:length(geno))
            geno[[i]]$data[geno[[i]]$data > 2] <- 1
    }
    else type <- "f2"

    totmar <- sum(sapply(geno,function(a) ncol(a$data)))
    cat("\t", totmar, "  total markers\n",sep="")

    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c(type,"cross")

    if(estimate.map) estmap <- TRUE
    else estmap <- FALSE

    # return cross + indicator of whether to run est.map
    list(cross,estmap)
}

# end of read.cross.qtx.R
