######################################################################
#
# write.cross.R
#
# copyright (c) 2001-2014, Karl W Broman and Hao Wu
# last modified June, 2014
# first written Feb, 2001
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
# Contains: write.cross, write.cross.mm, write.cross.csv,
#           write.cross.gary, write.cross.tidy, fixX4write
#           [See qtlcart_io.R for write.cross.qtlcart]
#           [write.cross.qtab in write.cross.qtab.R]
#
######################################################################


######################################################################
#
# write.cross: Wrapper for the other write.cross functions
#
######################################################################

write.cross <-
    function(cross, format=c("csv", "csvr", "csvs", "csvsr", "mm", "qtlcart",
                    "gary", "qtab", "mapqtl", "tidy"),
             filestem="data", chr, digits=NULL, descr)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    format <- match.arg(format)
    if(!missing(chr)) cross <- subset(cross,chr=chr)

    # revise X data
    chrtype <- sapply(cross$geno,class)
    crosstype <- class(cross)[1]
    if((crosstype=="bc" || crosstype=="f2") &&
       any(chrtype=="X")) {
        sexpgm <- getsex(cross)
        sex <- sexpgm$sex
        pgm <- sexpgm$pgm
        for(i in which(chrtype=="X"))
            cross$geno[[i]]$data <- fixX4write(cross$geno[[i]]$data,sex,pgm,crosstype)
    }

    if(crosstype == "bcsft") # convert BCsFt to intercross for writing
        class(cross)[1] <- "f2"

    if(format=="csv") write.cross.csv(cross,filestem,digits,FALSE,FALSE)
    else if(format=="csvr") write.cross.csv(cross,filestem,digits,TRUE,FALSE)
    else if(format=="csvs") write.cross.csv(cross,filestem,digits,FALSE,TRUE)
    else if(format=="csvsr") write.cross.csv(cross,filestem,digits,TRUE,TRUE)
    else if(format=="mm") write.cross.mm(cross,filestem,digits)
    else if(format=="qtlcart") write.cross.qtlcart(cross, filestem)
    else if(format=="gary") write.cross.gary(cross, digits)
    else if(format=="tidy") write.cross.tidy(cross, filestem, digits)
    else if(format=="qtab") {
        if(missing(descr)) descr <- paste(deparse(substitute(cross)), "from R/qtl")
        write.cross.qtab(cross, filestem, descr, verbose=FALSE)
    } else if(format == "mapqtl")
        write.cross.mq(cross, filestem, digits)
}




######################################################################
#
# write.cross.mm: Write data for an experimental cross in Mapmaker
#                 format
#
#           creates two files: "raw" file with geno & pheno data
#                              "prep" file with map information
#
######################################################################

write.cross.mm <-
    function(cross, filestem="data", digits=NULL)
{
    n.ind <- nind(cross)
    tot.mar <- totmar(cross)
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    type <- class(cross)[1]
    if(type=="riself" || type=="risib" || type=="dh" || type=="haploid") type <- "bc"
    if(type != "f2" && type != "bc")
        stop("write.cross.mm only works for intercross, backcross, doubled haploid and RI data.")

    # write genotype and phenotype data
    file <- paste(filestem, ".raw", sep="")

    # write experiment type
    if(type == "f2")
        write("data type f2 intercross", file, append=FALSE)
    else
        write("data type f2 backcross", file, append=FALSE)

    # write numbers of progeny, markers and phenotypes
    write(paste(n.ind, tot.mar, n.phe), file, append=TRUE)

    # max length of marker name
    mlmn <- max(nchar(unlist(lapply(cross$geno,function(a) colnames(a$data)))))+1

    # write marker data
    for(i in 1:n.chr) {
        for(j in 1:n.mar[i]) {
            mn <- paste("*", colnames(cross$geno[[i]]$data)[j], sep="")
            if(nchar(mn) < mlmn)
                mn <- paste(mn,paste(rep(" ", mlmn-nchar(mn)),collapse=""),sep="")
            g <- cross$geno[[i]]$data[,j]

            x <- rep("", n.ind)
            x[is.na(g)] <- "-"
            x[!is.na(g) & g==1] <- "A"
            x[!is.na(g) & g==2] <- "H"
            if(type == "f2") {
                x[!is.na(g) & g==3] <- "B"
                x[!is.na(g) & g==4] <- "D"
                x[!is.na(g) & g==5] <- "C"
            }

            if(n.ind < 60)
                write(paste(mn, paste(x,collapse="")), file, append=TRUE)
            else {
                lo <- seq(1,n.ind-1,by=60)
                hi <- c(lo[-1]-1,n.ind)
                for(k in seq(along=lo)) {
                    if(k==1) write(paste(mn,paste(x[lo[k]:hi[k]],collapse="")),file,append=TRUE)
                    else write(paste(paste(rep(" ", mlmn),collapse=""),
                                     paste(x[lo[k]:hi[k]],collapse="")),file,append=TRUE)
                }
            }

        }
    } # end writing marker data

    # max length of phenotype name
    mlpn <- max(nchar(colnames(cross$pheno)))+1

    # write phenotypes
    for(i in 1:n.phe) {
        pn <- paste("*",colnames(cross$pheno)[i],sep="")
        if(nchar(pn) < mlpn)
            pn <- paste(pn, paste(rep(" ", mlpn-nchar(pn)),collapse=""),sep="")

        if(!is.factor(cross$pheno[,i])) {
            if(is.null(digits)) x <- as.character(cross$pheno[,i])
            else x <- as.character(round(cross$pheno[,i],digits))
        }
        else
            x <- as.character(cross$pheno[,i])
        x[is.na(x)] <- "-"

        if(n.ind < 10)
            write(paste(pn, paste(x,collapse="")), file, append=TRUE)
        else {
            lo <- seq(1,n.ind-1,by=10)
            hi <- c(lo[-1]-1,n.ind)
            for(k in seq(along=lo)) {
                if(k==1) write(paste(pn,paste(x[lo[k]:hi[k]],collapse=" ")),file,append=TRUE)
                else write(paste(paste(rep(" ", mlpn),collapse=""),
                                 paste(x[lo[k]:hi[k]],collapse=" ")),file,append=TRUE)
            }
        }
    }


    # make "prep" file with map information
    file <- paste(filestem, ".prep", sep="")

    for(i in 1:n.chr) {
        cname <- paste("chr", names(cross$geno)[i], sep="")
        line <- paste("make chromosome", cname)
        if(i==1) write(line, file, append=FALSE)
        else write(line, file, append=TRUE)

        mn <- names(cross$geno[[i]]$map)
        #    dis <- round(diff(cross$geno[[i]]$map),2)
        #    dis <- paste("=", dis, sep="")
        #    write(paste(paste("sequence", mn[1]), paste(dis,mn[-1],collapse=" ")),
        #          file, append=TRUE)
        write(paste(paste("sequence", mn[1]), paste(mn[-1],collapse=" ")),
              file, append=TRUE)

        write(paste("anchor", cname), file, append=TRUE)
        write(paste("framework", cname), file, append=TRUE)
    }

}

######################################################################
#
# write.cross.csv: Write data for an experimental cross in
#                  comma-delimited format (the same format as is read
#                  by read.cross.csv)
#
######################################################################

write.cross.csv <-
    function(cross, filestem="data", digits=NULL, rotate=FALSE, split=FALSE)
{
    type <- class(cross)[1]
    if(type != "f2" && type != "bc" && type != "riself" && type != "risib" && type != "dh" && type != "haploid")
        stop("write.cross.csv only works for intercross, backcross, RI, doubled haploid, and haploid data.")

    if(!split)
        file <- paste(filestem, ".csv", sep="")
    else {
        genfile <- paste(filestem, "_gen.csv", sep="")
        phefile <- paste(filestem, "_phe.csv", sep="")
    }

    if(split) { # split files; need individual IDs
        id <- getid(cross)
        if(is.null(id)) {
            cross$pheno$id <- 1:nind(cross)
            id <- getid(cross)
        }
        id.col <- which(colnames(cross$pheno)==attr(id,"phenam"))
    }

    n.ind <- nind(cross)
    tot.mar <- totmar(cross)
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    geno <- matrix(ncol=tot.mar,nrow=n.ind)

    # allele codes to use
    if("alleles" %in% names(attributes(cross))) {
        alle <- attr(cross, "alleles")
        alleles <- c(paste(alle[1],alle[1],sep=""),
                     paste(alle[1],alle[2],sep=""),
                     paste(alle[2],alle[2],sep=""),
                     paste("not ", alle[2],alle[2],sep=""),
                     paste("not ", alle[1],alle[1],sep=""))
    }
    else alleles <- c("A","H","B","D","C")

    if(type=="dh" || type=="riself" || type=="risib") alleles[2:3] <- alleles[3:2]
    else if(type=="haploid") alleles <- alle

    firstmar <- 1
    for(i in 1:n.chr) {
        # replace allele numbers with
        geno[,firstmar:(firstmar+n.mar[i]-1)] <-
            alleles[match(cross$geno[[i]]$data,1:5)]
        firstmar <- firstmar + n.mar[i]
    }
    if(any(is.na(geno))) geno[is.na(geno)] <- "-"
    pheno <- cross$pheno
    for(i in 1:nphe(cross)) {
        if(is.factor(pheno[,i])) pheno[,i] <- as.character(pheno[,i])
        else if(is.numeric(pheno[,i])) {
            if(!is.null(digits)) pheno[,i] <- round(pheno[,i], digits)
            pheno[,i] <- as.character(pheno[,i])
        }
    }
    pheno <- matrix(unlist(pheno), nrow=n.ind)

    if(any(is.na(pheno))) pheno[is.na(pheno)] <- "-"
    thedata <- cbind(pheno,geno)
    colnames(thedata) <- c(colnames(cross$pheno),
                           unlist(lapply(cross$geno, function(a) colnames(a$data))))
    chr <- rep(names(cross$geno),n.mar)
    pos <- unlist(lapply(cross$geno,function(a) a$map))
    chr <- c(rep("",n.phe),chr)
    if(!is.null(digits))
        pos <- c(rep("",n.phe),as.character(round(pos,digits)))
    else
        pos <- c(rep("",n.phe),as.character(pos))

    # put it all together
    thenames <- colnames(thedata)
    thedata <- matrix(as.character(thedata), ncol=ncol(thedata))
    thedata <- rbind(thenames, chr, pos, thedata)

    if(!split) {
        if(!rotate)
            write.table(thedata, file, quote=FALSE, sep=",",
                        row.names=FALSE, col.names=FALSE)
        else
            write.table(t(thedata), file, quote=FALSE, sep=",",
                        row.names=FALSE, col.names=FALSE)
    }
    else { # split files: one for phenotypes and one for genotypes
        n.phe <- nphe(cross)
        phe <- thedata[-(2:3),1:n.phe]
        gen <- cbind(thedata[,id.col], thedata[,-(1:n.phe)])

        if(!rotate) {
            write.table(gen, genfile, quote=FALSE, sep=",",
                        row.names=FALSE, col.names=FALSE)
            write.table(phe, phefile, quote=FALSE, sep=",",
                        row.names=FALSE, col.names=FALSE)
        }
        else {
            write.table(t(gen), genfile, quote=FALSE, sep=",",
                        row.names=FALSE, col.names=FALSE)
            write.table(t(phe), phefile, quote=FALSE, sep=",",
                        row.names=FALSE, col.names=FALSE)
        }

    }


}


######################################################################
#
# write.cross.gary: Write data for an experimental cross in
# Gary's format. There will be 6 output files, they are:
#    chrid.dat - chromosome ids
#    markerpos.txt - marker position
#    mnames.txt - marker names
#    geno.data - genotypes
#    pheno.data - phenotypes
#    pnames.txt - phenotype names
#
######################################################################

write.cross.gary <-
    function(cross, digits=NULL)
{
    # local variables
    n.ind <- nind(cross)
    tot.mar <- totmar(cross)
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    # chromosome ids
    chrid <- NULL
    for(i in 1:n.chr) {
        # the name for this chromosome
        chrname <- names(cross$geno[i])
        # convert to number (why?)
        #    if(chrname=="X") chrname <- 20
        #    else chrname <- as.numeric(chrname)
        chrid <- c(chrid, rep(chrname, n.mar[i]))
    }
    write.table(chrid, file="chrid.dat", quote=FALSE, row.names=FALSE,
                col.names=FALSE)

    # marker position file
    markpos <- NULL
    for(i in 1:n.chr)
        markpos <- c(markpos, cross$geno[[i]]$map)
    write.table(markpos, file="markerpos.txt", quote=FALSE, sep="\t",
                row.names=TRUE, col.names=FALSE)

    # marker names
    mnames <- names(markpos)
    write.table(mnames, file="mnames.txt", quote=FALSE, row.names=FALSE,
                col.names=FALSE)

    # genotype
    geno <- NULL
    for(i in 1:n.chr)
        geno <- cbind(geno, cross$geno[[i]]$data)
    # note that gary's format codes genotype from 0
    # and 9 is for NA
    geno <- geno - 1 # note NA will still be NA
    write.table(geno, file="geno.dat", quote=FALSE, row.names=FALSE,
                col.names=FALSE, sep="\t", na="9")

    # phenotype
    pheno <- cross$pheno
    for(i in 1:nphe(cross)) {
        if(is.factor(pheno[,i])) pheno[,i] <- as.character(pheno[,i])
        else if(is.numeric(pheno[,i])) {
            if(is.null(digits)) pheno[,i] <- as.character(pheno[,i])
            else pheno[,i] <- as.character(round(pheno[,i], digits))
        }
    }
    pheno <- matrix(unlist(pheno), nrow=nrow(pheno))

    write.table(pheno, file="pheno.dat", quote=FALSE, row.names=FALSE,
                col.names=FALSE, sep="\t", na="-999")
    # phenotype names
    write.table(names(cross$pheno), file="pnames.txt", quote=FALSE,
                row.names=FALSE, col.names=FALSE, sep="\t", na="-999")

}


######################################################################
#
# write.cross.tidy: Write data for an experimental cross in tidy
# format. There will be 3 output files, they are:
#    geno.csv
#    pheno.csv
#    map.csv
#
######################################################################

write.cross.tidy <-
    function(cross, filestem="data", digits=NULL)
{
    genfile <- paste(filestem, "_gen.csv", sep="")
    phefile <- paste(filestem, "_phe.csv", sep="")
    mapfile <- paste(filestem, "_map.csv", sep = "")

    type <- class(cross)[1]

    id <- getid(cross)

    if(is.null(id)) {
        cross$pheno$id <- 1:nind(cross)
        id <- getid(cross)
    }
    id.col <- which(colnames(cross$pheno) == attr(id,"phenam"))

    # allele codes to use
    if("alleles" %in% names(attributes(cross))) {
        alle <- attr(cross, "alleles")
        alleles <- c(paste(alle[1],alle[1],sep=""),
                     paste(alle[1],alle[2],sep=""),
                     paste(alle[2],alle[2],sep=""),
                     paste("not ", alle[2],alle[2],sep=""),
                     paste("not ", alle[1],alle[1],sep=""))
    }
    else alleles <- c("A","H","B","D","C")

    if (type %in% c("dh", "riself", "risib")) alleles[2:3] <- alleles[3:2]
    else if (type == "haploid") alleles <- alle

    geno <- pull.geno(cross)
    geno <- matrix(alleles[geno], ncol = nind(cross), byrow = TRUE,
                   dimnames = list(markernames(cross), make.names(id)))

    if(any(is.na(geno))) geno[is.na(geno)] <- "-"

    pheno <- cross$pheno[-id.col]
    for(i in 1:ncol(pheno)) {
        if(is.factor(pheno[,i])) pheno[,i] <- as.character(pheno[,i])
        else if(is.numeric(pheno[,i])) {
            if(!is.null(digits)) pheno[,i] <- round(pheno[,i], digits)
            pheno[,i] <- as.character(pheno[,i])
        }
    }

    pheno <- matrix(unlist(pheno), nrow=nind(cross),
                    dimnames = list(make.names(id), phenames(cross)[-id.col]))
    pheno <- t(pheno)

    if(any(is.na(pheno))) pheno[is.na(pheno)] <- "-"

    map <- pull.map(cross, as.table = TRUE)

    if(!is.null(digits))
        map$pos <- as.character(round(map$pos, digits))
    else
        map$pos <- as.character(map$pos)

    write.table(geno,  genfile, quote = FALSE, sep = ",", col.names = NA)
    write.table(pheno, phefile, quote = FALSE, sep = ",", col.names = NA)
    write.table(map,   mapfile, quote = FALSE, sep = ",", col.names = NA)
}


######################################################################
# fixX4write
######################################################################

fixX4write <-
    function(geno,sex,pgm,crosstype)
{
    # males
    if(!is.null(sex) & any(sex==1)) {
        temp <- geno[sex==1,,drop=FALSE]
        temp[temp==2] <- 3
        geno[sex==1,] <- temp
    }

    if(crosstype == "f2") {
        # females
        if(!is.null(pgm)) {
            if(!is.null(sex)) {
                if(any(pgm==1 & sex==0)) {
                    temp <- geno[sex==0 & pgm==1,,drop=FALSE]
                    temp[temp==1] <- 3
                    geno[sex==0 & pgm==1,] <- temp
                }
            }
            else { # assume all females
                if(any(pgm==1)) {
                    temp <- geno[pgm==1,,drop=FALSE]
                    temp[temp==1] <- 3
                    geno[pgm==1,] <- temp
                }
            }
        }

    }

    geno
}

# end of write.cross.R
