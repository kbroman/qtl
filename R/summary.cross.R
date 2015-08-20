######################################################################
#
# summary.cross.R
#
# copyright (c) 2001-2014, Karl W Broman
# last modified Dec, 2014
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
# Contains: summary.cross, print.summary.cross, nind, nchr, nmar,
#           totmar, nphe, nmissing, ntyped, print.cross, chrlen
#
######################################################################

summary.cross <-
    function(object,...)
{
    if(!any(class(object) == "cross"))
        stop("Input should have class \"cross\".")

    n.ind <- nind(object)
    tot.mar <- totmar(object)
    n.phe <- nphe(object)
    n.chr <- nchr(object)
    n.mar <- nmar(object)
    type <- class(object)[1]

    if(!(type %in% c("f2", "bc", "4way", "riself", "risib", "dh", "haploid",
                     "ri4self", "ri4sib", "ri8self", "ri8selfIRIP1", "ri8sib", "bcsft", "bgmagic16")))
        stop("Cross type ", type, " is not supported.")

    # combine genotype data into one big matrix
    Geno <- pull.geno(object)

    # A and X-specific genotypes?
    chrtype <- sapply(object$geno, class)
    if(any(chrtype=="A")) {
        GenoA <- pull.geno(object, chr=chrnames(object)[chrtype=="A"])
    } else {
        GenoA <- NULL
    }
    if(any(chrtype=="X") && type %in% c("f2", "bc", "bcsft")) {
        GenoX <- pull.geno(object, chr=chrnames(object)[chrtype=="X"])
        GenoX <- reviseXdata(type, "full", getsex(object), geno=GenoX, cross.attr=attributes(object))
    } else{
        GenoX <- NULL
    }

    # proportion of missing genotype data
    missing.gen <- mean(is.na(Geno))

    # Get cross scheme for BCsFt.
    if(type == "bcsft") {
        cross.scheme <- attr(object, "scheme")
        is.bcs <- (cross.scheme[2] == 0)
    }
    else {
        cross.scheme <- rep(0,2)
        is.bcs <- FALSE
    }

    # table of genotype values
    typings <- typingsA <- typingsX <- NULL
    if(type %in% c("f2", "bcsft") & !is.bcs) {
        if(is.null(GenoX)) {
            typings <- table(factor(Geno[!is.na(Geno)], levels=1:5))
            temp <- getgenonames("f2", "A", cross.attr=attributes(object))
            names(typings) <- c(temp, paste("not", temp[c(3,1)]))
        } else {
            if(!is.null(GenoA)) {
                typingsA <- table(factor(GenoA[!is.na(GenoA)], levels=1:5))
                temp <- getgenonames("f2", "A", cross.attr=attributes(object))
                names(typingsA) <- c(temp, paste("not", temp[c(3,1)]))
            }
            if(!is.null(GenoX)) {
                gnames <- getgenonames("f2", "X", "full", getsex(object), cross.attr=attributes(object))
                typingsX <- table(factor(GenoX[!is.na(GenoX)], levels=1:length(gnames)))
                names(typingsX) <- gnames
            }
        }
    }
    else if(type %in% c("bc", "riself", "risib", "dh", "haploid", "bcsft")) {
        if(is.null(GenoX)) {
            typings <- table(factor(Geno[!is.na(Geno)], levels=1:2))
            temp <- getgenonames(type, "A", cross.attr=attributes(object))
            names(typings) <- temp
        } else {
            if(!is.null(GenoA)) {
                typingsA <- table(factor(GenoA[!is.na(GenoA)], levels=1:2))
                temp <- getgenonames(type, "A", cross.attr=attributes(object))
                names(typingsA) <- temp
            }
            if(!is.null(GenoX)) {
                gnames <- getgenonames(type, "X", "full", getsex(object), cross.attr=attributes(object))
                typingsX <- table(factor(GenoX[!is.na(GenoX)], levels=1:length(gnames)))
                names(typingsX) <- gnames
            }
        }
    }
    else if(type=="4way") {
        typings <- table(factor(Geno[!is.na(Geno)], levels=1:14))

        temp <- getgenonames("4way", "A", cross.attr=attributes(object))
        names(typings) <- c(temp,
                            paste(temp[c(1,3)], collapse="/"),
                            paste(temp[c(2,4)], collapse="/"),
                            paste(temp[c(1,2)], collapse="/"),
                            paste(temp[c(3,4)], collapse="/"),
                            paste(temp[c(1,4)], collapse="/"),
                            paste(temp[c(2,3)], collapse="/"),
                            paste("not", temp[1]),
                            paste("not", temp[2]),
                            paste("not", temp[3]),
                            paste("not", temp[4]))
    }
    else typings <- table(factor(Geno[!is.na(Geno)]))

    # turn into fractions
    if(!is.null(typings)) typings <- typings/sum(typings)
    if(!is.null(typingsA)) typingsA <- typingsA/sum(typingsA)
    if(!is.null(typingsX)) typingsX <- typingsX/sum(typingsX)


    # amount of missing phenotype data
    if(ncol(object$pheno) <= 30)
        missing.phe <- as.numeric(cbind(apply(object$pheno,2,function(a) mean(is.na(a)))))
    else
        missing.phe <- mean(as.numeric(is.na(object$pheno)))

    # check that, in the case of a "4way" cross, the genetic
    #     maps are matrices with 2 rows, and that for other crosses,
    #     the genetic maps are numeric vectors
    if(type=="4way") {
        if(any(!sapply(object$geno, function(a) (is.matrix(a$map) && nrow(a$map)==2))))
            warning("The genetic maps should all be matrices with two rows.")
    }
    else {
        if(any(sapply(object$geno, function(a) is.matrix(a$map))))
            warning("The genetic maps should all be numeric vectors rather than matrices.")
    }

    # check that object$geno[[i]]$data has colnames and that they match
    #     the names in object$geno[[i]]$map
    jitterwarning <- NULL
    for(i in 1:n.chr) {
        nam1 <- colnames(object$geno[[i]]$data)
        map <- object$geno[[i]]$map
        if(is.matrix(map)) nam2 <- colnames(map)
        else nam2 <- names(map)
        chr <- names(object$geno)[[i]]
        if(is.null(nam1)) {
            warn <- paste("The data matrix for chr", chr,
                          "lacks column names")
            warning(warn)
        }
        if(is.null(nam2)) {
            warn <- paste("The genetic map for chr", chr,
                          "lacks column names")
            warning(warn)
        }
        if(any(nam1 != nam2))
            stop("Marker names in the data matrix and genetic map\n",
                 "for chr ", chr, " do not match.")

        if((is.matrix(map) && (any(diff(map[1,])<0) || any(diff(map[2,])<0))) ||
           (!is.matrix(map) && any(diff(map)<0)))
            stop("Markers out of order on chr ", chr)


        # check that no two markers are on top of each other
        if(is.matrix(map)) { # sex-specific maps
            n <- ncol(map)
            if(n > 1) {
                d1 <- diff(map[1,])
                d2 <- diff(map[2,])
                if(any(d1 < 1e-14 & d2 < 1e-14)) {
                    if (is.null(jitterwarning)) jitterwarning<-list()
                    jitterwarning[[names(object$geno)[i]]]<-which(d1 < 1e-14 & d2 < 1e-14)
                }
            }
        }
        else {
            n <- length(map)
            if(n > 1) {
                d <- diff(map)
                if(any(d < 1e-14)) {
                    if (is.null(jitterwarning)) jitterwarning<-list()
                    jitterwarning[[names(object$geno)[i]]]<-which(d < 1e-14)
                }
            }
        }

    }

    if (!is.null(jitterwarning))
        warning("Some markers at the same position on chr ",
                paste(names(jitterwarning),collapse=",",sep=""),"; use jittermap().")

    if(!is.data.frame(object$pheno))
        warning("Phenotypes should be a data.frame.")

    if(is.null(colnames(object$pheno)))
        stop("Phenotype data needs column names")

    x <- table(colnames(object$pheno))

    if(any(x > 1))
        warning("Some phenotypes have the same name:\n",
                paste(names(x)[x>1], collapse="  "))

    # check genotype data
    if(type %in% c("bc", "riself", "risib", "dh", "haploid") | (type == "bcsft" & is.bcs)) {
        # Invalid genotypes?
        if(any(!is.na(Geno) & Geno != 1 & Geno != 2)) {
            u <- unique(as.numeric(Geno))
            u <- sort(u[!is.na(u)])
            warn <- paste("Invalid genotypes.",
                          "\n    Observed genotypes:",
                          paste(u, collapse=" "))
            warning(warn)
        }

        # Missing genotype category on autosomes?
        if(sum(!is.na(Geno) & Geno==2) == 0 ||
           sum(!is.na(Geno) & Geno==1) == 0) {
            warn <- paste("Strange genotype pattern on chr ", chr, ".", sep="")
            warning(warn)
        }
    }
    else if(type %in% c("f2","bcsft") & !is.bcs) {
        # invalid genotypes
        if(any(!is.na(Geno) & Geno!=1 & Geno!=2 & Geno!=3 &
               Geno!=4 & Geno!=5)) {
            u <- unique(as.numeric(Geno))
            u <- sort(u[!is.na(u)])
            warn <- paste("Invalid genotypes on chr", chr, ".",
                          "\n    Observed genotypes:",
                          paste(u, collapse=" "))
            warning(warn)
        }

        # X chromosome
        for(i in 1:n.chr) {
            if(class(object$geno[[i]]) == "X") {
                dat <- object$geno[[i]]$data
                if(any(!is.na(dat) & dat!=1 & dat!=2)) {
                    u <- unique(as.numeric(dat))
                    u <- sort(u[!is.na(u)])
                    warn <- paste("Invalid genotypes on X chromosome:",
                                  "\n    Observed genotypes:",
                                  paste(u, collapse=" "))
                    warning(warn)
                }
            }
        }

        # Missing genotype category on autosomes?
        dat <- NULL; flag <- 0
        for(i in 1:n.chr) {
            if(class(object$geno[[i]]) != "X") {
                dat <- cbind(dat,object$geno[[i]]$data)
                flag <- 1
            }
        }
        if(flag && (sum(!is.na(dat) & dat==2) == 0 ||
                    sum(!is.na(dat) & dat==1) == 0 ||
                    sum(!is.na(dat) & dat==3) == 0))
            warning("Strange genotype pattern.")
    }
    else if(type=="4way") {
        # Invalid genotypes?
        if(any(!is.na(Geno) & (Geno != round(Geno) | Geno < 1 | Geno > 14))) {
            u <- unique(as.numeric(Geno))
            u <- sort(u[!is.na(u)])
            warn <- paste("Invalid genotypes.",
                          "\n    Observed genotypes:",
                          paste(u, collapse=" "))
            warning(warn)
        }
    }
    else if(type %in% c("ri4sib", "ri4self", "ri8sib", "ri8self", "ri8selfIRIP1", "bgmagic16")) {
        if(type=="bgmagic16") n.str <- 16
        else n.str <- as.numeric(substr(type, 3, 3))
        if(any(!is.na(Geno) & (Geno != round(Geno) | Geno < 1 | Geno > 2^n.str-1))) {
            u <- unique(as.numeric(Geno))
            u <- sort(u[!is.na(u)])
            warn <- paste("Invalid genotypes.",
                          "\n    Observed genotypes:",
                          paste(u, collapse=" "))
            warning(warn)
        }
    }

    # Look for duplicate marker names
    mnames <- NULL
    for(i in 1:nchr(object))
        mnames <- c(mnames,colnames(object$geno[[i]]$data))

    o <- table(mnames)
    if(any(o > 1))
        warning("Duplicate markers [", paste(names(o)[o>1], collapse=", "), "]")

    # make sure the genotype data are matrices rather than data frames
    if(any(sapply(object$geno, function(a) is.data.frame(a$data))))
        warning("The $data objects should be simple matrices, not data frames.")

    # make sure each chromosome has class "A" or "X"
    chr.class <- sapply(object$geno, class)
    if(!all(chr.class == "A" | chr.class == "X"))
        warning("Each chromosome should have class \"A\" or \"X\".")
    chr.nam <- names(object$geno)
    if(is.null(chr.nam)) {
        warning("The chromosome names are missing.")
        chr.nam <- as.character(1:length(chr.class))
    }

    if(type != "riself" && any(chr.class=="A" & (chr.nam=="X" | chr.nam=="x"))) {
        wh <- which(chr.nam=="X" | chr.nam=="x")
        warning("Chromosome \"", chr.nam[wh], "\" has class \"A\" but probably ",
                "should have class \"X\".")
    }

    if(length(chr.nam) > length(unique(chr.nam))) {
        tab <- table(chr.nam)
        dups <- names(tab[tab > 1])
        warning("Duplicate chromosome names: ", paste(dups,sep=", "))
    }
    g <- grep("^-", chr.nam)
    if(length(g) > 0)
        warning("Chromosome names shouldn't start with '-': ", paste(chr.nam[g], sep=", "))

    # if more than one X chromosome, print a warning
    if(sum(chr.class=="X") > 1)
        warning("More than one X chromosome: [",
                paste(chr.nam[chr.class == "X"], collapse=", "), "]")

    if(any(chr.class=="A"))
        autosomes <- chr.nam[chr.class == "A"]
    else autosomes <- NULL
    if(any(chr.class=="X"))
        Xchr <- chr.nam[chr.class == "X"]
    else Xchr <- NULL

    # check individual IDs
    id <- getid(object)
    if(!is.null(id) && length(id) != length(unique(id)))
        warning("The individual IDs are not unique.")

    # check that chromosomes aren't too long
    mapsum <- summaryMap(object)
    if(ncol(mapsum)==7)  # sex-specific map
        maxlen <- max(mapsum[1:(nrow(mapsum)-1),2:3])
    else
        maxlen <- max(mapsum[1:(nrow(mapsum)-1),2])
    if(maxlen > 1000)
        warning(paste("Some chromosomes > 1000 cM in length; there may",
                      "be a problem with the genetic map.\n  (Perhaps it is in basepairs?)"))

    cross.summary <- list(type=type, n.ind = n.ind, n.phe=n.phe,
                          n.chr=n.chr, n.mar=n.mar,
                          missing.gen=missing.gen,
                          typing.freq=typings, typing.freq.A=typingsA, typing.freq.X=typingsX,
                          missing.phe=missing.phe,
                          autosomes=autosomes, Xchr=Xchr, cross.scheme=cross.scheme)
    class(cross.summary) <- "summary.cross"
    cross.summary

}


print.summary.cross <-
    function(x,...)
{
    print.genotypes <- TRUE

    if(x$type=="f2") cat("    F2 intercross\n\n")
    else if(x$type=="bc") cat("    Backcross\n\n")
    else if(x$type=="4way") cat("    4-way cross\n\n")
    else if(x$type=="riself") cat("    RI strains via selfing\n\n")
    else if(x$type=="risib") cat("    RI strains via sib matings\n\n")
    else if(x$type=="dh") cat("    Doubled haploids\n\n")
    else if(x$type=="haploid") cat("    Haploids\n\n")
    else if(x$type %in% c("ri4self", "ri4sib", "ri8self", "ri8selfIRIP1", "ri8sib")) {
        n.str <- substr(x$type, 3, 3)
        if(substr(x$type, 4, min(6, nchar(x$type)))=="sib") crosstype <- "sib-mating"
        else crosstype <- "selfing"
        print.genotypes <- FALSE
        cat("    ", n.str, "-way RIL by ", crosstype, "\n\n", sep="")
    }
    else if(x$type=="bcsft") cat(paste("    BC(", x$cross.scheme[1], ")F(", x$cross.scheme[2], ") cross\n\n", sep = ""))
    else if(x$type %in% c("bgmagic16")) {
        n.str <- 16
        print.genotypes <- FALSE
        cat("    ", n.str, "-way Biogemma MAGIC lines\n\n", sep="")
    }
    else cat("    cross", x$type, "\n\n",sep=" ")

    cat("    No. individuals:   ", x$n.ind,"\n\n")
    cat("    No. phenotypes:    ", x$n.phe,"\n")

    header <- "                       "
    width <- options("width")$width

    cat("    Percent phenotyped:")

    ######################################################################
    # function to print things nicely
    printnicely <-
        function(thetext, header, width, sep=" ")
        {
            nleft <- width - nchar(header)
            nsep <- nchar(sep)
            if(length(thetext) < 2) cat("", thetext, "\n", sep=sep)
            else {
                z <- paste("", thetext[1], sep=sep, collapse=sep)
                for(j in 2:length(thetext)) {
                    if(nchar(z) + nsep + nchar(thetext[j]) > nleft) {
                        cat(z, "\n")
                        nleft <- width
                        z <- paste(header, thetext[j], sep=sep)
                    }
                    else {
                        z <- paste(z, thetext[j], sep=sep)
                    }
                }
                cat(z, "\n")
            }
        }
    ######################################################################
    # function to pre-pad with spaces
    pad_w_spaces <-
        function(x, width, pre=TRUE)
        {
            padding <- vapply(nchar(x), function(n) paste(rep(" ", width-n), collapse=""), "")
            if(pre) return(paste0(padding, x))
            paste0(x, padding)
        }
    ######################################################################

    printnicely(round((1-x$missing.phe)*100,1), header, width)
    cat("\n")

    cat("    No. chromosomes:   ", x$n.chr,"\n")
    if(!is.null(x$autosomes)) {
        cat("        Autosomes:     ")
        printnicely(x$autosomes, header, width)
    }
    if(!is.null(x$Xchr)) {
        cat("        X chr:         ")
        printnicely(x$Xchr, header, width)
    }
    cat("\n")
    cat("    Total markers:     ", sum(x$n.mar), "\n")
    cat("    No. markers:       ")
    printnicely(x$n.mar, header, width)
    cat("    Percent genotyped: ", round((1-x$missing.gen)*100,1), "\n")
    if(print.genotypes) {
        cat("    Genotypes (%):    ")
        header <- "                      "
        if(!is.null(x$typing.freq.X)) {
            # contortions to line things up
            roundedX <- sprintf("%.1f", x$typing.freq.X*100)

            genoX <- paste(names(x$typing.freq.X),roundedX,sep=":")
            if(!is.null(x$typing.freq.A)) {
                roundedA <- sprintf("%.1f", x$typing.freq.A*100)
                genoA <- paste(names(x$typing.freq.A),roundedA,sep=":")

                # line up values
                maxchar <- max(nchar(c(roundedA, roundedX)))
                roundedA <- pad_w_spaces(roundedA, maxchar, pre=FALSE)
                roundedX <- pad_w_spaces(roundedX, maxchar, pre=FALSE)
                genoX <- paste(names(x$typing.freq.X),roundedX,sep=":")
                genoA <- paste(names(x$typing.freq.A),roundedA,sep=":")

                # line things up again
                maxchar <- max(nchar(c(genoA, genoX)))
                genoA <- pad_w_spaces(genoA, maxchar)
                genoX <- pad_w_spaces(genoX, maxchar)

                cat("\n          Autosomes:  ")
                printnicely(genoA, header, width, "  ")
                cat("       X chromosome:  ")
            }
            printnicely(genoX, header, width, "  ")
        } else {
            if(!is.null(x$typing.freq.A) && !is.null(x$typing.freq))
                x$typing.freq <- x$typing.freq.A
            geno <- paste(names(x$typing.freq),sprintf("%.1f", x$typing.freq*100), sep=":")
            printnicely(geno, header, width, "  ")
        }
    }
}



nind <-
    function(object)
{
    if(!any(class(object) == "cross"))
        stop("Input should have class \"cross\".")

    n.ind1 <- nrow(object$pheno)
    n.ind2 <- sapply(object$geno,function(x) nrow(x$data))
    if(any(n.ind2 != n.ind1))
        stop("Different numbers of individuals in genotypes and phenotypes.")
    n.ind1
}

nchr <-
    function(object)
{
    cl <- class(object)
    if(!("cross" %in% cl || "map" %in% cl))
        stop("Input should have class \"cross\" or \"map\".")

    if("map" %in% cl)
        return(length(object))

    length(object$geno)
}

nmar <-
    function(object)
{
    cl <- class(object)
    if(!("cross" %in% cl || "map" %in% cl))
        stop("Input should have class \"cross\" or \"map\".")

    if("map" %in% cl) {
        if(is.matrix(object[[1]]))
            return(sapply(object, function(x) ncol(x)))
        else return(sapply(object, function(x) length(x)))
    }

    if(length(object$geno) == 0)
        stop("There is no genotype data.")

    if(!is.matrix(object$geno[[1]]$map))
        n.mar1 <- sapply(object$geno, function(x) length(x$map))
    else # sex-specific maps
        n.mar1 <- sapply(object$geno, function(x) ncol(x$map))

    n.mar2 <- sapply(object$geno, function(x) ncol(x$data))
    if(any(n.mar1 != n.mar2))
        stop("Different numbers of markers in genotypes and maps.")
    n.mar1
}

totmar <-
    function(object)
{
    cl <- class(object)
    if(!("cross" %in% cl || "map" %in% cl))
        stop("Input should have class \"cross\" or \"map\".")

    if("map" %in% cl) {
        if(is.matrix(object[[1]]))
            return(sum(sapply(object, function(x) ncol(x))))
        else return(sum(sapply(object, function(x) length(x))))
    }

    if(length(object$geno) == 0)
        stop("There is no genotype data.")

    if(!is.matrix(object$geno[[1]]$map))
        totmar1 <- sum(sapply(object$geno, function(x) length(x$map)))
    else # sex-specific maps
        totmar1 <- sum(sapply(object$geno, function(x) ncol(x$map)))
    totmar2 <- sum(sapply(object$geno, function(x) ncol(x$data)))
    if(totmar1 != totmar2)
        stop("Different numbers of markers in genotypes and maps.")
    totmar1
}

nphe <-
    function(object)
{
    if(!any(class(object) == "cross"))
        stop("Input should have class \"cross\".")

    ncol(object$pheno)
}

# count number of missing genotypes for each individual or each marker
nmissing <-
    function(cross,what=c("ind","mar"))
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    what <- match.arg(what)

    if(what=="ind") {
        n.missing <- rep(0,nind(cross))
        for(i in 1:nchr(cross))
            n.missing <- n.missing +
                apply(cross$geno[[i]]$data,1,function(a) sum(is.na(a)))

        # individual IDs
        id <- getid(cross)
        if(!is.null(id)) names(n.missing) <- id
    }
    else {
        n.missing <- NULL
        for(i in 1:nchr(cross))
            n.missing <- c(n.missing,
                           apply(cross$geno[[i]]$data,2,function(a) sum(is.na(a))))
    }

    n.missing
}


# like nmissing, but for the opposite value
ntyped <-
    function(cross, what=c("ind","mar"))
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    what <- match.arg(what)

    if(what=="ind") n <- totmar(cross)
    else n <- nind(cross)

    n - nmissing(cross, what)
}

# "print" method for cross object
#
# to avoid ever printing the entire object, print just a little
#     warning message and then the summary
print.cross <-
    function(x, ...)
{
    cat("  This is an object of class \"cross\".\n")
    cat("  It is too complex to print, so we provide just this summary.\n")
    print(summary(x))
    return(summary(x))
}


# get chromosome lengths
chrlen <-
    function(object)
{
    if(!any(class(object) == "map") && !any(class(object) == "cross"))
        stop("Input should have class \"map\" or \"cross\".")

    if(!any(class(object) == "map"))
        x <- pull.map(object)
    else x <- object

    if(is.matrix(x[[1]]))
        return(sapply(x, apply, 1, function(a) diff(range(a))))

    sapply(x, function(a) diff(range(a)))
}


# end of summary.cross.R
