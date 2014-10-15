#####################################################################
#
# qtlcart_io.R
#
# copyright (c) 2002-2013, Brian S. Yandell
#          [with some modifications by Karl W. Broman and Hao Wu]
# last modified Mar, 2013
# first written Jun, 2002
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
# Contains: read.cross.qtlcart, read.cro.qtlcart, read.map.qtlcart,
#           write.cross.qtlcart
#           [See read.cross.R for the main read.cross function, and
#            write.cross.R for the main write.cross function.]
#
##############################################################################

######################################################################
# read.cross.qtlcart
#
# read QTL cross object in QTL cartographer format
######################################################################
read.cross.qtlcart <-
    function (dir, crofile, mapfile)
{
    if (missing(mapfile)) stop("Missing mapfile.")
    if (missing(crofile)) stop("Missing crofile.")

    if(!missing(dir) && dir != "") {
        mapfile <- file.path(dir, mapfile)
        crofile <- file.path(dir, crofile)
    }
    map <- read.map.qtlcart( mapfile )
    cro <- read.cro.qtlcart( crofile )

    cat(" --Read the following data:\n")
    cat("       Type of cross:         ", cro$cross.class, "\n")
    cat("       Number of individuals: ", nrow( cro$markers ), "\n")
    cat("       Number of markers:     ", ncol( cro$markers ), "\n")
    cat("       Number of phenotypes:  ", ncol( cro$traits ), "\n")

    maplen <- unlist(lapply(map,length))
    markers <- split( as.data.frame( t( cro$markers ), stringsAsFactors=TRUE),
                     ordered( rep(names( maplen ), maplen )))

    Geno <- list()
    for( i in names( map )) {
        name.markers <- names( map[[i]] )
        markers[[i]] <- t( markers[[i]] )
        colnames( markers[[i]] ) <- name.markers
        tmp <- list( data = markers[[i]], map = map[[i]] )

        # determine whether autosomal chromosome or X chromosome
        #     using the chromosome name
        class(tmp) <- ifelse(length(grep("[Xx]", i)), "X", "A")
        Geno[[i]] <- tmp
    }
    cross <- list(geno = Geno, pheno = cro$traits )
    cross$pheno <- as.data.frame(cross$pheno, stringsAsFactors=TRUE)
    class(cross) <- c( cro$cross.class, "cross")
    if(cro$cross.class == "bcsft")
        attr(cross, "scheme") <- cro$cross.scheme

    list(cross,FALSE)
}

######################################################################
# read.map.qtlcart
#
# read QTL Cartographer map file
######################################################################
read.map.qtlcart <-
    function (file)
{
    # only interested in chromosomes, marker IDs and positions
    f <- scan(file, what = "", blank.lines.skip = FALSE, sep = "\n",
              quiet = TRUE)
    ctrl <- seq(f)[substring(f, 1, 1) == "-"]
    getvalue <- function(s, f, ctrl) {
        tmp <- unlist(strsplit(f[ctrl[substring(f[ctrl], 2, 3) ==
                                      s]], " "))
        as.numeric(tmp["" != tmp][2])
    }
    nchrom <- getvalue("c ", f, ctrl)
    nmarkers <- getvalue("i ", f, ctrl)

    # marker positions
    tmp <- range(seq(f)[substring(f, 1, 3) == "-l "])
    s <- strsplit(f[tmp[1]], "")[[1]]
    b <- grep("\\|", s)
    s <- grep("0", s)
    s <- ceiling((s[length(s)] - b - 1)/nchrom)

    position <- scan(file, what=character(), sep="\n",
                     skip = tmp[1] - 1, n = tmp[2], quiet=TRUE)


    tmp <- grep("-b", f)
    if(length(tmp) < 1) stop("Marker names not found in map file\n")
    markers <- scan(file, list(1, 2, ""), skip = tmp[1], nlines = nmarkers,
                    blank.lines.skip = FALSE, quiet = TRUE)

    if(length(tmp) < 2) {
        warning("Chromosome names not found in map file\n")
        chroms <- as.character(1:nchrom)
    }
    else
        chroms <- scan(file, list(1, ""), skip = tmp[2], nlines = nchrom,
                       blank.lines.skip = FALSE, quiet = TRUE)[[2]]

    map <- list()
    n.markers <- table(markers[[1]])

    # make sense of "position"
    position <- strsplit(position, "\\s+")
    pos <- matrix(ncol=nchrom,nrow=max(n.markers))
    for(i in 1:max(n.markers)) {
        z <- position[[i]][-(1:3)]
        pos[i,which(n.markers >= i-1)] <- as.numeric(z)
    }

    for (i in seq(nchrom)) {
        tmp <- cumsum(pos[1:n.markers[i],i])
        names(tmp) <- markers[[3]][i == markers[[1]]]
        map[[chroms[i]]] <- tmp
    }
    map
}

######################################################################
# read.cro.qtlcart
#
# read QTL cartographer CRO file
######################################################################
read.cro.qtlcart <-
    function (file)
{
    # translation from cro to R/qtl (see read.cross)
    # -1  NA      missing data
    #  0  1       AA
    #  1  2       AB
    #  2  3       BB
    # 10  4       AA or AB
    # 12  5       AB or BB
    #
    f <- scan(file, what = "", blank.lines.skip = FALSE, sep = "\n",
              quiet = TRUE)
    ctrl <- seq(f)[substring(f, 1, 1) == "-"]
    s <- strsplit(f[ctrl], " ")
    ns <- character(length(ctrl))
    for (i in seq(ctrl)) {
        ns[i] <- substring(s[[i]][1], 2)
        s[[i]] <- s[[i]]["" != s[[i]]][-1]
    }
    names(s) <- ns
    size <- as.numeric(s$n[1])
    nmarkers <- as.numeric(s$p[1]) - 1
    ntraits <- as.numeric(s$traits[1])

    # cross.scheme (used for bcsft only)
    cross.scheme <- c(0,0)

    # cross type
    fix.bc1 <- fix.ridh <- FALSE # indicator of whether to fix genotypes
    cross <- tolower(s$cross[1])
    if (cross == "ri1" || cross=="riself") {
        cross <- "riself"
        fix.ridh <- TRUE
    }
    else if (cross == "ri2" || cross=="risib") {
        cross <- "risib"
        fix.ridh <- TRUE
    }
    else if (cross == "ri0" || cross=="dh") {
        cross <- "dh" # doubled haploid
        fix.ridh <- TRUE
    }
    else if (cross == "b1" || cross == "b2") {
        fix.bc1 <- cross == "b1"
        cross <- "bc"
    }
    else if (cross == "sf2" || cross == "rf2")
        cross <- "f2"
    else if (cross == "sf3") {
        cross <- "bcsft"
        cross.scheme <- c(0,3)
    }

    if(!cross %in% c("f2","bc","risib","riself","bcsft", "dh"))
        stop("Cross type ", cross, " not supported.")

    notraits <- as.numeric(s$otraits[1])
    skip <- ctrl["s" == ns]
    nlines <- ctrl["e" == ns] - skip - 1
    trait.names <- f[ctrl["Names" == ns][1] + 1:ntraits]
    if(notraits)
        trait.names <- c(trait.names, f[ctrl["Names" == ns][2] + 1:notraits] )
    ns <- strsplit(trait.names, " ")
    for (i in seq(ns)) ns[[i]] <- ns[[i]][length(ns[[i]])]
    trait.names <- unlist(ns)
    # kludge to handle factor phenos
    f <- matrix(scan(file, "", skip = skip, nlines = nlines, na.strings = ".",
                     blank.lines.skip = FALSE, quiet = TRUE), ncol = size)
    traits <- t(f[-(1:(2 + nmarkers)), ])
    traits = as.data.frame(traits, stringsAsFactors=TRUE)
    if (nrow(traits) == 1)
        traits <- as.data.frame(t(traits), stringsAsFactors=TRUE)
    colnames(traits) <- trait.names

    tmp = options(warn=-1)
    for(i in names(traits)){
        tmp1 = as.numeric(as.character(traits[[i]]))
        if(!all(is.na(tmp1))) traits[[i]] = tmp1
    }
    options(tmp)
    f <- t(f[3:(2 + nmarkers), ])

    # here is the translation
    f = array(as.numeric(f),dim(f))

    # omit negative genotypes (treat as missing)
    f[f<0] <- NA

    f[!is.na(f)] <- c(NA, 1:3, rep(NA, 7), 4, NA, 5)[2 + f[!is.na(f)]]
    if (fix.ridh && all(is.na(f) | f == 1 | f == 3))
        f[!is.na(f) & f == 3] <- 2
    if (fix.bc1) {
        f[!is.na(f) & f == 5] <- NA
        f[!is.na(f) & f == 2] <- 1
        f[!is.na(f) & f == 3] <- 2
    }
    list(traits = traits, markers = f, cross.class = cross, cross.scheme = cross.scheme)
}


######################################################################
# write.cross.qtlcart
#
# write a QTL cross object to files in QTL Cartographer format
######################################################################
write.cross.qtlcart <-
    function( cross, filestem="data")
{
    n.ind <- nind(cross)
    tot.mar <- totmar(cross)
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    type <- class(cross)[1]
    if(type=="bc") {
        type <- paste("B", 1 + match(1, names(table(c(pull.geno(cross)))),
                                     nomatch = 0), sep = "")
    }
    else if(type=="f2") type <- "RF2"
    else if(type=="riself") type <- "RI1"
    else if(type=="risib") type <- "RI2"
    else {
        warn <- paste("Cross type", type, "may not work with QTL Cartographer.")
        warning(warn)
    }

    # RIL data: convert genotypes to 1/3; later will be converted to 0/2
    if(type=="RI1" || type=="RI2") {
        for(i in 1:n.chr) {
            g <- cross$geno[[i]]$data
            g[!is.na(g) & g==2] <- 3
            cross$geno[[i]]$data <- g
        }
    }

    # write genotype and phenotype data
    file <- paste(filestem, ".cro", sep="")
    if( file.exists( file )) {
        warning( paste( "previous file", file, "moved to *.bak" ))
        file.rename( file, paste( file, "bak", sep = "." ))
    }
    write("#  123456789 -filetype Rcross.out", file, append=FALSE)

    # write numbers of progeny, markers and phenotypes
    write( paste( "-n   ", n.ind ), file, append=TRUE)
    write( paste( "-p   ", 1 + tot.mar ), file, append=TRUE)
    # write experiment type
    write( paste( "-cross", type ), file, append=TRUE)

    # write numbers of progeny, markers and phenotypes
    write( paste( "-traits   ", n.phe ), file, append=TRUE)
    write( "-Names of the traits...", file, append=TRUE)
    phe <- names( cross$pheno )
    for( i in seq( phe ))
        write( paste( i, phe[i] ), file, append=TRUE)
    write( paste( "-otraits   ", 0 ), file, append=TRUE)

    # write genotype and phenotype data by individual
    write( "-s", file, append=TRUE)
    for( ind in 1:n.ind ) {
        write( paste( ind, 1 ), file, append=TRUE)
        for(i in 1:n.chr) {
            g <- unlist( cross$geno[[i]]$data[ind,] )
            g[ is.na( g ) ] <- 0
            g <- c(-1,0,1,2,10,12)[ 1 + g ]

            #      if( length( g ) <= 40)
            write(paste( "      ", paste( g, collapse = " " )), file, append=TRUE)
            #      else {
            #        lo <- seq( 1, length(g), by=40)
            #        hi <- c( lo[-1]-1, length( g ))
            #        for(k in seq(along=lo)) {
            #          write( paste( "      ", paste( g[lo[k]:hi[k]], collapse = " " )),
            #                file, append=TRUE)
            #        }
            #      }
        } # end writing marker data
        p <- c( cross$pheno[ind,])
        tmp <- format( p )
        tmp[ is.na( p ) ] <- "."
        write( paste( "       ", tmp ), file, append = TRUE )
        # end of writing phenotype data
    }
    write( "-e", file, append = TRUE )
    write( "-q", file, append = TRUE )

    # make "prep" file with map information
    file <- paste(filestem, ".map", sep="")
    if( file.exists( file )) {
        warning( paste( "previous file", file, "moved to *.bak" ))
        file.rename( file, paste( file, "bak", sep = "." ))
    }
    write("#  123456789 -filetype Rmap.out", file, append=FALSE)

    # write numbers of progeny, markers and phenotypes
    write( "-s", file, append=TRUE)
    write( "-f 1", file, append=TRUE)
    write( "-p 0.0000", file, append=TRUE)
    write( "-u c", file, append=TRUE)
    write( "#", file, append=TRUE)

    write( paste( "-c", n.chr ), file, append=TRUE)
    write( paste( "-i", tot.mar ), file, append=TRUE)

    map <- lapply( cross$geno, function( x ) x$map )
    maplen <- unlist( lapply( map, length ))

    # mean and SD of number of markers
    write( paste( "-m", round( mean( maplen ), 3 )), file, append=TRUE)
    write( paste( "-vm", round( sqrt( var( maplen )), 3 )), file, append=TRUE)

    mapdif <- lapply( map, diff )
    # mean and SD of intermarker distances
    write( paste( "-d", round( mean( unlist( mapdif )), 3 )), file, append=TRUE)
    write( paste( "-vd", round( sqrt( var( unlist( mapdif ))), 3 )), file, append=TRUE)
    write( "-t 0.0000", file, append=TRUE)
    write( "#", file, append=TRUE)
    write( "          |   Chromosome----> ", file, append=TRUE)
    write( "--------------------------------------", file, append=TRUE)
    mapmat <- matrix( NA, 1 + max( maplen ), n.chr )
    mapmat[ 1, ] <- 0
    for( i in seq( along = maplen )) {
        tmp <- c( mapdif[[i]],0)
        mapmat[1 + seq( along = tmp ), i ] <- tmp
    }
    mapmat <- format(round(mapmat, 6))
    ncmap <- nchar( mapmat[1] )
    mapmat[ grep( "NA", mapmat ) ] <- paste( rep( " ", ncmap ), collapse = "" )
    tmp <- format( seq( n.chr ))
    write( paste( "Marker    |  ",
                 paste( paste(" ", tmp, sep=""), collapse =
                       paste( rep( " ", max( 1, ncmap - nchar( tmp ))), collapse = "" ))),
          file, append=TRUE)
    write( "--------------------------------------", file, append=TRUE)
    for( i in seq( nrow( mapmat ))) {
        j <- 5-nchar(i-1)
        if(j < 1) j <- 1
        write( paste( "-l", paste(rep(" ", j), collapse=""),
                     i - 1, "|",
                     paste( mapmat[i,], collapse = " " )),
              file, append=TRUE)
    }

    write( "---------------------------------------", file, append=TRUE)
    write( paste( "-Number   |", paste( maplen, collapse = "  " )),
          file, append=TRUE)

    write( "Names and positions of the markers", file, append=TRUE)
    write( "Chrom  Mark  Name", file, append=TRUE)
    write( "-b MarkerNames", file, append=TRUE)
    for( i in 1:n.chr )
        for( j in seq( along = map[[i]] ))
            write( paste( i, j, names( map[[i]] )[j] ), file, append=TRUE)
    write( "-e MarkerNames", file, append=TRUE)
    write( "Names of the Chromosomes", file, append=TRUE)
    write( "-b ChromosomeNames", file, append=TRUE)
    for( i in 1:n.chr )
        write( paste( i, names( map )[i] ), file, append=TRUE)
    write( "-e ChromosomeNames", file, append=TRUE)
}

# end of qtlcart_io.R
