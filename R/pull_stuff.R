#####################################################################
#
# pull_stuff.R
#
# copyright (c) 2001-2012, Karl W Broman
#     [find.pheno, find.flanking, and a modification to create.map
#      from Brian Yandell]
# last modified Mar, 2012
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
# Contains: pull.map, pull.geno, pull.pheno
#           pull.genoprob, pull.argmaxgeno
#
######################################################################

######################################################################
#
# pull.map
#
# pull out the map portion of a cross object, as a list
#
######################################################################

pull.map <-
    function(cross, chr, as.table=FALSE)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(!missing(chr)) cross <- subset(cross, chr=chr)
    if(!as.table) {
        result <- lapply(cross$geno,function(a) {
            b <- a$map
            class(b) <- as.character(class(a))
            b })
        class(result) <- "map"
        return(result)
    } else {
        return(map2table(pull.map(cross, as.table=FALSE)))
    }
}

map2table <-
    function(map, chr)
{
    if(!missing(chr)) {
        chr <- matchchr(chr, names(map))
        map <- map[chr]
    }

    if(is.matrix(map[[1]])) {
        map1 <- unlist(lapply(map, function(a) a[1,]))
        map2 <- unlist(lapply(map, function(a) a[2,]))
        result <- data.frame(chr=rep(names(map), vapply(map, ncol, 0)),
                             pos.female=map1, pos.male=map2, stringsAsFactors=FALSE)
        rownames(result) <- unlist(lapply(map, colnames))
    } else {
        result <- data.frame(chr=rep(names(map), vapply(map, length, 0)),
                             pos=unlist(map), stringsAsFactors=FALSE)
        rownames(result) <- unlist(lapply(map, names))
    }
    result[,1] <- factor(result[,1], levels=unique(result[,1]))

    result
}


######################################################################
# pull.geno
######################################################################
pull.geno <-
    function(cross, chr)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if(!missing(chr))
        cross <- subset(cross, chr=chr)

    X <- cross$geno[[1]]$data
    if(nchr(cross) > 1)
        for(i in 2:nchr(cross))
            X <- cbind(X, cross$geno[[i]]$data)
    X
}

######################################################################
# pull.pheno
######################################################################
pull.pheno <-
    function(cross, pheno.col)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    pheno <- cross$pheno

    if(!missing(pheno.col)) {

        if(is.character(pheno.col)) {
            m <- match(pheno.col, names(pheno))
            if(any(is.na(m))) {
                if(sum(is.na(m)) > 1)
                    warning("Phenotypes ", paste("\"", pheno.col[is.na(m)], "\"", sep="", collapse=" "), " not found.")
                else
                    warning("Phenotype ", paste("\"", pheno.col[is.na(m)], "\"", sep="", collapse=" "), " not found.")
            }
            if(all(is.na(m))) return(NULL)

            m <- m[!is.na(m)]
            pheno <- pheno[,m]
        }
        else if(is.logical(pheno.col)) {
            if(length(pheno.col) != ncol(pheno))
                stop("If pheno.col is logical, it should have length ", ncol(pheno))
            pheno <- pheno[,pheno.col]
        }
        else if(is.numeric(pheno.col)) {
            if(any(pheno.col > 0) && any(pheno.col < 0))
                stop("If pheno.col is numeric, values should be all > 0 or all < 0")
            if(any(pheno.col > 0) && (any(pheno.col < 1) || any(pheno.col > ncol(pheno))))
                stop("pheno.col values should be >= 1 and <= ", ncol(pheno))
            if(any(pheno.col < 0) && (any(pheno.col > -1) || any(pheno.col < -ncol(pheno))))
                stop("With negative pheno.col values, they should be between -", ncol(pheno), " and -1")
            pheno <- pheno[,pheno.col]
        }
    }

    if(is.data.frame(pheno) && ncol(pheno) == 1) pheno <- pheno[,1]

    pheno
}

######################################################################
# pull.genoprob
######################################################################
pull.genoprob <-
    function(cross, chr, omit.first.prob=FALSE, include.pos.info=FALSE, rotate=FALSE)
{
    if(!missing(chr))
        cross <- subset(cross, chr=chr)

    if(!("prob" %in% names(cross$geno[[1]])))
        stop("You must first run calc.genoprob.")

    if(include.pos.info && !rotate) {
        warning("If include.pos.info=TRUE, we assume rotate=TRUE as well.")
        rotate <- TRUE
    }

    pr <- lapply(cross$geno, function(a) a$prob)
    chrnames <- names(cross$geno)
    for(i in seq(along=pr)) {
        w <- colnames(pr[[i]])
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",chrnames[i],".",w[o],sep="")
        colnames(pr[[i]]) <- w
    }

    if(omit.first.prob)
        fullncol <- sum(sapply(pr, ncol))*(dim(pr[[1]])[3]-1)
    else
        fullncol <- sum(sapply(pr, ncol))*dim(pr[[1]])[3]

    fullpr <- matrix(nrow=nrow(pr[[1]]), ncol=fullncol)
    colnames(fullpr) <- 1:fullncol
    curcol <- 0
    thegen <- themarker <- rep(NA, fullncol)
    if(include.pos.info)
        thechr <- thepos <- rep(NA, fullncol)
    for(i in seq(along=pr)) {
        dim3 <- 1:dim(pr[[i]])[3]
        if(omit.first.prob) dim3 <- dim3[-1]
        for(j in seq(along=dim3)) {
            thecol <- curcol + ((1:ncol(pr[[i]]))-1)*length(dim3) + j
            fullpr[,thecol] <- pr[[i]][,,dim3[j]]
            thisgen <- dimnames(pr[[i]])[[3]][dim3[j]]
            thegen[thecol] <- rep(thisgen, length(thecol))
            themarker[thecol] <- colnames(pr[[i]])
            colnames(fullpr)[thecol] <- paste(colnames(pr[[i]]), thisgen, sep=":")
            if(include.pos.info) {
                thechr[thecol] <- rep(names(cross$geno)[i], length(thecol))
                thepos[thecol] <- attr(pr[[i]], "map")
            }
        }
        curcol <- curcol + ncol(pr[[i]])*length(dim3)
    }

    id <- getid(cross)
    if(is.null(id)) id <- paste("ind", 1:nrow(fullpr), sep="")
    rownames(fullpr) <- id

    if(rotate) {
        fullpr <- as.data.frame(t(fullpr))
        if(include.pos.info) {
            thechr <- factor(thechr, names(cross$geno))
            fullpr <- cbind(marker=themarker, gen=thegen, chr=thechr, pos=thepos, fullpr, stringsAsFactors=FALSE)
        }
    }

    fullpr
}


######################################################################
# pull.argmaxgeno
######################################################################
pull.argmaxgeno <-
    function(cross, chr, include.pos.info=FALSE, rotate=FALSE)
{
    if(!missing(chr))
        cross <- subset(cross, chr=chr)

    if(!("argmax" %in% names(cross$geno[[1]])))
        stop("You must first run argmax.geno.")

    if(include.pos.info && !rotate) {
        warning("If include.pos.info=TRUE, we assume rotate=TRUE as well.")
        rotate <- TRUE
    }

    am <- lapply(cross$geno, function(a) a$argmax)
    chrnames <- names(cross$geno)
    for(i in seq(along=am)) {
        w <- colnames(am[[i]])
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",chrnames[i],".",w[o],sep="")
        colnames(am[[i]]) <- w
    }

    fullncol <- sum(sapply(am, ncol))

    fullam <- matrix(nrow=nrow(am[[1]]), ncol=fullncol)
    colnames(fullam) <- 1:fullncol
    curcol <- 0
    if(include.pos.info)
        thechr <- thepos <- rep(NA, fullncol)
    for(i in seq(along=am)) {
        thecol <- curcol + 1:ncol(am[[i]])
        fullam[,thecol] <- am[[i]]
        colnames(fullam)[thecol] <- colnames(am[[i]])
        if(include.pos.info) {
            thechr[thecol] <- rep(names(cross$geno)[i], length(thecol))
            thepos[thecol] <- attr(am[[i]], "map")
        }
        curcol <- curcol + length(thecol)
    }

    id <- getid(cross)
    if(is.null(id)) id <- paste("ind", 1:nrow(fullam), sep="")
    rownames(fullam) <- id

    if(rotate) {
        fullam <- as.data.frame(t(fullam))
        if(include.pos.info) {
            thechr <- factor(thechr, names(cross$geno))
            fullam <- cbind(marker=rownames(fullam), chr=thechr, pos=thepos, fullam, stringsAsFactors=FALSE)
        }
    }

    fullam
}

######################################################################
# pull.draws
######################################################################
pull.draws <-
    function(cross, chr)
{
    if(!missing(chr))
        cross <- subset(cross, chr=chr)

    if(!("draws" %in% names(cross$geno[[1]])))
        stop("You must first run argmax.geno.")

    dr <- lapply(cross$geno, function(a) a$draws)
    chrnames <- names(cross$geno)
    for(i in seq(along=dr)) {
        w <- colnames(dr[[i]])
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",chrnames[i],".",w[o],sep="")
        colnames(dr[[i]]) <- w
    }

    fullncol <- sum(sapply(dr, ncol))
    d <- dim(dr[[1]])

    fulldr <- array(dim=c(d[1], fullncol, d[3]))
    colnames(fulldr) <- 1:fullncol
    curcol <- 0

    for(i in seq(along=dr)) {
        thecol <- curcol + 1:ncol(dr[[i]])
        fulldr[,thecol,] <- dr[[i]]
        colnames(fulldr)[thecol] <- colnames(dr[[i]])
        curcol <- curcol + length(thecol)
    }

    id <- getid(cross)
    if(is.null(id)) id <- paste("ind", 1:nrow(fulldr), sep="")
    rownames(fulldr) <- id

    fulldr
}

##############################
# table2map: create map object from a table
#
# rownames should be marker names
# first column chromosome
# second column position
##############################
table2map <-
    function(tab)
{
    mar <- rownames(tab)
    if(is.null(mar)) stop("marker names should be the row names")
    chr <- factor(tab[,1], levels=unique(tab[,1]))
    pos <- tab[,2]

    map <- split(pos, chr)
    mar <- split(mar, chr)
    for(i in seq(along=map))
        names(map[[i]]) <- mar[[i]]

    if(all(names(map) %in% c(1:20,"X"))) { # names are as in mouse
        for(i in seq(along=map))
            class(map[[i]]) <- ifelse(names(map)[i]=="X", "X", "A")
    }

    class(map) <- "map"
    map
}


# end of pull_stuff.R
