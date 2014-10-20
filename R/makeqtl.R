######################################################################
#
# makeqtl.R
#
# copyright (c) 2002-2013, Hao Wu and Karl W. Broman
# last modified Dec, 2013
# first written Apr, 2002
#
# Modified by Danny Arends
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
# Contains: makeqtl, replaceqtl, addtoqtl, dropfromqtl, locatemarker
#           print.qtl, summary.qtl, print.summary.qtl, reorderqtl
#           plot.qtl
#           print.compactqtl, summary.compactqtl, print.summary.compactqtl
#
######################################################################

######################################################################
#
# This is the function to construct an object of class "qtl"
# The phenotype data and genotype data for a given list of
# chromosome and locations will be extracted from the input
# "cross" object
#
######################################################################
makeqtl <-
    function(cross, chr, pos, qtl.name, what=c("draws", "prob"))
{
    if( !sum(class(cross) == "cross") )
        stop("The first input variable must be an object of class cross")

    # cross type
    type <- class(cross)[1]
    chrtype <- sapply(cross$geno, class)
    names(chrtype) <- names(cross$geno)
    sexpgm <- getsex(cross)

    what <- match.arg(what)

    themap <- pull.map(cross)

    # try to interpret chr argument
    if(!is.character(chr))
        chr <- as.character(chr)

    # chr, pos and qtl.name must have the same length
    if(length(chr) != length(pos))
        stop("Input chr and pos must have the same length.")
    else if( !missing(qtl.name) )
        if( length(chr) != length(qtl.name) )
            stop("Input chr and qtl.name must have the same length.")

    # local variables
    n.ind <- nrow(cross$pheno) # number of individuals
    n.pos <- length(chr) # number of selected markers
    n.gen <- NULL

    # initialize output object
    qtl <- NULL

    # take out the imputed genotypes and/or genoprobs for the
    # selected markers (if there are there)
    if(what == "draws") { # pull out draws
        if(!("draws" %in% names(cross$geno[[1]])))
            stop("You must first run sim.geno.")

        # take out imputed genotype data
        n.draws <- dim(cross$geno[[1]]$draws)[3] # number of draws

        # initialize geno matrix for selected markers
        geno <- array(rep(0, n.ind*n.pos*n.draws),
                      dim=c(n.ind, n.pos, n.draws))

        for(i in 1:n.pos) {
            # get the index for this chromosome
            i.chr <- which(chr[i]==names(cross$geno))
            if(length(i.chr) == 0) # no this chromosome in cross
                stop("There's no chromosome number ", chr[i], " in input cross object")
            i.pos <- pos[i] # marker position

            # make the genetic map for this chromosome
            if("map" %in% names(attributes(cross$geno[[i.chr]]$draws)))
                map <- attr(cross$geno[[i.chr]]$draws,"map")
            else {
                stp <- attr(cross$geno[[i.chr]]$draws, "step")
                oe <- attr(cross$geno[[i.chr]]$draws, "off.end")

                if("stepwidth" %in% names(attributes(cross$geno[[i.chr]]$draws)))
                    stpw <- attr(cross$geno[[i.chr]]$draws, "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[i.chr]]$map,stp,oe,stpw)
            }

            # pull out the female map if there are sex-specific maps
            if(is.matrix(map)) map <- map[1,]

            # locate this marker (given chromosome and position)
            marker.idx <- locatemarker(map, i.pos, i.chr, flag="draws")
            if(length(marker.idx) > 1)
                stop("Multiple markers at the same position; run jittermap.")

            # if everything is all right, take the genotype
            geno[,i,] <- cross$geno[[i.chr]]$draws[,marker.idx,]
            pos[i] <- map[marker.idx]

            # no. genotypes
            n.gen[i] <- length(getgenonames(type,chrtype[i.chr],"full",sexpgm, attributes(cross)))

            # Fix up X chromsome here
            if(chrtype[i.chr]=="X" && (type=="bc" || type=="f2"))
                geno[,i,] <- reviseXdata(type,"full",sexpgm,draws=geno[,i,,drop=FALSE],
                                         cross.attr=attributes(cross))
        }
        # give geno dimension names
        # the 2nd dimension called "Q1", "Q2", etc.
        dimnames(geno) <- list(NULL, paste("Q", 1:n.pos, sep=""), NULL)
        # output
        qtl$geno <- geno
    }
    else { # pull out probs

        if(!("prob" %in% names(cross$geno[[1]])))
            stop("You must first run calc.genoprob.")

        # initialize prob matrix
        prob <- vector("list",n.pos)

        # locate the marker
        for(i in 1:n.pos) {

            # get the index for this chromosome
            i.chr <- which(chr[i]==names(cross$geno))
            if(length(i.chr) == 0) # no this chromosome in cross
                stop("There's no chromosome number ", chr[i], " in input cross object")
            i.pos <- pos[i] # marker position


            if("map" %in% names(attributes(cross$geno[[i.chr]]$prob)))
                map <- attr(cross$geno[[i.chr]]$prob,"map")
            else {
                stp <- attr(cross$geno[[i.chr]]$prob, "step")
                oe <- attr(cross$geno[[i.chr]]$prob, "off.end")

                if("stepwidth" %in% names(attributes(cross$geno[[i.chr]]$prob)))
                    stpw <- attr(cross$geno[[i.chr]]$prob, "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[i.chr]]$map,stp,oe,stpw)
            }

            # pull out the female map if there are sex-specific maps
            if(is.matrix(map)) map <- map[1,]

            # locate this marker (given chromosome and position)
            marker.idx <- locatemarker(map, i.pos, i.chr, flag="prob")
            if(length(marker.idx) > 1)
                stop("Multiple markers at the same position; run jittermap.")

            # take genoprob
            if(chrtype[i.chr]=="X" && (type=="bc" || type=="f2")) { # fix X chromosome probs
                prob[[i]] <- reviseXdata(type, "full", sexpgm,
                                         prob=cross$geno[[i.chr]]$prob[,marker.idx,,drop=FALSE],
                                         cross.attr=attributes(cross))[,1,]
            }
            else
                prob[[i]] <- cross$geno[[i.chr]]$prob[,marker.idx,]

            pos[i] <- map[marker.idx]

            # no. genotypes
            n.gen[i] <- ncol(prob[[i]])
        }
        qtl$prob <- prob
    }

    if(missing(qtl.name))  { # no given qtl names
        dig <- 1

        if(what=="draws")
            step <- attr(cross$geno[[i.chr]]$draws, "step")
        else
            step <- attr(cross$geno[[i.chr]]$prob, "step")

        if(!is.null(step)) {
            if(step > 0) dig <- max(dig, -floor(log10(step)))
        }
        else {
            if(what=="draws")
                stepw <- attr(cross$geno[[i.chr]]$draws, "stepwidth")
            else
                stepw <- attr(cross$geno[[i.chr]]$prob, "stepwidth")

            if(!is.null(stepw) && stepw > 0) dig <- max(dig, -floor(log10(stepw)))
        }

        # make qtl names
        qtl.name <- paste( paste(chr,sep=""), charround(pos,dig), sep="@")
    }

    # output object
    qtl$name <- qtl.name
    qtl$altname <- paste("Q", 1:n.pos, sep="")
    qtl$chr <- chr
    qtl$pos <- pos
    qtl$n.qtl <- n.pos
    qtl$n.ind <- nind(cross)
    qtl$n.gen <- n.gen
    qtl$chrtype <- chrtype[qtl$chr]
    names(qtl$chrtype) <- NULL

    class(qtl) <- "qtl"
    attr(qtl, "map") <- themap

    qtl
}



######################################################################
#
# This is the function to replace one QTL by another.
#
######################################################################
replaceqtl <-
    function(cross, qtl, index, chr, pos, qtl.name, drop.lod.profile=TRUE)
{
    if(class(qtl) != "qtl")
        stop("qtl should have class \"qtl\".")

    if(any(index < 1 | index > qtl$n.qtl))
        stop("index should be between 1 and ", qtl$n.qtl)

    if(length(index) != length(chr) || length(index) != length(pos))
        stop("index, chr, and pos should all have the same length.")
    if(!missing(qtl.name) && length(index) != length(qtl.name))
        stop("index and qtl.name should have the same length.")

    if("geno" %in% names(qtl)) what <- "draws"
    else what <- "prob"

    if(missing(qtl.name))
        newqtl <- makeqtl(cross, chr, pos, what=what)
    else
        newqtl <- makeqtl(cross, chr, pos, qtl.name=qtl.name, what=what)

    if(what=="draws") {
        qtl$geno[,index,] <- newqtl$geno
    }
    else {
        qtl$prob[index] <- newqtl$prob
    }

    qtl$name[index] <- newqtl$name
    qtl$chr[index] <- newqtl$chr
    qtl$pos[index] <- newqtl$pos
    qtl$chrtype[index] <- newqtl$chrtype

    if(qtl$n.ind != newqtl$n.ind) stop("Mismatch in no. individuals")

    qtl$n.gen[index] <- newqtl$n.gen

    if(drop.lod.profile)
        attr(qtl, "lodprofile") <- NULL

    qtl
}


######################################################################
#
# This is the function to add a QTL to given qtl object
#
######################################################################

addtoqtl <-
    function(cross, qtl, chr, pos, qtl.name, drop.lod.profile=TRUE)
{
    if(class(qtl) != "qtl")
        stop("qtl should have class \"qtl\".")

    if("geno" %in% names(qtl)) what <- "draws"
    else what <- "prob"

    if(missing(qtl.name))
        newqtl <- makeqtl(cross, chr, pos, what=what)
    else
        newqtl <- makeqtl(cross, chr, pos, qtl.name=qtl.name, what=what)

    if(what=="draws") {
        do <- dim(qtl$geno)
        dn <- dim(newqtl$geno)
        if(do[1] != dn[1] || do[3] != dn[3])
            stop("Mismatch in number of individuals or number of imputations.")

        temp <- array(dim=c(do[1], do[2]+dn[2], do[3]))
        temp[,1:ncol(qtl$geno),] <- qtl$geno
        temp[,-(1:ncol(qtl$geno)),] <- newqtl$geno
        colnames(temp) <- paste("Q", 1:ncol(temp), sep="")
        qtl$geno <- temp
    }
    else {
        qtl$prob <- c(qtl$prob, newqtl$prob)
    }

    qtl$name <- c(qtl$name, newqtl$name)
    qtl$chr <- c(qtl$chr, newqtl$chr)
    qtl$pos <- c(qtl$pos, newqtl$pos)
    qtl$n.qtl <- qtl$n.qtl + newqtl$n.qtl
    qtl$altname <- paste("Q", 1:qtl$n.qtl, sep="")
    qtl$chrtype <- c(qtl$chrtype, newqtl$chrtype)
    if(qtl$n.ind != newqtl$n.ind)
        stop("Mismatch in no. individuals")
    qtl$n.gen <- c(qtl$n.gen, newqtl$n.gen)

    attr(qtl, "formula") <- NULL
    attr(qtl, "pLOD") <- NULL

    if(drop.lod.profile)
        attr(qtl, "lodprofile") <- NULL

    qtl
}

######################################################################
#
# This is the function to drop a QTL from a given qtl object
#
######################################################################
dropfromqtl <-
    function(qtl, index, chr, pos, qtl.name, drop.lod.profile=TRUE)
{
    if(class(qtl) != "qtl")
        stop("qtl should have class \"qtl\".")

    if(!missing(chr) || !missing(pos)) {
        if(missing(chr) || missing(pos))
            stop("Give both chr and pos, or give name, or give a numeric index")
        if(!missing(qtl.name) || !missing(index))
            stop("Give chr and pos or qtl.name or numeric index, but not multiple of these.")
        if(length(chr) != length(pos))
            stop("chr and pos must have the same lengths.")

        todrop <- NULL
        for(i in seq(along=chr)) {
            m <- which(qtl$chr == chr[i])

            if(length(m) < 1)
                stop("No QTL on chr ", chr[i], " in input qtl object.")

            for(j in seq(along=m)) {
                d <- abs(qtl$pos[m[j]] - pos[i])
                if(min(d) > 10) stop("No qtl near position ", pos[i], " on chr ", chr[i])

                wh <- m[d==min(d)]
                if(length(wh) > 1)
                    stop("Multiple QTL matching chr ", chr[i], " at pos ", pos[i])

                if(min(d) > 1)
                    warning("No QTL on chr ", chr[i], " exactly at ", pos[i],
                            "; dropping that at ", qtl$pos[wh])

                todrop <- c(todrop, wh)
            }
        }

        todrop <- unique(todrop)
    }
    else if(!missing(qtl.name)) {
        if(!missing(index))
            stop("Give chr and pos or qtl.name or numeric index, but not multiple of these.")

        m <- match(qtl.name, qtl$name)
        if(all(is.na(m))) # if no matches, try "altname"
            m <- match(qtl.name, qtl$altname)

        if(any(is.na(m)))
            warning("Didn't match QTL ", qtl.name[is.na(m)])

        todrop <- m[!is.na(m)]
    }
    else {
        if(missing(index))
            stop("Give chr and pos or qtl.name or numeric index, but not multiple of these.")
        if(any(index < 1 | index > qtl$n.qtl))
            stop("index should be between 1 and ", qtl$n.qtl)

        todrop <- index
    }

    # input drop is an integer index
    # get the index for exclusing drop QTL
    idx <- setdiff(1:qtl$n.qtl, todrop)

    # result object
    qtl$name <- qtl$name[idx]
    qtl$chr <- qtl$chr[idx]
    qtl$chrtype <- qtl$chrtype[idx]
    qtl$pos <- qtl$pos[idx]
    qtl$n.qtl <- qtl$n.qtl - length(todrop)
    qtl$altname <- paste("Q", 1:qtl$n.qtl, sep="")
    qtl$n.ind <- qtl$n.ind
    qtl$n.gen <- qtl$n.gen[idx]
    if("geno" %in% names(qtl)) {
        qtl$geno <- qtl$geno[,idx,,drop=FALSE]
        colnames(qtl$geno) <- paste("Q", 1:ncol(qtl$geno), sep="")
    }
    if("prob" %in% names(qtl))
        qtl$prob <- qtl$prob[idx]

    attr(qtl, "formula") <- NULL
    attr(qtl, "pLOD") <- NULL

    if(drop.lod.profile)
        attr(qtl, "lodprofile") <- NULL

    qtl
}


##################################################################
#
# locate the marker on a genetic map. Choose the nearest
# one if there's no marker or pseudomarker one the given
# location
#
# This is the internal function and not supposed to be used by user
#
###################################################################

locatemarker <-
    function(map, pos, chr, flag)
{
    marker.idx <- which(map == pos)
    if( length(marker.idx)==0 ) {
        # there's no this marker, take the nearest marker instead
        # if there's a tie, take the first nearst one
        m.tmp <- abs(pos-map)
        marker.idx <- which(m.tmp==min(m.tmp))[[1]]
    }

    if(length(marker.idx) > 1)
        marker.idx <- marker.idx[sample(length(marker.idx))]
    marker.idx
}

# print QTL object
print.qtl <-
    function(x, ...)
{
    print(summary(x))
}

# summary of QTL object
summary.qtl <-
    function(object, ...)
{
    if(is.null(object) || length(object) == 0 || length(object$chr)==0) {
        object <- numeric(0)
        class(object) <- "summary.qtl"
        return(object)
    }

    if("geno" %in% names(object)) {
        type <- "draws"
        n.draws <- dim(object$geno)[3]
    }
    else type <- "prob"

    output <- data.frame(name=object$name, chr=object$chr, pos=object$pos, n.gen=object$n.gen, stringsAsFactors=TRUE)
    rownames(output) <- object$altname

    attr(output, "type") <- type
    if(!is.null(attr(object,"mqm"))) attr(output, "mqm") <- attr(object,"mqm")
    if(type=="draws") attr(output, "n.draws") <- n.draws
    class(output) <- c("summary.qtl", "data.frame")

    if("formula" %in% names(attributes(object)))
        attr(output, "formula") <- attr(object, "formula")
    if("pLOD" %in% names(attributes(object)))
        attr(output, "pLOD") <- attr(object, "pLOD")

    output
}

# print summary of QTL object
print.summary.qtl <-
    function(x, ...)
{
    if(is.null(x) || length(x) == 0) {
        cat("  Null QTL model\n")
    }
    else {
        type <- attr(x, "type")
        if(type=="draws")
            thetext <- paste("imputed genotypes, with", attr(x, "n.draws"), "imputations.")
        else thetext <- "genotype probabilities."
        if(!is.null(attr(x,"mqm"))) thetext <- paste("model created by mqmscan")
        cat("  QTL object containing", thetext, "\n\n")

        print.data.frame(x, digits=5)
    }

    if("formula" %in% names(attributes(x))) {
        form <- attr(x, "formula")
        if(!is.character(form)) form <- deparseQTLformula(form)
        cat("\n  Formula:")
        w <- options("width")[[1]]
        printQTLformulanicely(form, "               ", w+5, w)
    }

    if("pLOD" %in% names(attributes(x)))
        cat("\n  pLOD: ", round(attr(x, "pLOD"),3), "\n")
}

######################################################################
# plot locations of QTLs on the genetic map
######################################################################
plot.qtl <-
    function(x, chr, horizontal=FALSE, shift=TRUE,
             show.marker.names=FALSE,  alternate.chrid=FALSE,
             justdots=FALSE, col="red", ...)
{
    if(!("qtl" %in% class(x)))
        stop("input should be a qtl object")

    if(length(x) == 0)
        stop("  There are no QTL to plot.")

    map <- attr(x, "map")
    if(is.null(map))
        stop("qtl object doesn't contain a genetic map.")

    if(missing(chr))
        chr <- names(map)
    else {
        chr <- matchchr(chr, names(map))
        map <- map[chr]
        class(map) <- "map"
    }

    if(horizontal)
        plotMap(map, horizontal=horizontal, shift=shift,
                show.marker.names=show.marker.names, alternate.chrid=alternate.chrid,
                ylim=c(length(map)+0.5, 0), ...)
    else
        plotMap(map, horizontal=horizontal, shift=shift,
                show.marker.names=show.marker.names, alternate.chrid=alternate.chrid,
                xlim=c(0.5,length(map)+1), ...)

    whchr <- match(x$chr, names(map))
    thepos <- x$pos
    thepos[is.na(whchr)] <- NA

    if(any(!is.na(thepos))) {
        whchr <- whchr[!is.na(whchr)]

        if(shift) thepos <- thepos - sapply(map[whchr], min)

        if(is.matrix(map[[1]])) whchr <- whchr - 0.3

        if(length(grep("^.+@[0-9\\.]+$", x$name)) == length(x$name))
            x$name <- x$altname

        if(!justdots) {
            if(horizontal) {
                arrows(thepos, whchr - 0.35, thepos, whchr, lwd=2, col=col, length=0.05)
                text(thepos, whchr-0.4, x$name, col=col, adj=c(0.5,0))
            }
            else {
                arrows(whchr + 0.35, thepos, whchr, thepos, lwd=2, col=col, length=0.05)
                text(whchr+0.4, thepos, x$name, col=col, adj=c(0,0.5))
            }
        }
        else {
            if(horizontal)
                points(thepos, whchr, pch=16, col=col)
            else {
                points(whchr, thepos, pch=16, col=col)
            }
        }
    }

    invisible()
}

######################################################################
#
# This is the function to reorder the QTL within a QTL object
#
######################################################################
reorderqtl <-
    function(qtl, neworder)
{
    if(class(qtl) != "qtl")
        stop("qtl should have class \"qtl\".")

    if(missing(neworder)) {
        if(!("map" %in% names(attributes(qtl))))
            stop("No map in the qtl object; you must provide argument 'neworder'.")
        chr <- names(attr(qtl, "map"))
        thechr <- match(qtl$chr, chr)
        if(any(is.na(thechr)))
            stop("Chr ", paste(qtl$chr[is.na(thechr)], " "), " not found.")
        neworder <- order(thechr, qtl$pos)
    }

    curorder <- seq(qtl$n.qtl)
    if(length(neworder) != qtl$n.qtl ||
       !all(curorder == sort(neworder)))
        stop("neworder should be an ordering of the integers from 1 to ", qtl$n.qtl)

    if(qtl$n.qtl == 1)
        stop("Nothing to do; just one qtl.")


    if("geno" %in% names(qtl))
        qtl$geno <- qtl$geno[,neworder,]
    else
        qtl$prob <- qtl$prob[neworder]

    qtl$name <- qtl$name[neworder]
    qtl$chr <- qtl$chr[neworder]
    qtl$pos <- qtl$pos[neworder]
    qtl$n.gen <- qtl$n.gen[neworder]
    qtl$chrtype <- qtl$chrtype[neworder]

    attr(qtl, "formula") <- NULL
    attr(qtl, "pLOD") <- NULL

    if("lodprofile" %in% names(attributes(qtl))) {
        lodprof <- attr(qtl, "lodprofile")
        if(length(lodprof) == length(neworder))
            attr(qtl, "lodprofile") <- lodprof[neworder]
    }

    qtl
}

# print compact version of QTL object
print.compactqtl <-
    function(x, ...)
{
    print(summary(x))
}

summary.compactqtl <-
    function(object, ...)
{
    class(object) <- c("summary.compactqtl", "list")
    object
}

print.summary.compactqtl <-
    function(x, ...)
{
    if(is.null(x) || length(x) == 0)
        cat("Null QTL model\n")
    else {
        temp <- as.data.frame(x, stringsAsFactors=TRUE)
        rownames(temp) <- paste("Q", 1:nrow(temp), sep="")
        print.data.frame(temp)
    }
    if("formula" %in% names(attributes(x))) {
        form <- attr(x, "formula")
        if(!is.character(form)) form <- deparseQTLformula(form)
        cat("  Formula:")
        w <- options("width")[[1]]
        printQTLformulanicely(form, "               ", w+5, w)
    }

    if("pLOD" %in% names(attributes(x)))
        cat("  pLOD: ", round(attr(x, "pLOD"),3), "\n")
}

# nqtl: number of qtl
nqtl <- function(qtl) ifelse(length(qtl)==0,0,length(qtl$chr))

# end of makeqtl.R
