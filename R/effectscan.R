######################################################################
#
# effectscan.R
#
# copyright (c) 2003-2013, Karl W. Broman
# [completely re-written in Sep, 2007, based partly on code from Hao Wu]
# last modified Sep, 2013
# first written Jan, 2003
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
# Contains: effectscan
#
######################################################################

effectscan <-
    function(cross, pheno.col=1, chr, get.se=FALSE, draw=TRUE,
             gap=25, ylim, mtick=c("line","triangle"),
             add.legend=TRUE, alternate.chrid=FALSE, ...)
{
    type <- class(cross)[1]
    mtick <- match.arg(mtick)
    if(type == "4way")
        stop("effect scan not working for 4-way cross yet.")

    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("effectscan can take just one phenotype; only the first will be used")
    }

    if(is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if(is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col, "\"")
        pheno.col <- num
    }

    if(pheno.col < 1 | pheno.col > nphe(cross))
        stop("pheno.col values should be between 1 and the no. phenotypes")

    if(!is.numeric(cross$pheno[,pheno.col]))
        stop("phenotype \"", colnames(cross$pheno)[pheno.col], "\" is not numeric.")

    pheno <- cross$pheno[,pheno.col]
    wh <- is.na(pheno)
    if(any(wh)) {
        pheno <- pheno[!wh]
        cross <- subset.cross(cross, ind=(!wh))
    }

    if(!missing(chr)) cross <- subset.cross(cross, chr=chr)

    chrtype <- sapply(cross$geno, class)

    n.ind <- length(pheno)

    results <- NULL
    for(i in 1:nchr(cross)) {
        if(!("draws" %in% names(cross$geno[[i]])))
            stop("You must first run sim.geno.")

        draws <- cross$geno[[i]]$draws

        # create map of positions
        if("map" %in% names(attributes(cross$geno[[i]]$draws)))
            map <- attr(cross$geno[[i]]$draws,"map")
        else {
            stp <- attr(cross$geno[[i]]$draws, "step")
            oe <- attr(cross$geno[[i]]$draws, "off.end")

            if("stepwidth" %in% names(attributes(cross$geno[[i]]$draws)))
                stpw <- attr(cross$geno[[i]]$draws, "stepwidth")
            else stpw <- "fixed"
            map <- create.map(cross$geno[[i]]$map,stp,oe,stpw)
        }

        if(is.matrix(map)) {
            marnam <- colnames(map)
            map <- map[1,]
        }
        else marnam <- names(map)

        if(type == "risib" || type=="riself" || type=="dh" || type=="haploid") {
            mapping <- rbind(c(+1, -1),
                             c(+1, +1))
            colnames(mapping) <- c("intercept","a")
            dropcol <- 1
        }
        else if(type=="bc") {
            if(chrtype[i] == "X") {
                sexpgm <- getsex(cross)
                draws <- reviseXdata(type, "full", sexpgm, draws=draws,
                                     cross.attr=attributes(cross))

                if(is.null(sexpgm$sex) || all(sexpgm$sex==0) || all(sexpgm$sex==1)) { # all one sex
                    mapping <- rbind(c(+1, -0.5),
                                     c(+1, +0.5))
                    colnames(mapping) <- c("intercept", "a")
                    dropcol <- 1
                }
                else { # some of each sex
                    mapping <- rbind(c(+1, 0,-0.5,   0),
                                     c(+1, 0,+0.5,   0),
                                     c(+1,+1, 0,  -0.5),
                                     c(+1,+1, 0,  +0.5))
                    colnames(mapping) <- c("intercept", "sex", "a.female", "a.male")
                    dropcol <- 1:2
                }
            } # end bc X chr
            else {
                mapping <- rbind(c(+1, -0.5),
                                 c(+1, +0.5))
                colnames(mapping) <- c("intercept", "a")
                dropcol <- 1
            } # end bc autosome
        } # end bc
        else { # intercross
            if(chrtype[i] == "X") {
                sexpgm <- getsex(cross)
                draws <- reviseXdata(type, "full", sexpgm, draws=draws,
                                     cross.attr=attributes(cross))


                if(is.null(sexpgm$pgm) || all(sexpgm$pgm==0) || all(sexpgm$pgm==1)) { # all one direction

                    if(is.null(sexpgm$sex) || all(sexpgm$sex==0) || all(sexpgm$sex==1)) { # all one sex
                        mapping <- rbind(c(+1, -0.5),
                                         c(+1, +0.5))
                        colnames(mapping) <- c("intercept", "a")
                        dropcol <- 1
                    }
                    else {
                        mapping <- rbind(c(+1, 0,-0.5,   0),
                                         c(+1, 0,+0.5,   0),
                                         c(+1,+1, 0,  -0.5),
                                         c(+1,+1, 0,  +0.5))
                        colnames(mapping) <- c("intercept", "sex", "a.female", "a.male")
                        dropcol <- 1:2
                    }
                }
                else { # some of each direction
                    if(is.null(sexpgm$sex) || all(sexpgm$sex==0)) { # all female
                        mapping <- rbind(c(+1, 0,-0.5,  0),
                                         c(+1, 0,+0.5,  0),
                                         c(+1,+1,   0,-0.5),
                                         c(+1,+1,   0,+0.5))
                        colnames(mapping) <- c("intercept","dir","a.forw","a.rev")
                        dropcol <- 1:2
                    }
                    else if(all(sexpgm$sex==1)) { # all male
                        mapping <- rbind(c(+1, -0.5),
                                         c(+1, +0.5))
                        colnames(mapping) <- c("intercept", "a")
                        dropcol <- 1
                    }
                    else { # some of each sex
                        mapping <- rbind(c(+1, 0, 0, -0.5,  0,   0),
                                         c(+1, 0, 0, +0.5,  0,   0),
                                         c(+1,+1, 0,   0,-0.5,   0),
                                         c(+1,+1, 0,   0,+0.5,   0),
                                         c(+1, 0,+1,   0,   0,-0.5),
                                         c(+1, 0,+1,   0,   0,+0.5))
                        colnames(mapping) <- c("intercept","dir","sex","a.femaleforw","a.femalerev","a.male")
                        dropcol <- 1:3
                    }
                }
            } # end f2 X chr
            else {
                mapping <- rbind(c(+1, -1,  0),
                                 c(+1,  0, +1),
                                 c(+1, +1,  0))
                colnames(mapping) <- c("intercept","a","d")
                dropcol <- 1
            } # f2 autosome
        } # end f2

        n.gen <- ncol(mapping)
        n.pos <- ncol(draws)
        n.imp <- dim(draws)[3]

        z <- .C("R_effectscan",
                as.integer(n.ind),
                as.integer(n.gen),
                as.integer(n.imp),
                as.integer(n.pos),
                as.integer(draws-1),
                as.double(pheno),
                as.double(mapping),
                beta=as.double(rep(0,n.pos*n.gen)),
                se=as.double(rep(0,n.pos*n.gen)),
                as.integer(get.se),
                PACKAGE="qtl")


        beta <- t(matrix(z$beta, ncol=n.pos))
        colnames(beta) <- colnames(mapping)

        if(get.se) {
            se <- t(matrix(z$se, ncol=n.pos))
            colnames(se) <- paste("se", colnames(mapping), sep=".")
            beta <- cbind(beta, se[,-dropcol,drop=FALSE])
        }
        z <- beta[,-dropcol,drop=FALSE]

        w <- marnam
        o <- grep("^loc-*[0-9]+",w)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
        rownames(z) <- w

        z <- as.data.frame(z, stringsAsFactors=TRUE)
        z <- cbind(chr=rep(names(cross$geno)[i],length(map)),
                   pos=as.numeric(map), z)
        rownames(z) <- w

        if(i==1) results <- z
        else {
            w <- match(colnames(z), colnames(results))

            if(any(is.na(w))) {
                curnam <- colnames(results)
                for(j in which(is.na(w)))
                    results <- cbind(results, rep(NA, nrow(results)))
                colnames(results) <- c(curnam, colnames(z)[is.na(w)])
            }

            w <- match(colnames(results), colnames(z))
            if(any(is.na(w))) {
                curnam <- colnames(z)
                for(j in which(is.na(w)))
                    z <- cbind(z, rep(NA, nrow(z)))
                colnames(z) <- c(curnam, colnames(results)[is.na(w)])
            }

            results <- rbind(results, z)

        }
    } # end  loop over chromosomes

    class(results) <- c("effectscan", "scanone", "data.frame")

    if(draw) { # make the figure
        if(missing(ylim))
            plot.effectscan(results, gap=gap, mtick=mtick, add.legend=add.legend,
                            alternate.chrid=alternate.chrid, ...)
        else
            plot.effectscan(results, gap=gap, mtick=mtick, add.legend=add.legend,
                            ylim=ylim, alternate.chrid=alternate.chrid, ...)
    }

    invisible(results)
}

# function to make the effectscan plot
plot.effectscan <-
    function(x, gap=25, ylim, mtick=c("line","triangle"),
             add.legend=TRUE, alternate.chrid=FALSE, ...)
{
    col <- c("blue","red","darkorange","darkgreen","purple")
    lightcol <- c("lightblue", "pink", "peachpuff1", "palegreen1", "thistle1")

    results <- x
    eff <- 3:ncol(results)

    if(length(grep("^se", colnames(results)))>0) get.se <- TRUE
    else get.se <- FALSE

    if(get.se) {
        se <- grep("^se", colnames(results)[eff])
        eff <- eff[-se]
        se <- se + 2

        lo <- as.matrix(results[,eff]) - as.matrix(results[,se])
        hi <- as.matrix(results[,eff]) + as.matrix(results[,se])

        yl <- range(c(lo,hi), na.rm=TRUE)
    }
    else yl <- range(results[,eff], na.rm=TRUE)

    if(!missing(ylim)) yl <- ylim

    plot.scanone(results, lodcolumn=1, ylim=yl, gap=gap, mtick=mtick, alternate.chrid=alternate.chrid,
                 col=col[1], ...)

    if(get.se) {
        begend <- matrix(unlist(tapply(results[,2],results[,1],range)),ncol=2,byrow=TRUE)
        rownames(begend) <- unique(results[,1])
        chr <- unique(as.character(results[,1]))

        begend <- begend[as.character(chr),,drop=FALSE]
        len <- begend[,2]-begend[,1]

        if(length(len)>1) start <- c(0,cumsum(len+gap))-c(begend[,1],0)
        else start <- 0

        x <- results[,2]
        for(i in seq(along=chr))
            x[results[,1]==chr[i]] <- results[results[,1]==chr[i],2]+start[i]

        for(i in seq(along=chr)) {
            wh <- results[,1]==chr[i]

            for(j in 1:ncol(lo)) {
                if(any(!is.na(lo[wh,j]))) {
                    xx <- c(x[wh], rev(x[wh]))
                    yy <- c(lo[wh,j], rev(hi[wh,j]))
                    polygon(xx, yy, col=lightcol[j], border=lightcol[j])
                }
            }
        }

        # go back and add lines at edges
        for(i in seq(along=chr)) {
            wh <- results[,1]==chr[i]

            for(j in 1:ncol(lo)) {
                if(any(!is.na(lo[wh,j]))) {
                    xx <- c(x[wh], rev(x[wh]))
                    yy <- c(lo[wh,j], rev(hi[wh,j]))
                    lines(xx, yy, col=lightcol[j])
                }
            }
        }

        plot.scanone(results, lodcolumn=1, add=TRUE, col=col[1])
    }

    if(length(eff) > 1) {
        for(i in seq(along=eff)[-1])
            plot.scanone(results, lodcolumn=eff[i]-2, gap=gap, add=TRUE, col=col[i])
    }

    if(add.legend)
        legend("top", legend=names(results)[eff], col=col[1:length(eff)], lwd=2)

    abline(h=0, lty=2)
}



# end of effectscan.R
