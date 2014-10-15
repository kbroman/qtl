######################################################################
#
# effectplot.R
#
# copyright (c) 2002-2012, Hao Wu and Karl W. Broman
# Last modified Sep, 2012
# first written Jul, 2002
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
# Modified by Hao Wu Feb 2005 for the following:
# 1. function will take marker, pseudomarker or phenotype as input;
# 2. separate functions to extract marker genodata given marker names
# and calculate means and ses;
#
# Part of the R/qtl package
# Contains: effectplot, effectplot.getmark, effectplot.calmeanse
#
######################################################################

effectplot <-
    function (cross, pheno.col = 1, mname1, mark1, geno1, mname2,
              mark2, geno2, main, ylim, xlab, ylab, col, add.legend = TRUE,
              legend.lab, draw=TRUE, var.flag=c("pooled","group"))
{
    if(!sum(class(cross) == "cross"))
        stop("The first input variable must be an object of class cross")

    if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }

    # if mname2 is given but not mname1, switch it around
    if((missing(mname1) && missing(mark1) && missing(geno1)) &&
       any(!missing(mname2) || missing(mark2) || missing(geno2))) {
        if(!missing(mname2)) { mname1 <- mname2; mname2 <- NULL }
        if(!missing(mark2)) { mark1 <- mark2; mark2 <- NULL }
        if(!missing(geno2)) { geno1 <- geno2; geno2 <- NULL }
    }
    if(missing(mark2)) mark2 <- NULL
    if(missing(mname2)) mname2 <- NULL

    if(length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("effectplot can take just one phenotype; only the first will be used")
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

    var.flag <- match.arg(var.flag)

    # local variables
    n.ind <- nind(cross)
    pheno <- cross$pheno[, pheno.col]
    type <- class(cross)[1]
    chrtype1 <- chrtype2 <- "A"
    gennames1 <- gennames2 <- NULL

    # If imputations are not available, create them
    if(!("draws" %in% names(cross$geno[[1]]))) {
        warning(" -Running sim.geno.")
        cross <- sim.geno(cross, n.draws=16)
    }

    ####################################################
    # get genotype data for markers given marker name
    # if marker genodata were given, this will be skipped
    ####################################################

    # used for alternative print name for pseudomarkers
    pm.pattern <- "^c.*\\.loc.*$" # pseudomarker names will be like "c1.loc10"
    dig <- 1
    step <- attr(cross$geno[[1]]$draws, "step")
    if(!is.null(step)) {
        if(step > 0) dig <- max(dig, -floor(log10(step)))
    }
    else {
        stepw <- attr(cross$geno[[1]]$draws, "stepwidth")
        if(!is.null(stepw) && stepw > 0) dig <- max(dig, -floor(log10(stepw)))
    }

    printname1 <- printname2 <- NULL
    # Get marker 1 genotype data
    if(missing(mark1)) { # no data given
        if(missing(mname1)) # no marker data or marker name, have to stop
            stop("Either mname1 or mark1 must be specified.")
        # get marker data according to marker name
        tmp <- effectplot.getmark(cross, mname1)
        tmptmp <- attr(tmp, "mname")
        if(!is.null(tmptmp)) mname1 <- tmptmp
        mark1 <- tmp$mark
        gennames1 <- tmp$genname

        # perhaps alternative print name
        if(length(grep(pm.pattern, mname1))>0) {
            tmp <- find.pseudomarkerpos(cross, mname1, "draws")
            printname1 <- paste(tmp[1,1], charround(tmp[1,2],dig), sep="@")
        }
    }
    else {
        # make mark1 a matrix if it's not
        if(class(mark1) != "matrix")
            mark1 <- matrix(mark1, ncol=1)
        if(dim(mark1)[1] != n.ind)
            stop("Marker 1 data hass the wrong dimension")
        if(missing(mname1))
            mname1 <- "Marker 1"
    }
    if(is.null(printname1)) printname1 <- mname1

    # Deal with marker 2
    if(!is.null(mname2) || !is.null(mark2)) {
        if(is.null(mark2)) {
            # get marker data according to marker name
            tmp <- effectplot.getmark(cross, mname2)
            tmptmp <- attr(tmp, "mname")
            if(!is.null(tmptmp)) mname2 <- tmptmp
            mark2 <- tmp$mark
            gennames2 <- tmp$genname

            # perhaps alternative print name
            if(length(grep(pm.pattern, mname2))>0) {
                tmp <- find.pseudomarkerpos(cross, mname2, "draws")
                printname2 <- paste(tmp[1,1], charround(tmp[1,2],dig), sep="@")
            }
        }
        else { # mark2 data is given
            # make mark2 a matrix if it's not
            if(class(mark2) != "matrix")
                mark2 <- matrix(mark2, ncol=1)
            if(dim(mark2)[1] != n.ind)
                stop("Marker 2 data has the wrong dimension")
            if(is.null(mname2))
                mname2 <- "Marker 2"
        }
        if(is.null(printname2)) printname2 <- mname2
    }
    else {
        mark2 <- NULL
        geno2 <- NULL
    }
    ### till now, mark1 and mark2 are genotype data in matrix

    ########################################################
    # deal with the data - if one of them is a pseudomarker,
    # make the other one the same number of draws
    ########################################################
    # determine number of draws - this part of codes works even if mark2 is NULL
    ndraws1 <- dim(mark1)[2]
    if(is.null(mark2))
        ndraws2 <- 1
    else
        ndraws2 <- dim(mark2)[2]

    # make them the same number of draws
    if( (ndraws1>1) && (ndraws2>1) ) {
        # two pseudomarkers, they must have the same number of draws
        if(ndraws1 != ndraws2)
            stop("Input two pseudomarkers have different number of draws.")
        else
            ndraws <- ndraws1
    }
    else if( (ndraws1>1) && (ndraws2==1) ) {
        # one pm and one typed marker
        if(!is.null(mark2))
            mark2 <- matrix(rep(mark2,ndraws1), ncol=ndraws1)
        ndraws <- ndraws1
    }
    else if( (ndraws1==1) && (ndraws2>1) ) {
        # one pm and one typed marker
        mark1 <- matrix(rep(mark1,ndraws2), ncol=ndraws2)
        ndraws <- ndraws2
    }
    else # they are all real markers
        ndraws <- 1

    # drop data for individuals with missing phenotypes or genotypes
    keepind <- !is.na(pheno)
    if(!is.null(mark1))
        keepind <- keepind & apply(mark1, 1, function(a) all(!is.na(a)))
    if(!is.null(mark2))
        keepind <- keepind & apply(mark2, 1, function(a) all(!is.na(a)))

    mark1 <- mark1[keepind,]
    mark2 <- mark2[keepind,]
    pheno <- pheno[keepind]

    ########################################################
    # 1. get level names - this part will be executed when
    # user only input mark without mname and geno
    # 2. adjust marker data if the input is not numeric
    ########################################################
    tmpf <- factor(mark1)
    if(!missing(geno1)) { # geno1 is given
        # check if it has the correct length
        if(length(geno1) < length(levels(tmpf)))
            stop("geno1 is too short.")
    }
    else {
        # geno1 is not given
        if(!is.null(gennames1)) # get it from genname1
            geno1 <- gennames1
        else if(is.factor(mark1)) { # or if it's factor, get it from level
            geno1 <- levels(mark1)
            mark1 <- as.numeric(mark1)
        }
        else if(!is.numeric(mark1)) {
            # if it's neither factor nor numeric, it must be a string vector
            # such like c("F","M","F")...
            geno1 <- levels(tmpf)
        }
        else { # otherwise, generate a standard one
            geno1 <- getgenonames(type, "A", cross.attr=attributes(cross))
            if(length(levels(tmpf)) > length(geno1))
                geno1 <- c(geno1, rep("?", length(levels(tmpf)) -
                                      length(geno1)))
        }
    }
    # adjust marker data - if the input is not numeric, convert them into numeric
    if(!is.numeric(mark1))
        mark1 <- matrix(as.numeric(tmpf, levels=sort(levels(tmpf))), ncol=ndraws)

    # Now work on mark2
    if(!is.null(mark2)) {
        tmpf <- factor(mark2)
        if(!missing(geno2)) { # geno2 is given
            # check if it has the correct length
            if(length(geno2) < length(levels(tmpf)))
                stop("geno2 is too short.")
        }
        else {
            # geno2 is not given
            if(!is.null(gennames2)) # get it from genname2
                geno2 <- gennames2
            else if(is.factor(mark2)) { # or if it's factor, get it from level
                geno2 <- levels(mark2)
                mark2 <- as.numeric(mark2)
            }
            else if(!is.numeric(mark2)) {
                # if it's neither factor nor numberic, it must be a string vector
                # such like c("F","M","F")...
                geno2 <- levels(tmpf)
            }
            else { # otherwise, generate a standard one
                geno2 <- getgenonames(type, "A", cross.attr=attributes(cross))
                if(length(levels(tmpf)) > length(geno2))
                    geno2 <- c(geno2, rep("?", length(levels(tmpf)) -
                                          length(geno2)))
            }
        }
        # adjust marker data - if the input is not numeric, convert them into numeric
        if(!is.numeric(mark2))
            mark2 <- matrix(as.numeric(tmpf, levels=sort(levels(tmpf))), ncol=ndraws)
    }

    # number of genotypes
    ngen1 <- length(geno1)
    if(!is.null(mark2))
        ngen2 <- length(geno2)

    # calculate means and ses
    # and make output object
    # the output will be a data frame. For two-marker case,
    # the rows corresponding to the first marker and the columns
    # corresponding to the second marker
    result <- effectplot.calmeanse(pheno, mark1, mark2, geno1, geno2, ndraws, var.flag)
    means <- result$Means
    ses <- result$SEs

    # assign column and row names
    if(is.null(mark2)) {
        if(length(means) != length(geno1)) {
            warning("Number of genotypes is different than length(geno1).")
            if(length(means) < length(geno1))
                geno1 <- geno1[1:length(means)]
            else geno1 <- c(geno1, rep("?", length(means) - length(geno1)))
            ngen1 <- length(geno1)
        }
        names(result$Means) <- paste(printname1, geno1, sep = ".")
        names(result$SEs) <- paste(printname1, geno1, sep = ".")
    }
    else {
        if(nrow(means) != length(geno1)) {
            warning("Number of genotypes in marker 1 is different than length(geno1).")
            if(nrow(means) < length(geno1))
                geno1 <- geno1[1:nrow(means)]
            else geno1 <- c(geno1, rep("?", nrow(means) - length(geno1)))
            ngen1 <- length(geno1)
        }
        if(ncol(means) != length(geno2)) {
            warning("Number of genotypes in marker 2 is different than length(geno2).")
            if(ncol(means) < length(geno2))
                geno2 <- geno2[1:ncol(means)]
            else geno2 <- c(geno2, rep("?", ncol(means) - length(geno2)))
            ngen2 <- length(geno2)
        }
        rownames(result$Means) <- paste(printname1, geno1, sep = ".")
        colnames(result$Means) <- paste(printname2, geno2, sep = ".")
        rownames(result$SEs) <- paste(printname1, geno1, sep = ".")
        colnames(result$SEs) <- paste(printname2, geno2, sep = ".")
    }

    # calculate lo's and hi's for plot
    lo <- means - ses
    hi <- means + ses

    ######### Draw the figure if requested ############
    if(draw) {
        # graphics parameters
        old.xpd <- par("xpd")
        old.las <- par("las")
        par(xpd = FALSE, las = 1)
        on.exit(par(xpd = old.xpd, las = old.las))

        # colors (for case of two markers)
        if(missing(col)) {
            if(ngen1 <= 5) {
                if(ngen1 == 1) int.color <- "black"
                else if(ngen1 == 2) int.color <- c("red", "blue")
                else int.color <- c("black", "red", "blue", "orange", "green")[1:ngen1]
            }
            else
                int.color <- c("black", rainbow(ngen1-1, start=0, end=2/3))
        }
        else int.color <- col

        # plot title
        if(missing(main)) {
            if(is.null(mark2))
                main <- paste("Effect plot for", printname1)
            else main <- paste("Interaction plot for", printname1, "and",
                               printname2)
        }

        # y axis limits
        if(missing(ylim)) {
            ylimits <- range(c(lo, means, hi), na.rm = TRUE)
            ylimits[2] <- ylimits[2] + diff(ylimits) * 0.1
        }
        else ylimits <- ylim

        # x axis limits
        if(is.null(mark2)) { # one marker
            u <- seq(along=geno1)
            d <- diff(u[1:2])
            xlimits <- c(min(u) - d/4, max(u) + d/4)
        }
        else { # two markers
            u <- seq(along=geno2)
            d <- diff(u[1:2])
            xlimits <- c(min(u) - d/4, max(u) + d/4)
        }

        ## fix of x limits
        d <- 1
        xlimits <- c(1 - d/4, length(u) + d/4)

        if(is.null(mark2)) { # single marker
            if(missing(xlab)) xlab <- printname1
            if(missing(ylab)) ylab <- names(cross$pheno)[pheno.col]
            if(missing(col)) col <- "black"

            # plot the means
            plot(1:ngen1, means, main = main, xlab = xlab, ylab = ylab,
                 pch = 1, col = col[1], ylim = ylimits, xaxt = "n",
                 type = "b", xlim = xlimits)
            # confidence limits
            for(i in 1:ngen1) {
                if(!is.na(lo[i]) && !is.na(hi[i]))
                    lines(c(i, i), c(lo[i], hi[i]), pch = 3, col = col[1],
                          type = "b", lty = 3)
            }

            # X-axis ticks
            a <- par("usr")
            ystart <- a[3]
            yend <- ystart - diff(a[3:4]) * 0.02
            ytext <- ystart - diff(a[3:4]) * 0.05
            #      for(i in 1:ngen1) {
            #        lines(x = c(i, i), y = c(ystart, yend), xpd = TRUE)
            #        text(i, ytext, geno1[i], xpd = TRUE)
            #      }
            axis(side=1, at=1:ngen1, labels=geno1)
        }
        else { # two markers
            if(missing(xlab)) xlab <- printname2
            if(missing(ylab)) ylab <- names(cross$pheno)[pheno.col]
            # plot the first genotype of marker 1
            plot(1:ngen2, means[1, ], main = main, xlab = xlab,
                 ylab = ylab, pch = 1, col = int.color[1],
                 ylim = ylimits, xaxt = "n", type = "b", xlim = xlimits)
            # confidence limits
            for(i in 1:ngen2) {
                if(!is.na(lo[1, i]) && !is.na(hi[1, i]))
                    lines(c(i, i), c(lo[1, i], hi[1, i]), pch = 3,
                          col = int.color[1], type = "b", lty = 3)
            }
            for(j in 2:ngen1) { # for the rest of genotypes for Marker 1
                lines(1:ngen2, means[j, ], col = int.color[j], pch = 1,
                      type = "b")
                # confidence limits
                for(i in 1:ngen2) {
                    if(!is.na(lo[j, i]) && !is.na(hi[j, i]))
                        lines(c(i, i), c(lo[j, i], hi[j, i]), pch = 3,
                              col = int.color[j], type = "b", lty = 3)
                }
            }

            # draw X-axis ticks
            a <- par("usr")
            ystart <- a[3]
            yend <- ystart - diff(a[3:4]) * 0.02
            ytext <- ystart - diff(a[3:4]) * 0.05
            #      for(i in 1:ngen2) {
            #        lines(x = c(i, i), y = c(ystart, yend), xpd = TRUE)
            #        text(i, ytext, geno2[i], xpd = TRUE)
            #      }
            axis(side=1, at=1:ngen2, labels=geno2)

            # add legend
            if(add.legend) {
                col <- int.color[1:ngen1]
                # legend position
                x.leg <- a[1]*0.25+a[2]*0.75
                y.leg <- a[4] - diff(a[3:4]) * 0.05
                y.leg2 <- a[4] - diff(a[3:4]) * 0.03
                legend(x.leg, y.leg, geno1, lty = 1, pch = 1, col = col,
                       cex = 1, xjust = 0.5)
                if(missing(legend.lab)) legend.lab <- printname1
                text(x.leg, y.leg2, legend.lab)
            }
        }
    }

    return(invisible(result))
}

##############################################
# function to get genotype data for a marker
# given marker name
##############################################

effectplot.getmark <-
    function (cross, mname)
{
    # cross type
    type <- class(cross)[1]
    # return variables
    mark <- NULL
    gennames <- NULL

    # for pseudomarkers refered to as "1@10.5":
    #    - check that it is not a phenotype or marker name
    #    - otherwise convert to the usual name via find.pseudomarker
    pmalt.pattern <- "@-*[0-9]+" # alternate way to refer to a pseudomarker ("1@10.5")
    if(length(grep(pmalt.pattern, mname))>0 &&
       !(mname %in% names(cross$pheno) || mname %in% colnames(pull.geno(cross)))) {
        ss <- unlist(strsplit(mname, "@"))
        if(!(ss[1] %in% names(cross$geno)))
            stop("Don't understand the marker name ", mname)
        mname <- find.pseudomarker(cross, ss[1], as.numeric(ss[2]), "draws")
    }

    # determine marker type - it could be a marker, a pseudomarker or a phenotype
    mar.type <- "none"
    # regular expression pattern for a pseudomarker names
    pm.pattern <- "^c.*\\.loc.*$" # pseudomarker names will be like "c1.loc10"
    if(mname %in% names(cross$pheno)) { # this is a phenotype
        mar.type <- "pheno"
        idx.pos <- which(mname==names(cross$pheno))
    }
    else if(length(grep(pm.pattern, mname)) > 0) { # like "c1.loc10", this is a pseudomarker
        # note that the column names for draws is like "loc10",
        # so I need to take the part after "." in mname
        tmp <- unlist(strsplit(mname, "loc"))
        chr <- substr(tmp[1],2,nchar(tmp[1])-1) # this will be like 1 or "X"
        if( !(chr %in% names(cross$geno)) )
            stop("Couldn't find marker ", mname)
        mar.type <- "pm"
        chrtype <- class(cross$geno[[chr]])
        pm.name <- paste("loc", tmp[2],sep="") # this will be like loc10
        idx.pos <- which(pm.name==colnames(cross$geno[[chr]]$draws))
        if(length(idx.pos) == 0)
            stop("Couldn't find marker ", mname)
        else if(length(idx.pos)>1) # take the first one for multiple markers with the same name
            idx.pos <- idx.pos[1]
    }
    else { # this is a real marker name but it could be a observed or imputed
        for(i in 1:length(cross$geno)) {
            if(mname %in% colnames(cross$geno[[i]]$draws)) { # this is a pseudomarker
                mar.type <- "pm"
                chr <- i
                chrtype <- class(cross$geno[[chr]])
                idx.pos <- which(mname == colnames(cross$geno[[i]]$draws))
                break
            }
            else if(mname %in% colnames(cross$geno[[i]]$data)) { # this is a typed marker
                mar.type <- "marker"
                chr <- i
                chrtype <- class(cross$geno[[i]])
                idx.pos <- which(mname == colnames(cross$geno[[i]]$data))
                break
            }
        }
    }

    # if didn't find this marker
    if(mar.type == "none")
        stop("Marker ", mname, " not found")

    # get data from typed marker, pseudomarker or phenotype
    if(mar.type == "pheno") { # this is a phenotype
        mark <- cross$pheno[,idx.pos]
        # the phenotype need to be categorical
        if(length(unique(mark)) > 5) { # I'm using arbitrary number here
            stop("The input phenotype ", mname, " is not a categorical trait")
        }
        gennames <- sort(unique(mark))
    }
    else if(mar.type=="marker") { # this is a real marker
        mark <- cross$geno[[chr]]$data[, idx.pos]
        # if X chr and backcross or intercross, get sex/dir data + revise data
        if(chrtype == "X" && (type %in% c("bc","f2","bcsft"))) {
            sexpgm <- getsex(cross)
            mark <- as.numeric(reviseXdata(type, "full", sexpgm,
                                           geno = as.matrix(mark),
                                           cross.attr=attributes(cross)))
            gennames <- getgenonames(type, chrtype, "full", sexpgm, attributes(cross))
        }
    }

    else if(mar.type=="pm") { # this is a pseudomarker
        # get the imputed genotype data for this marker
        mark <- cross$geno[[chr]]$draws[,idx.pos,,drop=FALSE]

        # if X chr and backcross or intercross, get sex/dir data + revise data
        if(chrtype == "X" && (type %in% c("bc","f2","bcsft"))) {
            sexpgm <- getsex(cross)
            mark <- reviseXdata(type, "full", sexpgm, draws=mark,
                                cross.attr=attributes(cross))[,1,]
            gennames <- getgenonames(type, chrtype, "full", sexpgm, attributes(cross))
        }
        else mark <- mark[,1,]
    }

    else  # none of the above
        stop("Couldn't find marker ", mname)

    # make mark a matrix if it's not one
    if(class(mark) != "matrix")
        mark <- matrix(mark, ncol=1)

    # return
    ret <- list(mark=mark, gennames=gennames)
    attr(ret, "mname") <- mname
    ret
}

##############################################
# function to calculate the means and ses
# if ndraws is 1, it's easy
# if ndraws > 1 (has pseudomarker),
# loop thru the draws
##############################################
effectplot.calmeanse <-
    function(pheno, mark1, mark2, geno1, geno2, ndraws, var.flag=c("pooled","group"))
{
    # local variables
    nind <- length(pheno)
    # method to calculate variances for estimated QTL effects
    var.flag <- match.arg(var.flag)

    result <- NULL
    nind <- sum(!is.na(pheno)) # number of individuals

    if(is.null(mark2)) { # if mark2 is missing
        if(ndraws > 1) { # more than one draws
            mark1.level <- seq(along=geno1) # level for mark1
            # init
            means.all <- matrix(NA, nrow=ndraws, ncol=length(mark1.level))
            colnames(means.all) <- mark1.level
            vars.all <- matrix(NA, nrow=ndraws, ncol=length(mark1.level))
            colnames(vars.all) <- mark1.level
            weight <- rep(0, ndraws) # weight for draws
            # loop thru draws
            for(i in 1:ndraws) {
                mark1.tmp <- mark1[,i] # data for current draw
                # fit a regression - this is used to calculate the weights
                mark1.factor <- factor(mark1.tmp, mark1.level)
                lm.tmp <- lm(pheno~mark1.factor-1)
                rss <- sum(lm.tmp$residuals^2)
                # compute the weight
                weight[i] <- (-nind/2)*log(rss)
                # group means
                means.tmp <- tapply(pheno, mark1.tmp, mean, na.rm = TRUE)
                # calculate group means and variances
                if(var.flag=="group") { # use variance in each group
                    vars.tmp <- tapply(pheno, mark1.tmp, function(a) var(a,na.rm = TRUE)/length(a))
                }
                else { # use pooled variance
                    vars.tmp <- tapply(mark1.tmp, mark1.tmp, function(a) rss/nind/length(a))
                }
                # note that there could be missing categories in draws
                means.all[i, names(means.tmp)] <- means.tmp
                vars.all[i, names(vars.tmp)] <- vars.tmp
            }
            # average across draws - for vars, it should be
            # mean of variance plus variance of means
            weight <- exp(weight-max(weight))
            means <- apply(means.all, 2, function(a) weighted.mean(a,weight,na.rm=TRUE))
            meanvar <- apply(vars.all, 2, function(a) weighted.mean(a,weight,na.rm=TRUE)) # mean of vars
            varmean <- apply(means.all, 2, function(a) weighted.mean((a-mean(a,na.rm=TRUE))^2,weight,na.rm=TRUE)) # var of means
            measured <- apply(means.all, 2, function(a) sum(!is.na(a)))
            means[measured==0] <- meanvar[measured==0] <- varmean[measured==0] <- NA
            # standard error
            ses <- sqrt(meanvar+varmean)
        }
        else { # ndraws is 1
            u <- sort(unique(mark1))
            if(any(!(u %in% seq(along=geno1)))) {
                newmark1 <- mark1
                for(i in seq(along=u))
                    newmark1[mark1==u[i]] <- i
                mark1 <- newmark1
            }
            means <- tapply(pheno, mark1, mean, na.rm = TRUE)

            if(var.flag == "group") { # use group variance
                ses <- tapply(pheno, mark1, function(a) sd(a, na.rm = TRUE)/sqrt(sum(!is.na(a))))
            }
            else { # use pooled variance
                mark1.factor <- factor(mark1, seq(along=geno1))
                lm.tmp <- lm(pheno~mark1.factor-1)
                rss <- sum(lm.tmp$residuals^2)
                ses <- tapply(mark1, mark1, function(a) sqrt(rss/nind/length(a)))
            }
        }
    }

    else { # with mark2
        if(ndraws > 1) {
            mark1.level <- seq(along=geno1) # level for mark1
            mark2.level <- seq(along=geno2) # level for mark2

            u <- sort(unique(as.numeric(mark1)))
            if(any(!(u %in% seq(along=geno1)))) {
                newmark1 <- mark1
                for(i in seq(along=u))
                    for(j in 1:ncol(mark2))
                        newmark1[mark1[,j]==u[i],j] <- i
                mark1 <- newmark1
            }

            u <- sort(unique(as.numeric(mark2)))
            if(any(!(u %in% seq(along=geno2)))) {
                newmark2 <- mark2
                for(i in seq(along=u))
                    for(j in 1:ncol(mark2))
                        newmark2[mark2[,j]==u[i],j] <- i
                mark2 <- newmark2
            }

            # init
            means.all <- array(NA, c(length(mark1.level), length(mark2.level), ndraws))
            dimnames(means.all) <- list(mark1.level, mark2.level, NULL)
            vars.all <- array(NA, c(length(mark1.level), length(mark2.level), ndraws))
            dimnames(vars.all) <- list(mark1.level, mark2.level, NULL)
            weight <- rep(0, ndraws) # weight for draws
            # loop thru draws
            for(i in 1:ndraws) {
                mark1.tmp <- mark1[,i] # data for current draw
                mark2.tmp <- mark2[,i]
                # fit a regression - this is used to calculate the weights
                mark1.factor <- factor(mark1.tmp, mark1.level)
                mark2.factor <- factor(mark2.tmp, mark2.level)
                lm.tmp <- lm(pheno~mark1.factor+mark2.factor+1)
                rss <- sum(lm.tmp$residuals^2)
                # compute the weight
                weight[i] <- (-nind/2)*log(rss)
                # group means
                means.tmp <- tapply(pheno, list(mark1.tmp, mark2.tmp), mean, na.rm = TRUE)
                # calculate group means and variances
                if(var.flag=="group") { # use variance in each group
                    vars.tmp <- tapply(pheno, list(mark1.tmp,mark2.tmp),
                                       function(a) var(a,na.rm = TRUE)/length(a))
                }
                else { # use pooled variance
                    vars.tmp <- tapply(mark1.tmp, list(mark1.tmp,mark2.tmp),
                                       function(a) rss/nind/length(a))
                }
                # note that there could be missing categories in draws
                means.all[dimnames(means.tmp)[[1]], dimnames(means.tmp)[[2]],i] <- means.tmp
                vars.all[dimnames(vars.tmp)[[1]], dimnames(vars.tmp)[[2]], i] <- vars.tmp
            }
            # average across draws - for vars, it should be
            # mean of variance plus variance of means
            weight <- exp(weight-max(weight))
            means <- apply(means.all, c(1,2), function(a) weighted.mean(a,weight,na.rm=TRUE))
            meanvar <- apply(vars.all, c(1,2), function(a) weighted.mean(a,weight,na.rm=TRUE))
            varmean <- apply(means.all, c(1,2), function(a) weighted.mean((a-mean(a,na.rm=TRUE))^2,weight,na.rm=TRUE)) # var of means
            measured <- apply(means.all, c(1,2), function(a) sum(!is.na(a)))
            means[measured==0] <- meanvar[measured==0] <- varmean[measured==0] <- NA
            # standard error
            ses <- sqrt(meanvar+varmean)
        }
        else { # ndraws is 1
            u <- sort(unique(mark1))
            if(any(!(u %in% seq(along=geno1)))) {
                newmark1 <- mark1
                for(i in seq(along=u))
                    newmark1[mark1==u[i]] <- i
                mark1 <- newmark1
            }

            u <- sort(unique(mark2))
            if(any(!(u %in% seq(along=geno2)))) {
                newmark2 <- mark2
                for(i in seq(along=u))
                    newmark2[mark2==u[i]] <- i
                mark2 <- newmark2
            }
            mark1 <- factor(mark1, seq(along=geno1))
            mark2 <- factor(mark2, seq(along=geno2))

            means <- tapply(pheno, list(mark1, mark2), mean, na.rm = TRUE)
            if(var.flag=="group") { # use group variance
                ses <- tapply(pheno, list(mark1, mark2), function(a) sd(a, na.rm = TRUE)/sqrt(sum(!is.na(a))))
            }
            else {# use pooled variance
                lm.tmp <- lm(pheno~mark1+mark2-1)
                rss <- sum(lm.tmp$residuals^2)
                ses <- tapply(mark1, list(mark1, mark2), function(a) sqrt(rss/nind/length(a)))
            }
        }
    }

    # result
    result$Means <- means
    result$SEs <- ses
    result
}



# end of effectplot.R
