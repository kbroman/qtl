######################################################################
#
# simulate.R
#
# copyright (c) 2001-2012, Karl W Broman
# last modified May, 2012
# first written Apr, 2001
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
# Contains: sim.map, sim.cross, sim.cross.bc, sim.cross.f2,
#           sim.cross.4way, sim.bcg
#
######################################################################

######################################################################
#
# sim.map: simulate a genetic map
#
######################################################################

sim.map <-
    function(len=rep(100,20), n.mar=10, anchor.tel=TRUE, include.x=TRUE,
             sex.sp=FALSE, eq.spacing=FALSE)
{
    if(length(len)!=length(n.mar) && length(len)!=1 && length(n.mar)!=1)
        stop("Lengths of vectors len and n.mar do not conform.")

    # make vectors the same length
    if(length(len) == 1) len <- rep(len,length(n.mar))
    else if(length(n.mar) == 1) n.mar <- rep(n.mar,length(len))

    n.chr <- length(n.mar)

    map <- vector("list",n.chr)
    names(map) <- as.character(1:n.chr)
    if(include.x) names(map)[n.chr] <- "X"

    for(i in 1:n.chr) {
        if(anchor.tel) {
            if(n.mar[i] < 2) n.mar[i] <- 2
            map[[i]] <- c(0,len[i])
            if(n.mar[i] > 2) {
                if(!eq.spacing)
                    map[[i]] <- sort(c(map[[i]],runif(n.mar[i]-2,0,len[i])))
                else # equal spacing
                    map[[i]] <- seq(0,len[i],length=n.mar[i])
            }
        }
        else {
            if(!eq.spacing) {
                map[[i]] <- sort(runif(n.mar[i],0,len[i]))
                map[[i]] <- map[[i]] - min(map[[i]])
            }
            else {  # equal spacing
                map[[i]] <- seq(0,len[i],length=n.mar[i]+1)
                map[[i]] <- map[[i]][-1] - map[[i]][2]/2
            }
        }
        names(map[[i]]) <- paste("D", names(map)[i], "M", 1:n.mar[i], sep="")
        class(map[[i]]) <- "A"
    }

    if(sex.sp) {
        for(i in 1:n.chr) {
            if(eq.spacing) tempmap <- map[[i]]
            else {
                if(anchor.tel) {
                    if(n.mar[i] < 2) n.mar[i] <- 2
                    tempmap <- c(0,len[i])
                    if(n.mar[i] > 2)
                        tempmap <- sort(c(tempmap,runif(n.mar[i]-2,0,len[i])))
                }
                else {
                    tempmap <- sort(runif(n.mar[i],0,len[i]))
                    tempmap <- tempmap - min(tempmap)
                }
            }
            map[[i]] <- rbind(map[[i]],tempmap)
            dimnames(map[[i]]) <- list(NULL,paste("D", names(map)[i], "M", 1:n.mar[i], sep=""))
            class(map[[i]]) <- "A"

            if(include.x && i==n.chr)  # if X chromosome, force no recombination in male
                map[[i]][2,] <- rep(0,ncol(map[[i]]))
        }
    }

    if(include.x) class(map[[n.chr]]) <- "X"

    class(map) <- "map"
    map
}

######################################################################
#
# sim.cross: Simulate an experimental cross
#
# Note: These functions are a bit of a mess.  I was in the "get it to
#       work without worrying about efficiency" mode while writing it.
#       Sorry!
#
######################################################################

sim.cross <-
    function(map, model=NULL, n.ind=100,
             type=c("f2", "bc", "4way", "risib", "riself",
             "ri4sib", "ri4self", "ri8sib", "ri8self","bcsft"),
             error.prob=0, missing.prob=0, partial.missing.prob=0,
             keep.qtlgeno=TRUE, keep.errorind=TRUE, m=0, p=0,
             map.function=c("haldane","kosambi","c-f","morgan"),
             founderGeno, random.cross=TRUE, ...)
{
    type <- match.arg(type)
    map.function <- match.arg(map.function)

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }

    # 2-way RIL by sibmating or selfing
    if(type=="risib" || type=="riself") {
        if(type=="risib") type <- "sibmating"
        else type <- "selfing"
        cross <- sim.ril(map, n.ind, type, "2", m=m, p=p,
                         error.prob=error.prob, missing.prob=missing.prob)
        cross$cross <- NULL

        return(cross)
    }
    # 4- or 8-way RIL by sibmating or selfing
    if(type=="ri4sib" || type=="ri4self" || type=="ri8sib" || type=="ri8self") {
        if(substr(type, 4, nchar(type))=="self") crosstype <- "selfing"
        else crosstype <- "sibmating"
        n.str <- substr(type, 3, 3)
        cross <- sim.ril(map, n.ind, crosstype, n.str, m=m, p=p,
                         random.cross=random.cross,
                         error.prob=0,
                         missing.prob=missing.prob)

        rcross <- convertMWril(cross, founderGeno, error.prob=error.prob)
        for(i in names(cross$geno))
            if(!("truegeno" %in% names(rcross$geno[[i]])))
                rcross$geno[[i]]$truegeno <- cross$geno[[i]]$data

        # remove "un" from cross type
        class(rcross)[1] <- substr(class(cross)[1], 1, nchar(class(cross)[1])-2)

        fg <- t(founderGeno[[1]])
        if(length(founderGeno)>1)
            for(i in 2:length(founderGeno))
                fg <- cbind(fg, t(founderGeno[[i]]))
        colnames(fg) <- markernames(rcross)
        rcross$founderGeno <- fg
        return(rcross)
    }

    # sort the model matrix
    if(!is.null(model) && is.matrix(model))
        model <- model[order(model[,1],model[,2]),]

    if(type=="bc")
        cross <- sim.cross.bc(map,model,n.ind,error.prob,missing.prob,
                              keep.errorind,m,p,map.function)
    else if(type=="f2")
        cross <- sim.cross.f2(map,model,n.ind,error.prob,missing.prob,
                              partial.missing.prob,keep.errorind,
                              m,p,map.function)
    else if(type=="bcsft")
        cross <- sim.cross.bcsft(map,model,n.ind,error.prob,missing.prob,
                                 partial.missing.prob,keep.errorind,
                                 m,p,map.function, ...)
    else
        cross <- sim.cross.4way(map,model,n.ind,error.prob,missing.prob,
                                partial.missing.prob,keep.errorind,
                                m,p,map.function)

    # remove QTL genotypes from data and, if keep.qtlgeno=TRUE,
    #     place them in cross$qtlgeno
    qtlgeno <- NULL
    for(i in 1:nchr(cross)) {
        o <- grep("^QTL[0-9]+", colnames(cross$geno[[i]]$data))
        if(length(o) != 0) {
            qtlgeno <- cbind(qtlgeno, cross$geno[[i]]$data[,o,drop=FALSE])
            cross$geno[[i]]$data <- cross$geno[[i]]$data[,-o,drop=FALSE]
            if(is.matrix(cross$geno[[i]]$map))
                cross$geno[[i]]$map <- cross$geno[[i]]$map[,-o,drop=FALSE]
            else
                cross$geno[[i]]$map <- cross$geno[[i]]$map[-o]
        }
    }
    if(keep.qtlgeno) cross$qtlgeno <- qtlgeno

    # store genotype data as integers
    for(i in 1:nchr(cross))
        storage.mode(cross$geno[[i]]$data) <- "integer"

    if(is.null(names(cross$geno)))
        names(cross$geno) <- 1:length(cross$geno)

    cross
}


######################################################################
#
# sim.cross.bc
#
######################################################################

sim.cross.bc <-
    function(map,model,n.ind,error.prob,missing.prob,
             keep.errorind,m,p,map.function)
{
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    if(any(sapply(map,is.matrix)))
        stop("Map must not be sex-specific.")

    chr.type <- sapply(map, function(a) ifelse(class(a)=="X","X","A"))

    n.chr <- length(map)

    if(is.null(model)) n.qtl <- 0
    else {
        if(!((!is.matrix(model) && length(model) == 3) ||
             (is.matrix(model) && ncol(model) == 3)))
            stop("Model must be a matrix with 3 columns (chr, pos and effect).")
        if(!is.matrix(model)) model <- rbind(model)
        n.qtl <- nrow(model)
        if(any(model[,1] < 0 | model[,1] > n.chr))
            stop("Chromosome indicators in model matrix out of range.")
        model[,2] <- model[,2]+1e-14 # so QTL not on top of marker
    }

    # if any QTLs, place qtls on map
    if(n.qtl > 0) {
        for(i in 1:n.qtl) {
            temp <- map[[model[i,1]]]
            if(model[i,2] < min(temp)) {
                temp <- c(model[i,2],temp)
                names(temp)[1] <- paste("QTL",i,sep="")
            }
            else if(model[i,2] > max(temp)) {
                temp <- c(temp,model[i,2])
                names(temp)[length(temp)] <- paste("QTL",i,sep="")
            }
            else {
                j <- max((seq(along=temp))[temp < model[i,2]])
                temp <- c(temp[1:j],model[i,2],temp[(j+1):length(temp)])
                names(temp)[j+1] <- paste("QTL",i,sep="")
            }
            map[[model[i,1]]] <- temp
        }
    }

    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    n.mar <- sapply(map,length)
    mar.names <- lapply(map,names)

    for(i in 1:n.chr) {
        # simulate genotype data
        thedata <- sim.bcg(n.ind, map[[i]], m, p, map.function)
        dimnames(thedata) <- list(NULL,mar.names[[i]])

        geno[[i]] <- list(data = thedata, map = map[[i]])
        class(geno[[i]]) <- chr.type[i]
        class(geno[[i]]$map) <- NULL

    } # end loop over chromosomes

    # simulate phenotypes
    pheno <- rnorm(n.ind,0,1)

    if(n.qtl > 0) {
        # find QTL positions in genotype data
        QTL.chr <- QTL.loc <- NULL
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0) {
                QTL.chr <- c(QTL.chr,rep(i,length(o)))
                QTL.loc <- c(QTL.loc,o)
            }
        }

        # incorporate QTL effects
        for(i in 1:n.qtl) {
            QTL.geno <- geno[[QTL.chr[i]]]$data[,QTL.loc[i]]
            pheno[QTL.geno==2] <- pheno[QTL.geno==2] + model[i,3]
        }

    } # end simulate phenotype

    n.mar <- sapply(geno, function(a) length(a$map))

    # add errors
    if(error.prob > 0) {
        for(i in 1:n.chr) {
            a <- sample(0:1,n.mar[i]*n.ind,replace=TRUE,
                        prob=c(1-error.prob,error.prob))
            geno[[i]]$data[a == 1] <- 3 - geno[[i]]$data[a == 1]
            if(keep.errorind) {
                errors <- matrix(0,n.ind,n.mar[i])
                errors[a==1] <- 1
                colnames(errors) <- colnames(geno[[i]]$data)
                geno[[i]]$errors <- errors
            }
        }
    }

    # add missing
    if(missing.prob > 0) {
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0)
                x <- geno[[i]]$data[,o]
            geno[[i]]$data[sample(c(TRUE,FALSE),n.mar[i]*n.ind,replace=TRUE,
                                  prob=c(missing.prob,1-missing.prob))] <- NA
            if(length(o)>0)
                geno[[i]]$data[,o] <- x
        }
    }

    pheno <- data.frame(phenotype=pheno, stringsAsFactors=TRUE)

    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c("bc","cross")

    cross
}

######################################################################
#
# sim.cross.f2
#
######################################################################

sim.cross.f2 <-
    function(map,model,n.ind,error.prob,missing.prob,partial.missing.prob,
             keep.errorind,m,p,map.function)
{
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    if(any(sapply(map,is.matrix)))
        stop("Map must not be sex-specific.")

    # chromosome types
    chr.type <- sapply(map,function(a) ifelse(class(a)=="X", "X", "A"))

    n.chr <- length(map)
    if(is.null(model)) n.qtl <- 0
    else {
        if(!((!is.matrix(model) && length(model) == 4) ||
             (is.matrix(model) && ncol(model) == 4))) {
            stop("Model must be a matrix with 4 columns (chr, pos and effects).")
        }
        if(!is.matrix(model)) model <- rbind(model)
        n.qtl <- nrow(model)
        if(any(model[,1] < 0 | model[,1] > n.chr))
            stop("Chromosome indicators in model matrix out of range.")
        model[,2] <- model[,2]+1e-14 # so QTL not on top of marker
    }

    # if any QTLs, place qtls on map
    if(n.qtl > 0) {
        for(i in 1:n.qtl) {
            temp <- map[[model[i,1]]]
            if(model[i,2] < min(temp)) {
                temp <- c(model[i,2],temp)
                names(temp)[1] <- paste("QTL",i,sep="")
            }
            else if(model[i,2] > max(temp)) {
                temp <- c(temp,model[i,2])
                names(temp)[length(temp)] <- paste("QTL",i,sep="")
            }
            else {
                j <- max((seq(along=temp))[temp < model[i,2]])
                temp <- c(temp[1:j],model[i,2],temp[(j+1):length(temp)])
                names(temp)[j+1] <- paste("QTL",i,sep="")
            }
            map[[model[i,1]]] <- temp
        }
    }

    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    n.mar <- sapply(map,length)
    mar.names <- lapply(map,names)

    for(i in 1:n.chr) {

        # simulate genotype data
        thedata <- sim.bcg(n.ind, map[[i]], m, p, map.function)
        dimnames(thedata) <- list(NULL,mar.names[[i]])

        if(chr.type[i] != "X")
            thedata <- thedata + sim.bcg(n.ind, map[[i]], m, p, map.function) - 1

        geno[[i]] <- list(data = thedata, map = map[[i]])
        class(geno[[i]]) <- chr.type[i]
        class(geno[[i]]$map) <- NULL

    } # end loop over chromosomes

    # simulate phenotypes
    pheno <- rnorm(n.ind,0,1)

    if(n.qtl > 0) {
        # find QTL positions in genotype data
        QTL.chr <- QTL.loc <- NULL
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0) {
                QTL.chr <- c(QTL.chr,rep(i,length(o)))
                QTL.loc <- c(QTL.loc,o)
            }
        }

        # incorporate QTL effects
        for(i in 1:n.qtl) {
            QTL.geno <- geno[[QTL.chr[i]]]$data[,QTL.loc[i]]
            pheno[QTL.geno==1] <- pheno[QTL.geno==1] - model[i,3]
            pheno[QTL.geno==2] <- pheno[QTL.geno==2] + model[i,4]
            pheno[QTL.geno==3] <- pheno[QTL.geno==3] + model[i,3]
        }

    } # end simulate phenotype

    n.mar <- sapply(geno, function(a) length(a$map))

    # add errors
    if(error.prob > 0) {
        for(i in 1:n.chr) {
            if(chr.type[i]=="X") {
                a <- sample(0:1,n.mar[i]*n.ind,replace=TRUE,
                            prob=c(1-error.prob,error.prob))
                geno[[i]]$data[a == 1] <- 3 - geno[[i]]$data[a == 1]
            }
            else {
                a <- sample(0:2,n.mar[i]*n.ind,replace=TRUE,
                            prob=c(1-error.prob,error.prob/2,error.prob/2))
                if(any(a>0 & geno[[i]]$data==1))
                    geno[[i]]$data[a>0 & geno[[i]]$data==1] <-
                        (geno[[i]]$data+a)[a>0 & geno[[i]]$data==1]
                if(any(a>0 & geno[[i]]$data==2)) {
                    geno[[i]]$data[a>0 & geno[[i]]$data==2] <-
                        (geno[[i]]$data+a)[a>0 & geno[[i]]$data==2]
                    geno[[i]]$data[geno[[i]]$data>3] <- 1
                }
                if(any(a>0 & geno[[i]]$data==3))
                    geno[[i]]$data[a>0 & geno[[i]]$data==3] <-
                        (geno[[i]]$data-a)[a>0 & geno[[i]]$data==3]
            }

            if(keep.errorind) {
                errors <- matrix(0,n.ind,n.mar[i])
                errors[a>0] <- 1
                colnames(errors) <- colnames(geno[[i]]$data)
                geno[[i]]$errors <- errors
            }

        } # end loop over chromosomes
    } # end simulate genotyping errors

    # add partial missing
    if(partial.missing.prob > 0) {
        for(i in 1:n.chr) {
            if(chr.type[i] != "X") {
                o <- sample(c(TRUE,FALSE),n.mar[i],replace=TRUE,
                            prob=c(partial.missing.prob,1-partial.missing.prob))
                if(any(o)) {
                    o2 <- grep("^QTL[0-9]+",mar.names[[i]])
                    if(length(o2)>0)
                        x <- geno[[i]]$data[,o2]
                    m <- (1:n.mar[i])[o]
                    for(j in m) {
                        if(runif(1) < 0.5)
                            geno[[i]]$data[geno[[i]]$data[,j]==1 | geno[[i]]$data[,j]==2,j] <- 4
                        else
                            geno[[i]]$data[geno[[i]]$data[,j]==3 | geno[[i]]$data[,j]==2,j] <- 5
                    }
                    if(length(o2)>0)
                        geno[[i]]$data[,o2] <- x
                }
            }

        } # end loop over chromosomes
    } # end simulate partially missing data

    # add missing
    if(missing.prob > 0) {
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0)
                x <- geno[[i]]$data[,o]
            geno[[i]]$data[sample(c(TRUE,FALSE),n.mar[i]*n.ind,replace=TRUE,
                                  prob=c(missing.prob,1-missing.prob))] <- NA
            if(length(o)>0)
                geno[[i]]$data[,o] <- x
        }
    }

    pheno <- data.frame(phenotype=pheno, stringsAsFactors=TRUE)

    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c("f2","cross")

    cross
}

######################################################################
#
# sim.cross.4way
#
######################################################################

sim.cross.4way <-
    function(map,model,n.ind,error.prob,missing.prob,partial.missing.prob,
             keep.errorind,m,p,map.function)
{
    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    if(!all(sapply(map,is.matrix)))
        stop("Map must be sex-specific.")

    n.chr <- length(map)
    if(is.null(model)) n.qtl <- 0
    else {
        if(!((!is.matrix(model) && length(model) == 5) ||
             (is.matrix(model) && ncol(model) == 5))) {
            stop("Model must be a matrix with 5 columns (chr, pos and effects).")
        }
        if(!is.matrix(model)) model <- rbind(model)
        n.qtl <- nrow(model)
        if(any(model[,1] < 0 | model[,1] > n.chr))
            stop("Chromosome indicators in model matrix out of range.")
        model[,2] <- model[,2]+1e-14 # so QTL not on top of marker
    }

    chr.type <- sapply(map,function(a) ifelse(class(a)=="X", "X", "A"))

    # if any QTLs, place qtls on map
    if(n.qtl > 0) {
        for(i in 1:n.qtl) {
            temp <- map[[model[i,1]]]
            temp1 <- temp[1,]
            temp2 <- temp[2,]
            qtlloc <- model[i,2]

            if(qtlloc < min(temp1)) {
                temp1 <- c(qtlloc,temp1)
                temp2 <- min(temp2) - (min(temp1)-qtlloc)/diff(range(temp1))*diff(range(temp2))
                temp1 <- temp1-min(temp1)
                temp2 <- temp2-min(temp2)
                n <- c(paste("QTL",i,sep=""),colnames(temp))
            }
            else if(qtlloc > max(temp1)) {
                temp1 <- c(temp1,qtlloc)
                temp2 <- (qtlloc-max(temp1))/diff(range(temp1))*diff(range(temp2))+max(temp2)
                n <- c(colnames(temp),paste("QTL",i,sep=""))
            }
            else {
                temp1 <- c(temp1,qtlloc)
                o <- order(temp1)
                wh <- (seq(along=temp1))[order(temp1)==length(temp1)]
                temp2 <- c(temp2[1:(wh-1)],NA,temp2[-(1:(wh-1))])
                temp2[wh] <- temp2[wh-1] + (temp1[wh]-temp1[wh-1])/(temp1[wh+1]-temp1[wh-1]) *
                    (temp2[wh+1]-temp2[wh-1])
                temp1 <- sort(temp1)
                n <- c(colnames(temp),paste("QTL",i,sep=""))[o]
            }
            map[[model[i,1]]] <- rbind(temp1,temp2)
            dimnames(map[[model[i,1]]]) <- list(NULL, n)
        }
    }

    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    n.mar <- sapply(map,ncol)
    mar.names <- lapply(map,function(a) colnames(a))

    for(i in 1:n.chr) {

        # simulate sex
        sex <- NULL
        if(chr.type[i]=="X")
            sex <- rep(0,n.ind)

        # simulate genotype data
        thedata <- sim.bcg(n.ind, map[[i]], m, p, map.function)
        dimnames(thedata) <- list(NULL,mar.names[[i]])

        if(chr.type[i] != "X")
            thedata <- thedata + 2*sim.bcg(n.ind, map[[i]][2:1,], m, p, map.function) - 2

        dimnames(thedata) <- list(NULL,mar.names[[i]])

        geno[[i]] <- list(data = thedata, map = map[[i]])

        class(geno[[i]]) <- chr.type[i]
        class(geno[[i]]$map) <- NULL
    } # end loop over chromosomes

    # simulate phenotypes
    pheno <- rnorm(n.ind,0,1)

    if(n.qtl > 0) {
        # find QTL positions
        QTL.chr <- QTL.loc <- NULL
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0) {
                QTL.chr <- c(QTL.chr,rep(i,length(o)))
                QTL.loc <- c(QTL.loc,o)
            }
        }

        # incorporate QTL effects
        for(i in 1:n.qtl) {
            QTL.geno <- geno[[QTL.chr[i]]]$data[,QTL.loc[i]]
            pheno[QTL.geno==1] <- pheno[QTL.geno==1] + model[i,3]
            pheno[QTL.geno==2] <- pheno[QTL.geno==2] + model[i,4]
            pheno[QTL.geno==3] <- pheno[QTL.geno==3] + model[i,5]
        }

    } # end simulate phenotype

    n.mar <- sapply(geno, function(a) ncol(a$map))

    # add errors
    if(error.prob > 0) {
        for(i in 1:n.chr) {
            if(chr.type[i] != "X") { # 4-way cross; autosomal
                a <- sample(0:3,n.mar[i]*n.ind,replace=TRUE,
                            prob=c(1-error.prob,rep(error.prob/3,3)))
                if(any(a>0 & geno[[i]]$data==1))
                    geno[[i]]$data[a>0 & geno[[i]]$data==1] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==1] + a[a>0 & geno[[i]]$data==1]
                if(any(a>0 & geno[[i]]$data==2))
                    geno[[i]]$data[a>0 & geno[[i]]$data==2] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==2] + c(-1,1,2)[a[a>0 & geno[[i]]$data==2]]
                if(any(a>0 & geno[[i]]$data==3))
                    geno[[i]]$data[a>0 & geno[[i]]$data==3] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==3] + c(-2,-1,1)[a[a>0 & geno[[i]]$data==3]]
                if(any(a>0 & geno[[i]]$data==4))
                    geno[[i]]$data[a>0 & geno[[i]]$data==4] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==4] - a[a>0 & geno[[i]]$data==4]
            }
            else {
                a <- sample(0:1,n.mar[i]*n.ind,replace=TRUE,
                            prob=c(1-error.prob,error.prob))
                if(any(a>0 & geno[[i]]$data==1))
                    geno[[i]]$data[a>0 & geno[[i]]$data==1] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==1] + 1
                if(any(a>0 & geno[[i]]$data==2))
                    geno[[i]]$data[a>0 & geno[[i]]$data==2] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==2] - 1
                if(any(a>0 & geno[[i]]$data==3))
                    geno[[i]]$data[a>0 & geno[[i]]$data==3] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==3] + 1
                if(any(a>0 & geno[[i]]$data==4))
                    geno[[i]]$data[a>0 & geno[[i]]$data==4] <-
                        geno[[i]]$data[a>0 & geno[[i]]$data==4] - 1
            }

            if(keep.errorind) {
                errors <- matrix(0,n.ind,n.mar[i])
                errors[a>0] <- 1
                colnames(errors) <- colnames(geno[[i]]$data)
                geno[[i]]$errors <- errors
            }

        } # end loop over chromosomes
    } # end simulate genotyping errors

    # add partial missing
    if(partial.missing.prob > 0) {
        for(i in 1:n.chr) {
            if(chr.type[i] != "X") {
                o <- sample(c(TRUE,FALSE),n.mar[i],replace=TRUE,
                            prob=c(partial.missing.prob,1-partial.missing.prob))

                if(any(o)) {
                    o2 <- grep("^QTL[0-9]+",mar.names[[i]])
                    if(length(o2)>0)
                        x <- geno[[i]]$data[,o2]
                    m <- (1:n.mar[i])[o]
                    for(j in m) {
                        a <- sample(1:4,1)
                        if(a==1) { # AB:AA marker
                            geno[[i]]$data[geno[[i]]$data[,j]==1 | geno[[i]]$data[,j]==3,j] <- 5
                            geno[[i]]$data[geno[[i]]$data[,j]==2 | geno[[i]]$data[,j]==4,j] <- 6
                        }
                        else if(a==2) { # AA:AB marker
                            geno[[i]]$data[geno[[i]]$data[,j]==1 | geno[[i]]$data[,j]==2,j] <- 7
                            geno[[i]]$data[geno[[i]]$data[,j]==3 | geno[[i]]$data[,j]==4,j] <- 8
                        }
                        else if(a==3)  # AB:AB marker
                            geno[[i]]$data[geno[[i]]$data[,j]==2 | geno[[i]]$data[,j]==3,j] <- 10
                        else  # AB:BA marker
                            geno[[i]]$data[geno[[i]]$data[,j]==1 | geno[[i]]$data[,j]==4,j] <- 9
                    }
                    if(length(o2) > 0)
                        geno[[i]]$data[,o2] <- x
                }
            }

        } # end loop over chromosomes
    } # end simulate partially missing data

    # add missing
    if(missing.prob > 0) {
        for(i in 1:n.chr) {
            o <- grep("^QTL[0-9]+",mar.names[[i]])
            if(length(o)>0)
                x <- geno[[i]]$data[,o]
            geno[[i]]$data[sample(c(TRUE,FALSE),n.mar[i]*n.ind,replace=TRUE,
                                  prob=c(missing.prob,1-missing.prob))] <- NA
            if(length(o)>0)
                geno[[i]]$data[,o] <- x
        }
    }

    if(!is.null(sex)) {
        pheno <- cbind(pheno,sex)
        dimnames(pheno) <- list(NULL, c("phenotype", "sex"))
    }
    else {
        pheno <- cbind(pheno)
        dimnames(pheno) <- list(NULL, "phenotype")
    }

    pheno <- as.data.frame(pheno, stringsAsFactors=TRUE)
    cross <- list(geno=geno,pheno=pheno)
    class(cross) <- c("4way","cross")

    cross
}


######################################################################
# sim.bcg
#
# simulate backcross genotype data for a single chromosome;
# output is a matrix of 1's and 0's
######################################################################

sim.bcg <-
    function(n.ind, map, m, p,
             map.function=c("haldane","kosambi","c-f","morgan"))
{
    map.function <- match.arg(map.function)

    if(map.function=="kosambi") mf <- mf.k
    else if(map.function=="c-f") mf <- mf.cf
    else if(map.function=="morgan") mf <- mf.m
    else mf <- mf.h

    if(m < 0 || p < 0 || p > 1)
        stop("Must have m >= 0 and 0 <= p <= 1")

    if(is.matrix(map)) map <- map[1,]
    map <- map-map[1]
    n.mar <- length(map)

    if(m==0 || p==1) { # no interference
        g <- .C("R_sim_bc_ni",
                as.integer(n.mar),
                as.integer(n.ind),
                as.double(mf(diff(map))),
                g=as.integer(rep(0, n.mar*n.ind)),
                PACKAGE="qtl")$g
    }
    else {
        g <- .C("R_sim_bc",
                as.integer(n.mar),
                as.integer(n.ind),
                as.double(map),
                as.integer(m),
                as.double(p),
                g=as.integer(rep(0, n.mar*n.ind)),
                PACKAGE="qtl")$g
    }
    matrix(g, ncol=n.mar)
}



# end of simulate.R
