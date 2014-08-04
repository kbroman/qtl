######################################################################
#
# read.cross.mq.R
#
# copyright (c) 2014, INRA (author: Timothée Flutre)
# last modified June, 2014
# first written May, 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License,
# version 3, as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose. See the GNU
# General Public License, version 3, for more details.
#
# A copy of the GNU General Public License, version 3, is available
# at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/qtl package
# Contains: read.cross.mq, read.cross.mq.loc, read.cross.mq.map,
#           read.cross.mq.qua, mq.rmv.comment
# [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.mq: read data from an experimental cross in MapQTL (and
# JoinMap) format.
#
# We need three files: a "loc" file containing the genotype data, a
# "map" file containing the linkage group assignments and map
# positions, and a "qua" file containing the phenotypes.
#
# File formats are described in the MapQTL manual available online at
# http://www.kyazma.nl/docs/MQ6Manual.pdf
# For the loc-file, each marker should be on a single line.
# Only 4-way crosses are supported ("CP" type in MapQTL/JoinMap).
#
######################################################################

read.cross.mq <-
function(dir, locfile, mapfile, quafile, estimate.map=TRUE)
{
    if(! missing(dir) && dir != "") {
        locfile <- file.path(dir, locfile)
        mapfile <- file.path(dir, mapfile)
        quafile <- file.path(dir, quafile)
    }

    loc <- read.cross.mq.loc(locfile)
    map <- read.cross.mq.map(mapfile)
    pheno <- read.cross.mq.qua(quafile)

    type <- loc$pop.type # only "4way" for the moment
    n.ind <- nrow(loc$genotypes)
    n.mar <- ncol(loc$genotypes)
    n.phe <- ncol(pheno)
    if(nrow(pheno) < n.ind){
        msg <- paste("qua-file should have at least the same number of",
                     " individuals than loc-file")
        stop(msg, call.=FALSE)
    }
    if(nrow(pheno) > n.ind){
        pheno <- pheno[-((n.ind+1):nrow(pheno)),]
    }

    cat(" --Read the following data:\n")
    cat("\tNumber of individuals: ", n.ind, "\n")
    cat("\tNumber of markers: ", n.mar, "\n")
    cat("\tNumber of phenotypes: ", n.phe, "\n")

    geno <- list()
    for(chr in levels(map$chr)){
        geno[[chr]] <- list(data=loc$genotypes[,map$marker[map$chr == chr]],
                            map=map$pos[map$chr == chr])
        names(geno[[chr]]$map) <- map$marker[map$chr == chr]
        if(type == "4way")
            geno[[chr]]$map <- rbind(geno[[chr]]$map,
                                     map$pos[map$chr == chr])
        if(chr %in% c("x", "X")){
            class(geno[[chr]]) <- "X"
        } else
            class(geno[[chr]]) <- "A"
    }

    cross <- list(geno=geno, pheno=pheno)
    class(cross) <- c(type, "cross")
    list(cross, estimate.map)
}

mq.rmv.comment <-
function(x, symbol=";")
{
    x.wo.comment <- x
    comments <- regexpr(pattern=symbol, text=x)
    if(comments[1] != -1)
        x.wo.comment <- substr(x=x, start=1, stop=comments[1]-1)
    return(x.wo.comment)
}

## each marker should be on a single line
## only 4-way crosses (CP type) are handled
read.cross.mq.loc <-
function(locfile){
    pop.name <- NULL
    pop.type <- NULL
    nb.loci <- NULL
    nb.inds <- NULL
    seg <- NULL
    phase <- NULL
    classif <- NULL
    genotypes <- NULL

    lines <- readLines(locfile)

    locus.id <- 1
    ind.id <- 1
    for(line.id in 1:length(lines)){

        ## skip empty lines
        tmp <- mq.rmv.comment(x=lines[line.id], symbol=";")
        tmp <- gsub(pattern=" + ", replacement=" ", x=tmp)
        if(tmp == "" || tmp == " ")
            next

        ## extract the population name
        if(grepl(pattern="name", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            pop.name <- strsplit(x=tmp, split="=")[[1]][2]
            next
        }

        ## extract the population type
        if(grepl(pattern="popt", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            pop.type <- strsplit(x=tmp, split="=")[[1]][2]
            pop.types <- c("BC1","F2","RIx","DH","DH1","DH2","HAP","HAP1","CP","BCpxFy",
                           "IMxFy")
            if(! pop.type %in% pop.types){
                msg <- paste("unknown population type", pop.type)
                stop(msg, call.=FALSE)
            }
            if(pop.type == "CP"){
                pop.type <- "4way"
            } else{
                msg <- paste("population type", pop.type,
                             "is not supported (yet)")
                stop(msg, call.=FALSE)
            }
            next
        }

        ## extract the number of loci
        if(grepl(pattern="nloc", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            nb.loci <- as.numeric(strsplit(x=tmp, split="=")[[1]][2])
            seg <- rep(NA, nb.loci)
            phase <- rep(NA, nb.loci)
            classif <- rep(NA, nb.loci)
            next
        }

        ## extract the number of individuals
        if(grepl(pattern="nind", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            nb.inds <- as.numeric(strsplit(x=tmp, split="=")[[1]][2])
            genotypes <- matrix(NA, nrow=nb.loci, ncol=nb.inds)
            rownames(genotypes) <- paste("loc", 1:nb.loci, sep=".")
            ## colnames(genotypes) <- paste("ind", 1:nb.inds, sep=".")
            next
        }

        tokens <- strsplit(x=tmp, split=" |\t")[[1]]
        if(length(tokens) > nb.inds + 4){
            msg <- paste("line", line.id, "should have a maximum of",
                         nb.inds + 4,
                         "columns separated by spaces or tabs")
            stop(msg, call.=FALSE)
        }
        rownames(genotypes)[locus.id] <- tokens[1]
        if(length(tokens) > nb.inds + 1){
            for(i in 2:(length(tokens)-nb.inds)){
                if(length(grep(pattern="\\{", x=tokens[i])) > 0){
                    phase[locus.id] <- tokens[i]
                } else if(length(grep(pattern="\\(", x=tokens[i])) > 0){
                    classif[locus.id] <- tokens[i]
                } else
                    seg[locus.id] <- tokens[i] # only for "4way" type
            }
        }
        genotypes[locus.id,] <- tokens[(length(tokens)-nb.inds+1):length(tokens)]
        locus.id <- locus.id + 1
    }

    if(locus.id > nb.loci + 1){
        msg <- paste("there seems to be more loci (", locus.id-1,
                             ") than indicated in the header (", nb.loci, ")")
        stop(msg, call.=FALSE)
    }

    genotypes <- t(genotypes) # individuals in rows, markers in columns

    ## convert all missing data to NA
    for(i in 1:length(genotypes)){
        if(grepl(pattern="\\.", x=genotypes[i]))
            genotypes[i] <- gsub(pattern="\\.", replacement="-",
                                 x=genotypes[i])
        if(grepl(pattern="u", x=genotypes[i]))
            genotypes[i] <- gsub(pattern="u", replacement="-",
                                 x=genotypes[i])
    }
    genotypes[which(genotypes == "--")] <- NA

    ## convert segregation types and genotypes to new MapQtl format
    convert.seg <- FALSE
    new.seg.types <- c("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>")
    for(seg.type in unique(seg))
        if(! seg.type %in% new.seg.types){ # e.g. <abxac>
            convert.seg <- TRUE
            break
        }
    if(convert.seg){
        if(pop.type != "4way"){
            msg <- paste("can't convert genotypes from old to new format",
                         "for 4way cross (yet)")
            stop(msg, call.=FALSE)
        }
        for(locus.id in 1:nb.loci){
            if(seg[locus.id] == "<abxcd>"){
                next
            } else if(seg[locus.id] == "<abxac>"){
                seg[locus.id] <- "<efxeg>"
                tmp <- which(genotypes[,locus.id] == "aa")
                genotypes[tmp, locus.id] <- "ee"
                tmp <- which(genotypes[,locus.id] == "ab")
                genotypes[tmp, locus.id] <- "ef"
                tmp <- which(genotypes[,locus.id] == "ac")
                genotypes[tmp, locus.id] <- "eg"
                tmp <- which(genotypes[,locus.id] == "bc")
                genotypes[tmp, locus.id] <- "fg"
            } else if(seg[locus.id] == "<abxab>"){
                seg[locus.id] <- "<hkxhk>"
                tmp <- which(genotypes[,locus.id] == "aa")
                genotypes[tmp, locus.id] <- "hh"
                tmp <- which(genotypes[,locus.id] == "ab")
                genotypes[tmp, locus.id] <- "hk"
                tmp <- which(genotypes[,locus.id] == "bb")
                genotypes[tmp, locus.id] <- "kk"
                tmp <- which(genotypes[,locus.id] == "a-")
                genotypes[tmp, locus.id] <- "h-"
                tmp <- which(genotypes[,locus.id] == "b-")
                genotypes[tmp, locus.id] <- "k-"
            } else if(seg[locus.id] == "<abxaa>"){
                seg[locus.id] <- "<lmxll>"
                tmp <- which(genotypes[,locus.id] == "aa")
                genotypes[tmp, locus.id] <- "ll"
                tmp <- which(genotypes[,locus.id] == "ab")
                genotypes[tmp, locus.id] <- "lm"
            } else if(seg[locus.id] == "<aaxab>"){
                seg[locus.id] <- "<nnxnp>"
                tmp <- which(genotypes[,locus.id] == "aa")
                genotypes[tmp, locus.id] <- "nn"
                tmp <- which(genotypes[,locus.id] == "ab")
                genotypes[tmp, locus.id] <- "np"
            } else{
                msg <- paste("unrecognized segregation type", seg[locus.id],
                             "at locus", locus.id)
                stop(msg, call.=FALSE)
            }
        }
    }

    ## check genotypes
    proper.genotypes <- c(NA, "ac", "ca", "ad", "da", "bc", "cb", "bd", "db",
                          "ee", "ef", "fe", "eg", "ge", "fg", "gf",
                          "hh", "hk", "kh", "kk", "h-", "k-",
                          "ll", "lm", "ml",
                          "nn", "np", "pn")
    for(genotype in unique(genotypes))
        if(! genotype %in% proper.genotypes){
            msg <- paste("unrecognized genotype", genotype)
            stop(msg, call.=FALSE)
        }

    ## replace all "ca" by "ac", etc -> speed-up next step?
    ## TODO

    ## convert genotypes to R/qtl code (mother=AB x father=CD)
    for(locus.id in 1:nb.loci){
        if(phase[locus.id] == "{0-}"){
            if(seg[locus.id] == "<lmxll>"){ # AB=lm ; CD=ll
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] == "ll"){ # AC or AD
                        genotypes[ind.id,locus.id] <- 5
                    } else if(genotypes[ind.id,locus.id] %in% c("lm","ml")){ # BC or BD
                        genotypes[ind.id,locus.id] <- 6
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id], "(should be <lmxll>)")
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{1-}"){
            if(seg[locus.id] == "<lmxll>"){ # AB=ml ; CD=ll
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] == "ll"){ # BC or BD
                        genotypes[ind.id,locus.id] <- 6
                    } else if(genotypes[ind.id,locus.id] %in% c("lm","ml")){ # AC or AD
                        genotypes[ind.id,locus.id] <- 5
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id], "(should be <lmxll>)")
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{-0}"){
            if(seg[locus.id] == "<nnxnp>"){ # AB=nn ; CD=np
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] == "nn"){ # AC or BC
                        genotypes[ind.id,locus.id] <- 7
                    } else if(genotypes[ind.id,locus.id] %in% c("np","pn")){ # AD or BD
                        genotypes[ind.id,locus.id] <- 8
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id], "(should be <nnxnp>)")
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{-1}"){
            if(seg[locus.id] == "<nnxnp>"){ # AB=nn ; CD=pn
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] == "nn"){ # AD or BD
                        genotypes[ind.id,locus.id] <- 8
                    } else if(genotypes[ind.id,locus.id] %in% c("np","pn")){ # AC or BC
                        genotypes[ind.id,locus.id] <- 7
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id], "(should be <nnxnp>)")
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{00}"){
            if(seg[locus.id] == "<abxcd>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("ac","ca")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("ad","da")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("bc","cb")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("bd","db")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<efxeg>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("ee")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("eg","ge")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("fe","ef")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("fg","gf")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<hkxhk>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("hh")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("hk","kh")){ # AD or BC
                        genotypes[ind.id,locus.id] <- 10
                    } else if(genotypes[ind.id,locus.id] %in% c("kk")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    } else if(genotypes[ind.id,locus.id] %in% c("h-","-h")){ # not BD
                        genotypes[ind.id,locus.id] <- 14
                    } else if(genotypes[ind.id,locus.id] %in% c("k-","-k")){ # not AC
                        genotypes[ind.id,locus.id] <- 11
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id])
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{01}"){
            if(seg[locus.id] == "<abxcd>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("ad","da")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("ac","ca")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("bd","db")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("bc","cb")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<efxeg>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("eg","ge")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("ee")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("fg","gf")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("fe","ef")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<hkxhk>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("hk","kh")){ # AC or BD
                        genotypes[ind.id,locus.id] <- 9
                    } else if(genotypes[ind.id,locus.id] %in% c("hh")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("kk")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("h-","-h")){ # not BC
                        genotypes[ind.id,locus.id] <- 12
                    } else if(genotypes[ind.id,locus.id] %in% c("k-","-k")){ # not AD
                        genotypes[ind.id,locus.id] <- 13
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id])
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{10}"){
            if(seg[locus.id] == "<abxcd>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("bc","cb")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("bd","db")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("ac","ca")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("ad","da")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<efxeg>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("fe","ef")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("fg","gf")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("ee")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("eg","ge")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<hkxhk>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("kh","hk")){ # AC or BD
                        genotypes[ind.id,locus.id] <- 9
                    } else if(genotypes[ind.id,locus.id] %in% c("kk")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("hh")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("h-","-h")){ # not AD
                        genotypes[ind.id,locus.id] <- 13
                    } else if(genotypes[ind.id,locus.id] %in% c("k-","-k")){ # not BC
                        genotypes[ind.id,locus.id] <- 12
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id])
                stop(msg, call.=FALSE)
            }
        } else if(phase[locus.id] == "{11}"){
            if(seg[locus.id] == "<abxcd>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("bd","db")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("bc","cb")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("ad","da")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("ac","ca")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<efxeg>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("fg","gf")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("fe","ef")){ # AD
                        genotypes[ind.id,locus.id] <- 3
                    } else if(genotypes[ind.id,locus.id] %in% c("eg","ge")){ # BC
                        genotypes[ind.id,locus.id] <- 2
                    } else if(genotypes[ind.id,locus.id] %in% c("ee")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    }
                }
            } else if(seg[locus.id] == "<hkxhk>"){
                for(ind.id in 1:nb.inds){
                    if(is.na(genotypes[ind.id,locus.id])){
                        next
                    } else if(genotypes[ind.id,locus.id] %in% c("kk")){ # AC
                        genotypes[ind.id,locus.id] <- 1
                    } else if(genotypes[ind.id,locus.id] %in% c("kh","hk")){ # AD or BC
                        genotypes[ind.id,locus.id] <- 10
                    } else if(genotypes[ind.id,locus.id] %in% c("hh")){ # BD
                        genotypes[ind.id,locus.id] <- 4
                    } else if(genotypes[ind.id,locus.id] %in% c("h-","-h")){ # not AC
                        genotypes[ind.id,locus.id] <- 11
                    } else if(genotypes[ind.id,locus.id] %in% c("k-","-k")){ # not BD
                        genotypes[ind.id,locus.id] <- 14
                    }
                }
            } else{
                msg <- paste("unrecognized segregation ", seg[locus.id],
                             "at locus", locus.id, "with phase",
                             phase[locus.id])
                stop(msg, call.=FALSE)
            }
        } else{
            msg <- paste("unrecognized phase", phase[locus.id],
                         "at locus", locus.id)
            stop(msg, call.=FALSE)
        }
    }
    storage.mode(genotypes) <- "numeric"

    list(pop.name=pop.name, pop.type=pop.type, genotypes=genotypes,
         seg=seg, phase=phase, classif=classif)
}

## returns a data.frame with 3 columns: chr (factor), marker (char), pos (num)
read.cross.mq.map <-
function(mapfile){
    lines <- readLines(mapfile)
    genmap <- data.frame(chr=rep(NA, length(lines)),
                         marker=NA, pos=NA)

    whole <- paste(lines, collapse="\n")
    groups <- strsplit(x=whole, split="group")[[1]]

    marker.id <- 1
    for(group.id in 1:length(groups)){
        tmp <- gsub(pattern=" + ", replacement=" ", x=groups[group.id])
        if(tmp == "" || tmp == " ") # empty groups
            next

        lines <- strsplit(tmp, split="\n")[[1]]
        group.name <- gsub(pattern=" ", replacement="", x=lines[1])

        for(line.id in 2:length(lines)){
            tmp <- mq.rmv.comment(x=lines[line.id], symbol=";")
            if(tmp == "" || tmp == " ") # empty lines
                next
            tokens <- strsplit(x=tmp, split=" |\t")[[1]]
            tmp <- gsub(pattern=" ", replacement="", x=tokens[2])
            genmap[marker.id, "chr"] <- group.name
            genmap[marker.id, "marker"] <- tokens[1]
            genmap[marker.id, "pos"] <- as.numeric(tokens[2])
            marker.id <- marker.id + 1
        }
    }

    genmap <- genmap[! is.na(genmap[,"chr"]), ]
    genmap[,"chr"] <- factor(genmap[,"chr"], levels=unique(genmap[,"chr"]))

    genmap
}

read.cross.mq.qua <-
function(quafile){
    nb.traits <- NULL
    nb.inds <- NULL
    miss <- NULL
    phenotypes <- NULL

    lines <- readLines(quafile)

    ind.id <- 1
    trait.id <- 1
    trait.names <- c()
    for(line.id in 1:length(lines)){

        ## skip empty lines
        tmp <- mq.rmv.comment(x=lines[line.id], symbol=";")
        tmp <- gsub(pattern=" + ", replacement=" ", x=tmp)
        if(tmp == "" || tmp == " ")
            next

        ## extract the number of traits
        if(grepl(pattern="ntrt", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            nb.traits <- as.numeric(strsplit(x=tmp, split="=")[[1]][2])
            next
        }

        ## extract the number of individuals
        if(grepl(pattern="nind", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            nb.inds <- as.numeric(strsplit(x=tmp, split="=")[[1]][2])
            next
        }

        ## extract the symbol for missing values
        if(grepl(pattern="miss", x=lines[line.id])){
            tmp <- gsub(pattern=" ", replacement="", x=tmp)
            miss <- strsplit(x=tmp, split="=")[[1]][2]
            next
        }

        tokens <- strsplit(x=tmp, split=" |\t")[[1]]

        if(trait.id <= nb.traits){
            if(length(tokens) == 1){ # one trait name per line
                trait.names <- c(trait.names, tokens[1])
                trait.id <- trait.id + 1
            } else{ # all trait names on the same line
                trait.names <- tokens
                trait.id <- nb.traits + 1
            }
            next
        } else if(trait.id > nb.traits && is.null(phenotypes)){
            if(length(trait.names) != nb.traits){
                msg <- paste("there seems to be fewer trait names (",
                             length(trait.names),
                             ") than indicated in the header (",
                             nb.traits, ")", sep="")
                stop(msg, call.=FALSE)
            }
            phenotypes <- matrix(NA, nrow=nb.inds, ncol=nb.traits)
        }

        if(length(tokens) != nb.traits){
            msg <- paste0("line ", line.id, " should have ", nb.traits,
                          " column", ifelse(nb.traits > 1, "s", ""),
                          " separated by spaces or tabs")
            stop(msg, call.=FALSE)
        }
        phenotypes[ind.id,] <- tokens
        ind.id <- ind.id + 1
    }

    phenotypes[which(phenotypes == miss)] <- NA
    phenotypes <- as.data.frame(phenotypes)

    for(j in 1:ncol(phenotypes))
        phenotypes[,j] <- tryCatch(expr=as.numeric(as.character(phenotypes[,j])),
                                   warning=function(w)
                                   as.factor(phenotypes[,j]))
    colnames(phenotypes) <- trait.names

    phenotypes
}

# end of read.cross.mq.R
