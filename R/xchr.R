#####################################################################
#
# xchr.R
#
# copyright (c) 2004-2013, Karl W Broman
# last modified Sep, 2013
# first written Apr, 2004
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
# Contains: Utilities for dealing with the X chromosome.
#           getsex, getgenonames, reviseXdata, scanoneXnull
#           revisecovar, dropXcol
#           [See also fixXgeno.bc & fixXgeno.f2 in read.cross.R]
#
######################################################################

# get sex and pgm columns from phenotype data
getsex <-
    function(cross)
{
    type <- class(cross)[1]
    if(type != "bc" && type != "f2" && type != "4way") return(list(sex=NULL, pgm=NULL))

    phe.names <- names(cross$pheno)

    sex.column <- grep("^[Ss][Ee][Xx]$", phe.names)
    pgm.column <- grep("^[Pp][Gg][Mm]$", phe.names)

    if(length(sex.column)==0) { # no sex included
        sex <- NULL
    }
    else {
        if(length(sex.column)>1)
            warning("'sex' included multiple times.  Using the first one.")
        temp <- cross$pheno[,sex.column[1]]
        if(is.numeric(temp)) {
            if(any(!is.na(temp) & temp != 0 & temp != 1)) {
                warning("Sex column should be coded as 0=female 1=male; sex ignored.")
                sex <- NULL
            }
            else sex <- temp
        }
        else {
            if(!is.factor(temp)) temp <- as.factor(temp)

            if(length(levels(temp)) == 1) {
                if(levels(temp) == "F" || levels(temp)=="f" ||
                   toupper(levels(temp)) == "FEMALE") sex <- rep(0,nind(cross))
                else if(levels(temp) == "M" || levels(temp)=="m" ||
                        toupper(levels(temp)) == "MALE") sex <- rep(1,nind(cross))
                else
                    warning("Sex column should be coded as 0=female 1=male; sex ignored.")
            }
            else if(length(levels(temp)) > 2) {
                warning("Sex column should be coded as a two-level factor; sex ignored.")
                sex <- NULL
            }
            else { # is a factor with two levels
                lev <- levels(temp)
                if(length(grep("^[Ff]",lev))>0 &&
                   length(males <- grep("^[Mm]",lev))>0) {
                    temp <- as.character(temp)
                    sex <- rep(0,length(temp))
                    sex[is.na(temp)] <- NA
                    sex[!is.na(temp) & temp==lev[males]] <- 1
                }
                else
                    warning("Don't understand levels in sex column; sex ignored.")
            }
        }
    }

    if(length(pgm.column)==0 || type=="4way") { # no pgm included
        pgm <- NULL
    }
    else {
        if(length(pgm.column)>1)
            warning("'pgm' included multiple times.  Using the first one.")
        temp <- cross$pheno[,pgm.column[1]]
        if(!is.numeric(temp))
            temp <- as.numeric(temp)-1
        if(any(!is.na(temp) & temp != 0 & temp != 1)) {
            warning("pgm column should be coded as 0/1; pgm ignored.")
            pgm <- NULL
        }
        else pgm <- temp
    }

    if(!is.null(sex) && any(is.na(sex))) {
        if(all(sex[!is.na(sex)]==1)) {
            warning(sum(is.na(sex)), " individuals with missing sex; assuming they're male like the others")
            sex[is.na(sex)] <- 1
        }
        else if(all(sex[!is.na(sex)]==0)) {
            warning(sum(is.na(sex)), " individuals with missing sex; assuming they're female like the others")
            sex[is.na(sex)] <- 0
        }
        else {
            warning(sum(is.na(sex)), " individuals with missing sex; assuming they're female")
            sex[is.na(sex)] <- 0
        }
    }

    if(!is.null(pgm) && any(is.na(pgm))) {
        if(all(pgm[!is.na(pgm)]==1)) {
            warning(sum(is.na(pgm)), " individuals with missing pgm; assuming pgm==1 like the others")
            pgm[is.na(pgm)] <- 1
        }
        else if(all(pgm[!is.na(pgm)]==0)) {
            warning(sum(is.na(pgm)), " individuals with missing pgm; assuming pgm==0 like the others")
            pgm[is.na(pgm)] <- 0
        }
        else {
            warning(sum(is.na(pgm)), " individuals with missing pgm; assuming pgm==0")
            pgm[is.na(pgm)] <- 0
        }
    }

    list(sex=sex,pgm=pgm)
}



# get names of genotypes
# used in discan, effectplot, plotPXG, scanone, scantwo, vbscan, reviseXdata
# cross.attr gives the cross attributes
getgenonames <-
    function(type=c("f2","bc","riself","risib","4way","dh","haploid","special","bcsft"),
             chrtype=c("A","X"), expandX=c("simple","standard","full"),
             sexpgm, cross.attr)
{
    type <- match.arg(type)
    chrtype <- match.arg(chrtype)
    expandX <- match.arg(expandX)

    ## Treat bcsft as bc if no intercross generations; otherwise as f2.
    if(type == "bcsft") {
        if(cross.attr$scheme[2] == 0)
            type <- "bc"
        else
            type <- "f2"
    }

    if(chrtype=="X") {
        sex <- sexpgm$sex
        pgm <- sexpgm$pgm
    }

    if(type=="special") return(cross.attr$genotypes)

    if(missing(cross.attr) || !("alleles" %in% names(cross.attr))) {
        if(type == "4way") alleles <- LETTERS[1:4]
        else alleles <- LETTERS[1:2]
    }
    else
        alleles <- cross.attr$alleles

    tempgn <- c(paste(rep(alleles[1],2),collapse=""),
                paste(alleles,collapse=""),
                paste(rep(alleles[2],2),collapse=""),
                paste(alleles[1],"Y",sep=""),
                paste(alleles[2],"Y",sep=""))

    # get rid of missing sex and pgm values, if there are any
    if(chrtype=="X") {
        if(length(sex)>0) sex <- sex[!is.na(sex)]
        if(length(pgm)>0) pgm <- pgm[!is.na(pgm)]
    }

    if(type=="riself" || type=="risib" || type=="dh")
        gen.names <- tempgn[c(1,3)]

    else if(type=="haploid")
        gen.names <- alleles

    else if(type == "4way") {
        if(chrtype=="A")
            gen.names <- c(paste(alleles[1],alleles[3],sep=""),
                           paste(alleles[2],alleles[3],sep=""),
                           paste(alleles[1],alleles[4],sep=""),
                           paste(alleles[2],alleles[4],sep=""))
        else
            gen.names <- c(paste(alleles[1],alleles[3],sep=""),
                           paste(alleles[2],alleles[3],sep=""),
                           paste(alleles[1],"Y",sep=""),
                           paste(alleles[2],"Y",sep=""))
    }

    else if(type == "bc") {

        if(chrtype=="A") # autosome
            gen.names <- tempgn[1:2] # AA AB

        else { # X chromosome

            #                 simple     standard       full
            #   -both sexes   A-/AB/BY   AA/AB/AY/BY    same as std
            #   -all females  AA/AB      same           same
            #   -all males    AY/BY      same           same

            if(length(sex)==0 || all(sex==0)) # all females
                gen.names <- tempgn[1:2] # AA AB
            else if(all(sex==1)) # all males
                gen.names <- tempgn[4:5] # AY BY
            else { # some of each
                if(expandX == "simple")
                    gen.names <- c(paste(alleles[1], "-", sep=""),
                                   tempgn[c(2,5)]) # A-, AB, BY
                else gen.names <- tempgn[c(1,2,4,5)]  # AA,AB,AY,BY
            }
        }
    }

    else { # intercross
        if(chrtype == "A")  # autosomal
            gen.names <- tempgn[1:3]
        else { # X chromsome

            # both crosses     simple     standard         full
            #   -both sexes   A-/AB/B-    AA/AB/BB/AY/BY   AA/AB1/AB2/BB/AY/BY
            #   -all females  AA/AB/BB    same as simple   AA/AB1/AB2/BB
            #   -all males    AY/BY       same             same
            # forw cross
            #   -both sexes   A-/AB/BY    AA/AB/AY/BY      same as std
            #   -all females  AA/AB       same             same
            #   -all males    AY/BY       same             same
            # backw cross
            #   -both sexes   B-/AB/AY    BB/AB/AY/BY      same as std
            #   -all females  BB/AB       same             same
            #   -all males    AY/BY       same             same

            if(length(sex)==0 || all(sex==0)) { # all females
                if(length(pgm)==0 || all(pgm==0)) # all forw dir
                    gen.names <- tempgn[1:2] # AA AB
                else if(all(pgm==1))  # all backw dir
                    gen.names <- tempgn[3:2] # BB AB
                else { # some of each direction
                    if(expandX=="full")
                        gen.names <- c(tempgn[1],
                                       paste(tempgn[2],c("f","r"), sep=""),
                                       tempgn[3])
                    else gen.names <- tempgn[1:3]
                }
            }
            else if(all(sex==1))  # all males
                gen.names <- tempgn[4:5]
            else { # some of each sex
                if(length(pgm)==0 || all(pgm==0)) { # all forw
                    if(expandX=="simple")
                        gen.names <- c(paste(alleles[1],"-", sep=""),
                                       tempgn[c(2,5)])
                    else gen.names <- tempgn[c(1,2,4,5)]
                }
                else if (all(pgm==1)) { # all backw
                    if(expandX=="simple")
                        gen.names <- c(paste(alleles[2], "-",sep=""),
                                       tempgn[c(2,4)])
                    else gen.names <- tempgn[c(3,2,4,5)]
                }
                else { # some of each dir
                    if(expandX=="simple")
                        gen.names <- c(paste(alleles[1],"-",sep=""),
                                       tempgn[2],
                                       paste(alleles[2],"-",sep=""))
                    else if(expandX=="standard")
                        gen.names <- tempgn
                    else
                        gen.names <- c(tempgn[1],
                                       paste(tempgn[2],c("f","r"),sep=""),
                                       tempgn[3:5])
                }
            }
        }
    }

    gen.names
}

# revise genotype data, probabilities or imputations for the X chromosome
reviseXdata <-
    function(type=c("f2","bc","bcsft"), expandX=c("simple","standard","full"),
             sexpgm, geno, prob, draws, pairprob, cross.attr, force=FALSE)
{
    type <- match.arg(type)
    expandX <- match.arg(expandX)

    ## Treat bcsft as bc if no intercross generations; otherwise as f2.
    if(type == "bcsft") {
        if(cross.attr$scheme[2] == 0)
            type <- "bc"
        else
            type <- "f2"
    }

    sex <- sexpgm$sex
    pgm <- sexpgm$pgm

    notmissing <- (!missing(geno)) + (!missing(prob)) + (!missing(draws)) +
        (!missing(pairprob))
    if(notmissing == 0)
        stop("Provide one of geno, prob, draws, pairprob.")
    if(notmissing > 1)
        stop("Provide just one of geno, prob, draws, pairprob.")

    # get genonames
    genonames <- getgenonames(type, "X", expandX, sexpgm, cross.attr)

    if(type == "bc") { # backcross

        if(length(sex)==0 || ((all(sex==0) || all(sex==1)) && !force)) { # all one sex
            # no changes necessary
            if(!missing(geno)) return(geno)
            else if(!missing(prob)) {
                dimnames(prob)[[3]] <- genonames
                return(prob)
            }
            else if(!missing(draws))
                return(draws)
            else # pairprob
                return(pairprob)
        }

        else { # both sexes

            if(!missing(geno)) {
                gmale <- geno[sex==1,]
                if(expandX=="simple")
                    gmale[!is.na(gmale) & gmale==2] <- 3
                else {
                    gmale[!is.na(gmale) & gmale==1] <- 3
                    gmale[!is.na(gmale) & gmale==2] <- 4
                }
                geno[sex==1,] <- gmale
                return(geno)
            }

            else if(!missing(draws)) {
                gmale <- draws[sex==1,,]
                if(expandX=="simple")
                    gmale[gmale==2] <- 3
                else {
                    gmale[gmale==1] <- 3
                    gmale[gmale==2] <- 4
                }
                draws[sex==1,,] <- gmale
                return(draws)
            }

            else if(!missing(prob)) {
                dimprob <- dim(prob)
                dimprob[3] <- length(genonames)
                newprob <- array(0,dim=dimprob)
                dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
                newprob[sex==0,,1:2] <- prob[sex==0,,1:2]

                if(expandX=="simple") {
                    newprob[sex==1,,1] <- prob[sex==1,,1]
                    newprob[sex==1,,3] <- prob[sex==1,,2]
                }
                else {
                    newprob[sex==1,,3] <- prob[sex==1,,1]
                    newprob[sex==1,,4] <- prob[sex==1,,2]
                }
                return(newprob)
            }

            else { # pairprob
                dimpairprob <- dim(pairprob)
                dimpairprob[3] <- dimpairprob[4] <- length(genonames)
                newpairprob <- array(0,dim=dimpairprob)
                newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]

                if(expandX=="simple") {
                    newpairprob[sex==1,,1,1] <- pairprob[sex==1,,1,1]
                    newpairprob[sex==1,,1,3] <- pairprob[sex==1,,1,2]
                    newpairprob[sex==1,,3,1] <- pairprob[sex==1,,2,1]
                    newpairprob[sex==1,,3,3] <- pairprob[sex==1,,2,2]
                }
                else {
                    newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
                    newpairprob[sex==1,,3,4] <- pairprob[sex==1,,1,2]
                    newpairprob[sex==1,,4,3] <- pairprob[sex==1,,2,1]
                    newpairprob[sex==1,,4,4] <- pairprob[sex==1,,2,2]
                }
                return(newpairprob)
            }

        } # end of "both sexes" / backcross

    } # end of backcross

    else { # intercross

        if(length(sex)==0 || all(sex==0)) { # all females

            if(length(pgm)==0 || ((all(pgm==0) || all(pgm==1)) && !force)) { # one dir, females
                if(!missing(geno)) return(geno)
                else if(!missing(draws)) return(draws)
                else if(!missing(pairprob)) return(pairprob)
                else {
                    dimnames(prob)[[3]] <- genonames
                    return(prob)
                }
            }

            else { # both dir, females
                if(!missing(geno)) {
                    gback <- geno[pgm==1,]
                    if(expandX!="full") {
                        gback[!is.na(gback) & gback==1] <- 3
                        geno[pgm==1,] <- gback
                    }
                    else {
                        gback[!is.na(gback) & gback==1] <- 4
                        gback[!is.na(gback) & gback==2] <- 3
                        geno[pgm==1,] <- gback
                    }
                    return(geno)
                }
                else if(!missing(draws)) {
                    gback <- draws[pgm==1,,]
                    if(expandX!="full") {
                        gback[!is.na(gback) & gback==1] <- 3
                    }
                    else {
                        gback[!is.na(gback) & gback==1] <- 4
                        gback[!is.na(gback) & gback==2] <- 3
                    }
                    draws[pgm==1,,] <- gback
                    return(draws)
                }
                else if(!missing(prob)) {
                    dimprob <- dim(prob)
                    dimprob[3] <- length(genonames)
                    newprob <- array(0,dim=dimprob)
                    dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
                    newprob[pgm==0,,1:2] <- prob[pgm==0,,1:2]

                    if(expandX!="full") { # simple/standard
                        newprob[pgm==1,,3] <- prob[pgm==1,,1]
                        newprob[pgm==1,,2] <- prob[pgm==1,,2]
                    }
                    else {
                        newprob[pgm==1,,4] <- prob[pgm==1,,1]
                        newprob[pgm==1,,3] <- prob[pgm==1,,2]
                    }
                    return(newprob)
                }
                else { # pairprob
                    dimpairprob <- dim(pairprob)
                    dimpairprob[3] <- dimpairprob[4] <- length(genonames)
                    newpairprob <- array(0,dim=dimpairprob)
                    newpairprob[pgm==0,,1:2,1:2] <- pairprob[pgm==0,,,]

                    if(expandX!="full") { # simple/standard
                        newpairprob[pgm==1,,3,3] <- pairprob[pgm==1,,1,1]
                        newpairprob[pgm==1,,3,2] <- pairprob[pgm==1,,1,2]
                        newpairprob[pgm==1,,2,3] <- pairprob[pgm==1,,2,1]
                        newpairprob[pgm==1,,2,2] <- pairprob[pgm==1,,2,2]
                    }
                    else {
                        newpairprob[pgm==1,,4,4] <- pairprob[pgm==1,,1,1]
                        newpairprob[pgm==1,,4,3] <- pairprob[pgm==1,,1,2]
                        newpairprob[pgm==1,,3,4] <- pairprob[pgm==1,,2,1]
                        newpairprob[pgm==1,,3,3] <- pairprob[pgm==1,,2,2]
                    }
                    return(newpairprob)
                }
            }
        }
        else if(all(sex==1) && !force)  { # all males
            if(!missing(geno)) return(geno)
            else if(!missing(draws)) return(draws)
            else if(!missing(pairprob)) return(pairprob)
            else {
                dimnames(prob)[[3]] <- genonames
                return(prob)
            }
        }

        else { # both sexes

            if(length(pgm)==0 || all(pgm==0)) { # both sexes, forw dir
                if(!missing(geno)) {
                    gmale <- geno[sex==1,]
                    if(expandX=="simple")
                        gmale[!is.na(gmale) & gmale==2] <- 3
                    else {
                        gmale[!is.na(gmale) & gmale==1] <- 3
                        gmale[!is.na(gmale) & gmale==2] <- 4
                    }
                    geno[sex==1,] <- gmale
                    return(geno)
                }

                else if(!missing(draws)) {
                    gmale <- draws[sex==1,,]
                    if(expandX=="simple")
                        gmale[gmale==2] <- 3
                    else {
                        gmale[gmale==1] <- 3
                        gmale[gmale==2] <- 4
                    }
                    draws[sex==1,,] <- gmale
                    return(draws)
                }

                else if(!missing(prob)) {
                    dimprob <- dim(prob)
                    dimprob[3] <- length(genonames)
                    newprob <- array(0,dim=dimprob)
                    dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
                    newprob[sex==0,,1:2] <- prob[sex==0,,1:2]

                    if(expandX=="simple") {
                        newprob[sex==1,,1] <- prob[sex==1,,1]
                        newprob[sex==1,,3] <- prob[sex==1,,2]
                    }
                    else {
                        newprob[sex==1,,3] <- prob[sex==1,,1]
                        newprob[sex==1,,4] <- prob[sex==1,,2]
                    }
                    return(newprob)
                }

                else { # pairprob
                    dimpairprob <- dim(pairprob)
                    dimpairprob[3] <- dimpairprob[4] <- length(genonames)
                    newpairprob <- array(0,dim=dimpairprob)
                    newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]

                    if(expandX=="simple") {
                        newpairprob[sex==1,,1,1] <- pairprob[sex==1,,1,1]
                        newpairprob[sex==1,,1,3] <- pairprob[sex==1,,1,2]
                        newpairprob[sex==1,,3,1] <- pairprob[sex==1,,2,1]
                        newpairprob[sex==1,,3,3] <- pairprob[sex==1,,2,2]
                    }
                    else {
                        newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
                        newpairprob[sex==1,,3,4] <- pairprob[sex==1,,1,2]
                        newpairprob[sex==1,,4,3] <- pairprob[sex==1,,2,1]
                        newpairprob[sex==1,,4,4] <- pairprob[sex==1,,2,2]
                    }
                    return(newpairprob)
                }
            } # both sexes, forw dir

            if(all(pgm==1) && !force) { # both sexes, backw dir
                if(!missing(geno)) {
                    gmale <- geno[sex==1,]
                    if(expandX!="full") {
                        gmale[!is.na(gmale) & gmale==1] <- 3
                        gmale[!is.na(gmale) & gmale==2] <- 1
                    }
                    else {
                        gmale[!is.na(gmale) & gmale==1] <- 3
                        gmale[!is.na(gmale) & gmale==2] <- 4
                    }
                    geno[sex==1,] <- gmale
                    return(geno)
                }

                else if(!missing(draws)) {
                    gmale <- draws[sex==1,,]
                    if(expandX!="full") {
                        gmale[gmale==1] <- 3
                        gmale[gmale==2] <- 1
                    }
                    else {
                        gmale[gmale==1] <- 3
                        gmale[gmale==2] <- 4
                    }
                    draws[sex==1,,] <- gmale
                    return(draws)
                }

                else if(!missing(prob)) {
                    dimprob <- dim(prob)
                    dimprob[3] <- length(genonames)
                    newprob <- array(0,dim=dimprob)
                    dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
                    newprob[sex==0,,1:2] <- prob[sex==0,,1:2]

                    if(expandX=="simple") {
                        newprob[sex==1,,3] <- prob[sex==1,,1]
                        newprob[sex==1,,1] <- prob[sex==1,,2]
                    }
                    else {
                        newprob[sex==1,,3] <- prob[sex==1,,1]
                        newprob[sex==1,,4] <- prob[sex==1,,2]
                    }
                    return(newprob)
                }

                else { # pairprob
                    dimpairprob <- dim(pairprob)
                    dimpairprob[3] <- dimpairprob[4] <- length(genonames)
                    newpairprob <- array(0,dim=dimpairprob)
                    newpairprob[sex==0,,1:2,1:2] <- pairprob[sex==0,,,]

                    if(expandX=="simple") {
                        newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
                        newpairprob[sex==1,,1,3] <- pairprob[sex==1,,2,1]
                        newpairprob[sex==1,,3,1] <- pairprob[sex==1,,1,2]
                        newpairprob[sex==1,,1,1] <- pairprob[sex==1,,2,2]
                    }
                    else {
                        newpairprob[sex==1,,3,3] <- pairprob[sex==1,,1,1]
                        newpairprob[sex==1,,3,4] <- pairprob[sex==1,,1,2]
                        newpairprob[sex==1,,4,3] <- pairprob[sex==1,,2,1]
                        newpairprob[sex==1,,4,4] <- pairprob[sex==1,,2,2]
                    }
                    return(newpairprob)
                }
            } # both sexes, backw dir

            else { # both dir, both sexes

                if(!missing(geno)) {
                    gmale <- geno[sex==1,]
                    gfemaler <- geno[sex==0 & pgm==1,]
                    if(expandX=="simple") {
                        gmale[!is.na(gmale) & gmale==2] <- 3
                        gfemaler[!is.na(gfemaler) & gfemaler==1] <- 3
                    }
                    else if(expandX=="standard") {
                        gmale[!is.na(gmale) & gmale==1] <- 4
                        gmale[!is.na(gmale) & gmale==2] <- 5
                        gfemaler[!is.na(gfemaler) & gfemaler==1] <- 3
                    }
                    else {
                        gmale[!is.na(gmale) & gmale==1] <- 5
                        gmale[!is.na(gmale) & gmale==2] <- 6
                        gfemaler[!is.na(gfemaler) & gfemaler==1] <- 4
                        gfemaler[!is.na(gfemaler) & gfemaler==2] <- 3
                    }
                    geno[sex==1,] <- gmale
                    geno[sex==0 & pgm==1,] <- gfemaler
                    return(geno)
                }

                else if(!missing(draws)) {
                    gmale <- draws[sex==1,,]
                    gfemaler <- draws[sex==0 & pgm==1,,]
                    if(expandX=="simple") {
                        gmale[gmale==2] <- 3
                        gfemaler[gfemaler==1] <- 3
                    }
                    else if(expandX=="standard") {
                        gmale[gmale==1] <- 4
                        gmale[gmale==2] <- 5
                        gfemaler[gfemaler==1] <- 3
                    }
                    else {
                        gmale[gmale==1] <- 5
                        gmale[gmale==2] <- 6
                        gfemaler[gfemaler==1] <- 4
                        gfemaler[gfemaler==2] <- 3
                    }
                    draws[sex==1,,] <- gmale
                    draws[sex==0 & pgm==1,,] <- gfemaler
                    return(draws)
                }

                else if(!missing(prob)) {
                    dimprob <- dim(prob)
                    dimprob[3] <- length(genonames)
                    newprob <- array(0,dim=dimprob)
                    dimnames(newprob) <- c(dimnames(prob)[1:2],list(genonames))
                    newprob[sex==0 & pgm==0,,1:2] <- prob[sex==0 & pgm==0,,1:2]

                    if(expandX=="simple") {
                        newprob[sex==1,,1] <- prob[sex==1,,1]
                        newprob[sex==1,,3] <- prob[sex==1,,2]
                        newprob[sex==0 & pgm==1,,3] <- prob[sex==0 & pgm==1,,1]
                        newprob[sex==0 & pgm==1,,2] <- prob[sex==0 & pgm==1,,2]
                    }
                    else if(expandX=="standard") {
                        newprob[sex==1,,4] <- prob[sex==1,,1]
                        newprob[sex==1,,5] <- prob[sex==1,,2]
                        newprob[sex==0 & pgm==1,,3] <- prob[sex==0 & pgm==1,,1]
                        newprob[sex==0 & pgm==1,,2] <- prob[sex==0 & pgm==1,,2]
                    }
                    else {
                        newprob[sex==1,,5] <- prob[sex==1,,1]
                        newprob[sex==1,,6] <- prob[sex==1,,2]
                        newprob[sex==0 & pgm==1,,4] <- prob[sex==0 & pgm==1,,1]
                        newprob[sex==0 & pgm==1,,3] <- prob[sex==0 & pgm==1,,2]
                    }
                    return(newprob)
                }

                else { # pairprob
                    dimpairprob <- dim(pairprob)
                    dimpairprob[3] <- dimpairprob[4] <- length(genonames)
                    newpairprob <- array(0,dim=dimpairprob)
                    newpairprob[sex==0 & pgm==0,,1:2,1:2] <- pairprob[sex==0 & pgm==0,,,]

                    male <- (sex==1)
                    femaler <- (sex==0) & (pgm==1)
                    if(expandX=="simple") {
                        newpairprob[male,,1,1] <- pairprob[male,,1,1]
                        newpairprob[male,,1,3] <- pairprob[male,,1,2]
                        newpairprob[male,,3,1] <- pairprob[male,,2,1]
                        newpairprob[male,,3,3] <- pairprob[male,,2,2]

                        newpairprob[femaler,,3,3] <- pairprob[femaler,,1,1]
                        newpairprob[femaler,,3,2] <- pairprob[femaler,,1,2]
                        newpairprob[femaler,,2,3] <- pairprob[femaler,,2,1]
                        newpairprob[femaler,,2,2] <- pairprob[femaler,,2,2]
                    }
                    else if(expandX=="standard") {
                        newpairprob[male,,4,4] <- pairprob[male,,1,1]
                        newpairprob[male,,4,5] <- pairprob[male,,1,2]
                        newpairprob[male,,5,4] <- pairprob[male,,2,1]
                        newpairprob[male,,5,5] <- pairprob[male,,2,2]

                        newpairprob[femaler,,3,3] <- pairprob[femaler,,1,1]
                        newpairprob[femaler,,3,2] <- pairprob[femaler,,1,2]
                        newpairprob[femaler,,2,3] <- pairprob[femaler,,2,1]
                        newpairprob[femaler,,2,2] <- pairprob[femaler,,2,2]
                    }
                    else {
                        newpairprob[male,,5,5] <- pairprob[male,,1,1]
                        newpairprob[male,,5,6] <- pairprob[male,,1,2]
                        newpairprob[male,,6,5] <- pairprob[male,,2,1]
                        newpairprob[male,,6,6] <- pairprob[male,,2,2]

                        newpairprob[femaler,,4,4] <- pairprob[femaler,,1,1]
                        newpairprob[femaler,,4,3] <- pairprob[femaler,,1,2]
                        newpairprob[femaler,,3,4] <- pairprob[femaler,,2,1]
                        newpairprob[femaler,,3,3] <- pairprob[femaler,,2,2]
                    }
                    return(newpairprob)
                }

            }
        }

    } # end of intercross

}

######################################################################
# scanoneXnull
#
# figure out null hypothesis business for scanone on X chromosome
######################################################################
scanoneXnull <-
    function(type, sexpgm, cross.attr)
{
    sex <- sexpgm$sex
    pgm <- sexpgm$pgm

    if(type == "risib" || type=="riself" || type=="dh" || type=="haploid") type <- "bc"

    if(type == "bcsft") {
        if(cross.attr$scheme[2] == 0)
            type <- "bc"
        else
            type <- "f2"
    }

    ### first figure out sex/pgm pattern

    # sex
    if(length(sex)==0 || all(sex==0)) { # all female
        onesex <- allfemale <- TRUE
    }
    else if(all(sex==1)) { # all male
        onesex <- TRUE
        allfemale <- FALSE
    }
    else { # both sexes
        onesex <- allfemale <- FALSE
    }
    # pgm
    if(length(pgm)==0 || all(pgm==0) || all(pgm==1)) # one direction
        onedir <- TRUE
    else onedir <- FALSE

    allmale <- onesex && !allfemale
    bothsex <- !onesex
    bothdir <- !onedir


    ### now figure out the null hypothesis and pull out appropriate
    ### covariates for the null

    # backcross, one sex
    # OR intercross, one dir and one sex
    # OR intercross, both dir and all male
    if((type=="bc" && onesex) ||
       (type=="f2" && ((onedir && onesex) || (bothdir && allmale)))) {
        adjustX <- FALSE
        parX0 <- 1
        sexpgmcovar <- sexpgmcovar.alt <- NULL
    }

    # backcross, both sexes
    # OR intercross, one direction and both sexes
    else if((type=="bc" && bothsex) ||
            (type=="f2" && onedir && bothsex)) {
        adjustX <- TRUE
        parX0 <- 2
        sexpgmcovar <- cbind(sex)
        sexpgmcovar.alt <- sex+1
    }

    # intercross, both dir and all female
    else if(type=="f2" && bothdir && allfemale) {
        adjustX <- TRUE
        parX0 <- 2
        sexpgmcovar <- cbind(pgm)
        sexpgmcovar.alt <- pgm+1
    }

    # intercross, both dir and both sexes
    else {
        adjustX <- TRUE
        parX0 <- 3
        sexpgmcovar <- cbind(sex,as.numeric(sex==0 & pgm==1))
        sexpgmcovar.alt <- rep(3,length(sex))
        sexpgmcovar.alt[sex==0 & pgm==0] <- 1
        sexpgmcovar.alt[sex==0 & pgm==1] <- 2
    }

    list(adjustX=adjustX, parX0=parX0, sexpgmcovar=sexpgmcovar,
         sexpgmcovar.alt=sexpgmcovar.alt)
}

######################################################################
# revisecovar
#
# Drop sex and pgm and their interxn as covariates for the X chr.
######################################################################
revisecovar <-
    function(sexpgm, covar)
{
    if(is.null(covar) || (is.null(sexpgm$sex) && is.null(sexpgm$pgm))) {
        if(!is.null(covar)) attr(covar, "n.dropped") <- 0
        return(covar)
    }

    covar <- as.matrix(covar)

    sex <- sexpgm$sex
    pgm <- sexpgm$pgm

    if(!is.null(pgm) && length(unique(pgm))==1) pgm <- NULL
    allfemale <- FALSE
    if(is.null(sex)) allfemale <- TRUE
    else {
        if(all(sex==0)) {
            allfemale <- TRUE
            sex <- NULL
        }
        else if(all(sex==1)) {
            allfemale <- FALSE
            sex <- NULL
        }
    }

    if(!is.null(pgm)) { # some of each direction
        if(!is.null(sex)) { # some of each sex
            femf <- as.numeric(pgm==0 & sex==0)
            femr <- as.numeric(pgm==1 & sex==0)
            mal <- sex
            X <- cbind(femf, femr, mal)
        }
        else { # all of one sex
            if(allfemale)
                X <- cbind(1-pgm, pgm)
            else
                X <- cbind(rep(1, nrow(covar)))
        }
    }
    else { # all of one direction
        if(!is.null(sex))  # some of each sex
            X <- cbind(sex, 1-sex)
        else X <- cbind(rep(1, nrow(covar)))
    }

    nc <- ncol(X)

    keep <- rep(TRUE,ncol(covar))
    for(i in 1:ncol(covar)) {
        if(qr(cbind(X,covar[,i]))$rank <= nc)
            keep[i] <- FALSE
    }
    if(!any(keep))
        covar <- numeric(0)
    else
        covar <- covar[,keep,drop=FALSE]

    attr(covar, "n.dropped") <- sum(!keep)
    covar
}

######################################################################
# dropXcol: for use with scantwo() for the X chromosome:
#           figure out what columns to drop...both for the full model
#           and for the additive model.
######################################################################
dropXcol <-
    function(type=c("f2","bc", "riself", "risib", "4way", "dh", "haploid", "special","bcsft"),
             sexpgm, cross.attr)
{
    type <- match.arg(type)

    ## Treat bcsft as bc if no intercross generations; otherwise as f2.
    if(type == "bcsft") {
        if(cross.attr$scheme[2] == 0)
            type <- "bc"
        else
            type <- "f2"
    }

    gn <- getgenonames(type, "X", "full", sexpgm, cross.attr)

    if(length(gn)==2) return(rep(0,4))

    if(length(gn)==4) return( c(0,0,0,0,0,1,0,  0,1,1,1,1,1,1,1,0) )

    if(length(gn)==6) {
        todrop <- c(rep(0,11), rep(1,25))
        todrop[c(8,10)] <- 1
        todrop[11+c(1,13,25)] <- 0
        return(todrop)
    }

    return(rep(0,length(gn)^2))
}


# end of xchr.R
