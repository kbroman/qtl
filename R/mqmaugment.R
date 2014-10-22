#####################################################################
#
# mqmaugment.R
#
# Copyright (c) 2009-2013, Danny Arends
#
# Modified by Pjotr Prins; slight modification by Karl Broman
#
#
# first written Februari 2009
# last modified Sep 2013
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
# Contains: mqmaugment
#
#
#####################################################################

######################################################################
#
# mqmaugment_on_cofactors: Data Augmentation routine for MQM only using the cofactors, the other markers are filled by fill.geno()
#
######################################################################

mqmaugment_on_cofactors <- function(cross, cofactors, maxaugind=82, minprob=0.1, strategy=c("default","impute","drop"), verbose=FALSE, ...){
    markernames <- as.character(unlist(lapply(pull.map(cross),names)))
    todrop <- markernames[which(cofactors==0)]
    toaugment <- drop.markers(cross,todrop)
    cross <- fill.geno(cross, ...)
    augmented <- mqmaugment(toaugment,maxaugind,minprob,strategy,verbose)
    newgenomatrix <- NULL
    for(marker in markernames){
        if(marker %in% todrop){
            newgenomatrix <- cbind(newgenomatrix,pull.geno(cross)[augmented$mqm[[3]],marker])
        }else{
            newgenomatrix <- cbind(newgenomatrix,pull.geno(augmented)[,marker])
        }
    }
    colnames(newgenomatrix) <- markernames
    augmented$geno <- vector("list",nchr(cross))
    for(chr in 1:length(cross$geno)){
        cat(chr,"\n")
        augmented$geno[[chr]]$map <- cross$geno[[chr]]$map
        augmented$geno[[chr]]$data <- newgenomatrix[,names(cross$geno[[chr]]$map)]
        names(augmented$geno) <- names(cross$geno)
        class(augmented$geno[[chr]]) <- as.character(class(cross$geno[[chr]]))
    }
    augmented
}

######################################################################
#
# mqmaugment: dataaugmentation routine for MQM
#
######################################################################

mqmaugment <- function(cross,maxaugind=82, minprob=0.1, strategy=c("default","impute","drop"), verbose=FALSE) {
    # ---- check input supplied by user, start a timer
    if(minprob <= 0 || minprob > 1){ stop("Error minprob should be a value between 0 and 1.") }

    starttime <- proc.time()
    maxaug    <- nind(cross) * maxaugind   # maxaug is the maximum of individuals to augment to
    supported <- c("default","impute","drop")
    strategy  <- pmatch(strategy, supported)

    # ---- check for supported crosses and set ctype

    isF2       = 1
    isBC       = 2
    isRIL      = 3

    isAA       = 1
    isAB       = 2
    isH        = 2
    isBB       = 3
    isNOTBB    = 4
    isNOTAA    = 5
    isMISSING  = 9

    crosstype <- class(cross)[1]

    if (crosstype == "f2") {
        ctype = isF2
    }
    else if (crosstype == "bc" || crosstype == "dh" || crosstype=="haploid") {
        ctype = isBC
    }
    else if (crosstype == "riself") {
        ctype = isRIL
    }
    else {
        stop("Currently only F2, BC, and selfed RIL crosses can be analyzed by MQM.")
    }

    if (verbose) cat("INFO: Received a valid cross file type:", crosstype,".\n")

    # ---- Check sex chromosome
    # check whether the X chromosome should be dropped
    # (backcross with one sex should be fine)
    chrtype <- sapply(cross$geno, class)
    # Drop the X chromosome in F2 and related crosses
    if (any(chrtype == "X") && (ctype == isF2 || length(getgenonames(crosstype, "X", "full", getsex(cross), attributes(cross))) != 2)) {
        warning("MQM not yet available for the X chromosome; omitting chr ",
                paste(names(cross$geno)[chrtype == "X"], collapse=" "))
        cross <- subset(cross, chr=(chrtype != "X"))
    }

    # ---- Count
    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    n.aug <- maxaug
    if (verbose) {
        cat("INFO: Number of individuals:",n.ind,".\n")
        cat("INFO: Number of chr:",n.chr,".\n")
    }

    # ---- Genotype
    geno <- pull.geno(cross)
    chr <- rep(1:nchr(cross), nmar(cross))
    dist <- unlist(pull.map(cross))
    # Create a fake phenotype for augmentation
    pheno <- rep(1:n.ind)
    n.mark <- ncol(geno)
    if (verbose) cat("INFO: Number of markers:",n.mark,".\n")

    # Check for NA genotypes and replace them with a 9
    geno[is.na(geno)] <- isMISSING
    if (ctype==isRIL) {
        nH = sum(geno==isH)
        if (nH>0) {
            #warning("RIL dataset contains ", nH," heterozygous genotypes")
            if (any(geno==isBB)) { # have 3/BB's, so replace 2/H's with missing values
                geno[geno==isH] <- isMISSING
                #warning("Removed heterozygous genotypes from RIL set")
            } else {
                #warning("Converting heterozygous genotypes to BB from RIL set")
                geno[geno==isH] <- isBB
            }
        }
    } # end if(RIL)

    # ---- Call data augmentation
    result <- .C("R_mqmaugment",
                 as.integer(geno),
                 as.double(dist),
                 as.double(pheno),
                 augGeno=as.integer(rep(0,n.mark*maxaug)),
                 augPheno=as.double(rep(0,maxaug)),
                 augIND=as.integer(rep(0,maxaugind*n.ind)),
                 nind=as.integer(n.ind),
                 naug=as.integer(n.aug),
                 as.integer(n.mark),
                 as.integer(1),    # 1 phenotype
                 as.integer(maxaug),
                 as.integer(maxaugind),
                 as.double(minprob),
                 as.integer(chr),
                 as.integer(ctype),
                 as.integer(strategy),
                 as.integer(verbose),
                 PACKAGE="qtl")

    n.indold = n.ind
    n.ind = result$nind
    n.aug = result$naug
    markONchr <- 0
    markdone <- 0
    pheno <- NULL
    pnames <- phenames(cross)
    oldpheno <- pull.pheno(cross)
    result$augIND <- result$augIND+1
    for(x in result$augPheno[1:n.aug]){
        if(nphe(cross)>1){
            pheno <- rbind(pheno,oldpheno[x,])
        }else{
            pheno <- c(pheno,oldpheno[x])
        }
    }
    for(c in 1:n.chr){
        #print(paste("Cromosome",c,"\n",sep=""))
        matri <- NULL
        markONchr <- dim(cross$geno[[c]]$data)[2]
        #print(paste("# markers",markONchr,"\n",sep=""))
        for(j in markdone:(markdone+markONchr-1)){
            #print(paste("Start",markdone,":End",(markdone+markONchr-1),"\n",sep=""))
            ind2 <- NULL
            ind2 <- result$augGeno[(1+(j*maxaug)):(n.aug+(j*maxaug))]
            matri <- rbind(matri,ind2)
        }
        matri <- t(matri)
        if(markdone==0){
            colnames(matri) <- colnames(geno)[markdone:(markdone+markONchr)]
        } else {
            #print(paste("Markdone",markdone,"End",(markdone+markONchr-1)))
            colnames(matri) <- colnames(geno)[(markdone+1):(markdone+markONchr)]
        }
        cross$geno[[c]]$data <- matri
        markdone <- (markdone+markONchr)
    }
    if(nphe(cross) > 1){ colnames(pheno) <- colnames(cross$pheno) }   # Add the phenotype names if we have multiple phenotypes
    cross$pheno <- as.data.frame(pheno, stringsAsFactors=TRUE)
    # Store extra information (needed by the MQM algorithm) which individual was which original etc..
    cross$mqm$Nind <- n.ind
    cross$mqm$Naug <- n.aug
    result$augIND <- result$augIND-1
    cross$mqm$augIND <- result$augIND[1:n.aug]
    # ---- RESULTS
    endtime <- proc.time()
    if(n.ind != n.indold){
        if(verbose) warning("SERIOUS WARNING: Dropped ",abs(n.ind - n.indold)," original individuals.\n  Information lost, please increase minprob.")
    }
    if(verbose) cat("INFO: DATA-Augmentation took: ",round((endtime-starttime)[3], digits=3)," seconds\n")

    colnames(cross$pheno) <- pnames # fix up phenotype names problem

    cross  # return cross type
}
