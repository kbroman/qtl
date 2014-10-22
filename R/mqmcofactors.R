#####################################################################
#
# mqmcofactors.R
#
# Copyright (c) 2009-2010, Danny Arends
#
# Modified by Karl Broman and Pjotr Prins
#
#
# first written Februari 2009
# last modified April 2010
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
# Contains: find.markerindex
#           mqmsetcofactors
#           mqmautocofactors
#           scoremissingmarkers
#           calculatedensity
#           mqmplot.cofactors
#           checkdistances
#
#####################################################################

######################################################################
#
# find.markerindex: Extracts the number of the marker when viewing the markers lineair
# mqmcofactors: Prepares a cofactor list to use with mqmscan
# mqmsetcofactors: Prepares a cofactor list to use with mqmscan
#
######################################################################


find.markerindex <- function(cross, name) {
    match(name, markernames(cross))
}

mqmsetcofactors <- function(cross,each = NULL,cofactors=NULL,sexfactors=NULL,verbose=FALSE){
    if(is.null(each) && is.null(cofactors))
        stop("Please set either the each parameter or the cofactors")
    if(missing(cross))
        stop("No cross file. Please supply a valid cross object.")
    individuals <- nind(cross)
    n.chr <- nchr(cross)
    n.mark <- sum(nmar(cross))
    cofactorlist <- rep(0,n.mark)
    if(!is.null(each)  && each > n.mark)
        stop("Not enough markers to place cofactors at this wide an interval.")


    if(verbose) {
        cat("INFO: Found",individuals,"individuals in the cross object.\n")
        cat("INFO: Mamimum amount of cofactors",(individuals-15)," (each =",ceiling(sum(n.mark)/(individuals-15)),") leaves 15 Degrees of Freedom (no Dominance).\n")
        cat("INFO: Mamimum amount of cofactors",(individuals-15)/2," (each =",ceiling(sum(n.mark)/(individuals-15))*2,") leaves 15 Degrees of Freedom (Dominance).\n")
    }

    if(is.null(cofactors)){
        cofactorlist <- rep(c(rep(0,each-1),1),(2*n.mark)/each)
        cofactorlist <- cofactorlist[1:n.mark]
    }else{
        if(max(cofactors) > n.mark)
            stop("Trying to set a non-existent marker as a cofactor.")
        if(min(cofactors) <= 0)
            stop("Trying to set a non-existent marker as a cofactor.")

        cofactorlist[cofactors]=1
        if(!is.null(sexfactors)){
            cofactorlist[sexfactors]=2
        }
    }
    if(sum(cofactorlist) > (individuals-15)){
        warning("Trying to set: ",ceiling(sum(n.mark)/each)," markers as cofactor. This leaves less than 15 Degrees of Freedom.\n")
    }
    cofactorlist
}

scoremissingmarkers <- function(cross){
    genotype <- pull.geno(cross)
    nind <- dim(genotype)[1]
    missing <- NULL
    for(x in 1:dim(genotype)[2]){
        missing <- c(missing,sum(is.na(genotype[,x]))/nind)
    }
    missing
}

calculatedensity <- function(cross,distance=30){
    genotype <- pull.geno(cross)
    densities <- NULL
    for(chr in 1:nchr(cross)){
        map <- pull.map(cross)[[chr]]
        for(x in 1:length(map)){
            densities <- c(densities,sum(map[which(map > map[x]-distance)] < map[x]+distance))
        }
    }
    densities
}

mqmautocofactors <- function(cross, num=50, distance=5,dominance=FALSE,plot=FALSE,verbose=FALSE){
    if(num > (nind(cross)-15) && !dominance){
        stop("Trying to set: ",num," markers as cofactor. This leaves less than 15 Degrees of Freedom.\n")
    }
    if(num > ((nind(cross)-15)/2) && dominance){
        stop("Trying to set: ",num," markers as cofactor. This leaves less than 15 Degrees of Freedom.\n")
    }
    if(distance < 0.1){
        distance <- 0.1
    }
    #  r <- scanone(cross)
    cofactors <- rep(0,sum(nmar(cross)))
    missing <- scoremissingmarkers(cross)
    densities <- calculatedensity(cross,distance*2)*missing
    cnt <- 0
    while(sum(cofactors) < num && cnt < num){
        lefttoset <- num - sum(cofactors)
        if(verbose) cat("Cofactors left",lefttoset,"/",num,"\n")
        possible <- which(max(densities)==densities)
        if(length(possible) > lefttoset){
            possible <- sample(possible,lefttoset)
        }
        cofactors[possible] <- 1
        densities[which(cofactors==1)] <- 0
        cofactors <- checkdistances(cross,cofactors,distance)
        cnt <- cnt+1
    }
    if(cnt==num && verbose) cat("Solution by iteration, there might be less cofactors then requested\n")
    if(plot) mqmplot.cofactors(cross,cofactors)
    cofactors
}

checkdistances <- function(cross,cofactors,dist=5){
    map <- unlist(pull.map(cross))
    newcofactors <- cofactors
    cnt_dropped <- 0
    for(x in which(cofactors==1)){
        for(y in which(cofactors==1)){
            if(x != y){
                chr_x <- strsplit(names(map[x]),'.',fixed=TRUE)[[1]][1]
                loc_x <- as.double(map[x])
                chr_y <- strsplit(names(map[y]),'.',fixed=TRUE)[[1]][1]
                loc_y <- as.double(map[y])
                if(chr_x==chr_y && abs(loc_x-loc_y) < dist){
                    newcofactors[y] <- 0
                    cnt_dropped <- cnt_dropped+1
                }
            }
        }
    }
    #cat("Dropped ",cnt_dropped," cofactors due to conflicting locations\n")
    newcofactors
}

mqmplot.cofactors <- function(cross,cofactors,...){
    map <- pull.map(cross)
    qc <- NULL
    qn <- NULL
    qp <- NULL
    mapnames <- NULL
    for(x in 1:nchr(cross)){
        mapnames <- c(mapnames,names(pull.map(cross)[[x]]))
    }
    chr <- 1
    genotype <- pull.geno(cross)
    for(x in 1:length(cofactors)){
        if(x > sum(nmar(cross)[1:chr])){
            chr <- chr+1
        }
        if(cofactors[x]>0){
            qn <- c(qn, mapnames[x])
            qc <- c(qc, as.character(names(map)[chr]))
            qp <- c(qp, as.double(unlist(map)[x]))
        }
    }
    plot(makeqtl(sim.geno(cross),qc,qp,qn),...)
}

# end of mqmcofactors.R
