#####################################################################
#
# mqmcofactors.R
#
# Copyright (c) 2009-2010, Danny Arends
#
# Modified by Karl Broman and unknown and Pjotr Prins
#
# 
# first written Februari 2009
# last modified March 2010
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
# Contains: which.marker
#           mqmcofactors
#           mqmsetcofactors
#           mqmautocofactors
#           
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

mqmcofactors <- function(cross,cofactors,sexfactors,verbose=FALSE){
	if(missing(cross))
		ourstop("No cross file. Please supply a valid cross object.")

	if(missing(cofactors))
		ourstop("Cofactors to set. Please supply a list of markers to serve as cofactors.")
		
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL
	individuals <- nind(cross)
        if(verbose) {
          cat("INFO: Found",individuals,"individuals in the cross object.\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-10),"leaves 10 Degrees of Freedom, (No Dominance).\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-10)/2,"leaves 10 Degrees of Freedom (Dominance).\n")
        }
	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	
	n.mark <- ncol(geno)

	if(max(cofactors) > n.mark){
		ourstop("Trying to set a non-existent marker as a cofactor.")	  
	}

	if(min(cofactors) <= 0){
		ourstop("Trying to set a non-existent marker as a cofactor.")	  
	}
	
	if(!missing(sexfactors)){
		if(max(sexfactors) > n.mark){
			ourstop("Trying to set a non-existent marker as a sexfactor.")
		}
		if(min(sexfactors) <= 0){
			ourstop("Trying to set a non-existent marker as a sexfactor.")
		}
	}else{
    sexfactors <- NULL
  }

  cofactorlist <- rep(0,n.mark)
	for(i in 1:length(cofactors)) {
	  cofactorlist[cofactors[i]]=1
	}
	if(!is.null(sexfactors)){
      for(i in 1:length(sexfactors)) {
	    cofactorlist[sexfactors[i]]=2
	  }
	}
  if(sum(cofactorlist) > (individuals-15)){
		stop("Trying to set: ",sum(cofactorlist)," markers as cofactor. This leaves less than 15 Degrees of Freedom.\n")
	}
  cofactorlist
}

mqmsetcofactors <- function(cross,each = 3,verbose=FALSE){
	if(missing(cross))
          ourstop("No cross file. Please supply a valid cross object.")

	individuals <- nind(cross)
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL

	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	n.mark <- ncol(geno)
	
        if(verbose) {
          cat("INFO: Found",individuals,"individuals in the cross object.\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-15)," (each =",ceiling(sum(n.mark)/(individuals-15)),") leaves 15 Degrees of Freedom (no Dominance).\n")
          cat("INFO: Mamimum amount of cofactors",(individuals-15)/2," (each =",ceiling(sum(n.mark)/(individuals-15))*2,") leaves 15 Degrees of Freedom (Dominance).\n")
        }

	if(each > n.mark){
      ourstop("Not enough markers to place cofactors at.")
	  return 
	}	
	
    cofactorlist <- rep(0,n.mark)
	for(i in 1:n.mark) {
		if(i%%as.integer(each)==0){
			cofactorlist[i] = 1
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

mqmautocofactors <- function(cross, num=50, distance=5,dominance=FALSE,plot=FALSE){
  if(num > (nind(cross)-15) && !dominance){
		stop("Trying to set: ",num," markers as cofactor. This leaves less than 15 Degrees of Freedom.\n")
	}
  if(num > ((nind(cross)-15)/2) && dominance){
		stop("Trying to set: ",num," markers as cofactor. This leaves less than 15 Degrees of Freedom.\n")
	}
  if(distance < 0.1){
    distance <- 0.1
  }
  r <- scanone(cross)
  cofactors <- rep(0,sum(nmar(cross)))
  missing <- scoremissingmarkers(cross)
  densities <- calculatedensity(cross,distance*2)*missing
  cnt <- 0
  while(sum(cofactors) < num && cnt < num){
    lefttoset <- num - sum(cofactors)
    cat("Cofactors left",lefttoset,"/",num,"\n")
    possible <- which(max(densities)==densities)
    if(length(possible) > lefttoset){
      possible <- sample(possible,lefttoset)
    }
    cofactors[possible] <- 1
    densities[which(cofactors==1)] <- 0
    cofactors <- checkdistances(cross,cofactors,distance)
    cnt <- cnt+1
  }
  if(cnt==num) cat("Solution by itteration, there might be less cofactors then requested\n")
  if(plot) plotcofactors(cross,cofactors)
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

plotcofactors <- function(cross,cofactors){
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
  plot(makeqtl(sim.geno(cross),qc,qp,qn),col="blue")
}

# end of mqmcofactors.R
