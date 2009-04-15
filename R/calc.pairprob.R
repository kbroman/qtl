######################################################################
#
# calc.pairprob.R
#
# copyright (c) 2001-8, Karl W Broman
# last modified Jun, 2008
# first written Nov, 2001
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
# Contains: calc.pairprob
#
######################################################################

######################################################################
#
# calc.pairprob: calculate joint genotype probabilities for all pairs
#                of putative QTLs, conditional on the observed marker
#                data
#
# This is an *internal* function, not to be called by the user.
#
# The input argument cross is assumed to have just one chromosome.
#
######################################################################

calc.pairprob <-
function(cross, step=0, off.end=0, error.prob=0.0001, 
         map.function=c("haldane","kosambi","c-f","morgan"),
         map, assumeCondIndep=FALSE)
{
  if(assumeCondIndep) { # assume conditional independence of QTL given markers
    if(!("prob" %in% names(cross$geno[[1]]))) {
      cross <- calc.genoprob(subset(cross, chr=1), step=step, off.end=off.end,
                             error.prob=error.prob, map.function=map.function)
    }
    prob <- cross$geno[[1]]$prob
    n.ind <- dim(prob)[1]
    n.pos <- dim(prob)[2]
    n.gen <- dim(prob)[3]

    if(n.pos < 2) stop("Must have > 1 position.")

    z <- .C("R_calc_pairprob_condindep",
            as.integer(n.ind),
            as.integer(n.pos),
            as.integer(n.gen),
            as.double(prob),
            pairprob=as.double(rep(0,n.ind*choose(n.pos, 2)*n.gen*n.gen)),
            PACKAGE="qtl")
    return(array(z$pairprob, dim=c(n.ind,n.pos*(n.pos-1)/2,n.gen,n.gen)))
  }

  if(step==0 && off.end > 0) step <- off.end*2

  # map function
  map.function <- match.arg(map.function)
  if(map.function=="kosambi") mf <- mf.k
  else if(map.function=="c-f") mf <- mf.cf
  else if(map.function=="morgan") mf <- mf.m
  else mf <- mf.h
 
  # don't let error.prob be exactly zero (or >1)
  if(error.prob < 1e-50) error.prob <- 1e-50
  if(error.prob > 1) {
    error.prob <- 1-1e-50
    warning("error.prob shouldn't be > 1!")
  }
  n.ind <- nind(cross)
  n.chr <- nchr(cross)

  # which type of cross is this?
  type <- class(cross)[1]
  if(type == "f2") {
    one.map <- TRUE
    if(class(cross$geno[[1]]) == "A") { # autosomal
      cfunc <- "calc_pairprob_f2"
      n.gen <- 3
      gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
    }
    else {                             # X chromsome 
      cfunc <- "calc_pairprob_bc"
      n.gen <- 2
      gen.names <- c("g1","g2")
    }
  }
  else if(type == "bc") {
    cfunc <- "calc_pairprob_bc"
    n.gen <- 2
    if(class(cross$geno[[1]]) == "A")
      gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
    else gen.names <- c("g1","g2")
    one.map <- TRUE
  }
  else if(type == "riself" || type=="risib" || type=="dh") {
    cfunc <- "calc_pairprob_bc"
    n.gen <- 2
    gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
    one.map <- TRUE
  }
  else if(type == "4way") {
    cfunc <- "calc_pairprob_4way"
    n.gen <- 4
    one.map <- FALSE
    gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
  }
  else 
    stop("calc.pairprob not available for cross type ", type, ".")

  # genotype data
  gen <- cross$geno[[1]]$data
  gen[is.na(gen)] <- 0
  
  # get recombination fractions
  if(one.map) {
#    map <- create.map(cross$geno[[1]]$map,step,off.end)
    rf <- mf(diff(map))
    if(type=="risib" || type=="riself")
      rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),class(cross$geno[[1]]))
    rf[rf < 1e-14] <- 1e-14
    
    # new genotype matrix with pseudomarkers filled in
    newgen <- matrix(ncol=length(map),nrow=nrow(gen))
    colnames(newgen) <- names(map)
    newgen[,colnames(gen)] <- gen
    newgen[is.na(newgen)] <- 0
    n.pos <- ncol(newgen)
    marnames <- names(map)
  }
  else {
#    map <- create.map(cross$geno[[1]]$map,step,off.end)
    rf <- mf(diff(map[1,]))
    rf[rf < 1e-14] <- 1e-14
    rf2 <- mf(diff(map[2,]))
    rf2[rf2 < 1e-14] <- 1e-14
    
    # new genotype matrix with pseudomarkers filled in
    newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
    colnames(newgen) <- colnames(map)
    newgen[,colnames(gen)] <- gen
    newgen[is.na(newgen)] <- 0
    n.pos <- ncol(newgen)
    marnames <- colnames(map)
  }
  
  if(n.pos < 2) return(NULL)

  # below: at least two positions
  # call the C function
  if(one.map) {
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.pos),         # number of markers
            as.integer(newgen),        # genotype data
            as.double(rf),             # recombination fractions
            as.double(error.prob),     # 
            as.double(rep(0,n.gen*n.ind*n.pos)),
            pairprob=as.double(rep(0,n.ind*n.pos*(n.pos-1)/2*n.gen^2)),
            PACKAGE="qtl")
  }
  else {
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.pos),         # number of markers
            as.integer(newgen),        # genotype data
            as.double(rf),             # recombination fractions
            as.double(rf2),            # recombination fractions
            as.double(error.prob),     # 
            as.double(rep(0,n.gen*n.ind*n.pos)),
            pairprob=as.double(rep(0,n.ind*n.pos*(n.pos-1)/2*n.gen^2)),
            PACKAGE="qtl")
  }
  
  array(z$pairprob, dim=c(n.ind,n.pos*(n.pos-1)/2,n.gen,n.gen))
}

# end of calc.pairprob.R
