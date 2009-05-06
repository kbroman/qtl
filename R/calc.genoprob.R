######################################################################
#
# calc.genoprob.R
#
# copyright (c) 2001-9, Karl W Broman
# last modified Apr, 2009
# first written Feb, 2001
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
# Contains: calc.genoprob
#
######################################################################

######################################################################
#
# calc.genoprob: calculate genotype probabilities conditional on 
#                observed marker genotypes
#
######################################################################

calc.genoprob <-
function(cross, step=0, off.end=0, error.prob=0.0001,
         map.function=c("haldane","kosambi","c-f","morgan"),
         stepwidth=c("fixed", "variable"))
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  # map function
  map.function <- match.arg(map.function)
  if(map.function=="kosambi") mf <- mf.k
  else if(map.function=="c-f") mf <- mf.cf
  else if(map.function=="morgan") mf <- mf.m
  else mf <- mf.h

  stepwidth <- match.arg(stepwidth)
 
  # don't let error.prob be exactly zero (or >1)
  if(error.prob < 1e-50) error.prob <- 1e-50
  if(error.prob > 1) {
    error.prob <- 1-1e-50
    warning("error.prob shouldn't be > 1!")
  }

  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)

  type <- class(cross)[1]

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {
    if(n.mar[i]==1) temp.offend <- max(c(off.end,5))
    else temp.offend <- off.end
    
    chrtype <- class(cross$geno[[i]])
    if(chrtype=="X") xchr <- TRUE
    else xchr <- FALSE

    # which type of cross is this?
    if(type == "f2") {
      one.map <- TRUE
      if(!xchr) { # autosomal
        cfunc <- "calc_genoprob_f2"
        n.gen <- 3
        gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
      }
      else {                             # X chromsome 
        cfunc <- "calc_genoprob_bc"
        n.gen <- 2
        gen.names <- c("g1","g2")
      }
    }
    else if(type == "bc") {
      cfunc <- "calc_genoprob_bc"
      n.gen <- 2
      if(!xchr)
        gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
      else gen.names <- c("g1","g2")
      one.map <- TRUE
    }
    else if(type == "riself" || type=="risib" || type=="dh") {
      cfunc <- "calc_genoprob_bc"
      n.gen <- 2
      gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
      one.map <- TRUE
    }
    else if(type == "4way") {
      cfunc <- "calc_genoprob_4way"
      n.gen <- 4
      one.map <- FALSE
      gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
    }
    else if(type=="ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self") {
      cfunc <- paste("calc_genoprob_", type, sep="")
      n.gen <- as.numeric(substr(type, 3, 3))
      one.map <- TRUE
      gen.names <- LETTERS[1:n.gen]
      if(xchr)
        warning("calc.genoprob not working properly for the X chromosome for 4- or 8-way RIL.")
    }
    else 
      stop("calc.genoprob not available for cross type ", type, ".")

    # genotype data
    gen <- cross$geno[[i]]$data
    gen[is.na(gen)] <- 0
    
    # recombination fractions
    if(one.map) {
      # recombination fractions
      map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
      rf <- mf(diff(map))
      if(type=="risib" || type=="riself")
        rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),chrtype)
      rf[rf < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=length(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,names(map))
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
      marnames <- names(map)
    }
    else {
      map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
      rf <- mf(diff(map[1,]))
      rf[rf < 1e-14] <- 1e-14
      rf2 <- mf(diff(map[2,]))
      rf2[rf2 < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,dimnames(map)[[2]])
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
      marnames <- colnames(map)
    }

    # call the C function
    if(one.map) {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(newgen),        # genotype data
              as.double(rf),             # recombination fractions
              as.double(error.prob),     # 
              genoprob=as.double(rep(0,n.gen*n.ind*n.pos)),
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
              genoprob=as.double(rep(0,n.gen*n.ind*n.pos)),
              PACKAGE="qtl")
    }

    # re-arrange marginal probabilites
    cross$geno[[i]]$prob <- array(z$genoprob,dim=c(n.ind,n.pos,n.gen))
    dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, gen.names)
    # attribute set to the error.prob value used, for later
    #     reference, especially by calc.errorlod()

    attr(cross$geno[[i]]$prob, "map") <- map
    attr(cross$geno[[i]]$prob,"error.prob") <- error.prob
    attr(cross$geno[[i]]$prob,"step") <- step
    attr(cross$geno[[i]]$prob,"off.end") <- temp.offend
    attr(cross$geno[[i]]$prob,"map.function") <- map.function
    attr(cross$geno[[i]]$prob,"stepwidth") <- stepwidth
  } # end loop over chromosomes

  # 4- and 8-way RIL: reorganize the results
  if(type=="ri4self" || type=="ri4sib" || type=="ri8self" || type=="ri8sib") 
    cross <- reorgRIgenoprob(cross)

  cross
}

######################################################################
#
# calc.genoprob.special: 
#    special version used by calc.errorlod
#    for each individual and marker, calculate probabilities allowing
#    that genotype to be in error but assuming that all other genotypes
#    are correct
#
######################################################################

calc.genoprob.special <-
function(cross, error.prob=0.0001,
         map.function=c("haldane","kosambi","c-f","morgan"))
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  step <- 0
  off.end <- 0
  stepwidth <- "fixed"

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
  n.mar <- nmar(cross)

  type <- class(cross)[1]

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {
    if(n.mar[i]==1) temp.offend <- max(c(off.end,5))
    else temp.offend <- off.end
    
    chrtype <- class(cross$geno[[i]])
    if(chrtype=="X") xchr <- TRUE
    else xchr <- FALSE

    # which type of cross is this?
    if(type == "f2") {
      one.map <- TRUE
      if(!xchr) { # autosomal
        cfunc <- "calc_genoprob_special_f2"
        n.gen <- 3
        gen.names <- getgenonames("f2", "A", cross.attr=attributes(cross))
      }
      else {                             # X chromsome 
        cfunc <- "calc_genoprob_special_bc"
        n.gen <- 2
        gen.names <- c("g1","g2")
      }
    }
    else if(type == "bc") {
      cfunc <- "calc_genoprob_special_bc"
      n.gen <- 2
      if(!xchr) 
        gen.names <- getgenonames("bc", "A", cross.attr=attributes(cross))
      else gen.names <- c("g1","g2")
      one.map <- TRUE
    }
    else if(type == "riself" || type=="risib" || type=="dh") {
      cfunc <- "calc_genoprob_special_bc"
      n.gen <- 2
      gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
      one.map <- TRUE
    }
    else if(type == "4way") {
      cfunc <- "calc_genoprob_special_4way"
      n.gen <- 4
      one.map <- FALSE
      gen.names <- getgenonames(type, "A", cross.attr=attributes(cross))
    }
    else if(type=="ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self") {
      cfunc <- paste("calc_genoprob_special_", type, sep="")
      n.gen <- as.numeric(substr(type, 3, 3))
      one.map <- TRUE
      gen.names <- LETTERS[1:n.gen]
      if(xchr)
        warning("calc.genoprob.special not working properly for the X chromosome for 4- or 8-way RIL.")
    }
    else 
      stop("calc.genoprob.special not available for cross type ", type, ".")

    # genotype data
    gen <- cross$geno[[i]]$data
    gen[is.na(gen)] <- 0
    
    # recombination fractions
    if(one.map) {
      # recombination fractions
      map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
      rf <- mf(diff(map))
      if(type=="risib" || type=="riself")
        rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),chrtype)
      rf[rf < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=length(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,names(map))
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
      marnames <- names(map)
    }
    else {
      map <- create.map(cross$geno[[i]]$map,step,temp.offend,stepwidth)
      rf <- mf(diff(map[1,]))
      rf[rf < 1e-14] <- 1e-14
      rf2 <- mf(diff(map[2,]))
      rf2[rf2 < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,dimnames(map)[[2]])
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
      marnames <- colnames(map)
    }

    # call the C function
    if(one.map) {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(newgen),        # genotype data
              as.double(rf),             # recombination fractions
              as.double(error.prob),     # 
              genoprob=as.double(rep(0,n.gen*n.ind*n.pos)),
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
              genoprob=as.double(rep(0,n.gen*n.ind*n.pos)),
              PACKAGE="qtl")
    }

    # re-arrange marginal probabilites
    cross$geno[[i]]$prob <- array(z$genoprob,dim=c(n.ind,n.pos,n.gen))
    dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, gen.names)
    # attribute set to the error.prob value used, for later
    #     reference, especially by calc.errorlod()

    attr(cross$geno[[i]]$prob, "map") <- map
    attr(cross$geno[[i]]$prob,"error.prob") <- error.prob
    attr(cross$geno[[i]]$prob,"step") <- step
    attr(cross$geno[[i]]$prob,"off.end") <- temp.offend
    attr(cross$geno[[i]]$prob,"map.function") <- map.function
    attr(cross$geno[[i]]$prob,"stepwidth") <- stepwidth
  } # end loop over chromosomes

  cross
}

# end of calc.genoprob.R
