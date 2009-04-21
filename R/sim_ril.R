#####################################################################
#
# sim_ril.R
#
# copyright (c) 2004-9, Karl W Broman
# last modified Apr, 2009
# first written May, 2004
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
# Contains: sim.ril, sim.cc
#
######################################################################

######################################################################
#
# sim.ril
#
# Simulate RILs by selfing or sibling mating from 2, 4, or 8
# parental strains
# map = map in the usual R/qtl map format
# m = interference parameter (0 is no interference)
######################################################################
sim.ril <-
function(map, n.ril=1, type=c("sibmating", "selfing"), n.str=c("2","4","8"),
         m=0, p=0, random.cross=FALSE)
{
  type <- match.arg(type)
  if(type=="sibmating") selfing <- 0
  else selfing <- 1
  n.str <- as.numeric(match.arg(n.str))
  n.chr <- length(map)
  n.mar <- sapply(map,length)
  tot.mar <- sum(n.mar)

  if(m < 0) stop("Must have m >= 0.")
  if(p < 0 || p > 1) stop("Must have 0 <= p <= 1.")
  if(p == 1) {
    p <- 0
    m <- 0
  }

  omap <- map
  map <- lapply(map, function(a) a-min(a))

  if(!selfing && class(omap[[length(omap)]])=="X")
    include.x <- TRUE
  else include.x <- FALSE

  x <- .C("R_sim_ril",
          as.integer(n.chr),
          as.integer(n.mar),
          as.integer(n.ril),
          as.double(unlist(map)),
          as.integer(n.str),
          as.integer(m),
          as.double(p),
          as.integer(include.x),
          as.integer(random.cross),
          as.integer(selfing),
          cross=as.integer(rep(0,n.ril*n.str)),
          res=as.integer(rep(0,tot.mar*n.ril)),
          PACKAGE="qtl")

  cross <- t(matrix(x$cross,ncol=n.ril,nrow=n.str))
  x <- t(matrix(x$res,nrow=tot.mar,ncol=n.ril))

  geno <- vector("list", n.chr)
  names(geno) <- names(map)
  cur <- 0
  for(i in 1:n.chr) {
    geno[[i]]$data <- x[,cur + 1:n.mar[i],drop=FALSE]
    colnames(geno[[i]]$data) <- names(map[[i]])
    geno[[i]]$map <- omap[[i]]
    cur <- cur + n.mar[i]
    class(geno[[i]]) <- class(map[[i]])
  }
  pheno <- data.frame(line=1:n.ril)
  x <- list(geno=geno,pheno=pheno,cross=cross)

  if(type=="sibmating") {
    if(n.str=="2")
      class(x) <- c("risib","cross")
    else
      class(x) <- c(paste("ri", n.str, "sib",sep=""),"cross")
  }
  else {
    if(n.str=="2")
      class(x) <- c("riself","cross")
    else
      class(x) <- c(paste("ri", n.str, "self",sep=""),"cross")
  }
  
  x
}

  
######################################################################
# sim.cc: Simulate the collaborative cross
#
# parents = Parental snp data, with genetic map
#           list with elements being chromosomes
#           each chromosome is a list with data=matrix n_mar x 8, 
#              map = vector of marker positions
#
# n.ril = number of lines to simulate
#
# error_prob = probability of genotyping error
# missing_prob = probability a genotype is missing
#
# m = interference parameter (0 is no interference)
#
# step = step size for intermediate loci to be simulated
######################################################################
sim.cc <-
function(parents, n.ril=1, type=c("sibmating", "selfing"),
         error.prob=0, missing.prob=0, m=0, p=0, step=0)
{
  type <- match.arg(type)
  map <- lapply(parents, function(a) a$map)
  markers <- vector("list",length(map))
  if(step<1e-8) {
    fmap <- map
    for(i in 1:length(map))
      markers[[i]] <- rep(TRUE,length(map[[i]]))
  }
  else {
    fmap <- vector("list",map)
    for(i in 1:length(map)) {
      fmap[[i]] <- create.map(map[[i]],step,0)
      class(fmap[[i]]) <- class(map[[i]])
      markers[[i]] <- map[[i]] %in% fmap[[i]]
      if(sum(markers[[i]]) != length(map[[i]]))
        warning("problem: screw up regarding create_map")
    }
  }

  if(m < 0) stop("Must have m >= 0.")
  if(p < 0 || p > 1) stop("Must have 0 <= p <= 1.")
  if(p == 1) {
    p <- 0
    m <- 0
  }

  cc <- sim.ril(fmap, n.ril, type, "8", m, p, TRUE)
  cc$truth <- cc$geno
  g <- pull.geno(cc)[,unlist(markers)]
  pg <- NULL
  for(i in 1:length(parents))
    pg <- rbind(pg,parents[[i]]$data)

  res <- .C("R_sim_cc",
            as.integer(n.ril), # no. ril
            as.integer(ncol(g)), # no. markers
            as.integer(pg), # SNP data on parents
            g=as.integer(g), # genotype data on rils
            as.double(error.prob), # error prob
            as.double(missing.prob), # missing data prob
            PACKAGE="qtl")$g

  n.mar <- sapply(map,length)

  g <- matrix(res,nrow=n.ril)
  
  # function for picking out the locations of breakpoints in the "truth"
  tempf <-
    function(a,b)  
      { 
        wh <- which(diff(a) != 0)
        x <- a[c(1,wh+1,length(a))]
        y <- c(b[1],(b[wh] + b[wh+1])/2,b[length(b)])
        y <- rbind(x,y)
        colnames(y) <- NULL
        y
      }

  # get genotype data back into cross
  cur <- 0
  for(i in seq(along=n.mar)) {
    cc$geno[[i]]$data <- g[,cur+(1:n.mar[i])]
    cc$geno[[i]]$map <- map[[i]]
    colnames(cc$geno[[i]]$data) <- names(map[[i]])
    cur <- cur + n.mar[i]

    # turn "truth" into a more compact form
    cc$truth[[i]] <- apply(cc$truth[[i]]$data,1,tempf,cc$truth[[i]]$map)
  }

  class(cc) <- c("cc","cross")
  cc
}

# end of sim_ril.R
