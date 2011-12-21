#####################################################################
#
# inferFounderHap.R
#
# copyright (c) 2011, Karl W Broman
# last modified Dec, 2011
# first written Dec, 2011
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
# Contains: inferFounderHap
#
# This is for reconstructing the founder haplotypes in inbred lines
# by a crude method using groups of adjacent SNPs
#
######################################################################

inferFounderHap <-
function(offspringGen, founderGen, max.n.marker=10, verbose=FALSE)
{
  n.mar <- ncol(founderGen)
  n.ind <- nrow(offspringGen)
  n.founders <- nrow(founderGen)
  if(n.mar != ncol(founderGen))
    stop("ncol(offspringGen) != ncol(founderGen)")
  if(any(!is.na(offspringGen) & offspringGen != 0 & offspringGen != 1))
    stop("offspringGen should be NA, 0 or 1")
  if(any(!is.na(founderGen) & founderGen != 0 & founderGen != 1))
    stop("founderGen should be NA, 0 or 1")
  
  nomissing <- apply(founderGen, 2, function(a) !any(is.na(a)))
  if(!any(nomissing))
    stop("No markers with complete founder genotypes")
  offspringGen <- offspringGen[,nomissing,drop=FALSE]
  founderGen <- founderGen[,nomissing,drop=FALSE]

  if(max.n.marker > n.mar) max.n.marker <- n.mar

  offspringGen[is.na(offspringGen)] <- -1

  z <- .C("R_inferFounderHap",
          as.integer(n.mar),
          as.integer(n.founders),
          as.integer(n.ind),
          as.integer(founderGen),
          as.integer(offspringGen),
          as.integer(max.n.marker),
          hap=as.integer(rep(0,n.mar * n.ind)),
          as.integer(verbose),
          PACKAGE="qtl")
  hap <- matrix(z$hap, ncol=n.mar, nrow=n.ind)
  hap[hap <= 0] <- NA

  hap
}


# end of inferFounderHap.R
