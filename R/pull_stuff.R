#####################################################################
#
# pull_stuff.R
#
# copyright (c) 2001-2011, Karl W Broman
#     [find.pheno, find.flanking, and a modification to create.map
#      from Brian Yandell]
# last modified Dec, 2011
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
# Contains: pull.map, pull.geno, pull.pheno
#
######################################################################

######################################################################
#
# pull.map
#
# pull out the map portion of a cross object, as a list
#
######################################################################

pull.map <-
function(cross, chr, as.table=FALSE)
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  if(!missing(chr)) cross <- subset(cross, chr=chr)
  if(!as.table) {
    a <- lapply(cross$geno,function(a) {
      b <- a$map
      class(b) <- as.character(class(a))
      b })
    class(a) <- "map"
  } else {
    themap <- pull.map(cross, as.table=FALSE)
    if(is.matrix(themap[[1]])) {
      themap1 <- unlist(lapply(themap, function(a) a[1,]))
      themap2 <- unlist(lapply(themap, function(a) a[2,]))
      a <- data.frame(chr=rep(names(cross$geno), nmar(cross)),
                      pos.female=themap1, pos.male=themap2, stringsAsFactors=TRUE)
    } else {
      a <- data.frame(chr=rep(names(cross$geno), nmar(cross)),
                      pos=unlist(themap), stringsAsFactors=TRUE)
    }
    rownames(a) <- markernames(cross)
  }

  a
}

######################################################################
# pull.geno
######################################################################
pull.geno <-
function(cross, chr)
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  if(!missing(chr))
    cross <- subset(cross, chr=chr)
  
  X <- cross$geno[[1]]$data
  if(nchr(cross) > 1)
    for(i in 2:nchr(cross))
      X <- cbind(X, cross$geno[[i]]$data)
  X
}

######################################################################
# pull.pheno
######################################################################
pull.pheno <-
function(cross, pheno.col)
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  pheno <- cross$pheno

  if(!missing(pheno.col)) {

    if(is.character(pheno.col)) {
      m <- match(pheno.col, names(pheno))
      if(any(is.na(m))) {
        if(sum(is.na(m)) > 1)
          warning("Phenotypes ", paste("\"", pheno.col[is.na(m)], "\"", sep="", collapse=" "), " not found.")
        else 
          warning("Phenotype ", paste("\"", pheno.col[is.na(m)], "\"", sep="", collapse=" "), " not found.")
      }
      if(all(is.na(m))) return(NULL)
      
      m <- m[!is.na(m)]
      pheno <- pheno[,m]
    }
    else if(is.logical(pheno.col)) {
      if(length(pheno.col) != ncol(pheno))
        stop("If pheno.col is logical, it should have length ", ncol(pheno))
      pheno <- pheno[,pheno.col]
    }
    else if(is.numeric(pheno.col)) {
      if(any(pheno.col > 0) && any(pheno.col < 0))
        stop("If pheno.col is numeric, values should be all > 0 or all < 0")
      if(any(pheno.col > 0) && (any(pheno.col < 1) || any(pheno.col > ncol(pheno))))
        stop("pheno.col values should be >= 1 and <= ", ncol(pheno))
      if(any(pheno.col < 0) && (any(pheno.col > -1) || any(pheno.col < -ncol(pheno))))
        stop("With negative pheno.col values, they should be between -", ncol(pheno), " and -1")
      pheno <- pheno[,pheno.col]
    }
  }

  if(is.data.frame(pheno) && ncol(pheno) == 1) pheno <- pheno[,1]

  pheno
}

# end of pull_stuff.R
