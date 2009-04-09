#####################################################################
#
# plotperm.R
#
# copyright (c) 2007-8, Karl W Broman
# last modified Jan, 2008
# first written Dec, 2007
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
# Contains: plot.scanoneperm, plot.scantwoperm, plot.scanoneboot 
#
######################################################################

############################################################
# plot.scanoneboot
#
# plot histogram of the results from scanoneboot
############################################################
plot.scanoneboot <-
function(x, ...)
{
  results <- attr(x, "results")

  markerpos <- results[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", rownames(results)),2]
  n.brk <- ceiling(2*sqrt(length(x)))
  xlim <- range(results[,2])

  dots <- list(...)
  main <- ""
  xlab <- "QTL position (cM)"

  if("breaks" %in% names(dots)) {
    if("main" %in% names(dots)) {
      if("xlab" %in% names(dots)) {
        hist(x, xlim=xlim, ...)
      }
      else {
        hist(x, xlim=xlim, xlab=xlab, ...)
      }
    }
    else {
      if("xlab" %in% names(dots)) {
        hist(x, xlim=xlim, main=main, ...)
      }
      else {
        hist(x, xlim=xlim, xlab=xlab, main=main, ...)
      }
    }
  }
  else {
    if("main" %in% names(dots)) {
      if("xlab" %in% names(dots)) {
        hist(x, breaks=n.brk, xlim=xlim, ...)
      }
      else {
        hist(x, breaks=n.brk, xlim=xlim, xlab=xlab, ...)
      }
    }
    else {
      if("xlab" %in% names(dots)) {
        hist(x, breaks=n.brk, xlim=xlim, main=main, ...)
      }
      else {
        hist(x, breaks=n.brk, xlim=xlim, xlab=xlab, main=main, ...)
      }
    }
  }
  
  rug(markerpos)
}

############################################################
# plot.scanoneperm
#
# plot histogram of the permutation results from scanone
############################################################
plot.scanoneperm <-
function(x, lodcolumn=1, ...)
{
  if(is.list(x)) { # separate X chr results
    if(lodcolumn < 1 || lodcolumn > ncol(x$A))
      stop("lodcolumn should be between 1 and ", ncol(x$A))
    A <- x$A[,lodcolumn]
    X <- x$X[,lodcolumn]
    
    n.brkA <- ceiling(2*sqrt(length(A)))
    n.brkX <- ceiling(2*sqrt(length(X)))
    xlim <- c(0,max(c(A,X)))

    dots <- list(...)
    xlab <- "maximum LOD score"

    mfrow <- par("mfrow")
    on.exit(par(mfrow=mfrow))
    par(mfrow=c(2,1))

    if("breaks" %in% names(dots)) {
      if("xlab" %in% names(dots)) {
        hist(A, xlim=xlim, main="Autosomes", ...)
        rug(A)
        hist(X, xlim=xlim, main="X chromosome", ...)
        rug(X)
      }
      else {
        hist(A, xlim=xlim, xlab=xlab, main="Autosomes", ...)
        rug(A)
        hist(X, xlim=xlim, xlab=xlab, main="X chromosome", ...)
        rug(X)
      }
    }
    else {
      if("xlab" %in% names(dots)) {
        hist(A, breaks=n.brkA, xlim=xlim, main="Autosomes", ...)
        rug(A)
        hist(X, breaks=n.brkX, xlim=xlim, main="X chromosome", ...)
        rug(X)
      }
      else {
        hist(A, breaks=n.brkA, xlim=xlim, xlab=xlab, main="Autosomes", ...)
        rug(A)
        hist(X, breaks=n.brkX, xlim=xlim, xlab=xlab, main="X chromosome", ...)
        rug(X)
      }
    }
  }
  else {  
    if(lodcolumn < 1 || lodcolumn > ncol(x))
      stop("lodcolumn should be between 1 and ", ncol(x$A))
    x <- x[,lodcolumn]
    n.brk <- ceiling(2*sqrt(length(x)))
    xlim <- c(0,max(as.numeric(x)))

    dots <- list(...)
    main <- ""
    xlab <- "maximum LOD score"

    if("breaks" %in% names(dots)) {
      if("main" %in% names(dots)) {
        if("xlab" %in% names(dots)) {
          hist(x, xlim=xlim, ...)
        }
        else {
          hist(x, xlim=xlim, xlab=xlab, ...)
        }
      }
      else {
        if("xlab" %in% names(dots)) {
          hist(x, xlim=xlim, main=main, ...)
        }
        else {
          hist(x, xlim=xlim, xlab=xlab, main=main, ...)
        }
      }
    }
    else {
      if("main" %in% names(dots)) {
        if("xlab" %in% names(dots)) {
          hist(x, breaks=n.brk, xlim=xlim, ...)
        }
        else {
          hist(x, breaks=n.brk, xlim=xlim, xlab=xlab, ...)
        }
      }
      else {
        if("xlab" %in% names(dots)) {
          hist(x, breaks=n.brk, xlim=xlim, main=main, ...)
        }
        else {
          hist(x, breaks=n.brk, xlim=xlim, xlab=xlab, main=main, ...)
        }
      }
    }
  
    rug(as.numeric(x))
  }
}


############################################################
# plot.scantwoperm
#
# plot histogram of the permutation results from scantwo
############################################################
plot.scantwoperm <-
function(x, lodcolumn=1, ...)
{
  if(lodcolumn < 1 || lodcolumn > ncol(x[[1]]))
    stop("lodcolumn should be between 1 and ", ncol(x[[1]]))
  x <- lapply(x, function(a,b) a[,b], lodcolumn)

  xlim <- c(0,max(unlist(x)))
  brks <- seq(0, max(unlist(x)), length=ceiling(4*sqrt(length(x[[1]])))+1)

  dots <- list(...)
  xlab <- "maximum LOD score"

  mfcol <- par("mfcol")
  on.exit(par(mfcol=mfcol))
  par(mfcol=c(3,2))

  if("breaks" %in% names(dots)) {
    if("xlab" %in% names(dots)) {
      for(i in seq(along=x)) {
        hist(x[[i]], xlim=xlim, main=names(x)[i], ...)
        rug(x[[i]])
      }
    }
    else {
      for(i in seq(along=x)) {
        hist(x[[i]], xlim=xlim, xlab=xlab, main=names(x)[i], ...)
        rug(x[[i]])
      }
    }
  }
  else {
    if("xlab" %in% names(dots)) {
      for(i in seq(along=x)) {
        hist(x[[i]], breaks=brks, xlim=xlim, main=names(x)[i], ...)
        rug(x[[i]])
      }
    }
    else {
      for(i in seq(along=x)) {
        hist(x[[i]], breaks=brks, xlim=xlim, xlab=xlab, main=names(x)[i], ...)
        rug(x[[i]])
      }
    }
  }
}

                                
# end of plotperm.R
