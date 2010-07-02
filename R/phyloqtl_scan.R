######################################################################
# phyloqtl_scan.R
#
# copyright (c) 2009-2010, Karl W Broman
# last modified Jul, 2010
# first written May, 2009
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
# Single-QTL scan to map QTL to a phylogenetic tree
#
# Part of the R/qtl package
# Contains: scanPhyloQTL, max.scanPhyloQTL, plot.scanPhyloQTL
#
######################################################################

scanPhyloQTL <-
function(crosses, partitions, chr, pheno.col=1,
         model=c("normal","binary"), method=c("em","imp","hk"),
         maxit=4000, tol=0.0001)
{  
  if(missing(chr)) chr <- names(crosses[[1]]$geno)
  model <- match.arg(model)
  method <- match.arg(method)

  thecrosses <- names(crosses)
  taxa <- sort(unique(unlist(strsplit(thecrosses, ""))))

  thecrosses <- names(crosses) <- checkPhyloCrosses(thecrosses, taxa)

  if(missing(partitions)) {
    # generate all partitions
    partitions <- genAllPartitions(length(taxa), taxa)
  }
  checkPhyloPartition(partitions, taxa)

  crossmat <- qtlByPartition(thecrosses, partitions)

  dimnames(crossmat) <- list(thecrosses, partitions)

  out <- vector("list", length(partitions))
  for(i in seq(along=partitions)) {
    cm <- crossmat[,i]

    # if all crosses need to be flipped, don't flip any of them
    if(!any(cm > 0)) cm <- -cm

    # pull out the crosses to use here
    x <- crosses[cm != 0]
    cm <- cm[cm != 0]

    # flip crosses if necessary
    if(any(cm < 0))
      for(j in which(cm < 0)) 
        x[[j]] <- flipcross(x[[j]])

    # combine the crosses
    xx <- x[[1]]
    n.phe <- nphe(x[[1]])
    if(length(x) > 1) {
      for(j in 2:length(x)) {
        xx <- c(xx, x[[j]])
        xx$pheno <- xx$pheno[,1:n.phe,drop=FALSE]
      }
    }

    # create cross indicators (as additive covariates)
    if(length(x)==1) addcovar<-NULL
    else {
      ni <- sapply(x, nind)
      addcovar <- matrix(0, ncol=length(x)-1, nrow=sum(ni))
      thestart <- cumsum(c(1,ni))
      end <- cumsum(ni)
      for(j in 2:length(x))
        addcovar[thestart[j]:end[j],j-1] <- 1
    }

    # do the scan
    out[[i]] <- scanone(xx, chr=chr, pheno.col=pheno.col, addcovar=addcovar,
                        model=model, method=method, maxit=maxit, tol=tol)
  }

  # just one partition
  if(length(out) == 1) return(out[[1]])

  # multiple partitions
  result <- out[[1]]
  for(j in 2:length(out)) 
    result <- c(result, out[[j]])
  colnames(result)[-(1:2)] <- partitions
  class(result) <- c("scanPhyloQTL", "scanone", "data.frame")
  result
}


max.scanPhyloQTL <-
function(object, format=c("lod", "postprob"), chr, na.rm=TRUE, ...)
{
  format <- match.arg(format)

  if(!missing(chr))
    object <- subset(object, chr=chr)

  # find the overall maximum
  themax <- apply(object[,-(1:2)], 2, max)
  wh <- which(themax==max(themax))
  if(length(wh) > 1) wh <- sample(wh, 1)

  # use max.scanone to create a filler for the output
  mx <- max.scanone(object, lodcolumn=wh, na.rm=TRUE)

  # chromosome where max occurred
  thechr <- as.character(mx$chr)

  # maximum LOD for each partition on that chromosome
  themax <- apply(object[object$chr==thechr,-(1:2)], 2, max)

  mx[-(1:2)] <- themax
  
  thelod <- -diff(sort(themax, decreasing=TRUE)[1:2])
  names(thelod) <- names(themax)[wh]
  attr(mx, "LOD") <- thelod

  if(format=="lod") 
    mx <- cbind(mx, inferred=colnames(object)[wh+2], loddif=thelod)
  else {
    mxmx <- apply(mx[,-(1:2)], 1, max)
    temp <- mx[,-(1:2)]
    mx[,-(1:2)] <- t(apply(temp, 1, function(a) 10^a/sum(10^a)))
    mx <- cbind(mx, inferred=colnames(object)[wh+2], maxlod=mxmx)
  }

  class(mx) <- c("summary.scanPhyloQTL", "summary.scanone", "data.frame")

  mx
}

summary.scanPhyloQTL <-
function(object, format=c("lod", "postprob"), threshold, na.rm=TRUE) 
{
  format <- match.arg(format)
  themax <- apply(object[,-(1:2)], 2, tapply, object[,1], max, na.rm=na.rm)
  wh <- apply(themax, 1, function(a) { a <- which(a==max(a)); if(length(a) > 1) a <- sample(a, 1); a })
  whpos <- rep(NA, length(wh))
  names(whpos) <- unique(object[,1])
  for(i in seq(along=whpos)) {
    temp <- object[object[,1]==names(whpos)[i],,drop=FALSE]
    temp2 <- which(temp[,wh[i]+2]==max(temp[,wh[i]+2],na.rm=na.rm))
    if(length(temp2) > 1) temp2 <- sample(temp2, 1)
    whpos[i] <- temp[temp2,2]
    names(whpos)[i] <- rownames(temp)[temp2]
  }
  if(format=="lod") {
    out <- data.frame(chr=unique(object[,1]), pos=whpos, themax,
                      inferred=colnames(object)[wh+2],
                      loddif=apply(themax, 1, function(a) -diff(sort(a, decreasing=TRUE)[1:2])))
  }
  else {
    out <- data.frame(chr=unique(object[,1]), pos=whpos, themax,
                      inferred=colnames(object)[wh+2],
                      maxlod=apply(themax, 1, max))
    temp <- out[,-c(1:2, ncol(out)-0:1)]
    out[,-c(1:2, ncol(out)-0:1)] <- t(apply(temp, 1, function(a) 10^a/sum(10^a)))
  }

  rownames(out) <- names(whpos)

  if(!missing(threshold)) 
    out <- out[out[ncol(out)] >= threshold,,drop=FALSE]

  class(out) <- c("summary.scanPhyloQTL", "summary.scanone", "data.frame")
  out
}

print.summary.scanPhyloQTL <- print.summary.scanone


# plot results of scanPhyloQTL
plot.scanPhyloQTL <-
function(x, chr, incl.markers=TRUE, col, xlim, ylim, lwd=2,
         gap=25, mtick=c("line", "triangle"),
         show.marker.names=FALSE, alternate.chrid=FALSE,
         legend=TRUE, ...)
{
  mtick <- match.arg(mtick)

  if(!missing(chr)) x <- subset(x, chr=chr)
  
  if(missing(col)) {
    col <- c("black","blue","red","green","orange","brown","gray","cyan","magenta")
    if(ncol(x)-2 > length(col))
      stop("Please give a list of colors")
    col <- col[1:(ncol(x)-2)]
  }

  if(missing(ylim))
    ylim <- c(0, max(apply(x[,-(1:2)], 2, max)))

  dots <- list(...)

  if(missing(xlim)) {
    if("ylab" %in% names(dots)) 
      plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                   ylim=ylim, lwd=lwd, gap=gap, mtick=mtick,
                   show.marker.names=show.marker.names, alternate.chrid=alternate.chrid, ...)
    else
      plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                   ylim=ylim, lwd=lwd, gap=gap, mtick=mtick,
                   show.marker.names=show.marker.names, alternate.chrid=alternate.chrid,
                   ylab="LOD score", ...)
  }
  else {
    if("ylab" %in% names(dots)) 
      plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                   ylim=ylim, xlim=xlim, lwd=lwd, gap=gap, mtick=mtick,
                   show.marker.names=show.marker.names, alternate.chrid=alternate.chrid, ...)
    else
      plot.scanone(x, incl.markers=incl.markers, col=col[1], lodcolumn=1,
                   ylim=ylim, xlim=xlim, lwd=lwd, gap=gap, mtick=mtick,
                   show.marker.names=show.marker.names, alternate.chrid=alternate.chrid,
                   ylab="LOD score", ...)
  }
    
  if(ncol(x) > 3)
    for(i in 2:(ncol(x)-2)) 
      plot.scanone(x, col=col[i], lodcolumn=i, add=TRUE, ...)

  if(is.character(legend) || legend) {
    if(is.character(legend))
      legend(legend, legend=colnames(x)[-(1:2)], col=col, lwd=lwd)
    else
      legend("topright", legend=colnames(x)[-(1:2)], col=col, lwd=lwd)
  }

  invisible()
}
    
  


# end of phyloqtl_scan.R
