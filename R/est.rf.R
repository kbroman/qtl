######################################################################
#
# est.rf.R
#
# copyright (c) 2001-9, Karl W Broman
# last modified Apr, 2009
# first written Apr, 2001
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
# Contains: est.rf, plot.rf, checkAlleles
#
######################################################################

######################################################################
#
# est.rf: Estimate sex-averaged recombination fractions between
#         all pairs of markers
#
######################################################################

est.rf <-
function(cross, maxit=10000, tol=1e-6) 
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  n.chr <- nchr(cross)
  n.mar <- totmar(cross)
  n.ind <- nind(cross)
  mar.names <- unlist(lapply(cross$geno,function(a) colnames(a$data)))
  
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno,class)

  xchrcol <- NULL
  fixX <- FALSE
  Geno <- NULL
  # create full genotype matrix
  for(i in 1:n.chr) {
    temp <- cross$geno[[i]]$data

    # treat X chromosome specially in an intercross
    if(type=="f2" && chrtype[i]=="X") {
      fixX <- TRUE
      if(i != 1) xchrcol <- c(xchrcol,ncol(Geno)+(1:ncol(cross$geno[[i]]$data)))
      else xchrcol <- 1:ncol(cross$geno[[i]]$data)
      xchr <- temp
      xchr[is.na(xchr)] <- 0
      temp <- reviseXdata(type,"simple",getsex(cross),geno=temp,
                          cross.attr=attributes(cross))
    }
    Geno <- cbind(Geno,temp)
  }

  # which type of cross is this?
  if(type == "f2")
    cfunc <- "est_rf_f2"
  else if(type == "bc" || type=="risib" || type=="riself" || type=="dh") 
    cfunc <- "est_rf_bc"
  else if(type == "4way") 
    cfunc <- "est_rf_4way"
  else if(type=="ri8sib" || type=="ri8self" || type=="ri4sib" || type=="ri4self") {
    cfunc <- paste("est_rf_", type, sep="")
    if(any(chrtype == "X"))
      warning("est.rf not working properly for the X chromosome for 4- or 8-way RIL.")
  }
  else 
    stop("est.rf not available for cross type ", type, ".")

  Geno[is.na(Geno)] <- 0
  
  if(type=="bc" || type=="risib" || type=="riself" || type=="dh")
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            rf = as.double(rep(0,n.mar*n.mar)),
            PACKAGE="qtl")
  else
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            rf = as.double(rep(0,n.mar*n.mar)),
            as.integer(maxit),
            as.double(tol),
            PACKAGE="qtl")

  cross$rf <- matrix(z$rf,ncol=n.mar)
  dimnames(cross$rf) <- list(mar.names,mar.names)

  if(fixX) {
    zz <- .C("est_rf_bc",
             as.integer(n.ind),
             as.integer(ncol(xchr)),
             as.integer(xchr),
             rf=as.double(rep(0,ncol(xchr)^2)),
             PACKAGE="qtl")
    zz <- matrix(zz$rf,ncol=ncol(xchr))
    cross$rf[xchrcol,xchrcol] <- zz
  }

  # check for alleles switches
  if(type == "risib" || type=="riself" || type=="f2" || type=="bc" || type=="dh") {
    out <- checkAlleles(cross, 5, FALSE)
    if(!is.null(out)) {
      out <- as.character(out[,1])
      warning("Alleles potentially switched at markers \n  ",
              paste(out, collapse=" "))
    }
  }

  cross
}

  

plot.rf <-
function(x, chr, what=c("both","lod","rf"),
         alternate.chrid=FALSE, zmax=12,
         mark.diagonal=FALSE, ...)
{
  if(!any(class(x) == "cross"))
    stop("Input should have class \"cross\".")

  what <- match.arg(what)
  
  if(!missing(chr)) x <- subset(x,chr=chr)
  
  if(!("rf" %in% names(x))) {
    warning("Running est.rf.")
    x <- est.rf(x)
  }
  g <- x$rf
  
  old.xpd <- par("xpd")
  old.las <- par("las")
  par(xpd=TRUE,las=1)
  on.exit(par(xpd=old.xpd,las=old.las))

  # if any of the rf's are NA (ie no data), put NAs in corresponding LODs
  if(any(is.na(g))) g[is.na(t(g))] <- NA

  # convert rf to -2*(log2(rf)+1); place zmax's on the diagonal;
  #    anything above zmax replaced by zmax;
  #    NA's replaced by -1
  g[row(g) > col(g) & g > 0.5] <- 0.5
  g[row(g) > col(g)] <- -4*(log2(g[row(g) > col(g)])+1)/12*zmax
  diag(g) <- zmax
  g[!is.na(g) & g>zmax] <- zmax
  
  g[is.na(g)] <- -1

  if(what=="lod") { # plot LOD scores 
    # copy upper triangle (LODs) to lower triangle (rec fracs)
    g[row(g) > col(g)] <- t(g)[row(g) > col(g)]
  }
  else if(what=="rf") { # plot recombination fractions
    # copy lower triangle (rec fracs) to upper triangle (LODs)
    g[row(g) < col(g)] <- t(g)[row(g) < col(g)]
  }
  br <- c(-1, seq(-1e-6, zmax, length=257))

  image(1:ncol(g),1:nrow(g),t(g),ylab="Markers",xlab="Markers",breaks=br,
        col=c("lightgray",rev(rainbow(256,start=0,end=2/3, gamma=0.6))))
  
  if(mark.diagonal) {
    for(i in 1:ncol(g))
      segments(i+c(-0.5, -0.5, -0.5, +0.5), i+c(-0.5, +0.5, -0.5, -0.5),
               i+c(-0.5, +0.5, +0.5, +0.5), i+c(+0.5, +0.5, -0.5, +0.5))
  }

  # plot lines at the chromosome boundaries
  n.mar <- nmar(x)
  n.chr <- nchr(x)
  a <- c(0.5,cumsum(n.mar)+0.5)
  abline(v=a,xpd=FALSE,col="white")
  abline(h=a,xpd=FALSE,col="white")

  # this line adds a line above the image
  #     (the image function leaves it out)
  abline(h=0.5+c(0,nrow(g)),xpd=FALSE)
  abline(v=0.5+c(0,nrow(g)),xpd=FALSE)

  # add chromosome numbers
  a <- par("usr")
  wh <- cumsum(c(0.5,n.mar))
  chrnam <- names(x$geno)
  chrpos <- (wh[-1] + wh[-length(wh)])/2
  if(!alternate.chrid || length(chrnam) < 2) {
    for(i in seq(along=chrpos)) {
      axis(side=3, at=chrpos[i], labels=chrnam[i], tick=FALSE, line=-0.8)
      axis(side=4, at=chrpos[i], labels=chrnam[i], tick=FALSE, line=-0.8)
    }
  }
  else {
    odd <- seq(1, length(chrpos), by=2)
    even <- seq(2, length(chrpos), by=2)
    for(i in odd) {
      axis(side=3, at=chrpos[i], labels=chrnam[i], line=-0.8, tick=FALSE)
      axis(side=4, at=chrpos[i], labels=chrnam[i], line=-0.8, tick=FALSE)
    }
    for(i in even) {
      axis(side=3, at=chrpos[i], labels=chrnam[i], line=0, tick=FALSE)
      axis(side=4, at=chrpos[i], labels=chrnam[i], line=0, tick=FALSE)
    }

  }


  dots <- list(...)
  if("main" %in% names(dots))
    title(main=dots$main)
  else {
    if(what=="lod") title(main="Pairwise LOD scores")
    else if(what=="rf") title(main="Recombination fractions")
    else title("Pairwise recombination fractions and LOD scores")
  }

  invisible()
}

######################################################################
# check for apparent errors in the recombination fractions
######################################################################
#checkrf <-
#function(cross, threshold=5)
#{
#  rf <- cross$rf
#  n.mar <- nmar(cross)
#  map <- pull.map(cross)
#  n <- ncol(rf)
#  mnam <- colnames(rf)
#  whpos <- unlist(lapply(map,function(a) 1:length(a)))
#  whchr <- rep(names(map),sapply(map,length))
#
#  # first check whether a locus has "significant" pairwise recombination
#  #     with rf > 0.5
#  for(i in 1:n) {
#    if(i == 1) {
#      lod <- rf[1,-1]
#      r <- rf[-1,1]
#    }
#    else if(i == n) {
#      lod <- rf[-n,n]
#      r <- rf[n,-n]
#    }
#    else {
#      lod <- c(rf[1:(i-1),i],rf[i,(i+1):n])
#      r <- c(rf[i,1:(i-1)],rf[(i+1):n,i])
#    }
#
#    # if rf > 1/2 and LOD > threshold for more than two other markers
#    if(sum(!is.na(lod) & !is.na(r) & lod > threshold & r > 0.5) >= 2)
#      warning("Genotypes potentially switched for marker ", mnam[i],
#          paste(" (",whpos[i],")",sep=""), " on chr ", whchr[i], "\n")
#    
#  }
#
#}


######################################################################
# checkAlleles()
#
# Function to find markers that may have alleles miscoded;
# we go through each marker, one at a time, swap alleles and
# then see what it does to pairwise linkage against all other
# markers
######################################################################

checkAlleles <-
function(cross, threshold=3, verbose=TRUE)
{
  if(!any(class(cross) == "cross"))
    stop("checkAlleles() only works for cross objects.")

  type <- class(cross)[1]
  if(type != "f2" && type != "bc" &&
     type != "risib" && type != "riself" && type != "dh")
    stop("checkAlleles not available for cross type ", type, ".")
    
  # drop X chromosome
  chrtype <- sapply(cross$geno,class)
  if(all(chrtype=="X")) {
    if(verbose) cat("checkAlleles() only works for autosomal data.\n")
    return(NULL)
  }

  cross <- subset(cross, chr = (chrtype != "X"))

  n.mar <- nmar(cross)
  mar.names <- unlist(lapply(cross$geno,function(a) colnames(a$data)))

  if(!("rf" %in% names(cross))) {
    warning("First running est.rf.")
    cross <- est.rf(cross) 
  }
  diag(cross$rf) <- 0
  lod <- rf <- cross$rf
  lod[lower.tri(lod)] <- t(lod)[lower.tri(lod)]
  rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]

  orig.lod <- rev.lod <- lod
  orig.lod[rf > 0.5] <- 0
  rev.lod[rf < 0.5] <- 0

  dif <- apply(rev.lod, 2, max, na.rm=TRUE) -
    apply(orig.lod, 2, max, na.rm=TRUE)

  results <- data.frame(marker=mar.names,
                        chr=rep(names(cross$geno), n.mar),
                        index=unlist(lapply(n.mar, function(a) 1:a)),
                        "diff in max LOD" = dif)
  rownames(results) <- 1:nrow(results)

  if(all(results[,4] < threshold)) {
    if(verbose) cat("No apparent problems.\n")
    return(invisible(NULL))
  }
  
  results[results[,4] >= threshold,]
}

# end of est.rf.R
