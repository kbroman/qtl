#####################################################################
#
# arithscan.R
#
# copyright (c) 2005-8, Karl W Broman
# last modified Apr, 2008
# first written Mar, 2005
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
# Contains: +.scanone, -.scanone, +.scanoneperm, -.scanoneperm
#           +.scantwo, -.scantwo, +.scantwoperm, -.scantwoperm
#
######################################################################


"-.scanone" <-
function(e1,e2)
{
  if(!any(class(e1) == "scanone"))
    stop("Input should have class \"scanone\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if(missing(e2)) {
    class(e1) <- "data.frame"
    e1[,-(1:2)] <- -e1[,-(1:2)]
    class(e1) <- c("scanone","data.frame")
    if(!is.null(df1))
      attr(e1, "df") <- -df1
    
    return(e1)
  }
  if(!any(class(e2) == "scanone"))
    stop("Input should have class \"scanone\".")

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  class(e1) <- class(e2) <- "data.frame"
  if(nrow(e1) != nrow(e2)) {
    u1 <- levels(e1[,1])
    u2 <- levels(e2[,2])
    u <- unique(c(u1,u2))
    if(length(u) == 0) stop("Can't subtract; no chromosomes in common.")
    e1 <- e1[!is.na(match(e1[,1],u)),]
    e2 <- e2[!is.na(match(e2[,1],u)),]
    if(nrow(e1) != nrow(e2) || any(e1[,1] != e2[,1]) || max(abs(e1[,2]-e2[,2])) < 0.01)
      stop("Can't subtract; arguments not compatible")
  }
  nc1 <- ncol(e1)
  nc2 <- ncol(e2)
  nc <- min(c(nc1,nc2))
  e1 <- e1[,1:nc]
  e1[,3:nc] <- e1[,3:nc] - e2[,3:nc]

  # zero out small stuff
  temp <- e1[,3:nc]
  temp[!is.na(temp) & abs(temp) < 1e-6] <- 0
  e1[,3:nc] <- temp
  
  class(e1) <- c("scanone","data.frame")

  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 - df2
  }

  e1
}

"+.scanone" <-
function(e1,e2)
{
  if(!any(class(e1) == "scanone"))
    stop("Input should have class \"scanone\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if(missing(e2)) return(e1)
  if(!any(class(e2) == "scanone"))
    stop("Input should have class \"scanone\".")

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  class(e1) <- class(e2) <- "data.frame"
  if(nrow(e1) != nrow(e2)) {
    u1 <- levels(e1[,1])
    u2 <- levels(e2[,2])
    u <- unique(c(u1,u2))
    if(length(u) == 0) stop("Can't subtract; no chromosomes in common.")
    e1 <- e1[!is.na(match(e1[,1],u)),]
    e2 <- e2[!is.na(match(e2[,1],u)),]
    if(nrow(e1) != nrow(e2) || any(e1[,1] != e2[,1]) || max(abs(e1[,2]-e2[,2])) < 0.01)
      stop("Can't subtract; arguments not compatible")
  }
  nc1 <- ncol(e1)
  nc2 <- ncol(e2)
  nc <- min(c(nc1,nc2))
  e1 <- e1[,1:nc]
  e1[,3:nc] <- e1[,3:nc] + e2[,3:nc]
  class(e1) <- c("scanone","data.frame")

  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 + df2
  }

  e1
}

"-.scanoneperm" <-
function(e1, e2)
{
  if(!any(class(e1) == "scanoneperm"))
    stop("Input should have class \"scanoneperm\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if(missing(e2)) {
    e1.x <- ("xchr" %in% names(attributes(e1)))
    if(e1.x) {
      e1$A <- -e1$A
      e1$X <- -e1$X
    }
    else {
      theclass <- class(e1)
      e1 <- -unclass(e1)
      class(e1) <- theclass
    }

    if(!is.null(df1))
      attr(e1, "df") <- -df1

    return(e1)
  }
  if(!any(class(e2) == "scanoneperm"))
    stop("Input should have class \"scanoneperm\".")
  
  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  # check input
  e1.x <- ("xchr" %in% names(attributes(e1)))
  e2.x <- ("xchr" %in% names(attributes(e2)))

  if(e1.x != e2.x)
    stop("Need both or neither input to be X-chr specific.\n")

  if((e1.x && (any(dim(e1$A)!=dim(e2$A)) || any(dim(e1$X)!=dim(e2$X)))) ||
     (!e1.x && any(dim(e1) != dim(e2))))
    stop("Need input to concern the same phenotypes and no. permutations.\n")
    
  if(e1.x) {
    e1$A <- e1$A - e2$A
    e1$X <- e1$X - e2$X

    # zero out small stuff
    e1$A[!is.na(e1$A) & abs(e1$A) < 1e-6] <- 0
    e1$X[!is.na(e1$X) & abs(e1$X) < 1e-6] <- 0
  }
  else {
    theclass <- class(e1)
    e1 <- unclass(e1) - unclass(e2)

    # zero out small stuff
    e1[!is.na(e1) & abs(e1) < 1e-6] <- 0

    class(e1) <- theclass
  }

  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 - df2
  }

  e1
}

"+.scanoneperm" <-
function(e1, e2)
{
  if(!any(class(e1) == "scanoneperm"))
    stop("Input should have class \"scanoneperm\".")
  if(missing(e2)) return(e1)
  if(!any(class(e2) == "scanoneperm"))
    stop("Input should have class \"scanoneperm\".")
  
  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  # check input
  e1.x <- ("xchr" %in% names(attributes(e1)))
  e2.x <- ("xchr" %in% names(attributes(e2)))

  if(e1.x != e2.x)
    stop("Need both or neither input to be X-chr specific.\n")

  if((e1.x && (any(dim(e1$A)!=dim(e2$A)) || any(dim(e1$X)!=dim(e2$X)))) ||
     (!e1.x && any(dim(e1) != dim(e2))))
    stop("Need input to concern the same phenotypes and no. permutations.\n")
    
  if(e1.x) {
    e1$A <- e1$A + e2$A
    e1$X <- e1$X + e2$X
  }
  else {
    theclass <- class(e1)
    e1 <- unclass(e1) + unclass(e2)
    class(e1) <- theclass
  }

  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 + df2
  }

  e1
}

######################################################################
# -.scantwo: subtract LOD scores in two scantwo results
######################################################################
"-.scantwo" <-
function(e1, e2)  
{
  if(!any(class(e1) == "scantwo"))
    stop("Input should have class \"scantwo\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if(missing(e2)) {
    e1$lod <- -e1$lod
    if("scanoneX" %in% names(e1))
      e1$scanoneX <- -e1$scanoneX
    if(!is.null(df1))
      attr(e1, "df") <- -df1
    
    return(e1)
  }
  if(!any(class(e2) == "scantwo"))
    stop("Input should have class \"scantwo\".")

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  e1x <- "scanoneX" %in% names(e1)
  e2x <- "scanoneX" %in% names(e2)

  if(any(dim(e1$map) != dim(e2$map)) ||
     length(dim(e1$lod)) != length(dim(e2$lod)) ||
     any(dim(e1$lod) != dim(e2$lod)) || e1x != e2x)
    stop("input arguments do not conform.")

  e1$lod <- e1$lod - e2$lod
  if(e1x) e1$scanoneX <- e1$scanoneX - e2$scanoneX
  
  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 - df2
  }

  e1
}


######################################################################
# +.scantwo: add LOD scores in two scantwo results
######################################################################
"+.scantwo" <-
function(e1, e2)
{
  if(!any(class(e1) == "scantwo"))
    stop("Input should have class \"scantwo\".")
  if(missing(e2)) return(e1)
  if(!any(class(e2) == "scantwo"))
    stop("Input should have class \"scantwo\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  e1x <- "scanoneX" %in% names(e1)
  e2x <- "scanoneX" %in% names(e2)

  if(any(dim(e1$map) != dim(e2$map)) ||
     length(dim(e1$lod)) != length(dim(e2$lod)) ||
     any(dim(e1$lod) != dim(e2$lod)) || e1x != e2x)
    stop("input arguments do not conform.")

  e1$lod <- e1$lod + e2$lod
  if(e1x) e1$scanoneX <- e1$scanoneX + e2$scanoneX
  
  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 + df2
  }

  e1
}  


######################################################################
# -.scantwoperm: subtract LOD scores in two scantwo permutation results
######################################################################
"-.scantwoperm" <-
function(e1, e2)  
{
  if(!any(class(e1) == "scantwoperm"))
    stop("Input should have class \"scantwoperm\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if(missing(e2)) {
    for(i in 1:length(e1))
      e1[[i]] <- -e1[[i]]
    if(!is.null(df1))
      attr(e1, "df") <- -df1
    
    return(e1)
  }

  if(!any(class(e2) == "scantwoperm"))
    stop("Input should have class \"scantwoperm\".")

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  dim1 <- sapply(e1, dim)
  dim2 <- sapply(e2, dim)
  if(any(dim1 != dim2)) 
    stop("Need input to concern the same phenotypes and no. permutations.\n")

  for(i in 1:length(e1))
    e1[[i]] <- e1[[i]] - e2[[i]]
  
  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 - df2
  }

  e1
}


######################################################################
# +.scantwoperm: add LOD scores in two scantwo permutation results 
######################################################################
"+.scantwoperm" <-
function(e1, e2)
{
  if(!any(class(e1) == "scantwoperm"))
    stop("Input should have class \"scantwoperm\".")
  if(missing(e2)) return(e1)
  if(!any(class(e2) == "scantwoperm"))
    stop("Input should have class \"scantwoperm\".")

  if("df" %in% names(attributes(e1)))
    df1 <- attr(e1, "df")
  else df1 <- NULL

  if("df" %in% names(attributes(e2)))
    df2 <- attr(e2, "df")
  else df2 <- NULL

  dim1 <- sapply(e1, dim)
  dim2 <- sapply(e2, dim)
  if(any(dim1 != dim2)) 
    stop("Need input to concern the same phenotypes and no. permutations.\n")

  for(i in 1:length(e1))
    e1[[i]] <- e1[[i]] + e2[[i]]
  
  if(!is.null(df1) && !is.null(df2)) {
    if(length(df1) != length(df2))
      warning("Dimensions of degrees of freedom don't match; this may indicate a problem.")
    else 
      attr(e1, "df") <- df1 + df2
  }

  e1
}  


# end of arithscan.R
