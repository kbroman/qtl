######################################################################
#
# tryallpositions.R
#
# copyright (c) 2007-8, Karl W Broman
# last modified Jun, 2008
# first written Oct, 2007
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/qtl package
# Contains: tryallpositions, markerloglik
#
######################################################################

######################################################################
# tryallpositions
#
# Place a given marker in all possible positions on a selected set of
# chromosomes, keeping the positions of all other markers fixed, and
# evaluate the likelihood and estimate the chromosome length
######################################################################

tryallpositions <-
function(cross, marker, chr, error.prob=0.0001,
         map.function=c("haldane","kosambi","c-f","morgan"),
         m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE, verbose=TRUE)
{
  map.function <- match.arg(map.function)
  if(missing(chr)) chr <- names(cross$geno)
  
  thechr <- find.markerpos(cross, marker)[1,1]
  if(is.na(thechr))
    stop("Marker ", marker, " not found.")

  markerll <- markerloglik(cross, marker, error.prob)

  allchr <- names(subset(cross, chr=chr)$geno)

  results <- NULL
  for(i in allchr) {
    if(i == thechr) { # marker already on this chromosome
      pos <- cross$geno[[i]]$map
      if(is.matrix(pos)) {
        pos <- pos[1,]
        matrixmap <- TRUE
      }
      else matrixmap <- FALSE
      
      n.mar <- ncol(cross$geno[[i]]$data)
      if(n.mar == 1) { # this is the only marker
        pos <- 0
        length <- length.male <- 0
        llik <- 0
        if(verbose) cat(i, pos, llik/log(10), "\n")
        interval <- "---"
      } # just this marker
      else if(n.mar == 2) { # just two markers
        pos <- 0
        themar <- colnames(cross$geno[[i]]$data)
        initialloglik <- sum(sapply(themar, function(a, b, c) markerloglik(b, a, c), cross, error.prob))
        
        nm <- est.map(subset(cross, chr=i), error.prob=error.prob,
                      map.function=map.function, m=m, p=p, maxit=maxit,
                      tol=tol, sex.sp=sex.sp, verbose=FALSE,
                      omit.noninformative=FALSE)[[1]]
        llik <- attr(nm, "loglik") - initialloglik
        if(verbose) cat(i, pos, llik/log(10), "\n")

        if(is.matrix(nm)) {
          length <- diff(range(nm[1,]))
          length.male <- diff(range(nm[2,]))
        }
        else length <- diff(range(nm))
        interval <- "---"
      } # just two markers
      else { # >2 markers
        temp <- drop.markers(cross, marker)

        pos <- temp$geno[[i]]$map
        if(is.matrix(pos)) {
          pos <- pos[1,]
          matrixmap <- TRUE
        }
        else matrixmap <- FALSE
      
        pos <- c(min(pos)-10, pos, max(pos)+10)
        pos <- (pos[-1] + pos[-length(pos)])/2

        themarkers <- colnames(temp$geno[[i]]$data)
        int2 <- c("*", match(themarkers, colnames(cross$geno[[i]]$data)), "*")
        themarkers <- c("pter", themarkers, "qter")
        interval <- paste(themarkers[-length(themarkers)], themarkers[-1],
                          sep="-")
        int2 <- paste("(", int2[-length(int2)], "-", int2[-1], ")", sep="")
        interval <- paste(interval, int2)

        initialmap <- est.map(subset(temp, chr=i),
                              error.prob=error.prob,
                              map.function=map.function, m=m, p=p, maxit=maxit,
                              tol=tol, sex.sp=sex.sp, verbose=FALSE,
                              omit.noninformative=FALSE)[[1]]
        initialloglik <- attr(initialmap, "loglik") + markerll

        llik <- length <- length.male <- rep(NA, length(pos))
        for(j in seq(along=pos)) {
          temp <- movemarker(cross, marker, i, pos[j])
          nm <- est.map(subset(temp, chr=i), error.prob=error.prob,
                        map.function=map.function, m=m, p=p, maxit=maxit,
                        tol=tol, sex.sp=sex.sp, verbose=FALSE,
                        omit.noninformative=FALSE)[[1]]
          llik[j] <- attr(nm, "loglik") - initialloglik
          if(verbose) cat(i, pos[j], llik[j]/log(10), "\n")

          if(is.matrix(nm)) {
            length[j] <- diff(range(nm[1,]))
            length.male[j] <- diff(range(nm[2,]))
          }
          else length[j] <- diff(range(nm))
        }
      } # >2 markers

    } 

    else { # marker not on this chromosome

      pos <- cross$geno[[i]]$map
      if(is.matrix(pos)) {
        pos <- pos[1,]
        matrixmap <- TRUE
      }
      else matrixmap <- FALSE

      n.mar <- ncol(cross$geno[[i]]$data)

      if(n.mar > 1) {
        pos <- c(min(pos)-10, pos, max(pos)+10)
        pos <- (pos[-1] + pos[-length(pos)])/2

        themarkers <- colnames(cross$geno[[i]]$data)
        int2 <- c("*", 1:length(themarkers), "*")
        themarkers <- c("pter", themarkers, "qter")
        interval <- paste(themarkers[-length(themarkers)], themarkers[-1],
                          sep="-")
        int2 <- paste("(", int2[-length(int2)], "-", int2[-1], ")", sep="")
        interval <- paste(interval, int2)

        initialmap <- est.map(subset(cross, chr=i), error.prob=error.prob,
                              map.function=map.function, m=m, p=p, maxit=maxit,
                              tol=tol, sex.sp=sex.sp, verbose=FALSE,
                              omit.noninformative=FALSE)[[1]]
        initialloglik <- attr(initialmap, "loglik") + markerll

        llik <- length <- length.male <- rep(NA, length(pos))
        for(j in seq(along=pos)) {
          temp <- movemarker(cross, marker, i, pos[j])

          nm <- est.map(subset(temp, chr=i), error.prob=error.prob,
                        map.function=map.function, m=m, p=p, maxit=maxit,
                        tol=tol, sex.sp=sex.sp, verbose=FALSE,
                        omit.noninformative=FALSE)[[1]]
          llik[j] <- attr(nm, "loglik") - initialloglik
          if(verbose) cat(i, pos[j], llik[j]/log(10), "\n")

          if(is.matrix(nm)) {
            length[j] <- diff(range(nm[1,]))
            length.male[j] <- diff(range(nm[2,]))
          }
          else length[j] <- diff(range(nm))
        }
      } # >1 marker on chromosome
      else {
        initialloglik <- markerloglik(cross, colnames(cross$geno[[i]]$data), error.prob) + markerll

        pos <- pos+10
        interval <- "---"

        temp <- movemarker(cross, marker, i, pos)

        nm <- est.map(subset(temp, chr=i), error.prob=error.prob,
                      map.function=map.function, m=m, p=p, maxit=maxit,
                      tol=tol, sex.sp=sex.sp, verbose=FALSE,
                      omit.noninformative=FALSE)[[1]]
        llik <- attr(nm, "loglik") - initialloglik
        if(verbose) cat(i, pos, llik/log(10), "\n")

        if(is.matrix(nm)) {
          length <- diff(range(nm[1,]))
          length.male <- diff(range(nm[2,]))
        }
        else length <- diff(range(nm))
        
      } # one marker on chromosome
    } # marker not on this chr

    if(matrixmap && sex.sp) 
      tempres <- data.frame(chr=rep(i, length(pos)),
                            pos=pos,
                            lod=llik/log(10),
                            length.female=length,
                            length.male=length.male,
                            interval=interval,
                            stringsAsFactors=FALSE)
    else
      tempres <- data.frame(chr=rep(i, length(pos)),
                            pos=pos,
                            lod=llik/log(10),
                            length=length,
                            interval=interval,
                            stringsAsFactors=FALSE)


    results <- rbind(results, tempres)
  } # loop over chromosomes

  rownames(results) <- 1:nrow(results)
  results
}

######################################################################
#
# markerloglik: Calculate log likelihood for a given marker
#
######################################################################

markerloglik <- 
function(cross, marker, error.prob=0.0001) 
{
  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  type <- class(cross)[1]

  # don't let error.prob be exactly zero (or >1)
  if(error.prob < 1e-50) error.prob <- 1e-50
  if(error.prob > 1) {
    error.prob <- 1-1e-50
    warning("error.prob shouldn't be > 1!")
  }

  n.ind <- nind(cross)

  thechr <- find.markerpos(cross, marker)[1,1]
  if(is.na(thechr))
    stop("Marker ", marker, " not found.")

  chrtype <- class(cross$geno[[thechr]])

  g <- pull.geno(cross, chr=thechr)
  m <- match(marker, colnames(g))
  if(is.na(m))
    stop("Marker ", marker, " not found.")
  g <- g[,m]
  g[is.na(g)] <- 0

  # which type of cross is this?
  if(type == "f2") {
    if(chrtype == "A") # autosomal
      cfunc <- "marker_loglik_f2"
    else                  # X chromsome 
      cfunc <- "marker_loglik_bc"
  }
  else if(type == "bc" || type=="riself" || type=="risib" || type=="dh") {
    cfunc <- "marker_loglik_bc"
  }
  else if(type == "4way") {
    cfunc <- "marker_loglik_4way"
  }
  else 
    stop("markerloglik not available for cross type ", type, ".")

  # call the C function
  z <- .C(cfunc,
          as.integer(n.ind),       # number of individuals
          as.integer(g),           # genotype data
          as.double(error.prob),     
          loglik=as.double(0),     # log likelihood
          PACKAGE="qtl")

  z$loglik
}

# end of tryallpositions.R
