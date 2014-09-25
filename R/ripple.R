######################################################################
#
# ripple.R
#
# copyright (c) 2001-2014, Karl W Broman
# last modified Aug, 2014
# first written Oct, 2001
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
# Contains: ripple, summary.ripple, print.summary.ripple
#           ripple.perm1, ripple.perm2, ripple.perm.sub
#
######################################################################

######################################################################
#
# ripple: Check marker orders for a given chromosome, comparing all
#         possible permutations of a sliding window of markers
#
######################################################################

ripple <-
    function(cross, chr, window=4, method=c("countxo","likelihood"),
             error.prob=0.0001, map.function=c("haldane","kosambi","c-f","morgan"),
             maxit=4000, tol=1e-6, sex.sp=TRUE, verbose=TRUE, n.cluster=1)
{
    if(!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    # pull out relevant chromosome
    if(missing(chr)) {
        chr <- names(cross$geno)[1]
        warning("chr argument not provided; assuming you want chr ", chr)
    }
    else {
        if(length(chr) > 1)
            stop("ripple only works for one chromosome at a time.")
        if(!testchr(chr, names(cross$geno)))
            stop("Chr ", chr, " not found.")
    }
    cross <- subset(cross,chr=chr)
    chr.name <- names(cross$geno)[1]

    if(nmar(cross)[1] < 3) {
        warning("Less than three markers.")
        return(NULL)
    }

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }

    # make sure window is an integer >= 2
    if(window < 2) {
        warning("The window argument must be > 1; using window=2.")
        window <- 2
    }
    window <- round(window)

    method <- match.arg(method)
    map.function <- match.arg(map.function)

    # get marker orders to test
    n.mar <- totmar(cross)
    if(n.mar <= window) # look at all possible orders
        orders <- ripple.perm2(n.mar)
    else {
        temp <- ripple.perm1(window)
        n <- nrow(temp)
        orders <- cbind(temp,matrix(rep((window+1):n.mar,n),
                                    byrow=TRUE,ncol=n.mar-window))
        for(i in 2:(n.mar-window+1)) {
            left <- matrix(rep(1:(i-1),n),byrow=TRUE,ncol=i-1)
            if(i < n.mar-window+1)
                right <- matrix(rep((i+window):n.mar,n),byrow=TRUE,ncol=n.mar-window-i+1)
            else
                right <- NULL
            orders <- rbind(orders,cbind(left,temp+i-1,right))
        }
        # keep only distinct orders
        orders <- as.numeric(unlist(strsplit(unique(apply(orders,1,paste,collapse=":")),":")))
        orders <- matrix(orders,ncol=n.mar,byrow=TRUE)
    }
    n.orders <- nrow(orders)


    # how often to print information about current order being considered
    if(n.orders > 49) print.by <- 10
    else if(n.orders > 14) print.by <- 5
    else print.by <- 2

    if(method=="likelihood") {
        # calculate log likelihoods (and est'd chr length) for each marker order
        loglik <- 1:n.orders
        chrlen <- 1:n.orders

        # create temporary cross
        m <- seq(0,by=5,length=n.mar)
        temcross <- cross
        if(is.matrix(cross$geno[[1]]$map))
            temcross$geno[[1]]$map <- rbind(m,m)
        else temcross$geno[[1]]$map <- m

        if(verbose) cat("  ", n.orders,"total orders\n")
        if(n.cluster > 1) {
            # parallelize
            if(n.orders <= n.cluster) n.cluster <- n.orders
            cl <- makeCluster(n.cluster)
            clusterStopped <- FALSE
            on.exit(if(!clusterStopped) stopCluster(cl))
            if(verbose) cat("   Running in", n.cluster, "clusters\n")
            clusterEvalQ(cl, library(qtl, quietly=TRUE))

            whclust <- sort(rep(1:n.cluster, ceiling(n.orders/n.cluster))[1:n.orders])
            order.split <- vector("list", n.cluster)
            for(i in 1:n.cluster)
                order.split[[i]] <- orders[whclust==i,,drop=FALSE]
            result <- parLapply(cl, order.split, rippleSnowLik, cross=temcross,
                                error.prob=error.prob, map.function=map.function, maxit=maxit, tol=tol,
                                sex.sp=sex.sp)
            loglik <- unlist(lapply(result, function(a) a$loglik))
            chrlen <- unlist(lapply(result, function(a) a$chrlen))
        }
        else {
            for(i in 1:n.orders) {
                if(verbose && (i %/% print.by)*print.by == i) cat("    --Order", i, "\n")
                temcross$geno[[1]]$data <- cross$geno[[1]]$data[,orders[i,]]
                newmap <- est.map(temcross, error.prob=error.prob, map.function=map.function,
                                  m=0, p=0, maxit=maxit, tol=tol, sex.sp=sex.sp, verbose=FALSE)
                loglik[i] <- attr(newmap[[1]],"loglik")
                chrlen[i] <- diff(range(newmap[[1]]))
            }
        }

        # re-scale log likelihoods and convert to lods
        loglik <- (loglik - loglik[1])/log(10)

        # sort orders by lod
        o <- order(loglik[-1], decreasing=TRUE)+1

        # create output
        orders <- cbind(orders,LOD=loglik,chrlen)[c(1,o),]
    }
    else { # count obligate crossovers for each order
        # which type of cross is this?
        type <- class(cross)[1]
        is.bcs <- type == "bcsft"
        if(is.bcs)
            is.bcs <- (attr(cross, "scheme")[2] == 0)

        if(type == "f2" || (type == "bcsft" && !is.bcs)) {
            if(class(cross$geno[[1]]) == "A") # autosomal
                func <- "R_ripple_f2"
            else func <- "R_ripple_bc"        # X chromsome
        }
        else if(type %in% c("bc", "riself", "risib", "dh", "haploid", "bcsft")) func <- "R_ripple_bc"
        else if(type == "4way") func <- "R_ripple_4way"
        else if(type=="ri4self" || type=="ri8self" || type=="ri4sib" || type=="ri8sib" || type=="bgmagic16")
            func <- "R_ripple_ril48"
        else
            stop("ripple not available for cross ", type)

        # data to be input
        genodat <- cross$geno[[1]]$data
        genodat[is.na(genodat)] <- 0
        n.ind <- nind(cross)

        if(verbose) cat("  ", n.orders,"total orders\n")
        if(n.cluster > 1) {
            # parallelize
            if(n.orders <= n.cluster) n.cluster <- n.orders

            cl <- makeCluster(n.cluster)
            clusterStopped <- FALSE
            on.exit(if(!clusterStopped) stopCluster(cl))
            if(verbose) cat("   Running in", n.cluster, "clusters\n")
            clusterEvalQ(cl, library(qtl, quietly=TRUE))

            whclust <- sort(rep(1:n.cluster, ceiling(n.orders/n.cluster))[1:n.orders])
            order.split <- vector("list", n.cluster)
            for(i in 1:n.cluster)
                order.split[[i]] <- orders[whclust==i,,drop=FALSE]
            oblxo <- unlist(parLapply(cl, order.split, rippleSnowCountxo, genodat=genodat, func=func))
            stopCluster(cl)
            clusterStopped <- TRUE
        }
        else {
            z <- .C(func,
                    as.integer(n.ind),
                    as.integer(n.mar),
                    as.integer(genodat),
                    as.integer(n.orders),
                    as.integer(orders-1),
                    oblxo=as.integer(rep(0,n.orders)),
                    as.integer(print.by),
                    PACKAGE="qtl")
            oblxo <- z$oblxo
        }
        # sort orders by lod
        o <- order(oblxo[-1])+1

        # create output
        orders <- cbind(orders,obligXO=oblxo)[c(1,o),]
    }

    rownames(orders) <- c("Initial", paste(1:(nrow(orders)-1)))
    class(orders) <- c("ripple","matrix")
    attr(orders,"chr") <- chr.name
    attr(orders,"window") <- window
    attr(orders,"error.prob") <- error.prob
    attr(orders,"method") <- method

    # make sure, for each order considered, that the proximal marker
    # (in the original order) is to the left of the distal marker
    # (in the original order)
    orders[,1:n.mar] <- t(apply(orders[,1:n.mar,drop=FALSE],1,
                                function(a) {
                                    n <- length(a)
                                    if((1:n)[a==1] > (1:n)[a==n]) return(rev(a))
                                    else return(a) }))

    orders
}

######################################################################
# function for method="likelihood", for parallel processing (formerly with snow pkg)
######################################################################
rippleSnowLik <-
    function(orders, cross, error.prob, map.function, maxit, tol, sex.sp)
{
    n.orders <- nrow(orders)
    temcross <- cross
    loglik <- chrlen <- rep(NA, n.orders)
    for(i in 1:n.orders) {

        temcross$geno[[1]]$data <- cross$geno[[1]]$data[,orders[i,]]
        newmap <- est.map(temcross, error.prob=error.prob, map.function=map.function,
                          m=0, p=0, maxit=maxit, tol=tol, sex.sp=sex.sp, verbose=FALSE)
        loglik[i] <- attr(newmap[[1]],"loglik")
        chrlen[i] <- diff(range(newmap[[1]]))
    }
    list(loglik=loglik, chrlen=chrlen)
}



######################################################################
# function for method="countxo", for parallel processing (formerly with snow pkg)
######################################################################
rippleSnowCountxo <-
    function(orders, genodat, func)
{
    func <- func # this avoids a Note from R CMD check
    .C(func,
       as.integer(nrow(genodat)),
       as.integer(ncol(genodat)),
       as.integer(genodat),
       as.integer(nrow(orders)),
       as.integer(orders-1),
       oblxo=as.integer(rep(0, nrow(orders))),
       as.integer(0),
       PACKAGE="qtl")$oblxo
}


######################################################################
#
# summary.ripple: print top results from ripple().  We do this so
#                 that we can return *all* results but allow easy
#                 view of only the important ones
#
######################################################################

summary.ripple <-
    function(object, lod.cutoff = -1, ...)
{
    if(!any(class(object) == "ripple"))
        stop("Input should have class \"ripple\".")

    n <- ncol(object)

    if("obligXO" %in% colnames(object)) # counts of crossovers
        o <- (object[-1,n] <= (object[1,n] - lod.cutoff*2))
    else o <- (object[-1,n-1] >= lod.cutoff) # likelihood analysis

    if(!any(o)) object <- object[1:2,,drop=FALSE]
    else  # make sure first row is included
        object <- object[c(TRUE,o),,drop=FALSE]

    rownames(object) <- c("Initial ", paste(1:(nrow(object)-1)))
    class(object) <- c("summary.ripple","matrix")
    object
}

######################################################################
#
# print.summary.ripple
#
######################################################################

print.summary.ripple <-
    function(x, ...)
{
    n <- ncol(x)
    x <- round(x,1)
    max.row <- 6

    if(!("obligXO" %in% colnames(x)))
        colnames(x)[n-1] <- "    LOD"

    class(x) <- "matrix"


    if(nrow(x) > max.row) {
        print(x[1:max.row,])
        cat("... [", nrow(x)-max.row, " additional rows] ...\n")
    }
    else print(x)
}

######################################################################
#
# ripple.perm1: Utility function for ripple().  Returns all possible
#               permutations of {1, 2, ..., n}
#
######################################################################

ripple.perm1 <-
    function(n)
{
    if(n == 1) return(rbind(1))
    o <- rbind(c(n-1,n),c(n,n-1))
    if(n > 2)
        for(i in (n-2):1)
            o <- ripple.perm.sub(i,o)
    dimnames(o) <- NULL
    o
}

######################################################################
#
# ripple.perm2: Utility function for ripple().  Returns all possible
#               permutations of {1, 2, ..., n}, up to orientation of
#               the entire group
#
######################################################################

ripple.perm2 <-
    function(n)
{
    if(n < 3) return(rbind(1:n))
    o <- rbind(c(n-2,n-1,n),c(n-1,n-2,n),c(n-1,n,n-2))
    if(n > 3)
        for(i in (n-3):1)
            o <- ripple.perm.sub(i,o)
    dimnames(o) <- NULL
    o
}

######################################################################
#
# ripple.perm.sub: Subroutine used for ripple().  I'm too tired to
#                  explain.
#
######################################################################

ripple.perm.sub <-
    function(x,mat)
{
    res <- cbind(x,mat)
    if(ncol(mat) > 1) {
        for(i in 1:ncol(mat))
            res <- rbind(res,cbind(mat[,1:i],x,mat[,-(1:i)]))
    }
    res
}

# end of ripple.R
