######################################################################
#
# est.map.R
#
# copyright (c) 2001-2014, Karl W Broman
# last modified Aug, 2014
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
# Contains: est.map
#
######################################################################

######################################################################
#
# est.map: re-estimate the genetic map for an experimental cross
#
######################################################################

est.map <-
    function(cross, chr, error.prob=0.0001, map.function=c("haldane","kosambi","c-f","morgan"),
             m=0, p=0, maxit=10000, tol=1e-6, sex.sp=TRUE, verbose=FALSE,
             omit.noninformative=TRUE, offset, n.cluster=1)
{
    if(!("cross" %in% class(cross)))
        stop("Input should have class \"cross\".")

    if(!missing(chr))
        cross <- subset(cross, chr=chr)

    type <- class(cross)[1]

    if(!missing(offset)) {
        if(length(offset)==1) offset <- rep(offset, nchr(cross))
        else if(length(offset) != nchr(cross))
            stop("offset must have length 1 or n.chr (", nchr(cross), ")")
    }

    if(m < 0 || p < 0 || p > 1)
        stop("Must have m >=0 and 0 <= p <= 1")

    if(m > 0 && p < 1 && type != "bc" && type != "f2") {
        warning("m and p currently used only for backcrosses and intercrosses.")
        m <- p <- 0
    }
    if(m > 0 && p < 1 && !missing(map.function))
        warning("Map function not used with interference model.")
    if(m > 0 && p < 1) interf.model <- TRUE
    else interf.model <- FALSE

    # map function
    map.function <- match.arg(map.function)
    if(map.function=="kosambi") {
        mf <- mf.k; imf <- imf.k
    }
    else if(map.function=="c-f") {
        mf <- mf.cf; imf <- imf.cf
    }
    else if(map.function=="morgan") {
        mf <- mf.m; imf <- imf.m
    }
    else {
        mf <- mf.h; imf <- imf.h
    }

    # don't let error.prob be exactly zero (or >1)
    if(error.prob < 1e-50) error.prob <- 1e-50
    if(error.prob > 1) {
        error.prob <- 1-1e-50
        warning("error.prob shouldn't be > 1!")
    }

    n.ind <- nind(cross)
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)

    newmap <- vector("list",n.chr)
    names(newmap) <- names(cross$geno)
    chrtype <- sapply(cross$geno, class)

    if(n.cluster > 1 && nchr(cross) > 1) {
        cat(" -Running est.map via a cluster of", n.cluster, "nodes.\n")
        cl <- makeCluster(n.cluster)
        clusterStopped <- FALSE
        on.exit(if(!clusterStopped) stopCluster(cl))
        clusterEvalQ(cl, library(qtl, quietly=TRUE))

        chr <- names(cross$geno)

        # temporary definition of est.map
        temp.est.map <- function(chr, cross, error.prob, map.function, m, p, maxit, tol,
                                 sex.sp, omit.noninformative)
            est.map(cross=cross, chr=chr, error.prob=error.prob, map.function=map.function,
                    m=m, p=p, maxit=maxit, tol=tol, sex.sp=sex.sp, omit.noninformative=omit.noninformative,
                    verbose=FALSE)#, n.cluster=1)

        newmap <- clusterApplyLB(cl, chr, temp.est.map, cross, error.prob, map.function, m, p,
                                 maxit, tol, sex.sp, omit.noninformative)
        for(i in seq(along=newmap)) {
            newmap[[i]] <- newmap[[i]][[1]]
            class(newmap[[i]]) <- class(cross$geno[[i]])
        }

        names(newmap) <- chr

        if(!missing(offset)) {  # shift map start positions
            for(i in seq(along=newmap))
                if(is.matrix(newmap[[i]])) {
                    for(j in 1:2)
                        newmap[[i]][j,] <- newmap[[i]][j,] - newmap[[i]][j,1] + offset[i]
                } else {
                    newmap[[i]] <- newmap[[i]] - newmap[[i]][1] + offset[i]
                }
        }

        class(newmap) <- "map"
        return(newmap)
    }

    # calculate genotype probabilities one chromosome at a time
    for(i in 1:n.chr) {

        if(n.mar[i] < 2) {
            newmap[[i]] <- cross$geno[[i]]$map
            next
        }

        # which type of cross is this?
        if(type == "f2") {
            one.map <- TRUE
            if(chrtype[i] != "X") # autosomal
                cfunc <- "est_map_f2"
            else                              # X chromsome
                cfunc <- "est_map_bc"
        }
        else if(type == "bc" || type=="riself" || type=="risib" || type=="dh" || type=="haploid") {
            one.map <- TRUE
            cfunc <- "est_map_bc"
        }
        else if(type == "4way") {
            one.map <- FALSE
            cfunc <- "est_map_4way"
        }
        else if(type=="ri8sib" || type=="ri4sib" || type=="ri8self" || type=="ri4self" || type=="bgmagic16") {
            cfunc <- paste("est_map_", type, sep="")
            one.map <- TRUE
            if(chrtype[i] == "X")
                warning("est.map not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else if(type == "bcsft") {
            one.map <- TRUE
            interf.model <- FALSE
            cfunc <- "est_map_bcsft"
            cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
            if(chrtype[i] == "X") { # X chromsome
                cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
                cross.scheme[2] <- 0
            }
            ## Tolerance: need two values.
            if(length(tol) == 1) {
                tol[2] <- 1e-6
            }
        }
        else
            stop("est.map not available for cross type ", type, ".")

        # genotype data
        gen <- cross$geno[[i]]$data
        gen[is.na(gen)] <- 0

        # remove individuals that have less than two typed markers
        if(omit.noninformative) {
            o <- apply(gen,1,function(a) sum(a!=0)>1)
            gen <- gen[o,,drop=FALSE]
        }

        # recombination fractions
        if(one.map) {
            # recombination fractions
            rf <- mf(diff(cross$geno[[i]]$map))
            if(type=="risib" || type=="riself")
                rf <- adjust.rf.ri(rf,substr(type,3,nchar(type)),chrtype[i])
            rf[rf < 1e-14] <- 1e-14
        }
        else {
            orig <- cross$geno[[i]]$map
            # randomize the maps a bit [we no longer do this]
            #      cross$geno[[i]]$map <- cross$geno[[i]]$map +
            #        runif(length(cross$geno[[i]]$map), -0.2, 0.2)

            rf <- mf(diff(cross$geno[[i]]$map[1,]))
            rf[rf < 1e-14] <- 1e-14
            rf2 <- mf(diff(cross$geno[[i]]$map[2,]))
            rf2[rf2 < 1e-14] <- 1e-14
            if(!sex.sp && chrtype[i]=="X")
                temp.sex.sp <- TRUE
            else temp.sex.sp <- sex.sp
        }

        if(interf.model)
            d <- diff(cross$geno[[i]]$map)

        if(verbose) cat(paste("Chr ", names(cross$geno)[i], ":\n",sep=""))

        # call the C function
        if(one.map && !interf.model) {
            ## Hide cross scheme in genoprob to pass to routine. BY
            temp <- 0
            if(type == "bcsft")
                temp[1] <- cross.scheme[1] * 1000 + cross.scheme[2]

            z <- .C(cfunc,
                    as.integer(nrow(gen)),         # number of individuals
                    as.integer(n.mar[i]),      # number of markers
                    as.integer(gen),           # genotype data
                    rf=as.double(rf),          # recombination fractions
                    as.double(error.prob),
                    loglik=as.double(temp),       # log likelihood
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(verbose),
                    PACKAGE="qtl")

            z$rf[z$rf < 1e-14] <- 1e-14
            if(type=="riself" || type=="risib")
                z$rf <- adjust.rf.ri(z$rf, substr(type, 3, nchar(type)),
                                     chrtype[i], expand=FALSE)
            newmap[[i]] <- cumsum(c(min(cross$geno[[i]]$map),imf(z$rf)))
            names(newmap[[i]]) <- names(cross$geno[[i]]$map)
            attr(newmap[[i]],"loglik") <- z$loglik
        }
        else if(interf.model) { # Chi-square / Stahl model
            if(type=="bc" || (type=="f2" && chrtype[i]=="X")) {
                z <- .C("R_est_map_bci",
                        as.integer(nrow(gen)),         # number of individuals
                        as.integer(n.mar[i]),      # number of markers
                        as.integer(gen),           # genotype data
                        d=as.double(d),          # cM distances
                        as.integer(m),
                        as.double(p),
                        as.double(error.prob),
                        loglik=as.double(0),       # log likelihood
                        as.integer(maxit),
                        as.double(tol),
                        as.integer(verbose),
                        PACKAGE="qtl")
            } else {
                z <- .C("R_est_map_f2i",
                        as.integer(nrow(gen)),         # number of individuals
                        as.integer(n.mar[i]),      # number of markers
                        as.integer(gen),           # genotype data
                        d=as.double(d),          # cM distances
                        as.integer(m),
                        as.double(p),
                        as.double(error.prob),
                        loglik=as.double(0),       # log likelihood
                        as.integer(maxit),
                        as.double(tol),
                        as.integer(verbose),
                        PACKAGE="qtl")
            }
            z$d[z$d < 1e-14] <- 1e-14
            newmap[[i]] <- cumsum(c(min(cross$geno[[i]]$map),z$d))
            names(newmap[[i]]) <- names(cross$geno[[i]]$map)
            attr(newmap[[i]], "loglik") <- z$loglik
            attr(newmap[[i]], "m") <- m
            attr(newmap[[i]], "p") <- p
        }
        else {
            z <- .C(cfunc,
                    as.integer(nrow(gen)),         # number of individuals
                    as.integer(n.mar[i]),      # number of markers
                    as.integer(gen),           # genotype data
                    rf=as.double(rf),          # recombination fractions
                    rf2=as.double(rf2),        # recombination fractions
                    as.double(error.prob),
                    loglik=as.double(0),       # log likelihood
                    as.integer(maxit),
                    as.double(tol),
                    as.integer(temp.sex.sp),
                    as.integer(verbose),
                    PACKAGE="qtl")

            z$rf[z$rf<1e-14] <- 1e-14
            z$rf2[z$rf2<1e-14] <- 1e-14

            if(!temp.sex.sp) z$rf2 <- z$rf

            newmap[[i]] <- rbind(cumsum(c(min(orig[1,]),imf(z$rf))),
                                 cumsum(c(min(orig[2,]),imf(z$rf2))))
            dimnames(newmap[[i]]) <- dimnames(cross$geno[[i]]$map)
            attr(newmap[[i]],"loglik") <- z$loglik
        }

        class(newmap[[i]]) <- chrtype[i]
    } # end loop over chromosomes


    if(!missing(offset)) {  # shift map start positions
        for(i in seq(along=newmap))
            if(is.matrix(newmap[[i]])) {
                for(j in 1:2)
                    newmap[[i]][j,] <- newmap[[i]][j,] - newmap[[i]][j,1] + offset[i]
            } else {
                newmap[[i]] <- newmap[[i]] - newmap[[i]][1] + offset[i]
            }
    }


    class(newmap) <- "map"
    newmap
}

# end of est.map.R
