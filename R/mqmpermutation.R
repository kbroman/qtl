#####################################################################
#
# mqmpermutation.R
#
# Copyright (c) 2009-2014, Danny Arends
#
# Modified by Pjotr Prins and Karl Broman
#
#
# first written Februari 2009
# last modified Jan 2014
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
# Contains: mqmscanfdr
#           mqmpermutation
#           mqmprocesspermutation
#
#
#####################################################################

mqmscanfdr <- function(cross, scanfunction=mqmscanall, thresholds=c(1,2,3,4,5,7,10,15,20), n.perm = 10, verbose=FALSE, ...){
    if(verbose){cat("Calculation of FDR estimate of threshold in multitrait analysis.\n")}
    results <- NULL
    above.in.real.res <- NULL
    res <- scanfunction(cross,...)
    for(threshold in thresholds){
        above.in.real <- 0
        for(x in 1:nphe(cross)){
            above.in.real = above.in.real + sum(res[[x]][,3] > threshold)
        }
        above.in.real.res <- c(above.in.real.res,above.in.real)
    }
    perm <- cross
    if(verbose){cat("QTL's above threshold:",above.in.real,"\n")}
    above.in.perm.res <- rep(0,length(thresholds))
    for(x in 1:n.perm){
        if(verbose){cat("Starting permutation",x,"\n")}
        perm.res <- NULL

        neworder <- sample(nind(cross))
        for(chr in 1:nchr(cross)){
            perm$geno[[1]]$data <- perm$geno[[1]]$data[neworder,]
        }
        res <- scanfunction(perm,...)
        for(threshold in thresholds){
            above.in.perm <- 0
            for(y in 1:nphe(cross)){
                above.in.perm = above.in.perm + sum(res[[y]][,3] > threshold)
            }
            perm.res <- c(perm.res,above.in.perm)
            #if(verbose){cat("Permutation",x,"QTL's above threshold:",above.in.perm,"\n")}
        }
        above.in.perm.res <- above.in.perm.res+perm.res
    }
    above.in.perm.res <- above.in.perm.res/n.perm
    results <- cbind(above.in.real.res,above.in.perm.res,above.in.perm.res/above.in.real.res)
    rownames(results) <- thresholds
    results
}
######################################################################
#
# mqmpermutation: Shuffles phenotype or does parametric bootstrapping of mqmscan
#
######################################################################

mqmpermutation <- function(cross,scanfunction=scanone,pheno.col=1,multicore=TRUE,n.perm=10,batchsize=10,file="MQM_output.txt",n.cluster=1,method=c("permutation","simulation"),cofactors=NULL,plot=FALSE,verbose=FALSE,...)
{
    bootmethod <- 0

    supported <- c("permutation","simulation")
    bootmethod <- pmatch(method, supported)[1]-1
    if(missing(cross))
        stop("No cross file. Please supply a valid cross object.")

    if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
        #Echo back the cross type
        if(verbose) {
            cat("------------------------------------------------------------------\n")
            cat("Starting permutation analysis\n")
            cat("Number of permutations:",n.perm,"\n")
            cat("Batchsize:",batchsize," & n.cluster:",n.cluster,"\n")
            cat("------------------------------------------------------------------\n")
            cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
        }
        b <- proc.time()
        if(!bootmethod){
            if(verbose) cat("INFO: Shuffleling traits between individuals.\n")
        }else{
            if(verbose) cat("INFO: Parametric permutation\nINFO: Calculating new traits for each individual.\n")
        }

        #Set the Phenotype under interest as the first
        cross$pheno[[1]] <- cross$pheno[[pheno.col]]
        names(cross$pheno)[1] <- names(cross$pheno)[pheno.col]

        if(n.cluster > batchsize){
            stop("Please have more items in a batch then clusters assigned per batch")
        }

        #Scan the original
        #cross <- fill.geno(cross) # <- this should be done outside of this function
        res0 <- lapply(1, FUN=snowCoreALL,all.data=cross,scanfunction=scanfunction,verbose=verbose,cofactors=cofactors,...)

        #Setup bootstraps by generating a list of random numbers to set as seed for each bootstrap
        bootstraps <- runif(n.perm)
        batches <- length(bootstraps) %/% batchsize
        last.batch.num <- length(bootstraps) %% batchsize
        results <- NULL
        if(last.batch.num > 0){
            batches = batches+1
        }
        SUM <- 0
        AVG <- 0
        LEFT <- 0
        if(multicore && n.cluster >1) {
            updateParallelRNG(n.cluster)
            if(verbose) cat("INFO: Using ",n.cluster," Cores/CPU's/PC's for calculation.\n")
            for(x in 1:(batches)){
                start <- proc.time()
                if(verbose) {
                    ourline()
                    cat("INFO: Starting with batch",x,"/",batches,"\n")
                    ourline()
                }
                if(x==batches && last.batch.num > 0){
                    boots <- bootstraps[((batchsize*(x-1))+1):((batchsize*(x-1))+last.batch.num)]
                }else{
                    boots <- bootstraps[((batchsize*(x-1))+1):(batchsize*(x-1)+batchsize)]
                }
                if(Sys.info()[1] == "Windows") { # Windows doesn't support mclapply, but it's faster if available
                    cl <- makeCluster(n.cluster)
                    on.exit(stopCluster(cl))
                    res <- clusterApply(cl, boots, snowCoreBOOT, all.data=cross, scanfunction=scanfunction, bootmethod=bootmethod,
                                        cofactors=cofactors, verbose=verbose, ...)
                }
                else {
                    res <- mclapply(boots, snowCoreBOOT, all.data=cross, scanfunction=scanfunction, bootmethod=bootmethod,
                                    cofactors=cofactors, verbose=verbose, mc.cores=n.cluster, ...)
                }
                results <- c(results,res)
                if(plot){
                    temp <- c(res0,results)
                    class(temp) <- c(class(temp),"mqmmulti")
                    mqmplot.permutations(temp)
                }
                end <- proc.time()
                SUM <- SUM + (end-start)[3]
                AVG <- SUM/x
                LEFT <- AVG*(batches-x)
                if(verbose) {
                    cat("INFO: Done with batch",x,"/",batches,"\n")
                    cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
                    cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
                    cat("INFO: Average time per batch:",round((AVG), digits=3)," per trait:",round((AVG %/% batchsize), digits=3),"seconds\n")
                    cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
                    ourline()
                }
            }
        }else{
            if(verbose) cat("INFO: Going into singlemode.\n")
            for(x in 1:(batches)){
                start <- proc.time()
                if(verbose) {
                    ourline()
                    cat("INFO: Starting with batch",x,"/",batches,"\n")
                    ourline()
                }
                if(x==batches && last.batch.num > 0){
                    boots <- bootstraps[((batchsize*(x-1))+1):((batchsize*(x-1))+last.batch.num)]
                }else{
                    boots <- bootstraps[((batchsize*(x-1))+1):(batchsize*(x-1)+batchsize)]
                }
                res <- lapply(boots, FUN=snowCoreBOOT,all.data=cross,scanfunction=scanfunction,bootmethod=bootmethod,cofactors=cofactors,verbose=verbose,...)
                results <- c(results,res)
                if(plot){
                    temp <- c(res0,results)
                    class(temp) <- c(class(temp),"mqmmulti")
                    mqmplot.permutations(temp)
                }
                end <- proc.time()
                SUM <- SUM + (end-start)[3]
                AVG <- SUM/x
                LEFT <- AVG*(batches-x)
                if(verbose) {
                    cat("INFO: Done with batch",x,"/",batches,"\n")
                    cat("INFO: Calculation of batch",x,"took:",round((end-start)[3], digits=3),"seconds\n")
                    cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
                    cat("INFO: Average time per batch:",round((AVG), digits=3),",per run:",round((AVG %/% batchsize), digits=3),"seconds\n")
                    cat("INFO: Estimated time left:",LEFT%/%3600,":",(LEFT%%3600)%/%60,":",round(LEFT%%60,digits=0),"(Hour:Min:Sec)\n")
                    ourline()
                }
            }
        }
        res <- c(res0,results)
        #Set the class of the result to mqmmulti (so we can use our plotting routines)
        class(res) <- c(class(res),"mqmmulti")
        e <- proc.time()
        SUM <- (e-b)[3]
        AVG <- SUM/(n.perm+1)
        if(verbose) {
            cat("INFO: Done with MQM permutation analysis\n")
            cat("------------------------------------------------------------------\n")
            cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
            cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
            cat("------------------------------------------------------------------\n")
        }
        res
    }else{
        stop("Currently only F2, BC, and selfed RIL crosses can be analyzed by MQM.")
    }
}

mqmprocesspermutation <- function(mqmpermutationresult = NULL){
    if(!is.null(mqmpermutationresult) && class(mqmpermutationresult)[2] == "mqmmulti"){
        result <- NULL
        result <- sapply(mqmpermutationresult[-1], function(a) max(a[,3], na.rm=TRUE))
        result <- as.matrix(result)
        colnames(result) <- colnames(mqmpermutationresult[[2]])[3]
        rownames(result) <- 1:(length(mqmpermutationresult)-1)
        class(result) <- c("scanoneperm",class(result))
        result
    }else{
        stop("Please supply a valid resultobject (mqmmulti).")
    }
}

# end of mqmpermutation.R
