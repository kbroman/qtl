#####################################################################
#
# mqmscanall.R
#
# Copyright (c) 2009-2014, Danny Arends
#
# Modified by Pjotr Prins and Karl Broman
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
# Contains: mqmscanall
#           scanall
#
#
#####################################################################

mqmscanall <- function(cross, multicore=TRUE, n.clusters=1, batchsize=10,cofactors=NULL, ...) {
    scanall(cross=cross, multicore=multicore, n.clusters=n.clusters, batchsize=batchsize,cofactors=cofactors, ..., scanfunction=mqmscan)
}

scanall <- function(cross, scanfunction=scanone, multicore=TRUE, n.clusters=1, batchsize=10, FF=0,cofactors=NULL, ..., plot=FALSE, verbose=FALSE)
{
    if(missing(cross)){
        ourstop("No cross file. Please supply a valid cross object.")
    }
    if(!(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"))
        stop("Currently only F2, BC, and selfed RIL crosses can be analyzed by MQM.")

    start <- proc.time()
    n.pheno <- nphe(cross)
    if(verbose) {
        ourline()
        cat("Starting R/QTL multitrait analysis\n")
        cat("Number of phenotypes:",n.pheno,"\n")
        cat("Batchsize:",batchsize," & n.clusters:",n.clusters,"\n")
        ourline()
    }

    result <- NULL      #BATCH result variable
    res <- NULL			#GLOBAL result variable
    all.data <- cross

    bootstraps <- 1:n.pheno
    batches <- length(bootstraps) %/% batchsize
    last.batch.num <- length(bootstraps) %% batchsize

    if(last.batch.num > 0){
        batches = batches+1
    }

    #INIT TIME VARS
    SUM <- 0
    AVG <- 0
    LEFT <- 0

    #TEST FOR SNOW CAPABILITIES
    if(multicore && n.clusters >1) {
        updateParallelRNG(n.clusters)

        if(verbose) cat("INFO: Using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
        for(x in 1:(batches)){
            start <- proc.time()
            if(verbose) cat("INFO: Starting with batch",x,"/",batches,"\n")
            if(x==batches && last.batch.num > 0){
                boots <- bootstraps[((batchsize*(x-1))+1):((batchsize*(x-1))+last.batch.num)]
            }else{
                boots <- bootstraps[((batchsize*(x-1))+1):(batchsize*(x-1)+batchsize)]
            }

            if(Sys.info()[1] == "Windows") { # Windows doesn't support mclapply, but it's faster if available
                cl <- makeCluster(n.clusters)
                on.exit(stopCluster(cl))
                result <- clusterApply(cl, boots, snowCoreALL, all.data=all.data, scanfunction=scanfunction, cofactors=cofactors,
                                       verbose=verbose, ...)
            }
            else {
                result <- mclapply(boots, snowCoreALL, all.data=all.data, scanfunction=scanfunction, cofactors=cofactors,
                                   verbose=verbose, mc.cores=n.clusters, ...)
            }

            if(plot){
                temp <- result
                class(temp) <- c(class(temp),"mqmmulti")
                mqmplot.multitrait(temp)
            }
            res <- c(res,result)
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
            if(verbose) cat("INFO: Starting with batch",x,"/",batches,"\n")
            if(x==batches && last.batch.num > 0){
                boots <- bootstraps[((batchsize*(x-1))+1):((batchsize*(x-1))+last.batch.num)]
            }else{
                boots <- bootstraps[((batchsize*(x-1))+1):(batchsize*(x-1)+batchsize)]
            }
            result <- lapply(boots, FUN=snowCoreALL,all.data=all.data,scanfunction=scanfunction,cofactors=cofactors,verbose=verbose,...)
            if(plot){
                temp <- result
                class(temp) <- c(class(temp),"mqmmulti")
                mqmplot.multitrait(temp)
            }
            res <- c(res,result)
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
    }
    if(FF){
        if(verbose) cat(rownames(res[[1]]),"\n",res[[1]][,1],"\n",res[[1]][,2],"\n",file="out.frank")
        for(i in 1:length(res)){
            if(verbose) cat("INFO: Saving trait",i,"in frankformat\n")
            qtl <- res[[i]]
            if(verbose) cat(colnames(qtl)[3],qtl[,3],"\n",file="out.frank",append = TRUE)
        }
    }
    #Return the results
    if(length(res) > 1){
        class(res) <- c(class(res),"mqmmulti")
    }else{
        class(res) <- c(class(res),"scanone")
    }
    #All done now plot the results
    end <- proc.time()
    SUM <- SUM + (end-start)[3]
    AVG <- SUM/n.pheno
    if(verbose) {
        cat("------------------------------------------------------------------\n")
        cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")
        cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
        cat("------------------------------------------------------------------\n")
    }
    res
}


# end of scanall.R
