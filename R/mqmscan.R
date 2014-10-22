#####################################################################
#
# mqmscan.R
#
# Copyright (c) 2009-2013, Danny Arends
#
# Modified by Pjotr Prins and Karl Broman
#
#
# first written Februari 2009
# last modified Sep 2013
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
# Contains: mqmscan
#
#
#####################################################################

######################################################################
#
# mqmscan: main scanning function to the mqmpackage
#
######################################################################

mqmscan <- function(cross,cofactors=NULL,pheno.col=1,model=c("additive","dominance"),forceML=FALSE,
                    cofactor.significance=0.02,em.iter=1000,window.size=25.0,step.size=5.0,logtransform = FALSE,
                    estimate.map = FALSE,plot=FALSE,verbose=FALSE, outputmarkers=TRUE, multicore=TRUE, batchsize=10, n.clusters=1,
                    test.normality=FALSE,off.end=0)
{

    start <- proc.time()
    model <- match.arg(model)

    #Because iirc we cannot pass booleans from R to C
    if(forceML){
        forceML <- 1           #1 -> Maximum Likelyhood
    }else{
        forceML <- 0           #0 -> Restricted Maximum Likelyhood
    }
    dominance <- 0            #We code :0 -> Additive model (no_dominance)
    if(model=="dominance"){
        dominance <- 1          #and 1 -> Dominance model
    }
    if(estimate.map){
        estimate.map <- 1
    }else{
        estimate.map <- 0
    }
    n.run <- 0
    if(is.null(cross)){
        stop("No cross file. Please supply a valid cross object.")
    }
    if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
        if(class(cross)[1] == "f2"){
            ctype = 1
        }
        if(class(cross)[1] == "bc" || class(cross)[1]=="dh" || class(cross)[1]=="haploid"){
            ctype = 2
        }
        if(class(cross)[1] == "riself") {
            ctype = 3

            # check genotypes
            g <- as.numeric(pull.geno(cross))
            g <- sort(unique(g[!is.na(g)]))
            if(max(g)==2) { # convert genotypes from 1/2 to 1/3
                for(i in seq(along=cross$geno))
                    cross$geno[[i]]$data[!is.na(cross$geno[[i]]$data) & cross$geno[[i]]$data==2] <- 3
            }
        }
        n.ind <- nind(cross)
        n.chr <- nchr(cross)
        if(verbose) {
            cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
            cat("INFO: Number of individuals: ",n.ind,"\n")
            cat("INFO: Number of chromosomes: ",n.chr,"\n")
        }
        savecross <- cross
        geno <- NULL
        chr <- NULL
        dist <- NULL
        newcmbase <- NULL
        out.qtl <- NULL
        for(i in 1:n.chr) {
            #We always shift the marker positions to 0
            geno <- cbind(geno,cross$geno[[i]]$data)
            chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
            newcmbase = c(newcmbase,min(cross$geno[[i]]$map))
            cross$geno[[i]]$map <- cross$geno[[i]]$map-(min(cross$geno[[i]]$map))
            dist <- c(dist,cross$geno[[i]]$map)
        }
        if(cofactor.significance <=0 || cofactor.significance >= 1){
            stop("cofactor.significance must be between 0 and 1.\n")
        }
        if(any(is.na(geno))){
            stop("Missing genotype information, please estimate unknown data, before running mqmscan.\n")
        }
        if(missing(cofactors)) cofactors <- rep(0,sum(nmar(cross)))
        numcofold <- sum(cofactors)
        cofactors <- checkdistances(cross,cofactors,1)
        numcofnew <- sum(cofactors)
        if(numcofold!=numcofnew){
            cat("INFO: Removed ",numcofold-numcofnew," cofactors that were close to eachother\n")
        }
        #Convert phenotypes from character to numeric
        pheno.col = stringPhenoToInt(cross,pheno.col)

        if (length(pheno.col) > 1){
            cross$pheno <- cross$pheno[,pheno.col]   #Scale down the triats
            result <- mqmscanall( cross,cofactors=cofactors,forceML=forceML,model=model,
                                 cofactor.significance=cofactor.significance,step.size=step.size,window.size=window.size,
                                 logtransform=logtransform, estimate.map = estimate.map,plot=plot, verbose=verbose,n.clusters=n.clusters,batchsize=batchsize)
            return(result)
        }
        if(pheno.col != 1){
            if(verbose) {
                cat("INFO: Selected phenotype ",pheno.col,".\n")
                cat("INFO: Number of phenotypes in object ",nphe(cross),".\n")
            }
            if(nphe(cross) < pheno.col || pheno.col < 1){
                stop("No such phenotype in cross object.\n")
            }
        }
        if(test.normality && !mqmtestnormal(cross, pheno.col, 0.01, FALSE)){
            warning("Trait might not be normal (Shapiro normality test)\n")
        }
        pheno <- cross$pheno[[pheno.col]]
        phenovar = var(pheno,na.rm = TRUE)
        if(phenovar > 1000){
            if(!logtransform){
                if(verbose) cat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"variation:",var(pheno,na.rm = TRUE),".\n")
                #warning(paste("WARNING: Set needs Log-transformation? (var=",phenovar,"), use mqmscan with logtransform=TRUE"))
            }
        }
        if(logtransform){
            #transform the cross file
            cross <- transformPheno(cross,pheno.col,transf=log)
            pheno <- cross$pheno[[pheno.col]]
        }
        n.mark <- ncol(geno)
        if(verbose) cat("INFO: Number of markers:",n.mark,"\n")

        #check for missing phenotypes
        dropped <- NULL
        droppedIND <- NULL
        for(i in 1:length(pheno)) {
            if(is.na(pheno[i]) || is.infinite(pheno[i])){
                if(verbose) cat("INFO: Dropped individual ",i," with missing phenotype.\n")
                dropped <- c(dropped,i)
                if(!is.null(cross$mqm)){
                    droppedIND <- c(droppedIND,cross$mqm$augIND[i])
                }
                n.ind = n.ind-1
            }
        }
        #Throw out missing phenotypes from phenotype vector and genotype matrix
        if(!is.null(dropped)){
            geno <- geno[-dropped,]
            pheno <- pheno[-dropped]
        }

        #CHECK for previously augmented dataset
        if(!is.null(cross$mqm)){
            augmentedNind <- cross$mqm$Nind
            augmentedInd <- cross$mqm$augIND
            if(verbose){
                cat("n.ind:",n.ind,"\n")
                cat("augmentedNind:",augmentedNind,"\n")
                cat("length(augmentedInd):",length(augmentedInd),"\n")
            }
            if(!is.null(dropped)){
                augmentedInd <- cross$mqm$augIND[-dropped]
                augmentedInd[1] <- 0
                for(x in 1:(length(augmentedInd)-1)){
                    if(augmentedInd[x+1] - 1 > augmentedInd[x] ){
                        for(y in (x+1):length(augmentedInd)){
                            augmentedInd[y] <- augmentedInd[y] - ((augmentedInd[x+1] - augmentedInd[x])-1)
                        }
                    }
                }
            }
            augmentedNind <- length(unique(augmentedInd))
            if(verbose){
                cat("New augmentedNind:",augmentedNind,"\n")
                cat("New length(augmentedInd):",length(augmentedInd),"\n")
            }
        }else{
            #No augmentation
            augmentedNind <- n.ind
            augmentedInd <- 0:n.ind
        }

        #CHECK if we have cofactors, so we can do backward elimination
        backward <- 0;
        if(missing(cofactors)){
            if(verbose) cat("INFO: No cofactors, setting cofactors to 0\n")
            cofactors = rep(0,n.mark)
        }else{
            if(length(cofactors) != n.mark){
                if(verbose) cat("ERROR: # Cofactors != # Markers\n")
                }else{
                    if(verbose) cat("INFO:",sum(cofactors!=0),"Cofactors received to be analyzed\n")
                    if((sum(cofactors) > n.ind-10 && dominance==0)){
                        stop("Cofactors don't look okay for use without dominance\n")
                    }
                    if((sum(cofactors)*2 > n.ind-10 && dominance==1)){
                        stop("Cofactors don't look okay for use with dominance\n")
                    }
                    if(sum(cofactors) > 0){
                        if(verbose) cat("INFO: Doing backward elimination of selected cofactors.\n")
                        backward <- 1;
                        n.run <- 0;
                    }else{
                        backward <- 0;
                        cofactors = rep(0,n.mark)
                    }
                }
        }
        step.min = -off.end;
        step.max = max(dist)+step.size+off.end;
        step.max <- as.integer(ceiling((step.max+step.size)/step.size)*step.size)
        if((step.min+step.size) > step.max){
            stop("step.max needs to be >= step.min + step.size")
        }
        #    cat("Step.min:",step.min," Step.max:",step.max,"\n")
        if(step.size < 1){
            stop("step.size needs to be >= 1")
        }
        qtlAchromo <- length(seq(step.min,step.max,step.size))
        if(verbose) cat("INFO: Number of locations per chromosome: ",qtlAchromo, "\n")
        end.1 <- proc.time()
        result <- .C("R_mqmscan",
                     as.integer(n.ind),
                     as.integer(n.mark),
                     as.integer(1),    # 1 phenotype
                     as.integer(geno),
                     as.integer(chr),
                     DIST=as.double(dist),
                     as.double(pheno),
                     COF=as.integer(cofactors),
                     as.integer(backward),
                     as.integer(forceML),
                     as.double(cofactor.significance),
                     as.integer(em.iter),
                     as.double(window.size),
                     as.double(step.size),
                     as.double(step.min),
                     as.double(step.max),
                     as.integer(n.run),
                     as.integer(augmentedNind),
                     as.integer(augmentedInd),
                     QTL=as.double(rep(0,2*n.chr*qtlAchromo)),
                     as.integer(estimate.map),
                     as.integer(ctype),
                     as.integer(dominance),
                     as.integer(verbose),
                     PACKAGE="qtl")
        end.2 <- proc.time()
        # initialize output object
        qtl <- NULL
        info <- NULL
        names <- NULL
        for(i in 1:(n.chr*qtlAchromo)) {
            #Store the result in the qtl object
            qtl <- rbind(qtl,c(ceiling(i/qtlAchromo),rep(seq(step.min,step.max,step.size),n.chr)[i],result$QTL[i]))
            info <- rbind(info,result$QTL[(n.chr*qtlAchromo)+i])
            #make names in the form: cX.locXX
            names <- c(names,paste("c",ceiling(i/qtlAchromo),".loc",rep(seq(step.min,step.max,step.size),n.chr)[i],sep=""))
        }
        if(estimate.map){
            new.map <- pull.map(cross)
            chrmarkers <- nmar(cross)
            sum <- 1
            for(i in 1:length(chrmarkers)) {
                for(j in 1:chrmarkers[[i]]) {
                    new.map[[i]][j] <- result$DIST[sum]
                    sum <- sum+1
                }
            }
        }
        if(plot){
            if(estimate.map && backward){
                op <- par(mfrow = c(3,1))
            }else{
                if(estimate.map || backward){
                    op <- par(mfrow = c(2,1))
                }else{
                    op <- par(mfrow = c(1,1))
                }
            }
            if(estimate.map){
                if(verbose) cat("INFO: Viewing the user supplied map versus genetic map used during analysis.\n")
                plotMap(pull.map(cross), new.map,main="Supplied map versus re-estimated map")
            }
        }
        if(backward){
            if(!estimate.map){
                new.map <- pull.map(cross)
            }
            chrmarkers <- nmar(cross)
            mapnames <- NULL
            for(x in 1:nchr(cross)){
                mapnames <- c(mapnames,names(pull.map(cross)[[x]]))
            }
            sum <- 1
            model.present <- 0
            qc <- NULL
            qp <- NULL
            qn <- NULL
            for(i in 1:length(chrmarkers)) {
                for(j in 1:chrmarkers[[i]]) {
                    #cat("INFO ",sum," ResultCOF:",result$COF[sum],"\n")
                    if(result$COF[sum] != 48){
                        if(verbose) cat("MODEL: Marker",sum,"named:", strsplit(names(unlist(new.map)),".",fixed=TRUE)[[sum]][2],"from model found, CHR=",i,",POSITION=",as.double(unlist(new.map)[sum])," cM\n")
                        qc <- c(qc, as.character(names(cross$geno)[i]))
                        qp <- c(qp, as.double(unlist(new.map)[sum]))
                        qn <- c(qn, mapnames[sum])
                        model.present <- 1
                    }
                    sum <- sum+1
                }
            }
            if(!is.null(qc) && model.present){
                why <- sim.geno(savecross,n.draws=1)
                QTLmodel <- makeqtl(why, qc, qp, qn, what="draws")
                attr(QTLmodel,"mqm") <- 1
                if(plot) plot(QTLmodel)
            }
        }
        rownames(qtl) <- names

        qtl <- cbind(qtl,1/(min(info))*(info-min(info)))
        qtl <- cbind(qtl,1/(min(info))*(info-min(info))*qtl[,3])
        colnames(qtl) = c("chr","pos (cM)",paste("LOD",colnames(cross$pheno)[pheno.col]),"info","LOD*info")
        #Convert to data/frame and scanone object so we can use the standard plotting routines
        qtl <- as.data.frame(qtl, stringsAsFactors=TRUE)
        if(backward && !is.null(qc) && model.present){
            attr(qtl,"mqmmodel") <- QTLmodel
            cimcovar <- as.data.frame(cbind(as.numeric(attr(qtl,"mqmmodel")[[4]]),as.data.frame(attr(qtl,"mqmmodel")[[5]])), stringsAsFactors=TRUE)
            rownames(cimcovar) <- attr(qtl,"mqmmodel")[[2]]
            colnames(cimcovar) <- c("chr","pos")
            attr(qtl, "marker.covar.pos") <- cimcovar
        }
        class(qtl) <- c("scanone",class(qtl))
        if(outputmarkers){
            #Remove pseudomarkers from the dataset and scale to the chromosome
            #Somewhat longer then off-end to be able to put back the original markers
            for( x in 1:nchr(cross)){
                to.remove <- NULL
                chr.length <- max(cross$geno[[x]]$map)
                markers.on.chr <- which(qtl[,1]==x)
                to.remove <- markers.on.chr[which(qtl[markers.on.chr,2] > chr.length+off.end+(2*step.size))]
                to.remove <- c(to.remove,markers.on.chr[which(qtl[markers.on.chr,2] < -off.end)])
                if(length(to.remove) > 0) qtl <- qtl[-to.remove,]
            }
            qtl <- addmarkerstointervalmap(cross,qtl)
            qtl <- as.data.frame(qtl, stringsAsFactors=TRUE)
            if(backward && !is.null(qc) && model.present){
                attr(qtl,"mqmmodel") <- QTLmodel
                cimcovar <- as.data.frame(cbind(as.numeric(attr(qtl,"mqmmodel")[[4]]),as.data.frame(attr(qtl,"mqmmodel")[[5]])), stringsAsFactors=TRUE)
                rownames(cimcovar) <- attr(qtl,"mqmmodel")[[2]]
                colnames(cimcovar) <- c("chr","pos")
                attr(qtl, "marker.covar.pos") <- cimcovar
            }
            class(qtl) <- c("scanone",class(qtl))
        }
        #Now we handle the off-end and any shifts comming from that (newcmbase)
        #We didn't thow away the old names, and could thus get duplicates
        qtlnew <- qtl
        rownames(qtlnew) <- NULL
        oldnames <- rownames(qtl)
        for(x in 1:n.chr){
            #Do per chromosome
            markers.on.chr <- which(qtl[,1]==x)
            if(newcmbase[x] !=0 && nrow(qtl[markers.on.chr,]) > 0){
                qtl[markers.on.chr,2] <- qtl[markers.on.chr,2]+newcmbase[x]
                #We end up with 3 posibilities
                for(n in 1:nrow(qtl[markers.on.chr,])){
                    name <- rownames(qtl[markers.on.chr,])[n]
                    id <- which(name==oldnames)
                    if(!is.na(name) && grepl(".loc",name,fixed=TRUE)){
                        #a marker we need to shift
                        rownames(qtlnew)[id] <- paste(strsplit(name,".loc",fixed=TRUE)[[1]][1],".loc",qtlnew[id,2],sep="")
                    }else{
                        #a marker with a user defined name, no need to shift
                        rownames(qtlnew)[id] <- name
                    }
                }
            }else{
                #No shift for the chromosome so we don't have to worry about shifting
                rownames(qtlnew)[markers.on.chr] <- rownames(qtl)[markers.on.chr]
            }
        }
        qtl <- qtlnew

        for( x in 1:nchr(cross)){
            #Remove the off ends that we left hanging a step before
            if(class(cross$geno[[x]])!="X"){
                to.remove <- NULL
                chr.length <- max(cross$geno[[x]]$map)
                markers.on.chr <- which(qtl[,1]==x)
                to.remove <- markers.on.chr[which(qtl[markers.on.chr,2] > chr.length+off.end)]
                to.remove <- c(to.remove,markers.on.chr[which(qtl[markers.on.chr,2] < -off.end)])
                qtl <- qtl[-to.remove,]
            }
        }

        qtl[,1] <- factor(names(cross$geno)[qtl[,1]], levels=names(cross$geno))

        #Plot the results if the user asked for it
        if(plot){
            info.c <- qtl
            #Check for errors in the information content IF err we can't do a second plot
            e <- 0
            for(i in 1:ncol(qtl)){
                if(is.na(info.c[i,5])){
                    e<- 1
                }
                if(is.infinite(info.c[i,5])){
                    e<- 1
                }
                if(is.null(info.c[i,5])){
                    e<- 1
                }
            }
            #No error do plot 2
            if(!e){
                mqmplot.singletrait(qtl,main=paste(colnames(cross$pheno)[pheno.col],"at significance=",cofactor.significance))
            }else{
                plot(qtl,main=paste(colnames(cross$pheno)[pheno.col],"at significance=",cofactor.significance),lwd=1)
                grid(max(qtl$chr),5)
                labels <- paste("QTL",colnames(cross$pheno)[pheno.col])
                legend("topright", labels,col=c("black"),lty=c(1))
            }
            op <- par(mfrow = c(1,1))
        }

        end.3 <- proc.time()
        if(verbose) cat("INFO: Calculation time (R->C,C,C-R): (",round((end.1-start)[3], digits=3), ",",round((end.2-end.1)[3], digits=3),",",round((end.3-end.2)[3], digits=3),") (in seconds)\n")

        qtl
    }else{
        stop("Currently only F2, BC, and selfed RIL crosses can be analyzed by MQM.")
    }
}

# end of mqmscan.R
