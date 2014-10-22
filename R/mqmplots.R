#####################################################################
#
# mqmplots.R
#
# Copyright (c) 2009-2011, Danny Arends
# Copyright polyplot routine (c) 2009, Rutger Brouwer
#
# Modified by Pjotr Prins and Karl Broman
#
#
# first written Februari 2009
# last modified Feb 2011
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
# Contains: mqmplot.directedqtl
#           mqmplot.cistrans
#           addloctocross
#           polyplot
#           getThird
#           getChr
#           mqmplot.multitrait
#           mqmplot.permutations
#           mqmplot.singletrait
#
#
#
#####################################################################

mqmplot.directedqtl <- function(cross, mqmresult, pheno.col=1, draw = TRUE)
{
    if(is.null(cross)){
        stop("No cross object. Please supply a valid cross object.")
    }
    if(!is.null(cross$mqm$Nind)){
        stop("Augmented crossobject. Please supply the original unaugmented dataset.")
    }
    if(is.null(mqmresult)){
        stop("No mqmresult object. Please supply a valid scanone object.")
    }
    if(!any(class(mqmresult)=="scanone")){
        stop("No mqmresult object. Please supply a valid scanone object.")
    }
    onlymarkers <- mqmextractmarkers(mqmresult)
    eff <- effectscan(sim.geno(cross),pheno.col=pheno.col,draw=FALSE)
    if(any(eff[,1]=="X")){
        eff <- eff[-which(eff[,1]=="X"),]
    }
    onlymarkers[,3] <- onlymarkers[,3]*(eff[,3]/abs(eff[,3]))
    if(draw) plot(ylim=c((min(onlymarkers[,3])*1.1),(max(onlymarkers[,3])*1.1)),onlymarkers)
    class(onlymarkers) <- c("scanone",class(onlymarkers))
    if(!is.null(attr(mqmresult,"mqmmodel"))) attr(onlymarkers,"mqmmodel") <- attr(mqmresult,"mqmmodel")
    invisible(onlymarkers)
}

mqmplot.heatmap <- function(cross, result, directed=TRUE, legend=FALSE, breaks = c(-100,-10,-3,0,3,10,100), col = c("darkblue","blue","lightblue","yellow","orange","red"), ...)
{
    if(is.null(cross)){
        stop("No cross object. Please supply a valid cross object.")
    }
    if(directed && !is.null(cross$mqm$Nind)){
        stop("Augmented crossobject. Please supply the original unaugmented dataset.")
    }
    if(is.null(result)){
        stop("No result object. Please supply a valid scanone object.")
    }
    if(!any(class(result)=="mqmmulti")){
        stop("Not a mqmmulti object. Please supply a valid scanone object.")
    }
    cross <- sim.geno(cross)
    names <- NULL
    for(x in 1:nphe(cross)){
        result[[x]] <- mqmextractpseudomarkers(result[[x]])
        if(directed){
            effect <- effectscan(sim.geno(cross,step=stepsize(result[[x]])), pheno.col=x, draw=FALSE)
            for(y in 1:nrow(result[[x]])){
                effectid <- which(rownames(effect)==rownames(result[[x]])[y])
                if(!is.na(effectid&&1)){
                    result[[x]][y,3]  <- result[[x]][y,3] *(effect[effectid,3]/abs(effect[effectid,3]))
                }
            }
        }
        names <- c(names,substring(colnames(result[[x]])[3],5))
    }
    chrs <- unique(lapply(result,getChr))
    data <- NULL
    for(x in 1:length(result)){
        data <- rbind(data,result[[x]][,3])
    }
    rownames(data) <- names
    if(nphe(cross) < 100){
        image(seq(0,nrow(result[[1]])),seq(0,nphe(cross)),t(data),xlab="Markers",ylab="Traits",breaks=breaks,col=col,yaxt="n",...)
        axis(2,at=seq(1,nphe(cross)),labels=colnames(pull.pheno(cross)),las=1)
    }else{
        image(seq(0,nrow(result[[1]])),seq(0,nphe(cross)),t(data),xlab="Markers",ylab="Traits",breaks=breaks,col=col,...)
    }
    abline(v=0)
    for(x in unique(chrs[[1]])){
        abline(v=sum(as.numeric(chrs[[1]])<=x))
    }
    for(x in 1:nphe(cross)){
        abline(h=x)
    }
    if(legend){
        leg <- NULL
        for(x in 2:length(breaks)){
            leg <- c(leg,paste("LOD",breaks[x-1],"to",breaks[x]))
        }
        legend("bottom",leg,col=col,lty=1,bg="white")
    }
    invisible(data)
}

mqmplot.clusteredheatmap <- function(cross, mqmresult, directed=TRUE, legend=FALSE, Colv=NA, scale="none", verbose=FALSE, breaks = c(-100,-10,-3,0,3,10,100), col = c("darkblue","blue","lightblue","yellow","orange","red"), ...)
{
    if(is.null(cross)){
        stop("No cross object. Please supply a valid cross object.")
    }
    if(directed && !is.null(cross$mqm$Nind)){
        stop("Augmented crossobject. Please supply the original unaugmented dataset.")
    }
    if(is.null(mqmresult)){
        stop("No mqmresult object. Please supply a valid mqmmulti object.")
    }
    if(!any(class(mqmresult)=="mqmmulti")){
        stop("Not a mqmmulti object. Please supply a valid mqmmulti object.")
    }
    cross <- sim.geno(cross)
    names <- NULL
    for(x in 1:nphe(cross)){
        mqmresult[[x]] <- mqmextractpseudomarkers(mqmresult[[x]])
        if(directed){
            effect <- effectscan(sim.geno(cross,step=stepsize(mqmresult[[x]])), pheno.col=x, draw=FALSE)
            if(verbose) cat(".")
            for(y in 1:nrow(mqmresult[[x]])){
                effectid <- which(rownames(effect)==rownames(mqmresult[[x]])[y])
                if(!is.na(effectid&&1)){
                    mqmresult[[x]][y,3]  <- mqmresult[[x]][y,3] *(effect[effectid,3]/abs(effect[effectid,3]))
                }
            }
        }
        names <- c(names,substring(colnames(mqmresult[[x]])[3],5))
    }
    if(verbose && directed) cat("\n")
    chrs <- unique(lapply(mqmresult,getChr))
    data <- NULL
    for(x in 1:length(mqmresult)){
        data <- rbind(data,mqmresult[[x]][,3])
    }
    colnames(data) <- rownames(mqmresult[[1]])
    rownames(data) <- names
    if(length(names) < 100){
        retresults <- heatmap(data,Colv=Colv,scale=scale, xlab="Markers",main="Clustered heatmap",breaks=breaks,col=col,keep.dendro = TRUE, ...)
    }else{
        retresults <- heatmap(data,Colv=Colv,scale=scale, xlab="Markers",main="Clustered heatmap",breaks=breaks,col=col,keep.dendro = TRUE,labRow=1:length(names), ...)
    }
    if(legend){
        leg <- NULL
        for(x in 2:length(breaks)){
            leg <- c(leg,paste("LOD",breaks[x-1],"to",breaks[x]))
        }
        legend("bottom",leg,col=col,lty=1,bg="white")
    }
    invisible(retresults)
}

mqmplot.cistrans <- function(result,cross,threshold=5,onlyPEAK=TRUE,highPEAK=FALSE,cisarea=10,pch=22,cex=0.5, verbose=FALSE, ...)
{
    if(is.null(cross)){
        stop("No cross object. Please supply a valid cross object.")
    }
    if(is.null(cross$locations)){
        stop("Please add trait locations to the cross file\n")
    }
    if(any(class(result) == "mqmmulti")){
        sum.map <- 0
        chr.breaks <- NULL
        for(j in 1:nchr(cross)){
            l.chr <- max(result[[1]][result[[1]][,1]==j,2])
            chr.breaks <- c(chr.breaks,sum.map)
            sum.map <- sum.map+l.chr
        }
        sum.map <- ceiling(sum.map)
        if(verbose) cat("Total maplength:",sum.map," cM in ",nchr(cross),"Chromosomes\nThe lengths are:",chr.breaks,"\n")
        locations <-  do.call(rbind,cross$locations)
        QTLs <- do.call(rbind,lapply(FUN=getThird,result))
        colnames(QTLs) <- rownames(result[[1]])
        bmatrix <- QTLs>threshold
        pmatrix <- NULL
        for(j in 1:nrow(QTLs)){
            if(verbose && (j%%1000 == 0)) cat("QTL row:",j,"\n")
            temp <- as.vector(bmatrix[j,])
            tempv <- QTLs[j,]
            value = 0
            index = -1
            for(l in 1:length(temp)){
                if(temp[l]){
                    if( tempv[l] > value){
                        #New highest marker set the old index to false
                        if(index != -1){
                            temp[index] <- FALSE
                        }
                        #and store the new one
                        value = tempv[l]
                        index <- l
                    }else{
                        #Lower marker
                        temp[l] <- FALSE
                    }
                }else{
                    index = -1
                    value = 0
                }
            }
            if(onlyPEAK){
                bmatrix[j,] <- temp
            }
            if(highPEAK){
                pmatrix <- rbind(pmatrix,temp)
            }
        }

        locz <- NULL
        for(marker in 1:ncol(QTLs)){
            if(verbose && (j%%1000 == 0)) cat("QTL col:",marker,"\n")
            pos <- find.markerpos(cross, colnames(QTLs)[marker])
            if(!is.na(pos[1,1])){
                locz <- c(locz,round(chr.breaks[as.numeric(pos[[1]])] + as.numeric(pos[[2]])))
            }else{
                mark <- colnames(QTLs)[marker]
                mchr <- substr(mark,sum(regexpr("c",mark)+attr(regexpr("c",mark),"match.length")),regexpr(".loc",mark)-1)
                mpos <- as.numeric(substr(mark,sum(regexpr("loc",mark)+attr(regexpr("loc",mark),"match.length")),nchar(mark)))
                locz <- c(locz,round(chr.breaks[as.numeric(mchr)] + as.numeric(mpos)))
            }
        }

        axi <- 1:sum.map
        plot(x=axi,y=axi,type="n",main=paste("Cis/Trans QTL plot at LOD",threshold),xlab="Markers (in cM)",ylab="Location of traits (in cM)",xaxt="n",yaxt="n")
        trait.locz <- NULL
        for(j in 1:nrow(QTLs)){
            if(verbose && (j%%10 == 0)) cat("QTL row:",j,"\n")
            values <- rep(NA,sum.map)
            aa <- locz[bmatrix[j,]]
            trait.locz <- c(trait.locz,chr.breaks[locations[j,1]] + locations[j,2])
            values[aa] = chr.breaks[locations[j,1]] + locations[j,2]
            if(!highPEAK){
                points(values,pch=pch,cex=cex)
            }else{
                points(values,pch=24,cex=1.25*cex,col="black",bg="red")
            }
        }
        points(axi,type="l")
        points(axi+(cisarea/2),type="l",col="green")
        points(axi-(cisarea/2),type="l",col="green")
        chr.breaks <- c(chr.breaks,sum.map)
        axis(1,at=chr.breaks,labels=FALSE)
        axis(2,at=chr.breaks,labels=FALSE)
        axis(1,at=locz,line=1,pch=24)
        axis(2,at=seq(0,sum.map,25),line=1)
    }else{
        stop("invalid object supplied\n")
    }
}

addloctocross <- function(cross,locations=NULL,locfile="locations.txt", verbose=FALSE)
{
    if(is.null(cross)){
        stop("No cross object. Please supply a valid cross object.")
    }
    if(is.null(locations)){
        locations <- read.table(locfile,row.names=1,header=TRUE)
    }
    if(verbose) {
        cat("Phenotypes in cross:",nphe(cross),"\n")
        cat("Phenotypes in file:",nrow(locations),"\n")
    }
    if(max(as.numeric(rownames(locations))) != nphe(cross)){
        stop("ID's of traits in file are larger than # of traits in crossfile.")
    }
    if(nphe(cross)==nrow(locations)){
        locs <- vector(mode = "list", length = nphe(cross))
        for(x in as.numeric(rownames(locations))){
            if(names(cross$pheno)[x] == locations[x,1]){
                locs[[x]] <- locations[x,2:3]
                rownames(locs[[x]]) <- locations[x,1]
            }else{
                warning("Mismatch between name of trait in cross & file.\n")
            }
        }
    }else{
        stop("Number of traits in cross & file don't match.")
    }
    cross$locations <- locs
    cross
}


polyplot <- function( x, type='b', legend=TRUE,legendloc=0, labels=NULL, cex = par("cex"), pch = 19, gpch = 21, bg = par("bg"), color = par("fg"), col=NULL, ylim=range(x[is.finite(x)]), xlim = NULL,
                     main = NULL, xlab = NULL, ylab = NULL, add=FALSE, ... )
{
    #Addition by Danny Arends
    if(legend){
        if(legendloc){
            op <- par(mfrow = c(1,2))
        }else{
            op <- par(mfrow = c(1,1))
        }
    }else{
        op <- par(mfrow = c(1,1))
    }
    #End of addition
    if ( is.vector(x) ) {
        x <- t( as.matrix(x) )
        if (is.null(labels) ) {
            rownames(x) <- c("unspecified gene")
        } else {
            rownames(x) <- labels
        }
    } else {
        x <- as.matrix(x)
    }
    if (is.null(labels) ) labels = rownames( x )
    if (is.null(col) )  col = rainbow( nrow(x),alpha=0.35 )
    if (is.null(xlab) ) xlab="Markers"
    if (is.null(ylab) ) ylab="QTL (LOD)"

    timepoints <- as.numeric( colnames(x) )
    tps        <- sort( unique( timepoints ) )

    if(is.null(xlim)) xlim = c(min(tps),max(tps))
    plot.new()															# make a new plot
    plot.window(xlim=xlim, ylim=ylim, log="")							# add the plot windows size
    #grid()
    for( k in 1:nrow( x ) ) {
        max.p  <- NULL			# the expression of the maximum
        min.p  <- NULL			#
        med.p  <- NULL			#
        for( i in 1:length(tps)) {	#
            idx <- ( 1:length( timepoints ) )[timepoints==tps[i] ]	# get the indeces of the
            work  <- x[k, idx]
            pmax  <- max(work[is.finite(work)], na.rm=TRUE)
            pmin  <- min(work[is.finite(work)], na.rm=TRUE)
            pmed  <- median(work[is.finite(work)], na.rm=TRUE)

            max.p <- append(max.p, pmax)	#
            min.p <- append(min.p, pmin)	#
            med.p <- append(med.p, pmed)	#
        }
        lines( x=tps, y=max.p, type='l', col=col[k] )		# add the lines if requested
        lines( x=tps, y=min.p, type='l', col=col[k] )		# add the lines if requested
        xp <- append(tps, rev(tps))
        yp <- append(max.p, rev(min.p) )

        polygon(xp, y=yp, col=col[k], border=FALSE)
        lines( x=tps, y=med.p, type='l', col=col[k] )		# add the lines if requested
    }


    axis(1)																# add the x axis
    axis(2)																# add the y axis

    title(main=main, xlab=xlab, ylab=ylab, ...)							# add the title and axis labels
    #Addition by Danny Arends
    if (legend){
        if(legendloc){
            plot.new()
        }
        legend("topright", labels, col=col, pch=pch)			# add a legend if requested
    }
    #End of addition
    op <- par(mfrow = c(1,1))
    invisible()															# return the plot
}

getThird <- function(x){
    x[,3]
}

getChr <- function(x){
    x[,1]
}

mqmplot.multitrait <- function(result, type=c("lines","image","contour","3Dplot"), group=NULL, meanprofile=c("none","mean","median"), theta=30, phi=15, ...)
{
    #Helperfunction to plot mqmmulti objects made by doing multiple mqmscan runs (in a LIST)
    type <- match.arg(type)
    meanprofile <- match.arg(meanprofile)
    if(class(result)[2] != "mqmmulti"){
        stop("Wrong type of result file, please supply a valid mqmmulti object.")
    }
    n.pheno <- length(result)
    temp <- lapply(result,getThird)
    chrs <- unique(lapply(result,getChr))
    qtldata <- do.call("rbind",temp)
    if(!is.null(group)){
        qtldata <- qtldata[group,]
        colors <- rep("blue",n.pheno)
    }else{
        group <- 1:n.pheno
        colors <- rainbow(n.pheno)
    }
    qtldata <- t(qtldata)
    if(type=="contour"){
        #Countour plot
        contour(x=seq(1,dim(qtldata)[1]),
                y=seq(1,dim(qtldata)[2]),
                qtldata,
                xlab="Markers",ylab="Trait", ...)
        for(x in unique(chrs[[1]])){
            abline(v=sum(as.numeric(chrs[[1]])<=x))
        }
    }
    if(type=="image"){
        #Image plot
        image(x=1:dim(qtldata)[1],y=1:dim(qtldata)[2],qtldata,xlab="Markers",ylab="Trait",...)
        for(x in unique(chrs[[1]])){
            abline(v=sum(as.numeric(chrs[[1]])<=x))
        }
    }
    if(type=="3Dplot"){
        #3D perspective plot
        persp(x=1:dim(qtldata)[1],y=1:dim(qtldata)[2],qtldata,
              theta = theta, phi = phi, expand = 1,
              col="gray", xlab = "Markers", ylab = "Traits", zlab = "LOD score")
    }
    if(type=="lines"){
        #Standard plotting option, Lineplot
        first <- TRUE
        for(i in group) {
            if(first){
                plot(result[[i]],ylim=c(0,max(qtldata)),col=colors[i],lwd=1,ylab="LOD score",xlab="Markers",main="Multiple profiles", ...)
                first <- FALSE
            }else{
                plot(result[[i]],add=TRUE,col=colors[i],lwd=1,...)
            }
        }
        if(meanprofile != "none"){
            temp <- result[[1]]
            if(meanprofile=="median"){
                temp[,3] <- apply(qtldata,1,median)
                legend("topright",c("QTL profiles","Median profile"),col=c("blue","black"),lwd=c(1,3))
            }
            if(meanprofile=="mean"){
                temp[,3] <- rowMeans(qtldata)
                legend("topright",c("QTL profiles","Mean profile"),col=c("blue","black"),lwd=c(1,3))
            }
            plot(temp,add=TRUE,col="black",lwd=3,...)
        }
    }
}

mqmplot.permutations <- function(permutationresult, ...)
{
    #Helperfunction to show mqmmulti objects made by doing multiple mqmscan runs (in a LIST)
    #This function should only be used for bootstrapped data
    matrix <- NULL
    row1 <- NULL
    row2 <- NULL
    i <- 1
    if(class(permutationresult)[2] != "mqmmulti")
        ourstop("Wrong type of result file, please supply a valid mqmmulti object.")

    for( j in 1:length( permutationresult[[i]][,3] ) ) {
        row1 <- NULL
        row2 <- NULL
        for(i in 1:length( permutationresult ) ) {
            if(i==1){
                row1 <- c(row1,rep(permutationresult[[i]][,3][j],(length( permutationresult )-1)))
                names(row1) <- rep(j,(length( permutationresult )-1))
            }else{
                row2 <- c(row2,permutationresult[[i]][,3][j])
            }
        }
        names(row2) <- rep(j,(length( permutationresult )-1))
        matrix <- cbind(matrix,rbind(row1,row2))
    }

    rownames(matrix) <- c("QTL trait",paste("# of bootstraps:",length(permutationresult)-1))

    #Because bootstrap only has 2 rows of data we can use black n blue
    polyplot(matrix,col=c(rgb(0,0,0,1),rgb(0,0,1,0.35)),...)
    #PLot some lines so we know what is significant
    perm.temp <- mqmprocesspermutation(permutationresult)			#Create a permutation object
    numresults <- dim(permutationresult[[1]])[1]
    lines(x=1:numresults,y=rep(summary(perm.temp)[1,1],numresults),col="green",lwd=2,lty=2)
    lines(x=1:numresults,y=rep(summary(perm.temp)[2,1],numresults),col="blue",lwd=2,lty=2)
    chrs <- unique(lapply(permutationresult,getChr))
    for(x in unique(chrs[[1]])){
        abline(v=sum(as.numeric(chrs[[1]])<=x),lty="dashed",col="gray",lwd=1)
    }
}

mqmplot.singletrait <- function(result, extended=0,...)
{
    #Helperfunction to show scanone objects made by doing mqmscan runs
    if(!("scanone" %in% class(result))){
        stop("Wrong type of result file, please supply a valid scanone (from MQM) object.")
    }
    if(is.null(result$"info")){
        stop("Wrong type of result file, please supply a valid scanone (from MQM) object.")
    }
    if(is.null(attr(result,"mqmmodel"))){
        op <- par(mfrow = c(1,1))
    }else{
        op <- par(mfrow = c(2,1))
    }
    info.l <- result
    info.l[,3] <- result[,4]
    if(extended){
        if(!is.null(attr(result,"mqmmodel"))){
            plot(attr(result,"mqmmodel"))
        }
        plot(result,lwd=1,col=c("black"),ylab="QTL (LOD)",...)
        par(new=TRUE)
        plot(info.l,lwd=1,col=c("red"),ylab="QTL (LOD)",yaxt="n",lty=1,...)
        grid(length(result$chr),5)
        labels <- c(colnames(result)[3],"Information Content")
        mtext("Information Content",side=4,col="red",line=4)
        axis(4, ylim=c(0,1), col="red",col.axis="red",las=1)

        legend("right", labels,col=c("black","red"),lty=c(1,1,1))
    }else{
        if(!is.null(attr(result,"mqmmodel"))){
            plot(attr(result,"mqmmodel"))
        }
        plot(result,lwd=1,ylab="QTL (LOD)",...)
        grid(length(result$chr),5)
        labels <- c(colnames(result)[3])
        legend("right", labels,col=c("black","blue"),lty=c(1,1))
    }
    op <- par(mfrow = c(1,1))
}

# end of mqmplots.R
