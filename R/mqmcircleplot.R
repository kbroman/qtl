#####################################################################
#
# mqmcircleplot.R
#
# Copyright (c) 2009-2011, Danny Arends
#
# Modified by Pjotr Prins and Karl Broman
#
#
# first written Februari 2009
# last modified May 2011
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
# Contains: mqmplot_c
#           mqmplotcircle
#           mqmplot_circle
#           circlelocations
#           drawspline
#           getchromosomelength
#           getgenomelength
#           locationtocircle
#           drawcirculargenome
#           loopthroughmulti
#
#
#####################################################################

mqmplot.circle <- function(cross, result, highlight=0, spacing=25, interactstrength=2, axis.legend = TRUE, col.legend=FALSE, verbose=FALSE, transparency=FALSE){
    if(is.null(cross)){
        stop("No cross object. Please supply a valid cross object.")
    }
    retresults <- NULL
    lod<-FALSE
    if(any(class(result)=="mqmmulti")){
        if(highlight > 0){
            templateresult <- mqmextractmarkers(result[[highlight]])
            lod<-TRUE
        }else{
            templateresult <- mqmextractmarkers(result[[1]])
            lod<-FALSE
        }
    }else{
        templateresult <- mqmextractmarkers(result)
        lod<-TRUE
    }
    if(!("scanone" %in% class(templateresult))){
        stop("Wrong type of result file, please supply a valid scanone object.")
    }
    if(transparency){
        colorz <- rainbow(length(result), alpha=0.8)
    }else{
        colorz <- rainbow(length(result))
    }
    if(!col.legend){
        totallength <- getgenomelength(templateresult)
        nchr <- length(unique(templateresult[,1]))
        cvalues <- circlelocations(totallength+(nchr*spacing))

        drawcirculargenome(templateresult,spacing=spacing,lodmarkers=lod)

        if(any(class(result)=="mqmmulti")){
            #multiple scan results, so plot em all
            todo <- 1:length(result)
            if(highlight!=0){
                #Unless we highlight only one of them
                todo <- highlight
            }
            for(x in todo){
                model <- mqmgetmodel(result[[x]])
                if(!is.null(model)){
                    if(verbose) cat("Model found for trait ",x,"\n")
                    if(!is.null(cross$locations)){
                        if(verbose) cat("Locations of traits available\n")
                        location <- as.numeric(cross$locations[[x]])
                        traitl <- locationtocircle(templateresult,location[1],location[2],spacing=spacing)
                    }else{
                        if(verbose) cat("Locations of traits not available\n")
                        traitl <- t(c(0,0+0.25*(0.5-(x/length(result)))))
                    }
                    if(highlight==0){
                        col=colorz[x]
                    }else{
                        col=rgb(0.1, 0.1, 0.1, 0.1)
                        if(highlight==x) title(main = paste("Circleplot of:",colnames(cross$pheno)[x]))
                    }
                    for(y in 1:length(model[[4]])){
                        qtll <- locationtocircle(templateresult,model[[4]][y],model[[5]][y],spacing=spacing)
                        if(highlight==x){
                            for(z in y:length(model[[4]])){
                                if(!z==y){
                                    cross <- sim.geno(cross)
                                    eff <- effectplot(cross,pheno.col=x,mname1=model$name[y],mname2=model$name[z],draw=FALSE)
                                    changeA <- (eff$Means[1,2]-eff$Means[1,1])
                                    changeB <- (eff$Means[2,2]-eff$Means[2,1])
                                    #interaction
                                    qtl2 <- locationtocircle(templateresult,model[[4]][z],model[[5]][z],spacing=spacing)
                                    if(!is.na(changeA) && !is.na(changeB) && !any(is.na(eff$SEs))){
                                        col <- "blue"
                                        if(changeA/abs(changeA) <  changeB/abs(changeB)) col <- "green"
                                        if(abs(abs(changeA)-abs(changeB)) > interactstrength*mean(eff$SEs)){
                                            retresults <- rbind(retresults,c(model$name[y],model$name[z],changeA,changeB,mean(eff$SEs)))
                                            drawspline(qtll,qtl2,lwd=2,col=col)
                                        }
                                    }
                                }
                            }
                        }
                        if(highlight==0){
                            points(traitl,col=col,pch=24,cex=1)
                            drawspline(traitl,qtll,col=col)
                            points(qtll*(1+0.1*((x/length(result)))),col=col,pch=19,cex=1)
                        }
                    }
                }else{
                    if(verbose) cat("Trait ",x," has no model\n")
                }
            }
            if(axis.legend) legend("topleft",c("Trait","QTL"),col=c("black","black"),pch=c(24,19),cex=1)
        }
        if(any(class(result)=="scanone") || highlight > 0){
            #single scan result or highlighting one of the multiple
            if(!any(class(result)=="scanone")){
                #Just highlight the template result
                result <- templateresult
            }
            model <- mqmgetmodel(result)
            if(!is.null(model)){
                for(y in 1:length(model[[4]])){
                    qtll <- locationtocircle(templateresult,model[[4]][y],model[[5]][y],spacing=spacing)
                    points(qtll,col="red",pch=19,cex=1)
                    text(qtll*1.15,model[[2]][y],col="red",cex=0.7)
                    if(!is.null(cross$locations)){
                        location <- as.numeric(cross$locations[[highlight]])
                        traitl <- locationtocircle(templateresult,location[1],location[2],spacing=spacing)
                        points(traitl,col="red",lwd=2,pch=24,cex=1.5)
                    }else{
                        traitl <- c(0,0)
                    }
                    if(!(highlight>0))drawspline(traitl,qtll,col="red")
                }
            }
            if(axis.legend) legend("topright",c("Selected cofactor(s)","Epistasis (+)","Epistasis (-)"), col=c("red", "blue", "green"), pch=19, lwd=c(0,1,2), cex=0.75)
            if(axis.legend) legend("bottomright",c("LOD 3","LOD 6","LOD 9","LOD 12"), pch=19, lwd=0, pt.cex = c(1, 2, 3, 4), cex=0.75)
            if(highlight==0) title(sub = "Single trait")
        }
    }else{
        plot(c(-1,1), c(-1, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
        title(main = "Legend to circular genome plot")
        legend("center", paste(colnames(cross$pheno)), col=colorz, pch=19, cex=0.75)
    }
    if(!is.null(retresults)){
        colnames(retresults) <- c("Marker","Marker","Change","Change","SEs")
        retresults <- as.data.frame(retresults, stringsAsFactors=TRUE)
        return(invisible(retresults))
    }
}

circlelocations <- function(nt){
    medpoints <- matrix(nrow = nt, ncol = 2)
    phi <- seq(0, 2 * pi, length = (nt + 1))
    complex.circle <- complex(modulus = 1, argument = phi)
    for (j in 1:nt) {
        medpoints[j, ] <- c(Im(complex.circle[j]), Re(complex.circle[j]))
    }
    medpoints
}

drawspline <- function (cn1, cn2, lwd = 1,col="blue",...){
    x <- cbind(cn1[1],0,cn2[1])
    y <- cbind(cn1[2],0,cn2[2])
    r <- xspline(x, y, lty=1, shape=1, lwd=lwd, border=col,...)
}

getchromosomelength <- function(result, chr){
    l <- ceiling(max(result[which(result[,1]==chr),2]))
    l
}

getgenomelength <- function(result){
    l <- 1
    for(x in unique(result[,1])){
        l <- l + getchromosomelength(result,x)
    }
    l
}

locationtocircle <- function(result, chr, loc, spacing=50, fixoutofbounds=TRUE, verbose=FALSE){
    templateresult <- result
    totallength <- getgenomelength(result)
    nchr <- length(unique(templateresult[,1]))
    cvalues <- circlelocations(totallength+(nchr*spacing))
    l <- 1
    for(x in unique(templateresult[,1])){
        if(x==chr){
            if(loc < getchromosomelength(result,x)){
                return(t(cvalues[(l+loc),]))
            }else{
                if(verbose) cat("Location out of chromosome bounds",loc," ",getchromosomelength(result,x),"\n")
                if(fixoutofbounds) return(t(cvalues[(l+getchromosomelength(result,x)),]))
                stop(paste("Location out of chromosome bounds",loc," ",getchromosomelength(result,x),"\n"))
            }
        }
        l <- l + getchromosomelength(result,x) + spacing
    }
    stop("No such chromosome")
}

drawcirculargenome <- function(result,lodmarkers=FALSE,spacing=50){
    result <- mqmextractmarkers(result)
    plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", axes = FALSE, xlab = "", ylab = "")
    totallength <- getgenomelength(result)
    nchr <- length(unique(result[,1]))
    cvalues <- circlelocations(totallength+(nchr*spacing))
    l <- 1
    for(x in unique(result[,1])){
        #Draw chromosomes
        nl <- l+getchromosomelength(result,x)
        lines(cvalues[l:nl,],cex=0.01)
        l <- nl + spacing
    }
    for(x in 1:nrow(result)){
        #Draw markers
        if(lodmarkers){
            points(locationtocircle(result,result[x,1],result[x,2],spacing=spacing),pch=20,cex=min(c((result[x,3]),4)))
        }else{
            points(locationtocircle(result,result[x,1],result[x,2],spacing=spacing),pch=20)
        }
    }
    for(x in 1:nchr){
        chrnumberloc <- locationtocircle(result,x,getchromosomelength(result,x)/2,spacing=spacing)
        points(t(c(-1.1, -1.15)))
        points(t(c(-0.9, -1.15)))
        points(t(c(-0.7, -1.15)))
        text(t(c(-0.9, -1.0)),paste("Distances in cM"),cex=0.8)
        text(t(c(-1.1, -1.1)),paste("0 cM"),cex=0.7)
        text(t(c(-0.9, -1.1)),paste(round((totallength+(nchr*spacing))*(0.2/(2*pi)),digits=1),"cM"),cex=0.7)
        text(t(c(-0.7, -1.1)),paste(round((totallength+(nchr*spacing))*(0.4/(2*pi)),digits=1),"cM"),cex=0.7)
        text(0.9*chrnumberloc,paste("Chr",x),cex=0.8)

    }
}

loopthroughmulti <- function(cross,result,save=FALSE,spacing=100){
    n <- 1
    while(n <= length(result)){
        if(save) png(filename=paste("circleplotT",n,".png",sep=""),width=1024,height=768)
        mqmplot.circle(cross,result,spacing=spacing,highlight=n)
        if(save) dev.off()
        n <- n+1
    }
}
