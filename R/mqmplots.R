#####################################################################
#
# MQMplots.R
#
# copyright polyplot (c) 2009, Rutger Brouwer
# copyright other functions (c) 2009, Danny Arends
# last modified Jun, 2009
# first written Feb, 2009
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
# Contains: Different plotting functions belonging to the MQM subpackage of R/QTL
#           polyplot, plotMQMall, plotMQMboot, plotMQMnice, plotMQMone
#           CisTransPlot, addloctocross 
######################################################################

CisTransPlot <- function(x,cross,threshold=5,onlyPEAK=TRUE,highPEAK=FALSE,cisarea=10,pch=22,cex=0.5, ...){
	if(is.null(cross$locations)){
		stop("Please add trait locations to the cross file\n")
	}
	locations <- NULL
	if(any(class(x) == "MQMmulti")){
		sum.map <- 0
		chr.breaks <- NULL
		for(j in 1:nchr(cross)){
			l.chr <- max(x[[1]][x[[1]][,1]==j,2])
			chr.breaks <- c(chr.breaks,sum.map)
			sum.map <- sum.map+l.chr
		}
		sum.map <- ceiling(sum.map)
		cat("Total maplength:",sum.map," cM in ",nchr(cross),"Chromosomes\nThe lengths are:",chr.breaks,"\n")		
		for( k in 1:length(x) ) {
			loc <- cross$locations[[k]]
			rownames(loc) <- k
			locations <- rbind(locations,loc)
		}
		QTLs <- NULL
		for(y in 1:nrow(locations)){
			qtl <- x[[y]][,3]
			QTLs <- rbind(QTLs,qtl)
		}
		colnames(QTLs) <- rownames(x[[1]])
		axi <- 1:sum.map
		plot(x=axi,y=axi,type="n",main="Cis/Trans QTLplot",sub=paste("QTLs above threshold:",threshold,"LOD"),xlab="Markers (in cM)",ylab="Location of traits (in cM)",xaxt="n",yaxt="n")
		bmatrix <- QTLs>threshold
		pmatrix <- NULL
		for(j in 1:nrow(QTLs)){
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
			pos <- find.markerpos(cross, colnames(QTLs)[marker])
			locz <- c(locz,round(chr.breaks[as.numeric(pos[[1]])] + as.numeric(pos[[2]])))
		}
		trait.locz <- NULL
		for(j in 1:nrow(QTLs)){
			values <- rep(NA,sum.map)
			aa <- locz[bmatrix[j,]]
			trait.locz <- c(trait.locz,chr.breaks[locations[j,1]] + locations[j,2])
			values[aa] = chr.breaks[locations[j,1]] + locations[j,2]
			points(values,pch=pch,cex=cex)
		}
		if(highPEAK){
		for(j in 1:nrow(QTLs)){
			values <- rep(NA,sum.map)
			aa <- locz[pmatrix[j,]]
			trait.locz <- c(trait.locz,chr.breaks[locations[j,1]] + locations[j,2])
			values[aa] = chr.breaks[locations[j,1]] + locations[j,2]
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

addloctocross <- function(cross,locations=NULL,locfile="locations.txt"){
	if(is.null(locations)){
		locations <- read.table(locfile,row.names=1,header=TRUE)
	}
	cat("Phenotypes in cross:",nphe(cross),"\n")
	cat("Phenotypes in file:",nrow(locations),"\n")
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
					  main = NULL, xlab = NULL, ylab = NULL, add=FALSE, ... ){
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
	grid()
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

plotMQMall <- function(result, type="C", theta=30, phi=15, ...){
	#Helperfunction to plot MQMmulti objects made by doing multiple mqmscan runs (in a LIST)
  if(class(result)[2] != "MQMmulti")
		stop("Wrong type of result file, please supply a valid MQMmulti object.") 

  if(type=="C"){
    #Countour plot
    temp <- lapply(result,getThird)
	c <- do.call("rbind",temp)
    c <- t(c)
    contour(
            x=seq(1,dim(c)[1]),
            y=seq(1,dim(c)[2]),
            c,
            xlab="Markers",ylab="Trait",
            col=rainbow((max(c)/5)+25,1,1.0,0.1),
            nlevels=(max(c)/5)
            )
  }
  if(type=="I"){
    #Image plot
    temp <- lapply(result,getThird)
	c <- do.call("rbind",temp)
	c <- t(c)
    image(x=1:dim(c)[1],y=1:dim(c)[2],c,
          xlab="Markers",ylab="Trait",
          col=rainbow((max(c)/5)+25,1,1.0,0.1),
          )
    
  }
  if(type=="D"){
    #3D perspective plot
    temp <- lapply(result,getThird)
	c <- do.call("rbind",temp)
    c <- t(c)
    persp(x=1:dim(c)[1],y=1:dim(c)[2],c,
          theta = theta, phi = phi, expand = 1,
          col="gray", xlab = "Markers", ylab = "Traits", zlab = "LOD score")
  }
  if(type=="P"){
    #Standard plotting option, Lineplot
    n.pheno <- length(result)
    colors <- rainbow(n.pheno)
    for(i in 1:n.pheno) {
      if(i !=1 ){
        plot(result[[i]],add=TRUE,col=colors[i],lwd=1,...)
      }else{
        plot(result[[i]],col="black",lwd=1,...)
      }
    }
  }

}

plotMQMboot <- function(result, ...){
	#Helperfunction to show MQMmulti objects made by doing multiple mqmscan runs (in a LIST)
	#This function should only be used for bootstrapped data
	matrix <- NULL
	row1 <- NULL
	row2 <- NULL
	i <- 1
	if(class(result)[2] != "MQMmulti")		
          ourstop("Wrong type of result file, please supply a valid MQMmulti object.") 

        for( j in 1:length( result[[i]][,3] ) ) {
          row1 <- NULL
          row2 <- NULL
          for(i in 1:length( result ) ) {
            if(i==1){
              row1 <- c(row1,rep(result[[i]][,3][j],(length( result )-1)))
              names(row1) <- rep(j,(length( result )-1))
            }else{
              row2 <- c(row2,result[[i]][,3][j])
            }
          }
          names(row2) <- rep(j,(length( result )-1))
          matrix <- cbind(matrix,rbind(row1,row2))
        }

        rownames(matrix) <- c("QTL trait",paste("# of bootstraps:",length(result)-1))
		
	#Because bootstrap only has 2 rows of data we can use black n blue
        polyplot(matrix,col=c(rgb(0,0,0,1),rgb(0,0,1,0.35)),...)
	#PLot some lines so we know what is significant
        perm.temp <- MQMpermObject(result)			#Create a permutation object
        numresults <- dim(result[[1]])[1]
        lines(x=1:numresults,y=rep(summary(perm.temp)[1,1],numresults),col="green",lwd=2,lty=2)
        lines(x=1:numresults,y=rep(summary(perm.temp)[2,1],numresults),col="blue",lwd=2,lty=2)			
}

plotMQMnice <- function(result, ...){
	#Helperfunction to show MQMmulti objects made by doing multiple mqmscan runs (in a LIST)
	matrix <- NULL
	names <- NULL
	i <- 1
	if(class(result)[2] != "MQMmulti")		
          ourstop("Wrong type of result file, please supply a valid MQMmulti object.") 

	for( j in 1:length( result[[i]][,3] ) ) {
		row <- NULL
		for( i in 1:length( result ) ) {
			row <- c(row,result[[i]][,3][j])
		}
		matrix <- rbind(matrix,row)
	}
	for( i in 1:length( result ) ) {
		if(colnames(result[[i]])[3] == "QTL phenotype"){
			if(i==1){
				names <- "QTL Phenotype"
			}else{
				names <- c(names,paste("Bootstrap run",(i-1)))
			}
		}else{
			names <- c(names,colnames(result[[i]])[3])
		}
	}
	matrix <- t(matrix)
	colnames(matrix) <- c(1:dim(matrix)[2])
	rownames(matrix) <- names
	polyplot(matrix,...)
}

plotMQMone <- function(result, result2, extended=0,...){
	#Helperfunction to show scanone objects made by doing mqmscan runs
	if(!("scanone" %in% class(result)))
          ourstop("Wrong type of result file, please supply a valid scanone (from MQM) object.") 

        info.c <- result
        info.c[,3]<- info.c[,5]
        if(extended){
          info.l <- result
          info.l[,3] <- result[,4]
          plot(result,info.c,info.l,lwd=1,col=c("black","blue","red"),ylab="QTL (LOD)",...)
          grid(max(result$chr),5)
          labels <- c(colnames(result)[3],colnames(result)[5],colnames(result)[4])
          legend("topright", labels,col=c("black","blue","red"),lty=c(1,1,1))		
        }else{
          if (!missing(result2) && !("scanone" %in% class(result2))) {
	    #MAX 3 scanone objects
            plot(result,info.c,result2,lwd=1,ylab="QTL (LOD)",...)
            grid(max(result$chr),5)
            labels <- c(colnames(result)[3],colnames(result)[5],colnames(result2)[3])
            legend("topright", labels,col=c("black","blue"),lty=c(1,1))
          }else{
	    #MAX 3 scanone objects (here we now have 2)
            plot(result,info.c,lwd=1,ylab="QTL (LOD)",...)
            grid(max(result$chr),5)
            labels <- c(colnames(result)[3],colnames(result)[5])
            legend("topright", labels,col=c("black","blue"),lty=c(1,1))			
          }
        }
}

# end of MQMplots.R
