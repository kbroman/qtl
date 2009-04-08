#####################################################################
#
# MQMplots.R
#
# copyright polyplot (c) 2009, Rutger Brouwer
# copyright other functions (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: scanMQM
#
######################################################################

polyplot <- function( x, type='b', legend=T,legendloc=0, labels=NULL, cex = par("cex"), pch = 19, gpch = 21, bg = par("bg"), color = par("fg"), col=NULL, ylim=range(x[is.finite(x)]), xlim = NULL, 
					  main = NULL, xlab = NULL, ylab = NULL, add=F, ... ){
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
		max_p  <- NULL			# the expression of the maximum
		min_p  <- NULL			#
		med_p  <- NULL			#
		for( i in 1:length(tps)) {	# 		
			idx <- ( 1:length( timepoints ) )[timepoints==tps[i] ]	# get the indeces of the
			work  <- x[k, idx]
			pmax  <- max(work[is.finite(work)], na.rm=T)
			pmin  <- min(work[is.finite(work)], na.rm=T)
			pmed  <- median(work[is.finite(work)], na.rm=T)
			
			max_p <- append(max_p, pmax)	# 
			min_p <- append(min_p, pmin)	#
			med_p <- append(med_p, pmed)	#
		}
		lines( x=tps, y=max_p, type='l', col=col[k] )		# add the lines if requested		
		lines( x=tps, y=min_p, type='l', col=col[k] )		# add the lines if requested
		xp <- append(tps, rev(tps)) 
		yp <- append(max_p, rev(min_p) )		
			
		polygon(xp, y=yp, col=col[k], border=F)
		lines( x=tps, y=med_p, type='l', col=col[k] )		# add the lines if requested
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

plot.MQMall <- function(result = NULL, type="C", theta=30, phi=15,...){
	#Helperfunction to plot MQMmulti objects made by doing multiple scanMQM runs (in a LIST)
	if(class(result)[2] == "MQMmulti"){
		if(type=="C"){
		#Countour plot
			c <- NULL
			for(i in 1:length(result)){
				#Collect all the "QTL PHENO_TYPE" colums of the result
				c <- rbind(c,result[[i]][,3])
			}
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
			c <- NULL
			for(i in 1:length(result)){
				c <- rbind(c,result[[i]][,3])
			}
			c <- t(c)
			image(x=1:dim(c)[1],y=1:dim(c)[2],c,
				  xlab="Markers",ylab="Trait",
				  col=rainbow((max(c)/5)+25,1,1.0,0.1),
			)
		
		}
		if(type=="D"){
		#3D perspective plot
			c <- NULL
			for(i in 1:length(result)){
				c <- rbind(c,result[[i]][,3])
			}
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
	}else{
		ourstop("Wrong type of result file, please supply a valid MQMmulti object.") 
	}
}

plot.MQMboot <- function(result = NULL,...){
	#Helperfunction to show MQMmulti objects made by doing multiple scanMQM runs (in a LIST)
	#This function should only be used for bootstrapped data
	matrix <- NULL
	row1 <- NULL
	row2 <- NULL
	i <- 1
	if(class(result)[2] == "MQMmulti"){		
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
		perm_temp <- MQMpermObject(result)			#Create a permutation object
		numresults <- dim(result[[1]])[1]
		lines(x=1:numresults,y=rep(summary(perm_temp)[1,1],numresults),col="green",lwd=2,lty=2)
		lines(x=1:numresults,y=rep(summary(perm_temp)[2,1],numresults),col="blue",lwd=2,lty=2)			
	}else{
		ourstop("Wrong type of result file, please supply a valid MQMmulti object.") 
	}	
}

plot.MQMnice <- function(result = NULL,...){
	#Helperfunction to show MQMmulti objects made by doing multiple scanMQM runs (in a LIST)
	matrix <- NULL
	names <- NULL
	i <- 1
	if(class(result)[2] == "MQMmulti"){		
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
	}else{
		ourstop("Wrong type of result file, please supply a valid MQMmulti object.") 
	}
}

plot.MQMone <- function(result = NULL,result2 = NULL, extended=0,...){
	#Helperfunction to show scanone objects made by doing scanMQM runs
	if(any(class(result) == "scanone")){
		info_c <- result
		info_c[,3]<- info_c[,5]
		if(extended){
			info_l <- result
			info_l[,3] <- result[,4]
			plot(result,info_c,info_l,lwd=1,col=c("black","blue","red"),ylab="QTL (LOD)",...)
			grid(max(result$chr),5)
			labels <- c(colnames(result)[3],colnames(result)[5],colnames(result)[4])
			legend("topright", labels,col=c("black","blue","red"),lty=c(1,1,1))		
		}else{
			if (any(class(result2) == "scanone")){
				#MAX 3 scanone objects
				plot(result,info_c,result2,lwd=1,ylab="QTL (LOD)",...)
				grid(max(result$chr),5)
				labels <- c(colnames(result)[3],colnames(result)[5],colnames(result2)[3])
				legend("topright", labels,col=c("black","blue"),lty=c(1,1))
			}else{
				#MAX 3 scanone objects (here we now have 2)
				plot(result,info_c,lwd=1,ylab="QTL (LOD)",...)
				grid(max(result$chr),5)
				labels <- c(colnames(result)[3],colnames(result)[5])
				legend("topright", labels,col=c("black","blue"),lty=c(1,1))			
			}
		}
	}else{
		ourstop("Wrong type of result file, please supply a valid scanone (from MQM) object.") 
	}
}

