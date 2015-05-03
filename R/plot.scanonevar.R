#####################################################################
#
# plot.scanonevar.R
#
# copyright (c) 2001-2014, Karl W Broman
# last modified Apr, 2014
# first written Feb, 2001
# modified by Robert Corty in March 2015
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
# Contains: plot.scanonevar,
#
######################################################################

######################################################################
#
# plot.scanone: plot output from scanone
#
######################################################################

plot.scanonevar <- function(x, chr, incl.markers = TRUE, xlim, ylim,
														lty=1, col=c("black", "blue", "red", "green"),
														lwd=2,gap=25,
														mtick=c("line", "triangle"), show.marker.names=FALSE,
														alternate.chrid=FALSE, bandcol=NULL, type="l", cex=1,
														pch=1, bg="transparent", bgrect=NULL, legend.pos = 'top',
														thresholds = NULL, ...)
{

	if(!any(class(x) == "scanonevar")) {
		stop("Input should have class \"scanonevar\".")
	}

	if (!is.factor(x$chr)) {
		x$chr <- factor(x$chr, levels=unique(x$chr))
	}

	dots <- list(...)

	# handle special arguments to be used in lines()
	if(length(col)==1) col <- rep(col, 3)
	if(length(type)==1) type <- rep(type, 3)
	if(length(cex)==1) cex <- rep(cex, 3)
	if(length(pch)==1) pch <- rep(pch, 3)
	if(length(bg)==1) bg <- rep(bg, 3)
	if(length(lty)==1) lty <- rep(lty, 3)
	if(length(lwd)==1) lwd <- rep(lwd, 3)

	mtick <- match.arg(mtick)

	if(length(dim(x))!=2)
		stop("Argument x must be a matrix or data.frame.")

	out <- x

	# pull out desired chromosomes
	if(missing(chr) || length(chr)==0)
		chr <- unique(as.character(out[,1]))
	else chr <- matchchr(chr, unique(out[,1]))

	out <- out[!is.na(match(out[,1],chr)),]

	# beginning and end of chromosomes
	temp <- out
	begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
	rownames(begend) <- unique(out[,1])
	begend <- begend[as.character(chr),,drop=FALSE]
	len <- begend[,2]-begend[,1]

	# locations to plot start of each chromosome
	start <- c(0,cumsum(len+gap))-c(begend[,1],0)

	maxx <- sum(len+gap)-gap

	# graphics parameters
	old.xpd <- par("xpd")
	old.las <- par("las")
	par(xpd=FALSE,las=1)
	on.exit(par(xpd=old.xpd,las=old.las))

	ys <- data.frame(out[,grepl(pattern = 'lod', x = names(out))])
	if(length(thresholds)==3) {
		ys <- ys - rep(thresholds, each = nrow(ys))
		thresholds <- rep(1, 3)
	}
	if(length(thresholds)==6) {
		ys <- ys - rep(thresholds[c(1,3,5)], each = nrow(ys))
		thresholds <- c(1, thresholds[2] - thresholds[1], 1, thresholds[4] - thresholds[3], 1, thresholds[6] - thresholds[5])
	}

	# make frame of plot
	if (missing(ylim)) {
		if ('lod' %in% names(ys)) {	ylim.set <- c(min(ys, thresholds), 1.05*max(ys, thresholds))	} else { ylim.set <- range(ys) }
	} else {
		ylim.set <- ylim
	}
	if(missing(xlim)) {
		xlim <- c(-gap/2,maxx+gap/2)
	}

	plot(0,0,ylim=ylim.set,xlim=xlim,type="n",xaxt="n",
			 xlab="Chromosome",ylab=NA, xaxs="i",	 ...)

	if(!is.null(bgrect)) {
		u <- par("usr")
		rect(u[1], u[3], u[2], u[4], col=bgrect)
	}
	if(!is.null(bandcol)) {
		u <- par("usr")
		for(i in seq(2, by=2, length(chr))) {
			rect(min(out[out[,1]==chr[i],2]) + start[i]-gap/2, u[3],
					 max(out[out[,1]==chr[i],2]) + start[i]+gap/2, u[4], border=bandcol, col=bandcol)
		}
		abline(h=u[3:4])
	}

	# initialize xtick and xtickmark
	xtick <- NULL
	xticklabel <- NULL
	for(i in 1:length(chr)) {
		# plot first out
		x <- out[out[,1]==chr[i], 2]+start[i]

		for(j in 1:ncol(ys)) {

			y <- ys[out[,1]==chr[i], j]

			lines(x,y,lwd=lwd[j],lty=lty[j],col=col[j], type=type[j], cex=cex[j], pch=pch[j], bg=bg[j])

		}

		# plot chromosome number
		tloc <- mean(c(min(x),max(x)))
		xtick <- c(xtick, tloc)
		xticklabel <- c(xticklabel, as.character(chr[i]))

		# plot lines or triangles at marker positions
		nam <- dimnames(out)[[1]][out[,1]==chr[i]]
		wh.genoprob <- grep("^c.+\\.loc-*[0-9]+", nam)
		if(length(wh.genoprob)==0) wh.genoprob <- seq(along=nam)
		else wh.genoprob <- (seq(along=nam))[-wh.genoprob]
		pos <- out[out[,1]==chr[i],2][wh.genoprob]+start[i]

		if(incl.markers) {
			if(mtick=="line")
				rug(pos, 0.02, quiet=TRUE)
			else {
				a <- par("usr")
				points(pos, rep(a[3]+diff(a[3:4])*0.01, length(pos)), pch=17, cex=1.5)
			}
		}

		if(show.marker.names) {
			a <- par("usr")
			text(pos, rep(a[3]+diff(a[3:4])*0.03, length(pos)), nam[wh.genoprob],
					 srt=90, adj=c(0,0.5))
		}

	}

	# draw the thresholds
	if (length(thresholds) == 3) {
		abline(h = 1)
	}
	if (length(thresholds) == 6) {
		abline(h = c(0, thresholds[2], thresholds[4], thresholds[6]),
					 col = c('black', 'black', 'blue', 'red'),
					 lty = c(1, 2, 2, 2))
	}

	# draw the legend
	legend(x = legend.pos, legend = names(ys),
				 fill = col, cex = 0.8, bty = 'n',
				 x.intersp = 0.3, y.intersp = 0.8, xjust = 0.5, yjust = 0)

	# draw the axis
	if(!alternate.chrid || length(xtick) < 2) {
		for(i in seq(along=xtick))
			axis(side=1, at=xtick[i], labels=xticklabel[i])
	}
	else {
		odd <- seq(1, length(xtick), by=2)
		even <- seq(2, length(xtick), by=2)
		for(i in odd) {
			axis(side=1, at=xtick[i], labels="")
			axis(side=1, at=xtick[i], labels=xticklabel[i], line=-0.4, tick=FALSE)
		}
		for(i in even) {
			axis(side=1, at=xtick[i], labels="")
			axis(side=1, at=xtick[i], labels=xticklabel[i], line=+0.4, tick=FALSE)
		}
	}

	if(!is.null(bgrect)) {
		u <- par("usr")
		rect(u[1], u[3], u[2], u[4])
	}

	mtext(text = attr(x = out, 'pheno'), side = 3, outer = TRUE, line = 1)

	invisible()
}

# end of plot.scanone.R
