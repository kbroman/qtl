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
												 lty=1, col=c("black","blue","red"), lwd=2,gap=25,
												 mtick=c("line", "triangle"), show.marker.names=FALSE,
												 alternate.chrid=FALSE, bandcol=NULL, type="l", cex=1,
												 pch=1, bg="transparent", bgrect=NULL, legend.pos = 'top', ...)
{

	if(!any(class(x) == "scanone")) {
		stop("Input should have class \"scanone\".")
	}

	if (!is.factor(x$chr)) {
		x$chr <- factor(x$chr, levels=unique(x$chr))
	}

	dots <- list(...)

	num.plots <- 2
	if (attr(x, 'dom') == TRUE) { num.plots <- 4 }

	# handle special arguments to be used in lines()
	if(length(type)==1) type <- rep(type, 2)
	if(length(cex)==1) cex <- rep(cex, 2)
	if(length(pch)==1) pch <- rep(pch, 2)
	if(length(bg)==1) bg <- rep(bg, 2)
	if(length(lty)==1) lty <- rep(lty, 2)
	if(length(lwd)==1) lwd <- rep(lwd, 2)

	if(length(legend.pos) == 1) { legend.pos <- rep(legend.pos, num.plots)}

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

	par(mfrow=c(2,num.plots/2),
			oma = c(1,1,3,1),
			mar = c(2,3,0,1))

	for (plot.num in 1:num.plots) {

		y1.idx <- (plot.num-1)*2+5
		y2.idx <- (plot.num-1)*2+6
		y1 <- out[, y1.idx]
		y2 <- out[, y2.idx]

		# make frame of plot
		if (missing(ylim)) {
			if (grepl(pattern = 'neglogP', x = names(out)[y1.idx])) {
				ylim.set <- c(0, max(y1, y2))
			} else if (grepl(pattern = 'effect', x =  names(out)[y1.idx])) {
				ylim.set <- range(c(y1, y2))
			} else {
				stop("Input object must have columns 5-12 as neglogPs and effects")
			}
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
			y1 <- out[out[,1]==chr[i], y1.idx]
			lines(x,y1,lwd=lwd[1],lty=lty[1],col=col[1], type=type[1], cex=cex[1], pch=pch[1], bg=bg[1])
			y2 <- out[out[,1]==chr[i], y2.idx]
			lines(x,y2,lwd=lwd[2],lty=lty[2],col=col[2], type=type[2], cex=cex[2], pch=pch[2], bg=bg[2])

			legend(x = legend.pos[plot.num], legend = names(out)[c(y1.idx, y2.idx)],
						 fill = col[1:2], cex = 0.8, bty = 'n',
						 x.intersp = 0.3, y.intersp = 0.8, xjust = 0.5, yjust = 0)

			# plot chromosome number
			if (plot.num > num.plots /2) {
				tloc <- mean(c(min(x),max(x)))
				xtick <- c(xtick, tloc)
				xticklabel <- c(xticklabel, as.character(chr[i]))
			}


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
	}

	mtext(text = attr(x = out, 'pheno'), side = 3, outer = TRUE, line = 1)

	invisible()
}

# end of plot.scanone.R
