#'  @title plot.scanonevar
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{plot.scanonevar} implements the plot generic for objects of class 'scanonevar'.
#'    Because scanonevar objects can be viewed in terms of LODs or empirical p-values, 
#'    this plotting function checks the 'units' attribute to determine which to plot.
#'  
#'  @param scanonevar the \code{scanonevar} object to be plotted
#'  @param chrs optionally, the subset of the chromosomes to plot
#'  @param col optionally, a vector specifying the colors of the scan lines.  Defaults to \code{c("black", "blue", "red", "darkgreen")}.
#'  @param bandcol optionally, a background color for the even-index chromosomes in this scan.
#'  
#'  @return Returns nothing.  Only makes the plot.
#'    
#'  @details If such a strong signal was observed that the empirical p-value underflows R's 
#'    float type, this function produces an error.  The author is open to suggestions on how
#'    to deal with this situation better.
#'    
#'  @seealso  \code{\link{scanonevar}}, \code{\link{convert.scanonevar.to.empirical.ps}}

plot.scanonevar <- function(scanonevar,
														chrs = unique(scanonevar$chr),
														col = c("black", "blue", "red", "darkgreen"),
														bandcol = 'lightgray',
														legends = c('mean or var', 'mean', 'var'),
														legend.pos = 'top',
														gap = 25,
														incl.markers = TRUE,
														main = attr(scanonevar, 'pheno'),
														ylim = c(0, 1.05*max(coords.y.locus, na.rm = TRUE)),
														scanone.for.comparison = NULL,
														show.equations = length(chrs) != 1)
{

  # validate scanonevar object
  if (!is.scanonevar(scanonevar)) {
    stop(paste('scanonevar argument is not a valid scanonevar object:', attr(is.scanonevar(scanonevar), 'why.not')))
  }
  
	# store current graphical parameters and customize them for this plot
	start.pars <- par(no.readonly = TRUE)
	if (show.equations) {
	  par(mar = c(3,4,5,1))
	} else {
	  par(mar = c(3,4,2,1))
	}

	# convert scanone.for.comparison to tbl_df if needed
	if (!any(is.null(scanone.for.comparison), is.tbl(scanone.for.comparison))) {
	  scanone.for.comparison <- tbl_df(scanone.for.comparison)
	}
	
	# subset scanonevar to necessary chrs only
	if (!identical(chrs, unique(scanonevar$chr))) {
	  
	  temp <- dplyr::filter(scanonevar, chr %in% chrs)
	  scanone.for.comparison <- dplyr::filter(scanone.for.comparison, chr %in% chrs)

		class(temp) <- class(scanonevar)
		attr(temp, 'method') <- attr(scanonevar, 'method')
		attr(temp, 'type') <- attr(scanonevar, 'type')
		attr(temp, 'model') <- attr(scanonevar, 'model')
		attr(temp, 'mean.null.formula') <- attr(scanonevar, 'mean.null.formula')
		attr(temp, 'mean.alt.formula') <- attr(scanonevar, 'mean.alt.formula')
		attr(temp, 'var.null.formula') <- attr(scanonevar, 'var.null.formula')
		attr(temp, 'var.alt.formula') <- attr(scanonevar, 'var.alt.formula')
		attr(temp, 'pheno') <- attr(scanonevar, 'pheno')
		attr(temp, 'null.effects') <- attr(scanonevar, 'null.effects')
		attr(temp, 'units') <- attr(scanonevar, 'units')
		attr(temp, 'null.fit') <- attr(scanonevar, 'null.fit')
		temp$chr <- factor(temp$chr)

		scanonevar <- temp
	}

	# x coordinates for plotting
	scanonevar$chr <- factor(scanonevar$chr)
	levels(scanonevar$chr) <- mixedsort(levels(scanonevar$chr))
	coords.x.chr <- scanonevar %>%
		group_by(chr) %>%
		summarise(len = max(pos) - min(pos)) %>%
		mutate(start = cumsum(dplyr::lag(len + gap, default = 0))) %>%
		mutate(end = start + len) %>%
		mutate(middle = (start + end)/2)

	coords.x.locus <- left_join(coords.x.chr, scanonevar, by = 'chr') %>%
		mutate(coord.x = start + pos)

	# y coordinates for plotting
	if (attr(scanonevar, 'units') == 'lods') {
		coords.y.locus <- scanonevar %>% select(matches('lod'))
		ylab <- 'LOD'
	}
	if (attr(scanonevar, 'units') == 'emp.ps') {
		coords.y.locus <- -log10(select(scanonevar, matches('emp.p')))
		ylab = '-log10(empirical p)'
	}

	# make plotting area
	xlim <- c(-gap/2, max(coords.x.chr$end) + gap/2)
	plot(-42, -42, xlim = xlim, ylim = ylim,
			 type = 'n', xaxs = 'i',
			 xlab = NA, ylab = NA, axes = FALSE)

	# shade in bg for every other chr
	if (!is.null(bandcol) & length(chrs) > 1) {
		names.chrs.even <- chrs[seq(from = 2, to = length(chrs), by = 2)]
		xs.chrs.even <- dplyr::filter(coords.x.chr, chr %in% names.chrs.even)

		rect(xleft = xs.chrs.even$start - gap/2,
				 xright = xs.chrs.even$end + gap/2,
				 ybottom = ylim[1],
				 ytop = ylim[2],
				 border = bandcol, col = bandcol)
	}

	# draw x axis and label chrs
	segments(x0 = xlim[1], x1 = xlim[2],
					 y0 = 0, y1 = 0)
	mtext(side = 1, text = chrs, at = coords.x.chr$middle, line = 0)
	mtext(side = 1, text = 'Chromosome', at = mean(xlim), line = 1)

	# draw y axis and label chrs
	axis(side = 2)
	mtext(side = 2, text = ylab, line = 2, at = mean(ylim))

	# plot lines
	for (test.idx in 1:ncol(coords.y.locus)) {
		test.name <- names(coords.y.locus)[test.idx]
		segments(x0 = coords.x.locus$coord.x,
						 x1 = lead(coords.x.locus$coord.x),
						 y0 = coords.y.locus[[test.name]],
						 y1 = lead(coords.y.locus[[test.name]]),
						 lwd = 2*with(coords.x.locus, pos != len),
						 col = col[test.idx])
	}

	# plot scanone for comparison
	if (attr(scanonevar, 'units') == 'lods') {
	  if (!is.null(scanone.for.comparison)) {
	    segments(x0 = coords.x.locus$coord.x,
	             x1 = lead(coords.x.locus$coord.x),
	             y0 = scanone.for.comparison$lod,
	             y1 = lead(scanone.for.comparison$lod),
	             lwd = 2*with(coords.x.locus, pos != len),
	             col = col[length(col)])
	  }
	}
	if (attr(scanonevar, 'units') == 'emp.ps') {
	  if (!is.null(scanone.for.comparison)) {
	    segments(x0 = coords.x.locus$coord.x,
	             x1 = lead(coords.x.locus$coord.x),
	             y0 = -log10(scanone.for.comparison$emp.p.scanone),
	             y1 = -log10(lead(scanone.for.comparison$emp.p.scanone)),
	             lwd = 2*with(coords.x.locus, pos != len),
	             col = col[length(col)])
	  }
	}
	
	# note: I got rid of thresholds on the LOD scale, because it is visually confusing to
	# look at 3 or 4 different alpha05 thresholds and/or 3 or 4 different alpha 01's  -RWC

	# add thresholds on p-value scale
	if (attr(scanonevar, 'units') == 'emp.ps') {
	  abline(h = -log10(c(0.05, 0.01)), lty = c(1, 2))
	  text(x = coords.x.locus$coord.x[1], y = -log10(0.05), labels = 'alpha=0.05', adj = c(0.5, -0.2))
	}
	
	# plot lines at marker positions (rug plot)
	if (incl.markers) {
		marker.idxs <- !grepl(pattern = 'loc', x = scanonevar$marker.name)
		segments(x0 = coords.x.locus$coord.x[marker.idxs],
						 x1 = coords.x.locus$coord.x[marker.idxs],
						 y0 = 0,
						 y1 = ylim[2]*0.03,
						 col = 'gray50')
	}

	# draw the legend
	if (!is.null(scanone.for.comparison)) { 
	  legends <- c(legends, 'scanone for comparison')
	}
	legend(x = legend.pos, 
	       legend = legends,
				 fill = col, cex = 0.8, bty = 'n',
				 x.intersp = 0.3, y.intersp = 1, xjust = 0.5, yjust = 0)

	# add the title
	if (show.equations) {
	  title <- paste(attr(x = scanonevar, 'pheno'),
	                 '\n', 'mean null:',
	                 paste(as.character(attr(scanonevar, 'mean.null.formula'))[c(2,1,3)], collapse = ' '),
	                 '\n', 'mean alt:',
	                 paste(as.character(attr(scanonevar, 'mean.alt.formula'))[c(2,1,3)], collapse = ' '),
	                 '\n', 'var null:',
	                 paste(as.character(attr(scanonevar, 'var.null.formula')), collapse = ' '),
	                 '\n', 'var alt:',
	                 paste(as.character(attr(scanonevar, 'var.alt.formula')), collapse = ' '))
	  mtext(text = title, side = 3, line = 0)
	} else {
	  mtext(text = attr(scanonevar, 'pheno'), side = 3, line = 0)
	}
	
	# reset graphical parameteers to how they were on start
	par(start.pars)

	# return nothing
	invisible()
}
