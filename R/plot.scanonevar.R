###############################################################################
#
# plot.scanonevar.R
#
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
################################################################################

################################################################################
#
# plot.scanone: plot output from scanone
#
################################################################################

plot.scanonevar <- function(varscan,
														chrs = unique(varscan$chr),
														col = c("black", "blue", "red", "green"),
														bandcol = 'lightgray',
														legend.pos = 'top',
														gap = 25,
														incl.markers = TRUE,
														main = attr(varscan, 'pheno'),
														ylim = c(0, 1.05*max(coords.y.locus, na.rm = TRUE)))
{

	# store current graphical parameters and customize them for this plot
	start.pars <- par(no.readonly = TRUE)
	par(mar = c(3,4,6,1))

	# subset varscan to necessary chrs only
	if (!identical(chrs, unique(varscan$chr))) {
		temp <- dplyr::filter(varscan, chr %in% chrs)

		class(temp) <- c("scanonevar", "tbl_df", "data.frame")
		attr(temp, 'method') <- attr(varscan, 'method')
		attr(temp, 'type') <- attr(varscan, 'type')
		attr(temp, 'model') <- attr(varscan, 'model')
		attr(temp, 'mean.null.formula') <- attr(varscan, 'mean.null.formula')
		attr(temp, 'mean.alt.formula') <- attr(varscan, 'mean.alt.formula')
		attr(temp, 'var.null.formula') <- attr(varscan, 'var.null.formula')
		attr(temp, 'var.alt.formula') <- attr(varscan, 'var.alt.formula')
		attr(temp, 'pheno') <- attr(varscan, 'pheno')
		attr(temp, 'null.effects') <- attr(varscan, 'null.effects')
		attr(temp, 'units') <- attr(varscan, 'units')
		temp$chr <- factor(temp$chr)

		varscan <- temp
	}

	# x coordinates for plotting
	levels(varscan$chr) <- mixedsort(levels(varscan$chr))
	coords.x.chr <- varscan %>%
		group_by(chr) %>%
		summarise(len = max(pos) - min(pos)) %>%
		mutate(start = cumsum(lag(len + gap, default = 0))) %>%
		mutate(end = start + len) %>%
		mutate(middle = (start + end)/2)

	coords.x.locus <- left_join(coords.x.chr, varscan, by = 'chr') %>%
		mutate(coord.x = start + pos)

	# y coordinates for plotting
	if (attr(varscan, 'units') == 'lods') {
		coords.y.locus <- select(varscan, matches('lod'))
		ylab <- 'LOD'
	}
	if (attr(varscan, 'units') == 'emp.ps') {
		coords.y.locus <- -log10(select(varscan, matches('emp.p')))
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

	# plot lines at marker positions (rug plot)
	if (incl.markers) {
		marker.idxs <- !grepl(pattern = 'loc', x = varscan$marker.name)
		segments(x0 = coords.x.locus$coord.x[marker.idxs],
						 x1 = coords.x.locus$coord.x[marker.idxs],
						 y0 = 0,
						 y1 = ylim[2]*0.03,
						 col = 'gray50')
	}

	# draw the legend
	legend(x = legend.pos, legend = names(coords.y.locus),
				 fill = col, cex = 0.8, bty = 'n',
				 x.intersp = 0.3, y.intersp = 0.8, xjust = 0.5, yjust = 0)

	# add the title
	title <- paste(attr(x = varscan, 'pheno'),
	               '\n', 'mean null:',
	               paste(as.character(attr(varscan, 'mean.null.formula'))[c(2,1,3)], collapse = ' '),
	               '\n', 'mean alt:',
	               paste(as.character(attr(varscan, 'mean.alt.formula'))[c(2,1,3)], collapse = ' '),
	               '\n', 'var null:',
	               paste(as.character(attr(varscan, 'var.null.formula')), collapse = ' '),
	               '\n', 'var alt:',
	               paste(as.character(attr(varscan, 'var.alt.formula')), collapse = ' '))
	mtext(text = title, side = 3, line = 0)

	# draw the alpha = 0.05 and 0.01 lines
	# iff we are looking at empirical ps
	if (units(varscan) == 'emp.ps') {
		abline(h = -log10(c(0.05, 0.01)), lty = c(1, 2))
	}

	# reset graphical parameteers to how they were on start
	par(start.pars)

	# return nothing
	invisible()
}
