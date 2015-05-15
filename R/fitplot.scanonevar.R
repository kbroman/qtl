#####################################################################
#
# fitplot.scanonevar.R
#
# copyright (c) 2001-2014, Karl W Broman
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
# Contains: fitplot.scanonevar,
#
######################################################################

######################################################################
#
# fitplot.scanonevar:
# Plot to show the fitted model that results from scanonevar() at
# one specific locus of interest.
# Useful to illustrate what data pattern drives a significant LOD score.
#
######################################################################

fitplot.scanonevar <- function (cross, marker.name, scanonevar, ...) {

# 	par(mfrow = c(1,1), mar = c(4, 4, 1, 1))
	set.seed(27599)

	phen.name <- attr(x = scanonevar, which = 'pheno')
	var.scan <- scanonevar$scan
	phenotypes <- cross$pheno[[phen.name]]

	chr.with.marker <- which(sapply(X = cross$geno,
																	FUN = function(x, mn) { mn %in% colnames(x$data)},
																	marker.name))
	if (length(chr.with.marker) == 0) {	stop('Marker not found.')	}
	if (length(chr.with.marker) >= 2) {	stop('Marker found 2+ times.')	}

	genotypes <- cross$geno[[chr.with.marker]]$data[,marker.name]

	# plot phenotypes, broken out by genotype at this marker
	plot(x = genotypes + runif(n = nind(cross), min = -0.2, max = 0.2),
			 y = phenotypes,
			 xlab = NA,
			 ylab = NA,
			 xaxt = 'n',
			 xlim = c(0.5, 3.5),
			 ylim = range(phenotypes)*c(0.95, 1.05),
			 col = 'gray',
			 ...)
	axis(side = 1, at = 1:3, labels = c('AA', 'AB', 'BB'), line = 0)
	mtext(text = 'Genotype', side = 1, line = 2.2, cex = 1.2)
	mtext(text = phen.name, side = 2, line = 2.2, cex = 1.2, las = 0)


	# put I-bar to the right of each sample based on scanonevar fit
	params <- var.scan[marker.name,]

	aa.mean <- params$mean.baseline - params$mean_add_effect
	ab.mean <- params$mean.baseline
	bb.mean <- params$mean.baseline + params$mean_add_effect

	aa.sd <- params$disp.baseline - params$disp_add_effect
	ab.sd <- params$disp.baseline
	bb.sd <- params$disp.baseline + params$disp_add_effect

	if (attr(scanonevar, 'dom')) {

		ab.mean <- ab.mean + params$mean_dom_effect
		ab.sd <- ab.sd + params$disp_dom_effect
	}

	aa.sd <- sqrt(exp(aa.sd))
	ab.sd <- sqrt(exp(ab.sd))
	bb.sd <- sqrt(exp(bb.sd))

	means <- c(aa.mean, ab.mean, bb.mean)
	sds <- c(aa.sd, ab.sd, bb.sd)

	# vertical lines
	segments(x0 = (1:3),
					 y0 = means - sds,
					 x1 = (1:3),
					 y1 = means + sds,
					 lwd = 3,
					 col = 'black')

	# fitten mean lines
	segments(x0 = (1:3)-0.05,
					 x1 = (1:3)+0.05,
					 y0 = means,
					 y1 = means,
					 lwd = 3,
					 col = 'black')

	# 1 SD away from mean lines
	segments(x0 = rep((1:3)-0.05, each = 2),
					 x1 = rep((1:3)+0.05, each = 2),
					 y0 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
					 y1 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
					 lwd = 3,
					 col = 'black')

	# horizontalish connectors
	segments(x0 = rep(1:2, 3),
					 x1 = rep(2:3, 3),
					 y0 = c(means[-3]+sds[-3], means[-3], means[-3]-sds[-3]),
					 y1 = c(means[-1]+sds[-1], means[-1], means[-1]-sds[-1]),
					 col = c(rep('red', 2), rep('blue', 2), rep('red', 2)))

}
