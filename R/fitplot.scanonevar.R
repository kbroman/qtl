#####################################################################
#
# fitplot.scanonevar.R
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

fitplot.scanonevar <- function(cross, 
                               name.of.marker.to.plot,
                               varscan,
                               ...) {

	# store current graphical parameters and customize them for this plot
	start.pars <- par(no.readonly = TRUE)
	par(mar = c(4,4,4,2))
	set.seed(27599)

	# pull out phenotype values
	phen.name <- attr(x = varscan, which = 'pheno')
	phen <- cross$pheno[[phen.name]]
	if (is.null(phen)) { stop('Phenotype of scan not found in cross')}

	# pull out genotypes at the marker to plot
	orig.name <- strsplit(x = name.of.marker.to.plot, split = '_', fixed = TRUE)[[1]][-1]
	genotypes <- GetGenoprobsByMarkerName(cross = cross,
																				marker.name = orig.name)
  genotypes <- apply(X = genotypes, MARGIN = 1, FUN = which.max)
  
	# check that varscan included chr where marker is
	mles <- dplyr::filter(varscan, marker.name == name.of.marker.to.plot)
	if (nrow(mles) == 0) { stop('Desired marker is not in scan.')}

	# find all the peaks on this chromosome
	# name.of.chr.with.marker <- mles$chr
# 	chr.peaks <- GetPeaksFromVarscan(vs = dplyr::filter(varscan, chr == name.of.chr.with.marker))
# 
# 	# distance of desired marker to each peak
# 	full.peaks <- dplyr::filter(chr.peaks, full.peak == TRUE)
# 	dists <- abs(mles$pos - (full.peaks)$pos)
# 	best.full.locus <- slice(full.peaks, which.min(dists))
# 
# 	mean.peaks <- dplyr::filter(chr.peaks, mean.peak == TRUE)
# 	dists <- abs(mles$pos - (mean.peaks)$pos)
# 	best.mean.locus <- slice(mean.peaks, which.min(dists))
# 
# 	var.peaks <- dplyr::filter(chr.peaks, var.peak == TRUE)
# 	dists <- abs(mles$pos - (var.peaks)$pos)
# 	best.var.locus <- slice(var.peaks, which.min(dists))


	# plot phenotypes, broken out by genotype at this marker
	plot(x = genotypes + runif(n = nind(cross), min = -0.2, max = 0.2),
			 y = phen,
			 xlab = NA,
			 ylab = NA,
			 xaxt = 'n',
			 xlim = c(0.5, 3.5),
			 ylim = range(phen)*c(0.95, 1.05),
			 col = 'black', 
			 main = 'Chromosome 2, marker at 48.3 cM',
			 ...)
	mtext(text = 'Genotype', side = 1, line = 2.2, cex = 1.2)
	mtext(text = phen.name, side = 2, line = 2.2, cex = 1.2, las = 0)

# 	if (units(varscan) == 'lods') {
# 		title(paste('Full Model: LOD here',
# 								round(mles$lod.full, 2),
# 								'LOD at closest full peak:',
# 								round(best.full.locus$lod.full, 2),
# 								'\nMean Model: LOD here',
# 								round(mles$lod.mean, 2),
# 								'LOD at closest mean peak:',
# 								round(best.mean.locus$lod.mean, 2),
# 								'\nVar Model: LOD here',
# 								round(mles$lod.var, 2),
# 								'LOD at closest var peak:',
# 								round(best.var.locus$lod.var, 2)))
# 	}
# 	if (units(varscan) == 'emp.ps') {
# 		title(paste('Full Model: -log10p here',
# 								round(-log10(mles$emp.p.full), 2),
# 								'-log10p at closest full peak:',
# 								round(-log10(best.full.locus$emp.p.full), 2),
# 								'\nMean Model: -log10p here',
# 								round(-log10(mles$emp.p.mean), 2),
# 								'-log10p at closest mean peak:',
# 								round(-log10(best.mean.locus$emp.p.mean), 2),
# 								'\nVar Model: -log10p here',
# 								round(-log10(mles$emp.p.var), 2),
# 								'-log10p at closest var peak:',
# 								round(-log10(best.var.locus$emp.p.var), 2)))
# 	}

	# autosome-specific stuff
	# if (mles$chrtype == 'A') {

		axis(side = 1, at = 1:3,
				 labels = c('AA', 'AB', 'BB'),
				 line = 0)

# 		# compute predictions for each group based on
# 		# parameters fitted in scanonevar
# 		aa.mean <- mles$mean.baseline 
# 		ab.mean <- mles$mean.baseline
# 		bb.mean <- mles$mean.baseline
# 
# 		if ('mean.add' %in% names(mles)) {
# 		  aa.mean <- aa.mean - mles$mean.add
# 		  bb.mean <- bb.mean + mles$mean.add
# 		}
# 		
# 		if ('mean.dom' %in% names(mles)) {
# 		  ab.mean <- ab.mean + mles$mean.dom
# 		}
# 		
# 		aa.var <- mles$var.baseline 
# 		ab.var <- mles$var.baseline
# 		bb.var <- mles$var.baseline
# 		
# 		if ('var.add' %in% names(mles)) {
# 		  aa.var <- aa.var - mles$var.add
# 		  bb.var <- bb.var + mles$var.add
# 		}
# 		
# 		if ('var.dom' %in% names(mles)) {
# 		  ab.var <- ab.var + mles$var.dom
# 		}
	# }

# 	# X chromosome-specific stuff
# 	if (mles$chrtype == 'X') {
# 
# 		axis(side = 1, at = 1:3,
# 				 labels = c('A- or AA', 'AB', 'B- or BB'),
# 				 line = 0)
# 
# 		# compute predictions for each group based on
# 		# parameters fitted in scanonevar (additive only here)
# 	  aa.mean <- mles$mean.baseline 
# 	  ab.mean <- mles$mean.baseline
# 	  bb.mean <- mles$mean.baseline
# 	  
# 	  if ('mean.add' %in% names(mles)) {
# 	    ab.mean <- ab.mean + 0.5*mles$mean.add.effect
# 	    bb.mean <- bb.mean + mles$mean.add.effect
# 	  }
# 	  
# 	  aa.var <- mles$var.baseline 
# 	  ab.var <- mles$var.baseline
# 	  bb.var <- mles$var.baseline
# 	  
# 	  if ('var.add' %in% names(mles)) {
# 	    ab.var <- ab.var + 0.5*mles$var.add.effect
# 	    bb.var <- bb.var + mles$var.add.effect
# 	  }
# 	}


# 	aa.sd <- sqrt(exp(aa.var))
# 	ab.sd <- sqrt(exp(ab.var))
# 	bb.sd <- sqrt(exp(bb.var))
# 
# 	means <- c(aa.mean, ab.mean, bb.mean)
# 	sds <- c(aa.sd, ab.sd, bb.sd)
# 
# 	# vertical lines
# 	segments(x0 = (1:3),
# 					 y0 = means - sds,
# 					 x1 = (1:3),
# 					 y1 = means + sds,
# 					 lwd = 3,
# 					 col = 'black')
# 
# 	# fitted mean lines
# 	segments(x0 = (1:3) - 0.05,
# 					 x1 = (1:3) + 0.05,
# 					 y0 = means,
# 					 y1 = means,
# 					 lwd = 3,
# 					 col = 'black')
# 
# 	# 1 SD away from mean lines
# 	segments(x0 = rep((1:3) - 0.05, each = 2),
# 					 x1 = rep((1:3) + 0.05, each = 2),
# 					 y0 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
# 					 y1 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
# 					 lwd = 3,
# 					 col = 'black')
# 
# 	# horizontalish connectors
# 	segments(x0 = rep(1:2, 3),
# 					 x1 = rep(2:3, 3),
# 					 y0 = c(means[-3] + sds[-3], means[-3], means[-3] - sds[-3]),
# 					 y1 = c(means[-1] + sds[-1], means[-1], means[-1] - sds[-1]),
# 					 col = c(rep('red', 2), rep('blue', 2), rep('red', 2)))

	par(start.pars)
}
