######################################################################
#
# fitplot.scanonevar.R
#
#
# Part of the R/qtl package
#
######################################################################

fitplot.scanonevar <- function (cross, marker.name, var.scan) {

	par(mfrow = c(1,1), mar = c(4, 4, 1, 1))

	phen.name <- attr(x = var.scan, which = 'pheno')
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
			 ylim = range(phenotypes)*c(0.9, 1.1),
			 col = 'gray')
	axis(side = 1, at = 1:3, labels = c('AA', 'AB', 'BB'), line = 0)
	mtext(text = 'Genotype', side = 1, line = 2.2, cex = 1.2)
	mtext(text = phen.name, side = 2, line = 2.2, cex = 1.2, las = 0)


	# put I-bar based on lm fit
	g <- (genotypes - 1)/2
	fit <- lm(phenotypes ~ g)
	baseline <- coef(fit)[1]
	add.effect <- coef(fit)[2]
	std <- summary(fit)$sigma

	aa.bot <- baseline - std
	ab.bot <- baseline + add.effect/2 - std
	bb.bot <- baseline + add.effect - std

	aa.top <- baseline + std
	ab.top <- baseline + add.effect/2 + std
	bb.top <- baseline + add.effect + std

	segments(x0 = (1:3)-0.2,
					 y0 = c(aa.bot, ab.bot, bb.bot),
					 x1 = (1:3)-0.2,
					 y1 = c(aa.top, ab.top, bb.top),
					 lwd = 3,
					 col = 'darkred')
	segments(x0 = rep((1:3)-0.25, each = 2),
					 x1 = rep((1:3)-0.15, each = 2),
					 y0 = c(aa.bot, aa.top, ab.bot, ab.top, bb.bot, bb.top),
					 y1 = c(aa.bot, aa.top, ab.bot, ab.top, bb.bot, bb.top),
					 lwd = 3,
					 col = 'darkred')



	# put I-bar on each sample based on sample mean and sd
	aa.mean <- mean(phenotypes[genotypes == 1], na.rm = TRUE)
	ab.mean <- mean(phenotypes[genotypes == 2], na.rm = TRUE)
	bb.mean <- mean(phenotypes[genotypes == 3], na.rm = TRUE)
	aa.sd <- sd(phenotypes[genotypes == 1], na.rm = TRUE)
	ab.sd <- sd(phenotypes[genotypes == 2], na.rm = TRUE)
	bb.sd <- sd(phenotypes[genotypes == 3], na.rm = TRUE)

	aa.top <- aa.mean + aa.sd
	aa.bot <- aa.mean - aa.sd
	ab.top <- ab.mean + ab.sd
	ab.bot <- ab.mean - ab.sd
	bb.top <- bb.mean + bb.sd
	bb.bot <- bb.mean - bb.sd

	segments(x0 = 1:3,
					 y0 = c(aa.bot, ab.bot, bb.bot),
					 x1 = 1:3,
					 y1 = c(aa.top, ab.top, bb.top),
					 lwd = 3)
	segments(x0 = rep((1:3)-0.05, each = 2),
					 x1 = rep((1:3)+0.05, each = 2),
					 y0 = c(aa.bot, aa.top, ab.bot, ab.top, bb.bot, bb.top),
					 y1 = c(aa.bot, aa.top, ab.bot, ab.top, bb.bot, bb.top),
					 lwd = 3)


	# put I-bar to the right of each sample based on scanonevar fit
	params <- var.scan[marker.name,]

	aa.mean <- params$mean.baseline - params$mean_add_effect
	ab.mean <- params$mean.baseline
	bb.mean <- params$mean.baseline + params$mean_add_effect

	aa.sd <- params$disp.baseline - params$disp_add_effect
	ab.sd <- params$disp.baseline
	bb.sd <- params$disp.baseline + params$disp_add_effect

	if (attr(var.scan, 'dom')) {

		ab.mean <- ab.mean + params$mean_dom_effect

		ab.sd <- ab.sd + params$disp_dom_effect
	}

	aa.sd <- sqrt(exp(aa.sd))
	ab.sd <- sqrt(exp(ab.sd))
	bb.sd <- sqrt(exp(bb.sd))

	means <- c(aa.mean, ab.mean, bb.mean)
	sds <- c(aa.sd, ab.sd, bb.sd)

	segments(x0 = (1:3)+0.2,
					 y0 = means - sds,
					 x1 = (1:3)+0.2,
					 y1 = means + sds,
					 lwd = 3,
					 col = 'darkgreen')
	segments(x0 = rep((1:3)+0.15, each = 2),
					 x1 = rep((1:3)+0.25, each = 2),
					 y0 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
					 y1 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
					 lwd = 3,
					 col = 'darkgreen')

}
