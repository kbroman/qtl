####################################################################
#
# scanonevar_.R
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
# Contains: scanonevar_,
#
######################################################################

######################################################################
#
# scanonevar_:
# workhorse function called by user-facing function scanonevar()
#
######################################################################


scanonevar_ <- function(cross,
												chrs = 1:length(cross$geno),
												dom,
												X,
												quiet = TRUE,
												mean.null.formula,
												var.null.formula,
												mean.alt.formula,
												var.alt.formula,
												return.effects = TRUE,
												mean.perm = seq(1, nind(cross)),
												var.perm = seq(1, nind(cross)),
												calc.lod.full = TRUE,
												calc.lod.mean = TRUE,
												calc.lod.var = TRUE) {

	crosstype <- class(cross)[1]

	# in case user specifies mean.perm or var.perm as NULL
	# same effect as omitting them, but a bit more explicit in the call 
	# as to what will be done
	if (is.null(mean.perm)) { mean.perm <- seq(1, nind(cross)) }
	if (is.null(var.perm)) { var.perm <- seq(1, nind(cross)) }

	# fit covariate-only model
	# covariates have mean and varianc effect, genetic markers have none
	# only need to do this once per scan
	d.fit.null <- dglm(formula = mean.null.formula,
										 dformula = var.null.formula,
										 data = X)
	ln.lik.null <- -0.5*d.fit.null$m2loglik
	log10.lik.null <- ln.lik.null / log(10)
	null.effects <- list(meancovs = d.fit.null$coef[-1],
											 varcovs = d.fit.null$disp$coef[-1])

	# scan by looping over the chromosomes
	scan <- NULL
	for(chr.idx in chrs) {

		if(!quiet) message(" - Chr ", chrnames(cross)[chr.idx])
		this.chr <- cross$geno[[chr.idx]]

		# error out if chr is not autosome and not X
		if (!(class(this.chr)[1] %in% c('A', 'X'))) {
			stop(paste('scanonevar not implemented for chr',
								 chr.idx,
								 'of type',
								 class(this.chr)))
		}

		# set up genotype design matrix for autosome
		if (class(this.chr)[1] == 'A') {
			if (crosstype=="f2") {
				g11 <- this.chr$prob[,,1]
				g12 <- this.chr$prob[,,2]
				g13 <- this.chr$prob[,,3]
				a1  <- -g11 + g13
				d1 <-  g12
				chr.dom <- dom
			}
			else {
				a1 <- cross$geno[[chr.idx]]$prob[,,2]
			}
		}

		# set up genotype design matrix for X chromosome
		if (class(this.chr)[1] == 'X') {
			chr.dom <- FALSE
			if (crosstype=="f2") {
				a1 <- this.chr$prob[,,2]
			}
			else {
				print("scanonevar not implemented for X chromosome of non-F2 crosses")
				next
			}
		}


		# set up containers for this chr
		n.loci <- ncol(a1)
		lod.full <- lod.mean <- lod.var <- numeric(n.loci)
		mean.baseline <- var.baseline <- numeric(n.loci)
		mean.add.effect <- var.add.effect <- numeric(n.loci)
		if (chr.dom) { mean.dom.effect <- var.dom.effect <- numeric(n.loci) }

		for(i in 1:n.loci) { # loop over positions within chromosome

			# fill in genotype probs for this locus
			# mean.perm and var.perm are just 1:N by default so don't actually do any permuting
			X$mean.add.genet <- a1[mean.perm, i]
			X$var.add.genet <- a1[var.perm, i]
			if (chr.dom) {
				X$mean.dom.genet <- d1[mean.perm, i]
				X$var.dom.genet <- d1[var.perm, i]
			}

			# fit full model
			d.fit.full <- dglm(formula = mean.alt.formula,
												 dformula = var.alt.formula,
												 data = X)
			ln.lik.full <- -0.5*d.fit.full$m2loglik
			log10.lik.full <- ln.lik.full / log(10)

			# fit variance-only model if needed (for lod.mean)
			# covariates, but not the genetic marker, have mean effects
			# both covariates and genetic marker have variance effects
			if (calc.lod.mean) {
				d.fit.nomean <- dglm(formula = mean.null.formula,
														 dformula = var.alt.formula,
														 data = X)
				ln.lik.nomean <- -0.5*d.fit.nomean$m2loglik
				log10.lik.nomean <- ln.lik.nomean / log(10)
			}

			# fit mean-only model if needed (for lod.var)
			# covariates, but not the genetic marker, have variance effects
			# both covariates and genetic marker have mean effects
			if (calc.lod.var) {
				d.fit.novar <- dglm(formula = mean.alt.formula,
														dformula = var.null.formula,
														data = X)
				ln.lik.novar <- -0.5*d.fit.novar$m2loglik
				log10.lik.novar <- ln.lik.novar / log(10)
			}


			# save LOD's and fitted effects to chr-level container
			if (calc.lod.full) { lod.full[i] <- log10.lik.full - log10.lik.null }
			if (calc.lod.mean) { lod.mean[i] <- log10.lik.full - log10.lik.nomean }
			if (calc.lod.var) { lod.var[i] <- log10.lik.full - log10.lik.novar }

			if (return.effects) {

				mean.baseline[i] <- coef(d.fit.full)[1]
				var.baseline[i] <- coef(d.fit.full$dispersion.fit)[1]
				mean.add.effect[i] <- coef(d.fit.full)[2]
				var.add.effect[i] <- coef(d.fit.full$dispersion.fit)[2]

				if (chr.dom) {
					mean.dom.effect[i] <- coef(d.fit.full)[3]
					var.dom.effect[i] <- coef(d.fit.full$dispersion.fit)[3]
				}
			}

		}

		# set up the output
		map <- attr(cross$geno[[chr.idx]]$prob, "map")
		this.chr.scan <- tbl_df(data.frame(chrtype = class(this.chr)[1],
																			 chr = chrnames(cross)[chr.idx],
																			 pos = as.vector(map),
																			 marker.name = names(map),
																			 stringsAsFactors = FALSE,
																			 row.names = NULL))

		if (calc.lod.full) { this.chr.scan <- cbind(this.chr.scan, lod.full) }
		if (calc.lod.mean) { this.chr.scan <- cbind(this.chr.scan, lod.mean) }
		if (calc.lod.var) { this.chr.scan <- cbind(this.chr.scan, lod.var) }

		if (return.effects) {

			this.chr.scan <- cbind(this.chr.scan,
														 data.frame(mean.baseline = mean.baseline,
														 					 var.baseline = var.baseline,
														 					 mean.add.effect = mean.add.effect,
														 					 var.add.effect = var.add.effect))

			# if dom is used in this chr, save those effects to output
			# if dom is used in this scan but not in this chr, put NA so
			# that the df of this chr conforms with the others
			if (dom) {
				if (chr.dom) {
					this.chr.scan$mean.dom.effect <- mean.dom.effect
					this.chr.scan$var.dom.effect <- var.dom.effect
				} else {
					this.chr.scan$mean.dom.effect <- this.chr.scan$var.dom.effect <- rep(NA, n.loci)
				}
			}

		}

		if(is.null(scan)) {
			scan <- this.chr.scan
		}	else {
			scan <- rbind(scan, this.chr.scan)
		}
	}

	# formalize output
	scan$chr <- factor(scan$chr)

	class(scan) <- c("scanonevar", "tbl_df", "data.frame")
	attr(scan, 'method') <- 'dglm package'
	attr(scan, 'type') <- crosstype
	attr(scan, 'model') <- 'normal dglm'
	attr(scan, 'dom') <- dom
	attr(scan, 'pheno') <- names(X)[1]
	attr(scan, 'null.effects') <- null.effects
	attr(scan, 'units') <- 'lods'

	return(scan)
}
