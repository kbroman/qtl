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
												X,
												mean.null.formula,
												var.null.formula,
												mean.alt.formula,
												var.alt.formula,
												quiet = TRUE,
												mean.add = TRUE,
												mean.dom = TRUE,
												var.add = TRUE,
												var.dom = TRUE,
												return.effects = FALSE,
												mean.perm = seq(1, nind(cross)),
												var.perm = seq(1, nind(cross)),
												calc.lod.full = (mean.add | mean.dom | var.add | var.dom) ,
												calc.lod.mean = (mean.add | mean.dom) ,
												calc.lod.var = (var.add | var.dom) ) {

	crosstype <- class(cross)[1]

	# fit covariate-only model
	# covariates have mean and varianc effect, genetic markers have none
	# only need to do this once per scan
	# and only if lod.full desired
	if (calc.lod.full) {
	  d.fit.null <- dglm(formula = mean.null.formula,
	                     dformula = var.null.formula,
	                     data = X)
	  ln.lik.null <- -0.5*d.fit.null$m2loglik
	  log10.lik.null <- ln.lik.null / log(10)
	  null.effects <- list(meancovs = d.fit.null$coef[-1],
	                       varcovs = d.fit.null$disp$coef[-1])
	}
	
	# loop over chromosomes
	scan <- NULL
	for (chr.idx in chrs) {

		if (!quiet) message(" - Chr ", chrnames(cross)[chr.idx])
		this.chr <- cross$geno[[chr.idx]]
    this.chr.type <- class(this.chr)[1]
    
		# error if chr is not autosome and not X
		if (!(class(this.chr)[1] %in% c('A', 'X'))) {
			stop(paste('scanonevar not implemented for chr',
								 chr.idx,
								 'of type',
								 class(this.chr)))
		}

		# pull out genotypes for this chr if it's an autosome
		if (this.chr.type == 'A') {
			if (crosstype == "f2") {
				g11 <- this.chr$prob[,,1]
				g12 <- this.chr$prob[,,2]
				g13 <- this.chr$prob[,,3]
				add.design  <- -g11 + g13
				if (or(mean.dom, var.dom)) { dom.design <-  g12 }
				this.chr.mean.add <- mean.add
				this.chr.mean.dom <- mean.dom
				this.chr.var.add <- var.add
				this.chr.var.dom <- var.dom
			}
			else {
				add.design <- cross$geno[[chr.idx]]$prob[,,2]
			}
		}

		# pull out genotypes for this chr if it's X chr
		# note that with read.cross(..., convertX = FALSE),
		# prob[,,1] = P(A) and prob[,,2] = P(B)
		if (this.chr.type == 'X') {
			if (crosstype == "f2") {
				add.design <- this.chr$prob[,,2]
				this.chr.mean.add <- mean.add
				this.chr.mean.dom <- mean.dom
				this.chr.mean.dom <- this.chr.var.dom <- FALSE
			}
			else {
				print("scanonevar not implemented for X chromosome of non-F2 crosses")
				next
			}
		}

		# set up containers for this chromosome
		# todo each covar needs to be separate for mean and var
		n.loci <- ncol(add.design)
		
		# 4 cols for metadata: chrtype, chr, pos, marker.name
		map <- attr(cross$geno[[chr.idx]]$prob, "map")
		result.this.chr <- data.frame(chrtype = class(this.chr)[1],
		                              chr = chrnames(cross)[chr.idx],
		                              pos = as.vector(map),
		                              marker.name = names(map),
		                              stringsAsFactors = FALSE,
		                              row.names = NULL)
		
		# up to 3 cols for lod scores
		if (calc.lod.full) { result.this.chr$lod.full <- NA }
		if (calc.lod.mean) { result.this.chr$lod.mean <- NA }
		if (calc.lod.var) { result.this.chr$lod.var <- NA }
		
		# cols for genetic and covariate effects added later...

		# loop over loci in this chromosome
		for (i in 1:n.loci) {

		  # remove genetic effects where they are identical to covariates
		  locus.name <- colnames(add.design)[i]
		  if (locus.name %in% all.vars(mean.alt.formula)) {
		    this.chr.mean.add <- FALSE
		    this.chr.mean.dom <- FALSE
		  }
		  if (locus.name %in% all.vars(var.alt.formula)) {
		    this.chr.var.add <- FALSE
		    this.chr.var.dom <- FALSE
		  }
			# fill in genotype probs for this locus
			# mean.perm and var.perm are 1:N by default
		  # setting covars to all zero makes dglm() to ignore them
		  X$mean.add <- X$var.add <- 0
		  X$mean.dom <- X$var.dom <- 0
		  if (this.chr.mean.add) { X$mean.add <- add.design[mean.perm, i] }
			if (this.chr.var.add) { X$var.add <- add.design[var.perm, i] }
			if (this.chr.mean.dom) { X$mean.dom <- dom.design[mean.perm, i] }
			if (this.chr.var.dom) { X$var.dom <- dom.design[var.perm, i] }

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
			{
			if (calc.lod.full) { 
			  result.this.chr$lod.full[i] <- log10.lik.full - log10.lik.null
			}
			if (calc.lod.mean) { 
			  result.this.chr$lod.mean[i] <- log10.lik.full - log10.lik.nomean
			}
			if (calc.lod.var) {
			  result.this.chr$lod.var[i] <- log10.lik.full - log10.lik.novar
			}
			}
			
			if (return.effects) {

			  # baselines
				result.this.chr$mean.baseline[i] <- coef(d.fit.full)[1]
				result.this.chr$var.baseline[i] <- coef(d.fit.full$dispersion.fit)[1]
				
				# qtl and covariate effects on mean
				# I suspect there is a better way to do this....
				for (effect.name in names(coef(d.fit.full))[-1]) {
				  result.this.chr[[effect.name]][i] <- coef(d.fit.full)[[effect.name]]
				}
				
				# qtl and covariate effects on mean
				# I suspect there is a better way to do this....
				for (effect.name in names(coef(d.fit.full$dispersion.fit))[-1]) {
				  result.this.chr[[effect.name]][i] <- 
				    coef(d.fit.full$dispersion.fit)[[effect.name]]
				}
				
			}

		}
		
		if (is.null(scan)) {
		  scan <- result.this.chr
		}	else {
		  scan <- rbind(scan, result.this.chr)
		}
	}

	
	
	
	# formalize output
	scan$chr <- factor(scan$chr,
	                   levels = mixedsort(unique(scan$chr)),
	                   ordered = TRUE)
	
	class(scan) <- c("scanonevar", "tbl_df", "data.frame")
	attr(scan, 'method') <- 'dglm package'
	attr(scan, 'type') <- crosstype
	attr(scan, 'model') <- 'normal dglm'
	attr(scan, 'pheno') <- names(X)[1]
	if (calc.lod.full) { attr(scan, 'null.effects') <- null.effects }
	attr(scan, 'units') <- 'lods'
	attr(scan, 'mean.null.formula') <- mean.null.formula
	attr(scan, 'var.null.formula') <- var.null.formula
	attr(scan, 'mean.alt.formula') <- mean.alt.formula
	attr(scan, 'var.alt.formula') <- var.alt.formula

	return(scan)
}
