#####################################################################
#
# scanonevar.R
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
# Contains: scanonevar,
#
######################################################################

######################################################################
#
# scanonevar:
# single-QTL genome scan for QTL affecting trait mean or variance
#
######################################################################

scanonevar <-
	function(cross, pheno.col=1, mean_covar = NULL, var_covar = NULL,
					 dom = TRUE, maxit = 25 , tol=1e-6, quiet=TRUE, chrs)
	{

		# check input
		crosstype <- class(cross)[1]
		if(!(crosstype %in% c("bc", "dh", "f2", "haploid", "risib", "riself"))) {
			stop('scanonevar not implemented for cross type "', crosstype, '"')
		}

# 		chrtype <- sapply(cross$geno, class)
# 		if(any(chrtype=="X")) {
# 			warning("Analysis of X chromosome not implemented for scanonevar; omitted.")
# 			cross <- subset(cross, chr=(chrtype != "X"))
# 		}

		# grab phenotype
		if(LikePheVector(pheno.col, nind(cross), nphe(cross))) {
			cross$pheno <- cbind(pheno.col, cross$pheno)
			pheno.col <- 1
		}
		if(is.character(pheno.col)) {
			num <- find.pheno(cross, pheno.col)
			if(any(is.na(num))) {
				if(sum(is.na(num)) > 1)
					stop("Couldn't identify phenotypes ", paste(paste("\"", pheno.col[is.na(num)], "\"", sep=""),
																											collapse=" "))
				else
					stop("Couldn't identify phenotype \"", pheno.col[is.na(num)], "\"")
			}
			pheno.col <- num
		}
		if(any(pheno.col < 1 | pheno.col > nphe(cross))) {
			stop("pheno.col values should be between 1 and the no. phenotypes")
		}
		pheno <- cross$pheno[,pheno.col]
		if(is.matrix(pheno) && ncol(pheno) > 1) {
			pheno <- pheno[,1]
			warning('scanonevar requires a single phenotype; all but "', phenames(cross)[pheno.col[1]], '" omitted.')
		}

		N <- length(pheno) # No. individuals
		n.chr <- nchr(cross) #No. chromosomes
		chr.names <- chrnames(cross)

		# need to run calc.genoprob?
		if(!("prob" %in% names(cross$geno[[1]]))) {
			warning("First running calc.genoprob")
			cross <- calc.genoprob(cross)
		}

		scan.logPm <- scan.logPd <- chr.names.out <- NULL

		# set up data and formulas
		X <- cbind(pheno = pheno, add = rep(0, length(pheno)))
		mean_formula <- var_formula <- "pheno ~ add"
		mean_null_formula <- var_null_formula <- "pheno ~ 1"

		if (dom) {
			X <- cbind(X, dom = rep(0, length(pheno)))
			mean_formula <- var_formula <- "pheno ~ add + dom"
		}

		# todo: give the same treatment to covariates as we gave to phenotype earlier
		# i.e., do some searching to figure out if the user provided a name of a pheno column
		if(!is.null(mean_covar)) {
			ncolX <- ncol(X)
			X <- cbind(X, mean_covar)
			meancovarnames <- paste0("meancov", 1:(ncol(X)-ncolX))
			colnames(X)[-(1:ncolX)] <- meancovarnames
			mean_formula <- paste(mean_formula, "+",
														paste(meancovarnames, collapse="+"))
			mean_null_formula <- paste(mean_null_formula, "+",
																 paste(meancovarnames, collapse="+"))

		}

		if(!is.null(var_covar)) {
			ncolX <- ncol(X)
			X <- cbind(X, var_covar)
			varcovarnames <- paste0("varcov", 1:(ncol(X)-ncolX))
			colnames(X)[-(1:ncolX)] <- varcovarnames
			var_formula <- paste(var_formula, "+",
													 paste(varcovarnames, collapse="+"))
			var_null_formula <- paste(var_null_formula, "+",
																paste(varcovarnames, collapse="+"))
		}

		X <- as.data.frame(X)
		mean_null_formula <- as.formula(mean_null_formula)
		var_null_formula <- as.formula(var_null_formula)
		mean_formula <- as.formula(mean_formula)
		var_formula <- as.formula(var_formula)

		# fit covariate-only model (covariates, but not genetic marker, may have mean and varianc effect)
		# only need to do this once per scan
		d.fit.null <- dglm(formula = mean_null_formula,
											 dformula = var_null_formula,
											 data = X)
		ln.lik.null <- -0.5*d.fit.null$m2loglik
		log10.lik.null <- ln.lik.null / log(10)
		null.effects <- list(meancovs = d.fit.null$coef[-1],
												 varcovs = d.fit.null$disp$coef[-1])

		result <- NULL
		if (missing(chrs)) { chrs <- 1:length(cross$geno) }

		for(j in chrs) { # loop over chromosomes

			if(!quiet) message(" - Chr ", chr.names[j])

			this.chr <- cross$geno[[j]]

			if (class(this.chr) == 'A') {
				if (crosstype=="f2") {
					g11 <- this.chr$prob[,,1]
					g12 <- this.chr$prob[,,2]
					g13 <- this.chr$prob[,,3]
					a1  <- -g11 + g13
					d1 <-  g12
					chr.dom <- dom
				}
				else {
					a1 <- cross$geno[[j]]$prob[,,1]
				}
			}

			if (class(this.chr) == 'X') {
				chr.dom <- FALSE
				if (crosstype=="f2") {
					a1 <- this.chr$prob[,,1]
				}
				else {
					print("X not implemented for non-F2 crosses")
					next
				}
			}

			n.loci <- dim(a1)[2]

			lod.full <- lod.mean <- lod.disp <- numeric(n.loci)

			mean.baseline <- disp.baseline <- numeric(n.loci)
			mean.add.effect <- disp.add.effect <- numeric(n.loci)
			if (chr.dom) { mean.dom.effect <- disp.dom.effect <- numeric(n.loci) }

			for(i in 1:n.loci) { # loop over positions within chromosome

				# fill in genotype probs for this locus
				X[,2] <- a1[,i]
				if (chr.dom) { X[,3] <- d1[,i] }

				# fit full model
				d.fit.full <- dglm(formula = mean_formula,
													 dformula = var_formula,
													 data = X)
				ln.lik.full <- -0.5*d.fit.full$m2loglik
				log10.lik.full <- ln.lik.full / log(10)

				# fit variance-only model (covariates, but not the genetic marker, may have mean effects)
				d.fit.nomean <- dglm(formula = mean_null_formula,
														 dformula = var_formula,
														 data = X)
				ln.lik.nomean <- -0.5*d.fit.nomean$m2loglik
				log10.lik.nomean <- ln.lik.nomean / log(10)

				# fit mean-only model (covariates, but not the genetic marker, may have variance effects)
				d.fit.nodisp <- dglm(formula = mean_formula,
														 dformula = var_null_formula,
														 data = X)
				ln.lik.nodisp <- -0.5*d.fit.nodisp$m2loglik
				log10.lik.nodisp <- ln.lik.nodisp / log(10)

				lod.full[i] <- log10.lik.full - log10.lik.null
				lod.mean[i] <- log10.lik.full - log10.lik.nomean
				lod.disp[i] <- log10.lik.full - log10.lik.nodisp
				mean.baseline[i] <- coef(d.fit.full)[1]
				disp.baseline[i] <- coef(d.fit.full$dispersion.fit)[1]
				mean.add.effect[i] <- coef(d.fit.full)[2]
				disp.add.effect[i] <- coef(d.fit.full$dispersion.fit)[2]

				if (chr.dom) {
					mean.dom.effect[i] <- coef(d.fit.full)[3]
					disp.dom.effect[i] <- coef(d.fit.full$dispersion.fit)[3]
				}

			}

			# set up the output
			map <- attr(cross$geno[[j]]$prob,"map")
			w <- names(map)
			o <- grep("^loc-*[0-9]+",w)
			if(length(o) > 0)  { # inter-marker locations cited as "c*.loc*"
				w[o] <- paste("c",chr.names[j],".",w[o],sep="")
			}
			thischr <- data.frame(chr = rep(chr.names[j], length(w)),
														pos = unclass(map),
														mean.baseline = mean.baseline,
														disp.baseline = disp.baseline,
														lod.full = lod.full,
														lod.mean = lod.mean,
														lod.disp = lod.disp,
														mean_add_effect = mean.add.effect,
														disp_add_effect = disp.add.effect,
														stringsAsFactors = FALSE)
			if (chr.dom) {
				thischr$mean_dom_effect = mean.dom.effect
				thischr$disp_dom_effect = disp.dom.effect
			}
			if (dom & !chr.dom) {
				thischr$mean_dom_effect <- thischr$disp_dom_effect <- rep(NA, n.loci)
			}

			rownames(thischr) <- w

			if(is.null(result)) result <- thischr
			else result <- rbind(result, thischr)
		}

		toret <- list(scan = result, null.effects = null.effects)
		class(toret) <- c("scanonevar", "scanone", "data.frame")
		attr(toret, "method") <- "scanonevar"
		attr(toret, 'dom') <- dom
		attr(toret, 'pheno') <- names(cross$pheno)[pheno.col]

		return(toret)
	}

# DGLM_norm <- function(m.form, d.form, indata, maxiter=20, conv=1e-6) {
#     X.mean <- model.matrix(m.form, data = indata)
#     X.disp <- model.matrix(d.form, data = indata)
#     y.name <- all.vars(m.form)[1]
#     y <- indata[,y.name]
#     w <- rep(1, nrow(indata))
#     convergence <- 1
#     iter <- 0
#     while (convergence > conv & iter < maxiter) {
#         iter <- iter + 1
#         w.old <- w
#         glm1 <- lm(y~.-1, weights=w, data=data.frame(X.mean))
#         res <- resid(glm1)
#         q <- hatvalues(glm1)
#         y2 <- res^2/(1-q)
#         glm2 <- glm(y2~.-1, family=Gamma(link=log), weights=(1-q)/2, data=data.frame(X.disp))
#         w <- 1/fitted(glm2)
#         convergence <- (max(abs(w.old-w)) + (summary(glm1)$sigma-1) )
#     }
#     return(list(mean=glm1, disp=glm2, iter=iter))
# }
