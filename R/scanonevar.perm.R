#####################################################################
#
# scanonevar.perm.R
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
# Contains: scanonevar.perm,
#
######################################################################

######################################################################
#
# scanonevar.perm:
# permutations for single-QTL genome scan affecting trait mean or variance
#
######################################################################

scanonevar.perm <- function(cross, pheno.col=1, mean_covar = NULL, var_covar = NULL,
														dom = FALSE, maxit = 25 , tol=1e-6, quiet=TRUE, chrs,
														job.num = 1, num.perms = 5,
														mean.perm = TRUE, var.perm = TRUE, meanvar.perm = TRUE)
{

	set.seed(job.num)

	# check input
	crosstype <- class(cross)[1]
	if(!(crosstype %in% c("bc", "dh", "f2", "haploid", "risib", "riself"))) {
		stop('scanonevar not implemented for cross type "', crosstype, '"')
	}

	chrtype <- sapply(cross$geno, class)
	if(any(chrtype=="X")) {
		warning("Analysis of X chromosome not implemented for scanonevar; omitted.")
		cross <- subset(cross, chr=(chrtype != "X"))
	}

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

	# 		# check that we have at least one dglm-fitting method to use
	# 		if(sum(use.dglm.package, use.custom.em) == 0) {
	# 			stop("Need at least one of 'use.dglm.package' and 'use.custom.em' arguments to be TRUE")
	# 		}

	scan.logPm <- scan.logPd <- chr.names.out <- NULL

	# set up data and formulas
	X <- cbind(pheno = pheno,
						 mean_add = rep(0, length(pheno)),
						 var_add = rep(0, length(pheno)))
	mean_formula <- 'pheno ~ mean_add'
	var_formula <- 'pheno ~ var_add'
	mean_null_formula <- var_null_formula <- 'pheno ~ 1'

	if (dom) {
		X <- cbind(X, mean_dom = rep(0, length(pheno)), var_dom = rep(0, length(pheno)))
		mean_formula <- paste0(mean_formula, ' + mean_dom')
		var_formula <- paste0(var_formula, ' + var_dom')
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
															paste(meancovarnames, collapse="+"))
	}

	# X used for meanvar perms, Y for mean perms, Z for var perms
	X <- Y <- Z <- as.data.frame(X)
	mean_null_formula <- as.formula(mean_null_formula)
	var_null_formula <- as.formula(var_null_formula)
	mean_formula <- as.formula(mean_formula)
	var_formula <- as.formula(var_formula)

	all.results <- vector(mode = 'list', length = num.perms)

	for (perm.num in 1:num.perms) {

		print(paste('Starting permutation number', perm.num, '...'))
		perm1 <- sample(N)

		result <- NULL
		if (missing(chrs)) { chrs <- 1:length(cross$geno) }
		for(j in chrs) { # loop over chromosomes
			if(!quiet) message(" - Chr ", chr.names[j])

			if (crosstype=="f2") {
				g11 <- cross$geno[[j]]$prob[,,1]
				g12 <- cross$geno[[j]]$prob[,,2]
				g13 <- cross$geno[[j]]$prob[,,3]
				a1  <- -g11 + g13
				d1 <-  g12
			}
			else {
				a1 <- cross$geno[[j]]$prob[,,1]
			}

			n.loci <- dim(a1)[2]

			lod.full <- lod.mean <- lod.disp <- numeric(n.loci)

			mean.baseline <- disp.baseline <- numeric(n.loci)
			mean.add.effect <- disp.add.effect <- numeric(n.loci)
			if (dom) { mean.dom.effect <- disp.dom.effect <- numeric(n.loci) }

			for(i in 1:n.loci) { # loop over positions within chromosome

				# fill in genotype probs for this locus
				# TODO: use $mean_add naming for setting columns of X, Y, and Z for clarity
				X$mean_add <- a1[perm1, i]
				Y$mean_add <- a1[perm1, i]
				Z$mean_add <- a1[,i]

				X$var_add <- a1[perm1, i]
				Y$var_add <- a1[,i]
				Z$var_add <- a1[perm1, i]

				if (dom) {
					X$mean_dom <- d1[perm1, i]
					Y$mean_dom <- a1[perm1, i]
					Z$mean_dom <- a1[,i]

					X$var_dom <- d1[perm1, i]
					Y$var_dom <- a1[,i]
					Z$var_dom <- a1[perm1, i]
				}

				if (meanvar.perm) {

					d.fit.full <- dglm(formula = mean_formula,
														 dformula = var_formula,
														 family = gaussian,
														 data = X)
					ln.lik.full <- -0.5*d.fit.full$m2loglik
					log10.lik.full <- ln.lik.full / log(10)

					d.fit.null <- dglm(formula = mean_null_formula,
														 dformula = var_null_formula,
														 data = X)
					ln.lik.null <- -0.5*d.fit.null$m2loglik
					log10.lik.null <- ln.lik.null / log(10)

					lod.full[i] <- log10.lik.full - log10.lik.null
				}

				if (mean.perm) {

					d.fit.full <- dglm(formula = mean_formula,
														 dformula = var_formula,
														 family = gaussian,
														 data = Y)
					ln.lik.full <- -0.5*d.fit.full$m2loglik
					log10.lik.full <- ln.lik.full / log(10)

					d.fit.mean <- dglm(formula = mean_null_formula,
														 dformula = var_formula,
														 family = gaussian,
														 data = Y)
					ln.lik.mean <- -0.5*d.fit.mean$m2loglik
					log10.lik.mean <- ln.lik.mean / log(10)

					lod.mean[i] <- log10.lik.full - log10.lik.mean
				}

				if (var.perm) {

					d.fit.full <- dglm(formula = mean_formula,
														 dformula = var_formula,
														 family = gaussian,
														 data = Z)
					ln.lik.full <- -0.5*d.fit.full$m2loglik
					log10.lik.full <- ln.lik.full / log(10)

					d.fit.disp <- dglm(formula = mean_formula,
														 dformula = var_null_formula,
														 family = gaussian,
														 data = Z)
					ln.lik.disp <- -0.5*d.fit.disp$m2loglik
					log10.lik.disp <- ln.lik.disp / log(10)

					lod.disp[i] <- log10.lik.full - log10.lik.disp
				}


			}

			# set up the output
			map <- attr(cross$geno[[j]]$prob,"map")
			w <- names(map)
			o <- grep("^loc-*[0-9]+",w)
			perm.idx <- (job.num - 1)*num.perms + perm.num
			if(length(o) > 0)  { # inter-marker locations cited as "c*.loc*"
				w[o] <- paste("c",chr.names[j],".",w[o],sep="")
			}
			thischr <- data.frame(perm.idx = perm.idx,
														chr = rep(chr.names[j], length(w)),
														pos = unclass(map),
														stringsAsFactors = FALSE)
			if (meanvar.perm) { thischr$lod.full = lod.full }
			if (mean.perm) { thischr$lod.mean = lod.mean }
			if (var.perm) { thischr$lod.disp = lod.disp }

			rownames(thischr) <- w

			if(is.null(result)) result <- thischr
			else result <- rbind(result, thischr)
		}

		class(result) <- c("scanonevar", "scanone", "data.frame")
		attr(result, "method") <- "scanonevar"
		attr(result, 'dom') <- dom
		attr(result, 'pheno') <- names(cross$pheno)[pheno.col]

		all.results[[perm.num]] <- result

		print(paste('Done with permutation number', perm.num, '.'))

	}

	return(all.results)
}
