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

scanonevar.perm <- function(cross, pheno.name, chrs,
														mean.covar.names = NULL,
														var.covar.names = NULL,
														dom = TRUE,
														quiet = TRUE,
														seed = 27599,
														num.perms = 5,
														mean.perm = TRUE,
														var.perm = TRUE,
														meanvar.perm = TRUE)
{

	set.seed(seed)

	# check that pheno.name and covar.names exist in the cross and that cross is of an acceptable type
	validated.input <- validate.input.scanonevar(cross = cross,
																							 pheno.name = pheno.name,
																							 mean.covar.names = mean.covar.names,
																							 var.covar.names = var.covar.names)
	pheno <- validated.input[[1]]
	mean.covars <- validated.input[[2]]
	var.covars <- validated.input[[3]]

	# set up design matrix (X) and formulae for use in dglm
	df.and.formulae <- assemble.df.and.formulae(pheno = pheno,
																							dom = dom,
																							mean.covars = mean.covars,
																							var.covars = var.covars)
	X <- df.and.formulae[[1]]


	# iterate through num.perms, calling the workhorse of scanonevar.perm() each time
	scans <- vector(mode = 'list', length = num.perms)
	for (perm.num in 1:num.perms) {

		if (!quiet) { message(paste('Starting permutation number', perm.num, '...')) }
		perm <- sample(nrow(pheno))

		meanvarperm <- scanonevar_(cross = cross, chrs = chrs, dom = dom, X = X,
															 null.formulae = df.and.formulae[c(2,3)],
															 alt.formulae = df.and.formulae[c(4,5)],
															 mean.perm = perm, var.perm = perm,
															 return.effects = FALSE,
															 calc.lod.mean = FALSE, calc.lod.var = FALSE)

		meanperm <- scanonevar_(cross = cross, chrs = chrs, dom = dom, X = X,
														null.formulae = df.and.formulae[c(2,5)],
														alt.formulae = df.and.formulae[c(4,5)],
														mean.perm = perm,	return.effects = FALSE,
														calc.lod.full = FALSE, calc.lod.var = FALSE)

		varperm <- scanonevar_(cross = cross, chrs = chrs, dom = dom, X = X,
													null.formulae = df.and.formulae[c(4,3)],
													alt.formulae = df.and.formulae[c(4,5)],
													var.perm = perm, return.effects = FALSE,
													calc.lod.full = FALSE, calc.lod.mean = FALSE)

		n.loci <- nrow(meanvarperm)
		result <- cbind(perm.idx = rep(perm.num, n.loci),
										as.data.frame(meanvarperm[,c('chrtype', 'chr', 'pos', 'lod.full')]),
										as.data.frame(meanperm[, 'lod.mean', drop = FALSE]),
										as.data.frame(varperm[, 'lod.var', drop = FALSE]))

		scans[[perm.num]] <- tbl_df(result) %>%
			group_by(perm.idx, chrtype) %>%
			summarize(max.lod.full = max(lod.full),
								max.lod.mean = max(lod.mean),
								max.lod.var = max(lod.var))
		if (!quiet) { message(paste('Done with permutation number', perm.num, '...')) }
	}

	return(rbind_all(scans))
}
