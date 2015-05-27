ValidateScanonevarInput <- function(cross,
																		pheno.name,
																		mean.covar.names,
																		var.covar.names)
{

	# check cross type
	crosstype <- class(cross)[1]
	if(!(crosstype %in% c("bc", "dh", "f2", "haploid", "risib", "riself"))) {
		stop('scanonevar not implemented for cross type "', crosstype, '"')
	}

	# grab phenotype
	if (!(pheno.name %in% names(cross$pheno))) {
		stop(paste(pheno.name, 'not found in cross'))
	}
	message(paste('Setting up scanonevar on phenotype:', pheno.name))
	pheno <- subset(cross$pheno, select = pheno.name)

	# grab names of all markers (in case a covariate is a marker)
	marker.names <- unlist(lapply(X = cross$geno,
																FUN = function(chr) { colnames(chr$data)}))

	# grab names of all phenotypes (in case covariate is a phenotype)
	phen.names <- names(cross$pheno)

	# split covariate names into marker names and phenotype names
	mean.covar.marker.idxs <- which(mean.covar.names %in% marker.names)
	mean.covar.phen.idxs <- which(mean.covar.names %in% phen.names)

	# if any covariate name is in neither list, error
	if (any(and(mean.covar.marker.idxs, mean.covar.phen.idxs))) {
		ambig.covar.namez <- which(and(mean.covar.marker.idxs,
																	 mean.covar.phen.idxs))
		error(paste('Ambiguity:  The following covariate names could refer',
								'to a phenotype or a marker',
								ambig.covar.names))
	}

	# if any covariate name is in both lists, error
	if (any(!or(mean.covar.marker.idxs, mean.covar.phen.idxs))) {
		missing.covar.namez <- which(!or(mean.covar.marker.idxs,
																		 mean.covar.phen.idxs))
		error(paste('The following covariate names could not be found',
								'in the list of phenotypes and markers',
								missing.covar.names))
	}

	# at this point each mean covar name is in exactly one list
	mean.covar.marker.names <- mean.covar.names[mean.covar.marker.idxs]
	mean.covar.phen.names <- mean.covar.names[mean.covar.phen.idxs]

	# grab mean covars
	mean.covars <- NULL
	if (!is.null(mean.covar.names)) {
		if (!all(mean.covar.names %in% names(cross$pheno))) {
			stop(paste('at least one mean covariate is not found in the cross'))
		}
		mean.covars <- select(cross$pheno, mean.covar.names)
	}

	# grab var covars
	var.covars <- NULL
	if (!is.null(var.covar.names)) {
		if (!all(var.covar.names %in% names(cross$pheno))) {
			stop(paste('at least one var covariate is not found in the cross'))
		}
		var.covars <- subset(cross$pheno, select = var.covar.names)
	}

	# calc genotype probabilities if needed
	if(!("prob" %in% names(cross$geno[[1]]))) {
		warning("First running calc.genoprob")
		cross <- calc.genoprob(cross)
	}

	return(list(pheno, mean.covars, var.covars))
}
