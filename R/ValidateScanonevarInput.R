ValidateScanonevarInput <- function(cross,
																		pheno.name,
																		mean.covar.names,
																		var.covar.names)
{

	# check cross type
	crosstype <- class(cross)[1]
	if (!(crosstype %in% c("bc", "dh", "f2", "haploid", "risib", "riself"))) {
		stop('scanonevar not implemented for cross type "', crosstype, '"')
	}

	# calc genotype probabilities if needed
	if (!("prob" %in% names(cross$geno[[1]]))) {
		message("Setup: Running calc.genoprob()")
		cross <- calc.genoprob(cross, step = 2.0)
	}

	# grab phenotype
	phen.names <- names(cross$pheno)
	if (!(pheno.name %in% phen.names)) {
		stop(paste(pheno.name, 'not found in cross'))
	}
	if (length(pheno.name) != 1) {
		stop('Scanonevar takes exactly one phenotype')
	}
	message(paste('Setting up scanonevar on phenotype:', pheno.name))
	pheno <- subset(cross$pheno, select = pheno.name)



	# store all marker names and ensure they are valid R variable names
	marker.names <- unlist(lapply(X = cross$geno,
																FUN = function(chr) { colnames(chr$data)}))
  stopifnot(all(marker.names == make.names(marker.names)))
  
	# split covar names into phenotypes and markers
	mean.covar.names <- SplitVecByMembership(v = mean.covar.names,
																					 a = phen.names,
																					 b = marker.names)
	mean.covar.phen.names <- mean.covar.names[[1]]
	mean.covar.marker.names <- mean.covar.names[[2]]


	var.covar.names <- SplitVecByMembership(v = var.covar.names,
																					a = phen.names,
																					b = marker.names)
	var.covar.phen.names <- var.covar.names[[1]]
	var.covar.marker.names <- var.covar.names[[2]]

	# pull the data specified by the covariate names
	mean.phen.covars <- data.frame(row.names = 1:nrow(pheno))
	if (length(mean.covar.phen.names)) {
		mean.phen.covars <- cross$pheno[, mean.covar.phen.names, drop = FALSE]
	}

	mean.marker.covars <- data.frame(row.names = 1:nrow(pheno))
	if (length(mean.covar.marker.names)) {

		for (mean.covar.marker.name in mean.covar.marker.names) {

			mean.marker.covars <-
				cbind(mean.marker.covars,
							GetGenoprobsByMarkerName(cross = cross,
																			 marker.name = mean.covar.marker.name)[,-1])
		}
	}

	var.phen.covars <- data.frame(row.names = 1:nrow(pheno))
	if (length(var.covar.phen.names)) {
		var.phen.covars <- cross$pheno[, var.covar.phen.names, drop = FALSE]
	}

	var.marker.covars <- data.frame(row.names = 1:nrow(pheno))
	if (length(var.covar.marker.names)) {
		for (var.covar.marker.name in var.covar.marker.names) {
			var.marker.covars <-
				cbind(var.marker.covars,
							GetGenoprobsByMarkerName(cross = cross,
																			 marker.name = var.covar.marker.name)[,-1])
		}
	}

	mean.covars <- cbind(mean.phen.covars, mean.marker.covars)
	var.covars <- cbind(var.phen.covars, var.marker.covars)

	return(list(cross = cross,
							pheno = pheno,
							mean.covars = mean.covars,
							var.covars = var.covars))
}
