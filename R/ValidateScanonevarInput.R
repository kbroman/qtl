ValidateScanonevarInput <- function(cross,
																		pheno.name,
																		mean.add,
																		mean.dom,
																		var.add,
																		var.dom,
																		mean.covar.names,
																		var.covar.names)
{

  # VALIDATE THE CROSS
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

	
	# VALIDATE THE PHENOTYPE
	# grab phenotype
	phen.names <- names(cross$pheno)
	if (!(pheno.name %in% phen.names)) {
		stop(paste(pheno.name, 'not found in cross'))
	}
	if (length(pheno.name) != 1) {
		stop('Scanonevar takes exactly one phenotype')
	}
	pheno <- subset(cross$pheno, select = pheno.name)


  # VALIDATE THE GENETIC DESIGN
	if (!mean.add & mean.dom) {
	  stop('Doesnt make sense to fit dominance effect without additive effect.  You tried to do this with regard to phenotype mean.')
	}
	if (!var.add & var.dom) {
	  stop('Doesnt make sense to fit dominance effect without additive effect.  You tried to do this with regard to phenotype variance')
	}

	# VALIDATE AND PULL COVARIATES
	# store all marker names and ensure they are valid R variable names
	marker.names <- unlist(lapply(X = cross$geno,
																FUN = function(chr) { colnames(chr$data)}))
  stopifnot(all(marker.names == make.names(marker.names)))
  
	# split mean covar names into phenotypes and markers
	mean.covar.names <- SplitVecByMembership(v = mean.covar.names,
																					 a = phen.names,
																					 b = marker.names)
	mean.covar.phen.names <- mean.covar.names[[1]]
	mean.covar.marker.names <- mean.covar.names[[2]]

	# pull the data for phenotype covariates
	mean.phen.covars <- data.frame(row.names = 1:nrow(pheno))
	if (length(mean.covar.phen.names)) {
		mean.phen.covars <- cross$pheno[, mean.covar.phen.names, drop = FALSE]
	}
	
	# pull the data for marker covariates
	mean.marker.covars <- data.frame(row.names = 1:nrow(pheno))
	if (length(mean.covar.marker.names)) {

		for (mean.covar.marker.name in mean.covar.marker.names) {

			mean.marker.covars <-
				cbind(mean.marker.covars,
							GetGenoprobsByMarkerName(cross = cross,
																			 marker.name = mean.covar.marker.name)[,-1])
		}
	}

	
	# same logic for variance covariates
	var.covar.names <- SplitVecByMembership(v = var.covar.names,
	                                        a = phen.names,
	                                        b = marker.names)
	var.covar.phen.names <- var.covar.names[[1]]
	var.covar.marker.names <- var.covar.names[[2]]
	
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

	# compile all covariates
	mean.covars <- cbind(mean.phen.covars, mean.marker.covars)
	names(mean.covars) <- paste0('mean.', names(mean.covars))
	var.covars <- cbind(var.phen.covars, var.marker.covars)
	names(var.covars) <- paste0('var.', names(var.covars))

	return(list(cross = cross,
							pheno = pheno,
							mean.covars = mean.covars,
							var.covars = var.covars))
}
