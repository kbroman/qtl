AssembleDesignMatAndFormulae <- function(pheno, 
                                         mean.add,
                                         mean.dom,
                                         var.add,
                                         var.dom,
                                         mean.covars = NULL, 
                                         var.covars = NULL) {

	# set up data and formulas
  n <- length(pheno)
	X <- data.frame(pheno = pheno)
	
	mean.formula <- mean.null.formula <- paste(names(pheno), '~ 1')
	var.formula <- var.null.formula <- ' ~ 1'
	
	if (mean.add) { 
	  X$mean.add <- rep(NA, n) 
	  mean.formula <- paste(names(pheno), '~ mean.add')
	}
	if (mean.dom) { 
	  X$mean.dom <- rep(NA, n)
	  mean.formula <- paste(mean.formula, '+ mean.dom')
	}
	if (var.add) { 
	  X$var.add <- rep(NA, n)
	  var.formula <- '~ var.add'
	}
	if (var.dom) {
	  X$var.dom <- rep(NA, n)
	  var.formula <- paste(var.formula, '+ var.dom')
	}

	# add covariates to X and mean formulae
	if (length(mean.covars)) {
		ncolX <- ncol(X)
		X <- cbind(X, mean.covars)
		colnames(X)[-(1:ncolX)] <- names(mean.covars)
		mean.formula <- paste(mean.formula, "+",
													paste(names(mean.covars), collapse = " + "))
		mean.null.formula <- paste(mean.null.formula, "+",
															 paste(names(mean.covars), collapse = " + "))
	}

	# add covariates to X and var formulae
	if (length(var.covars)) {
		ncolX <- ncol(X)
		X <- cbind(X, var.covars)
		colnames(X)[-(1:ncolX)] <- names(var.covars)
		var.formula <- paste(var.formula, "+",
												 paste(names(var.covars), collapse = " + "))
		var.null.formula <- paste(var.null.formula, "+",
															paste(names(var.covars), collapse = " + "))
	}

	# get everything in the right type
	X <- as.data.frame(X)
	X <- X[, !duplicated(names(X))]
	mean.null.formula <- as.formula(mean.null.formula)
	var.null.formula <- as.formula(var.null.formula)
	mean.formula <- as.formula(mean.formula)
	var.formula <- as.formula(var.formula)

	return(list(X = X,
							mean.null.formula = mean.null.formula,
							var.null.formula = var.null.formula,
							mean.alt.formula = mean.formula,
							var.alt.formula = var.formula))
}
