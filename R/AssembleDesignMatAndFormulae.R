assemble.df.and.formulae <- function(pheno, dom, mean.covars = NULL, var.covars = NULL) {

	# set up data and formulas
	X <- cbind(pheno = pheno,
						 mean.add.genet = rep(0, length(pheno)),
						 var.add.genet = rep(0, length(pheno)))

	mean.formula <- paste(names(pheno), '~ mean.add.genet')
	var.formula <- paste(names(pheno), '~ var.add.genet')
	mean.null.formula <- var.null.formula <- paste(names(pheno), '~ 1')

	# add two more columns to X and change alternative formulae if dom
	if (dom) {
		X <- cbind(X,
							 mean.dom.genet = rep(0, length(pheno)),
							 var.dom.genet = rep(0, length(pheno)))
		mean.formula <- paste(names(pheno), '~ mean.add.genet + mean.dom.genet')
		var.formula <- paste(names(pheno), '~ var.add.genet + var.dom.genet')
	}

	# add covariates to X and mean formulae
	if(!is.null(mean.covars)) {
		ncolX <- ncol(X)
		X <- cbind(X, mean.covars)
		colnames(X)[-(1:ncolX)] <- names(mean.covars)
		mean.formula <- paste(mean.formula, "+", paste(names(mean.covars), collapse=" + "))
		mean.null.formula <- paste(mean.null.formula, "+", paste(names(mean.covars), collapse=" + "))
	}

	# add covariates to X and var formulae
	if(!is.null(var.covars)) {
		ncolX <- ncol(X)
		X <- cbind(X, var.covars)
		colnames(X)[-(1:ncolX)] <- names(var.covars)
		var.formula <- paste(var.formula, "+", paste(names(var.covars), collapse=" + "))
		var.null.formula <- paste(var.null.formula, "+", paste(names(var.covars), collapse=" + "))
	}

	# get everything in the right type
	X <- as.data.frame(X)
	mean.null.formula <- as.formula(mean.null.formula)
	var.null.formula <- as.formula(var.null.formula)
	mean.formula <- as.formula(mean.formula)
	var.formula <- as.formula(var.formula)

	return(list(X, mean.null.formula, var.null.formula, mean.formula, var.formula))
}
