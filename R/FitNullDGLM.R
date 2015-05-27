fit.null.dglm <- function(X, null.formulae) {

	mean.null.formula <- null.formulae[[1]]
	var.null.formula <- null.formulae[[2]]

	# fit covariate-only model (covariates, but not genetic marker, may have mean and varianc effect)
	# only need to do this once per scan
	d.fit.null <- dglm(formula = mean.null.formula,
										 dformula = var.null.formula,
										 data = X)
	ln.lik.null <- -0.5*d.fit.null$m2loglik
	log10.lik.null <- ln.lik.null / log(10)
	null.effects <- list(meancovs = d.fit.null$coef[-1],
											 varcovs = d.fit.null$disp$coef[-1])

	return(list(log10.lik.null, null.effects))
}
