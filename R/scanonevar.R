# scanonevar
# single-QTL genome scan for QTL affecting variance
# with code from Lars Ronnegard

scanonevar <-
	function(cross, pheno.col=1, mean_covar = NULL, var_covar = NULL,
					 #use.dglm.package = TRUE, use.custom.em = FALSE,
					 dom = FALSE, maxit = 25 , tol=1e-6, quiet=TRUE, chrs)
	{

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
		X <- cbind(pheno = pheno, add = rep(0, length(pheno)))
		mean_formula <- var_formula <- "pheno ~ add"

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
		}

		if(!is.null(var_covar)) {
			ncolX <- ncol(X)
			X <- cbind(X, var_covar)
			varcovarnames <- paste0("varcov", 1:(ncol(X)-ncolX))
			colnames(X)[-(1:ncolX)] <- varcovarnames
			var_formula <- paste(var_formula, "+",
													 paste(varcovarnames, collapse="+"))
		}

		X <- as.data.frame(X)
		mean_formula <- as.formula(mean_formula)
		var_formula <- as.formula(var_formula)

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

			logP.mean.add <- logP.disp.add <- numeric(n.loci)
			if (dom) { logP.mean.dom <- logP.disp.dom <- numeric(n.loci) }

			mean.baseline <- disp.baseline <- numeric(n.loci)
			mean.add.effect <- disp.add.effect <- numeric(n.loci)
			if (dom) { mean.dom.effect <- disp.dom.effect <- numeric(n.loci) }

			for(i in 1:n.loci) { # loop over positions within chromosome

				# fill in genotype probs for this locus
				X[,2] <- a1[,i]
				if (dom) { X[,3] <- d1[,i] }

# 				if (use.dglm.package) {
#
# 					d.fit <- dglm(formula = mean_formula,
# 												dformula = var_formula,
# 												data = X)
#
# 					mean.baseline[i] <- coef(d.fit)[1]
# 					disp.baseline[i] <- coef(d.fit$dispersion.fit)[1]
# 					logP.mean.add[i] <- -log10(summary(d.fit)$coefficients[2,4])
# 					logP.disp.add[i] <- -log10(summary(d.fit)$dispersion.summary$coefficients[2,4])
# 					mean.add.effect[i] <- coef(d.fit)[2]
# 					disp.add.effect[i] <- coef(d.fit$dispersion.fit)[2]
#
# 					if (dom) {
# 						mean.dom.effect[i] <- coef(d.fit)[3]
# 						disp.dom.effect[i] <- coef(d.fit$dispersion.fit)[3]
# 						logP.mean.dom[i] <- -log10(summary(d.fit)$coefficients[3,4])
# 						logP.disp.dom[i] <- -log10(summary(d.fit)$dispersion.summary$coefficients[3,4])
# 					}
# 				}

				#if (use.custom.em) {

					d.fit <- DGLM_norm(m.form  = mean_formula,
														 d.form  = var_formula,
														 indata  = X,
														 maxiter = maxit,
														 conv    = tol)

					mean.baseline[i] <- summary(d.fit$mean)$coef[1,1]
					disp.baseline[i] <- summary(d.fit$disp)$coef[1,1]
					logP.mean.add[i] <- -log10(summary(d.fit$mean)$coef[2,4])
					logP.disp.add[i] <- -log10(summary(d.fit$disp)$coef[2,4])
					mean.add.effect[i] <- summary(d.fit$mean)$coef[2,1]
					disp.add.effect[i] <- summary(d.fit$disp)$coef[2,1]

					if (dom) {
						mean.dom.effect[i] <- summary(d.fit$mean)$coef[3,1]
						disp.dom.effect[i] <- summary(d.fit$disp)$coef[3,1]
						logP.mean.dom[i] <- -log10(summary(d.fit$mean)$coef[3, 4])
						logP.disp.dom[i] <- -log10(summary(d.fit$disp)$coef[3, 4])
					}

					if (d.fit$iter >= maxit) {
						logP.disp.add[i] <-  0
						if (dom) {
							logP.disp.dom[i] <- 0
						}
						warning("dglm did not converge on chr", chr.names[j], " position ", i)
					}

				#}

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
														neglogP_mean_add = logP.mean.add,
														neglogP_disp_add = logP.disp.add,
														mean_add_effect = mean.add.effect,
														disp_add_effect = disp.add.effect,
														stringsAsFactors = FALSE)
			if (dom) {
				thischr$neglogP_mean_dom = logP.mean.dom
				thischr$neglogP_disp_dom = logP.disp.dom
				thischr$mean_dom_effect = mean.dom.effect
				thischr$disp_dom_effect = disp.dom.effect
			}

			rownames(thischr) <- w

			if(is.null(result)) result <- thischr
			else result <- rbind(result, thischr)
		}

		class(result) <- c("scanonevar", "scanone", "data.frame")
		attr(result, "method") <- "scanonevar"
		attr(result, 'dom') <- dom
		attr(result, 'pheno') <- names(cross$pheno)[pheno.col]

		result
	}

DGLM_norm <- function(m.form, d.form, indata, maxiter=20, conv=1e-6) {
    X.mean <- model.matrix(m.form, data = indata)
    X.disp <- model.matrix(d.form, data = indata)
    y.name <- all.vars(m.form)[1]
    y <- indata[,y.name]
    w <- rep(1, nrow(indata))
    convergence <- 1
    iter <- 0
    while (convergence > conv & iter < maxiter) {
        iter <- iter + 1
        w.old <- w
        glm1 <- lm(y~.-1, weights=w, data=data.frame(X.mean))
        res <- resid(glm1)
        q <- hatvalues(glm1)
        y2 <- res^2/(1-q)
        glm2 <- glm(y2~.-1, family=Gamma(link=log), weights=(1-q)/2, data=data.frame(X.disp))
        w <- 1/fitted(glm2)
        convergence <- (max(abs(w.old-w)) + (summary(glm1)$sigma-1) )
    }
    return(list(mean=glm1, disp=glm2, iter=iter))
}
