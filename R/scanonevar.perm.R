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

scanonevar.perm <- function(cross,
														pheno.name,
														chrs = 1:length(cross$geno),
														mean.covar.names = NULL,
														var.covar.names = NULL,
														mean.add = TRUE,
														mean.dom = TRUE,
														var.add = TRUE,
														var.dom = TRUE,
														quiet = TRUE,
														seed = 27599,
														num.perms = 1,
														mean.perm = TRUE,
														var.perm = TRUE,
														meanvar.perm = TRUE,
														num.processes = 1)
{

	set.seed(seed)

 	# check that pheno.name and covar.names exist in the cross and that cross is of an acceptable type
	validated.input <- ValidateScanonevarInput(cross = cross,
																						 pheno.name = pheno.name,
																						 mean.add = mean.add,
																						 mean.dom = mean.dom,
																						 var.add = var.add,
																						 var.dom = var.dom,
																						 mean.covar.names = mean.covar.names,
																						 var.covar.names = var.covar.names)
	cross <- validated.input$cross
	pheno <- validated.input$pheno
	mean.covars <- validated.input$mean.covars
	var.covars <- validated.input$var.covars

	# set up design matrix (X) and formulae for use in dglm
	X.and.formulae <- AssembleDesignMatAndFormulae(pheno = pheno,
	                                               mean.add = mean.add,
	                                               mean.dom = mean.dom,
	                                               var.add = var.add,
	                                               var.dom = var.dom,
	                                               mean.covars = mean.covars,
	                                               var.covars = var.covars)
	X <- X.and.formulae$X
	mean.null.formula <- X.and.formulae$mean.null.formula
	var.null.formula <- X.and.formulae$var.null.formula
	mean.alt.formula <- X.and.formulae$mean.alt.formula
	var.alt.formula <- X.and.formulae$var.alt.formula

	# iterate through num.perms, calling the workhorse of scanonevar.perm() each time
	# do this in parallel if num.cores > 1, and in series otherwise
 	if (num.processes == 1) {
 	  scans <- vector(mode = 'list', length = num.perms)
 	  for (perm.num in 1:num.perms) {
 	    
 	    if (!quiet) { print(paste('Starting permutation number', perm.num, '...')) }
 	    perm <- sample(nrow(pheno))
 	    
 	    meanvarperm <- scanonevar_(cross = cross, 
 	                               chrs = chrs, 
 	                               mean.add = mean.add,
 	                               mean.dom = mean.dom,
 	                               var.add = var.add,
 	                               var.dom = var.dom,
 	                               X = X,
 	                               mean.null.formula = mean.null.formula,
 	                               var.null.formula = var.null.formula,
 	                               mean.alt.formula = mean.alt.formula,
 	                               var.alt.formula = var.alt.formula,
 	                               mean.perm = perm, 
 	                               var.perm = perm,
 	                               return.effects = FALSE,
 	                               calc.lod.full = TRUE,
 	                               calc.lod.mean = FALSE, 
 	                               calc.lod.var = FALSE)
 	    
 	    meanperm <- scanonevar_(cross = cross, 
 	                            chrs = chrs, 
 	                            mean.add = mean.add,
 	                            mean.dom = mean.dom,
 	                            var.add = var.add,
 	                            var.dom = var.dom,
 	                            X = X,
 	                            mean.null.formula = mean.null.formula,
 	                            var.null.formula = var.alt.formula,
 	                            mean.alt.formula = mean.alt.formula,
 	                            var.alt.formula = var.alt.formula,
 	                            mean.perm = perm, 
 	                            return.effects = FALSE,
 	                            calc.lod.full = FALSE,
 	                            calc.lod.mean = TRUE, 
 	                            calc.lod.var = FALSE)
 	    
 	    
 	    varperm <- scanonevar_(cross = cross, 
 	                           chrs = chrs, 
 	                           mean.add = mean.add,
 	                           mean.dom = mean.dom,
 	                           var.add = var.add,
 	                           var.dom = var.dom,
 	                           X = X,
 	                           mean.null.formula = mean.alt.formula,
 	                           var.null.formula = var.null.formula,
 	                           mean.alt.formula = mean.alt.formula,
 	                           var.alt.formula = var.alt.formula,
 	                           var.perm = perm,
 	                           return.effects = FALSE,
 	                           calc.lod.full = FALSE,
 	                           calc.lod.mean = FALSE, 
 	                           calc.lod.var = TRUE)
 	    
 	    n.loci <- nrow(meanvarperm)
 	    result <- cbind(perm.idx = rep(perm.num, n.loci),
 	                    as.data.frame(meanvarperm[,c('chrtype', 'chr', 'pos', 'lod.full')]),
 	                    as.data.frame(meanperm[, 'lod.mean', drop = FALSE]),
 	                    as.data.frame(varperm[, 'lod.var', drop = FALSE]))
 	    
 	    result <- tbl_df(result) %>%
 	      group_by(perm.idx, chrtype) %>%
 	      summarize(max.lod.full = max(lod.full),
 	                max.lod.mean = max(lod.mean),
 	                max.lod.var = max(lod.var))
 	    scans[[perm.num]] <- result
 	  }
 	}
	
	if (num.processes > 1) {
	  
	  # set up parllel backend
	  cl <- makeCluster(num.processes)
	  registerDoParallel(cl)
	  
	  scans <- foreach(perm.num = 1:num.perms) %dopar% {
	    
	    if (!quiet) { print(paste('Starting permutation number', perm.num, '...')) }
	    perm <- sample(nrow(pheno))
	    
	    meanvarperm <- scanonevar_(cross = cross, 
	                               chrs = chrs, 
	                               mean.add = mean.add,
	                               mean.dom = mean.dom,
	                               var.add = var.add,
	                               var.dom = var.dom,
	                               X = X,
	                               mean.null.formula = mean.null.formula,
	                               var.null.formula = var.null.formula,
	                               mean.alt.formula = mean.alt.formula,
	                               var.alt.formula = var.alt.formula,
	                               mean.perm = perm, 
	                               var.perm = perm,
	                               return.effects = FALSE,
	                               calc.lod.full = TRUE,
	                               calc.lod.mean = FALSE, 
	                               calc.lod.var = FALSE)
	    
	    meanperm <- scanonevar_(cross = cross, 
	                            chrs = chrs, 
	                            mean.add = mean.add,
	                            mean.dom = mean.dom,
	                            var.add = var.add,
	                            var.dom = var.dom,
	                            X = X,
	                            mean.null.formula = mean.null.formula,
	                            var.null.formula = var.alt.formula,
	                            mean.alt.formula = mean.alt.formula,
	                            var.alt.formula = var.alt.formula,
	                            mean.perm = perm, 
	                            return.effects = FALSE,
	                            calc.lod.full = FALSE,
	                            calc.lod.mean = TRUE, 
	                            calc.lod.var = FALSE)
	    
	    
	    varperm <- scanonevar_(cross = cross, 
	                           chrs = chrs, 
	                           mean.add = mean.add,
	                           mean.dom = mean.dom,
	                           var.add = var.add,
	                           var.dom = var.dom,
	                           X = X,
	                           mean.null.formula = mean.alt.formula,
	                           var.null.formula = var.null.formula,
	                           mean.alt.formula = mean.alt.formula,
	                           var.alt.formula = var.alt.formula,
	                           var.perm = perm,
	                           return.effects = FALSE,
	                           calc.lod.full = FALSE,
	                           calc.lod.mean = FALSE, 
	                           calc.lod.var = TRUE)
	    
	    n.loci <- nrow(meanvarperm)
	    result <- cbind(perm.idx = rep(perm.num, n.loci),
	                    as.data.frame(meanvarperm[,c('chrtype', 'chr', 'pos', 'lod.full')]),
	                    as.data.frame(meanperm[, 'lod.mean', drop = FALSE]),
	                    as.data.frame(varperm[, 'lod.var', drop = FALSE]))
	    
	    result <- tbl_df(result) %>%
	      group_by(perm.idx, chrtype) %>%
	      summarize(max.lod.full = max(lod.full),
	                max.lod.mean = max(lod.mean),
	                max.lod.var = max(lod.var))
		
		return(result)
	  }
	  # close down parallel backend
	  stopCluster(cl)
	}

	return(rbind_all(scans))
}
