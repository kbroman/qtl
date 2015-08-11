#####################################################################
#
# scanonevar.R
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
# Contains: scanonevar,
#
######################################################################

######################################################################
#
# scanonevar:
# single-QTL genome scan for QTL affecting trait mean or variance
#
######################################################################

scanonevar <- function(cross,
											 pheno.name,
											 mean.add.covar.names = NULL,
											 var.covar.names = NULL,
											 mean.int.covar.names = NULL,
											 var.int.covar.names = NULL,
											 chrs = 1:length(cross$geno),
											 mean.add = TRUE,
											 mean.dom = TRUE,
											 var.add = TRUE,
											 var.dom = TRUE,
											 return.effects = TRUE,
											 quiet = TRUE)
{

	# check that pheno.name and covar.names exist in the cross
	# and that cross is of an acceptable type
	# return the phenotype and covariate data (not names)
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

	# workhorse of scanonevar()
	# applies dglm to each locus and compares to various null fits to return LOD for each locus
	scan <- scanonevar_(cross = cross,
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
											return.effects = return.effects,
											quiet = quiet)

	return(scan)
}
