/**********************************************************************
 * 
 * fitqtl_imp.h
 *
 * copyright (c) 2002-2013, Hao Wu
 *     Modified by Karl W. Broman to get estimates of QTL effects
 *
 * last modified Feb, 2013
 * first written Apr, 2002
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 *
 * These functions are for fitting a fixed multiple-QTL model by 
 * imputation.
 *
 * Contains: R_fitqtl_imp, fitqtl_imp, nullRss0, galtRss
 *
 **********************************************************************/

void R_fitqtl_imp(int *n_ind, int *n_qtl, int *n_gen, int *n_draws,
		  int *draws, int *n_cov, double *cov, int *model, 
		  int *n_int, double *pheno, int *get_ests,
		  /* return variables */
		  double *lod, int *df, double *ests, double *ests_covar,
		  double *design_mat, int *matrix_rank);
/**********************************************************************
 * 
 * fitqtl_imp
 *
 * Fits a fixed multiple-QTL model by multiple imputation.
 * 
 * n_ind        Number of individuals
 *
 * n_qtl        Number of QTLs in the model 
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of impiutations
 *
 * Draws        Array of genotype imputations, indexed as 
 *              Draws[draw][mar][ind]
 *
 * Cov          covariates matrix, Cov[mar][ind]
 *
 * n_cov        Number of covariates
 *
 * model        Model matrix
 *
 * n_int        Number of interactions in the model
 *
 * pheno        Phenotype data, as a vector
 *
 * get_ests     0/1: If 1, return estimated effects and their variances
 *
 * lod          Return LOD score
 *
 * df           Return degree of freedom
 *
 * ests         Return ests (vector of length sizefull)
 *
 * ests_covar   Return covariance matrix of ests (sizefull^2 matrix)
 *
 * matrix_rank  Return min (across imputations) of rank of design matrix
 *
 **********************************************************************/

void fitqtl_imp(int n_ind, int n_qtl, int *n_gen, int n_draws, 
		int ***Draws, double **Cov, int n_cov, 
		int *model, int n_int, double *pheno, int get_ests,
		double *lod, int *df, double *ests, double *ests_covar,
		double *design_mat, int *matrix_rank);

/* function to calculate the null model RSS. This function is different
   from the function used in scanone_imp and scantwo_imp, which contain 
   the additive covariate in the null model */
double nullRss0(double *pheno, int n_ind);

/* galtRss - calculate RSS for full model in general scan */
double galtRss(double *pheno, int n_ind, int *n_gen, int n_qtl, 
	       int **Draws, double **Cov, int n_cov, int *model, 
	       int n_int, double *dwork, int *iwork, int sizefull,
	       int get_ests, double *ests, double **Ests_covar,
	       int save_design, double *designmat, int *matrix_rank);

/* end of fitqtl_imp.h */
