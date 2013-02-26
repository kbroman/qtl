/**********************************************************************
 * 
 * fitqtl_hk.h
 *
 * copyright (c) 2007-2013, Karl W Broman
 *
 * last modified Feb, 2013
 * first written Nov, 2007
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
 * Haley-Knott regression.
 *
 * Contains: R_fitqtl_hk, fitqtl_hk, galtRssHK
 *
 **********************************************************************/
void R_fitqtl_hk(int *n_ind, int *n_qtl, int *n_gen, 
		 double *genoprob, int *n_cov, double *cov, int *model, 
		 int *n_int, double *pheno, int *get_ests,
		  /* return variables */
		 double *lod, int *df, double *ests, double *ests_covar,
		 double *design_mat, int *matrix_rank);

/**********************************************************************
 * 
 * fitqtl_hk
 *
 * Fits a fixed multiple-QTL model by Haley-Knott regression.
 * 
 * n_ind        Number of individuals
 *
 * n_qtl        Number of QTLs in the model 
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     QTL genotype probabilities
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
 * matrix_rank  On return, rank of design matrix
 *
 **********************************************************************/

void fitqtl_hk(int n_ind, int n_qtl, int *n_gen, double ***Genoprob,
	       double **Cov, int n_cov, 
	       int *model, int n_int, double *pheno, int get_ests,
	       double *lod, int *df, double *ests, double *ests_covar,
	       double *design_mat, int *matrix_rank); 

/* galtRssHK - calculate RSS for full model by Haley-Knott regression */
double galtRssHK(double *pheno, int n_ind, int *n_gen, int n_qtl, 
		 double ***Genoprob, double **Cov, int n_cov, int *model, 
		 int n_int, double *dwork, int *iwork, int sizefull,
		 int get_ests, double *ests, double **Ests_covar,
		 double *designmat, int *matrix_rank);

/* end of fitqtl_hk.h */
