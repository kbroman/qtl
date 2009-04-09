/**********************************************************************
 * 
 * scanone_ehk.h
 *
 * copyright (c) 2006, Karl W Broman
 *
 * last modified Jul, 2006
 * first written Jul, 2006
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
 * These functions are for performing a genome scan with a 
 * single QTL model by the extended Haley-Knott method
 *
 * Contains: R_scanone_ehk, scanone_ehk, calc_mvz
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_ehk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_ehk.
 * 
 **********************************************************************/

void R_scanone_ehk(int *n_ind, int *n_pos, int *n_gen,
		   double *genoprob, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno,
		   double *weights, double *result, int *maxit,
		   double *tol);

/**********************************************************************
 * 
 * scanone_ehk
 *
 * Performs genome scan using the extended Haley-Knott method
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              Indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates: Addcov[cov][ind]
 * 
 * n_addcov     Number of columns of Addcov
 *
 * Intcov       Number of interactive covariates: Intcov[cov][ind]
 *
 * n_intcov     Number of columns of Intcov
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * result       Vector of length n_pos, to contain the LOD scores
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 **********************************************************************/

void scanone_ehk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		 double **Addcov, int n_addcov, double **Intcov, 
		 int n_intcov, double *pheno, double *weights, 
		 double *result, int maxit, double tol);

/* calc_mvz */
void calc_mvz(int n_ind, int curpos, int n_gen, double ***Genoprob, 
	      double **Addcov, int n_addcov, double **Intcov, int n_intcov, 
	      double *pheno, double *weights, double *coef, double sigmasq,
	      double *m, double *v, double *z);

/* end of scanone_ehk.h */
