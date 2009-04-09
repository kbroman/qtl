/**********************************************************************
 * 
 * discan_covar.h
 *
 * copyright (c) 2004, Karl W Broman
 *
 * last modified Dec, 2004
 * first written Dec, 2004
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
 * These functions are for performing a genome scan with a binary 
 * trait and a single QTL model in the presence of covariates.
 *
 * Contains: discan_covar, discan_covar_em, discan_covar_loglik,
 *           R_discan_covar
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * discan_covar
 *
 * Performs genome scan using interval mapping in the presence of
 * covariates.  (The multipoint genotype probabilities have already 
 * been calculated in calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
 *              indexed as Genoprob[gen][pos][ind]
 *
 * Addcov       Matrix of additive covariates indexed as 
 *              Addcov[cov][ind]
 *
 * n_addcov     Number of columns in Addcov
 *
 * Intcov       Matrix of interactive covariates indexed as 
 *              Intcov[cov][ind]
 *
 * n_intcov     Number of columns in Intcov
 *
 * pheno        Phenotype data (0/1), as a vector
 *
 * start        Starting values; vector of length 
 *              n_gen + n_addcov + n_intcov*(n_gen-1)
 *
 * result       Result vector of length n_pos; upon return, contains 
 *              the LOD scores.
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * verbose      If 1, print out log likelihood at each iteration
 *
 **********************************************************************/

void discan_covar(int n_ind, int n_pos, int n_gen, 
		  double ***Genoprob, double **Addcov, int n_addcov,
		  double **Intcov, int n_intcov, int *pheno, 
		  double *start, double *result, int maxit, double tol, 
		  int verbose);

void R_discan_covar(int *n_ind, int *n_pos, int *n_gen, 
		    double *genoprob, double *addcov, int *n_addcov,
		    double *intcov, int *n_intcov, int *pheno, 
		    double *start, double *result, int *maxit, double *tol, 
		    int *verbose);

/**********************************************************************
 * 
 * discan_covar_em
 *
 **********************************************************************/

double discan_covar_em(int n_ind, int pos, int n_gen, int n_par,
		       double ***Genoprob, double **Addcov, int n_addcov,
		       double **Intcov, int n_intcov, int *pheno, 
		       double *start, int maxit, double tol, int verbose);


/**********************************************************************
 * 
 * discan_covar_loglik
 *
 **********************************************************************/

double discan_covar_loglik(int n_ind, int pos, int n_gen, int n_par,
			   double *par, 
			   double ***Genoprob, double **Addcov, int n_addcov,
			   double **Intcov, int n_intcov, int *pheno);

/* end of discan_covar.h */

