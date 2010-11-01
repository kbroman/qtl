/**********************************************************************
 * 
 * scanone_hk_binary.h
 *
 * copyright (c) 2010, Karl W Broman
 *
 * last modified Jul, 2010
 * first written Jun, 2010
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
 * single QTL model by Haley-Knott regression
 *
 * Contains: R_scanone_hk_binary, scanone_hk_binary
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_hk_binary
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk_binary(int *n_ind, int *n_pos, int *n_gen,
			 double *genoprob, double *addcov, int *n_addcov, 
			 double *intcov, int *n_intcov, double *pheno,
			 double *result, double *tol, int *maxit, 
			 int *verbose, int *ind_noqtl);

/**********************************************************************
 * 
 * scanone_hk_binary
 *
 * Performs genome scan using the Haley-Knott regression method
 * for a binary trait
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
 * result       vector of length n_ind, to contain the log10 likelihood values
 *
 * tol          tolerance for convergence
 *
 * maxit        maximum number of iterations
 *
 * verbose      if TRUE, give some output
 *
 * ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *              from QTL effects.  
 *
 **********************************************************************/

void scanone_hk_binary(int n_ind, int n_pos, int n_gen, double ***Genoprob,
		       double **Addcov, int n_addcov, double **Intcov, 
		       int n_intcov, double *pheno, 
		       double *result, double tol, int maxit, int verbose,
		       int *ind_noqtl);

/* end of scanone_hk_binary.h */
