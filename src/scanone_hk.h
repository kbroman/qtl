/**********************************************************************
 * 
 * scanone_hk.h
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Feb, 2006
 * first written Nov, 2001
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
 * Contains: R_scanone_hk, scanone_hk
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_hk
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_hk.
 * 
 **********************************************************************/

void R_scanone_hk(int *n_ind, int *n_pos, int *n_gen,
		  double *genoprob, double *addcov, int *n_addcov, 
                  double *intcov, int *n_intcov, double *pheno, int *nphe,
		  double *weights, double *result);

/**********************************************************************
 * 
 * scanone_hk
 *
 * Performs genome scan using the Haley-Knott regression method
 * (regressing phenotypes on conditional genotype probabilities; the
 * multipoint genotype probabilities have already been calculated in
 * calc.genoprob)
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Genoprob     Array of conditional genotype probabilities
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
 * nphe         Number of phenotypes
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Result matrix of size [n_pos x (nphe)] containing the
 *              LOD scores for each phenotype
 *
 **********************************************************************/

void scanone_hk(int n_ind, int n_pos, int n_gen, double ***Genoprob,
                double **Addcov, int n_addcov, double **Intcov, 
		int n_intcov, double *pheno, int nphe, double *weights,
		double **Result);

/* end of scanone_hk.h */
