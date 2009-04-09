/**********************************************************************
 * 
 * scanone_imp.h
 *
 * copyright (c) 2001-6, Karl W Broman and Hao Wu
 *
 * This file is written by Hao Wu
 * with slight modifications by Karl Broman.
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
 * single QTL model by imputation.  
 *
 * Contains: R_scanone_imp, scanone_imp, nullRss, altRss
 *
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanone_imp
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls scanone_imp.
 * 
 **********************************************************************/

void R_scanone_imp(int *n_ind, int *n_pos, int *n_gen, int *n_draws, 
		   int *draws, double *addcov, int *n_addcov, 
		   double *intcov, int *n_intcov, double *pheno, 
		   int *nphe, double *weights,
		   double *result);

/**********************************************************************
 * 
 * scanone_imp
 *
 * Performs genome scan using the pseudomarker algorithm (imputation) 
 * method of Sen and Churchill (2001).
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * n_draws      Number of impiutations
 *
 * Draws        Array of genotype imputations, indexed as 
 *              Draws[repl][mar][ind]
 *
 * Addcov	Additive covariates matrix, Addcov[mar][ind]
 *
 * n_addcov     Number of additive covariates
 *
 * Intcov	Interacting covariates matrix, Intcov[mar][ind]
 *
 * n_intcov     Number of interacting covariates
 *
 * pheno        Phenotype data, as a vector
 *
 * weights      Vector of positive weights, of length n_ind
 *
 * Result       Matrix of size [n_pos x nphe]; upon return, contains
 *              the "LPD" (log posterior distribution of QTL location).
 * 
 **********************************************************************/

void scanone_imp(int n_ind, int n_pos, int n_gen, int n_draws, 
		 int ***Draws, double **Addcov, int n_addcov, 
		 double **Intcov, int n_intcov, double *pheno, 
		 int nphe, double *weights,
		 double **Result);

/* function to calculate the null model RSS for scanone_imp */
void nullRss(double *tmppheno, double *pheno, int nphe, int n_ind,
             double **Addcov, int n_addcov, double *dwork_null,
             int multivar, double *rss0, double *weights);

/* function to calculate the alternative model RSS. 
   This function is called by scanone_imp */
void altRss1(double *tmppheno, double *pheno, int nphe, int n_ind, int n_gen,
	     int *Draws, double **Addcov, int n_addcov, double **Intcov,
	     int n_intcov, double *dwork, int multivar, double *rss, 
	     double *weights);

/* end of scanone_imp.h */
