/**********************************************************************
 * 
 * discan.h
 *
 * copyright (c) 2001-6, Karl W Broman
 *
 * last modified Feb, 2006
 * first written Oct, 2001
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
 * single QTL model
 *
 * Contains: R_discan_mr, discan_mr,
 *           R_discan_im, discan_im
 *  
 **********************************************************************/

/**********************************************************************
 * 
 * R_discan_mr
 *
 * Wrapper for call from R; reorganizes genotype and result matrix
 * and calls discan_mr.
 * 
 **********************************************************************/

void R_discan_mr(int *n_ind, int *n_pos, int *n_gen,
		    int *geno, int *pheno, double *result);

/**********************************************************************
 * 
 * R_discan_im
 *
 * Wrapper for call from R; reorganizes genotype prob and result matrix
 * and calls discan_im.
 * 
 **********************************************************************/

void R_discan_im(int *n_ind, int *n_pos, int *n_gen, 
		 double *genoprob, int *pheno, double *result, 
		 int *maxit, double *tol);

/**********************************************************************
 * 
 * discan_mr
 *
 * Performs genome scan using marker regression for a dichotomous trait
 * 
 * n_ind        Number of individuals
 *
 * n_pos        Number of marker positions
 *
 * n_gen        Number of different genotypes
 *
 * Geno         Genotype matrix
 *
 * pheno        Phenotype data, as a vector
 *
 * result       Upon return, to contain the log10 likelihoods
 *
 * means        Space for the phenotype means for each genotype
 *
 **********************************************************************/

void discan_mr(int n_ind, int n_pos, int n_gen, int **Geno, 
		  int *pheno, double *result, double *means);

/**********************************************************************
 * 
 * discan_im
 *
 * Performs genome scan using interval mapping.  (The multipoint
 * genotype probabilities have already been calculated in 
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
 * pheno        Phenotype data, as a vector
 *
 * result       Upon return, to contain the log10 likelihoods
 *
 * maxit        Maximum number of iterations in the EM algorithm
 *
 * tol          Tolerance for determining convergence in EM
 *
 * work         Workspace of length n_gen
 *
 * means        Space for the phenotype means for each genotype
 *
 **********************************************************************/

void discan_im(int n_ind, int n_pos, int n_gen, double ***Genoprob,
	       int *pheno, double *result, 
	       int maxit, double tol, double **work, double *means);

/* end of discan.h */

